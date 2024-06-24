import functools
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import os
import pandas as pd

from analysispkg import pkg_analyze_cast
from analysispkg import pkg_utils

"""
Score comparisons

The scoring comparison data is high-dimensional. For instance, water level score comparisons are a function of
(station,period,scoretype,model,wltype) which is 5D, and for CTD we have (station,period,scoretype,model,variable,region)
which is 6D. Meanwhile, tables are 2D, so the data needs to be "flattened" to be presentable. Three strategies seem to work:

S1/ We can arrange 3D data into a series of 2D tables by 'looping' over a particular dimension. So 3D data with
    dimensions (A,B,C) data becomes C tables of shape A,B

S2/ We can "pivot" / "stack" dimensions into columns to effectively collapse a dimension. The (A,B,C) example would become
              c1                c2
         b1   b2   b3      b1   b2   b3
    a1  v111 v121 v131    v112 v122 v132
    a2  v211 v221 v231    ...

    This one works fairly well when the B or C dimensions are fairly short. For ex, model dimension as "B" here works as
    it's 2-4 in length. Similarly "variable" works here, as it may be just 2 (temp and salt). You can apply this multiple
    times to stack higher, and there is a need to choose the ordering of the stacking.

S3/ We can aggregate a dimension to eliminate it, such as by averaging. This only works for some dimensions -- we can
    average over stations, but not over variables/scores/models.


The code below uses a combination of these strategies to make useful tables. It attempts to be systematic:

  1) Assemble all scoring data into a single pandas dataframe, and reduce it to a shorter set of stations to present.
   This can be (a) a strict set-intersection across models, or (b) that set expanded a bit, where there will be gaps
   gaps for missing stations. See reduce_to_common_stations() below and pkg_utils.get_common_stations() for details.

  2) Apply strategy 1 to produce the "full" tables of scores by looping over period first, then region or scoretype etc,
     as needed to get the full tables.

  3) Produce some pivot tables by applying strategy 2, with careful choices about which dimensions to stack. Here we use
     the pandas pivot() function, and choose the index dimension (rows) and the column stacking order to get a useful table.
     Sorting is alphabetical by default, and we ensure "casename" columns appear in the order provided in the comparison YAML.

  4) Produce some aggregated pivot tables. Same as above but we average a dimension (typically Station) and then make
     choices about index and stacking columns. This uses the pandas pivot_table() function.

Lastly, we save scores to CSV, TEX, and XLSX spreadsheet files (see default_comparison.yaml for switches to control this).
The spreadsheets are a means of assembling strategy 1 tables into a single file by putting each table on it's own sheet.
"""


def plot_score(df, opt, mods, typ, period, score):
    """ Plot of the scores recorded in score comparison CSV files.

    Parameters
    ----------
    df : Pandas table
        Scores.
    opt : dict
        Options.
    mods: list
        List of model options dicts
    typ :
        Instrument type.
    period : str
        Label for the analyzed period.
    score : str
        Score to process.
    """

    # generate series names: shorten the original ones so xlabels are not too long
    def make_xlabel(stn):
        if typ == 'TG':
            r = stn
        else:
            tmp = stn.split('_' + typ + '_')         # split stn name and the rest, skip instr type
            r = tmp[0] + ' ' + tmp[1].split('_')[0]  # from the second part, take start date
        return r
    stn2xlabel = {stn:make_xlabel(stn) for stn in set(df['Station name'])}

    nstn = len(selectdrop(df, {'casename':mods[0]['casename']}))
    figwidth = max(8.5, nstn * 0.2)  # big enough to resolve stn names along the axis, but not too small
    fig, ax = plt.subplots(figsize=(figwidth, 0.7 * figwidth))

    for mod,modopt in zip(mods,opt['mods']):
        s = selectdrop(df, {'casename':mod['casename']})
        stns,scor,dep = s['Station name'], s['score'], s['Depth']
        xlabels = [stn2xlabel[stn] + ' ' + str(round(z)) + 'm' for stn, z in zip(stns, df['Depth'])]
        plt.plot(xlabels,scor, '.-', color=modopt['color'], label=mod['casename'])
    ax.set_xticks(np.arange(nstn))
    ax.set_xticklabels(xlabels, rotation=270)
    ax.grid()
    ax.legend()
    ax.set_ylabel(score)  #TODO Add units
    ax.set_title('{} {} {}'.format(typ, period, score))
    fig.tight_layout()
    figname = os.path.join(opt['dir_plots'], typ, period, '_'.join(('score', score + '.png')))
    fig.savefig(figname, dpi=100)
    plt.close(fig=fig)


def selectdrop(df, conditions):
    """
    Helper function to select exact matches for conditions (dict), and we return the results
    with those columns dropped
    """
    # create masks for each condition
    masks = [df[key] == val for key,val in conditions.items()]
    # logical-and them all together to get final mask
    mask = functools.reduce(lambda a,b: a & b, masks)
    # use the mask to select rows and then drop the matching-condition columns
    r = df[mask].drop(columns=conditions.keys())
    return r


def reduce_to_common_stations(opt, df, column='Station name', searchcolumn='casename'):
    """
    Helper function to reduce stations (see pkg_utils.get_common_stations for details about intersection & union)
    """
    if len(df) == 0:
        return df
    # Determine unique values in searchcolumn
    uniq = set(df[searchcolumn].to_list())
    # Get data from column for each unique value
    cols = [set(df[df[searchcolumn] == x][column].to_list()) for x in uniq]
    # Get the common stations
    r = pkg_utils.get_common_stations(opt, uniq, cols)
    # Retain only these common stations
    df = df[df[column].isin(r)]
    return df


def savecsv(df, filename, float_format='%.03f'):
    """
    Wrapper to reduce boilerplate and specify some defaults
    """
    dirname, fname = os.path.split(filename)
    os.makedirs(dirname, exist_ok=True)
    df.to_csv(filename, float_format=float_format)


def savetex(df, filename, float_format='%.03f'):
    """
    Wrapper to reduce boilerplate and specify some defaults
    """
    dirname, fname = os.path.split(filename)
    os.makedirs(dirname, exist_ok=True)
    with open(filename,'w') as f:
        df.to_latex(buf=f, float_format=float_format)


def savexls(df, filename, sheet_name, float_format='%.03f', formatter=None):
    """
    Wrapper to save a dataframe as a new sheet in a spreadsheet file
    """
    if len(df) == 0:
        print("savexls: received empty df, not writing sheet {} to {}".format(sheet_name,filename))
        return
    def columns_best_fit(ws):
        # Based on https://stackoverflow.com/questions/13197574/openpyxl-adjust-column-width-size/57634315#57634315
        def as_text(value):
            return "" if value is None else str(value)
        for column_cells in ws.columns:
            new_column_length = max(len(as_text(cell.value)) for cell in column_cells)
            new_column_letter = (openpyxl.utils.get_column_letter(column_cells[0].column))
            if new_column_length > 0:
                ws.column_dimensions[new_column_letter].width = new_column_length + 4

    dirname, fname = os.path.split(filename)
    os.makedirs(dirname, exist_ok=True)
    if os.path.exists(filename):
        writer = pd.ExcelWriter(filename, mode='a', engine='openpyxl', if_sheet_exists='replace')
    else:
        writer = pd.ExcelWriter(filename, mode='w', engine='openpyxl')
    df.to_excel(writer, sheet_name=sheet_name, float_format=float_format)
    worksheet = writer.sheets[sheet_name]
    columns_best_fit(worksheet)
    if formatter is not None:
        formatter(worksheet)
    writer.close()


def savescores(opt, df, out_dir, filename_part1, filename_part2, float_format='%.03f', formatter=None):
    """
    Helper function to write the scores
    For CSV and TEX, filename_part1 and filename_part2 get combined to form the filename
    For XLSX, filename_part1 is the filename and filename_part2 is the sheet name
    formatter : function, optional
        Custom function for xlsx formatting
    """
    if opt['score']['save_csv']:
        csvfile = os.path.join(out_dir, filename_part1 + "_" + filename_part2 + ".csv")
        savecsv(df, csvfile, float_format=float_format)

    if opt['score']['save_tex']:
        texfile = os.path.join(out_dir, filename_part1 + "_" + filename_part2 + ".tex")
        # savetex(df, texfile, float_format=float_format)
        try:
            savetex(df, texfile, float_format=float_format)
        except:
            print(f'***** Failed to save tex file {texfile}')
            # TODO Figure out the problem. Occurs in currents_residual
            #   MD:  some debugging: if I change the type of the data from 'object' to 'float' it fixes it, sort of
            #   eg, df['crmse'] = df['crmse'].astype(float)
            #   however that doesn't work for the non-float complex numbers, so that's not a solution that we can just
            #   stick in the data loading part...

    if opt['score']['save_xlsx']:
        xlsfile = os.path.join(out_dir, filename_part1 + ".xlsx")
        savexls(df, xlsfile, sheet_name=filename_part2, float_format=float_format, formatter=formatter)

def get_case_order(df, mods):
    casenames = [mod['casename'] for mod in mods]
    casenames_present = set(df['casename'].to_list())
    case_order = [x for x in casenames if x in casenames_present]
    return case_order

def ts_timeseries(opt, mods, typ):
    """ Output comparison tables for T & S time series for LH and SST.

    Data here is 5D: (station,period,scoretype,model,variable); for SST the variable dimension has length 1

    typ : {'LH' | 'SST'}
    """
    scores = sorted(['rmse', 'mad', 'mae', 'gamma2', 'bias', 'pearson', 'skill1981', 'crmse'])
    periods = opt['analysis']['periods'].keys()
    if typ in ['LH']:
        vars=['temperature', 'salinity']
    elif typ in ['SST']:
        vars = ['temperature']
    else:
        print("pkg_csv_comparison.ts_timeseries: unrecognized type: ", typ)
        return

    def loadscores():
        # Data loading loop
        df = pd.DataFrame()
        for period in periods:
            for score in scores:
                for var in vars:
                    for mod in mods:
                        indir = pkg_utils.get_plot_dir(mod, typ, period)
                        fn = os.path.join(indir, pkg_utils.csv_file(typ, period, 'scores_' + var, score, mod))
                        if not os.path.exists(fn):
                            print("Can not find {}".format(fn))
                            continue
                        dfi = pd.read_csv(fn, usecols=['Station name', score]).rename(columns={score:'score'})
                        dfi['period'] = period
                        dfi['var'] = var
                        dfi['scoretype'] = score
                        dfi['casename'] = mod['casename']
                        df = pd.concat([df, dfi], ignore_index=True)
        df = reduce_to_common_stations(opt, df, 'Station name', 'casename')
        return df

    df = loadscores()
    case_order = get_case_order(df, mods)
    if len(df) == 0:
        return

    # Loop over periods and scores
    for period in periods:
        out_dir = pkg_utils.get_plot_dir(opt, typ, period)
        for score in scores:
            if typ in ['SST']:
                s = selectdrop(df, {'period':period,'scoretype':score,'var':'temperature'}).rename(columns={'score':score})
                ss = s.pivot(index='Station name', columns=['casename']).reindex(case_order, axis=1, level=1) # (station, casename)
            else:
                s = selectdrop(df, {'period':period,'scoretype':score}).rename(columns={'score':score})
                ss = s.pivot(index='Station name', columns=['var','casename']).reindex(case_order, axis=1, level=2) # (station, variable&casename)
            savescores(opt, ss, out_dir, typ.lower() + '_' + period, score)


    print('ts_timeseries for ' + typ + ' comparison scores written.')


def ferry(opt, mods):
    """ Output comparison tables for FERRY.

    Reads the CSVs generated from pkg_csv and uses pandas for merging
    """
    typ = "FERRY"
    domains = set.union(*[ set(mod['analysis']['domains']) for mod in mods])
    scores = sorted(['rmse', 'mad', 'mae', 'gamma2', 'bias', 'pearson', 'skill1981', 'crmse'])
    periods = opt['analysis']['periods'].keys()
    vars = ['T', 'S']

    # Data loading loop
    def loaddata():
        df = pd.DataFrame()
        for period in periods:
            for score in scores:
                for var in vars:
                    for mod in mods:
                        indir = pkg_utils.get_plot_dir(mod, typ, period)
                        fn = os.path.join(indir, pkg_utils.csv_file(typ, period, 'scores_' + var, score, mod))
                        if not os.path.exists(fn):
                            print("Can not find {}".format(fn))
                            continue
                        dfi = pd.read_csv(fn, usecols=['Ferry name', score, 'Domain']).rename(columns={score: 'score'})
                        dfi['period'] = period
                        dfi['var'] = var
                        dfi['scoretype'] = score
                        dfi['casename'] = mod['casename']
                        df = pd.concat([df, dfi], ignore_index=True)

        df = reduce_to_common_stations(opt, df, 'Ferry name', 'casename')
        return df

    df = loaddata()

    # Loop over periods, domains, and scores
    for period in periods:
        out_dir = pkg_utils.get_plot_dir(opt, typ, period)
        for domain in domains:
            tmp = selectdrop(df, {'period': period, 'Domain': domain})
            if np.any(tmp['score']): # proceed only if we have scores for this domain
                for score in scores:
                    s = selectdrop(df, {'period':period,'Domain':domain,'scoretype':score}).rename(columns={'score':score})
                    ss = s.pivot(index='Ferry name', columns=['casename','var']) # (station, casename&variable)
                    savescores(opt, ss, out_dir, typ.lower() + '_' + period + '_' + domain, score)

    # One big score table per domain
    out_dir = pkg_utils.get_plot_dir(opt, typ)
    for domain in domains:
        s = selectdrop(df, {'Domain': domain})
        if np.any(s['score']): # proceed only if we have scores for this domain
            ss = s.pivot(index=['scoretype', 'Ferry name', 'period'], columns=['var', 'casename'])
            savescores(opt, ss, out_dir, typ.lower() + '_' + domain, "allscores")

    # Per ferry score table for ACOM
    scores = ['bias', 'crmse']
    df = loaddata()
    ferries = set(df['Ferry name'].to_list())
    for domain in domains:
        for ferry in ferries:
            s = selectdrop(df, {'Domain': domain, 'Ferry name': ferry})
            if np.any(s['score']): # proceed only if we have scores for this domain
                ss = s.pivot(index=['scoretype', 'period'], columns=['var', 'casename'])
                savescores(opt, ss, out_dir, 'ACOM_' + typ.lower() + '_' + domain, ferry)

    print('FERRY comparison scores written.')


def ctd(opt, mods):
    """ Output comparison tables for CTD.

    Reads the CSVs generated from pkg_csv and uses pandas for merging
    The data here is 6D: (station,period,scoretype,variable,region,casename)
    """
    typ = "CTD"
    periods = opt['analysis']['periods'].keys()
    scores = ['bias', 'rmse', 'crmse', 'skill1981']
    vars = ['T', 'S']

    # Data loading loop
    def loadscores():
        df = pd.DataFrame()
        for period in periods:
            for score in scores:
                for var in vars:
                    for mod in mods:
                        indir = pkg_utils.get_plot_dir(mod, typ, period)
                        fn = os.path.join(indir, pkg_utils.csv_file(typ, period, 'scores', score + '_' + var, mod))
                        if not os.path.exists(fn):
                            print("Can not find {}, done".format(fn))
                            continue
                        dfi = pd.read_csv(fn, usecols=['Station name', 'Region', score]).rename(columns={score:'score'})
                        dfi['period'] = period
                        dfi['var'] = var
                        dfi['scoretype'] = score
                        dfi['casename'] = mod['casename']
                        df = pd.concat([df, dfi], ignore_index=True)

        df = reduce_to_common_stations(opt, df, 'Station name', 'casename')
        return df

    df = loadscores()
    regions = set(df['Region'])

    # Loop over periods, regions and scores
    for period in periods:
        out_dir = pkg_utils.get_plot_dir(opt, typ, period)
        for region in regions:
            for score in scores:
                s = selectdrop(df, {'period': period, 'scoretype': score, 'Region':region}).rename(columns={'score': score})
                ss = s.pivot(index='Station name', columns=['var', 'casename'])
                savescores(opt, ss, out_dir, 'CTD_' + period + '_' + region, score)

    # ACOM scores
    scores = ['bias', 'crmse']
    df = loadscores()
    case_order = get_case_order(df, mods)
    _, region_defs = pkg_analyze_cast.read_CTD_domain(opt['CTD']['domain_file'])
    region_order = region_defs.keys()
    # Loop over periods, regions and scores
    for period in periods:
        print(period)
        out_dir = pkg_utils.get_plot_dir(opt, typ, period)
        s0 = selectdrop(df, {'period':period}).rename(columns={'scoretype':'score','score':'data'})
        s1 = s0.pivot_table(index=['Region','score'], columns=['var','casename'], values='data', aggfunc='mean')
        s2 = s1.reindex(vars, axis=1, level=0).reindex(case_order, axis=1, level=1).reindex(region_order, axis=0, level=0)
        savescores(opt, s2, out_dir, 'ACOM', period)

        ncasts = []
        for region in region_order:
            tmp = selectdrop(df, {'period': period, 'Region':region}).rename(columns={'scoretype': 'score', 'score': 'data'})
            ncasts += [len(set(tmp['Station name']))]
        s0 = pd.DataFrame(data={'Region':region_order, 'Number of casts': ncasts})
        savescores(opt, s0, out_dir, 'ACOM_casts_per_region', period)

    # TODO: mean is not appropriate with intersection+union, needs to be reworked to just intersectio if we want to keep this
    # # Loop over region, mean-over-stations
    # out_dir = pkg_utils.get_plot_dir(opt, typ)
    # for region in regions:
    #     s = selectdrop(df, {'Region':region})
    #
    #     ss = s.pivot_table(index='period', columns=['scoretype', 'var', 'casename'], aggfunc='mean')
    #     savescores(opt, ss, out_dir, 'CTD_' + region, "by_period")
    #
    #     ss = s.pivot_table(index='casename', columns=['scoretype', 'var', 'period'], aggfunc='mean')
    #     savescores(opt, ss, out_dir, 'CTD_' + region, "by_model")

    print('CTD comparison scores written.')




