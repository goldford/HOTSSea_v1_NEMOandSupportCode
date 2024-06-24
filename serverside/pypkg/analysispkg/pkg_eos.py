import gsw
"""
    Handle the TEOS-10 <-> EOS-80 conversions
"""

def CTSA_to_PTSP(CT,SA,dep,lon,lat):
    p = gsw.p_from_z(-dep, lat)           # depth to pressure
    SP = gsw.SP_from_SA(SA, p, lon, lat)  # Practical Salinity (PSS-78) from Absolute Salinity, g/kg
    PT = gsw.pt_from_CT(SA, CT)           # Potential Temperature from Conservative Temperature
    return PT, SP

def CTSR_to_PTSP(CT,SR,dep,lon,lat):
    p = gsw.p_from_z(-dep, lat)            # depth to pressure
    SP = gsw.SP_from_SR(SR)                # Practical Salinity (PSS-78) from Reference Salinity, g/kg
    SA = gsw.SA_from_SP(SP, p, lon, lat)   # Absolute Salinity, g/kg from Practical Salinity (PSS-78)
    PT = gsw.pt_from_CT(SA, CT)            # Potential Temperature from Conservative Temperature
    return PT, SP

def eoscheck(ncf,temp_name,sal_name):
    # Detect if the model has saved TEOS-10 temp & salinity
    flag_CT, flag_SR, flag_SA = False, False, False
    ncT = ncf[temp_name]
    ncS = ncf[sal_name]
    if 'standard_name' in ncT.ncattrs():
        if 'conservative' in ncT.standard_name:
            flag_CT = True
    if 'standard_name' in ncS.ncattrs():
        if 'reference' in ncS.standard_name:
            flag_SR = True
        elif 'absolute' in ncS.standard_name:
            flag_SA = True
    return flag_CT, flag_SR, flag_SA
