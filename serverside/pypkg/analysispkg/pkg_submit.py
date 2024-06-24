# Functionality for running tasks as submitted jobs
import os
import socket
import subprocess
import time

from analysispkg import pkg_utils


def get_systems():
    systems_gpsc = {
        'gpsc7': {
            'cores_per_node': 64,
            'memorymb': 460000,
            'hours': 6,
            'submitter': 'edc',
            'host': 'gpsc7.science.gc.ca',
        }
    }
    systems_sharcnet = {
        'sharcnet': {
            'cores_per_node': 32,
            'memorymb': 120000,
            'hours': 6,
            'submitter': 'slurm'
        }
    }

    hostname=socket.getfqdn()
    if "science.gc.ca" in hostname:  # Assume GPSC
        return systems_gpsc
    elif "sharcnet" in hostname:  # SHARCNET
        return systems_sharcnet
    else: # no systems to submit to
        return {}

def get_sysnames():
    systems = get_systems() # dictionary
    sysnames = list(systems.keys()) # list
    return sysnames


def submit(sysname, myargs, logdir, cores=None, nodes=None, hours=None, jobnumbers=None, mpi=False, ppn=None):
    """
    Submits a task to the job queue.

    - sysname is the name of the system to submit to
    - myargs is the command to run in list-of-strings form, eg:
     ['extract.py','configfile.yaml','-i','CTD','-j','44']
    - cores, nodes, hours are job resource parameters, provided as integers
    - jobnumbers is a list-of-strings that this job should wait for before running, eg
     ['473762','473763']
    - mpi is a flag to specify that this job is an MPI job
    - ppn is the number of processes per node for mpi job (set to 1 for tidal diagnostic, don't use otherwise)

    This function returns the job number to facilitate job chaining.
    """
    system = get_systems()[sysname]
    submitter = system['submitter']
    cores_per_node = system['cores_per_node']
    memorymb = system['memorymb']

    # Allow for partial node use (eg for plots or scores that doesn't need a full node)
    if nodes is not None and 0 < nodes < 1:
        frac = nodes
        nodes = 1
        if cores is None:
            cores = max(int(cores_per_node*frac),1)
        memorymb = int(memorymb*frac)

    if hours is None: # Use default runtime limit
        hours = system['hours']

    if cores is None: # If cores not provided, we use full nodes
        cores = cores_per_node

    if nodes is None: # If nodes not provided, we use one node, and clip #cores at node size
        nodes = 1
        cores = min(cores, cores_per_node)

    # This is needed to ensure we override whatever is in the YAML
    if mpi is False:
        myargs += ['-j',str(cores)]

    print("Submission params: nodes={}, cores={}, memory={}mb, hours={}".format(nodes,cores,memorymb,hours))

    datestamp = time.strftime('%Y-%m-%dH%H-%M-%S')
    os.makedirs(logdir, exist_ok=True)

    _, JOBNAME = os.path.split('_'.join(myargs))
    CORESPERNODE = str(cores)
    NNODES = str(nodes)
    MEMORYMB = str(memorymb)
    HOURS = str(hours)
    STDOUTLOG = os.path.join(logdir, datestamp + '_' + JOBNAME + '.log')

    PKGBIN = os.path.normpath(os.path.dirname(__file__) + "/../bin")
    CONFIG_WD  = os.getcwd()

    # if "analyze.py" in ''.join(myargs) and "TDL" in myargs:
    #     # Special case for analyzing TDL
    #     # t_tide is threaded via numpy, so it uses 44 cores efficiently with one process
    #     # So we can use mpirun here with one process per node to make full use of >1 nodes

    if mpi:
        if ppn is None:
            PPN = CORERSPERNODE
        else:
            PPN = str(ppn)
        if submitter == 'gpsc':
            JOBCMD = ' '.join(['rumpirun', '-np', NNODES, "-ppn", PPN] + myargs)
        else:
            JOBCMD = ' '.join(['mpirun', '-np', NNODES, "-ppn", PPN] + myargs)
    else:
        JOBCMD = ' '.join(myargs)

    if submitter == 'gpsc':
        PROJECT = pkg_utils.getuservar('project')
        CONDA_SCRIPT = pkg_utils.getuservar('CONDA_SCRIPT')
        templatefile = os.path.join(PKGBIN, "template_gpsc.sh")
    elif submitter == 'edc':
        PROJECT = pkg_utils.getuservar('project').replace('-','_')
        CONDA_SCRIPT = pkg_utils.getuservar('CONDA_SCRIPT')
        templatefile = os.path.join(PKGBIN, "template_edc.sh")
    elif submitter == 'slurm':
        PROJECT = 'none'
        CONDA_SCRIPT = 'none'
        templatefile = os.path.join(PKGBIN, "template_slurm.sh")
    else:
        print("Unknown submission system, unable to choose job submission template, quitting")
        raise

    with open(templatefile) as f:
        script = f.read() \
                    .replace("JOBNAME", JOBNAME) \
                    .replace("PROJECT", PROJECT) \
                    .replace("NNODES", NNODES) \
                    .replace("CORESPERNODE", CORESPERNODE) \
                    .replace("MEMORYMB", MEMORYMB) \
                    .replace("HOURS", HOURS) \
                    .replace("STDOUTLOG", STDOUTLOG) \
                    .replace("CONDA_SCRIPT", CONDA_SCRIPT) \
                    .replace("CONFIG_WD", CONFIG_WD) \
                    .replace("JOBCMD", JOBCMD)

    scriptfile = os.path.join(logdir, datestamp + '_' + JOBNAME + ".sh")
    with open(scriptfile,'w') as f:
        f.write(script)

    if submitter == 'gpsc':
        # Base submit command
        subcmd = ['jobsub', '-c', system['host'], scriptfile]

        # If jobnumbers provided, we should hold for their completion
        if jobnumbers is not None:
            subcmd += ['--', '-hold_jid', ','.join(jobnumbers)]

        print("Submitting job: ", ' '.join(subcmd))
        result = subprocess.run(subcmd,capture_output=True)

        # Extract job number
        jobnumber = str(result.stdout).split("Your job ")[1].split()[0]

    elif submitter == 'edc':
        # check whether it's a vis node or interactive
        isvis = subprocess.run(['uname', '-n'], capture_output=True)
        if 'vis' in isvis.stdout.decode('utf-8'):
            subcmd = ['ssh', 'inter-dfo-ubuntu1804.science.gc.ca']
        else: 
            subcmd = []

        if jobnumbers is not None:
            subcmd.extend(['sbatch', '-d', ','.join(jobnumbers), scriptfile])
        else:
            subcmd.extend(['sbatch', scriptfile])

        print("Submitting job: ", ' '.join(subcmd))
        result = subprocess.run(subcmd, capture_output=True)

        # Extract job number
        jobnumber = result.stdout.decode('utf-8').split()[-1]


    elif submitter == 'slurm' : #or submitter == 'edc':
        if jobnumbers is not None:
            subcmd = ['sbatch', '-d', ','.join(jobnumbers), scriptfile]
        else:
            subcmd = ['sbatch', scriptfile]

        print("Submitting job: ", ' '.join(subcmd))
        result = subprocess.run(subcmd, capture_output=True)

        # Extract job number
        jobnumber = result.stdout.decode('utf-8').split()[-1]

    return [jobnumber]
