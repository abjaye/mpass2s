#!/usr/bin/env python
import os, shutil, sys, argparse

from argparse              import RawTextHelpFormatter    
from datetime import timedelta, date, datetime
from ecflow import Defs, Suite, Task, Edit, Trigger, Family

#from standard_script_setup import *

def parse_command_line(args, description):
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    #CIME.utils.setup_standard_logging_options(parser)
    parser.add_argument("--date",
                        help="Specify a start Date")

    args = parser.parse_args()

    if args.date:
        try:
            date = datetime.strptime(args.date, '%Y-%m-%d')
        except ValueError as verr:
            raise ValueError("Incorrect data format, should be YYYY-MM-DD") from verr
    elif cdate:
        date = datetime.strptime(cdate, '%Y-%m-%d')
    else:
        date = datetime.today() - timedelta(days=1)

    return date.strftime("%Y_%m_%d")

def workflow(description):
    print("Creating suite definition")
    workflow = os.getenv("MPAS_WORKFLOW")
    home = os.path.join(os.getenv("HOME"),workflow,"ecflow")
    machine = os.getenv("NCAR_HOST")
    user = os.getenv("USER")
    #expect(home and workflow and machine and user,f"Missing required env variable: home={home} MPAS_WORKFLOW={workflow} machine={machine} user={user}")
    
    # user changes here
    ensemble_start = 0
    ensemble_end = 10
    project="CESM0020"

    #fcstdate="2023_09_25"
    fcstdate = parse_command_line(sys.argv, __doc__)
    
    print(f"running for {fcstdate}")
    fcsthome=os.path.join(os.getenv("HOME"),workflow)
    fcstwork=os.path.join("/glade/work/",user,machine,"cases",workflow+"_"+fcstdate)
    
    print(f"workflow={workflow} fcsthome={fcsthome} fcstwork={fcstwork}")

    workdir = os.path.join(home,workflow+"_"+fcstdate)
    if os.path.isdir(workdir):
        txt = input(f"Direcotry {workdir} exists, delete? (y or n) ")
        if txt == 'y':
            shutil.rmtree(workdir)
        else:
            sys.exit(-1)
            
    shutil.copytree(os.path.join(home,"template"), workdir)
    workflow_date = workflow+"_"+fcstdate
    init_member = os.path.join(home,workflow_date,"init_family","init_member.ecf")
    init_noah_member = os.path.join(home,workflow_date,"init_noah_family","init_noah_member.ecf")
    run_member = os.path.join(home,workflow_date,"run_family","run_member.ecf")
    pp_member =  os.path.join(home,workflow_date,"postprocess_family","postprocess_member.ecf")
    build_ens = os.path.join(home,workflow_date,"build_ens.ecf")
    for i in range(ensemble_start, ensemble_end+1):
        new_init_member = os.path.join(home,workflow_date,"init_family",f"init_member_{i:02d}.ecf")
        shutil.copy(init_member, new_init_member)
        new_init_noah_member = os.path.join(home,workflow_date,"init_noah_family",f"init_noah_member_{i:02d}.ecf")
        shutil.copy(init_noah_member, new_init_noah_member)
        new_member = os.path.join(home,workflow_date,"run_family",f"run_member_{i:02d}.ecf")
        shutil.copy(run_member, new_member)
        new_pp = os.path.join(home,workflow_date,"postprocess_family",f"postprocess_member_{i:02d}.ecf")
        shutil.copy(pp_member, new_pp)
    logdir = os.path.join(home,workflow_date,"log")
    os.mkdir(logdir)
    defs = Defs(Suite(workflow_date,
                      Edit(PROJECT=project,
                           ECF_JOB_CMD="qsub %ECF_JOB% 1>%ECF_JOBOUT% 2>&1",
                           MPAS_WORKFLOW=workflow,
                           FCSTDATE=fcstdate.replace("_","-"),
                           ENSEMBLE_START=ensemble_start,
                           ENSEMBLE_END=ensemble_end,
                           FCST_HOME=fcsthome,
                           FCST_WORK=fcstwork,
                           ECF_INCLUDE=os.path.join(home,"include"),
                           ECF_PORT=40179,
                           ECF_HOST="derecho7.hpc.ucar.edu",
                           ECF_HOME=home,
                           LOGDIR=logdir),
                      Task(f"build_ens"),
                      Family("init_family").add(
                          [Task(f"init_member_{i:02d}",
                                Trigger("/"+workflow_date + "/build_ens == complete"))
                for i in range(ensemble_start, ensemble_end+1)]),
                      Family("init_noah_family").add(
                          [Task(f"init_noah_member_{i:02d}",
                                Trigger("/" + workflow_date + f"/init_family/init_member_{i:02d} == complete"))
                           for i in range(ensemble_start, ensemble_end+1)
                for i in range(ensemble_start, ensemble_end+1)]),
                      Family("run_family").add(
                          [Task(f"run_member_{i:02d}",
                                Trigger("/" + workflow_date + f"/init_noah_family/init_noah_member_{i:02d} == complete")
)
                           for i in range(ensemble_start, ensemble_end+1)
                for i in range(ensemble_start, ensemble_end+1)]),
                      Family("postprocess_family").add(
                          [Task(f"postprocess_member_{i:02d}",
                                Trigger("/" + workflow_date + f"/run_family/run_member_{i:02d} == complete"))
                           for i in range(ensemble_start, ensemble_end+1)]),
                      ))

    print(defs)

    print("Checking job creation: .ecf -> .job0")
    print(defs.check_job_creation())

    print(f"Saving definition to file '{workflow_date}.def'")
    defs.save_as_defs(f"{workflow_date}.def")

# To restore the definition from file 'test.def' we can use:
# restored_defs = ecflow.Defs("test.def")

if __name__ == "__main__":
    workflow(__doc__)
