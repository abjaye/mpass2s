#!/usr/bin/env python3
import os, sys
import glob, shutil
from datetime import timedelta, datetime
import argparse
import subprocess
import warnings
warnings.filterwarnings("ignore")

def parse_command_line(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--date",
                        help="Specify a start Date")
    parser.add_argument("--ensemble-start",default=0,
                        help="Specify the first ensemble member")
    parser.add_argument("--ensemble-end",default=0,
                        help="Specify the last ensemble member")

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

    return date.strftime("%Y-%m-%d"),int(args.ensemble_start),int(args.ensemble_end)

if __name__ == '__main__':
    date, ensemble_start, ensemble_end = parse_command_line(sys.argv)

    fcsthome = os.getenv("FCST_HOME")
    workflow = os.getenv("MPAS_WORKFLOW")
    baseroot = os.getenv("FCST_WORK")
    caseroot = os.path.join(baseroot,f"{workflow}_{date}.00")
    outroot = f"/glade/derecho/scratch/espstoch/{workflow}/"
    outcaseroot = os.path.join(outroot,f"{workflow}_{date}.00")

    if workflow == "mpass2s_120km":
        cells = "40962"
        dt = "720.0"
        len_disp = "120000.0"
    elif workflow == "mpass2s":
        cells = "163842"
        dt = "420.0"
        len_disp = "60000.0"

    basemonth = int(date[5:7])
    baseyear = int(date[0:4])
    baseday = int(date[8:9])

    base_date_str = date  # starting date
    num_days = 45  # number of days

    startday = datetime.strptime(base_date_str,"%Y-%m-%d").date()
    endday = startday + timedelta(days=num_days)
    startday_str = startday.strftime("%Y-%m-%d")
    endday_str = endday.strftime("%Y-%m-%d")

    startval = "00"
    nint = len(startval)

    for i in range(ensemble_start, ensemble_end+1):
        member_string = "{{0:0{0:d}d}}".format(nint).format(i)
        caseroot = caseroot[:-nint] + member_string
        outcaseroot = outcaseroot[:-nint] + member_string
        print(caseroot)

        os.system(f"rm -rf {caseroot}")
        os.system(f"rm -rf {outcaseroot}")

        if not os.path.isdir(baseroot):
            os.mkdir(baseroot)
        if not os.path.isdir(caseroot):
            os.mkdir(caseroot)
        if not os.path.isdir(outroot):
            os.mkdir(outroot)
        if not os.path.isdir(outcaseroot):
            os.mkdir(outcaseroot)

        if i==0:
            print("using ERA5 reanalysis!")
            subprocess.run(["python",os.path.join(fcsthome,"bin","era5_to_int","era5_to_int.py"), \
                            "--po",outcaseroot,"-i",f"{date}_00"])
        else:
            print("using ERA5 ensemble "+str(i-1))
            subprocess.run(["python",os.path.join(fcsthome,"bin","era5_to_int","era5_to_int_ens.py"), \
                            "--path","/glade/campaign/mmm/c3we/jaye/era5_ens/","--po",outcaseroot,"-i", \
                            "-e",str(i-1),f"{date}_00"])


        #now need to link the SST update, static and grid files to each ens directory

        sst_update = f"/glade/work/espstoch/mpas_grid_input/x1.{cells}.sfc_update.nc"
        os.system(f"ln -sf {sst_update} {caseroot}/.")
        grid_file = f"/glade/work/espstoch/mpas_grid_input/x1.{cells}.grid.nc"
        os.system(f"ln -sf {grid_file} {caseroot}/.")
        static_file = f"/glade/work/espstoch/mpas_grid_input/x1.{cells}.static.nc"
        os.system(f"ln -sf {static_file} {caseroot}/.")
        graph_file = f"/glade/work/espstoch/mpas_grid_input/x1.{cells}.graph.info.*"
        os.system(f"cp {graph_file} {caseroot}/.")

        # copy model run dir to ens dirs

        model_dir = "/glade/work/espstoch/mpas_workflow/mpas_base/"
        os.system(f"cp {model_dir}* {caseroot}/.")

        # change dates in namelists

        with open(f"{caseroot}/namelist.init_atmosphere_temp") as fin, open(f"{caseroot}/namelist.init_atmosphere","w") as fout:
            input_lines = fin.readlines()
            for line in input_lines:
                if "    config_start_time = '" in line:
                    fout.write("    config_start_time = '{}_00:00:00'\n".format(startday_str))
                elif "    config_met_prefix = 'ERA5'" in line:
                    fout.write("    config_met_prefix = '{}/ERA5'\n".format(outcaseroot))
                elif "    config_block_decomp_file_prefix = 'x1.163842.graph.info.part.'" in line:
                    fout.write("    config_block_decomp_file_prefix = 'x1.{}.graph.info.part.'\n".format(cells))
                else:
                    fout.write(line)

        with open(f"{caseroot}/namelist.atmosphere_temp") as fin, open(f"{caseroot}/namelist.atmosphere","w") as fout:
            input_lines = fin.readlines()
            for line in input_lines:
                if "    config_dt = 360.0" in line:
                    fout.write("    config_dt = {}\n".format(dt))
                elif "    config_start_time = '" in line:
                    fout.write("    config_start_time = '{}_00:00:00'\n".format(startday_str))
                elif "    config_len_disp = 60000.0" in line:
                    fout.write("    config_len_disp = {}\n".format(len_disp))
                elif "    config_block_decomp_file_prefix = 'x1.163842.graph.info.part.'" in line:
                    fout.write("    config_block_decomp_file_prefix = 'x1.{}.graph.info.part.'\n".format(cells))
                else:
                    fout.write(line)

        # Change out directories in streams

        with open(f"{caseroot}/streams.init_atmosphere_temp") as fin, open(f"{caseroot}/streams.init_atmosphere","w") as fout:
            input_lines = fin.readlines()
            for line in input_lines:
                if "                  filename_template='x1.163842.static.nc'" in line:
                    fout.write("                  filename_template='x1.{}.static.nc'\n".format(cells))
                elif "                  filename_template='x1.163842.init.nc'" in line:
                    fout.write("                  filename_template='{}/x1.{}.init.nc'\n".format(outcaseroot,cells))
                elif "                  filename_template='x1.163842.sfc_update.nc'" in line:
                    fout.write("                  filename_template='x1.{}.sfc_update.nc'\n".format(cells))
                else:
                    fout.write(line)

        with open(f"{caseroot}/streams.atmosphere_temp") as fin, open(f"{caseroot}/streams.atmosphere","w") as fout:
            input_lines = fin.readlines()
            for line in input_lines:
                if "                  filename_template ='x1.163842.init.nc_spinup'" in line:
                    fout.write("                  filename_template ='{}/x1.{}.init.nc_spinup'\n".format(outcaseroot,cells))
                elif "                  filename_template ='restart.$Y-$M-$D_$h.$m.$s.nc'" in line:
                    fout.write("                  filename_template ='{}/restart.$Y-$M-$D_$h.$m.$s.nc'\n".format(outcaseroot))
                elif "        filename_template ='history.$Y-$M-$D_$h.$m.$s.nc'" in line:
                    fout.write("        filename_template ='{}/history.$Y-$M-$D_$h.$m.$s.nc'\n".format(outcaseroot))
                elif "        filename_template ='diag_p1p2.$Y-$M-$D_$h.$m.$s.nc'" in line:
                    fout.write("        filename_template ='{}/diag_p1p2.$Y-$M-$D_$h.$m.$s.nc'\n".format(outcaseroot))
                elif "        filename_template ='diag_sfc.$Y-$M-$D_$h.$m.$s.nc'" in line:
                    fout.write("        filename_template ='{}/diag_sfc.$Y-$M-$D_$h.$m.$s.nc'\n".format(outcaseroot))
                elif "        filename_template='x1.163842.sfc_update.nc'" in line:
                    fout.write("        filename_template='x1.{}.sfc_update.nc'\n".format(cells))
                else:
                    fout.write(line)
