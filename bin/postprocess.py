#!/usr/bin/env python3
import os, sys
import glob, shutil
from pathlib import Path
from datetime import timedelta, datetime
import argparse
import subprocess
import xarray as xr
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

def parse_command_line(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--date",
                        help="Specify a start Date")
    parser.add_argument("--ens",default=0,
                        help="Specify the ensemble member")

    args = parser.parse_args()

    if args.date:
        try:
            date = datetime.strptime(args.date, '%Y-%m-%d')
        except ValueError as verr:
            raise ValueError("Incorrect data format, should be YYYY-MM-DD or YYYY-MM") from verr
    elif cdate:
        date = datetime.strptime(cdate, '%Y-%m-%d')
    else:
        date = datetime.today() - timedelta(days=1)

    return date.strftime("%Y-%m-%d"),args.ens

if __name__ == '__main__':
    date, ens = parse_command_line(sys.argv)

    mem = f"{int(ens):02d}"

    fcsthome = os.getenv("FCST_HOME")
    workflow = os.getenv("MPAS_WORKFLOW")
    baseroot = os.getenv("FCST_WORK")
    caseroot = os.path.join(baseroot, f"{workflow}_{date}.{mem}")
    outroot = f"/glade/derecho/scratch/espstoch/{workflow}/"
    outcaseroot = os.path.join(outroot, f"{workflow}_{date}.{mem}")
    #path_tmp = "/glade/derecho/scratch/espstoch/tmp/"

    if workflow == "mpass2s_120km":
        cells = "40962"
    elif workflow == "mpass2s":
        cells = "163842"

    month_abbr = ["","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]
    yyyy = date.split("-")[0]; mm = date.split("-")[1]; dd = date.split("-")[2]
    outdate = dd + month_abbr[int(mm)] + yyyy

    files = f"{outcaseroot}/diag_sfc*.nc"
    files2 = f"{outcaseroot}/diag_p1p2*.nc"
    grid_path = f"{caseroot}/x1.{cells}.grid.nc"
    uxds = xr.open_mfdataset(files,combine="nested",concat_dim="Time").astype("float32")
    uxds2 = xr.open_mfdataset(files2,combine="nested",concat_dim="Time").astype("float32")
    uxds["Time"] = pd.to_datetime(uxds["Time"].values).round("3H")
    uxds2["Time"] = pd.to_datetime(uxds2["Time"].values).round("6H")

    comp = dict(zlib=True, complevel=4, dtype="float32")

    # Make p1 calculations
    prec_c = uxds["prec_acc_c"].resample(Time="D",closed="right",label="right").sum()
    prec_nc = uxds["prec_acc_nc"].resample(Time="D",closed="right",label="right").sum()
    pr_sfc = (prec_c + prec_nc)/86400
    pr_sfc.attrs ={"units": "kg m^-2 s^-1", "long_name": "Total (convective and large-scale) precipitation rate (liq + ice)"}
    pr_sfc = pr_sfc.to_dataset(name="pr_sfc")

    rlut = uxds2["olrtoa"].resample(Time="D").mean().to_dataset(name="rlut")
    tas_2m = uxds["t2m"].resample(Time="D").mean().to_dataset(name="tas_2m")
    ts = uxds["skintemp"].resample(Time="D").mean().to_dataset(name="ts")

    zg_200 = uxds2["height_200hPa"].resample(Time="D").mean().to_dataset(name="zg_200")
    zg_500 = uxds2["height_500hPa"].resample(Time="D").mean().to_dataset(name="zg_500")

    ua_200 = uxds2["uzonal_200hPa"].resample(Time="D").mean().to_dataset(name="ua_200")
    ua_850 = uxds2["uzonal_850hPa"].resample(Time="D").mean().to_dataset(name="ua_500")
    va_200 = uxds2["umeridional_200hPa"].resample(Time="D").mean().to_dataset(name="ua_200")
    va_850 = uxds2["umeridional_850hPa"].resample(Time="D").mean().to_dataset(name="ua_850")

    datasetsp1 = [rlut,pr_sfc,tas_2m,ts,zg_200,zg_500,ua_200,ua_850,va_200,va_850]
    varsp1 = ["rlut","pr_sfc","tas_2m","ts","zg_200","zg_500","ua_200","ua_850","va_200","va_850"]

    for name, ds in zip(varsp1, datasetsp1):
        print(name)
        path1 = f"{outroot}/{workflow}_native/p1/{name}/{yyyy}/{mm}/"
        if not os.path.exists(path1):
            os.system(f"mkdir -p {path1}")
        path2 = f"{outroot}/{workflow}_latlon/p1/{name}/{yyyy}/{mm}/"
        if not os.path.exists(path2):
            os.system(f"mkdir -p {path2}")
        filename = f"{path1}{name}_{workflow}_{outdate}_00z_unst_d46_m{mem}.nc"
        filename_regrid = f"{path2}{name}_{workflow}_{outdate}_00z_d01_d46_m{mem}.nc"
        encoding = {var: comp for var in ds.data_vars}
        ds = ds.chunk({"Time": 1, "nCells": int(cells)})
        ds.to_netcdf(filename,unlimited_dims="Time",encoding=encoding,mode="w")
        os.system(f"CDO_COMPRESS=4 cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:{grid_path} {filename} {filename_regrid}")

    # Make p2 calculations

    tc850 = uxds2["temperature_850hPa"].resample(Time="D").mean() - 273.15
    rh850 = uxds2["relhum_850hPa"].resample(Time="D").mean()
    pressure = 850.
    e_s = 6.112 * np.exp((17.67 * tc850) / (tc850 + 243.5))
    e = rh850 / 100 * e_s
    q850 = (0.622 * e) / (pressure - (0.378 * e))
    q850.attrs ={"units": "kg/kg", "long_name": "Specific Humidity at 850 mbar pressure surface"}
    huss_850 = q850.to_dataset(name="huss_850")

    q2m = uxds["q2"].resample(Time="D").mean()
    sfcp = uxds["surface_pressure"].resample(Time="D").mean()/100
    tdc = q2m*sfcp/(.622 + q2m)
    tdc = tdc.where(tdc>0.001, 0.001)
    tdps = ((243.5*np.log(tdc) - 440.8)/(19.48 - np.log(tdc))) + 273.15
    tdps.attrs ={"units": "K", "long_name": "Dewpoint temperature at 2 meters"}
    tdps = tdps.to_dataset(name="tdps")

    swdnb = uxds["swdnb"].resample(Time="D").mean()
    swupb = uxds["swupb"].resample(Time="D").mean()
    lwdnb = uxds["lwdnb"].resample(Time="D").mean()
    lwupb = uxds["lwupb"].resample(Time="D").mean()
    rad_sfc = lwupb-lwdnb+swupb-swdnb
    rad_sfc.attrs ={"units": "W/m2", "long_name": "Net surface radiation (lwupb-lwdnb+swupb-swdnb)"}
    rad_sfc = rad_sfc.to_dataset(name="rad_sfc")

    smois = uxds["smois"].resample(Time="D").mean()
    mrso = smois.sum(dim="nSoilLevels")
    mrso.attrs ={"units": "kg/m2", "long_name": "Vertically integrated soil moisture"}
    mrso = mrso.to_dataset(name="mrso")

    snow = uxds2["snow"].resample(Time="D").mean()
    snowh = uxds2["snowh"].resample(Time="D").mean()
    snd = snow/snowh
    snd.attrs ={"units": "kg/m3", "long_name": "Snow Density (snow/snowh)"}
    snd = snd.to_dataset(name="snd")

    snc = uxds2["snowc"].resample(Time="D").mean()*100
    snc.attrs ={"units": "%", "long_name": "Snow Cover"}
    snc = snc.to_dataset(name="snc")
    sic = uxds2["xice"].resample(Time="D").mean()*100
    sic.attrs ={"units": "%", "long_name": "Sea Ice Concentration"}
    sic = sic.to_dataset(name="sic")

    w500 = uxds2["w_500hPa"].resample(Time="D").mean()
    tk500 = uxds2["temperature_500hPa"].resample(Time="D").mean()
    wap_500 = -((50000.*9.81)/(287.04*tk500))*w500
    wap_500.attrs ={"units": "Pa/s", "long_name": "Vertical Velocity (omega)"}
    wap_500 = wap_500.to_dataset(name="wap_500")

    tasmax_2m = uxds2["tasmax"].resample(Time="D",closed="right",label="right").max().to_dataset(name="tasmax_2m")
    tasmin_2m = uxds2["tasmin"].resample(Time="D",closed="right",label="right").min().to_dataset(name="tasmin_2m")
    hfss_sfc = uxds["hfx"].resample(Time="D").mean().to_dataset(name="hfss_sfc")
    hfls_sfc = uxds["lh"].resample(Time="D").mean().to_dataset(name="hfls_sfc")
    ua_100 = uxds2["uzonal_100hPa"].resample(Time="D").mean().to_dataset(name="ua_100")
    va_100 = uxds2["umeridional_100hPa"].resample(Time="D").mean().to_dataset(name="va_100")
    uas = uxds["u10"].resample(Time="D").mean().to_dataset(name="uas")
    vas = uxds["v10"].resample(Time="D").mean().to_dataset(name="vas")
    psl = uxds["mslp"].resample(Time="D").mean().to_dataset(name="psl")
    swe = uxds["snow_acc_nc"].resample(Time="D",closed="right",label="right").sum().to_dataset(name="swe")
    cape = uxds["cape"].resample(Time="D").mean().to_dataset(name="cape")

    datasetsp2 = [huss_850,tasmax_2m,tasmin_2m,hfss_sfc,hfls_sfc,wap_500,ua_100,va_100,uas,vas,tdps,psl,swe,rad_sfc,snd,snc,mrso,sic,cape]
    varsp2 = ["huss_850","tasmax_2m","tasmin_2m","hfss_sfc","hfls_sfc","wap_500","ua_100","va_100","uas","vas","tdps","psl","swe","rad_sfc","snd","snc","mrso","sic","cape"]


    for name, ds in zip(varsp2, datasetsp2):
        print(name)
        path1 = f"{outroot}/{workflow}_native/p2/{name}/{yyyy}/{mm}/"
        if not os.path.exists(path1):
            os.system(f"mkdir -p {path1}")
        path2 = f"{outroot}/{workflow}_latlon/p2/{name}/{yyyy}/{mm}/"
        if not os.path.exists(path2):
            os.system(f"mkdir -p {path2}")
        filename = f"{path1}{name}_{workflow}_{outdate}_00z_unst_d46_m{mem}.nc"
        filename_regrid = f"{path2}{name}_{workflow}_{outdate}_00z_d01_d46_m{mem}.nc"
        encoding = {var: comp for var in ds.data_vars}
        ds = ds.chunk({"Time": 1, "nCells": int(cells)})
        ds.to_netcdf(filename,unlimited_dims="Time",encoding=encoding,mode="w")
        os.system(f"CDO_COMPRESS=4 cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:{grid_path} {filename} {filename_regrid}")


    # Now for daily plevs for OMEGA, RELHUM, T, U, V, Z
    levs = [50,100,200,250,500,700,850,925]
    omega_lev = {}; relhum_lev = {}; t_lev = {}; u_lev = {}; v_lev = {}; z_lev = {}
    omega = {}; relhum = {}; t = {}; u = {}; v = {}; z = {}
    for lev in levs:
        wlev = f"w_{str(lev)}hPa"; tlev = f"temperature_{str(lev)}hPa";
        ulev = f"uzonal_{str(lev)}hPa"; vlev = f"umeridional_{str(lev)}hPa";
        rhlev = f"relhum_{str(lev)}hPa"; zlev = f"height_{str(lev)}hPa";

        w = uxds2[wlev].resample(Time="D").mean()
        t_lev[lev] = uxds2[tlev].resample(Time="D").mean()
        omega_lev[lev] = -(((lev*100.)*9.81)/(287.04*t_lev[lev]))*w
        u_lev[lev] = uxds2[ulev].resample(Time="D").mean()
        v_lev[lev] = uxds2[vlev].resample(Time="D").mean()
        relhum_lev[lev] = uxds2[rhlev].resample(Time="D").mean()
        z_lev[lev] = uxds2[zlev].resample(Time="D").mean()

    t = xr.concat([t_lev[l] for l in levs],dim="lev_p").assign_coords(lev_p=levs)
    t.attrs ={"units": "K", "long_name": "Temperature", "cell_methods": "time: mean"}
    t = t.to_dataset(name="T"); t["T"] = t["T"].transpose("Time","lev_p","nCells")
    omega = xr.concat([omega_lev[l] for l in levs],dim="lev_p").assign_coords(lev_p=levs)
    omega.attrs ={"units": "Pa/s", "long_name": "Vertical Velocity (pressure)", "cell_methods": "time: mean"}
    omega = omega.to_dataset(name="OMEGA"); omega["OMEGA"] = omega["OMEGA"].transpose("Time","lev_p","nCells")
    u = xr.concat([u_lev[l] for l in levs],dim="lev_p").assign_coords(lev_p=levs)
    u.attrs ={"units": "m/s", "long_name": "Zonal wind", "cell_methods": "time: mean"}
    u = u.to_dataset(name="U"); u["U"] = u["U"].transpose("Time","lev_p","nCells")
    v = xr.concat([v_lev[l] for l in levs],dim="lev_p").assign_coords(lev_p=levs)
    v.attrs ={"units": "m/s", "long_name": "Meridional wind", "cell_methods": "time: mean"}
    v = v.to_dataset(name="V"); v["V"] = v["V"].transpose("Time","lev_p","nCells")
    z = xr.concat([z_lev[l] for l in levs],dim="lev_p").assign_coords(lev_p=levs)
    z.attrs ={"units": "m", "long_name": "Geopotential height", "cell_methods": "time: mean"}
    z = z.to_dataset(name="Z"); z["Z"] = z["Z"].transpose("Time","lev_p","nCells")
    relhum = xr.concat([relhum_lev[l] for l in levs],dim="lev_p").assign_coords(lev_p=levs)
    relhum.attrs ={"units": "percent", "long_name": "Relative humidity", "cell_methods": "time: mean"}
    relhum = relhum.to_dataset(name="RELHUM"); relhum["RELHUM"] = relhum["RELHUM"].transpose("Time","lev_p","nCells")

    datasetslev = [t,omega,u,v,z,relhum]
    varlev = ["T","OMEGA","U","V","Z","RELHUM"]
    for name, ds in zip(varlev, datasetslev):
        print(name)
        path1 = f"{outroot}/{workflow}_native/{name}/{yyyy}/{mm}/"
        if not os.path.exists(path1):
            os.system(f"mkdir -p {path1}")
        #path2 = f"{outroot}/{workflow}_latlon/{name}/{yyyy}/{mm}/"
        #if not os.path.exists(path2):
        #    os.system(f"mkdir -p {path2}")
        filename = f"{path1}{name}_{workflow}_{outdate}_00z_unst_d46_m{mem}.nc"
        #filename_regrid = f"{path2}{name}_{workflow}_{outdate}_00z_d01_d46_m{mem}.nc"
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(filename,unlimited_dims="Time",encoding=encoding,mode="w")
        #for i, lev in enumerate(ds.lev_p.values):
        #    sub = ds.isel(lev_p=i)
        #    filename_sub = f"{path_tmp}{name}_{int(lev)}hPa_{workflow}_{outdate}_m{mem}.nc"
        #    filename_sub_regrid = f"{path_tmp}{name}_{int(lev)}hPa_{workflow}_{outdate}_m{mem}_regrid.nc"
        #    sub.to_netcdf(filename_sub,unlimited_dims="Time",encoding={name: comp},mode="w")
        #    os.system(f"CDO_COMPRESS=4 cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:{grid_path} {filename_sub} {filename_sub_regrid}")
        #merged_ds = xr.concat([xr.open_dataset(f"{path_tmp}{name}_{lev}hPa_{workflow}_{outdate}_m{mem}_regrid.nc")[name].assign_coords(lev_p=lev) for lev in levs],dim="lev_p").to_dataset(name=name).transpose("Time","lev_p","lat","lon")
        #merged_ds.to_netcdf(filename_regrid,unlimited_dims="Time",encoding={name: comp},mode="w")

    # Write out, compress, and regrid 3hourly data
    path1 = f"{outroot}/{workflow}_native/3hourly/{yyyy}/{mm}/"
    if not os.path.exists(path1):
        os.system(f"mkdir -p {path1}")
    #path2 = f"{outroot}/{workflow}_latlon/3hourly/{yyyy}/{mm}/"
    #if not os.path.exists(path2):
    #    os.system(f"mkdir -p {path2}")
    filename = f"{path1}diag_sfc_{workflow}_{outdate}_00z_unst_d46_m{mem}.nc"
    #filename_2d = f"{path_tmp}diag_sfc_{workflow}_{outdate}_m{mem}_2d.nc"
    #filename_regrid_2d = f"{path_tmp}diag_sfc_{workflow}_{outdate}_m{mem}_2d_regrid.nc"
    #filename_regrid = f"{path2}diag_sfc_{workflow}_{outdate}_00z_d01_d46_m{mem}.nc"
    encoding = {var: comp for var in uxds.data_vars}
    uxds.to_netcdf(filename,unlimited_dims="Time",encoding=encoding,mode="w")
    #vars_3d = [v for v in uxds.data_vars if "nSoilLevels" in uxds[v].dims]
    #vars_2d = [v for v in uxds.data_vars if "nSoilLevels" not in uxds[v].dims]
    #uxds_2d = uxds[vars_2d]
    #uxds_2d.to_netcdf(filename_2d,unlimited_dims="Time",encoding={v: comp for v in vars_2d},mode="w")
    #os.system(f"CDO_COMPRESS=4 cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:{grid_path} {filename_2d} {filename_regrid_2d}")
    #for var in vars_3d:
    #    print(f"Regridding 3D var: {var}")
    #    lev_values = uxds[var].coords["nSoilLevels"].values if "nSoilLevels" in uxds[var].coords else range(uxds[var].shape[-1])
    #    files = []
    #    for i, lev in enumerate(lev_values):
    #        sub = uxds[var].isel(nSoilLevels=i).to_dataset(name=var)
    #        filename_sub = f"{path_tmp}{var}_{int(lev)}_{workflow}_{outdate}_m{mem}.nc"
    #        filename_sub_regrid = f"{path_tmp}{var}_{int(lev)}_{workflow}_{outdate}_m{mem}_regrid.nc"
    #        sub.to_netcdf(filename_sub, unlimited_dims="Time", encoding={var: comp},mode="w")
    #        os.system(f"CDO_COMPRESS=4 cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:{grid_path} {filename_sub} {filename_sub_regrid}")
    #        ds_regridded = xr.open_dataset(filename_sub_regrid)[var].assign_coords(nSoilLevels=lev)
    #        files.append(ds_regridded)
    #    ds_out = xr.concat(files, dim="nSoilLevels").to_dataset(name=var)
    #    filename_regrid_3d_combine = f"{path_tmp}{var}_{workflow}_{outdate}_m{mem}_regrid.nc"
    #    ds_out.transpose("Time","nSoilLevels","lat","lon").to_netcdf(filename_regrid_3d_combine,unlimited_dims="Time",encoding={var: comp},mode="w")
    #ds_3d_regridded = xr.merge([xr.open_dataset(f"{path_tmp}{var}_{workflow}_{outdate}_m{mem}_regrid.nc") for var in vars_3d])
    #ds_2d_regridded = xr.open_dataset(filename_regrid_2d)
    #ds_full = xr.merge([ds_2d_regridded,ds_3d_regridded])
    #ds_full.to_netcdf(filename_regrid, unlimited_dims="Time",encoding={var: comp for var in ds_full.data_vars},mode="w")


    # Write out, compress, and regrid 6hourly data
    path1 = f"{outroot}/{workflow}_native/6hourly/{yyyy}/{mm}/"
    if not os.path.exists(path1):
        os.system(f"mkdir -p {path1}")
    #path2 = f"{outroot}/{workflow}_latlon/6hourly/{yyyy}/{mm}/"
    #if not os.path.exists(path2):
    #    os.system(f"mkdir -p {path2}")
    filename = f"{path1}diag_p1p2_{workflow}_{outdate}_00z_unst_d46_m{mem}.nc"
    #filename_regrid = f"{path2}diag_p1p2_{workflow}_{outdate}_00z_d01_d46_m{mem}.nc"
    encoding = {var: comp for var in uxds2.data_vars}
    uxds2.to_netcdf(filename,unlimited_dims="Time",encoding=encoding,mode="w")
    #os.system(f"CDO_COMPRESS=4 cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:{grid_path} {filename} {filename_regrid}")

    # Delete temp files
    #temp_file = glob.glob(f"{path_tmp}*{workflow}_{outdate}_m{mem}*")
    #for f in temp_file:
    #    path = Path(f)
    #    if path.exists():
    #        path.unlink()
