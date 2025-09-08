#!/bin/bash
# Abby Jaye, April 25, 2024

cd /glade/u/home/espstoch/mpass2s/ecflow
export ECF_PORT=40179
export ECF_HOST=derecho7
export CESM_WORKFLOW=mpass2s
export PROJECT=CESM0020

d=2018-01-01 # specify monday start date (the INCLUDED monday)
while [ "$d" != 2018-01-08 ]; do # specify monday end date (the EXCLUDED monday)
	year=$(date -d "$d" +%Y)
	month=$(date -d "$d" +%m)
	day=$(date -d "$d" +%d)
	if (( ${month#0} < 04 || ${month#0} > 10 )); then
	echo ${d}
	echo ${year}_${month}_${day}
	python workflow.py --date ${d}
	python mpass2s_${year}_${month}_${day}/client.py
	fi
	echo	
	d=$(date -I -d "$d + 7 day")
done

