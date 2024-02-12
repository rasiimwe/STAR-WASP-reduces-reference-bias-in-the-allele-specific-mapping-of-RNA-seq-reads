#!/home/asiimwe/anaconda3/bin/python
#Extracting summary mapping statistics of all runs from respective Log.final.out files per sample and run
import sys
import csv
import os
import string
import subprocess
from itertools import chain
import pandas as pd
import re

STAR_path = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs/"
STAR_WASP_path = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/"
WASP_path = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/WASP_Runs/"

os.remove("/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/dataExtractions/final_log_results_all_runs.txt") 

subprocess.call(["touch", "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/dataExtractions/final_log_results_all_runs.txt"])

outputfile = open("/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/dataExtractions/final_log_results_all_runs.txt", "a")


wrt_header1= "Sample" + '|' + "Run" + '|' + "Thread" + '|' + "Param"  + '|' + "Value" + "\n"
outputfile.write(wrt_header1)


paths = (STAR_path, STAR_WASP_path, WASP_path)

for path, dirs, files in chain.from_iterable(os.walk(path) for path in paths):
	for file in files:
		if file.endswith("Log.final.out"):
			pathx = os.path.join(path, file)
			sample_name = "/".join(pathx.split("/")[-3:-2])
			run = "/".join(pathx.split("/")[-4:-3])
			threads = "/".join(pathx.split("/")[-2:-1])
			f1 = open(pathx, "r")
			reader = f1.readlines()
			for i in reader[1:]:
				param = ':'.join(i.split(':')[:1])
				value = ':'.join(i.split(':')[1:])
				line_out = sample_name + "|" + run + "|" + threads + "|" + param +  "%s" %value  + "\n"
				outputfile.write(line_out)


