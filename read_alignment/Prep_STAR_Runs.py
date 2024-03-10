#!/home/asiimwe/anaconda3/bin/python
#Etracting all STAR-run commands per sample and passing lines to Run_STAR.sh for sequential execution 
import sys
import csv
import os
import string
import subprocess
from itertools import chain
import pandas as pd
from os import mkdir
import shutil
import glob
import fileinput
import stat


samples_path = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs/" 

STARdir = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/"

os.remove("/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/Run_STAR.sh") #this file needs to be in the directory prior to running this python file
subprocess.call(["touch", "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/Run_STAR.sh"])

for path, dirs, files in os.walk(samples_path):
	for file in files:
		if file.endswith(".sh"):
			#print(path) #/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs/NA19238/16threads
			os.chdir(path) #setting working dir to specific thread dir per sample
			bashfile = os.path.join(path, file)
			#print(bashfile) #/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs/HG00514/16threads/BaseCode_STAR_Runs.sh
			time_command = "/usr/bin/time -v -o "
			threads = path.split("/")[9]
			#print(threads)
			resource_log = threads + "_resource_log.txt"
			#print(resource_log)
			resource_log_path = os.path.join(path, resource_log)
			#print(resource_log_path)
			line = time_command + resource_log_path + " " + bashfile 
			#print(line)
			line2 = time_command + resource_log_path + " " + bashfile  + "\n"
			#print(line2)
			#subprocess.call(line, executable='/bin/bash', shell = True) #subprocess.call runs jobs in parallel which we don't want. Each line to be executed will be written out to a bash file and run in sequence
			batchrunfile = STARdir + "Run_STAR.sh"
			st = os.stat(batchrunfile)
			os.chmod(batchrunfile, st.st_mode | stat.S_IEXEC) 
			with open(batchrunfile, "a") as file:
				file.write(line2)
				
	
