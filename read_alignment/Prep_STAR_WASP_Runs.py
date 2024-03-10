#!/home/asiimwe/anaconda3/bin/python
#Etracting all STAR+WASP-run commands per sample and passing lines to Run_STAR_WASP.sh for sequential execution (suprocess.call() runs jobs in parallel which is not desirable for our runs)
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


samples_path = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/" 

STAR_WASPdir = "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/"

os.remove("/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/Run_STAR_WASP.sh") #this file needs to be in the directory prior to running this python file
subprocess.call(["touch", "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/Run_STAR_WASP.sh"])

for path, dirs, files in os.walk(samples_path):
	#for dir in dirs:
	 #      print(dir)
	for file in files:
		if file.endswith(".sh"):
			#print(path) #/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/NA12878_PolyA/8threads
			os.chdir(path) #setting working dir to specific thread dir per sample
			bashfile = os.path.join(path, file)
			#print(bashfile) #/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/HG00732/16threads/BaseCode_STAR_WASP_Runs.sh
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
			batchrunfile = STAR_WASPdir + "Run_STAR_WASP.sh"
			st = os.stat(batchrunfile)
			os.chmod(batchrunfile, st.st_mode | stat.S_IEXEC) 
			with open(batchrunfile, "a") as file:
				file.write(line2)
				
	
