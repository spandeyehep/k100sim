import os, sys
import uproot, awkward
#import ROOT
import numpy as np

def getNcapIndex(list_):
	for i,it in enumerate(list_):
		if(it == 1):
			return i
	return -1

def GetNcapEvents(nCap):
	nevents = 0
	for list__ in nCap:
		ncap_index = getNcapIndex(list__)
		if(ncap_index >= 0):
			nevents += 1
	return nevents

def generate_nrCascadeSim(events,outfilename):
	executable = "/Users/shubhampandey/work/NuclearCascade/nrCascadeSim_build/realizeCascades"
	print ("] Using executable = %s"%(executable))
	if(not os.path.exists(executable)):
		print ("] Executable not found. Please make sure file name and path is correct.")
		print ("] Exiting.")
		sys.exit(0)
	cmd = executable + " -n " + str(events) + " -o " + outfilename + " /Users/shubhampandey/work/NuclearCascade/nrCascadeSim/levelfiles/Si28_ngam_all_cascades_rfmt_sorted.txt"
	print ("\n \n Running : "+ cmd + "\n \n" )
	os.system(cmd)
	print ("] File generated: %s \n"%(outfilename))
	if(not os.path.exists(outfilename)):
		print ("something went wrong in nrCascadeSim file generation.")
		sys.exit(0)

