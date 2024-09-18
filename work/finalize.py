#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os
import sys

#_____________________________________________________________________________
if __name__ == "__main__":

    #get the name of root file to process
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./finalize.py
            
    basedir = "/gpfs01/star/pwg/truhlar/Run17_P20ic/"
    if len(args) == 1:
        tag = str(args.pop(0))
        if tag == "l":
            infile = "AnalysisOutput.root"
            location = "/star/u/truhlar/star-upcDst/work/"
        else:
            infile = "StRP_production_0000.root"
            location = basedir + tag + "/merged/"
    else: 
        #get outputs directory from submit used for production
        for line in open("submit.xml").xreadlines(  ): 
            if "Location" in line:
                basedir = line.lstrip().split(">")[1].split("/sched")[0]
        infile = "StRP_production_0000.root"
        location = basedir + "/merged/"


    print "  Input file:     ", infile
    print "  File location:     ", location
    
    #root_cmd = "root -l -b -q /star/u/truhlar/star-upcDst/work/finalize.C+(\"" + location + infile + "\")"
    root_cmd = "root -l -b -q /star/u/truhlar/star-upcDst/work/plotManager.C+(\"" + location + infile + "\")"
    root = Popen(root_cmd.split(), stdout=PIPE, stderr=PIPE).communicate()

    print root[0], root[1]
    print "All done."
