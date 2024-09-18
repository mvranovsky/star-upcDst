#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os
import sys
import re
import shutil

#_____________________________________________________________________________
if __name__ == "__main__":

    #get the name of root file to process
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./runPlotManager.py
            
    #currDir = os.getcwd()
    
    basedir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/"
    if len(args) == 1:
        tag = str(args.pop(0))
        if tag == "l":
            infile = "AnalysisOutput.root"
            location = "/star/u/mvranovsk/star-upcDst/work/"
        else:
            infile = "StRP_production_0000.root"
            location = basedir + tag + "/merged/"
    else: 
        #get outputs directory from submit used for production
        for line in open("submit.xml").xreadlines(  ): 
            if "Location" in line:
                basedir = line.lstrip().split(">")[1].split("/sched")[0]
                tag = basedir.rsplit("/",1)[1]
        infile = "StRP_production_0000.root"
        location = basedir + "/merged/"

    #Copy libstar-upc.so
    shutil.copyfile("/star/u/mvranovsk/star-upcDst/build/libstar-upc.so", "libstar-upc.so")

    print "  Input file:     ", infile
    print "  File location:     ", location
    #print "  Tag :      ", tag
    if tag != "l":
        shutil.copyfile("/star/u/mvranovsk/star-upcDst/work/include/RunDef.h", "/star/u/mvranovsk/star-upcDst/work/include/RunDef_backup.h")
        shutil.copyfile("/star/u/mvranovsk/star-upcDst/work/config/" + tag + ".h", "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h")

    make = Popen("make", stdout=PIPE, stderr=PIPE).communicate()
    if len(make[1]) != 0:
        print make[1]
    #exit()



    defaultDir = "/gpfs01/star/pwg/mvranovsk/Embedding/"            
    for line in open("anaEmbed/submit.xml").xreadlines(  ): 
        if "Location" in line:
            embedDir = line.lstrip().split(">")[1].split("/sched")[0]

    embedFile = embedDir + "/merged/StRP_production_0000.root"
    root_cmd = "./PlotManager " + location + infile #"\",\"" + embedFile + "\")"
    root = Popen(root_cmd.split(), stdout=PIPE, stderr=PIPE).communicate()

    print root[0], root[1]
    if tag != "l":
        shutil.copyfile("/star/u/mvranovsk/star-upcDst/work/include/RunDef_backup.h", "/star/u/mvranovsk/star-upcDst/work/include/RunDef.h")
    print "All done."
