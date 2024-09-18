#!/usr/bin/python

#--------------------------------------------------------------
# macro to submit rpDst analysis 
#
# usage:
# ./SubmitPlugin.py outputDir inputSource
# ./SubmitPlugin.py outputDir
#
#--------------------------------------------------------------

import os
from subprocess import Popen
import re
import sys
import time


#_____________________________________________________________________________
if __name__ == "__main__":

    #name of config file from command line argument
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./SubmitPlugin
    if len(args) == 2:
        outputDir = args.pop(0) # read third input from terminal = outputDir
        inputSource = args.pop(0) # read second input from terminal = inputSource
    elif len(args) == 1:
        outputDir = args.pop(0) # read third input from terminal = outputDir
        #inputSource = "/star/u/mvranovsk/star-upcDst/build/goodRun17.list"
        inputSource = "/star/u/mvranovsk/star-upcDst/build/inputcheck.list"

    else:
        print "Wrong input arguments"
        exit()

    exe = "Ana"

    top = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic"
    srcDir = os.getcwd()
    basedir = top + "/" + outputDir
    inputSource = "<input URL=\"filelist:" + inputSource + "\"/>"
    print("Input source: "+inputSource)
    print("Output dir: "+basedir)
    print("Source dir: "+srcDir)
    #create jobs output folders
    Popen("mkdir -p "+basedir, shell=True).communicate()
    #proceed with scheduler submit
    Popen("mkdir "+basedir+"/logs", shell=True).communicate()
    Popen("mkdir "+basedir+"/sched", shell=True).communicate()
    #create xml
    xmlname = "submit.xml"
    xml = open(xmlname, "w")
    for line in open("scheduler_template.xml", "r"):
        line = re.sub(r"__BASEDIR__", basedir, line)
        line = re.sub(r"__INPUT__", inputSource, line)
        line = re.sub(r"__SRCDIR__", srcDir, line)
        line = re.sub(r"__EXE__", exe, line)
        xml.write(line)
    xml.close()
    #submit for a given catalog query
    Popen("star-submit "+xmlname, shell=True).communicate()
    print "All done."




















