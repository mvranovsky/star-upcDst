#!/usr/bin/python

# Macro to create lists for new upcDst production 
# from merging the old upcDst and picoDst output
# It takes direction to picoDst output and create a list containg the output and old upcDst 

import os
from subprocess import Popen, PIPE
from os.path import exists
import re
import sys
import time


#_____________________________________________________________________________
if __name__ == "__main__":

    # Setting the input list of upcDst files
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./SubmitPlugin
    if len(args) == 1:
        inputSource = args.pop(0) # read second input from terminal = inputSource
    	listName = "createdList.list"
    elif len(args) == 2:
    	listName = args.pop(0)
    	inputSource = args.pop(0)
    else:
        print "missing input. First argument should be the directory with the files to create a list of. Second argument is the name of the created list."
        exit()


    cmd = "ls " + inputSource +"*.root"
    out = Popen(cmd.split(), stdout=PIPE).communicate()[0].split("\n")

    with open(listName, "w") as dest_file:
    	for line in out:
    		file = line
    		line.replace(".root", "")
    		destination = inputSource + line + "/"+file + "\n"
    		print destination
    		dest_file.write(destination)
    	dest_file.close()


    	
