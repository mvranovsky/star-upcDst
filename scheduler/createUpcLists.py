#!/usr/bin/python

# Macro to create lists for new upcDst production 
# from merging the old upcDst and picoDst output
# It takes direction to picoDst output and create a list containg the output and old upcDst 

import os
from subprocess import Popen
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
    else:
        print "missing input"
        exit()

    upcDds_location = "/gpfs01/star/pwg_tasks/upc03/pp17/UPC_P20ic/"

    with open("config_from2upcDst_template.in", 'r') as source_file, open("config_from2upcDst.in", 'w') as dest_file:
        for line in source_file:
            if 'add_input' in line:
                for index_file in open(inputSource, "r"):
                    run_number = index_file[index_file.rfind("/") - 8:index_file.rfind("/")]
                    main_file = upcDds_location + run_number + ".root\n"
                    file_list = open("lists/" + run_number + ".list", "w")
                    file_list.write(main_file)
                    file_list.write(index_file)
                    input_line = "add_input " + run_number + " " + os.path.abspath(file_list.name) + "\n"
                    file_list.close()
                    dest_file.write(input_line)

            else:
                dest_file.write(line)  # Write the original line
















