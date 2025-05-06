#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os
import sys


#_____________________________________________________________________________
if __name__ == "__main__":

    #merge output files by chunks of a given size
    #get outputs directory from submit used for production

    args = sys.argv
    args.pop(0) # cut first input from terminal = ./"macro".py
    deleteFiles = False
    for arg in args:
        if arg == "-del":
            deleteFiles = True


            
    defaultdir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/"

    for line in open("submit.xml").xreadlines(  ): 
        if "Location" in line:
            basedir = line.lstrip().split(">")[1].split("/sched")[0]

    if len(args) == 1:
        basedir = defaultdir + str(args.pop(0))

    outfile = "StRP_production.root"
    outlist = "StRP_production.list"

    #basedir = "/gpfs01/star/pwg/truhlar/Run17_P20ic/set3a"
    print "  Basedir:     ", basedir 
    print "  Output file:     ", outfile
    print "  Output list:     ", outlist

    #size of one output chunk
    chunksiz = int(5e6)  # approx kB

    #create output directory for merged files
    merge_path = basedir+"/merged"
    if os.access(merge_path, os.F_OK) == False:
        os.makedirs(merge_path)

    #output files list
    chunk_list = open(merge_path+"/"+outlist, "w")

    #list done jobs
    totsiz = 0    # total chunk size
    ichunk = 0    # chunk index
    nchunk = 0    # number of files in chunk
    tmpnam = "files.tmp"   # temporary to list files for a given chunk
    tmp = open(tmpnam, "w")
    cmd = "ls -s " + basedir + "/*.root"
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")

    #remove empty elements (last was empty string)
    out = [i for i in out if i]

    print "Number of all outputs:", len(out)
    print
    #list loop
    for iline in xrange(len(out)):
        fline = out[iline]
        size = 0
        if len(fline) > 0:
            siznam = fline.lstrip().split(" ")
            size = int(siznam[0])
            name = siznam[1]
        if size > 0:
            totsiz += size
            tmp.write(name+"\n")
            nchunk += 1
        if totsiz >= chunksiz or iline == len(out)-1:
            tmp.close()
            print "Chunk:", ichunk
            print "Num files:", nchunk
            print "Size:", totsiz
            iform = "_{0:04d}".format(ichunk)
            chunkout = merge_path+"/"+outfile.split(".root")[0]+iform+".root"
            merge_cmd = "root -l -b -q MergeFiles.C(\"" + tmpnam + "\",\"" + chunkout + "\")"
            merg = Popen(merge_cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            print merg[0], merg[1]
            chunk_list.write(chunkout+"\n")
            #show the merged file
            os.system("ls -alh "+chunkout)
            #delete the original files
            if deleteFiles:
                print "Deleting original files..."
                tmp = open(tmpnam, "r")
                lines = tmp.xreadlines()
                for line in lines:
                    del_cmd = "rm " + line
                    del_out = Popen(del_cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            #reset the temporary and indices
            totsiz = 0
            nchunk = 0
            tmp = open(tmpnam, "w")
            #increment the chunk index
            ichunk += 1
    #list loop end

    #clean-up the list temporary file
    tmp.close()
    #os.remove(tmpnam)

    print "All done."
