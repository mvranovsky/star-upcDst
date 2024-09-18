#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os
import sys

#_____________________________________________________________________________
if __name__ == "__main__":

    runList = "/star/u/truhlar/star-upcDst/work/lists/goodRuns1109.list"  
    outDir = "/star/u/truhlar/star-upcDst/work/OffSetsRun17.list"  
    out = open(outDir, 'w')
    cptRuns = []
    for line in open(runList).readlines(  ):
        run = line.lstrip().split(".")[0].split("/")[6]
        cptRuns.append( run )

    with open("pp2ppOffsets_v2017.0.2.txt", 'r+') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 5):
            line = [j for j in lines[i].lstrip().split(" ") if j] #remove empty strings
            runnumber = line[4]
            if runnumber not in cptRuns:
                continue
            line = [j for j in lines[i+1].lstrip().split(" ") if j] #remove empty strings
            E1U_X = (float(line[2]) + float(line[6]))/2000 #convert from mm to m
            E1U_Y = (float(line[0]) + float(line[4]))/2000
            E1D_X = (float(line[10]) + float(line[14]))/2000
            E1D_Y = (float(line[8]) + float(line[12]))/2000
            line = [j for j in lines[i+2].lstrip().split(" ") if j] #remove empty strings
            E2U_X = (float(line[2]) + float(line[6]))/2000
            E2U_Y = (float(line[0]) + float(line[4]))/2000
            E2D_X = (float(line[10]) + float(line[14]))/2000
            E2D_Y = (float(line[8]) + float(line[12]))/2000
            line = [j for j in lines[i+3].lstrip().split(" ") if j] #remove empty strings
            W1U_X = (float(line[2]) + float(line[6]))/2000
            W1U_Y = (float(line[0]) + float(line[4]))/2000
            W1D_X = (float(line[10]) + float(line[14]))/2000
            W1D_Y = (float(line[8]) + float(line[12]))/2000
            line = [j for j in lines[i+4].lstrip().split(" ") if j] #remove empty strings
            W2U_X = (float(line[2]) + float(line[6]))/2000
            W2U_Y = (float(line[0]) + float(line[4]))/2000
            W2D_X = (float(line[10]) + float(line[14]))/2000
            W2D_Y = (float(line[8]) + float(line[12]))/2000
            out.write(runnumber + " {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(E1U_X, E1U_Y, E1D_X, E1D_Y, E2U_X, E2U_Y, E2D_X, E2D_Y, W1U_X, W1U_Y, W1D_X, W1D_Y, W2U_X, W2U_Y, W2D_X, W2D_Y) + '\n')
