#!/usr/bin/python

import random
from subprocess import Popen, PIPE
import subprocess
import os
import sys

printOut = True

filename = "BEMCParameters.txt"

gpfs_path = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/EmbeddingJPsi_BemcSysStudy/merged/"

out_path = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/EmbeddingJPsi_BemcSysStudy/results/"

def writeToFile(eps0, n, pTThr, sigma):
    stringToWrite = "%.3f %.3f %.3f %.3f" % (eps0, n, pTThr, sigma)
    with open(filename, 'w') as f:
        f.write(stringToWrite)


def randomNumber(min, max):
    return random.uniform(min, max)

#_____________________________________________________________________________
if __name__ == "__main__":



    for i in range(100):
        eps0 = 0
        n = randomNumber(0.423-0.005, 0.423+0.005)
        pTThr = randomNumber(0.42-0.05, 0.42 + 0.05)
        sigma = randomNumber(0.39 - 0.05, 0.39 + 0.05)
        writeToFile(eps0, n, pTThr, sigma)

        print("epsilon_0: ", eps0)
        print("n: ", n)
        print("pTThr: ", pTThr)
        print("sigma: ", sigma)

        path_to_embedding = "lists/embeddingJPsi.upcDst.root"
        cmd = "./AnalysisManager " + path_to_embedding
        if printOut: print(cmd)
        out = Popen(cmd.split(), stdout=PIPE).communicate()[0].split("\n")
        if printOut: print(out)

        # move the output to gpfs
        cmd = "cp AnalysisOutput.root " + gpfs_path + "StRP_production_0000.root"
        if printOut: print(cmd)
        out = Popen(cmd.split(), stdout=PIPE).communicate()[0]
        if printOut: print(out)

        os.chdir("../../plots/build")

        # run the plotting script
        cmd = "./RunPlot.py EmbeddingJPsi_BemcSysStudy"
        if printOut: print(cmd)
        out = Popen(cmd.split(), stdout=PIPE).communicate()[0]
        if printOut: print(out)

        os.chdir("../../star-upcDst/work")

        # copy the plots to an out dir
        cmd = "cp " + gpfs_path + "AnalysisOutput.root " + out_path + "AnalysisOutput" + str(i) + ".root"
        if printOut: print(cmd)
        out = Popen(cmd.split(), stdout=PIPE).communicate()[0]
        if printOut: print(out)

        '''
        dir_name = "pTThr_" + str(pTThr) + "_sigma_" + str(sigma)
        # Run the ROOT macro with the directory name
        cmd = 'root -l -b -q \'SaveResults.C("{}")\''.format(dir_name)
        out = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()[0]
        print(out)
        '''
        
    


