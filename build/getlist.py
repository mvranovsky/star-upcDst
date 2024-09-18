#!/usr/bin/python


import glob
import sys


args = sys.argv
args.pop(0)
inputDirectory = args.pop(0)

print inputDirectory


out = glob.glob(inputDirectory + "/*/*.root")


for i in out:
	print i