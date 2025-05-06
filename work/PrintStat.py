#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import sys
import os


#_____________________________________________________________________________
if __name__ == "__main__":

    #command line argumets
    resubmit = False

    #get the name of root file to process
    args = sys.argv
    args.pop(0) # cut first input from terminal = ./"macro".py
            
    defaultdir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/"

    for line in open("submit.xml").xreadlines(  ): 
        if "Location" in line:
            basedir = line.lstrip().split(">")[1].split("/sched")[0]
            location = basedir + "/merged/"

    if len(args) == 1:
        tag = str(args.pop(0))
        if tag == "r": 
            resubmit = True
        else:
            basedir = defaultdir + tag
            location = basedir + "/merged/"

    if len(args) == 2:
        basedir = defaultdir + str(args.pop(1))
        location = basedir + "/merged/"
        resubmit = True

    print basedir 
    #submitted jobs
    joblist = []
    for job in glob(basedir +"/sched/*_*.csh"):
        joblist.append( job.split("sched/sched")[1].split(".csh")[0] )

    print "Submitted:", len(joblist)

    #running jobs
    running = []
    cmd = "condor_q "+os.getlogin()
    out = Popen(cmd.split(), stdout=PIPE).communicate()[0].split("\n")
    for i in out:
        i1 = i.find("sched/sched")+len("sched/sched")
        i2 = i.find(".csh")
        if i2 < 0: continue
        jobid = i[i1:i2]
        #select jobs belonging to this production
        if jobid not in joblist: continue
        running.append(jobid)

    print "Running:", len(running)

    #done jobs
    donelist = []
    totsiz = 0
    #list all root files with size
    cmd = "ls -s " + basedir + "/*.root"
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")
    for fline in out:
        if len(fline) == 0: continue
        #get name and size for each file
        siznam = fline.lstrip().split(" ")
        size = int(siznam[0])
        name = siznam[1]
        #test for zero-sized outputs
        if size == 0:
            continue
        totsiz += size
        #job id from file name
        donelist.append(name.split("/")[-1][:-len(".root")])

    #missing jobs
    missing = []
    for job in joblist:
        if job in running or job in donelist:
            continue
        missing.append(job)
	#print job

    print "Errors:", len(missing)
    #print missing
    print "Done:", len(donelist)
    print "Output total size:", totsiz

    #resubmit missing jobs if requested
    if resubmit is True and len(missing) > 0:
        print "Resubmitting the missing jobs"
        ijob = 0
        os.chdir("sched")
        for job in missing:
            print "Clearing the log files"
            log_path = basedir + "/logs/"
            cmd = "rm " + log_path + job + ".out"
            out = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            cmd = "rm " + log_path + job + ".err"
            out = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            print "Job to resubmit:", job, "idx in list:", ijob
            session_num = job.split("_")
            cmd = "star-submit -r " + session_num[1] + " sched" + session_num[0] + ".session.xml"
            out = Popen(cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            print out[0]
            print out[1]
            ijob += 1
        print "Resubmittion done."



















