#!/usr/bin/python

from subprocess import check_output
from math import ceil
from subprocess import Popen
import os

query = [
    #"production=P16id,trgsetupname=AuAu_200_production_2014,filetype=daq_reco_MuDst,filename~st_upc,storage=",
    #"production=P16id,trgsetupname=AuAu_200_production_low_2014,filetype=daq_reco_MuDst,filename~st_upc,storage=",
    #"production=P16id,trgsetupname=AuAu_200_production_mid_2014,filetype=daq_reco_MuDst,filename~st_upc,storage=",
    #"production=P16id,trgsetupname=AuAu_200_production_high_2014,filetype=daq_reco_MuDst,filename~st_upc,storage="
    #"production=P18ib,trgsetupname=pp500_production_2017,filetype=daq_reco_MuDst,filename~st_rp,storage="
    "production=P20ic,trgsetupname=pp500_production_2017,filetype=daq_reco_picoDst,filename~st_rp,storage!="
]

#storage="local"
#storage="nfs"
storage="hpss"

get_list = ["get_file_list.pl", "-keys", "path,filename,size", "-limit", "0", "-cond"]

out = []
for i in query:
    out += check_output(get_list + [i+storage]).split("\n")



#remove empty entries (last one is empty)
out = [i for i in out if i]



qlist = []

for i in out:
    i = i.replace("_adc", "")
    q = i.split("::")[1]
    q = q.split("/")[-1]
    qlist.append(q.split("_")[2])


qlist = list(set(qlist))


directory = "/gpfs01/star/pwg/mvranovsk/run17PicoOut/"
done_list = os.listdir(directory)

done_list.pop(0)

for num in done_list:
    num.replace("/","")
    

qlist = [item for item in qlist if item not in done_list]


for q in qlist:
    q = "add_input " + q + " production=P20ic,trgsetupname=pp500_production_2017,filetype=daq_reco_picoDst,filename~st_rp,storage!=hpss,runnumber="
    print q
with open("config_picoMissing.in") as file:
    file.write(qlist)



'''

print len(out_done)
out_done = [i for i in out_done if i]
out_done.pop(0)
for i in out_done:
    print i


#------------------------------------------------------
flist_all = []
totsiz = 0.
for i in out:
    il = i.split("::")
    i.replace('::', "/")
    #i.replace('localhost/', 'root://xrdstar.rcf.bnl.gov:1095')
    file = il[0] + "/" + il[1]
    print file
    flist_all.append( file )
    totsiz += float(il[2])

#remove dupplicate files
flist = list( dict.fromkeys(flist_all) )

nall = len(flist_all)
nf = len(flist)

print "All files:", nall
print "After removing duplicates:", nf
print "Duplicate files removed:", nall-nf

print "Total size:", ceil(totsiz/1024**4), "T"

#print out
#print flist

#for i in flist: print i

'''
