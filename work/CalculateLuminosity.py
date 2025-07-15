# macro to study luminosity in run 17 for CPT2noBBCL
# assume a root file per run

import math
import ROOT
import sys
from array import array

# Assume these are defined above or passed in
# path_to_rootfiles = ...
# lumi_map = ...
# total_lumi_from_Jaime = ...
# total_events_from_Jaime = ...

lumi_ZB_path = "/gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/luminosityForZB.list"

lumi_path = "/star/u/mvranovsk/star-upcDst/work/lists/lum_perrun_JPsi_HTTP.txt"  
file_path = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaGoodRun_6.5.25/merged/StRP_production_0000.root" 
#From the table: https://www.star.bnl.gov/protected/common/common2017/trigger2017/lumipp500GeV/lum_pertriggerid_pp2017_500GeV.html

total_lumi_from_Jaime = 330.220
total_events_from_Jaime = 140.084e+6


#parameters for final cross section
reco_efficiency = 0.1082
reco_error = 0.0012
trigger_efficiency = 0.42
RP_efficiency = 0.88 
vertexReco_efficiency = 0.9

#yields
yield_withRP = 110
err_withRP = 10
yield_withoutRP = 1580
err_withoutRP = 40

delta_rapidity = 2.0

def errors_in_cross_section(Y,y_err, lumi, lumi_err, eff, eff_err):

    err1 = (y_err/Y)/(lumi*eff*delta_rapidity) 
    err2 = -(Y*lumi_err/lumi)/(lumi*lumi*eff*delta_rapidity)
    err3 = -(Y*eff_err/eff)/(lumi*eff*eff*delta_rapidity)

    print "y_err: ", err1, " lumi_err: ", err2, " eff_err: ", err3

    return math.sqrt(y_err**2 + lumi_err**2 + eff_err**2)

def veto_correction(lumi):
    A = 1.04
    B = -0.011
    return A*math.exp(B*lumi)

def read_file_and_store_values(file_path):
    data_map = {}

    with open(file_path, 'r') as file:
        for line in file:
            values = line.split()
            run = int(values[0])
            lumi = float(values[4]) # [pb-1]
            events = int(float(values[11]))

            time0 = float(values[1])
            time1 = float(values[2])
            prescale = float(values[5])
            livetime = float(values[6])
            timeDiff = time1-time0
            inst_lumi = 0
            if livetime*timeDiff != 0:
                inst_lumi = prescale*lumi*1000000/(livetime*timeDiff); # convert from [pb-1] pico to [mub] mikro

            data_map[run] = (lumi, events, inst_lumi)

    return data_map


def load_events_from_rootfile(file_path):
    events_map = {}

    # Open the ROOT file
    root_file = ROOT.TFile(file_path)
    # Get the tree from the file
    tree = root_file.Get("RunInfo")  # Replace "tree" with the actual name of your tree
    # Loop over the entries in the tree 
    #load the branches to variables for ROOT 5

    mRunNumber = array('i', [0])
    nEventsAll = array('i', [0])
    nEventsPassed = array('i', [0])
    luminosity = array('d', [0])
    nEventsJPsi = array('i', [0])


    tree.SetBranchAddress("RunNumber", mRunNumber)
    tree.SetBranchAddress("nEventsAll", nEventsAll)
    tree.SetBranchAddress("nEventsPassed", nEventsPassed)
    tree.SetBranchAddress("luminosity", luminosity)
    tree.SetBranchAddress("nEventsJPsi", nEventsJPsi)

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        events_map[mRunNumber[0]] = ( nEventsPassed[0], nEventsAll[0], luminosity[0], nEventsJPsi[0] )

    root_file.Close()

    return events_map


def calculateCS(lumi_sum, yield_JPsi, withRP=True):

    eff = 1
    if withRP:
        eff = reco_efficiency * trigger_efficiency * RP_efficiency * vertexReco_efficiency
    else:
        eff = reco_efficiency * trigger_efficiency * vertexReco_efficiency
    result = yield_JPsi/(lumi_sum * eff * delta_rapidity)

    return result

total_integrated_luminosity = 0.0 # lumi from Jaime * events in root / events from Jaime
total_int_lum_perrun = 0.0 # run by run method
total_cor_int_lum_perrun = 0.0 # corrected by veto eff
total_cor_int_lum_perrun_min = 0.0 # corrected by veto eff
total_cor_int_lum_perrun_max = 0.0 # corrected by veto eff
total_int_lum_perrun = 0.0
total_uncorrected_lumi = 0.0
total_events = 0.0


if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) != 1:
        print "Usage: python3 CalculateLuminosity.py <run_number>"
        sys.exit(1)
    
    input_list = args[0]
    with open(input_list, 'r') as f:
        lines = f.readlines()
    
    lumi_map = read_file_and_store_values(lumi_path)
    print "lumi_map size: ", len(lumi_map)

    events_map = load_events_from_rootfile(file_path)

    for line in lines:
        run_number = int(line.split("/")[-1].replace(".root", ""))


        #check if the run number is in the map
        if run_number not in lumi_map:
            print "Run number not found in the map: ", run_number
            continue
        lumi, events, inst_lumi = lumi_map[run_number]

        if run_number not in events_map:
            print "Run number not found in the events map: ", run_number
            continue
        
        events_passed, events_all, luminosity, nEventsJPsi = events_map[run_number]
    
        if events == 0:
            print "No events in the run: ", run_number, ". nEventsJPsi: ", nEventsJPsi
            continue


        #inst_Lumi == luminosity
        if nEventsJPsi >1000000:
            print "A lot of events in the run: ", run_number, ". nEventsJPsi: ", nEventsJPsi
            continue


        #print "Run number: ", run_number, " Lumi: ", lumi, " Events: ", events, " Inst Lumi: ", inst_lumi, " luminosity: ", luminosity
        total_events += nEventsJPsi
        total_int_lum_perrun += lumi*nEventsJPsi/events
        total_cor_int_lum_perrun += lumi*veto_correction(inst_lumi)*nEventsJPsi/events
        total_cor_int_lum_perrun_min += lumi*(veto_correction(inst_lumi)-0.01)*nEventsJPsi/events
        total_cor_int_lum_perrun_max += lumi*(veto_correction(inst_lumi)+0.01)*nEventsJPsi/events
        total_uncorrected_lumi += lumi


    CS_events, CS, CS_min, CS_max = 0, 0, 0, 0 

    #error = 0.0
    eff_combined = 0
    if "noRP" in input_list:
        CS = calculateCS(total_cor_int_lum_perrun, yield_withoutRP, False)
        CS_min = calculateCS(total_cor_int_lum_perrun_min, yield_withoutRP, False)
        CS_max = calculateCS(total_cor_int_lum_perrun_max, yield_withoutRP, False)
        eff_combined = reco_efficiency * trigger_efficiency * vertexReco_efficiency
        #error = errors_in_cross_section(yield_withoutRP, err_withoutRP, total_cor_int_lum_perrun, 0.0, reco_efficiency*trigger_efficiency*vertexReco_efficiency, reco_error)
    else:
        CS = calculateCS(total_cor_int_lum_perrun, yield_withRP, True)
        CS_min = calculateCS(total_cor_int_lum_perrun_min, yield_withRP, True)
        CS_max = calculateCS(total_cor_int_lum_perrun_max, yield_withRP, True)
        eff_combined = reco_efficiency * trigger_efficiency * RP_efficiency * vertexReco_efficiency
        #error = errors_in_cross_section(yield_withRP, err_withRP, total_cor_int_lum_perrun, 0.0, reco_efficiency*trigger_efficiency*RP_efficiency*vertexReco_efficiency, reco_error)


    #print the results
    print "-------------------------------------------------------"
    #print "Total uncorected integrated luminosity: ", total_uncorrected_lumi , " pb^-1"
    print "total integrated luminosity per run: ", total_int_lum_perrun , " pb^-1"
    print "Total integrated luminosity from Jaime:" , total_lumi_from_Jaime*(total_events/total_events_from_Jaime), " pb^-1"
    print "Total events: ", total_events

    print "-------------------------------------------------------"
    print "Combined efficiency: ", eff_combined*100 , " %"  
    print "Total integrated luminosity corrected by veto correction: ", total_cor_int_lum_perrun , " pb^-1"
    print "Final cross section:", CS, " pb"
    print "-------------------------------------------------------"
    print "Final integrated luminosity min:", total_cor_int_lum_perrun_min, " pb^-1"
    print "Final cross section min:", CS_min , " pb"
    print "-------------------------------------------------------"
    print "Final integrated luminosity max:", total_cor_int_lum_perrun_max, " pb^-1"
    print "Final cross section max:", CS_max , " pb"
    print "-------------------------------------------------------"






    

    











