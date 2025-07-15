#!/usr/bin/python

from subprocess import Popen, PIPE
import os

if __name__ == "__main__":
    # Define the command to run
    
    
    top = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaGoodRun_6.5.25"
    

    dirs = [dir for dir in os.listdir(top) if dir.endswith(".root")]
    for dir in dirs:
        # Construct the command
        print(top + "/" + dir)
        command = "root -l -b -q " + top  +"/" + dir 
        
        # Execute the command and capture the output
        process = Popen(command, shell=True, stdout=PIPE, stderr=PIPE).communicate()
        if "Error" in process[1]:
            print(process[1])

            out_path = top + "/logs/" + dir.replace(".root", ".out")
            print("out_path: " + out_path)

            with open(out_path, "r") as out_file:
                lines = out_file.readlines()
                for line in lines:
                    if "Running: cp " in line:
                        input_file = line.split("cp ")[1].split(" ")[0]
                        print("input_file: " + input_file)
                        break
            
            cmd = "./AnalysisManager " + input_file
            print("cmd: " + cmd)
            out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()
            cmd = "mv AnalysisOutput.root " + top + "/" + dir
            print("cmd: " + cmd)
            out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()
        else:
            print("No error in " + dir)
            continue






    

    
