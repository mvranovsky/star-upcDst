# CEP Analysis Repository

This repository contains code for measuring Time of Flight detector efficiency using tag/probe method. The code is based on [*star-upcDst*](https://github.com/adamjaro/star-upcDst), a new framework mainly developed by [Jarda Adam](https://github.com/adamjaro) and Tomas Truhlar to simplify analysis that are related to forward and UPC physics. 

The code consists of two parts:
- [*star-upcDst*](https://github.com/adamjaro/star-upcDst): Use MuDst and PicoDst files of one's interest as input, run so-called upcDstMaker to produce a upcDst file with information that are needed for one's analysis only. 
- work directory: Use the result of part 1, the upcDst, to build your personal analysis tree or histograms. 

## [*star-upcDst*](https://github.com/adamjaro/star-upcDst):

- Make a clean area on RACF and checkout the repository:

<pre><code> git clone https://github.com/mvranovsky/star-upcDst.git </pre></code>

- Go to the main directory star-upcDst:

<pre><code> cd star-upcDst </pre></code>

- Load a specific branch:
<pre><code> git checkout cleanCodes </pre></code>

- Setup the StRoot and build ( already include compiling ) by doing:

<pre><code> ./build-maker.sh </pre></code>

- To perform analysis on upcDst files, one needs to compile the package and setup the link to the library. This can be simply done in a few steps:

<pre><code> mkdir build </pre></code>
<pre><code> cd build </code></pre>
<pre><code> cmake ../ </code></pre>
<pre><code> make </code></pre>

- For the upcDst production, check the original [*star-upcDst*](https://github.com/adamjaro/star-upcDst) repository. But be aware that some changes were made. For example an option to merge old upcDst with PicoDst (for analyzing global tracks) was added. Also some scripts were modified.  

## work directory:
- This directory contains the code to fully reproduce V0(K0S, Lambda) analysis and calculate Time of Flight detector efficiency. This is done by using oppositely charged pion pairs for K0S reconstruction. 
- The analysis is set in "include/RunDef.h". It consists of several boolean operators based on which certain analysis is run. The default analysis is set to TofEffMult with which Time of Flight detector efficiency was measured. The difference between TofEff and TofEffMult is that in Mult, one looks for multiple pairs(5) of pions in a single event as compared to single pair in TofEff.
- Firstly, one has to go to the work dir and compile the code:

<pre><code> cd ../work </code></pre>
<pre><code> cmake .</code></pre>
<pre><code> make </code></pre>

### Run the code locally:
- To run the code locally on a single file: 
<pre><code> AnalysisManager /star/data01/pwg_tasks/upc03/pp17/ExtendedUpcDst/18100001.root </code></pre>
- To run the code locally on a list of files: 
<pre><code> AnalysisManager /gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/extendedUpc.list </code></pre>
- Or to run on the i-th file from list of files: 
<pre><code> AnalysisManager /gpfs01/star/pwg_tasks/upc02/CEP_input_files/lists/extendedUpc.list 111 </code></pre>
- The output will be stored in the work dir as "AnalysisOutput.root"

### Submitting jobs on condor:
- Before the jobs can be submitted, one has to set the output directory. The directory has to be set in "SubmitPlugin.py":
<pre><code> top = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic" </pre></code>
- And in "PrintStat.py" and "RunMerge.py":
<pre><code> defaultdir = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic" </pre></code>

- Also one has to create sched directory for the scheduler output
<pre><code> mdkir sched </pre></code>

- To submit the jobs on condor, run:
<pre><code> SubmitPlugin.py tag /gpfs01/star/u/mvranovsk/star-upcDst/work/lists/extendedUpc.list </pre></code>
- If the list is not given, then the default one from SubmitPlugin.py is used.
- The first argument (tag) is the name of the analysis. It should match string name of the analysis one is running, therefore when running TofEffMult, one could use: *TofEffMult_4.3.2025* with date. The name of the specific analysis is important for creating plots.

- One can check the job status using "PrintStat.py":
<pre><code> PrintStat.py tag </pre></code>
- If tag is not given, then the last submitted jobs are checked.
- If there are any failed jobs, one can resubmit them by:
<pre><code> PrintStat.py -r tag </pre></code>


- Each jobs produce an individual output. One can merge the output files by:
<pre><code> RunMerge.py tag </pre></code>
- Again, if tag is not given, then the last submitted jobs are merged.

### Create final plots:
- For creation of final plots, one has to use a different program. The code is available on GitHub:

<pre><code> https://github.com/mvranovsky/plots.git </pre></code>

## Contact:
Michal Vranovsky: <vranomic@fjfi.cvut.cz>








