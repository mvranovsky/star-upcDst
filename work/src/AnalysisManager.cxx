// Run by: ./AnalysisManager file.list
// e.g. ./AnalysisManager /gpfs01/star/pwg/truhlar/Final/CPtrig/merge_files/StUPCRP_production.list
// or you can open just root file
// ./AnalysisManager /star/data01/pwg_tasks/upc02/Part9/18143045/18143045.root
// or you can open n-th root file in file.list
// ./AnalysisManager file.list index 

#include "AnalysisManager.h"


//_____________________________________________________________________________
int main(int argc, char** argv) 
{
   cout<<"Starting the analysis..."<<endl;
   std::cout << "C++ version: " << __cplusplus << std::endl;
   //connect input file(s)
   if(!ConnectInput(argc, argv))
   {
      cout << "Wrong input parameters..." << endl; 
      return 1;
   }

   string inputPath = argv[1];
   //cout << "Input path: " << inputPath << endl;
   //open output file
   outFile = CreateOutputFile("AnalysisOutput.root"); 
   if(!outFile) 
   {
      cout << "Can not open output file." << endl; 
      return 2;
   }


   cout<<"Output file created..."<<endl;

   //if( runMAINANA )
      //mAnaVector.push_back(new MainAna(outFile));
   if( runAnaBP ){
      mAnaVector.push_back(new AnaBP(outFile));
      cout << "Will run analysis AnaBP." << endl;
   }
   if( runMCAna && runAnaV0 ){
      mAnaVector.push_back(new AnaV0(outFile));
      cout << "Will run V0 analysis of MC simulation." << endl;
   }
   if( !runMCAna && runAnaV0 ){
      mAnaVector.push_back(new AnaV0(outFile));
      cout << "Will run analysis AnaV0." << endl;
   }
   if( runAnaV0Mult ){
      mAnaVector.push_back(new AnaV0Mult(outFile));
      cout << "Will run analysis AnaV0SingleState." << endl;
   }
   if( runTofEff ){
      mAnaVector.push_back(new TofEff(outFile));
      cout << "Will run analysis ToF efficiency." << endl; 
   }
   if( !runMCAna && runTofEffMult ){
      mAnaVector.push_back(new TofEffMult(outFile));
      cout << "Will run analysis ToF efficiency with multiple states." << endl;       
   }
   if( runMCAna && runTofEffMult ){
      mAnaVector.push_back(new TofEffMult(outFile));
      cout << "Will run analysis MC ToF efficiency with multiple states." << endl;       
   }
   if( runAnaJPsi ){
      mAnaVector.push_back(new AnaJPsi(outFile));
      cout << "Will run analysis of JPsi." << endl;       
   }
   if( runEmbeddingJPsi ){
      mAnaVector.push_back(new EmbeddingJPsi(outFile));
      cout << "Will run embedding of JPsi." << endl;
   }
   if( runAnaJPSI ){
      mAnaVector.push_back(new AnaJPSI(outFile));
      cout << "Will run analysis of JPSI tryout." << endl;       
   }
   if( runAnaGoodRun ){

      mAnaVector.push_back(new AnaGoodRun(outFile));
      cout << "Will run analysis of good runs." << endl;
   }

   
   // Load RP off-sets with off-set corrections 
   if( !LoadOffsetFile(offsetFilePath, mOffSet) )
      return 3;
   #if !defined ALIGNMENT
      if( !LoadOffsetFile(offsetCorrectionsFilePath, mCorrection) )
         return 3;
   #endif
   cout<<"RP offset loaded..."<<endl;
   // Initiate histograms
   
   Init();


   //ask for number of events
   nEvents = upcTree->GetEntries();
   if (DEBUG)
      nEvents = 30000; // use for debugging and testing
   cout<<"Processing "<<nEvents<<" events"<<endl;

   //check if the run number is in the list of bad runs
   for(Long64_t iev=0; iev<nEvents; ++iev) 
   { //get the event
      //if( iev%1000 == 0)
      //   cout<<"Analyzing "<<iev<<". event "<<endl;
      upcTree->GetEntry(iev);
      
      if(iev == 0){
         cout << "This is MC: " << upcEvt->getIsMC() << endl;
      }
      
      mRunNumber = upcEvt->getRunNumber();
      // check if the run number is in the list of bad runs
      //check if a number is in std::vector<int>
      
      if(!(runMCAna || runEmbeddingJPsi) ){
         SetRpEvent();
      }
      
      Make();


      ReleaseTheMemoryOfCorrRpEvt();
   }
   
   
   if( runAnaGoodRun ){
      analysisOfTracks(mRunNumber, nEvents,inputPath);
      AnaGoodRun* anaGoodRun = dynamic_cast<AnaGoodRun*>(mAnaVector[0]);
      if( anaGoodRun ){
         fillRunTree(anaGoodRun, inputPath);
      }
      

   }
   
   //close the outputs
   //CleanMemory();
   outFile->Write(0, TObject::kOverwrite);
   outFile->Close();

   cout<<"Ending Analysis... GOOD BYE!"<<endl;
   return 0;
}//main

void Make()
{
   for (unsigned int i = 0; i < mAnaVector.size(); ++i)
   {
      mAnaVector[i]->SetEvent(upcEvt, rpEvt, mcEvt);
      mAnaVector[i]->Make();
   }
}


void Init()
{
   mUtil = new Util();
   for (unsigned int i = 0; i < mAnaVector.size(); ++i){
      mAnaVector[i]->Init();
   }

}  

void analysisOfTracks(int RUNNUMBER, int nEvents, string inputPath){

   ofstream runInfoFile("RunInfo.txt");
   AnaGoodRun* anaGoodRun = dynamic_cast<AnaGoodRun*>(mAnaVector[0]);
   if( anaGoodRun ){

      if(nEvents == 0){
         cout << "No events were processed!" << endl;
         // load the single argument from the command line, it is the path to the file. It looks like /path/to/file/RUNNUMBER.root

         int runNumber = stoi(inputPath.substr(inputPath.find_last_of("/")+1, inputPath.find_last_of(".") - inputPath.find_last_of("/") - 1));
         runInfoFile << runNumber << endl;  //means no events were processed
         runInfoFile << "-1" << endl;  
         runInfoFile.close();
         return;
      }else{
         runInfoFile << RUNNUMBER << endl;
      }

      bool RPsAreClose = areRPsCloseEnough(RUNNUMBER);
      if(RPsAreClose){
         runInfoFile << "1" << endl;  //means the run is good (RPs close)
      }else{
         runInfoFile << "0" << endl;  //means the run is bad (RPs far)
      }

      runInfoFile << fixed << setprecision(6) << anaGoodRun->isJPsiTrigger() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageBemcTracks() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageBemcClusters() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageTpcTracks() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageTOFTracks() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageVertices() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageTpcEta() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageBemcEta() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageTpcPhi() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getAverageBemcPhi() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getNEventsPassed() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getNEventsAll() << endl;
      runInfoFile << fixed << setprecision(6) << anaGoodRun->getLuminosity(RUNNUMBER) << endl;

      cout << "Finished filling up RunInfo.txt" << endl;

      runInfoFile.close();
      return;
   }else{
      cout << "No AnaGoodRun object found!" << endl;
      runInfoFile << "No AnaGoodRun object found!" << endl;
      runInfoFile.close();
      return;
   }

   
}

void fillRunTree(AnaGoodRun *Ana, string inputPath){

   int mRunNumber = 0;
   int atLeast1JPsiTrigger, RPsClose;
   double nTracksBemc, nClustersBEMC, nTracksTPC, nTracksTOF, nVertices;
   double tpcEtaAverage, bemcEtaAverage, tpcPhiAverage, bemcPhiAverage;
   int nEventsAll, nEventsPassed, nEventsJPsi;
   double luminosity, luminosityError;

   TTree *runTree = new TTree("RunInfo", "RunInfo");
   runTree->Branch("RunNumber", &mRunNumber, "RunNumber/I");
   runTree->Branch("AtLeast1JPsiTrigger", &atLeast1JPsiTrigger, "AtLeast1JPsiTrigger/I");
   runTree->Branch("RPsClose", &RPsClose, "RPsClose/I");
   runTree->Branch("nEventsAll", &nEventsAll, "nEventsAll/I");
   runTree->Branch("nEventsPassed", &nEventsPassed, "nEventsPassed/I");
   runTree->Branch("nEventsJPsi", &nEventsJPsi, "nEventsJPsi/I");
   runTree->Branch("luminosity", &luminosity, "luminosity/D");
   runTree->Branch("luminosityError", &luminosityError, "luminosityError/D");
   runTree->Branch("nTracksBEMC", &nTracksBemc, "nTracksBEMC/D");
   runTree->Branch("nClustersBEMC", &nClustersBEMC, "nClustersBEMC/D");
   runTree->Branch("nTracksTPC", &nTracksTPC, "nTracksTPC/D");
   runTree->Branch("nTracksTOF", &nTracksTOF, "nTracksTOF/D");
   runTree->Branch("nVertices", &nVertices, "nVertices/D");
   runTree->Branch("tpcEtaAverage", &tpcEtaAverage, "tpcEtaAverage/D");
   runTree->Branch("bemcEtaAverage", &bemcEtaAverage, "bemcEtaAverage/D");
   runTree->Branch("tpcPhiAverage", &tpcPhiAverage, "tpcPhiAverage/D");
   runTree->Branch("bemcPhiAverage", &bemcPhiAverage, "bemcPhiAverage/D");

   if( upcTree->GetEntries() <= 0 ){
      // convert inputPath to int
      // inputPath looks like /path/to/file/RUNNUMBER.root
      // get the last part of the path, which is the file name
      // get the run number from the file name
      mRunNumber = stoi(inputPath.substr(inputPath.find_last_of("/")+1, inputPath.find_last_of(".") - inputPath.find_last_of("/") - 1));  //means no events were processed

      atLeast1JPsiTrigger = 0;
      RPsClose = 0;
      nEventsAll = 0;
      nEventsPassed = 0;
      nEventsJPsi = 0;
      luminosity = 0.0;
      luminosityError = 0.0; // not used yet
      nTracksBemc = 0.0;
      nClustersBEMC = 0.0;
      nTracksTPC = 0.0;
      nTracksTOF = 0.0;
      nVertices = 0.0;
      tpcEtaAverage = 0.0;
      bemcEtaAverage = 0.0;
      tpcPhiAverage = 0.0;
      bemcPhiAverage = 0.0;
      runTree->Fill();

   }else{
      mRunNumber = upcEvt->getRunNumber();
      if(Ana->isJPsiTrigger()){
         atLeast1JPsiTrigger = 1;
      }else{
         atLeast1JPsiTrigger = 0;
      }
      if(areRPsCloseEnough(mRunNumber)){
         RPsClose = 1;
      }else{
         RPsClose = 0;
      }
      nEventsAll = Ana->getNEventsAll();
      nEventsPassed = Ana->getNEventsPassed();
      nEventsJPsi = Ana->getJPsiTriggerEvents();
      luminosity = Ana->getLuminosity(mRunNumber);
      luminosityError = 0.0; // not used yet
      nTracksBemc = Ana->getAverageBemcTracks();
      nClustersBEMC = Ana->getAverageBemcClusters();
      nTracksTPC = Ana->getAverageTpcTracks();
      nTracksTOF = Ana->getAverageTOFTracks();
      nVertices = Ana->getAverageVertices();
      tpcEtaAverage = Ana->getAverageTpcEta();
      bemcEtaAverage = Ana->getAverageBemcEta();
      tpcPhiAverage = Ana->getAverageTpcPhi();
      bemcPhiAverage = Ana->getAverageBemcPhi();
      runTree->Fill();
      
   }
   
   outFile->cd();
   runTree->Write(0, TObject::kOverwrite);

}

