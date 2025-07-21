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
   if( runSysStudyEmbedding){
      mAnaVector.push_back(new EmbeddingJPsi(outFile));
      cout << "Will run systematics study of JPsi." << endl;
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
   if( runSysStudy ){
      mAnaVector.push_back(new AnaJPsi(outFile));
      cout << "Will run systematics study." << endl;
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

      if(!(runMCAna || runEmbeddingJPsi || runSysStudyEmbedding) ){
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


   outFile->cd();
   RecTree* mRecTree = new RecTree(nameOfAnaGoodRunTree, AnaGoodRunTreeBits, false);


   if( upcTree->GetEntries() <= 0 ){
      // convert inputPath to int
      // inputPath looks like /path/to/file/RUNNUMBER.root
      // get the last part of the path, which is the file name
      // get the run number from the file name
      mRunNumber = stoi(inputPath.substr(inputPath.find_last_of("/")+1, inputPath.find_last_of(".") - inputPath.find_last_of("/") - 1));  //means no events were processed

      mRecTree->setRunNumber(mRunNumber);
      mRecTree->setAtLeast1JPsiTrigger(0);
      mRecTree->setRPsClose(0);
      mRecTree->setNEventsAll(0);
      mRecTree->setNEventsPassed(0);
      mRecTree->setNEventsJPsi(0);
      mRecTree->setLuminosity(0.0);
      mRecTree->setLuminosityError(0.0);
      mRecTree->setNTracksBemc(0.0);
      mRecTree->setNClustersBemc(0.0);
      mRecTree->setNTracksTpc(0.0);
      mRecTree->setNTracksTof(0.0);
      mRecTree->setNVertices(0.0);
      mRecTree->setTpcEtaAverage(0.0);
      mRecTree->setBemcEtaAverage(0.0);
      mRecTree->setTpcPhiAverage(0.0);
      mRecTree->setBemcPhiAverage(0.0);


   }else{
      mRunNumber = upcEvt->getRunNumber();
      if(Ana->isJPsiTrigger()){
         mRecTree->setAtLeast1JPsiTrigger(1);
      }else{
         mRecTree->setAtLeast1JPsiTrigger(0);
      }
      if(areRPsCloseEnough(mRunNumber)){
         mRecTree->setRPsClose(1);
      }else{
         mRecTree->setRPsClose(0);
      }
      mRecTree->setNEventsAll(Ana->getNEventsAll());
      mRecTree->setNEventsPassed(Ana->getNEventsPassed());
      mRecTree->setNEventsJPsi(Ana->getJPsiTriggerEvents());
      mRecTree->setLuminosity(Ana->getLuminosity(mRunNumber));
      mRecTree->setLuminosityError(0.0); // not used yet
      mRecTree->setNTracksBemc(Ana->getAverageBemcTracks());
      mRecTree->setNClustersBemc(Ana->getAverageBemcClusters());
      mRecTree->setNTracksTpc(Ana->getAverageTpcTracks());
      mRecTree->setNTracksTof(Ana->getAverageTOFTracks());
      mRecTree->setNVertices(Ana->getAverageVertices());
      mRecTree->setTpcEtaAverage(Ana->getAverageTpcEta());
      mRecTree->setBemcEtaAverage(Ana->getAverageBemcEta());
      mRecTree->setTpcPhiAverage(Ana->getAverageTpcPhi());
      mRecTree->setBemcPhiAverage(Ana->getAverageBemcPhi());
      
   }

   mRecTree->FillRecTree();


}

