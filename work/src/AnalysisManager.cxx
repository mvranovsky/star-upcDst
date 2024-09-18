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
   if( runMCAna && runAnaV0SingleState ){
      mAnaVector.push_back(new AnaV0SingleState(outFile));
      cout << "Will run V0 analysis of MC simulation." << endl;
   }
   if( runAnaV0){
      mAnaVector.push_back(new AnaV0(outFile));
      cout << "Will run analysis AnaV0." << endl;
   }
   if( runAnaV0Control){
      mAnaVector.push_back(new AnaV0Control(outFile));
      cout << "Will run analysis AnaV0Control." << endl;
   }
   if( runAnaV0SingleState && !runMCAna){
      mAnaVector.push_back(new AnaV0SingleState(outFile));
      cout << "Will run analysis AnaV0SingleState." << endl;
   }


   // Load RP off-sets with off-set corrections 
   if( !LoadOffsetFile(nameOfOffSetFile, mOffSet) )
      return 3;
   #if !defined ALIGNMENT
      if( !LoadOffsetFile(nameOfOffSetCorrectionFile, mCorrection) )
         return 3;
   #endif
   cout<<"RP offset loaded..."<<endl;
   // Initiate histograms
   Init();


   //ask for number of events
   nEvents = upcTree->GetEntries();
   if (DEBUG)
      nEvents = 30000; // use for debugging and testing
   cout<<"Proccesing "<<nEvents<<" events"<<endl;

   for(Long64_t iev=0; iev<nEvents; ++iev) 
   { //get the event
      if( iev%1000 == 0)
         cout<<"Analyzing "<<iev<<". event "<<endl;
      upcTree->GetEntry(iev);
      mRunNumber = upcEvt->getRunNumber();
      if(!runMCAna)
         SetRpEvent();

      Make();

      //if( runRPMCANA )
         //RunRpMcAna();
      ReleaseTheMemoryOfCorrRpEvt();
   }
   /*
   if( runALIGNMENT ){
      RunAlignment();
      SaveAlignment(); 
   }*/
   
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
   
   //if( runRPMCANA )
      //mRpMCAna->Init();
}  
/*
void RunRpMcAna()
{
   if(!mRpMCAna->CheckTriggers(&ZBtriggers, upcEvt, nullptr))
      return;

   // Set RP correction
   TVector3 corr[nRomanPots];
   TVector3 offset[nRomanPots];
   for (int i = 0; i < nRomanPots; ++i){
      corr[i] = mCorrection[i][mRunNumber];
      offset[i] = mOffSet[i][mRunNumber];
   }
   mRpMCAna->SetRpPosition(corr, offset); 
   for (int iEmbed = 0; iEmbed < nMcEventsPerZbEvent; ++iEmbed)
   {
      if(iMcEvnt == nMcEvents)
         iMcEvnt = 0;

      mcTree->GetEntry(iMcEvnt++);

      mRpMCAna->SetEvent(upcEvt, rpEvt, mcEvt);
      mRpMCAna->SetMCInfo( mc_vrtx, mc_p);
      if(!mRpMCAna->AreTracksInRp())
         continue;

      mRpMCAna->Make();
   }

}

void RunAlignment()
{
   //The offset corrections are set at zero at the beginning. 
   //Thus, the afterburner for the first iteration keeps tracks unchanged.
   //This function is designed to be run on single run i.e. 1 job for 1 run
   ElasticAna *aligAna =  new ElasticAna(outFile);
   aligAna->Init();
   for (unsigned int iter = 0; iter < nAligIteration; ++iter)
   {
      for(Long64_t iev=0; iev<nEvents; ++iev) 
      { //get the event
         upcTree->GetEntry(iev);
         mRunNumber = upcEvt->getRunNumber();
         for (unsigned int iRP = 0; iRP < nRomanPots; ++iRP){
            mCorrection[iRP][mRunNumber][Z] = 0;
         }
         SetRpEvent(); 
         aligAna->SetEvent(upcEvt, rpEvt, mcEvt);
         aligAna->Make();
         ReleaseTheMemoryOfCorrRpEvt();
      }          
      for (unsigned int iRP = 0; iRP < nRomanPots; ++iRP){
         mCorrection[iRP][mRunNumber] += aligAna->CalculateOffsetCorrForRP(iRP, iter);
      }

   }
   delete aligAna;   
}

void SaveAlignment()
{

   TTree *aligTree = new TTree("aligTree", "aligTree");

   aligTree->Branch("runNumber", &mRunNumber);

   // RP event info
   for (int i = 0; i < nRomanPots; ++i)
   {
      aligTree->Branch("events_" + mUtil->rpName(i), &mCorrection[i][mRunNumber][Z]);
      aligTree->Branch("X_" + mUtil->rpName(i), &mCorrection[i][mRunNumber][X]);
      aligTree->Branch("Y_" + mUtil->rpName(i), &mCorrection[i][mRunNumber][Y]);
   }

   aligTree->Fill();
}
*/