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
      //if( iev%1000 == 0)
      //   cout<<"Analyzing "<<iev<<". event "<<endl;
      upcTree->GetEntry(iev);
      mRunNumber = upcEvt->getRunNumber();
      if(!runMCAna)
         SetRpEvent();

      Make();

      ReleaseTheMemoryOfCorrRpEvt();
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
