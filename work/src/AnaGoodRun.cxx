#include "AnaGoodRun.h"

//_____________________________________________________________________________
AnaGoodRun::AnaGoodRun(TFile *outFile): Ana(outFile){}

AnaGoodRun::~AnaGoodRun(){
   //if(mUtil) delete mUtil;
}


void AnaGoodRun::Make(){

   tracksBEMC.clear();
   tracksTPC.clear();


   //trigger
   hAnalysisFlow->Fill(GRALL);
   // just fill hTriggerBits
   CheckTriggers(&JPSItriggers, mUpcEvt, hTriggerBits);
   bool atLeast1Trigger = false;
   if(mUpcEvt->isTrigger(JPSItriggers[0])){
      mRecTree->setJPsiTrigger1(1);
      atLeast1Trigger = true;
      hAnalysisFlow->Fill(GRTRIGGER);
   }
   if(mUpcEvt->isTrigger(JPSItriggers[1])){
      mRecTree->setJPsiTrigger2(1);
      atLeast1Trigger = true;
      hAnalysisFlow->Fill(GRTRIGGER);
   }
   if(mUpcEvt->isTrigger(JPSItriggers[2])){
      mRecTree->setJPsiTrigger3(1);
      hAnalysisFlow->Fill(GRTRIGGER);
      atLeast1Trigger = true;
   }

   SaveEventInfo(mUpcEvt);
   
   //if(!atLeast1Trigger){
      //   mRecTree->FillRecTree();
      //cout << "No JPsi trigger found!" << endl;
      //   return;
      //}
      
   int RUNNUMBER = mUpcEvt->getRunNumber();
   mRecTree->setRpOk(1);
   // check if offsets are loaded
   if(offsets[0].find(RUNNUMBER) == offsets[0].end()){
      mRecTree->setRpOk(-1);
      //mRecTree->FillRecTree();
      //cout << "RunNumber: " << RUNNUMBER << " not found in offsets file!" << endl;
      //return;
   }

   // check if RPs are close enough
   for (int iRP = 0; iRP < nRomanPots; ++iRP){
      double offsetX = offsets[iRP][RUNNUMBER].X() + corrections[iRP][RUNNUMBER].X();
      double offsetY = offsets[iRP][RUNNUMBER].Y() + corrections[iRP][RUNNUMBER].Y();

      if( offsetY > 0.04){
         mRecTree->setRpOk(-1);
         //mRecTree->FillRecTree();
         //cout << "Offsets for RP" << iRP << " are too big!" << endl;
         //return;
      }
   }



   // track loop
   int nTofTracks = 0;
   vector<Double_t> tpcEta, tpcPhi, tpcNHitsFit, tpcNHitsDEdx, tpcNSigmaP, tpcNSigmaPi, tpcNSigmaK, tpcNSigmaE, bemcEta, bemcPhi;
   for (int iTrk = 0; iTrk < mUpcEvt->getNumberOfTracks(); ++iTrk){
      const StUPCTrack *trk = mUpcEvt->getTrack(iTrk);

      if(!trk->getFlag(StUPCTrack::kPrimary)) // original star-upcDst
         continue;

      if(!goodQualityTrack(trk))
         continue;
      
      tracksTPC.push_back(iTrk);
      hAnalysisFlow->Fill(GRGOODTRACKTPC);

      tpcEta.push_back(trk->getEta());
      tpcPhi.push_back(trk->getPhi());
      tpcNHitsFit.push_back(trk->getNhitsFit());
      tpcNHitsDEdx.push_back(trk->getNhitsDEdx());
      tpcNSigmaP.push_back(trk->getNSigmasTPCProton());
      tpcNSigmaPi.push_back(trk->getNSigmasTPCPion());
      tpcNSigmaK.push_back(trk->getNSigmasTPCKaon());
      tpcNSigmaE.push_back(trk->getNSigmasTPCElectron());


      if( !trk->getFlag(StUPCTrack::kTof) )
         nTofTracks++;
      if( !trk->getFlag(StUPCTrack::kBemc) )
         continue;

      hAnalysisFlow->Fill(GRGOODTRACKBEMC);
      tracksBEMC.push_back(iTrk);
      bemcEta.push_back(trk->getBemcEta());
      bemcPhi.push_back(trk->getBemcPhi());
   }

   //if( tracksTPC.size() == 0 ){
   //   mRecTree->FillRecTree();
      //cout << "No good TPC tracks found!" << endl;
   //   return;
   //}

   //fill tree with info
   mRecTree->setNTracksBemc(tracksBEMC.size());
   mRecTree->setNTracksTof(nTofTracks);
   mRecTree->setNClustersBemc(mUpcEvt->getNumberOfClusters());
   mRecTree->setTpcTrackPhi(tpcPhi);
   mRecTree->setTpcTrackEta(tpcEta);
   mRecTree->setBemcTrackPhi(bemcPhi);
   mRecTree->setBemcTrackEta(bemcEta);
   mRecTree->setTpcNHitsFit(tpcNHitsFit);
   mRecTree->setTpcNHitsDEdx(tpcNHitsDEdx);
   mRecTree->setTpcNSigmaProton(tpcNSigmaP);
   mRecTree->setTpcNSigmaPion(tpcNSigmaPi);
   mRecTree->setTpcNSigmaElectron(tpcNSigmaE);
   mRecTree->setTpcNSigmaKaon(tpcNSigmaK);

   mRecTree->FillRecTree();


   if(DEBUG){
      cout << "Finished AnaGoodRun::Make()" << endl;
   }
}



void AnaGoodRun::Init(){

   mUtil = new Util();


   if( DEBUG )
     cout<<"AnaGoodRun::Init() called"<<endl;

   mOutFile->cd();
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nGRCuts-1, 1, nGRCuts);
   for(int tb=1; tb<nGRCuts; ++tb) {
     hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisGoodRun(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
     TString label; label.Form("%d",triggerID[tb]);
     hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }

   mRecTree = new RecTree(nameOfAnaGoodRunTree, AnaGoodRunTreeBits, false); 


   // load offset file
   if(!LoadOffsetFile(offsetFilePath, offsets )){
      cout << "Could not load file with offsets!" << endl;
      return;
   }

   // load offset correction file 
   if(!LoadOffsetFile(offsetCorrectionsFilePath, corrections ) ){
      cout << "Could not load file with offset corrections!" << endl;
      return;
   }



   cout << "Finished AnaGoodRun::Init()" << endl;
}

bool AnaGoodRun::goodQualityTrack(const StUPCTrack *trk){

   
   if(!(abs(trk->getBemcEta()) < maxEta))  // eta
      return false;

   if( !(trk->getDcaZ() > minDcaZ && trk->getDcaZ() < maxDcaZ) )  //DCA z
      return false;

   if( !(trk->getDcaXY() < maxDcaXY) ) //DCA xy
      return false;

   if( !(trk->getNhitsFit() > minNHitsFit) )  //NhitsFit
      return false;

   if( !(trk->getNhitsDEdx() > minNHitsDEdx) ) //NhitsdEdx
      return false;

   return true;

}


bool AnaGoodRun::LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]){

   string line;
   ifstream file( fileName );
   if (!file.is_open() )
   {
      cout << "\n ERROR in LoadOffsets(): Problems with opening a file: "<< fileName<< endl;
      return false;
   }   
   int RUNNUMBER;
   double var[2*nRomanPots]; // placeholder for offsets or offsets correction in x and y
   while ( getline(file, line) )
   {
      stringstream ss(line);
      ss >> RUNNUMBER;
      for (int i = 0; i < 2*nRomanPots; ++i)
         ss >> var[i]; 

      for (int iRP = 0; iRP < nRomanPots; ++iRP)
         offsets[iRP][RUNNUMBER].SetXYZ(var[2*iRP], var[2*iRP + 1], 0.0); 
   }

   file.close();
   return true;
}
