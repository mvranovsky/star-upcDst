#include "AnaZeroBias.h"

//_____________________________________________________________________________
AnaZeroBias::AnaZeroBias(TFile *outFile): Ana(outFile){}

AnaZeroBias::~AnaZeroBias(){
   //if(mUtil) delete mUtil;
}


void AnaZeroBias::Make(){


   nEventsAll++;

   //----------------------------------------------------------------
   // here begins analysis of zerobias data
   hAnalysisFlow->Fill(ZBALL);
   if(!CheckTriggers(&ZBtriggers,mUpcEvt,hTriggerBits)) {
      return;
   }
   hAnalysisFlow->Fill(ZBTRIGGER);

   addZBEventAll();

   SaveEventInfo(mUpcEvt);
   SaveTriggerInfo(mUpcEvt, mRpEvt);

   hNVertices->Fill(mUpcEvt->getNumberOfVertices());
   
   if(mUpcEvt->getNumberOfVertices() != 1){
      return;
   }

   StUPCVertex* vtx = mUpcEvt->getVertex(0);
   if(!vtx) return;

   SaveVertexInfo(vtx, 0);
   // find at least 2 ToF tracks
   vector<int> tofTracks, bemcTracks, tpcTracks;

   for(int i = 0; i < mUpcEvt->getNumberOfTracks(); ++i){ // outer loop
      const StUPCTrack *trk1 = mUpcEvt->getTrack(i);
      if(!trk1->getFlag(StUPCTrack::kPrimary)) continue;

      if(trk1->getVertexId() != 0) continue;

      if(!goodQualityTrack(trk1)) continue;

      tpcTracks.push_back(i);

      if(trk1->getFlag(StUPCTrack::kBemc)) {
         bemcTracks.push_back(i);
      }

      if(trk1->getFlag(StUPCTrack::kTof)) {
         tofTracks.push_back(i);
      }
   }

   hNTracksTof->Fill(tofTracks.size());
   hNTracksBemc->Fill(bemcTracks.size());
   
   mRecTree->setNTracksTof(tofTracks.size());
   mRecTree->setNTracksBemc(bemcTracks.size());
   mRecTree->setNTracksTpc(tpcTracks.size());

   if(bemcTracks.size() < 2) return;

   hAnalysisFlow->Fill(ZBBEMCTRACKS);
   
   //back to back tracks
   pair<int,int> pair = make_pair(-1, -1);
   for(unsigned int i = 0; i < bemcTracks.size(); ++i){
      const StUPCTrack *trk1 = mUpcEvt->getTrack(bemcTracks[i]);
      if(!trk1) continue;
      for(unsigned int j = i + 1; j < bemcTracks.size(); ++j){
         const StUPCTrack *trk2 = mUpcEvt->getTrack(bemcTracks[j]);
         if(!trk2) continue;

         if(backToBack(trk1, trk2, false)) {
            pair = make_pair(bemcTracks[i], bemcTracks[j]);
            break;
         } 
      } // inner loop
      if(pair.first != -1 && pair.second != -1) break;
   } // outer loop

   if( pair.first == -1 || pair.second == -1 ) return;

   const StUPCTrack *track1 = mUpcEvt->getTrack(pair.first);
   const StUPCTrack *track2 = mUpcEvt->getTrack(pair.second);

   if(!track1 || !track2){
      cout << "Error: tracks not found in the event!" << endl;
      return;
   }
   hAnalysisFlow->Fill(ZBBACKTOBACK);

   //PID
   if(!goodPID(track1, track2)) return;

   hAnalysisFlow->Fill(ZBPID);
   SaveChiSquareInfo(track1, track2);

   // eta vtxz cut
   if( abs(vtx->getPosZ()) > vertexRange) return;

   if( !IsGoodEtaTrack(track1, 0) || !IsGoodEtaTrack(track2, 0)) return; 

   hAnalysisFlow->Fill(ZBVTXZETA);

   // save info about tracks
   TLorentzVector electron1, electron2, state;
   track1->getLorentzVector(electron1, mUtil->mass(ELECTRON));
   track2->getLorentzVector(electron2, mUtil->mass(ELECTRON));
   state = electron1 + electron2;
   int totQ = track1->getCharge() + track2->getCharge();
   SaveTrackInfo(track1, electron1, 0);
   SaveTrackInfo(track2, electron2, 1);
   SaveStateInfo(state, totQ, 0);


   // check if zerobias satisfies back to back signals in BEMC
   UShort_t dsm3 = mUpcEvt->getLastDSM3();
   Bool_t emc_http = dsm3 & (1 << 3);

   
   if(emc_http == true){
      hAnalysisFlow->Fill(ZBTRIGGERCONDITION);
      mRecTree->setIsTopology(1);
   }else{
      mRecTree->setIsTopology(-1);
   }
   
   mRecTree->FillRecTree();


   if(DEBUG){
      cout << "Finished AnaZeroBias::Make()" << endl;
   }
}



void AnaZeroBias::Init(){

   mUtil = new Util();
   
   
   if( DEBUG )
   cout<<"AnaZeroBias::Init() called"<<endl;
   
   mOutFile->cd();
   mRecTree = new RecTree(nameOfAnaZeroBiasTree, ZeroBiasTreeBits, false);
   mOutFile->mkdir(nameOfAnaZeroBiasDir);
   mOutFile->cd(nameOfAnaZeroBiasDir);

   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nZBCuts-1, 1, nZBCuts);
   for(int tb=1; tb<nZBCuts; ++tb) {
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisZeroBias(tb));
   }
   
   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   
   nEventsZBAll = 0;
   nEventsAll = 0;
   nEventsJPsi = 0; 

   hNTracksTof = new TH1D("hNTracksTof", "Number of ToF Tracks", 11, -0.5, 10.5);
   hNTracksBemc = new TH1D("hNTracksBemc", "Number of BEMC Tracks", 6, -0.5, 5.5);
   hNVertices = new TH1D("hNVertices", "Number of Vertices", 6, -0.5, 5.5);

   cout << "Finished AnaZeroBias::Init()" << endl;

}


bool AnaZeroBias::goodQualityTrack(const StUPCTrack *trk){

   
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


bool AnaZeroBias::goodPID(const StUPCTrack *trk1, const StUPCTrack *trk2){

   double chiSquare = pow(trk1->getNSigmasTPCElectron(),2);
   chiSquare += pow(trk2->getNSigmasTPCElectron(),2);

   if( chiSquare > maxPidChiEE ) {
      return false;
   }else{
      return true;
   }

}
