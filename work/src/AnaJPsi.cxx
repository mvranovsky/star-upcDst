#include "AnaJPsi.h"

//_____________________________________________________________________________
AnaJPsi::AnaJPsi(TFile *outFile): Ana(outFile){}

AnaJPsi::~AnaJPsi(){
   //if(mUtil) delete mUtil;
}

void AnaJPsi::Make(){

   //trigger
   hAnalysisFlow->Fill(JPSIALL);
   if( !CheckTriggers(&JPSItriggers, mUpcEvt, hTriggerBits) ) 
      return;
   hAnalysisFlow->Fill(JPSITRIG);

   SaveEventInfo(mUpcEvt);

   SaveTriggerInfo(mUpcEvt, mRpEvt);

   // 1 vertex
   hNVertices->Fill(mUpcEvt->getNumberOfVertices());
   if(mUpcEvt->getNumberOfVertices() != 1){
      return;
   }

   hAnalysisFlow->Fill(JPSI1VTX);

   // fill vertex histogram by fill for systematic study
   int fillNum = mUpcEvt->getFillNumber();
   hVtxZByFillNum->SetName(Form("hVtxZFillNum_%d", fillNum));
   hVtxZByFillNum->Fill(mUpcEvt->getVertex(0)->getPosZ());

   const StUPCVertex *vtx = mUpcEvt->getVertex(0);

   if(!vtx){
      cout << "Could not load vertex. Leaving event." << endl;
      return;
   }


   unsigned int vertexID = vtx->getId();

   tpcCounter = 0;
   tracksBEMC.clear();
   //central system good quality tracks + BEMC
   for (int iTrk = 0; iTrk < mUpcEvt->getNumberOfTracks(); ++iTrk){
      hTrackQualityFlow->Fill(1);
      const StUPCTrack *trk = mUpcEvt->getTrack(iTrk);

      if(!trk->getFlag(StUPCTrack::kPrimary)) // original star-upcDst
         continue;
      if(trk->getVertexId() != vertexID){ // only belonging to the 1 vertex
         cout << "Found a track that has a different vertex id than the single vertex. Vertex id: "<< vertexID << "track vertex id: " << trk->getVertexId() << endl;
         continue;
      }

      fillTrackQualityCuts(trk);
      
      if(!goodQualityTrack(trk))  continue;
      
      fillTrackQualityCutsAfter(trk);

      tpcCounter += 1;      
      
      if( !trk->getFlag(StUPCTrack::kBemc) )   continue;

      hTrackQualityFlow->Fill(7);
      tracksBEMC.push_back(iTrk);
   }
   hNTracksBEMC->Fill( tracksBEMC.size() );
   hNTracksTpc->Fill(tpcCounter);

   if(tracksBEMC.size() != 2)   return;
   
   hAnalysisFlow->Fill(JPSIBEMC);
   const StUPCTrack* track1 = mUpcEvt->getTrack( tracksBEMC[0] );
   const StUPCTrack* track2 = mUpcEvt->getTrack( tracksBEMC[1] );

   if(!track1 || !track2)   return;

   SaveBemcInfo(track1, 0);
   SaveBemcInfo(track2, 1);
   SaveVertexInfo(vtx,0);

   fillEtaVtxPlotsBefore(track1, track2, vtx->getPosZ());

   if(abs(vtx->getPosZ()) > VERTEXZRANGE)  return;

   if(abs(track1->getEta()) > MAXETA || abs(track2->getEta()) > MAXETA) return;

   if(!IsGoodEtaTrack(track1, 0) || !IsGoodEtaTrack(track2, 0))  return;

   fillEtaVtxPlotsAfter(track1, track2, vtx->getPosZ());
   hAnalysisFlow->Fill(JPSIVTXZETA);


   // back to back
   if(!backToBack(track1, track2))   return;

   hAnalysisFlow->Fill(JPSIBACKTOBACK);

   //PID
   fillNSigmaPlots(track1);
   fillNSigmaPlots(track2);
   if(!chiSquarePID(track1,track2))   return;

   hAnalysisFlow->Fill(JPSIPID);


   // save info about tracks
   TLorentzVector electron1, electron2, state;
   track1->getLorentzVector(electron1, mUtil->mass(ELECTRON));
   track2->getLorentzVector(electron2, mUtil->mass(ELECTRON));
   state = electron1 + electron2;
   SaveStateInfo(state,track1->getCharge() + track2->getCharge(),0 );
   SaveChiSquareInfo(track1, track2);
   SaveTrackInfo(track1,electron1, 0 );
   SaveTrackInfo(track2,electron2, 1);

   int side = 0;
   if(analysisWithRPs){
      
      if(!exactly1RPTrack(side))   return;

      const StUPCRpsTrack *trackRP = mRpEvt->getTrack(0);
      if(!trackRP){
         cout << "No RP track found. Leaving event." << endl;
         return;
      }
      
      hAnalysisFlow->Fill(JPSI1RP);

      SaveRPinfo(trackRP,side);

      if(!fiducialVolume(trackRP, side)){   
         return;
      }

      hAnalysisFlow->Fill(JPSIRPFIDCUT);
      
      SaveMissingPtInfo(track1, track2, trackRP);

   }else{
      hAnalysisFlow->Fill(JPSI1RP);
      hAnalysisFlow->Fill(JPSIRPFIDCUT);
   }
   
   
   //Qtot
   double totQ = track1->getCharge() + track2->getCharge();
   hTotQ->Fill(totQ);

   
   if(totQ == 0){
      mRecTree->FillRecTree();
      hAnalysisFlow->Fill(JPSIQTOT);
   }else{
      mRecTree->FillBcgTree();
   }

   if(DEBUG){
      cout << "Finished AnaJPsi::Make()" << endl;
   }
}



void AnaJPsi::Init(){

   mUtil = new Util();


   if( DEBUG )
     cout<<"AnaJPsi::Init() called"<<endl;

   mOutFile->cd();
   if(runSysStudy){
      mRecTree = new RecTree(nameOfSysStudyTree, SysStudyTreeBits, true);
      loadCuts(runSysStudy);
      mOutFile->mkdir(nameOfSysStudyDir);
      mOutFile->cd(nameOfSysStudyDir);
   }else {
      mRecTree = new RecTree(nameOfAnaJPsiTree, AnaJPsiTreeBits, true);
      loadCuts(runSysStudy);
      mOutFile->mkdir(nameOfAnaJPsiDir);
      mOutFile->cd(nameOfAnaJPsiDir);
   }
   
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nJPSISelectionCuts-1, 1, nJPSISelectionCuts);
   for(int tb=1; tb<nJPSISelectionCuts; ++tb) {
     hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisJPSI(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
     TString label; label.Form("%d",triggerID[tb]);
     hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }

   hTrackQualityFlow = new TH1D("hTrackQualityFlow", "hTrackQualityFlow", 6,1,7);
   hTrackQualityFlow->GetXaxis()->SetBinLabel(1, TString("all"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(2, TString::Format("|DCA_{Z}| < %.1f cm", MAXDCAZ));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(3, TString::Format("DCA_{XY} < %.1f cm", MAXDCAXY));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(4, TString::Format("N^{fit}_{hits} > %d", MINNHITSFIT));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(5, TString::Format("N^{dE/dx}_{hits} > %d", MINNHITSDEDX));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(6, TString("BEMC match"));
   

   hEta = new TH1D("hEta", "Pseudorapidity; #eta_{e}; counts", 60, -2, 2);
   hEtaCut = new TH1D("hEtaCut", "Pseudorapidity; #eta_{e} [-]; counts", 60, -2, 2);

   hEtaBemc = new TH1D("hEtaBemc", "Pseudorapidity BEMC; #eta^{e}_{BEMC}; counts", 60, -2, 2);
   hEtaBemcCut = new TH1D("hEtaBemcCut", "Pseudorapidity BEMC; #eta^{e}_{BEMC} [-]; counts", 60, -2, 2);

   hEtaPhi = new TH2D("hEtaPhi","Phi vs eta of TPC tracks; #eta_{e}; #varphi",100,-2,2,100,-4,4);
   hEtaPhiCut = new TH2D("hEtaPhiCut","Phi vs eta of TPC tracks; #eta_{e} ; #varphi",100,-2,2,100,-4,4);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 
   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 

   hEtaVtxZ = new TH2D("hEtaVtxZ", "hEtaVtxZ; #eta_{e} [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);
   hEtaVtxZCut = new TH2D("hEtaVtxZCut", "hEtaVtxZCut; #eta_{e} [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);

   hNSigmaPiPcorr = new TH2F("hNSigmaPiPcorr","n_{#sigma} pions against protons;n#sigma_{#pi};n#sigma_{p}", 100, -25, 25, 100, -25, 25);
   hNSigmaPiKcorr = new TH2F("hNSigmaPiKcorr","n_{#sigma} pions against kaons;n#sigma_{#pi};n#sigma_{K}", 100, -25, 25, 100, -25, 25);
   hNSigmaPiecorr = new TH2F("hNSigmaPiecorr","n_{#sigma} pions against electrons;n#sigma_{#pi};n#sigma_{e}", 100, -25, 25, 100, -25, 25);
   hNSigmaPPicorr = new TH2F("hNSigmaPPicorr","n_{#sigma} protons against pions;n#sigma_{p};n#sigma_{#pi}", 100, -25, 25, 100, -25, 25);
   hNSigmaPKcorr = new TH2F("hNSigmaPKcorr","n_{#sigma} protons against kaons;n#sigma_{p};n#sigma_{K}", 100, -25, 25, 100, -25, 25);
   hNSigmaPecorr = new TH2F("hNSigmaPecorr","n_{#sigma} protons against electrons;n#sigma_{p};n#sigma_{e}", 100, -25, 25, 100, -25, 25);
   hNSigmaKPicorr = new TH2F("hNSigmaKPicorr","n_{#sigma} kaons against pions;n#sigma_{K};n#sigma_{#pi}", 100, -25, 25, 100, -25, 25);
   hNSigmaKPcorr = new TH2F("hNSigmaKPcorr","n_{#sigma} kaons against protons;n#sigma_{K};n#sigma_{p}", 100, -25, 25, 100, -25, 25);
   hNSigmaKecorr = new TH2F("hNSigmaKecorr","n_{#sigma} kaons against electrons;n#sigma_{K};n#sigma_{e}", 100, -25, 25, 100, -25, 25);
   hNSigmaPi = new TH1D("hNSigmaPi", "n_{#sigma} pions; n_{#sigma} #pi [-]; counts", 100, -10, 10 );
   hNSigmaP = new TH1D("hNSigmaP", "n_{#sigma} protons; n_{#sigma} p [-]; counts", 100, -10, 10 );
   hNSigmaK = new TH1D("hNSigmaK", "n_{#sigma} kaons; n_{#sigma} K [-]; counts", 100, -10, 10 );
   hDEdxSignal = new TH1D("hDEdxSignal", "dE/dx signal; ln#frac{dE}{dx} [MeV/cm];counts", 100, 0, 0.1 );

   hPIDChiee = new TH1D("hPIDChiee", "hPIDChiee; #chi_{ee} [-]; counts", 50, 0, 50);
   hPIDChipp = new TH1D("hPIDChipp", "hPIDChipp; #chi_{pp} [-]; counts", 50, 0, 50);
   hPIDChipipi = new TH1D("hPIDChipipi", "hPIDChipipi; #chi_{#pi #pi} [-]; counts", 50, 0, 50);
   hPIDChikk = new TH1D("hPIDChikk", "hPIDChikk; #chi_{kk} [-]; counts", 50, 0, 50);

   hNSigmaEE1 = new TH2F("hNSigmaEE1","n_{#sigma} e against e;n#sigma_{e} [-];n#sigma_{e} [-]", 40, -10, 10, 40, -10, 10);
   hNSigmaEE2 = new TH2F("hNSigmaEE2","n_{#sigma} e against e;n#sigma_{e} [-];n#sigma_{e} [-]", 40, -10, 10, 40, -10, 10);

   hNSigmaPP1 = new TH2F("hNSigmaPP1","n_{#sigma} p against p;n#sigma_{p} [-];n#sigma_{p} [-]", 40, -10, 10, 40, -10, 10);
   hNSigmaPP2 = new TH2F("hNSigmaPP2","n_{#sigma} p against p;n#sigma_{p} [-];n#sigma_{p} [-]", 40, -10, 10, 40, -10, 10);

   hNSigmaKK1 = new TH2F("hNSigmaKK1","n_{#sigma} K against K;n#sigma_{K} [-];n#sigma_{K} [-]", 40, -10, 10, 40, -10, 10);
   hNSigmaKK2 = new TH2F("hNSigmaKK2","n_{#sigma} K against K;n#sigma_{K} [-];n#sigma_{K} [-]", 40, -10, 10, 40, -10, 10);

   hNSigmaPiPi1 = new TH2F("hNSigmaPiPi1","n_{#sigma} #pi against #pi;n#sigma_{#pi} [-];n#sigma_{#pi} [-]", 40, -10, 10, 40, -10, 10);
   hNSigmaPiPi2 = new TH2F("hNSigmaPiPi2","n_{#sigma} #pi against #pi;n#sigma_{#pi} [-];n#sigma_{#pi} [-]", 40, -10, 10, 40, -10, 10);

   //initialize 2D graphs of before and after fiducial RP cut
   hRPcorr[0] = new TH2F("hRPcorr","p_{y} vs p_{x} of protons in RP;p_{x} [GeV];p_{y} [GeV]", 200, -0.7, 0.8, 150, -1, 1);

   hRPcorr[1] = new TH2F("hRPcorr_fid","p_{y} vs p_{x} of proton in RP fiducial region;p_{x} [GeV];p_{y} [GeV]", 200, -0.7, 0.8, 150, -1, 1);

   hRPcorrWest[0] = new TH2F("hRPcorrWest","p_{y} vs p_{x} of protons in RP (west);p_{x} [GeV];p_{y} [GeV]", 200, -0.7, 0.8, 150, -1, 1);

   hRPcorrWest[1] = new TH2F("hRPcorrWest_fid","p_{y} vs p_{x} of proton in RP fiducial region (west);p_{x} [GeV];p_{y} [GeV]", 200, -0.7, 0.8, 150, -1, 1);

   hRPcorrEast[0] = new TH2F("hRPcorrEast","p_{y} vs p_{x} of protons in RP (east);p_{x} [GeV];p_{y} [GeV]", 200, -0.7, 0.8, 150, -1, 1);

   hRPcorrEast[1] = new TH2F("hRPcorrEast_fid","p_{y} vs p_{x} of proton in RP fiducial region (east);p_{x} [GeV];p_{y} [GeV]",  200, -0.7, 0.8, 150, -1, 1);


   hNTracksRP = new TH1D("hNTracksRP", "Number of Tracks in RPs per event; n_{tracks}^{RP} [-]; counts", 30, -0.5, 29.5);

   hBranchRP = new TH1D("hBranchRP", "hBranchRP; detector branch; counts", 4,-0.5,3.5);
   //detectors branch, EU=0, ED=1, WU=2, WD=3
   hBranchRP->GetXaxis()->SetBinLabel(1, "East Up");
   hBranchRP->GetXaxis()->SetBinLabel(2, "East Down");
   hBranchRP->GetXaxis()->SetBinLabel(3, "West Up");
   hBranchRP->GetXaxis()->SetBinLabel(4, "West Down");

   hNTracksBEMC = new TH1D("hNTracksBEMC", "Number of Tracks in BEMC per event; Number of tracks in BEMC; counts", 21, -0.5, 20.5);

   hPt = new TH1D("hPt", "Transverse momentum of hadrons; p^{e}_{T} [GeV/c^{2}]; counts", 30, 0, 3);
   hPtCut = new TH1D("hPtCut", "hPtCut;p^{e}_{T} [GeV/c]; counts", 30, 0, 3);

   hDcaZ = new TH1D("hDcaZ", "hDcaZ; DCA_{Z} [cm]; counts", 40, -2, 2);
   hDcaZCut = new TH1D("hDcaZCut", "hDcaZCut; DCA_{Z} [cm]; counts", 40, -2, 2);

   hDcaXY = new TH1D("hDcaXY" , "hDcaXY; DCA_{XY}; counts",  40,0,4);
   hDcaXYCut = new TH1D("hDcaXYCut" , "hDcaXYCut; DCA_{XY}; counts", 40,0,4);

   hNfitHits = new TH1D("hNfitHits", "NfitHits; N^{fit}_{hits} [-]; counts", 50, 0, 50);
   hNfitHitsCut = new TH1D("hNfitHitsCut", "NfitHitsCut; N^{fit}_{hits} [-]; counts", 50, 0, 50);

   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);
   hNhitsDEdxCut = new TH1D("hNhitsDEdxCut", "NhitsDEdxCut; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);

   hNVertices = new TH1D("hNVertices", "Number of vertices before cut; n_{vertices} [-]; counts", 11, -0.5, 10.5);

   hTotQ = new TH1D("hTotQ", "Total charge of pair; Q_{tot} [-]; counts", 3, -1.5, 1.5);

   hNTracksTof = new TH1D("hNTracksTof", "hNTracksTof; Number of ToF tracks [-]; counts", 11, -0.5, 10.5);

   hNTracksTpc = new TH1D("hNTracksTpc", "hNTracksTpc; Number of TPC tracks [-]; counts", 26, -0.5, 25.5);

   hVtxZByFillNum = new TH1D("hVtxZByFillNum", "hVtxZByFillNum; V_{Z} [cm]; counts", 40, -200, 200);
  
   cout << "Finished AnaJPsi::Init()" << endl;
}

void AnaJPsi::loadCuts(bool isSysStudy){

   if(isSysStudy){
      VERTEXZRANGE = vertexRangeLoose;
      MINNHITSFIT = minNHitsFitLoose;
      MINNHITSDEDX = minNHitsDEdxLoose;
      MINDCAZ = minDcaZLoose;
      MAXDCAZ = maxDcaZLoose;
      MAXDCAXY = maxDcaXYLoose;
      MAXETA = maxEtaLoose;
      MINPIDCHIPP = minPidChiPPLoose;
      MINPIDCHIPIPI = minPidChiPiPiLoose;
      MINPIDCHIKK = minPidChiKKLoose;
      MAXPIDCHIEE = maxPidChiEELoose; // no cut
   }else{
      VERTEXZRANGE = vertexRange;
      MINNHITSFIT = minNHitsFit;
      MINNHITSDEDX = minNHitsDEdx;
      MINDCAZ = minDcaZ;
      MAXDCAZ = maxDcaZ;
      MAXDCAXY = maxDcaXY;
      MAXETA = maxEta;
      MINPIDCHIPP = minPidChiPP;
      MINPIDCHIPIPI = minPidChiPiPi;
      MINPIDCHIKK = minPidChiKK;
      MAXPIDCHIEE = maxPidChiEE; // no cut
   }
}


bool AnaJPsi::goodQualityTrack(const StUPCTrack *trk){


   hTrackQualityFlow->Fill(2);
   if( !(trk->getDcaZ() > MINDCAZ && trk->getDcaZ() < MAXDCAZ) )  //DCA z
      return false;
   hTrackQualityFlow->Fill(3);
   if( !(trk->getDcaXY() < MAXDCAXY) ) //DCA xy
      return false;
   hTrackQualityFlow->Fill(4);
   if( !(trk->getNhitsFit() > MINNHITSFIT) )  //NhitsFit
      return false;
   hTrackQualityFlow->Fill(5);
   if( !(trk->getNhitsDEdx() > MINNHITSDEDX) ) //NhitsdEdx
      return false;
   hTrackQualityFlow->Fill(6);

   return true;

}

bool AnaJPsi::sameVertex(const StUPCTrack *trk1,const StUPCTrack *trk2){
   int vtx1ID, vtx2ID;
   vtx1ID = trk1->getVertexId();
   vtx2ID = trk2->getVertexId();

   const StUPCVertex *vtx = mUpcEvt->getVertex(vtx1ID);

   if(!vtx)
      return false;

   if(vtx1ID != vtx2ID)
      return false;

   // additional condition on DCA XY
   if(trk1->getDcaXY() > maxDcaXY || trk2->getDcaXY() > maxDcaXY)
      return false;
   
   // additional condition on DCA Z
   if((trk1->getDcaZ() < minDcaZ || trk1->getDcaZ() > maxDcaZ) || (trk2->getDcaZ() < minDcaZ || trk2->getDcaZ() > maxDcaZ) )
      return false;

   return true;
}


void AnaJPsi::fillTrackQualityCuts(const StUPCTrack* trk){

   hEtaBemc->Fill( trk->getBemcEta() );
   hEta->Fill(trk->getEta() );
   hEtaPhi->Fill(trk->getEta(), trk->getPhi());
   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
   hDcaZ->Fill(trk->getDcaZ());
   hDcaXY->Fill(trk->getDcaXY());
}

void AnaJPsi::fillTrackQualityCutsAfter(const StUPCTrack* trk){

   hEtaBemcCut->Fill( trk->getBemcEta() );
   hEtaCut->Fill(trk->getEta() );
   hEtaPhiCut->Fill(trk->getEta(), trk->getPhi());
   hNfitHitsCut->Fill(trk->getNhitsFit() );
   hNhitsDEdxCut->Fill(trk->getNhitsDEdx() );
   hPtCut->Fill(trk->getPt());
   hDcaZCut->Fill(trk->getDcaZ());
   hDcaXYCut->Fill(trk->getDcaXY());

}

bool AnaJPsi::chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2){

   Double_t chi_pp = pow(trk1->getNSigmasTPCProton(),2) + pow(trk2->getNSigmasTPCProton(),2);
   Double_t chi_ee = pow(trk1->getNSigmasTPCElectron(),2) + pow(trk2->getNSigmasTPCElectron(),2);
   Double_t chi_kk = pow(trk1->getNSigmasTPCKaon(),2) + pow(trk2->getNSigmasTPCKaon(),2);
   Double_t chi_pipi = pow(trk1->getNSigmasTPCPion(),2) + pow(trk2->getNSigmasTPCPion(),2);
   
   hPIDChipipi->Fill(chi_pipi);
   hPIDChiee->Fill(chi_ee);
   hPIDChipp->Fill(chi_pp);
   hPIDChikk->Fill(chi_kk);
   
   hNSigmaEE1->Fill(trk1->getNSigmasTPCElectron(), trk2->getNSigmasTPCElectron());
   hNSigmaPP1->Fill(trk1->getNSigmasTPCProton(), trk2->getNSigmasTPCProton());
   hNSigmaKK1->Fill(trk1->getNSigmasTPCKaon(), trk2->getNSigmasTPCKaon());
   hNSigmaPiPi1->Fill(trk1->getNSigmasTPCPion(), trk2->getNSigmasTPCPion());

   if(chi_pp < MINPIDCHIPP)  return false;

   if(chi_pipi < MINPIDCHIPIPI)  return false;

   if(chi_kk < MINPIDCHIKK)  return false;

   hNSigmaEE2->Fill(trk1->getNSigmasTPCElectron(), trk2->getNSigmasTPCElectron());
   hNSigmaPP2->Fill(trk1->getNSigmasTPCProton(), trk2->getNSigmasTPCProton());
   hNSigmaKK2->Fill(trk1->getNSigmasTPCKaon(), trk2->getNSigmasTPCKaon());
   hNSigmaPiPi2->Fill(trk1->getNSigmasTPCPion(), trk2->getNSigmasTPCPion());

   if(chi_ee > MAXPIDCHIEE)  return false;

   return true;
}

void AnaJPsi::fillNSigmaPlots(const StUPCTrack *trk){
    Double_t nSigmaPion, nSigmaProton, nSigmaKaon, nSigmaElectron;
    nSigmaPion = trk->getNSigmasTPCPion();
    nSigmaProton = trk->getNSigmasTPCProton();
    nSigmaKaon = trk->getNSigmasTPCKaon();
    nSigmaElectron = trk->getNSigmasTPCElectron();

    hNSigmaPi->Fill(nSigmaPion);
    hNSigmaP->Fill(nSigmaProton);
    hNSigmaK->Fill(nSigmaKaon);
    hDEdxSignal->Fill( trk->getDEdxSignal() );
    hNSigmaPiPcorr->Fill(nSigmaPion, nSigmaProton);
    hNSigmaPiKcorr->Fill(nSigmaPion, nSigmaKaon);
    hNSigmaPiecorr->Fill(nSigmaPion, nSigmaElectron);
    hNSigmaPKcorr->Fill(nSigmaProton, nSigmaKaon);
    hNSigmaPecorr->Fill(nSigmaProton, nSigmaElectron);
    hNSigmaPPicorr->Fill(nSigmaProton, nSigmaPion);
    hNSigmaKecorr->Fill(nSigmaKaon, nSigmaElectron);
    hNSigmaKPcorr->Fill(nSigmaKaon, nSigmaProton);
    hNSigmaKPicorr->Fill(nSigmaKaon, nSigmaPion);

}


bool AnaJPsi::exactly1RPTrack(int &side){
   //1 RP track condition  
   AnaRpTracks(mRpEvt);

   unsigned int nRpTracksTotal = 0;
   side = 0;
   for (unsigned int iSide = 0; iSide < nSides; ++iSide)
   {
      for (unsigned int iTrck = 0; iTrck < mRpTrackIdVec_perSide[iSide].size(); ++iTrck)
      {
         StUPCRpsTrack* trkRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][iTrck]);
         if(!trkRP){
            cout << "Incorrect RP track. Leaving this event." << endl;
            return false;
         }
         nRpTracksTotal++;
         side = iSide;
      }
   }

   hNTracksRP->Fill(nRpTracksTotal);
   if(nRpTracksTotal == 1){
      return true;
   }else {
      return false;
   }

}

bool AnaJPsi::fiducialVolume(const StUPCRpsTrack* trackRP, int side){
   hRPcorr[0]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   if(side == East){
      hRPcorrEast[0]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   }else{
      hRPcorrWest[0]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   }

   //fiducial volume condition
   if( !RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y()) )
      return false;
   
   hRPcorr[1]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   if(side == East){
      hRPcorrEast[1]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   }else{
      hRPcorrWest[1]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
   }
   hBranchRP->Fill( trackRP->branch() );

   return true;

}


void AnaJPsi::fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ){

   Double_t eta2 = trk2->getEta();
   Double_t phi2 = trk2->getPhi();

   Double_t eta1 = trk1->getEta();
   Double_t phi1 = trk1->getPhi();

   hPosZ->Fill(posZ);

   hEtaVtxZ->Fill(eta1 , posZ );
   hEtaVtxZ->Fill(eta2 , posZ );

   hEtaPhi->Fill(eta1, phi1);
   hEtaPhi->Fill(eta2, phi2);

   hEta->Fill(eta1);
   hEta->Fill(eta2);

}

void AnaJPsi::fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ){

   Double_t eta2 = trk2->getEta();
   Double_t phi2 = trk2->getPhi();

   Double_t eta1 = trk1->getEta();
   Double_t phi1 = trk1->getPhi();

   hPosZCut->Fill(posZ);

   hEtaVtxZCut->Fill(eta1 , posZ );
   hEtaVtxZCut->Fill(eta2 , posZ );

   hEtaPhiCut->Fill(eta1, phi1);
   hEtaPhiCut->Fill(eta2, phi2);

   hEtaCut->Fill(eta1);
   hEtaCut->Fill(eta2);
}
