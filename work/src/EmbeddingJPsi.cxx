#include "EmbeddingJPsi.h"

//_____________________________________________________________________________
EmbeddingJPsi::EmbeddingJPsi(TFile *outFile): Ana(outFile){}

EmbeddingJPsi::~EmbeddingJPsi(){
   //if(mUtil) delete mUtil;
}


void EmbeddingJPsi::Make(){

   //trigger
   hAnalysisFlow->Fill(JPSIALL);

   hAnalysisFlow->Fill(JPSITRIG);

   SaveEventInfo(mUpcEvt);
   // 1 vertex
   hNVertices->Fill(mUpcEvt->getNumberOfVertices());
   if(mUpcEvt->getNumberOfVertices() != 1){
      return;
   }
   hAnalysisFlow->Fill(JPSI1VTX);

   const StUPCVertex *vtx = mUpcEvt->getVertex(0);

   if(!vtx){
      cout << "Could not load vertex. Leaving event." << endl;
      return;
   }
   hPosZ->Fill(vtx->getPosZ());
   if(abs(vtx->getPosZ()) > vertexRange)
      return;
   hPosZCut->Fill(vtx->getPosZ());


   hAnalysisFlow->Fill(JPSIVTXZ);
   unsigned int vertexID = vtx->getId();

   tpcCounter = 0;
   tracksBEMC.clear();
   hNClustersBEMC->Fill( mUpcEvt->getNumberOfClusters() );

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

      hTrackQualityFlow->Fill(2);
      if( !trk->getFlag(StUPCTrack::kBemc) )
         continue;
      if( !trk->getFlag(StUPCTrack::kTof) )
         continue;

      fillTrackQualityCuts(trk);
      hEtaVtxZ->Fill(trk->getEta(), vtx->getPosZ());

      if(!goodQualityTrack(trk))
         continue;
      
      fillTrackQualityCutsAfter(trk);
      hEtaVtxZCut->Fill(trk->getEta(), vtx->getPosZ());

      hNTracksTpc->Fill( tpcCounter );

      tracksBEMC.push_back(iTrk);
   }

   hNTracksBEMC->Fill( tracksBEMC.size() );

   if(tracksBEMC.size() != 2){
      cout << "Found " << tracksBEMC.size() << " tracks in BEMC. Leaving event." << endl;
      return;
   }
   hAnalysisFlow->Fill(JPSIBEMC);
   const StUPCTrack* track1 = mUpcEvt->getTrack( tracksBEMC[0] );
   const StUPCTrack* track2 = mUpcEvt->getTrack( tracksBEMC[1] );

   if(!track1 || !track2)
      return;

   //back to back condition
   int sectionTrk1, sectionTrk2;
   for (int i = 0; i < 6; ++i){
      double lower, upper;
      lower = -TMath::Pi() + i*TMath::Pi()/3;
      upper = -TMath::Pi() + (i+1)*TMath::Pi()/3;

      if(track1->getBemcPhi() >= lower && track1->getBemcPhi() < upper){
         sectionTrk1 = i;
      }
      if(track2->getBemcPhi() >= lower && track2->getBemcPhi() < upper){
         sectionTrk2 = i;
      }
   }
   double phiDelta = abs(track1->getPhi() - track2->getPhi());
   if(phiDelta <= 2.6)
      return;
   int deltaSections = abs(sectionTrk1 - sectionTrk2);
   if(deltaSections != 3)
      return;
   hAnalysisFlow->Fill(JPSIBACKTOBACK);


   //PID
   fillNSigmaPlots(track1);
   fillNSigmaPlots(track2);
   if(!chiSquarePID(track1,track2))
      return;

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

   // just fill the ana flows for RP cuts even though they are not used
   hAnalysisFlow->Fill(JPSI1RP);
   hAnalysisFlow->Fill(JPSIRPFIDCUT);

   //Qtot
   double totQ = track1->getCharge() + track2->getCharge();
   //cout << "Charge: " << totQ << endl;
   hTotQ->Fill(totQ);


   if(totQ == 0){
      hAnalysisFlow->Fill(JPSIQTOT);
      mRecTree->FillRecTree();
      hInvMassJPsi->Fill(state.M());
   }else{
      mRecTree->FillBcgTree();
      hInvMassJPsiBcg->Fill(state.M());
   }


   if(DEBUG){
      cout << "Finished EmbeddingJPsi::Make()" << endl;
   }
}



void EmbeddingJPsi::Init(){

   mUtil = new Util();


   if( DEBUG )
     cout<<"EmbeddingJPsi::Init() called"<<endl;

   mOutFile->cd();
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nJPSISelectionCuts-1, 1, nJPSISelectionCuts);
   for(int tb=1; tb<nJPSISelectionCuts; ++tb) {
     hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisJPSI(tb));
   }

   hTrackQualityFlow = new TH1D("hTrackQualityFlow", "hTrackQualityFlow", 8,1,9);
   hTrackQualityFlow->GetXaxis()->SetBinLabel(1, TString("all"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(2, TString("primary vertex"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(3, TString("BEMC hit"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(4, TString::Format("|#eta_{BEMC}| < %.1f", maxEta));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(5, TString::Format("|DCA_{Z}| < %.1f cm", maxDcaZ));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(6, TString::Format("DCA_{XY} < %.1f cm",maxDcaXY ));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(7, TString::Format("N^{fit}_{hits} > %d", minNHitsFit));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(8, TString::Format("N^{dE/dx}_{hits} > %d", minNHitsDEdx));
   
   mRecTree = new RecTree(nameOfAnaJPsiTree, AnaJPsiTreeBits, true); 

   hEta = new TH1D("hEta", "Pseudorapidity; #eta; counts", 60, -2, 2);
   hEtaCut = new TH1D("hEtaCut", "Pseudorapidity; #eta [-]; counts", 60, -2, 2);

   hEtaBemc = new TH1D("hEtaBemc", "Pseudorapidity BEMC; #eta_{BEMC}; counts", 60, -2, 2);
   hEtaBemcCut = new TH1D("hEtaBemcCut", "Pseudorapidity BEMC; #eta_{BEMC} [-]; counts", 60, -2, 2);


   hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks; #eta; #varphi",100,-2,2,100,-4,4);
   hEtaPhiCut = new TH2F("hEtaPhiCut","Phi vs eta of TOF tracks; #eta ; #varphi",100,-2,2,100,-4,4);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 
   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 

   hEtaVtxZ = new TH2F("hEtaVtxZ", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);
   hEtaVtxZCut = new TH2F("hEtaVtxZCut", "hEtaVtxZCut; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);

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

   hPIDChiep = new TH2F("hPIDChiep", "hPIDChiep; #chi_{ee} [-]; #chi_{pp}", 50, 0, 50, 50, 0, 50);
   hPIDChiek = new TH2F("hPIDChiek", "hPIDChiek; #chi_{ee} [-]; #chi_{kk}", 50, 0, 50, 50, 0, 50);
   hPIDChiepi = new TH2F("hPIDChiepi", "hPIDChiepi; #chi_{ee} [-]; #chi_{#pi #pi}", 50, 0, 50, 50, 0, 50);
   hPIDChipip = new TH2F("hPIDChipip", "hPIDChipip; #chi_{#pi #pi} [-]; #chi_{pp}", 50, 0, 50, 50, 0, 50);


   //initialize 2D graphs of before and after fiducial RP cut
   hRPcorr[0] = new TH2F("hRPcorr","p_{y} vs p_{x} of protons in RP", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorr[0]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorr[0]->GetYaxis()->SetTitle("p_{y} [GeV]");
   hRPcorr[1] = new TH2F("hRPcorr_fid","p_{y} vs p_{x} of proton in RP fiducial region", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorr[1]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorr[1]->GetYaxis()->SetTitle("p_{y} [GeV]");

   hRPcorrWest[0] = new TH2F("hRPcorrWest","p_{y} vs p_{x} of protons in RP (west)", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorrWest[0]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrWest[0]->GetYaxis()->SetTitle("p_{y} [GeV]");
   hRPcorrWest[1] = new TH2F("hRPcorrWest_fid","p_{y} vs p_{x} of proton in RP fiducial region (west)", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorrWest[1]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrWest[1]->GetYaxis()->SetTitle("p_{y} [GeV]");

   hRPcorrEast[0] = new TH2F("hRPcorrEast","p_{y} vs p_{x} of protons in RP (east)", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorrEast[0]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrEast[0]->GetYaxis()->SetTitle("p_{y} [GeV]");
   hRPcorrEast[1] = new TH2F("hRPcorrEast_fid","p_{y} vs p_{x} of proton in RP fiducial region (east)",  200, -0.7, 0.8, 150, -1, 1);
   hRPcorrEast[1]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrEast[1]->GetYaxis()->SetTitle("p_{y} [GeV]");

   hNTracksRP = new TH1D("hNTracksRP", "Number of Tracks in RPs per event", 400, 0, 40);
   hNTracksRP->GetXaxis()->SetTitle("Number of tracks in RPs");
   hNTracksRP->GetYaxis()->SetTitle(YAxisDescription);
   hNTracksRP->SetTitle("Number of tracks in RPs per event");

   hBranchRP = new TH1D("hBranchRP", "hBranchRP; detector branch; counts", 4,-0.5,3.5);
   //detectors branch, EU=0, ED=1, WU=2, WD=3
   hBranchRP->GetXaxis()->SetBinLabel(1, "East Up");
   hBranchRP->GetXaxis()->SetBinLabel(2, "East Down");
   hBranchRP->GetXaxis()->SetBinLabel(3, "West Up");
   hBranchRP->GetXaxis()->SetBinLabel(4, "West Down");

   hNTracksBEMC = new TH1D("hNTracksBEMC", "Number of Tracks in BEMC per event; Number of tracks in BEMC; counts", 21, -0.5, 20.5);

   hNClustersBEMC = new TH1D("hNClustersBEMC", "Number of Clusters in BEMC per event; Number of clusters in BEMC; counts", 21, -0.5, 20.5);


   hPt = new TH1D("hPt", "Transverse momentum of hadrons", 30, 0, 3);
   hPt->SetTitle("Distribution of p_{T}");
   hPt->GetXaxis()->SetTitle("p_{T} [GeV]");
   hPt->GetYaxis()->SetTitle(YAxisDescription);

   hPtCut = new TH1D("hPtCut", "hPtCut;p_{T} [GeV/c]; counts", 30, 0, 3);


   hDcaZ = new TH1D("hDcaZ", "hDcaZ; DCA_{Z} [cm]; counts", 40, -2, 2);
   hDcaZCut = new TH1D("hDcaZCut", "hDcaZCut; DCA_{Z} [cm]; counts", 40, -2, 2);

   hDcaXY = new TH1D("hDcaXY" , "hDcaXY; DCA_{XY}; counts",  40,0,4);
   hDcaXYCut = new TH1D("hDcaXYCut" , "hDcaXYCut; DCA_{XY}; counts", 40,0,4);

   hNfitHits = new TH1D("hNfitHits", "NfitHits; N^{fit}_{hits} [-]; counts", 50, 0, 50);
   hNfitHitsCut = new TH1D("hNfitHitsCut", "NfitHitsCut; N^{fit}_{hits} [-]; counts", 50, 0, 50);

   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);
   hNhitsDEdxCut = new TH1D("hNhitsDEdxCut", "NhitsDEdxCut; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);

   hNVertices = new TH1D("hNVertices", "Number of vertices before cut", 11, -0.5, 10.5);
   hNVertices->SetTitle("Number of vertices");
   hNVertices->GetXaxis()->SetTitle("Number of vertices");
   hNVertices->GetYaxis()->SetTitle(YAxisDescription);

   hTotQ = new TH1D("hTotQ", "Total charge of pair", 5, -2.5, 2.5);
   hTotQ->SetTitle("Total charge of pair");
   hTotQ->GetXaxis()->SetTitle("Charge of pair");
   hTotQ->GetYaxis()->SetTitle(YAxisDescription);

   hNTracksTof = new TH1D("hNTracksTof", "hNTracksTof; Number of ToF tracks [-]; counts", 11, -0.5, 10.5);

   hNTracksTpc = new TH1D("hNTracksTpc", "hNTracksTpc; Number of TPC tracks [-]; counts", 26, -0.5, 25.5);
  
   hInvMassJPsi = new TH1D("hInvMassJPsi", "hInvMassJPsi; m_{e e} [GeV/c^{2}]; counts", 40,2, 4);
   hInvMassJPsiBcg = new TH1D("hInvMassJPsiBcg", "hInvMassJPsiBcg; m_{e e} [GeV/c^{2}]; counts", 40,2, 4);


   cout << "Finished EmbeddingJPsi::Init()" << endl;
}

bool EmbeddingJPsi::goodQualityTrack(const StUPCTrack *trk){

   
   if(trk->getPt() < minPt)  //pT
      return false;

   hTrackQualityFlow->Fill(3);
   if(!(abs(trk->getBemcEta()) < maxEta))  // eta
      return false;
   hTrackQualityFlow->Fill(4);
   if( !(trk->getDcaZ() > minDcaZ && trk->getDcaZ() < maxDcaZ) )  //DCA z
      return false;
   hTrackQualityFlow->Fill(5);   
   if( !(trk->getDcaXY() < maxDcaXY) ) //DCA xy
      return false;
   hTrackQualityFlow->Fill(6);
   if( !(trk->getNhitsFit() > minNHitsFit) )  //NhitsFit
      return false;
   hTrackQualityFlow->Fill(7);
   //if( !(trk->getNhitsDEdx() > minNHitsDEdx) ) //NhitsdEdx
   //   return false;
   hTrackQualityFlow->Fill(8);
   tpcCounter += 1;

   return true;

}

bool EmbeddingJPsi::sameVertex(const StUPCTrack *trk1,const StUPCTrack *trk2){
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


void EmbeddingJPsi::fillTrackQualityCuts(const StUPCTrack* trk){

   //hEtaBemc->Fill( trk->getBemcEta() );
   hEta->Fill(trk->getEta() );
   hEtaPhi->Fill(trk->getEta(), trk->getPhi());
   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
   hDcaZ->Fill(trk->getDcaZ());
   hDcaXY->Fill(trk->getDcaXY());
}

void EmbeddingJPsi::fillTrackQualityCutsAfter(const StUPCTrack* trk){

   //hEtaBemcCut->Fill( trk->getBemcEta() );
   hEtaCut->Fill(trk->getEta() );
   hEtaPhiCut->Fill(trk->getEta(), trk->getPhi());
   hNfitHitsCut->Fill(trk->getNhitsFit() );
   hNhitsDEdxCut->Fill(trk->getNhitsDEdx() );
   hPtCut->Fill(trk->getPt());
   hDcaZCut->Fill(trk->getDcaZ());
   hDcaXYCut->Fill(trk->getDcaXY());

}
bool EmbeddingJPsi::chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2){

   Double_t chi_ee = pow(trk1->getNSigmasTPCElectron(),2) + pow(trk2->getNSigmasTPCElectron(),2);
   Double_t chi_pp = pow(trk1->getNSigmasTPCProton(),2) + pow(trk2->getNSigmasTPCProton(),2);
   Double_t chi_kk = pow(trk1->getNSigmasTPCKaon(),2) + pow(trk2->getNSigmasTPCKaon(),2);
   Double_t chi_pipi = pow(trk1->getNSigmasTPCPion(),2) + pow(trk2->getNSigmasTPCPion(),2);

   hPIDChiee->Fill(chi_ee);
   hPIDChipp->Fill(chi_pp);
   hPIDChipipi->Fill(chi_pipi);
   hPIDChikk->Fill(chi_kk);
   hPIDChiep->Fill(chi_ee, chi_pp);
   hPIDChiek->Fill(chi_ee, chi_kk);
   hPIDChiepi->Fill(chi_ee, chi_pipi);
   hPIDChipip->Fill(chi_pipi, chi_pp);

   if(chi_ee > maxPIDChiEE)
      return false;
   if(chi_pp < minPIDChiPP)
      return false;
   if(chi_pipi < minPIDChiPiPi)
      return false;
   if(chi_kk < minPIDChiKK)
      return false;

   return true;
}

void EmbeddingJPsi::fillNSigmaPlots(const StUPCTrack *trk){
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
