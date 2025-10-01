#include "../include/EmbeddingJPsi.h"

//_____________________________________________________________________________
EmbeddingJPsi::EmbeddingJPsi(TFile *outFile): Ana(outFile){}

EmbeddingJPsi::~EmbeddingJPsi(){
   //if(mUtil) delete mUtil;
}


void EmbeddingJPsi::Make(){

   //trigger
   hAnalysisFlow->Fill(EMBEDDINGALL);

   trueMCPeak();


   SaveEventInfo(mUpcEvt);
   // 1 vertex

   tpcCounter = 0;
   tracksBEMC.clear();
   hNClustersBEMC->Fill( mUpcEvt->getNumberOfClusters() );
   
   //central system good quality tracks + BEMC
   fillBemcInfoAll();
   for (int iTrk = 0; iTrk < mUpcEvt->getNumberOfTracks(); ++iTrk){
      const StUPCTrack *trk = mUpcEvt->getTrack(iTrk);
      
      if(!trk)  continue;

      ++tpcCounter;

      // vertexing does not work in embedding for pp, i have to use the truth information otherwise i would have to go over global tracks
      if(!(trk->getIdTruth() == 1 || trk->getIdTruth() == 2) )  continue;
      
      hTrackQualityFlow->Fill(1);
      
      if( trk->getFlag(StUPCTrack::kBemc) )  fillBemcInfo(trk);  // before, here was continue
      
      if(runCustomBemcSimulator){  // custom simulator for BEMC efficiency
         if( isBemcHit(trk) ){
            continue;
         }
      }else{  // standard way from embedding
         if( !trk->getFlag(StUPCTrack::kBemc) )  {
            continue;
         }
      }

      hTrackQualityFlow->Fill(2);

      fillTrackQualityCuts(trk);

      if(!goodQualityTrack(trk))  continue;
      
      fillTrackQualityCutsAfter(trk);
      
      tracksBEMC.push_back(iTrk);
   }
   
   hNTracksTpc->Fill( tpcCounter );
   
   hNTracksBEMC->Fill( tracksBEMC.size() );

   if(tracksBEMC.size() != 2)  return;
   
   hAnalysisFlow->Fill(EMBEDDING2BEMC);
   
   const StUPCTrack *track1 = mUpcEvt->getTrack(tracksBEMC[0]);

   const StUPCTrack *track2 = mUpcEvt->getTrack(tracksBEMC[1]);
   
   SaveBemcInfo(track1, 0);
   SaveBemcInfo(track2, 1);

   // eta condition
   if( !(abs(track1->getEta()) < MAXETA && abs(track2->getEta()) < MAXETA ) )   return;

   hAnalysisFlow->Fill(EMBEDDINGETA);

   if(!backToBack(track1, track2))  return;

   mRecTree->setIsBackToBack(1,0);

   hAnalysisFlow->Fill(EMBEDDINGBACKTOBACK);

   //PID
   fillNSigmaPlots(track1);
   fillNSigmaPlots(track2);

   if(!chiSquarePID(track1,track2))  return;
   
   hAnalysisFlow->Fill(EMBEDDINGPID);
   

   TLorentzVector electron1, electron2, state;   
   track1->getLorentzVector(electron1, mUtil->mass(ELECTRON));
   track2->getLorentzVector(electron2, mUtil->mass(ELECTRON));
   state = electron1 + electron2;

   hPtPair->Fill(state.Pt());

   SaveStateInfo(state,track1->getCharge() + track2->getCharge(),0 );
   SaveChiSquareInfo(track1, track2);
   SaveTrackInfo(track1,electron1, 0 );
   SaveTrackInfo(track2,electron2, 1);
   
   double totQ = track1->getCharge() + track2->getCharge();

   hTotQ->Fill(totQ);
   
   if(totQ == 0){
      mRecTree->FillRecTree();
      hAnalysisFlow->Fill(EMBEDDINGQTOT);
   }else{
      mRecTree->FillBcgTree();
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
     mRecTree = new RecTree(nameOfEmbeddingJPsiTree, EmbeddingJPsiTreeBits, true);
     mOutFile->mkdir(nameOfEmbeddingJPsiDir);
     mOutFile->cd(nameOfEmbeddingJPsiDir);
   loadCuts();  // if running sys study, will include loose cuts 
   
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nEmbeddingCuts-1, 1, nEmbeddingCuts);
   for(int tb=1; tb<nEmbeddingCuts; ++tb) {
     hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->embeddingName(tb));
   }


   hTrackQualityFlow = new TH1D("hTrackQualityFlow", "hTrackQualityFlow", 5,1,6);
   hTrackQualityFlow->GetXaxis()->SetBinLabel(1, TString("all"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(2, TString("BEMC hit"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(3, TString::Format("|#eta_{BEMC}| < %f", maxEta));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(4, TString::Format("N^{fit}_{hits} > %d", minNHitsFit));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(5, TString::Format("N^{dE/dx}_{hits} > %d", minNHitsDEdx));
   
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
   hDEdxSignal = new TH1D("hDEdxSignal", "dE/dx signal; ln#frac{dE}{dx} [MeV/cm];counts", 100, 0, 1 );

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


   hNTracksBEMC = new TH1D("hNTracksBEMC", "Number of Tracks in BEMC per event; Number of tracks in BEMC; counts", 21, -0.5, 20.5);

   hNClustersBEMC = new TH1D("hNClustersBEMC", "Number of Clusters in BEMC per event; Number of clusters in BEMC; counts", 21, -0.5, 20.5);

   hEClusters = new TH1D("hEClusters", "Energy of clusters in BEMC; E_{cluster} [GeV]; counts", 200, 0, 2);
   hEClusterTrack = new TH2D("hEClusterTrack", "Energy of clusters in BEMC; E_{cluster} [GeV]; E_{track}", 200, 0, 2, 200, 0, 2);
   hBemcTrackPhi = new TH1D("hBemcTrackPhi", "BEMC track phi; #varphi_{BEMC} [rad]; counts", 100, -4, 4);
   hBemcClusterPhi = new TH1D("hBemcClusterPhi", "BEMC cluster phi; #varphi_{BEMC} [rad]; counts", 100, -4, 4);
   hBemcEtaPhiTrack = new TH2D("hBemcEtaPhiTrack", "BEMC track eta phi; #eta_{BEMC}; #varphi_{BEMC}", 100, -1, 1, 100, -4, 4);
   hBemcEtaPhiCluster = new TH2D("hBemcEtaPhiCluster", "BEMC cluster eta phi; #eta_{BEMC}; #varphi_{BEMC}", 100, -1, 1, 100, -4, 4);
   hEClustersAll = new TH1D("hEClustersAll", "Energy of clusters in BEMC; E_{cluster} [GeV]; counts", 200, 0, 2);
   hBemcEtaPhiClusterAll = new TH2D("hBemcEtaPhiClusterAll", "BEMC cluster eta phi; #eta_{BEMC}; #varphi_{BEMC}", 100, -1, 1, 100, -4, 4);
   hBemcClusterPhiAll = new TH1D("hBemcClusterPhiAll", "BEMC cluster phi; #varphi_{BEMC} [rad]; counts", 100, -4, 4);
   hClusterMatched = new TH1D("hClusterMatched", "Cluster matched; Cluster matching (0 == all, 1 == matched with track); counts", 2, -0.5, 1.5);

   hPt = new TH1D("hPt", "Transverse momentum of hadrons", 120, 0, 3);
   hPt->SetTitle("Distribution of p_{T}");
   hPt->GetXaxis()->SetTitle("p_{T} [GeV]");
   hPt->GetYaxis()->SetTitle(YAxisDescription);

   hPtCut = new TH1D("hPtCut", "hPtCut;p_{T} [GeV/c]; counts", 120, 0, 3);

   hPtPair = new TH1D("hPtPair","transverse momentum of pair of electrons", 120, 0, 3);


   hNVertices = new TH1D("hNVertices", "Number of vertices before cut; n_{vertices} [-]; counts", 11, -0.5, 10.5);

   hDcaZ = new TH1D("hDcaZ", "hDcaZ; DCA_{Z} [cm]; counts", 40, -2, 2);
   hDcaZCut = new TH1D("hDcaZCut", "hDcaZCut; DCA_{Z} [cm]; counts", 40, -2, 2);

   hDcaXY = new TH1D("hDcaXY" , "hDcaXY; DCA_{XY}; counts",  40,0,4);
   hDcaXYCut = new TH1D("hDcaXYCut" , "hDcaXYCut; DCA_{XY}; counts", 40,0,4);

   hNfitHits = new TH1D("hNfitHits", "NfitHits; N^{fit}_{hits} [-]; counts", 50, 0, 50);
   hNfitHitsCut = new TH1D("hNfitHitsCut", "NfitHitsCut; N^{fit}_{hits} [-]; counts", 50, 0, 50);

   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);
   hNhitsDEdxCut = new TH1D("hNhitsDEdxCut", "NhitsDEdxCut; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);

   hTotQ = new TH1D("hTotQ", "Total charge of pair", 5, -2.5, 2.5);
   hTotQ->SetTitle("Total charge of pair");
   hTotQ->GetXaxis()->SetTitle("Charge of pair");
   hTotQ->GetYaxis()->SetTitle(YAxisDescription);

   hNTracksTpc = new TH1D("hNTracksTpc", "hNTracksTpc; Number of TPC tracks [-]; counts", 26, -0.5, 25.5);
  
   hInvMassJPsi = new TH1D("hInvMassTrueMC", "hInvMassJPsi; m_{ee} [GeV/c^{2}]; counts", 225, 0.5, 5);

   gBEControl = new TGraph(); // graph for control of BEMC efficiency
   gBEControl->SetName("gBEControl");
   gBEControl->SetTitle("BEMC efficiency used in MC simulation; p_{T} [GeV/c]; efficiency");
   runControlOfBemcEfficiency();
  
   bField = -4.991; // guess based on real runs
   beamline[0] = 0;
   beamline[2] = 0;
   beamline[1] = 0;
   beamline[3] = 0;

   cout << "Finished EmbeddingJPsi::Init()" << endl;
}


void EmbeddingJPsi::loadCuts(){

   if(runSysStudy){
      MINNHITSFIT = minNHitsFitLoose;
      MINNHITSDEDX = minNHitsDEdxLoose;
      MAXETA = maxEtaLoose;
      MINPIDCHIPP = minPidChiPP;
      MINPIDCHIPIPI = minPidChiPiPi;
      MINPIDCHIKK = minPidChiKK;
      MAXPIDCHIEE = maxPidChiEELoose; // no cut
   }else{
      MINNHITSFIT = minNHitsFit;
      MINNHITSDEDX = minNHitsDEdx;
      MAXETA = maxEta;
      MINPIDCHIPP = minPidChiPP;
      MINPIDCHIPIPI = minPidChiPiPi;
      MINPIDCHIKK = minPidChiKK;
      MAXPIDCHIEE = maxPidChiEE; // no cut
   }
}

void EmbeddingJPsi::trueMCPeak(){


   TParticle *mc1, *mc2;
   for(int iMcTrk = 0;iMcTrk < mUpcEvt->getNumberOfMCParticles(); ++iMcTrk){
      TParticle *mcTrk = mUpcEvt->getMCParticle(iMcTrk);
      if(!mcTrk)
         continue;
      if(mcTrk->GetPdgCode() == 11 ){
         mc1 = mcTrk;
      }
      else if(mcTrk->GetPdgCode() == -11){
         mc2 = mcTrk;
      }

   }

   if(!mc1 || !mc2){
      return;
   }

   TLorentzVector stateMC(mc1->Px() + mc2->Px(), mc1->Py() + mc2->Py(), mc1->Pz() + mc2->Pz(), mc1->Energy() + mc2->Energy());

   hInvMassJPsi->Fill(stateMC.M());


}

void EmbeddingJPsi::fillBemcInfo(const StUPCTrack *trk){

   if(!trk){
      return;
   }

   StUPCBemcCluster *cluster = trk->getBemcCluster();
   hEClusters->Fill(cluster->getEnergy());
   hEClusterTrack->Fill(cluster->getEnergy(), trk->getBemcHitE());
   hBemcTrackPhi->Fill(trk->getBemcPhi());
   hBemcClusterPhi->Fill(cluster->getPhi());
   hBemcEtaPhiTrack->Fill(trk->getBemcEta(), trk->getBemcPhi());
   hBemcEtaPhiCluster->Fill(cluster->getEta(), cluster->getPhi());
   hClusterMatched->Fill(1);

}

void EmbeddingJPsi::fillBemcInfoAll(){

   for(int i = 0; i < mUpcEvt->getNumberOfClusters(); ++i){
      StUPCBemcCluster *cluster = mUpcEvt->getCluster(i);
      if(!cluster)
         continue;
      hEClustersAll->Fill(cluster->getEnergy());
      hBemcClusterPhiAll->Fill(cluster->getPhi());
      hBemcEtaPhiClusterAll->Fill(cluster->getEta(), cluster->getPhi());
      hClusterMatched->Fill(0);
   }

}


bool EmbeddingJPsi::goodQualityTrack(const StUPCTrack *trk){


   if( !(abs(trk->getEta()) < MAXETA))  // eta
      return false;
   hTrackQualityFlow->Fill(3);   
   if( !(trk->getNhitsFit() > MINNHITSFIT) )  //NhitsFit
      return false;
   hTrackQualityFlow->Fill(4);
   if( !(trk->getNhitsDEdx() > MINNHITSDEDX) ) //NhitsdEdx
      return false;
   hTrackQualityFlow->Fill(5);

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
   //if(trk1->getDcaXY() > maxDcaXY || trk2->getDcaXY() > maxDcaXY)
   //   return false;
   
   // additional condition on DCA Z
   if((trk1->getDcaZ() < minDcaZ || trk1->getDcaZ() > maxDcaZ) || (trk2->getDcaZ() < minDcaZ || trk2->getDcaZ() > maxDcaZ) )
      return false;

   return true;
}


void EmbeddingJPsi::fillTrackQualityCuts(const StUPCTrack* trk){

   if(trk->getFlag(StUPCTrack::kBemc))
      hEtaBemc->Fill( trk->getBemcEta() );
   hEta->Fill(trk->getEta() );
   hEtaPhi->Fill(trk->getEta(), trk->getPhi());
   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
   hDcaZ->Fill(trk->getDcaZ());
   hDcaXY->Fill(trk->getDcaXY());
}

void EmbeddingJPsi::fillTrackQualityCutsAfter(const StUPCTrack* trk){

   if(trk->getFlag(StUPCTrack::kBemc))
      hEtaBemcCut->Fill( trk->getBemcEta() );
   hEtaCut->Fill(trk->getEta() );
   hEtaPhiCut->Fill(trk->getEta(), trk->getPhi());
   hNfitHitsCut->Fill(trk->getNhitsFit() );
   hNhitsDEdxCut->Fill(trk->getNhitsDEdx() );
   hPtCut->Fill(trk->getPt());
   hDcaZCut->Fill(trk->getDcaZ());
   hDcaXYCut->Fill(trk->getDcaXY());

}
bool EmbeddingJPsi::chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2){

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

void EmbeddingJPsi::fillNSigmaPlots(const StUPCTrack *trk){
    Double_t nSigmaPion, nSigmaProton, nSigmaKaon, nSigmaElectron;
    nSigmaPion = trk->getNSigmasTPCPion();
    nSigmaProton = trk->getNSigmasTPCProton();
    nSigmaKaon = trk->getNSigmasTPCKaon();
    nSigmaElectron = trk->getNSigmasTPCElectron();

    hNSigmaPi->Fill(nSigmaPion);
    hNSigmaP->Fill(nSigmaProton);
    hNSigmaK->Fill(nSigmaKaon);
    hDEdxSignal->Fill( trk->getDEdxSignal()*1000000 );
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


bool EmbeddingJPsi::isBemcHit(const StUPCTrack *trk){ // MC simulator of BEMC efficiency

   TRandom3 randGen;
   double random = randGen.Uniform(0,1);

   double val = bemcEfficiency( trk->getPt() );

   if(random < val) {
      return true;
   }else{
      return false;   
   }
   
}


double EmbeddingJPsi::bemcEfficiency(double pT){ //function which returns bemc efficiency based on fit to real data
   // for now, the results of the fit are here hard coded
   double eps0 = -1.303;
   double n = 1.0758;
   double pTThr = 0.0105;
   double sigma = 0.524;


   double res = eps0 + n*( 1 + TMath::Erf( (pT - pTThr)/(sqrt(2)*sigma ) ) );

   return res;
}


void EmbeddingJPsi::runControlOfBemcEfficiency(){

   int j = 0;
   for(int i = 50; i < 500; ++i){
      gBEControl->SetPoint(j , i*0.01, bemcEfficiency(i*0.01) );
      j++;
   }

   mOutFile->cd();
   gBEControl->Write();

}