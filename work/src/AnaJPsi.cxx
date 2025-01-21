#include "AnaJPsi.h"
#include "RunDef.h"

//_____________________________________________________________________________
AnaJPsi::AnaJPsi(TFile *outFile): Ana(outFile){}

AnaJPsi::~AnaJPsi(){
   //if(mUtil) delete mUtil;
}


void AnaJPsi::Make(){

   //trigger
   hAnalysisFlow->Fill(JPSIALL);
   if(!CheckTriggers(&triggerID, mUpcEvt, hTriggerBits))
      return;
   hAnalysisFlow->Fill(JPSITRIG);

   SaveEventInfo(mUpcEvt);
   SaveTriggerInfo(mUpcEvt, mRpEvt);

   tpcCounter = 0;
   //central system good quality tracks + BEMC
   for (int iTrk = 0; iTrk < mUpcEvt->getNumberOfTracks(); ++iTrk){
      const StUPCTrack *trk = mUpcEvt->getTrack(iTrk);
      if(!trk)
         continue;
      if(!trk->getFlag(StUPCTrack::kPrimary))
         continue;

      fillTrackQualityCuts(trk);
      if(!goodQualityTrack(trk))
         continue;
      fillTrackQualityCutsAfter(trk);
   hNTracksTpc->Fill( tpcCounter );

      if( !trk->getFlag(StUPCTrack::kBemc) )
         continue;
      hTrackQualityFlow->Fill(5);

      tracksBEMC.push_back(iTrk);
   }
   hAnalysisFlow->Fill(JPSIBEMC);
   hNTracksBEMC->Fill( tracksBEMC.size() );

   if(tracksBEMC.size() == 0){
      return;
   }

   //vertex + parovanie(back to back condition)
   // outer track loop
   const StUPCTrack *track1, *track2;
   for (unsigned int iTrk = 0; iTrk < tracksBEMC.size() ; ++iTrk){
      const StUPCTrack* trk1 = mUpcEvt->getTrack( tracksBEMC[iTrk] );
      if(!trk1)
         continue;

      // inner track loop
      for (unsigned int jTrk = iTrk; jTrk < tracksBEMC.size(); ++jTrk){
         if(tracksBEMC[iTrk] == tracksBEMC[jTrk])
            continue;
         const StUPCTrack* trk2 = mUpcEvt->getTrack( tracksBEMC[jTrk] );
         if(!trk2)
            continue;
         if(!sameVertex(trk1, trk2))
            continue;
         hAnalysisFlow->Fill(JPSI1VTX);


         //VtxZ - eta cut
         const StUPCVertex *vtx = mUpcEvt->getVertex(trk1->getVertexId());
         if(!vtx)
            continue;

         fillEtaVtxPlotsBefore(trk1, trk2, vtx->getPosZ());
         if(abs(vtx->getPosZ()) > vertexRangeForEVz )
            continue;
         SaveVertexInfo(vtx, 0);

         if( !( IsGoodEtaTrack(trk1, 0) && IsGoodEtaTrack(trk2, 0 ) ) )
            continue;

         fillEtaVtxPlotsAfter(trk1, trk2, vtx->getPosZ());
         hAnalysisFlow->Fill(JPSIETAVTXZ);


         //PID
         fillNSigmaPlots(trk1);
         fillNSigmaPlots(trk2);
         if(!chiSquarePID(trk1,trk2))
            continue;
         hAnalysisFlow->Fill(JPSIPID);

         // back to back
         Double_t deltaPhi = trk1->getBemcPhi() - trk2->getBemcPhi();
         if(deltaPhi > minBEMCPhi && deltaPhi < maxBEMCPhi ){
            track1 = trk1;
            track2 = trk2;
            hAnalysisFlow->Fill(JPSIBACKTOBACK);
            break;
         }else continue;
      }
      if(track1 && track2)
         break;
   }
   if(!track1 || !track2)
      return;

   // save info about tracks
   TLorentzVector hadron1, hadron2, state;
   track1->getLorentzVector(hadron1, mUtil->mass(ELECTRON));
   track2->getLorentzVector(hadron2, mUtil->mass(ELECTRON));

   SaveStateInfo(state,track1->getCharge() + track2->getCharge(),0 );
   SaveTrackInfo(track1,hadron1, 0 );
   SaveTrackInfo(track2,hadron2, 1);


   //1 RP track condition
   AnaRpTracks(mRpEvt);

   unsigned int nRpTracksTotal = 0;
   StUPCRpsTrack *trackRP;
   int side = 0;
   for (unsigned int iSide = 0; iSide < nSides; ++iSide)
   {
      for (unsigned int iTrck = 0; iTrck < mRpTrackIdVec_perSide[iSide].size(); ++iTrck)
      {
         StUPCRpsTrack* trkRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][iTrck]);
         if(!trkRP){
            cout << "Incorrect RP track. Leaving this event." << endl;
            return;
         }
         nRpTracksTotal++;
         trackRP = trkRP;
         side = iSide;
      }
   }
   if(nRpTracksTotal != 1)
      return;
   
   hAnalysisFlow->Fill(JPSI1RP);
   SaveRPinfo(trackRP,side);
   

   //fiducial volume condition
   if( !RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y()) )
      return;
   hAnalysisFlow->Fill(JPSIRPFIDCUT);


   //Qtot
   double totQ = track1->getCharge() + track2->getCharge();
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
      cout << "Finished AnaJPsi::Make()" << endl;
   }
}



void AnaJPsi::Init(){

   mUtil = new Util();


   if( DEBUG )
     cout<<"AnaJPsi::Init() called"<<endl;

   mOutFile->cd();
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nJPSISelectionCuts-1, 1, nJPSISelectionCuts);
   for(int tb=1; tb<nJPSISelectionCuts; ++tb) {
     hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisJPSI(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
     TString label; label.Form("%d",triggerID[tb]);
     hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }

   hTrackQualityFlow = new TH1D("hTrackQualityFlow", "hTrackQualityFlow", 5,1,6);
   hTrackQualityFlow->GetXaxis()->SetBinLabel(1, TString("all"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(2, TString("p_{T} > 0.2 GeV/c"));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(3, TString::Format("N^{fit}_{hits} > %d", minNHitsFit));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(4, TString::Format("N^{dE/dx}_{hits} > %d", minNHitsDEdx));
   hTrackQualityFlow->GetXaxis()->SetBinLabel(5, TString("BEMC"));

   mRecTree = new RecTree(nameOfAnaJPsiTree, AnaJPsiTreeBits, true); 

   hEta = new TH1D("hEta", "Pseudorapidity; #eta; counts", 60, -2, 2);
   hEtaCut = new TH1D("hEtaCut", "Pseudorapidity; #eta [-]; counts", 60, -2, 2);

   hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks; #eta; #varphi",100,-2,2,100,-4,4);
   hEtaPhiCut = new TH2F("hEtaPhiCut","Phi vs eta of TOF tracks; #eta ; #varphi",100,-2,2,100,-4,4);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 
   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 

   hEtaVtxZ = new TH2F("hEtaVtxZ", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);
   hEtaVtxZCut = new TH2F("hEtaVtxZCut", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);

   hNSigmaPiPcorr = new TH2F("hNSigmaPiPcorr","n_{#sigma} pions against protons;n#sigma_{#pi};n#sigma_{p}", 100, -25, 25, 100, -25, 25);

   hNSigmaPiKcorr = new TH2F("hNSigmaPiKcorr","n_{#sigma} pions against kaons;n#sigma_{#pi};n#sigma_{K}", 100, -25, 25, 100, -25, 25);

   hNSigmaPiecorr = new TH2F("hNSigmaPiecorr","n_{#sigma} pions against electrons;n#sigma_{#pi};n#sigma_{e}", 100, -25, 25, 100, -25, 25);

   hNSigmaPPicorr = new TH2F("hNSigmaPPicorr","n_{#sigma} protons against pions;n#sigma_{p};n#sigma_{#pi}", 100, -25, 25, 100, -25, 25);

   hNSigmaPKcorr = new TH2F("hNSigmaPKcorr","n_{#sigma} protons against kaons;n#sigma_{p};n#sigma_{K}", 100, -25, 25, 100, -25, 25);

   hNSigmaPecorr = new TH2F("hNSigmaPecorr","n_{#sigma} protons against electrons;n#sigma_{p};n#sigma_{e}", 100, -25, 25, 100, -25, 25);

   hNSigmaKPicorr = new TH2F("hNSigmaKPicorr","n_{#sigma} kaons against pions;n#sigma_{K};n#sigma_{#pi}", 100, -25, 25, 100, -25, 25);

   hNSigmaKPcorr = new TH2F("hNSigmaKPcorr","n_{#sigma} kaons against protons;n#sigma_{K};n#sigma_{p}", 100, -25, 25, 100, -25, 25);

   hNSigmaKecorr = new TH2F("hNSigmaKecorr","n_{#sigma} kaons against electrons;n#sigma_{K};n#sigma_{e}", 100, -25, 25, 100, -25, 25);

   hNSigmaPi = new TH1D("hNSigmaPi", "n_{#sigma} pions; n_{#sigma} #pi [-]; counts", 200, -10, 10 );

   hNSigmaP = new TH1D("hNSigmaP", "n_{#sigma} protons; n_{#sigma} p [-]; counts", 200, -10, 10 );

   hNSigmaK = new TH1D("hNSigmaK", "n_{#sigma} kaons; n_{#sigma} K [-]; counts", 200, -10, 10 );

   hDEdxSignal = new TH1D("hDEdxSignal", "dE/dx signal; ln#frac{dE}{dx} [MeV/cm];counts", 200, 0, 10 );

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

   hNTracksBEMC = new TH1D("hNTracksBEMC", "Number of Tracks in BEMC per event; Number of tracks in BEMC; counts", 21, -0.5, 20.5);

   hPt = new TH1D("hPt", "Transverse momentum of hadrons", 30, 0, 3);
   hPt->SetTitle("Distribution of p_{T}");
   hPt->GetXaxis()->SetTitle("p_{T} [GeV]");
   hPt->GetYaxis()->SetTitle(YAxisDescription);


   hNfitHits = new TH1D("hNfitHits", "NfitHits; N^{fit}_{hits} [-]; counts", 50, 0, 50);
   hNfitHitsCut = new TH1D("hNfitHitsCut", "NfitHitsCut; N^{fit}_{hits} [-]; counts", 50, 0, 50);

   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);
   hNhitsDEdxCut = new TH1D("hNhitsDEdxCut", "NhitsDEdxCut; N^{dEdx}_{hits} [-]; counts", 50, 0, 50);

   hNVertices = new TH1D("hNVertices", "Number of vertices before cut", 100, 0, 10);
   hNVertices->SetTitle("Number of vertices");
   hNVertices->GetXaxis()->SetTitle("Number of vertices");
   hNVertices->GetYaxis()->SetTitle(YAxisDescription);

   hTotQ = new TH1D("hTotQ", "Total charge of pair", 3, -1.5, 1.5);
   hTotQ->SetTitle("Total charge of pair");
   hTotQ->GetXaxis()->SetTitle("Charge of pair");
   hTotQ->GetYaxis()->SetTitle(YAxisDescription);

   hInvMassEta = new TH2D("hInvMassEta", "hInvMassEta; m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; #eta [-]", 80, 2,4, 200, -1.,1. );
   hInvMassEta->SetTitle("correlation plot of invMass and eta");

   hNTracksTof = new TH1D("hNTracksTof", "hNTracksTof; Number of ToF tracks [-]; counts", 11, -0.5, 10.5);

   hNTracksTpc = new TH1D("hNTracksTpc", "hNTracksTpc; Number of TPC tracks [-]; counts", 26, -0.5, 25.5);
  
   hInvMassJPsi = new TH1D("hInvMassJPsi", "hInvMassJPsi; m_{e e} [GeV/c^{2}]; counts", 40,2, 4);
   hInvMassJPsiBcg = new TH1D("hInvMassJPsiBcg", "hInvMassJPsiBcg; m_{e e} [GeV/c^{2}]; counts", 40,2, 4);


}

bool AnaJPsi::goodQualityTrack(const StUPCTrack *trk){

   //pT
   hTrackQualityFlow->Fill(1);
   if(trk->getPt() < minPt)
      return false;
   hTrackQualityFlow->Fill(2);
   if(trk->getNhitsFit() < minNHitsFit)
      return false;
   hTrackQualityFlow->Fill(3);
   if(trk->getNhitsDEdx() < minNHitsDEdx)
      return false;
   hTrackQualityFlow->Fill(4);
   tpcCounter += 1;

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

   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
}

void AnaJPsi::fillTrackQualityCutsAfter(const StUPCTrack* trk){

   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
}
bool AnaJPsi::chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2){

   Double_t chi_ee = pow(trk1->getNSigmasTPCElectron(),2) + pow(trk2->getNSigmasTPCElectron(),2);
   Double_t chi_pp = pow(trk1->getNSigmasTPCProton(),2) + pow(trk2->getNSigmasTPCProton(),2);
   Double_t chi_kk = pow(trk1->getNSigmasTPCKaon(),2) + pow(trk2->getNSigmasTPCKaon(),2);
   Double_t chi_pipi = pow(trk1->getNSigmasTPCPion(),2) + pow(trk2->getNSigmasTPCPion(),2);


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


