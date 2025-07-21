#include "AnaBP.h"

AnaBP::AnaBP(TFile *outFile): Ana(outFile){}

void AnaBP::Make()
{
   //cerr<< "okeeej" << endl; 
   hAnalysisFlow->Fill(ALL);
   if(!CheckTriggers(&triggerID, mUpcEvt, hTriggerBits))
      return;
   hAnalysisFlow->Fill(TRIG);

   /*
   //what are tac and adc values in RPs?
   for (int iRp = 0; iRp < nRomanPots; ++iRp)
      for (int iPmt = 0; iPmt < 2; ++iPmt)
      {
         hRpAdc[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
         hRpTac[2*iRp+iPmt]->Fill(mRpEvt->tac(iRp, iPmt));      
         if(mRpEvt->tac(iRp, iPmt) > 200 && mRpEvt->tac(iRp, iPmt) < 1750)
            hRpAdcInWindow[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
      }
    */ 


    int numberOfTracksPerBranch[nBranches] = {0, 0, 0, 0};
    vector<int> rpTrackIdVec_perBranch[nBranches];
    vector<int> rpTrackIdVec_perSide[nSides];
    
   
   //clear the values, loop over all tracks and clusters, save to mRPTrackIdVec, ..._perBranch, ...perSide 
   AnaRpTracks(mRpEvt);

   //fill number of RP tracks alltogether
   Int_t RPtracksTot = rpTrackIdVec_perSide[West].size() + rpTrackIdVec_perSide[East].size();
   hNumberRPTracks->Fill(RPtracksTot);
    
   //condition for one RP track on each side 
   if( !(mRpTrackIdVec_perSide[E].size()==1 && mRpTrackIdVec_perSide[W].size()==1))
      return;
   hAnalysisFlow->Fill(TWORPTRKS);

   //condition for RP fiducial region
   bool protonInRange[nSides];
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      StUPCRpsTrack* trackRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][0]);
      if(!trackRP) return;
      hRPcorr[0]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
      if(iSide == 0)
        hRPcorrEast[0]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
      else
        hRPcorrWest[0]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());

      protonInRange[iSide] = RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y());
      SaveRPinfo(trackRP, iSide);

      if (protonInRange[0] && protonInRange[1]){
          hRPcorr[1]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
          if(iSide == 0)
            hRPcorrEast[1]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
          else
            hRPcorrWest[1]->Fill(trackRP->pVec().X(), trackRP->pVec().Y());
      }

   }

   if(!(protonInRange[0] && protonInRange[1])) 
      return;

   hAnalysisFlow->Fill(INFID);

   //momentum sum of protons in RPs
   mRPpTBalance = mRpEvt->getTrack( mRpTrackIdVec_perSide[E][0] )->pVec() + mRpEvt->getTrack( mRpTrackIdVec_perSide[W][0] )->pVec();

   //condition for number of vertices
   Int_t nVertices = mUpcEvt->getNumberOfVertices();
   hNVertices->Fill(nVertices);
   if( nVertices != 1) return;
   hAnalysisFlow->Fill(ONETOFVX);

   // save event info
   mRunNumber = mUpcEvt->getRunNumber();
   SaveEventInfo(mUpcEvt);
   TVector3 vtx;
   vtx.SetXYZ(mUpcEvt->getVertex(0)->getPosX(), mUpcEvt->getVertex(0)->getPosY(), mUpcEvt->getVertex(0)->getPosZ());
   SaveVertexInfo(vtx, 0);

   //condition for z-position of vertex 
   hPosZ->Fill(mRecTree->getVertexZInCm(0));
   if( abs(mRecTree->getVertexZInCm(0)) > vertexRange ) return;
   hAnalysisFlow->Fill(ZVERTEX);

   unsigned int vertexID = mUpcEvt->getVertex(0)->getId();
   int totalCharge = 0;
   unsigned int nTpcGoodTracks = 0;
   vector<int> hadronId;

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // Or they are originating from different vertex than selected
      // In general, you want to skip all bad tracks

      fillGoodTrackCuts(trk);
      if( !trk->getFlag(StUPCTrack::kPrimary) || trk->getVertexId() != vertexID || !IsGoodTrack(trk)) 
         continue;

      // count number of good quality TPC tracks in the vertex
      nTpcGoodTracks++;

      if( !trk->getFlag(StUPCTrack::kTof) || !IsGoodTofTrack(trk)) 
         continue;

      hadronId.push_back(trackID);
      totalCharge += static_cast<int>( trk->getCharge() );
   } 
   mRecTree->setNGoodTpcTrks(nTpcGoodTracks);

   //condition on number of TOF tracks(usually 2)
   hTOFTracks->Fill(hadronId.size());
   if(hadronId.size() != nSigns)
      return;

   hAnalysisFlow->Fill(TWOTOFTRKS);

   //pseudorapidity cut for both tracks
   if(mUpcEvt->getTrack(hadronId[0])->getEta() > maxEta ||  mUpcEvt->getTrack(hadronId[1])->getEta() > maxEta )
      return;
   hAnalysisFlow->Fill(ETA);


   TLorentzVector hadron[nParticles], state[nNeutralParticles];
   for (unsigned int id = 0; id < hadronId.size(); ++id){
      const StUPCTrack* trk = mUpcEvt->getTrack(hadronId[id]);
      SaveTrackInfo(trk, id);
      saveNSigmaCorr(trk);

   }
   //identify if pairs are pion pairs, or proton pion pairs

   SaveTriggerInfo(mUpcEvt, mRpEvt);

   const StUPCTrack* hadron1 = mUpcEvt->getTrack(hadronId[0]);
   const StUPCTrack* hadron2 = mUpcEvt->getTrack(hadronId[1]);
   //important part of code which differentiates between like and unlike signs and specific particle pairs
   if(oppositePair(hadron1, hadron2)){
       hAnalysisFlow->Fill(OPPOSITE);
       if (is2pions(hadron1, hadron2)){
           hAnalysisFlow->Fill(PIPI);
           mRecTree->setPairID(K0S,0);
           hadron1->getLorentzVector(hadron[0], mUtil->mass(PION));
           hadron2->getLorentzVector(hadron[1], mUtil->mass(PION));
           state[K0S] = hadron[0] + hadron[1];

       }else if(isProtonPion(hadron1, hadron2) && lambda(hadron1, hadron2)){
           hAnalysisFlow->Fill(PPI); 
           mRecTree->setPairID(LAMBDA,0);
           if(hasGoodTPCnSigma(hadron1) == PROTON){
               hadron1->getLorentzVector(hadron[0], mUtil->mass(PROTON));
               hadron2->getLorentzVector(hadron[1], mUtil->mass(PION));
           }else{
               hadron1->getLorentzVector(hadron[0], mUtil->mass(PION));
               hadron2->getLorentzVector(hadron[1], mUtil->mass(PROTON));
           }
           state[LAMBDA] = hadron[0] + hadron[1];
       }else if(isProtonPion(hadron1, hadron2) && !lambda(hadron1, hadron2)){
           hAnalysisFlow->Fill(PIPBAR); 
           mRecTree->setPairID(LAMBDABAR,0);
           if(hasGoodTPCnSigma(hadron1) == PROTON){
               hadron1->getLorentzVector(hadron[0], mUtil->mass(PROTON));
               hadron2->getLorentzVector(hadron[1], mUtil->mass(PION));
           }else{
               hadron1->getLorentzVector(hadron[0], mUtil->mass(PION));
               hadron2->getLorentzVector(hadron[1], mUtil->mass(PROTON));
           }
           state[LAMBDABAR] = hadron[0] + hadron[1];
       }else return;
   } else{
       mRecTree->FillBcgTree();
       return;
   }

   if( mRecTree->getPairID(0) == -1)
      return;

   SaveStateInfo(state[mRecTree->getPairID(0)], totalCharge,0);
   TVector3 missingMomenta = state[mRecTree->getPairID(0)].Vect() + mRPpTBalance;
   SaveMissingMomenta(missingMomenta);
   mRecTree->FillRecTree();

}

void AnaBP::Init()
{
   mUtil = new Util();
 
   if( DEBUG )
   cout<<"AnaBP::Init() called"<<endl;
   
   mOutFile->cd();
   mOutFile->mkdir(nameOfAnaBPDir);
   mOutFile->cd(nameOfAnaBPDir);

   hAnalysisFlow = new TH1D("AnalysisFlow", "CutsFlow", nAnalysisCuts-1, 1, nAnalysisCuts);
   
   for(int tb=1; tb<nAnalysisCuts; ++tb) {
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisCutName(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   mRecTree = new RecTree(nameOfAnaBPTree, AnaBPTreeBits, true); 
   

   hNSigmaPiPcorr = new TH2F("hNSigmaPiPcorr","n_{#sigma} pions against protons", 100, -25, 25, 100, -25, 25);
   hNSigmaPiPcorr->GetXaxis()->SetTitle("n#sigma_{#pi}");
   hNSigmaPiPcorr->GetYaxis()->SetTitle("n#sigma_{p}");

   hNSigmaPiKcorr = new TH2F("hNSigmaPiKcorr","n_{#sigma} pions against kaons", 100, -25, 25, 100, -25, 25);
   hNSigmaPiKcorr->GetXaxis()->SetTitle("n#sigma_{#pi}");
   hNSigmaPiKcorr->GetYaxis()->SetTitle("n#sigma_{K}");

   hNSigmaPiecorr = new TH2F("hNSigmaPiecorr","n_{#sigma} pions against electrons", 100, -25, 25, 100, -25, 25);
   hNSigmaPiecorr->GetXaxis()->SetTitle("n#sigma_{#pi}");
   hNSigmaPiecorr->GetYaxis()->SetTitle("n#sigma_{e}");

   hNSigmaPUPCrr = new TH2F("hNSigmaPUPCrr","n_{#sigma} protons against pions", 100, -25, 25, 100, -25, 25);
   hNSigmaPUPCrr->GetXaxis()->SetTitle("n#sigma_{p}");
   hNSigmaPUPCrr->GetYaxis()->SetTitle("n#sigma_{#pi}");

   hNSigmaPKcorr = new TH2F("hNSigmaPKcorr","n_{#sigma} protons against kaons", 100, -25, 25, 100, -25, 25);
   hNSigmaPKcorr->GetXaxis()->SetTitle("n#sigma_{p}");
   hNSigmaPKcorr->GetYaxis()->SetTitle("n#sigma_{K}");

   hNSigmaPecorr = new TH2F("hNSigmaPecorr","n_{#sigma} protons against electrons", 100, -25, 25, 100, -25, 25);
   hNSigmaPecorr->GetXaxis()->SetTitle("n#sigma_{p}");
   hNSigmaPecorr->GetYaxis()->SetTitle("n#sigma_{e}");

   hNSigmaKUPCrr = new TH2F("hNSigmaKUPCrr","n_{#sigma} kaons against pions", 100, -25, 25, 100, -25, 25);
   hNSigmaKUPCrr->GetXaxis()->SetTitle("n#sigma_{K}");
   hNSigmaKUPCrr->GetYaxis()->SetTitle("n#sigma_{#pi}");

   hNSigmaKPcorr = new TH2F("hNSigmaKPcorr","n_{#sigma} kaons against protons", 100, -25, 25, 100, -25, 25);
   hNSigmaKPcorr->GetXaxis()->SetTitle("n#sigma_{K}");
   hNSigmaKPcorr->GetYaxis()->SetTitle("n#sigma_{p}");

   hNSigmaKecorr = new TH2F("hNSigmaKecorr","n_{#sigma} kaons against electrons", 100, -25, 25, 100, -25, 25);
   hNSigmaKecorr->GetXaxis()->SetTitle("n#sigma_{K}");
   hNSigmaKecorr->GetYaxis()->SetTitle("n#sigma_{e}");


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


   hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks",100,-2,2,100,-4,4);
   hEtaPhi->GetXaxis()->SetTitle("#eta");
   hEtaPhi->GetYaxis()->SetTitle("#varphi");
   hEtaPhi->SetTitle("Pseudorapidity and azimuthal angle distribution");

   hNumberRPTracks = new TH1D("NumberRPTracks", "Number of Tracks in RPs per event", 400, 0, 40);
   hNumberRPTracks->GetXaxis()->SetTitle("Number of tracks in RPs");
   hNumberRPTracks->GetYaxis()->SetTitle(YAxisDescription);
   hNumberRPTracks->SetTitle("Number of tracks in RPs per event");

   hGlobalTracks = new TH1D("NumberGlobalTracks", "Number of Tracks in RPs per event;n^{global}_tracks [-]", 50, 0, 50);
   hGlobalTracks->GetYaxis()->SetTitle(YAxisDescription);

   hPt = new TH1D("hPt", "Transverse momentum of hadrons", 30, 0, 3);
   hPt->SetTitle("Distribution of p_{T}");
   hPt->GetXaxis()->SetTitle("p_{T} [GeV]");
   hPt->GetYaxis()->SetTitle(YAxisDescription);

   hEta = new TH1D("hEta", "Pseudorapidity", 60, -2, 2);
   hEta->SetTitle("Distribution of pseudorapidity ");
   hEta->GetXaxis()->SetTitle("#eta");
   hEta->GetYaxis()->SetTitle(YAxisDescription);

   hDcaZ =  new TH1D("hDcaZ", "DcaZ", 100, -3,3); 
   hDcaZ->SetTitle("Distribution of DCA_{z} ");
   hDcaZ->GetXaxis()->SetTitle("DCA_{z} [cm]");
   hDcaZ->GetYaxis()->SetTitle(YAxisDescription);

   hDcaXY =  new TH1D("hDcaXY", "DcaXY", 60, -0.5,5.5); 
   hDcaXY->SetTitle("Distribution of DCA_{xy} ");
   hDcaXY->GetXaxis()->SetTitle("DCA_{xy} [cm]");
   hDcaXY->GetYaxis()->SetTitle(YAxisDescription);
 
   hNfitHits = new TH1D("hNfitHits", "NfitHits", 500, -0.5, 49.5);
   hNfitHits->SetTitle("Distribution of number of TPC hits");
   hNfitHits->GetXaxis()->SetTitle("hits");
   hNfitHits->GetYaxis()->SetTitle(YAxisDescription);
 
   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx", 500, -0.5, 49.5);
   hNhitsDEdx->SetTitle("Distribution of number of hits dE/dx ");
   hNhitsDEdx->GetXaxis()->SetTitle("Number of hits dE/dx");
   hNhitsDEdx->GetYaxis()->SetTitle(YAxisDescription);
   
   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}", 300, -300, 300); 
   hPosZ->SetTitle("Distribution of position of z_{vertex}");
   hPosZ->GetXaxis()->SetTitle("Z_{vertex} [cm]");
   hPosZ->GetYaxis()->SetTitle(YAxisDescription);

   hTOFTracks = new TH1D("hTOFtracks", "Number of tracks in TOF before cut", 100, 0, 10);
   hTOFTracks->SetTitle("Number of tracks in TOF");
   hTOFTracks->GetXaxis()->SetTitle("Number of tracks");
   hTOFTracks->GetYaxis()->SetTitle(YAxisDescription);
 
   hNVertices = new TH1D("hNVertices", "Number of vertices before cut", 100, 0, 10);
   hNVertices->SetTitle("Number of vertices");
   hNVertices->GetXaxis()->SetTitle("Number of vertices");
   hNVertices->GetYaxis()->SetTitle(YAxisDescription);

   hTotQ = new TH1D("TotQ", "Total charge of pair", 1000, -10, 10);
   hTotQ->SetTitle("Total charge of pair");
   hTotQ->GetXaxis()->SetTitle("Charge of pair");
   hTotQ->GetYaxis()->SetTitle(YAxisDescription);


}


bool AnaBP::lambda(const StUPCTrack *trk1,const StUPCTrack *trk2){
    if (hasGoodTPCnSigma(trk1) == PROTON && trk1->getCharge() > 0)
        return true;
    else if(hasGoodTPCnSigma(trk2) == PROTON && trk2->getCharge() > 0)
        return true;
    else 
        return false;
}



bool AnaBP::oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2){
    Double_t charge;
    charge = trk1->getCharge() + trk2->getCharge();
    if(charge == 0)
        return true;
    else
        return false;
}



bool AnaBP::is2pions(const StUPCTrack *trk1,const StUPCTrack *trk2){

    if((hasGoodTPCnSigma(trk1) == PION) && (hasGoodTPCnSigma(trk2) == PION)){
        return true;
    }
    else
        return false;

}
bool AnaBP::isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2){

   if( hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PROTON ){
       return true;
   }else if(hasGoodTPCnSigma(trk1) == PROTON && hasGoodTPCnSigma(trk2) == PION){
       return true;     
   }
   else
       return false;

}


int AnaBP::hasGoodTPCnSigma(const StUPCTrack *trk){
  //check for good nSigma of hadron in TPC, 10 means unidentified
    if((trk->getNSigmasTPCProton() < 3) && (trk->getNSigmasTPCKaon() > 3) && (trk->getNSigmasTPCPion()>3))
        return PROTON;
    else if((trk->getNSigmasTPCProton() > 3) && (trk->getNSigmasTPCKaon() < 3) && (trk->getNSigmasTPCPion()>3))
        return KAON;
    else if(trk->getNSigmasTPCPion() < 3)
        return PION;
    else
        return 10;

}



void AnaBP::fillGoodTrackCuts(const StUPCTrack* trk){
    hDcaZ->Fill(trk->getDcaZ());
    hDcaXY->Fill(trk->getDcaXY());
    hNfitHits->Fill(trk->getNhitsFit());
    hNhitsDEdx->Fill(trk->getNhitsDEdx());
}


void AnaBP::saveNSigmaCorr(const StUPCTrack *trk){
    Double_t nSigmaPion, nSigmaProton, nSigmaKaon, nSigmaElectron;
    nSigmaPion = trk->getNSigmasTPCPion();
    nSigmaProton = trk->getNSigmasTPCProton();
    nSigmaKaon = trk->getNSigmasTPCKaon();
    nSigmaElectron = trk->getNSigmasTPCElectron();
    hNSigmaPiPcorr->Fill(nSigmaPion, nSigmaProton);
    hNSigmaPiKcorr->Fill(nSigmaPion, nSigmaKaon);
    hNSigmaPiecorr->Fill(nSigmaPion, nSigmaElectron);
    hNSigmaPKcorr->Fill(nSigmaProton, nSigmaKaon);
    hNSigmaPecorr->Fill(nSigmaProton, nSigmaElectron);
    hNSigmaPUPCrr->Fill(nSigmaProton, nSigmaPion);
    hNSigmaKecorr->Fill(nSigmaKaon, nSigmaElectron);
    hNSigmaKPcorr->Fill(nSigmaKaon, nSigmaProton);
    hNSigmaKUPCrr->Fill(nSigmaKaon, nSigmaPion);
}

/*
void AnaBP::CalculateTOFEff(unsigned int tagID)
{
   // eta cut and PID not implemented
   unsigned int vertexID = mUpcEvt->getVertex(0)->getId();
   vector<int> hadronId;

   TLorentzVector tag, probe, state;
   const StUPCTrack* tagTrk = mUpcEvt->getTrack(tagID);
   int tagCharge = static_cast<int>( tagTrk->getCharge() );
   tagTrk->getLorentzVector(tag, mUtil->mass(PION));
   if( abs(tagTrk->getNSigmasTPC(StUPCTrack::kPion)) > 3)
      return;
   if( abs(tagTrk->getEta()) > 0.7 )
      return;

   for(unsigned int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      if( trackID == tagID)
         continue;

      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      if( !trk->getFlag(StUPCTrack::kPrimary) || trk->getVertexId() != vertexID || !IsGoodTrack(trk)) 
         continue;

      if( abs(trk->getNSigmasTPC(StUPCTrack::kPion)) > 3)
         continue;

      if( static_cast<int>( trk->getCharge() ) == -tagCharge  );
         hadronId.push_back(trackID);
   } 

   double missPt, minMissPt;
   unsigned int probeID;
   minMissPt = 999;
   nProbes = hadronId.size();
   for (unsigned int id = 0; id < hadronId.size(); ++id)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(hadronId[id]);

      trk->getLorentzVector(probe, mUtil->mass(PION));
      state = tag + probe;

      missPt = (state.Vect() + mRPpTBalance).Pt();
      if( missPt < minMissPt){
         minMissPt = missPt;
         probeID = hadronId[id];
      }
   }
   
   if( minMissPt > 2 )
      return;

   const StUPCTrack* probeTrk = mUpcEvt->getTrack(probeID);
   probeTrk->getLorentzVector(probe, mUtil->mass(PION));
   state = tag + probe;

   deltaPhiProtons = abs(mRecTree->getPhiRp(East) - mRecTree->getPhiRp(West)); 
   pTState = state.Pt();
   phiState = abs(state.Phi());
   double deltaPhiPions = abs(tag.Phi() - probe.Phi());
   double piHalf = mUtil->pi() / 2;

   
   if( deltaPhiProtons > piHalf){ // is elastic
      if(pTState > 0.7)
         return;
      if( deltaPhiPions < 2.36)
         return;
   }else{ // is inelastic
      if( pTState < 0.7)
         return;
      if( phiState < 1.14 || phiState > 2.0)
         return;
      if( deltaPhiPions > 2.36)
         return;
   }

   // Fill sample total
   pTMiss = minMissPt;
   invMass = state.M();

   tagPt = tagTrk->getPt();
   tagEta = tagTrk->getEta();
   tagPhi = tagTrk->getPhi();

   probePt = probeTrk->getPt();
   probeEta = probeTrk->getEta();
   probePhi = probeTrk->getPhi();
   probeTof = probeTrk->getFlag(StUPCTrack::kTof);
   probeTofHit = IsGoodTofTrack(probeTrk);

   tofEffTree->Fill();
}

void AnaBP::CalculatePID()
{

   Int_t pairID = -1;
   double chiPair[nParticles];
   for(int iPart = 0; iPart < nParticles; ++iPart)
      chiPair[iPart] = mRecTree->getNSigmaTPC(0, iPart)*mRecTree->getNSigmaTPC(0, iPart) + 
                     mRecTree->getNSigmaTPC(1, iPart)*mRecTree->getNSigmaTPC(1, iPart);

   double deltaTOF = mRecTree->getTofTimeInNs(1) - mRecTree->getTofTimeInNs(0);
   double speedOfLight2 = mUtil->c()*mUtil->c();
   double speedOfLight4 = speedOfLight2*speedOfLight2; 
   double length1Squared = mRecTree->getTofLengthInCm(0)*mRecTree->getTofLengthInCm(0)/(100*100); // convert TOFlength from cm to m
   double length2Squared = mRecTree->getTofLengthInCm(1)*mRecTree->getTofLengthInCm(1)/(100*100); // convert TOFlength from cm to m
   double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
   double deltaTime4 = deltaTime2*deltaTime2;
   double oneOverMomentum1sq = 1/(mRecTree->getMomentumInGeV(0)*mRecTree->getMomentumInGeV(0));
   double oneOverMomentum2sq = 1/(mRecTree->getMomentumInGeV(1)*mRecTree->getMomentumInGeV(1));
   double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
   double bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
   double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
   double mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);
   double pT[] = { mRecTree->getPtInGev(0), mRecTree->getPtInGev(1)};

   if(chiPair[PION] > 9 && chiPair[KAON] > 9 && chiPair[PROTON] < 9 && mSquared > 0.6) // it is... proton!
   {
      if(pT[0] > 0.4 && pT[1] > 0.4 && (pT[0] < 1.1 || pT[1] < 1.1) )
         pairID = PROTON;
   }
   else if(chiPair[PION] > 9 && chiPair[KAON] < 9 && chiPair[PROTON] > 9 && mSquared > 0.15) // it is... kaon!
   {
      if(pT[0] > 0.3 && pT[1] > 0.3 && (pT[0] < 0.7 || pT[1] < 0.7) )
         pairID = KAON;
   }
   else if( chiPair[PION] < 12) // it is... pion!
   {
      pairID = PION;
   }
   
   mRecTree->setPairID(pairID);
   mRecTree->setMSquared(mSquared);
}

Int_t AnaBP::PIDStartegy(int strategy)
{
   // 0 - the main PID stragy without mSquared
   // 1 - the main PID stragy without mSquared and without pT cuts
   // else return -1
   if( strategy < 0 || strategy > 1)
      return -1;

   Int_t pairID = -1;
   double chiPair[nParticles];
   for(int iPart = 0; iPart < nParticles; ++iPart)
      chiPair[iPart] = mRecTree->getNSigmaTPC(0, iPart)*mRecTree->getNSigmaTPC(0, iPart) + 
                     mRecTree->getNSigmaTPC(1, iPart)*mRecTree->getNSigmaTPC(1, iPart);
   double pT[] = { mRecTree->getPtInGev(0), mRecTree->getPtInGev(1)};

   if(chiPair[PION] > 9 && chiPair[KAON] > 9 && chiPair[PROTON] < 9) // it is... proton!
   {
      if(pT[0] > 0.4 && pT[1] > 0.4 && (pT[0] < 1.1 || pT[1] < 1.1) && strategy == 0)
         pairID = PROTON;
      if( strategy == 1)
         pairID = PROTON;
   }
   else if(chiPair[PION] > 9 && chiPair[KAON] < 9 && chiPair[PROTON] > 9) // it is... kaon!
   {
      if(pT[0] > 0.3 && pT[1] > 0.3 && (pT[0] < 0.7 || pT[1] < 0.7) && strategy == 0 )
         pairID = KAON;
      if( strategy == 1)
         pairID = KAON;
   }
   else if( chiPair[PION] < 12) // it is... pion!
   {
      pairID = PION;
   }
   
   return pairID;
}

bool AnaBP::IsPairOf(int type)
{
   // 0 - PION
   // 1 - KAON
   // 2 - PROTON
   if( type < 0 || type > 2)
      return false;

   return (mRecTree->getNSigmaTPC(0, type) < 2 && mRecTree->getNSigmaTPC(1, type) < 2);

}
*/
void AnaBP::SaveMissingMomenta(TVector3 missP)
{
   mRecTree->setPtMissing( missP.Pt() );
   mRecTree->setPxMissing( missP.X() );
   mRecTree->setPyMissing( missP.Y() );   
}

/*
void AnaBP::FillMSquared()
{
   int totCharge = mRecTree->getCharge(0) + mRecTree->getCharge(1);
   double pTMissing = mRecTree->getPtMissing();
   double mSquared = mRecTree->getMSquared();
   
   // hMSquared[strategy][state][particle]
   // particle: All Pion Kaon Proton
   // state: before tot charge | after tot charge | exclusive
   // strategy: standard PID | standard PID w/O pT | selective dEdx 

   // fill hMSquared for "before tot charge cut"
   for (int iState = 0; iState < 3; ++iState)
   {
      // iState == 0 -> fill hMSquared for "before tot charge cut"
      // iState == 1 -> fill hMSquared for "after tot charge cut"
      // iState == 2 -> fill hMSquared for "after exclusive cut"

      if( iState == 1 && totCharge != 0 )
         continue;

      if( iState == 2 && pTMissing > exclusivityCut )
         continue;

      for (int iStrategy = 0; iStrategy < 3; ++iStrategy) // for each strategy
         hMSquared[iStrategy][iState][0]->Fill(mSquared);

      int partTypeStrategyOne = PIDStartegy(0) + 1;
      int partTypeStrategyTwo = PIDStartegy(1) + 1;

      if( partTypeStrategyOne > 0 )
         hMSquared[0][iState][partTypeStrategyOne]->Fill(mSquared);

      if( partTypeStrategyTwo > 0 )
         hMSquared[1][iState][partTypeStrategyTwo]->Fill(mSquared);

      for(int iPart = 0; iPart < nParticles; ++iPart)
         if( IsPairOf(iPart))
            hMSquared[2][iState][iPart+1]->Fill(mSquared);
   }
}
*/