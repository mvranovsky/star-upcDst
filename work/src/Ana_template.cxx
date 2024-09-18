#include "Ana_template.h"

Ana_template::Ana_template(TFile *outFile): Ana(outFile){}

void Ana_template::Make()
{
   hAnalysisFlow->Fill(ALL);
   if(!CheckTriggers(&CEPtriggers, mUpcEvt, hTriggerBits))
      return;
   hAnalysisFlow->Fill(TRIG);

   for (int iRp = 0; iRp < nRomanPots; ++iRp)
      for (int iPmt = 0; iPmt < 2; ++iPmt)
      {
         hRpAdc[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
         hRpTac[2*iRp+iPmt]->Fill(mRpEvt->tac(iRp, iPmt));      
         if(mRpEvt->tac(iRp, iPmt) > 200 && mRpEvt->tac(iRp, iPmt) < 1750)
            hRpAdcInWindow[2*iRp+iPmt]->Fill(mRpEvt->adc(iRp, iPmt));
      }

   AnaRpTracks(mRpEvt);

   if( !(mRpTrackIdVec_perSide[E].size()==1 && mRpTrackIdVec_perSide[W].size()==1))
      return;
   hAnalysisFlow->Fill(TWORPTRKS);

   bool protonInRange[nSides];
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      StUPCRpsTrack* trackRP = mRpEvt->getTrack(mRpTrackIdVec_perSide[iSide][0]);
      if(!trackRP) return;
      protonInRange[iSide] = RPInFidRange(trackRP->pVec().X(), trackRP->pVec().Y());
      SaveRPinfo(trackRP, iSide);
   }

   if(!(protonInRange[0] && protonInRange[1])) 
      return;
   /*
// temporary code
   SaveEventInfo(mUpcEvt);
   mRecTree->FillRecTree(); // Fill analysis (reco) Tree
   return;
// tmp code ends here
   */
   hAnalysisFlow->Fill(INFID);
   mRPpTBalance = mRpEvt->getTrack( mRpTrackIdVec_perSide[E][0] )->pVec() + mRpEvt->getTrack( mRpTrackIdVec_perSide[W][0] )->pVec();

   // Skip all events with more vertecies than 1
   if( mUpcEvt->getNumberOfVertices() != 1) return;
   hAnalysisFlow->Fill(ONETOFVX);

   // save event info
   mRunNumber = mUpcEvt->getRunNumber();
   SaveEventInfo(mUpcEvt);
   SaveVertexInfo(mUpcEvt->getVertex(0));

   if( abs(mRecTree->getVertexZInCm()) > vertexRange ) return;
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

/*
   const int nVertecies = mUpcEvt->getNumberOfVertices();

   vector<int> tracksInVertex[nVertecies]; // save the tracks Id of good TOF tracks in that vertex
   int chargeInVertex[nVertecies]; 
   for (int iVertex = 0; iVertex < nVertecies; ++iVertex)
      chargeInVertex[iVertex] = 0;

   for(int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      // Skip all tracks that are not primary or they are not matched with TOF
      // In general, you want to skip all bad tracks
      if( !trk->getFlag(StUPCTrack::kPrimary) || !trk->getFlag(StUPCTrack::kTof) || !mUtil->IsGoodTrack(trk)) 
         continue;

      tracksInVertex[trk->getVertexId()].push_back(trackID);
      chargeInVertex[trk->getVertexId()] += static_cast<int>( trk->getCharge() );
   } 

   for (int iVertex = 0; iVertex < nVertecies; ++iVertex)
   {
      if( tracksInVertex[iVertex].size() != 2)
         continue;
      if( chargeInVertex[iVertex] )
         continue;
      
   }
*/
   /*
   if( (hadronId.size() == 1 || hadronId.size() == 2) )//&& nTpcGoodTracks < 5)
      CalculateTOFEff( hadronId[0] );
   */
   // Skip events with other than 2 TOF tracks
   if(hadronId.size() != nSigns)
      return;

   hAnalysisFlow->Fill(TWOTOFTRKS);

   TLorentzVector hadron[nParticles], state[nParticles];
   for (unsigned int id = 0; id < hadronId.size(); ++id)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(hadronId[id]);
      SaveTrackInfo(trk, id);
      for (int iPart = 0; iPart < nParticles; ++iPart){
         trk->getLorentzVector(hadron[iPart], mUtil->mass(iPart));
         state[iPart] += hadron[iPart];
      }
   }

   if( abs(mRecTree->getEta(0)) > maxEta || abs(mRecTree->getEta(1)) > maxEta)
      return;

   hAnalysisFlow->Fill(ETA);
   SaveTriggerInfo(mUpcEvt, mRpEvt);
   CalculatePID();

   hPIDStats[0]->Fill(mRecTree->getPairID());

   TVector3 missingMomenta = state[mRecTree->getPairID()].Vect() + mRPpTBalance;
   SaveMissingMomenta(missingMomenta);
   FillMSquared();

   if( mRecTree->getPairID() == -1)
      return;
   // Fill Ana flow
   SaveStateInfo(state[mRecTree->getPairID()]);

   if(totalCharge) // Total charge is not zero => background event
   {
      mRecTree->FillBcgTree(); // Fill background Tree
      return;
   }
   hAnalysisFlow->Fill(OPPOSITE);
   hPIDStats[1]->Fill(mRecTree->getPairID());
   if( missingMomenta.Pt() > exclusivityCut) 
      return;
   
   hAnalysisFlow->Fill(EXCLUSIVE);
   mRecTree->FillRecTree(); // Fill analysis (reco) Tree

   if( mRecTree->getPairID() == PION)
      hAnalysisFlow->Fill(PIPI);
   if( mRecTree->getPairID() == KAON)
      hAnalysisFlow->Fill(KK);
   if( mRecTree->getPairID() == PROTON)
      hAnalysisFlow->Fill(PPBAR);

}

void Ana_template::Init()
{
   if( DEBUG )
      cout<<"Ana_template::Init() called"<<endl;
   mOutFile->cd();
   hAnalysisFlow = new TH1D("AnalysisFlow", "CutsFlow", nAnalysisCuts-1, 1, nAnalysisCuts);
   for(int tb=1; tb<nAnalysisCuts; ++tb) 
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisCutName(tb));

   mOutFile->mkdir("PID")->cd();
   hPIDStats[0] = new TH1D("PIDStatsTotCh", "PID stats before tot. charge cut", 4, -1.5, 2.5);
   hPIDStats[1] = new TH1D("PIDStatsExlsv", "PID stats before pT exclusive cut", 4, -1.5, 2.5);
   for (int pid = 0; pid < 3; ++pid)
      for (int stat = 0; stat < 3; ++stat)
         for (int part = 0; part < 4; ++part)
            hMSquared[pid][stat][part] = new TH1D( Form("hMSquared_%i_%i_%i",pid, stat, part), Form("hMSquared_%i_%i_%i",pid, stat, part), 200, -0.5, 1.5);
   mOutFile->cd();

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   mRecTree = new RecTree(nameOfAna_templateTree, Ana_templateTreeBits, true); 
   mOutFile->mkdir("CPT2noBBCL")->cd();
   for (int iRp = 0; iRp < 2*nRomanPots; ++iRp)
   {
      hRpAdc[iRp]= new TH1D( mUtil->rpName(iRp/2) + Form("_%i_ADC",iRp%2), "ADC", 100, 0, 600);
      hRpAdcInWindow[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_ADCinTAC",iRp%2), "ADC in TAC window", 100, 0, 600);
      hRpTac[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_TAC",iRp%2), "TAC", 100, 0, 2000);
   }

   mOutFile->cd();
   tofEffTree = new TTree("tofEffTree", "tofEffTree");
   tofEffTree->Branch("probePt", &probePt);
   tofEffTree->Branch("probeEta", &probeEta);
   tofEffTree->Branch("probePhi", &probePhi);
   tofEffTree->Branch("probeTof", &probeTof);
   tofEffTree->Branch("probeTofHit", &probeTofHit);

   tofEffTree->Branch("tagPt", &tagPt);
   tofEffTree->Branch("tagEta", &tagEta);
   tofEffTree->Branch("tagPhi", &tagPhi);

   tofEffTree->Branch("pTMiss", &pTMiss); 
   tofEffTree->Branch("invMass", &invMass);
   tofEffTree->Branch("nProbes", &nProbes);
   tofEffTree->Branch("deltaPhiProtons", &deltaPhiProtons); 
   tofEffTree->Branch("pTState", &pTState);
   tofEffTree->Branch("phiState", &phiState); 
}


void Ana_template::CalculateTOFEff(unsigned int tagID)
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

void Ana_template::CalculatePID()
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

Int_t Ana_template::PIDStartegy(int strategy)
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

bool Ana_template::IsPairOf(int type)
{
   // 0 - PION
   // 1 - KAON
   // 2 - PROTON
   if( type < 0 || type > 2)
      return false;

   return (mRecTree->getNSigmaTPC(0, type) < 2 && mRecTree->getNSigmaTPC(1, type) < 2);

}

void Ana_template::SaveMissingMomenta(TVector3 missP)
{
   mRecTree->setPtMissing( missP.Pt() );
   mRecTree->setPxMissing( missP.X() );
   mRecTree->setPyMissing( missP.Y() );   
}

void Ana_template::FillMSquared()
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
