#include "../include/RecTree.h"

//_____________________________________________________________________________
RecTree::RecTree(TString treeName, bitset<16> treeVersion, bool isBcgTree) {
   // treeVersion bits: 11111111 = all information is stored
   // 0000 0001 - event info
   // 0000 0010 - Vertex info
   // ...

   mUtil = new Util();
   //standard reconstructed tree
   mRecTree = new TTree(treeName, treeName);

   // event info
   mRecTree->Branch("runNumber", &mRunNumber);

   if( treeVersion.test(0) )
   {
      mRecTree->Branch("nVertecies", &mNVertecies);
      mRecTree->Branch("nGoodTpcTrks", &mNGoodTpcTrks);
      mRecTree->Branch("tofMult", &mTofMult);
      mRecTree->Branch("eventNumber", &mEventNumber);
      mRecTree->Branch("fillNumber", &mFillNumber);
      mRecTree->Branch("bunchCrossId", &mBunchCrossId);
      mRecTree->Branch("bunchCrossId7bit", &mBunchCrossId7bit);
   }
   // Vertex info
   if( treeVersion.test(1) )
   {
      for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
      {
         TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
         mRecTree->Branch("vertexZInCm" + tag, &mVertexZInCm[0][iTag]);
         mRecTree->Branch("vertexYInCm" + tag, &mVertexYInCm[0][iTag]);
         mRecTree->Branch("vertexXInCm" + tag, &mVertexXInCm[0][iTag]);
         if( iTag == TRUEMC )
            mRecTree->Branch("vertexIdTruth" + tag, &mVertexIdTruth);
      }
   }

   // State info
   if( treeVersion.test(2) )
   {
      mRecTree->Branch("invMass", &mInvMass[0]);
      mRecTree->Branch("theta", &mTheta[0]);
      mRecTree->Branch("phi", &mPhi[0]);
      mRecTree->Branch("p", &mP[0]);
      mRecTree->Branch("pt", &mPt[0]);
      mRecTree->Branch("pairRapidity", &mRap[0]);
      mRecTree->Branch("mSquared", &mMSquared[0]);
      mRecTree->Branch("pairID", &mPairID[0]);
      mRecTree->Branch("totQ", &mTotQ[0]);
      if( treeVersion.test(11) ){ // 11 == add chi square for pid
         for (int iPart = 0; iPart < nParticles; ++iPart)
            mRecTree->Branch("chiSquare" + mUtil->particleName(iPart), &mChiSquare[iPart]);
      }
   }
  // Central hadrons info
   if( treeVersion.test(3) )
   {
      for (int i = 0; i < nSigns; ++i)
      { 
         for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
         {
            TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
            mRecTree->Branch(Form("momentumInGev%i",i) + tag, &mMomentumInGev[i][iTag]);
            mRecTree->Branch(Form("pTInGev%i",i) + tag, &mPtInGev[i][iTag]);
            mRecTree->Branch(Form("pXInGev%i",i) + tag, &mPxInGev[i][iTag]);
            mRecTree->Branch(Form("pYInGev%i",i) + tag, &mPyInGev[i][iTag]);
            mRecTree->Branch(Form("pZInGev%i",i) + tag, &mPzInGev[i][iTag]);
            mRecTree->Branch(Form("charge%i",i) + tag, &mCharge[i][iTag]);
            mRecTree->Branch(Form("etaHadron%i",i) + tag, &mEtaHadrons[i][iTag]);
            mRecTree->Branch(Form("phiHadron%i",i) + tag, &mPhiHadrons[i][iTag]);
            mRecTree->Branch(Form("tofHit%i",i)+ tag, &mTofHit[i]);
            if( iTag == TRUEMC )
               mRecTree->Branch(Form("QAHadron%i",i) + tag, &mQATruth[i]);

            if( treeVersion.test(12) ){ // 12 == BEMC info
               mRecTree->Branch(Form("bemcEta%i", i), &mBemcEta[i]);
               mRecTree->Branch(Form("bemcPt%i", i), &mBemcPt[i]);
               mRecTree->Branch(Form("bemcPhi%i", i), &mBemcPhi[i]);
               mRecTree->Branch(Form("bemcE%i", i), &mBemcE[i]);
            }
         } 

         mRecTree->Branch(Form("dEdxInKevCm%i",i), &mDEdxInKevCm[i]);
         mRecTree->Branch(Form("tofTimeInNs%i",i), &mTofTimeInNs[i]);
         mRecTree->Branch(Form("tofLengthInCm%i",i), &mTofLengthInCm[i]);
         mRecTree->Branch(Form("dcaXYInCm%i",i), &mDcaXYInCm[i]);
         mRecTree->Branch(Form("dcaZInCm%i",i), &mDcaZInCm[i]);
         mRecTree->Branch(Form("nHitsFit%i",i), &mNHitsFit[i]);
         mRecTree->Branch(Form("nHitsDEdx%i",i), &mNHitsDEdx[i]);
         for (int iPart = 0; iPart < nParticles; ++iPart)
            mRecTree->Branch("nSigmaTPC" + mUtil->particleName(iPart) + mUtil->signName(i), &mNSigmaTPC[i][iPart]);
      }      
   }

   // RP track info 
   if( treeVersion.test(4) )
   {
      for (int i = 0; i < nSides; ++i)
      {
         mRecTree->Branch("rpX" + mUtil->sideName(i), &mRpX[i]);
         mRecTree->Branch("rpY" + mUtil->sideName(i), &mRpY[i]);
         mRecTree->Branch("rpZ" + mUtil->sideName(i), &mRpZ[i]);
         mRecTree->Branch("thetaRp" + mUtil->sideName(i), &mThetaRp[i]);
         mRecTree->Branch("phiRp" + mUtil->sideName(i), &mPhiRp[i]);
         mRecTree->Branch("timeRp" + mUtil->sideName(i), &mTimeRp[i]);
         mRecTree->Branch("pRp" + mUtil->sideName(i), &mPRp[i]);
         mRecTree->Branch("ptRp" + mUtil->sideName(i), &mPtRp[i]);
         mRecTree->Branch("etaRp" + mUtil->sideName(i), &mEtaRp[i]);
         mRecTree->Branch("pXRp" + mUtil->sideName(i), &mPx[i]);
         mRecTree->Branch("pYRp" + mUtil->sideName(i), &mPy[i]);
         mRecTree->Branch("pZRp" + mUtil->sideName(i), &mPz[i]);
         mRecTree->Branch("t" + mUtil->sideName(i), &mT[i]);
         mRecTree->Branch("xi" + mUtil->sideName(i), &mXi[i]);
      }
   }


   // pT Missing info (if cental hadrons and RP track info)
   if( treeVersion.test(3) && treeVersion.test(4) )
   {
      mRecTree->Branch("pTMissing", &mPtMissing);
      mRecTree->Branch("pXMissing", &mPxMissing);
      mRecTree->Branch("pYMissing", &mPyMissing);
   }


   
   // ZDC and BBC info
   if( treeVersion.test(5) )
   {
      mRecTree->Branch("BbcSmallEast", &mBbcSmallEast);
      mRecTree->Branch("BbcSmallWest", &mBbcSmallWest);
      mRecTree->Branch("BbcLargeEast", &mBbcLargeEast);
      mRecTree->Branch("BbcLargeWest", &mBbcLargeWest);
      for (unsigned int iPmt = 0; iPmt < 3; ++iPmt)
      {
         mRecTree->Branch(Form("ZdcAdcEastPmt%i",iPmt), &mZdcAdcEastPmt[iPmt]);
         mRecTree->Branch(Form("ZdcAdcWestPmt%i",iPmt), &mZdcAdcWestPmt[iPmt]);
      }
      mRecTree->Branch("ZdcTdcEast", &mZdcTdcEast);
      mRecTree->Branch("ZdcTdcWest", &mZdcTdcWest);
      mRecTree->Branch("ZdcTimeDiff", &mZdcTimeDiff);
      mRecTree->Branch("ZdcVertexZ", &mZdcVertexZ);
      mRecTree->Branch("ZdcEastRate", &mZdcEastRate); 
      mRecTree->Branch("ZdcWestRate", &mZdcWestRate); 
      mRecTree->Branch("ZdcEastUA", &mZdcEastUA); 
      mRecTree->Branch("ZdcWestUA", &mZdcWestUA); 
   }
   // trigger efficiency info
   if( treeVersion.test(6) )
   {
      mRecTree->Branch("BbcDsmBit", &mBbcDsmBit);  
      mRecTree->Branch("BbcSmallEDsmBit", &mBbcSmallEDsmBit);
      mRecTree->Branch("BbcSmallWDsmBit", &mBbcSmallWDsmBit);
      mRecTree->Branch("BbcLargeEDsmBit", &mBbcLargeEDsmBit);
      mRecTree->Branch("BbcLargeWDsmBit", &mBbcLargeWDsmBit);
      mRecTree->Branch("ZdcDsmBit", &mZdcDsmBit);
      mRecTree->Branch("ZdcEDsmBit", &mZdcEDsmBit);
      mRecTree->Branch("ZdcWDsmBit", &mZdcWDsmBit);

      mRecTree->Branch("BbcTrigBit", &mBbcTrigBit);  
      mRecTree->Branch("ZdcTrigBit", &mZdcTrigBit);
      mRecTree->Branch("ZdcETrigBit", &mZdcETrigBit);
      mRecTree->Branch("ZdcWTrigBit", &mZdcWTrigBit);  

      mRecTree->Branch("tofTrigBit", &mTofTrigBit);
      mRecTree->Branch("tofDsmBit", &mTofDsmBit);
      mRecTree->Branch("tofDsmABit", &mTofDsmABit);
      mRecTree->Branch("tofDsmBBit", &mTofDsmBBit);

      mRecTree->Branch("RpTrigBit", &mRpTrigBit);  
      mRecTree->Branch("rpEtTrigBit", &mRpEtTrigBit);  
      mRecTree->Branch("rpItTrigBit", &mRpItTrigBit);  
      mRecTree->Branch("RpDsmBit", &mRpDsmBit);  
      mRecTree->Branch("rpEtDsmBit", &mRpEtDsmBit);
      mRecTree->Branch("rpItDsmBit", &mRpItDsmBit); 
      for (int iRp = 0; iRp < nRomanPots; ++iRp)
      {
         mRecTree->Branch("rp" + mUtil->rpName(iRp) + "DsmBit", &mRpTrigBits[iRp]); 
      }
   }
   


   // Vertex info
   if( treeVersion.test(8) )
   {
      for (int i = 0; i < nStates; ++i){

         for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
         {
            TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
            mRecTree->Branch(Form("vertexZInCm%i",i) + tag, &mVertexZInCm[i][iTag]);
            mRecTree->Branch(Form("vertexYInCm%i",i) + tag, &mVertexYInCm[i][iTag]);
            mRecTree->Branch(Form("vertexXInCm%i",i) + tag, &mVertexXInCm[i][iTag]);
            mRecTree->Branch(Form("dcaDaughtersInCm%i",i) + tag, &mDcaDaughtersInCm[i][iTag]);
            mRecTree->Branch(Form("dcaBeamlineInCm%i",i) + tag, &mDcaBeamlineInCm[i][iTag]);
            mRecTree->Branch(Form("PointingAngle%i",i) + tag, &mPointingAngle[i][iTag]);
            mRecTree->Branch(Form("decayLengthInCm%i",i) + tag, &mDecayLengthInCm[i][iTag]);
            mRecTree->Branch(Form("vertexDiffInCm%i",i) + tag, &mVertexDiffInCm[i][iTag]);
            if( iTag == TRUEMC )
               mRecTree->Branch("vertexIdTruth" + tag, &mVertexIdTruth);
         }   
      }
   }

   // State info
   if( treeVersion.test(9) )
   {
      for (int i = 0; i < nStates; ++i){
      
         mRecTree->Branch(Form("invMass%i",i), &mInvMass[i]);
         mRecTree->Branch(Form("theta%i",i), &mTheta[i]);
         mRecTree->Branch(Form("phi%i",i), &mPhi[i]);
         mRecTree->Branch(Form("p%i",i), &mP[i]);
         mRecTree->Branch(Form("pt%i",i), &mPt[i]);
         mRecTree->Branch(Form("pairRapidity%i",i), &mRap[i]);
         mRecTree->Branch(Form("mSquared%i",i), &mMSquared[i]);
         mRecTree->Branch(Form("pairID%i",i), &mPairID[i]);
         mRecTree->Branch(Form("totQ%i",i), &mTotQ[i]);
      }
   }

   // Central hadrons info
   if( treeVersion.test(10) )
   {
      int j;
      for (int i = 0; i < nHadrons; ++i)
      { 
         for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
         {
            TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
            mRecTree->Branch(Form("momentumInGev%i",i) + tag, &mMomentumInGev[i][iTag]);
            mRecTree->Branch(Form("pTInGev%i",i) + tag, &mPtInGev[i][iTag]);
            mRecTree->Branch(Form("pXInGev%i",i) + tag, &mPxInGev[i][iTag]);
            mRecTree->Branch(Form("pYInGev%i",i) + tag, &mPyInGev[i][iTag]);
            mRecTree->Branch(Form("pZInGev%i",i) + tag, &mPzInGev[i][iTag]);
            mRecTree->Branch(Form("charge%i",i) + tag, &mCharge[i][iTag]);
            mRecTree->Branch(Form("etaHadron%i",i) + tag, &mEtaHadrons[i][iTag]);
            mRecTree->Branch(Form("phiHadron%i",i) + tag, &mPhiHadrons[i][iTag]);
            if( iTag == TRUEMC )
               mRecTree->Branch(Form("QAHadron%i",i) + tag, &mQATruth[i]);
         } 

         mRecTree->Branch(Form("tofHit%i",i), &mTofHit[i]);
         mRecTree->Branch(Form("dEdxInKevCm%i",i), &mDEdxInKevCm[i]);
         mRecTree->Branch(Form("tofTimeInNs%i",i), &mTofTimeInNs[i]);
         mRecTree->Branch(Form("tofLengthInCm%i",i), &mTofLengthInCm[i]);
         mRecTree->Branch(Form("dcaXYInCm%i",i), &mDcaXYInCm[i]);
         mRecTree->Branch(Form("dcaZInCm%i",i), &mDcaZInCm[i]);
         mRecTree->Branch(Form("nHitsFit%i",i), &mNHitsFit[i]);
         mRecTree->Branch(Form("nHitsDEdx%i",i), &mNHitsDEdx[i]);
         j = i % 2 == 0 ? PLUS : MINUS;
         for (int iPart = 0; iPart < nParticles; ++iPart)
            mRecTree->Branch("nSigmaTPC" + mUtil->particleName(iPart) + mUtil->signName(j), &mNSigmaTPC[i][iPart]);
      }      
   }


   // Setting background Tree
   if(isBcgTree){
      mBcgTree = mRecTree->CloneTree(0);
      mBcgTree->SetName(treeName + "_Bcg");
   }
}//RecTree::CreateRecTree


//_____________________________________________________________________________
RecTree::RecTree(TTree* tree, bitset<16> treeVersion) {

   mUtil = new Util();
   // event info
   if( treeVersion.test(0) )
   {
      tree->SetBranchAddress("nVertecies", &mNVertecies);
      tree->SetBranchAddress("nGoodTpcTrks", &mNGoodTpcTrks);
      tree->SetBranchAddress("tofMult", &mTofMult);
      tree->SetBranchAddress("eventNumber", &mEventNumber);
      tree->SetBranchAddress("fillNumber", &mFillNumber);
      tree->SetBranchAddress("bunchCrossId", &mBunchCrossId);
      tree->SetBranchAddress("bunchCrossId7bit", &mBunchCrossId7bit);
   }

   // Vertex info
   if( treeVersion.test(1) )
   {
      for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
      {
         TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
         tree->SetBranchAddress("vertexZInCm" + tag, &mVertexZInCm[iTag]);
         tree->SetBranchAddress("vertexYInCm" + tag, &mVertexYInCm[iTag]);
         tree->SetBranchAddress("vertexXInCm" + tag, &mVertexXInCm[iTag]);
         if( iTag == TRUEMC )
            tree->SetBranchAddress("vertexIdTruth" + tag, &mVertexIdTruth);
      }
   }


   // State info
   if( treeVersion.test(2) )
   {
      tree->SetBranchAddress("invMass", &mInvMass);
      tree->SetBranchAddress("theta", &mTheta);
      tree->SetBranchAddress("phi", &mPhi);
      tree->SetBranchAddress("p", &mP);
      tree->SetBranchAddress("pt", &mPt);
      tree->SetBranchAddress("pairRapidity", &mRap);
      tree->SetBranchAddress("mSquared", &mMSquared);
      tree->SetBranchAddress("pairID", &mPairID);

   }

   // Central hadrons info
   if( treeVersion.test(3) )
   {
      for (int i = 0; i < nStates; ++i)
      { 
         for (int iTag = treeVersion.test(7) ? TRUEMC : RECO; iTag < nDataTag; ++iTag)
         {
            TString tag = iTag == TRUEMC ? "_" + mUtil->dataTagName(iTag) : "";
            tree->SetBranchAddress(Form("momentumInGev%i",i) + tag, &mMomentumInGev[i][iTag]);
            tree->SetBranchAddress(Form("pTInGev%i",i) + tag, &mPtInGev[i][iTag]);
            tree->SetBranchAddress(Form("pXInGev%i",i) + tag, &mPxInGev[i][iTag]);
            tree->SetBranchAddress(Form("pYInGev%i",i) + tag, &mPyInGev[i][iTag]);
            tree->SetBranchAddress(Form("pZInGev%i",i) + tag, &mPzInGev[i][iTag]);
            tree->SetBranchAddress(Form("charge%i",i) + tag, &mCharge[i][iTag]);
            tree->SetBranchAddress(Form("etaHadron%i",i) + tag, &mEtaHadrons[i][iTag]);
            tree->SetBranchAddress(Form("phiHadron%i",i) + tag, &mPhiHadrons[i][iTag]);
            if( iTag == TRUEMC )
               tree->SetBranchAddress(Form("QAHadron%i",i) + tag, &mQATruth[i]);
         } 

         tree->SetBranchAddress(Form("dEdxInKevCm%i",i), &mDEdxInKevCm[i]);
         tree->SetBranchAddress(Form("tofTimeInNs%i",i), &mTofTimeInNs[i]);
         tree->SetBranchAddress(Form("tofLengthInCm%i",i), &mTofLengthInCm[i]);
         tree->SetBranchAddress(Form("dcaXYInCm%i",i), &mDcaXYInCm[i]);
         tree->SetBranchAddress(Form("dcaZInCm%i",i), &mDcaZInCm[i]);
         tree->SetBranchAddress(Form("nHitsFit%i",i), &mNHitsFit[i]);
         tree->SetBranchAddress(Form("nHitsDEdx%i",i), &mNHitsDEdx[i]);

         for (int iPart = 0; iPart < nParticles; ++iPart)
            tree->SetBranchAddress("nSigmaTPC" + mUtil->particleName(iPart) + mUtil->signName(i), &mNSigmaTPC[i][iPart]);
         
         if( treeVersion.test(12) ){ //BEMC info
            tree->SetBranchAddress(Form("bemcEta%i", i), &mBemcEta[i]);
            tree->SetBranchAddress(Form("bemcPt%i", i), &mBemcPt[i]);
            tree->SetBranchAddress(Form("bemcPhi%i", i), &mBemcPhi[i]);
            tree->SetBranchAddress(Form("bemcE%i", i), &mBemcE[i]);
         }
      }
   }

   // RP track info 
   if( treeVersion.test(4) )
   {
      for (int i = 0; i < nSides; ++i)
      {
         tree->SetBranchAddress("rpX" + mUtil->sideName(i), &mRpX[i]);
         tree->SetBranchAddress("rpY" + mUtil->sideName(i), &mRpY[i]);
         tree->SetBranchAddress("rpZ" + mUtil->sideName(i), &mRpZ[i]);
         tree->SetBranchAddress("thetaRp" + mUtil->sideName(i), &mThetaRp[i]);
         tree->SetBranchAddress("phiRp" + mUtil->sideName(i), &mPhiRp[i]);
         tree->SetBranchAddress("timeRp" + mUtil->sideName(i), &mTimeRp[i]);
         tree->SetBranchAddress("pRp" + mUtil->sideName(i), &mPRp[i]);
         tree->SetBranchAddress("ptRp" + mUtil->sideName(i), &mPtRp[i]);
         tree->SetBranchAddress("etaRp" + mUtil->sideName(i), &mEtaRp[i]);
         tree->SetBranchAddress("pXRp" + mUtil->sideName(i), &mPx[i]);
         tree->SetBranchAddress("pYRp" + mUtil->sideName(i), &mPy[i]);
         tree->SetBranchAddress("pZRp" + mUtil->sideName(i), &mPz[i]);
         tree->SetBranchAddress("t" + mUtil->sideName(i), &mT[i]);
         tree->SetBranchAddress("xi" + mUtil->sideName(i), &mXi[i]);
      }
   }

   // pT Missing info (if cental hadrons and RP track info)
   if( treeVersion.test(3) && treeVersion.test(4))
   {
      tree->SetBranchAddress("pTMissing", &mPtMissing);
      tree->SetBranchAddress("pXMissing", &mPxMissing);
      tree->SetBranchAddress("pYMissing", &mPyMissing);
   }

   // ZDC and BBC info
   if( treeVersion.test(5) )
   {
      tree->SetBranchAddress("BbcSmallEast", &mBbcSmallEast);
      tree->SetBranchAddress("BbcSmallWest", &mBbcSmallWest);
      tree->SetBranchAddress("BbcLargeEast", &mBbcLargeEast);
      tree->SetBranchAddress("BbcLargeWest", &mBbcLargeWest);
      for (unsigned int iPmt = 0; iPmt < 3; ++iPmt)
      {
         tree->SetBranchAddress(Form("ZdcAdcEastPmt%i",iPmt), &mZdcAdcEastPmt[iPmt]);
         tree->SetBranchAddress(Form("ZdcAdcWestPmt%i",iPmt), &mZdcAdcWestPmt[iPmt]);
      }
      tree->SetBranchAddress("ZdcTdcEast", &mZdcTdcEast);
      tree->SetBranchAddress("ZdcTdcWest", &mZdcTdcWest);
      tree->SetBranchAddress("ZdcTimeDiff", &mZdcTimeDiff);
      tree->SetBranchAddress("ZdcVertexZ", &mZdcVertexZ);
      tree->SetBranchAddress("ZdcEastRate", &mZdcEastRate); 
      tree->SetBranchAddress("ZdcWestRate", &mZdcWestRate); 
      tree->SetBranchAddress("ZdcEastUA", &mZdcEastUA); 
      tree->SetBranchAddress("ZdcWestUA", &mZdcWestUA); 
   }

   // trigger efficiency info
   if( treeVersion.test(6) )
   {
      tree->SetBranchAddress("BbcDsmBit", &mBbcDsmBit);  
      tree->SetBranchAddress("BbcSmallEDsmBit", &mBbcSmallEDsmBit);
      tree->SetBranchAddress("BbcSmallWDsmBit", &mBbcSmallWDsmBit);
      tree->SetBranchAddress("BbcLargeEDsmBit", &mBbcLargeEDsmBit);
      tree->SetBranchAddress("BbcLargeWDsmBit", &mBbcLargeWDsmBit);
      tree->SetBranchAddress("ZdcDsmBit", &mZdcDsmBit);
      tree->SetBranchAddress("ZdcEDsmBit", &mZdcEDsmBit);
      tree->SetBranchAddress("ZdcWDsmBit", &mZdcWDsmBit);

      tree->SetBranchAddress("BbcTrigBit", &mBbcTrigBit);  
      tree->SetBranchAddress("ZdcTrigBit", &mZdcTrigBit);
      tree->SetBranchAddress("ZdcETrigBit", &mZdcETrigBit);
      tree->SetBranchAddress("ZdcWTrigBit", &mZdcWTrigBit);  

      tree->SetBranchAddress("tofTrigBit", &mTofTrigBit);
      tree->SetBranchAddress("tofDsmBit", &mTofDsmBit);
      tree->SetBranchAddress("tofDsmABit", &mTofDsmABit);
      tree->SetBranchAddress("tofDsmBBit", &mTofDsmBBit);

      tree->SetBranchAddress("RpTrigBit", &mRpTrigBit);  
      tree->SetBranchAddress("rpEtTrigBit", &mRpEtTrigBit);  
      tree->SetBranchAddress("rpItTrigBit", &mRpItTrigBit);  
      tree->SetBranchAddress("RpDsmBit", &mRpDsmBit);  
      tree->SetBranchAddress("rpEtDsmBit", &mRpEtDsmBit);
      tree->SetBranchAddress("rpItDsmBit", &mRpItDsmBit); 
      for (int iRp = 0; iRp < nRomanPots; ++iRp)
      {
         tree->SetBranchAddress("rp" + mUtil->rpName(iRp) + "DsmBit", &mRpTrigBits[iRp]); 
      }
   }

   mRecTree = tree;
}//RecTree::LoadTTree

RecTree::~RecTree(){
   if(mUtil) delete mUtil;
}



void RecTree::InitRPMCInfo()
{
   mRecTree->Branch("rpTrueVertexX", &mTrueVertexX);
   mRecTree->Branch("rpTrueVertexY", &mTrueVertexY);
   mRecTree->Branch("rpTrueVertexZ", &mTrueVertexZ);

   for (int i = 0; i < nSides; ++i)
   {
      mRecTree->Branch("rpTruePx" + mUtil->sideName(i), &mTruePx[i]);
      mRecTree->Branch("rpTruePy" + mUtil->sideName(i), &mTruePy[i]);
      mRecTree->Branch("rpTruePz" + mUtil->sideName(i), &mTruePz[i]);
   }
}

void RecTree::InitRPMCInfo(TTree* tree)
{
   tree->SetBranchAddress("rpTrueVertexX", &mTrueVertexX);
   tree->SetBranchAddress("rpTrueVertexY", &mTrueVertexY);
   tree->SetBranchAddress("rpTrueVertexZ", &mTrueVertexZ);

   for (int i = 0; i < nSides; ++i)
   {
      tree->SetBranchAddress("rpTruePx" + mUtil->sideName(i), &mTruePx[i]);
      tree->SetBranchAddress("rpTruePy" + mUtil->sideName(i), &mTruePy[i]);
      tree->SetBranchAddress("rpTruePz" + mUtil->sideName(i), &mTruePz[i]);
   }
}

void RecTree::InitVertexRecoStudy()
{
   mRecTree->Branch("dcaParticles", &mDcaParticles);
   mRecTree->Branch("dcaBeamline", &mDcaBeamline);
   mRecTree->Branch("primary", &mPrimary);
   mRecTree->Branch("sameVertex", &mSameVertex);
}

void RecTree::InitVertexRecoStudy(TTree* tree)
{
   tree->SetBranchAddress("dcaParticles", &mDcaParticles);
   tree->SetBranchAddress("dcaBeamline", &mDcaBeamline);
   tree->SetBranchAddress("primary", &mPrimary);
   tree->SetBranchAddress("sameVertex", &mSameVertex);
}
/*
void RecTree::CalculatePID()
{

   Int_t pairID = -1;
   double chiPair[nParticles];
   for(int iPart = 0; iPart < nParticles; ++iPart)
      chiPair[iPart] = this->getNSigmaTPC(0, iPart)*this->getNSigmaTPC(0, iPart) + 
                     this->getNSigmaTPC(1, iPart)*this->getNSigmaTPC(1, iPart);

   double deltaTOF = this->getTofTimeInNs(1) - this->getTofTimeInNs(0);
   double speedOfLight2 = mUtil->c()*mUtil->c();
   double speedOfLight4 = speedOfLight2*speedOfLight2; 
   double length1Squared = this->getTofLengthInCm(0)*this->getTofLengthInCm(0)/(100*100); // convert TOFlength from cm to m
   double length2Squared = this->getTofLengthInCm(1)*this->getTofLengthInCm(1)/(100*100); // convert TOFlength from cm to m
   double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
   double deltaTime4 = deltaTime2*deltaTime2;
   double oneOverMomentum1sq = 1/(this->getMomentumInGeV(0)*this->getMomentumInGeV(0));
   double oneOverMomentum2sq = 1/(this->getMomentumInGeV(1)*this->getMomentumInGeV(1));
   double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
   double bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
   double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
   double mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);
   double pT[] = { this->getPtInGev(0), this->getPtInGev(1)};

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
   
   this->setPairID(pairID);
   this->setMSquared(mSquared);
}
*/