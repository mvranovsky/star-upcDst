#ifndef RecTree_h
#define RecTree_h

#include "Libreries.h"
#include "Util.h"

using namespace std;
using namespace UTIL;


class RecTree{
   public: 
      RecTree(TString treeName, bitset<16> treeVersion, bool isBcgTree);
      RecTree(TTree* tree, bitset<16> treeVersion);
      ~RecTree(); 

      void FillRecTree(){ mRecTree->Fill(); }
      void FillBcgTree(){ mBcgTree->Fill(); }

      void InitRPMCInfo();
      void InitRPMCInfo(TTree* tree);

      void InitVertexRecoStudy();
      void InitVertexRecoStudy(TTree* tree);

      //getters
      TTree* getTTree() const { return mRecTree; }
      UInt_t getRunNumber() const { return mRunNumber; }
      UInt_t getEventNumber() const { return mEventNumber; }
      UInt_t getFillNumber() const { return mFillNumber; }
      UInt_t getBunchCrossId() const { return mBunchCrossId; }
      UInt_t getBunchCrossId7bit() const { return mBunchCrossId7bit; }
      UInt_t getNVertecies() const { return mNVertecies; }
      UInt_t getNGoodTpcTrks() const { return mNGoodTpcTrks; }
      UShort_t getTofMult() const { return mTofMult; }
      UInt_t getBbcSmallEast() const { return mBbcSmallEast; }
      UInt_t getBbcSmallWest() const { return mBbcSmallWest; }
      UInt_t getBbcLargeEast() const { return mBbcLargeEast; }
      UInt_t getBbcLargeWest() const { return mBbcLargeWest; }
      UInt_t getZdcAdcEastPmt(UInt_t pmt) const { return mZdcAdcEastPmt[pmt]; } // pmt must be: 0, 1 or 2
      UInt_t getZdcAdcWestPmt(UInt_t pmt) const { return mZdcAdcWestPmt[pmt]; } // pmt must be: 0, 1 or 2
      UShort_t getZdcTdcEast() const { return mZdcTdcEast; }
      UShort_t getZdcTdcWest() const { return mZdcTdcWest; }
      UShort_t getZdcTimeDiff() const { return mZdcTimeDiff; }
      Double_t getZdcVertexZ() const { return mZdcVertexZ; }
      Double_t getZdcEastRate() const { return mZdcEastRate; }
      Double_t getZdcWestRate() const { return mZdcWestRate; }
      UShort_t getZdcUnAttEast() const { return mZdcEastUA; }
      UShort_t getZdcUnAttWest() const { return mZdcWestUA; }
      bool IsVertexSet(unsigned int s, UInt_t dataTag = RECO) const { return (mVertexZInCm[s][dataTag] != -9999); }
      Double_t getVertexZInCm(unsigned int s,UInt_t dataTag = RECO) const { return mVertexZInCm[s][dataTag]; }
      Double_t getVertexYInCm(unsigned int s,UInt_t dataTag = RECO) const { return mVertexYInCm[s][dataTag]; }
      Double_t getVertexXInCm(unsigned int s,UInt_t dataTag = RECO) const { return mVertexXInCm[s][dataTag]; }
      Double_t getDcaDaughters(unsigned int s, UInt_t dataTag = RECO) const { return mDcaDaughtersInCm[s][dataTag]; }
      Double_t getDcaBeamline(unsigned int s, UInt_t dataTag = RECO) const { return mDcaBeamlineInCm[s][dataTag]; }
      Double_t getPointingAngle(unsigned int s, UInt_t dataTag = RECO) const { return mPointingAngle[s][dataTag]; }
      Double_t getDecayLength(unsigned int s, UInt_t dataTag = RECO) const { return mDecayLengthInCm[s][dataTag]; }
      Double_t getVertexDiff(unsigned int s, UInt_t dataTag = RECO) const { return mVertexDiffInCm[s][dataTag]; }
      UInt_t getVertexIdTruth() const { return mVertexIdTruth; }
      Double_t getInvMass(unsigned int s) const { return mInvMass[s]; }
      Double_t getTheta(unsigned int s) const { return mTheta[s]; }
      Double_t getPhi(unsigned int s) const { return mPhi[s]; }
      Double_t getP(unsigned int s) const { return mP[s]; }
      Double_t getPt(unsigned int s) const { return mPt[s]; }
      Double_t getRap(unsigned int s) const { return mRap[s]; }
      Double_t getMSquared(unsigned int s) const { return mMSquared[s]; }
      Int_t getPairID(unsigned int s) const { return mPairID[s]; }
      Int_t getTotQ(unsigned int s) const { return mTotQ[s]; }
      Double_t getPIDChiSquare(unsigned int s) { return mChiSquare[s]; }
      Double_t getMomentumInGeV(unsigned int s, UInt_t dataTag = RECO) const { return mMomentumInGev[s][dataTag]; }   
      Double_t getPtInGev(unsigned int s, UInt_t dataTag = RECO) const { return mPtInGev[s][dataTag]; }  
      Double_t getPxInGev(unsigned int s, UInt_t dataTag = RECO) const { return mPxInGev[s][dataTag]; }
      Double_t getPyInGev(unsigned int s, UInt_t dataTag = RECO)const { return mPyInGev[s][dataTag]; }
      Double_t getPzInGev(unsigned int s, UInt_t dataTag = RECO) const { return mPzInGev[s][dataTag]; }
      Double_t getEta(unsigned int s, UInt_t dataTag = RECO) const { return mEtaHadrons[s][dataTag]; } 
      Double_t getPhi(unsigned int s, UInt_t dataTag = RECO) const { return mPhiHadrons[s][dataTag]; } 
      Double_t getCharge(unsigned int s, UInt_t dataTag = RECO) const { return mCharge[s][dataTag]; } 
      Double_t getDEdxInKevCm(unsigned int s) const { return mDEdxInKevCm[s]; } 
      Int_t getTofHit(unsigned int s) const { return mTofHit[s]; }
      Double_t getTofTimeInNs(unsigned int s) const { return mTofTimeInNs[s]; }   
      Double_t getTofLengthInCm(unsigned int s) const { return mTofLengthInCm[s]; }  
      Double_t getDcaXYInCm(unsigned int s) const { return mDcaXYInCm[s]; } 
      Double_t getDcaZInCm(unsigned int s) const { return mDcaZInCm[s]; }   
      Double_t getNHitsFit(unsigned int s) const { return mNHitsFit[s]; }   
      Double_t getNHitsDEdx(unsigned int s) const { return mNHitsDEdx[s]; } 
      Double_t getNSigmaTPC(unsigned int s, int p) const { return mNSigmaTPC[s][p]; }
      UInt_t   getQATruth(unsigned int s) const { return mQATruth[s]; }
      Double_t getThetaRp(unsigned int s) const { return mThetaRp[s]; }
      Double_t getPhiRp(unsigned int s) const { return mPhiRp[s]; }
      Double_t getTimeRp(unsigned int s) const { return mTimeRp[s]; }
      Double_t getT(unsigned int s) const { return mT[s]; }
      Double_t getPRp(unsigned int s) const { return mPRp[s]; }
      Double_t getPtRp(unsigned int s) const { return mPtRp[s]; }
      Double_t getEtaRp(unsigned int s) const { return mEtaRp[s]; }
      Double_t getRpX(unsigned int s) const { return mRpX[s]; }
      Double_t getRpZ(unsigned int s) const { return mRpZ[s]; }
      Double_t getRpY(unsigned int s) const { return mRpY[s]; }
      Double_t getXi(unsigned int s) const { return mXi[s]; }
      Double_t getPx(unsigned int s) const { return mPx[s]; }
      Double_t getPy(unsigned int s) const { return mPy[s]; }
      Double_t getPz(unsigned int s) const { return mPz[s]; }
      Double_t getTrueRPVertexZInCm() { return mTrueVertexZ; }
      Double_t getTrueRPVertexXInCm() { return mTrueVertexY; }
      Double_t getTrueRPVertexYInCm() { return mTrueVertexX; }
      Double_t getRPTruePx(unsigned int s) { return mTruePx[s]; }
      Double_t getRPTruePy(unsigned int s) { return mTruePy[s]; }
      Double_t getRPTruePz(unsigned int s) { return mTruePz[s]; }
      Double_t getPtMissing() const { return mPtMissing; };
      Double_t getPxMissing() const { return mPxMissing; };
      Double_t getPyMissing() const { return mPyMissing; }; 
      bool getTofTrigBit() const { return mTofTrigBit; }
      bool getRpTrigBit() const { return mRpTrigBit; }
      bool getRpItTrigBit() const { return mRpItTrigBit; }
      bool getRpEtTrigBit() const { return mRpEtTrigBit; }
      bool getBbcTrigBit() const { return mBbcTrigBit; }
      bool getZdcTrigBit() const { return mZdcTrigBit; }
      bool getZdcETrigBit() const { return mZdcETrigBit; }
      bool getZdcWTrigBit() const { return mZdcWTrigBit; }
      bool getTofDsmBit() const { return mTofDsmBit; }
      bool getTofDsmABit() const { return mTofDsmABit; }
      bool getTofDsmBBit() const { return mTofDsmBBit; }
      bool getRpDsmBit() const { return mRpDsmBit; }
      bool getRpItDsmBit() const { return mRpItDsmBit; }
      bool getRpEtDsmBit() const { return mRpEtDsmBit; }
      bool getBbcDsmBit() const { return mBbcDsmBit; }
      bool getBbcSmallEDsmBit() const { return mBbcSmallEDsmBit; }
      bool getBbcSmallWDsmBit() const { return mBbcSmallWDsmBit; }
      bool getBbcLargeEDsmBit() const { return mBbcLargeEDsmBit; }
      bool getBbcLargeWDsmBit() const { return mBbcLargeWDsmBit; }
      bool getZdcDsmBit() const { return mZdcDsmBit; }
      bool getZdcEDsmBit() const { return mZdcEDsmBit; }
      bool getZdcWDsmBit() const { return mZdcWDsmBit; }
      bool getRpTrigBits(unsigned int rp) const { return mRpTrigBits[rp]; }
      bool getVertexStudyPrimary() const { return mPrimary; }
      bool getVertexStudySameVertex() const { return mSameVertex; }
      Double_t getVertexStudyDcaParticles() const { return mDcaParticles; }
      Double_t getVertexStudyDcaBeamline() const { return mDcaBeamline; }
      Double_t getBemcTrackE(unsigned int s) { return mBemcE[s]; }
      Double_t getBemcTrackEta(unsigned int s) { return mBemcEta[s]; }
      Double_t getBemcTrackPhi(unsigned int s) { return mBemcPhi[s]; }
      Double_t getBemcTrackPt(unsigned int s) { return mBemcPt[s]; }
      Double_t getBemcClusterE(unsigned int s) { return mBemcClusterE[s]; }
      Double_t getBemcClusterEta(unsigned int s) { return mBemcClusterEta[s]; }
      Double_t getBemcClusterPhi(unsigned int s) { return mBemcClusterPhi[s]; }
      Double_t getBemcClusterHTE(unsigned int s) { return mBemcClusterHTE[s]; }
      Double_t getBemcClusterSigmaEta(unsigned int s) { return mBemcClusterSigmaEta[s]; }
      Double_t getBemcClusterSigmaPhi(unsigned int s) { return mBemcClusterSigmaPhi[s]; }

      // setters for good run list study
      Int_t getRpOk() const { return mIsRpOk; }
      Int_t getJPsiTrigger1() const { return mIsJPsiTrigger1; }
      Int_t getJPsiTrigger2() const { return mIsJPsiTrigger2; }
      Int_t getJPsiTrigger3() const { return mIsJPsiTrigger3; }
      vector<Double_t> getTpcTrackPhi() const { return mTpcTrack_phi; }
      vector<Double_t> getTpcTrackEta() const { return mTpcTrack_eta; }
      vector<Double_t> getBemcTrackPhi() const { return mBemcTrack_phi; }
      vector<Double_t> getBemcTrackEta() const { return mBemcTrack_eta; }
      vector<Double_t> getTpcNHitsFit() const { return mTpcNHitsFit; }
      vector<Double_t> getTpcNHitsDEdx() const { return mTpcNHitsDEdx; }
      vector<Double_t> getTpcNSigmaElectron() const {return mTpcNSigmaElectron; }
      vector<Double_t> getTpcNSigmaProton() const {return mTpcNSigmaProton; }
      vector<Double_t> getTpcNSigmaKaon() const {return mTpcNSigmaKaon;}
      vector<Double_t> getTpcNSigmaPion() const {return mTpcNSigmaPion; }
      Int_t getNTracksBemc() { return nTracksBemc; }
      Int_t getNTracksTof() { return nTracksTof; }
      Int_t getNClustersBemc() { return nClustersBemc; }


      //setters
      void setRunNumber(UInt_t var) { mRunNumber = var; }
      void setEventNumber(UInt_t var) { mEventNumber = var; }
      void setFillNumber(UInt_t var) { mFillNumber = var; }
      void setBunchCrossId(UInt_t var) { mBunchCrossId = var; }
      void setBunchCrossId7bit(UInt_t var) { mBunchCrossId7bit = var; }
      void setNVertecies(UInt_t var) { mNVertecies = var; }
      void setNGoodTpcTrks(UInt_t var) { mNGoodTpcTrks = var; }
      void setTofMult(UShort_t var) { mTofMult = var; }
      void setBbcSmallEast(UInt_t var) { mBbcSmallEast = var; }
      void setBbcSmallWest(UInt_t var) { mBbcSmallWest = var; }
      void setBbcLargeEast(UInt_t var) { mBbcLargeEast = var; }
      void setBbcLargeWest(UInt_t var) { mBbcLargeWest = var; }
      void setZdcAdcEastPmt(UInt_t var, UInt_t pmt) { mZdcAdcEastPmt[pmt] = var; } // pmt must be: 0, 1 or 2
      void setZdcAdcWestPmt(UInt_t var, UInt_t pmt) { mZdcAdcWestPmt[pmt] = var; } // pmt must be: 0, 1 or 2
      void setZdcEastRate(Double_t var) { mZdcEastRate = var; }
      void setZdcWestRate(Double_t var) { mZdcWestRate = var; }
      void setZdcUnAttEast(UShort_t var) { mZdcEastUA = var; }
      void setZdcUnAttWest(UShort_t var) { mZdcWestUA = var; }
      void setZdcTdcEast(UShort_t var) { mZdcTdcEast = var; }
      void setZdcTdcWest(UShort_t var) { mZdcTdcWest = var; }
      void setZdcTimeDiff(UShort_t var) { mZdcTimeDiff = var; }
      void setZdcVertexZ(Double_t var) { mZdcVertexZ = var; }
      void setVertexZInCm(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mVertexZInCm[s][dataTag] = var; }
      void setVertexYInCm(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mVertexYInCm[s][dataTag] = var; }
      void setVertexXInCm(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mVertexXInCm[s][dataTag] = var; }
      void setDcaDaughters(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mDcaDaughtersInCm[s][dataTag] = var; }
      void setDcaBeamline(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mDcaBeamlineInCm[s][dataTag] = var; }
      void setPointingAngle(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mPointingAngle[s][dataTag] = var; }
      void setDecayLength(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mDecayLengthInCm[s][dataTag] = var; }
      void setVertexDiff(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mVertexDiffInCm[s][dataTag] = var; }      
      void setVertexIdTruth(UInt_t var) { mVertexIdTruth = var; }
      void setInvMass(Double_t var, unsigned int s) { mInvMass[s] = var; }
      void setTheta(Double_t var, unsigned int s) { mTheta[s] = var; }
      void setPhi(Double_t var, unsigned int s) { mPhi[s] = var; }
      void setP(Double_t var, unsigned int s) { mP[s] = var; }
      void setPt(Double_t var, unsigned int s) { mPt[s] = var; }
      void setRap(Double_t var, unsigned int s) { mRap[s] = var; }
      void setMSquared(Double_t var, unsigned int s) { mMSquared[s] = var; }
      void setPairID(Int_t var, unsigned int s) { mPairID[s] = var; }
      void setTotQ(Int_t var, unsigned int s) { mTotQ[s] = var; }
      void setPIDChiSquare(double var, unsigned int s) {mChiSquare[s] = var; }
      void setMomentumInGev(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mMomentumInGev[s][dataTag] = var; }   
      void setPtInGev(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mPtInGev[s][dataTag] = var; }  
      void setPxInGev(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mPxInGev[s][dataTag] = var; }
      void setPyInGev(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mPyInGev[s][dataTag] = var; }
      void setPzInGev(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mPzInGev[s][dataTag] = var; }
      void setEta(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mEtaHadrons[s][dataTag] = var; } 
      void setPhiAngle(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mPhiHadrons[s][dataTag] = var; } 
      void setCharge(Double_t var, unsigned int s, UInt_t dataTag = RECO) { mCharge[s][dataTag] = var; } 
      void setDEdxInKevCm(Double_t var, unsigned int s) { mDEdxInKevCm[s] = var; }
      void setTofHit(Double_t var, unsigned int s) { mTofHit[s] = var; }
      void setTofTimeInNs(Double_t var, unsigned int s) { mTofTimeInNs[s] = var; }   
      void setTofLengthInCm(Double_t var, unsigned int s) { mTofLengthInCm[s] = var; }  
      void setDcaXYInCm(Double_t var, unsigned int s) { mDcaXYInCm[s] = var; } 
      void setDcaZInCm(Double_t var, unsigned int s) { mDcaZInCm[s] = var; }   
      void setNHitsFit(Double_t var, unsigned int s) { mNHitsFit[s] = var; }   
      void setNHitsDEdx(Double_t var, unsigned int s) { mNHitsDEdx[s] = var; } 
      void setNSigmaTPC(Double_t var, unsigned int s, int p) { mNSigmaTPC[s][p] = var; }
      void setQATruth(UInt_t var, unsigned int s) { mQATruth[s] = var; }
      void setThetaRp(Double_t var, unsigned int s) { mThetaRp[s] = var; }
      void setPhiRp(Double_t var, unsigned int s) { mPhiRp[s] = var; }
      void setTimeRp(Double_t var, unsigned int s) { mTimeRp[s] = var; }
      void setT(Double_t var, unsigned int s) { mT[s] = var; }
      void setPRp(Double_t var, unsigned int s) { mPRp[s] = var; }
      void setPtRp(Double_t var, unsigned int s) { mPtRp[s] = var; }
      void setEtaRp(Double_t var, unsigned int s) { mEtaRp[s] = var; }
      void setRpX(Double_t var, unsigned int s) { mRpX[s] = var; }
      void setRpZ(Double_t var, unsigned int s) { mRpZ[s] = var; }
      void setRpY(Double_t var, unsigned int s) { mRpY[s] = var; }
      void setXi(Double_t var, unsigned int s) { mXi[s] = var; }
      void setPx(Double_t var, unsigned int s) { mPx[s] = var; }
      void setPy(Double_t var, unsigned int s) { mPy[s] = var; }
      void setPz(Double_t var, unsigned int s) { mPz[s] = var; }
      void setTrueRPVertexZInCm(Double_t var) { mTrueVertexZ = var; }
      void setTrueRPVertexXInCm(Double_t var) { mTrueVertexY = var; }
      void setTrueRPVertexYInCm(Double_t var) { mTrueVertexX = var; }
      void setRPTruePx(Double_t var, unsigned int s) { mTruePx[s] = var; }
      void setRPTruePy(Double_t var, unsigned int s) { mTruePy[s] = var; }
      void setRPTruePz(Double_t var, unsigned int s) { mTruePz[s] = var; }
      void setPtMissing(Double_t var) { mPtMissing = var; };
      void setPxMissing(Double_t var) { mPxMissing = var; };
      void setPyMissing(Double_t var) { mPyMissing = var; }; 
      void setTofTrigBit(bool var) { mTofTrigBit = var; }
      void setRpTrigBit(bool var) { mRpTrigBit = var; }
      void setRpItTrigBit(bool var) { mRpItTrigBit = var; }
      void setRpEtTrigBit(bool var) { mRpEtTrigBit = var; }
      void setBbcTrigBit(bool var) { mBbcTrigBit = var; }
      void setZdcTrigBit(bool var) { mZdcTrigBit = var; }
      void setZdcETrigBit(bool var) { mZdcETrigBit = var; }
      void setZdcWTrigBit(bool var) { mZdcWTrigBit = var; }
      void setTofDsmABit(bool var) { mTofDsmABit = var; }
      void setTofDsmBBit(bool var) { mTofDsmBBit = var; }
      void setTofDsmBit(bool var) { mTofDsmBit = var; }
      void setRpDsmBit(bool var) { mRpDsmBit = var; }
      void setRpItDsmBit(bool var) { mRpItDsmBit = var; }
      void setRpEtDsmBit(bool var) { mRpEtDsmBit = var; }
      void setBbcDsmBit(bool var) { mBbcDsmBit = var; }
      void setBbcSmallEDsmBit(bool var) { mBbcSmallEDsmBit = var; }
      void setBbcSmallWDsmBit(bool var) { mBbcSmallWDsmBit = var; }
      void setBbcLargeEDsmBit(bool var) { mBbcLargeEDsmBit = var; }
      void setBbcLargeWDsmBit(bool var) { mBbcLargeWDsmBit = var; }
      void setZdcDsmBit(bool var) { mZdcDsmBit = var; }
      void setZdcEDsmBit(bool var) { mZdcEDsmBit = var; }
      void setZdcWDsmBit(bool var) { mZdcWDsmBit = var; }
      void setRpTrigBits(bool var, unsigned int rp) { mRpTrigBits[rp] = var; }
      void setBemcTrackE(Double_t var, unsigned int s) { mBemcE[s] = var; }
      void setBemcTrackEta(Double_t var, unsigned int s) { mBemcEta[s] = var; }
      void setBemcTrackPhi(Double_t var, unsigned int s) { mBemcPhi[s] = var; }
      void setBemcTrackPt(Double_t var, unsigned int s) { mBemcPt[s] = var; }
      void setBemcClusterE(Double_t var, unsigned int s) { mBemcClusterE[s] = var; }
      void setBemcClusterEta(Double_t var, unsigned int s) { mBemcClusterEta[s] = var; }
      void setBemcClusterPhi(Double_t var, unsigned int s) { mBemcClusterPhi[s] = var; }
      void setBemcClusterHTE(Double_t var, unsigned int s) { mBemcClusterHTE[s] = var; }
      void setBemcClusterSigmaEta(Double_t var, unsigned int s) { mBemcClusterSigmaEta[s] = var; }
      void setBemcClusterSigmaPhi(Double_t var, unsigned int s) { mBemcClusterSigmaPhi[s] = var; }

      void setVertexStudyPrimary(bool var) { mPrimary = var; }
      void setVertexStudySameVertex(bool var) { mSameVertex = var; }
      void setVertexStudyDcaParticles(Double_t var) { mDcaParticles = var; }
      void setVertexStudyDcaBeamline(Double_t var) { mDcaBeamline = var; }

      // setters for good run list study
      void setRpOk(Int_t var) {mIsRpOk = var; }
      void setJPsiTrigger1(Int_t var) { mIsJPsiTrigger1 = var; }
      void setJPsiTrigger2(Int_t var) { mIsJPsiTrigger2 = var; }
      void setJPsiTrigger3(Int_t var) { mIsJPsiTrigger3 = var; }
      void setTpcTrackPhi(vector<Double_t> var) { mTpcTrack_phi = var; }
      void setTpcTrackEta(vector<Double_t> var) { mTpcTrack_eta = var; }
      void setBemcTrackPhi(vector<Double_t> var) { mBemcTrack_phi = var; }
      void setBemcTrackEta(vector<Double_t> var) { mBemcTrack_eta = var; }
      void setTpcNHitsFit(vector<Double_t> var) { mTpcNHitsFit = var; }
      void setTpcNHitsDEdx(vector<Double_t> var) { mTpcNHitsDEdx = var; }
      void setTpcNSigmaElectron(vector<Double_t> var) { mTpcNSigmaElectron = var; }
      void setTpcNSigmaProton(vector<Double_t> var) { mTpcNSigmaProton = var; }
      void setTpcNSigmaKaon(vector<Double_t> var) { mTpcNSigmaKaon = var;}
      void setTpcNSigmaPion(vector<Double_t> var) { mTpcNSigmaPion = var; }
      void setNTracksBemc(Int_t var) { nTracksBemc = var; }
      void setNTracksTof(Int_t var) { nTracksTof = var; }
      void setNClustersBemc(Int_t var) { nClustersBemc = var; }

      //void CalculatePID();

   private:
      Util* mUtil;
      TTree *mRecTree, *mBcgTree;
      // event info
      Double_t mRunNumber;
      UInt_t mEventNumber, mFillNumber, mBunchCrossId, mBunchCrossId7bit;
      UShort_t mTofMult;
      UInt_t mNVertecies, mNGoodTpcTrks;
      
      // BBC and ZDC info
      UInt_t mBbcSmallEast, mBbcSmallWest, mBbcLargeEast, mBbcLargeWest;
      UInt_t mZdcAdcEastPmt[3], mZdcAdcWestPmt[3];
      Double_t mZdcEastRate, mZdcWestRate;
      UShort_t mZdcEastUA, mZdcWestUA;
      UShort_t mZdcTdcEast, mZdcTdcWest, mZdcTimeDiff;
      Double_t mZdcVertexZ;

      // Vertex info
      Double_t mVertexZInCm[nStates][nDataTag], mVertexXInCm[nStates][nDataTag], mVertexYInCm[nStates][nDataTag];
      UInt_t mVertexIdTruth; // for MC vertices
      Double_t mVertexDiffInCm[nStates][nDataTag], mDcaDaughtersInCm[nStates][nDataTag], mDcaBeamlineInCm[nStates][nDataTag], mPointingAngle[nStates][nDataTag], mDecayLengthInCm[nStates][nDataTag];
      
      // State info
      Double_t mInvMass[nStates], mMSquared[nStates];
      Double_t mTheta[nStates], mPhi[nStates], mP[nStates], mPt[nStates], mRap[nStates];
      Int_t mPairID[nStates], mTotQ[nStates];

      // State PID chi square info
      Double_t mChiSquare[nParticles];

      // Central hadrons info
      Double_t mEtaHadrons[nHadrons][nDataTag], mPhiHadrons[nHadrons][nDataTag], mPtInGev[nHadrons][nDataTag];
      Double_t mCharge[nHadrons][nDataTag], mMomentumInGev[nHadrons][nDataTag];
      Double_t mDEdxInKevCm[nHadrons], mTofTimeInNs[nHadrons], mTofLengthInCm[nHadrons];
      Double_t mPxInGev[nHadrons][nDataTag], mPyInGev[nHadrons][nDataTag], mPzInGev[nHadrons][nDataTag];
      Double_t mNSigmaTPC[nHadrons][nParticles]; 
      Double_t mDcaXYInCm[nHadrons], mDcaZInCm[nHadrons], mNHitsFit[nHadrons], mNHitsDEdx[nHadrons];
      UInt_t mQATruth[nHadrons]; // for true MC hadrons
      Int_t mTofHit[nHadrons];

      // RP track info
      Double_t mThetaRp[nSides], mPhiRp[nSides], mTimeRp[nSides], mT[nSides];
      Double_t mPRp[nSides], mPtRp[nSides], mEtaRp[nSides], mRpX[nSides], mRpZ[nSides], mRpY[nSides];
      Double_t mXi[nSides], mPx[nSides], mPy[nSides], mPz[nSides];
      // RP MC info
      Double_t mTruePx[nSides], mTruePy[nSides], mTruePz[nSides];
      Double_t mTrueVertexX, mTrueVertexY, mTrueVertexZ;

      // pT missing info
      Double_t mPtMissing, mPxMissing, mPyMissing;

      // triger-bits info
      bool mTofTrigBit, mBbcTrigBit, mZdcTrigBit, mZdcETrigBit, mZdcWTrigBit, mBbcDsmBit, mZdcDsmBit, mRpDsmBit, mRpTrigBit;  
      bool mBbcSmallEDsmBit, mBbcSmallWDsmBit, mBbcLargeEDsmBit, mBbcLargeWDsmBit, mZdcEDsmBit, mZdcWDsmBit;
      bool mTofDsmBit, mTofDsmABit, mTofDsmBBit;
      bool mRpEtTrigBit, mRpItTrigBit, mRpEtDsmBit, mRpItDsmBit;
      bool mRpTrigBits[nRomanPots];

      // vertex eff part
      Double_t mDcaParticles, mDcaBeamline;
      bool mPrimary, mSameVertex;

      // BEMC info
      Double_t mBemcEta[nSigns], mBemcPhi[nSigns], mBemcPt[nSigns], mBemcE[nSigns];
      Double_t mBemcClusterEta[nSigns], mBemcClusterPhi[nSigns], mBemcClusterE[nSigns];
      Double_t mBemcClusterHTE[nSigns], mBemcClusterSigmaEta[nSigns], mBemcClusterSigmaPhi[nSigns];
      Int_t  nTracksBemc, nClustersBemc;

      // info about good run list
      Int_t mIsRpOk , mIsJPsiTrigger1, mIsJPsiTrigger2, mIsJPsiTrigger3, nTracksTof;
      vector<Double_t> mTpcTrack_phi, mTpcTrack_eta, mBemcTrack_eta, mBemcTrack_phi, mTpcNHitsFit, mTpcNHitsDEdx, mTpcNSigmaElectron, mTpcNSigmaPion, mTpcNSigmaProton, mTpcNSigmaKaon;



};

#endif

      

