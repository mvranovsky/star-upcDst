#ifndef Util_h
#define Util_h

//_____________________________________________________________________________
//    Class to handle all useful utils
//    Author: Truhlar Tomas
//_____________________________________________________________________________

#include "Libreries.h"

namespace UTIL {
   class Util; // just tell the compiler to expect a class def

   // Enumerations - very helpful and convenient !
   enum STUDYMAP { kMAINANA = 0, kEMBEDQA, kVERTEXSTUDY, kTOFQA, kTRIGEFF, kFULLZB, kELASTICANA, kRPMCANA, kALIGNMENT, nStudies };
   enum SIDE { E=0, East=0, W=1, West=1, nSides };
   enum RPPORIENTATIONS { Up=0, Down=1, nRpOrientations };
   enum RPCONFIGURATION { IT=0, ET=1, nRpConfigurations};
   enum XY_COORDINATE { X = 0, Y, Z, nCoordinates, nXYCoordinates = Z};
   enum TRIGGER_ID { UPC_v1, UPC_v2, SDT, ET_v1, ET_v2, CPT2_v1, CPT2_v2, CPT2noBBCL, Zerobias, JPsiHTTP_v1, 
   JPsiHTTP_v2, JPsiHTTP_v3, SDT_RHICf, ET_RHICf, CPT2_RHICf, CPT2noBBCL_RHICf_v1, CPT2noBBCL_RHICf_v2, nTriggers };
   enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots };
   enum PLANE_ID {A, B, C, D, nPlanes };
   enum BRANCH_ID { EU, ED, WU, WD, nBranches };
   enum BRANCHES_CONFIGURATION_ID { CONF_EU_WU=0, CONF_ED_WD, CONF_EU_WD, CONF_ED_WU, nBranchesConfigurations };
   enum ARM_ID { EU_WD, ED_WU, nArms };
   enum STATION_ID { E1, E2, W1, W2, nStations };
   enum STATION_ORDER { RP1, RP2, nStationPerSide, nRpPerStation = nStationPerSide};
   enum PARTICLE { PION = 0, KAON, PROTON, ELECTRON,nParticles };
   enum NEUTRALPARTICLE { K0S = 0, LAMBDA, LAMBDABAR , nNeutralParticles};
   enum TPC_TRACK_TYPE { GLO, PRI, TOF, QUA, nTpcTrkTypes }; // GLO=global(all), PRI=primary, TOF=PRI&TofMatched, QUA=TOF&QualityCuts
   enum BUNCH_CROSSING { CB, AG, nBnchXngsTypes }; // CB=colliding bunches, AG=abort gaps
   enum QSUM_2TRKS { OPPO, SAME, nCharges2Trks };
   enum SIGN { PLUS, MINUS, nSigns };
   enum LIST_OF_EFF_CORRECTIONS { RPACC, TPCRECOEFF, TOFMATCHEFF, nEffCorrections };
   enum ANALYSIS_CUT { ALL = 1, TRIG, TWORPTRKS, INFID, ONETOFVX, ZVERTEX, TWOTOFTRKS, ETA, OPPOSITE, PIPI, PPI, PIPBAR, nAnalysisCuts };
   enum V0SELECTION_CUT {V0ALL = 1, V0TRIG,V0FLAG, V0PID, V0PAIR, V0ETAVTXZ, V0OPPOSITE, V0PIPI, V0PPI, V0PIPBAR , nV0SelectionCuts };
   enum JPSISELECTION_CUT {JPSIALL = 1, JPSITRIG, JPSI1VTX, JPSIVTXZ,JPSIBEMC , JPSIBACKTOBACK,JPSIPID, JPSI1RP, JPSIRPFIDCUT,JPSIQTOT, nJPSISelectionCuts };
   enum JPSISELECTION_CUT2 {JPSI2ALL = 1, JPSI2TRIG,JPSI2BEMC, JPSI2SAMEVTX,JPSI2VTXZ, JPSI2BACKTOBACK,JPSI2PID, JPSI21RP, JPSI2RPFIDCUT,JPSI2QTOT,nJPSI2SelectionCuts };
   enum EMBEDDING{ EMBEDDINGALL = 1, EMBEDDING2BEMC, EMBEDDINGBACKTOBACK, EMBEDDINGPID, EMBEDDINGQTOT, nEmbeddingCuts };
   enum GOODRUN_CUT {GRALL = 1, GRTRIGGER , GRRP, GRGOODTRACKTPC, GRGOODTRACKBEMC, nGRCuts };
   enum TOFEFF_CUT {TOFALL = 1, TOFTRIG, TOFTRACKQUALITY, TOFETAVTXZ, TOFPAIR, TOFOPPOSITE , nTOFEFFCuts};
   enum RANGE_LIMIT { MIN, MAX };
   enum DATASET { MC = 0, MCZB, DATA, nDataSets };
   enum DATATAG { TRUEMC = 0, RECO, nDataTag };
   enum NUMBEREDHADRONS {PLUS0 = 0,MINUS0, PLUS1, MINUS1, PLUS2, MINUS2, PLUS3, MINUS3, PLUS4, MINUS4,/*PLUS5,MINUS5, PLUS6, MINUS6, PLUS7, MINUS7, PLUS8, MINUS8,PLUS9, MINU9,*/ nHadrons};
   enum NUMBEREDSTATES {STATE0 = 0, STATE1, STATE2, STATE3, STATE4,/*STATE5,STATE6, STATE7,STATE8,STATE9,*/ nStates};
}

using namespace std;
using namespace UTIL;

class UTIL::Util{
   public:
      Util();
      ~Util();
      inline TString sideName(UInt_t id) const { if(id<nSides) return mSideName[id]; else{ std::cerr << "ERROR in Util::sideName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString coordinateName(UInt_t id) const { if(id<nCoordinates) return mCoordinateName[id]; else{ std::cerr << "ERROR in Util::coordinateName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString triggerName(UInt_t id) const { if(id<nTriggers) return mTriggerName[id]; else{ std::cerr << "ERROR in Util::triggerName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString armName(UInt_t id) const { if(id<nArms) return mArmName[id]; else{ std::cerr << "ERROR in Util::armName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString branchName(UInt_t id) const { if(id<nBranches) return mBranchName[id]; else{ std::cerr << "ERROR in Util::branchName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString branchesConfigurationName(UInt_t id) const { if(id<nBranchesConfigurations) return mBranchesConfigurationName[id]; else{ std::cerr << "ERROR in Util::branchesConfigurationName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString rpName(UInt_t id) const { if(id<nRomanPots) return mRpName[id]; else{ std::cerr << "ERROR in Util::rpName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString rpConfigName(UInt_t id) const { if(id<nRpConfigurations) return mRpConfigName[id]; else{ std::cerr << "ERROR in Util::rpConfigName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString planeName(UInt_t id) const { if(id<nPlanes) return mPlaneName[id]; else{ std::cerr << "ERROR in Util::planeName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString stationName(UInt_t id) const { if(id<nStations) return mStationName[id]; else{ std::cerr << "ERROR in Util::stationName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString particleName(UInt_t id) const { if(id<nParticles) return mParticleName[id]; else{ std::cerr << "ERROR in Util::particleName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString neutralParticleName(UInt_t id) const { if(id<nNeutralParticles) return mNeutralParticleName[id]; else{ std::cerr << "ERROR in Util::particleName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString tpcTrackTypeName(UInt_t id) const { if(id<nTpcTrkTypes) return mTpcTrackTypeName[id]; else{ std::cerr << "ERROR in Util::tpcTrackTypeName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString bunchXngTypeName(UInt_t id) const { if(id<nBnchXngsTypes) return mBunchCrossingTypeName[id]; else{ std::cerr << "ERROR in Util::bunchXngTypeName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString qSum2TrksName(UInt_t id) const { if(id<nCharges2Trks) return mChargeSum2TrksName[id]; else{ std::cerr << "ERROR in Util::qSum2TrksName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString signName(UInt_t id) const { if(id<nSigns) return mSignName[id]; else{ std::cerr << "ERROR in Util::signName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString efficiencyName(UInt_t id) const { if(id<nEffCorrections) return mEfficiencyName[id]; else{ std::cerr << "ERROR in Util::efficiencyName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisCutName(UInt_t id) const { if(id<nAnalysisCuts) return mCutName[id]; else{ std::cerr << "ERROR in Util::analysisCutName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisV0SelectionName(UInt_t id) const { if(id<nV0SelectionCuts) return mV0CutName[id]; else{ std::cerr << "ERROR in Util::AnaV0CutName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisJPSI(UInt_t id) const { if(id<nJPSISelectionCuts) return mJPSICutName[id]; else{ std::cerr << "ERROR in Util::JPsiCutName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisJPSI2(UInt_t id) const { if(id<nJPSI2SelectionCuts) return mJPSI2CutName[id]; else{ std::cerr << "ERROR in Util::JPSI2CutName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString embeddingName(UInt_t id) const { if(id<nEmbeddingCuts) return mEmbeddingName[id]; else{ std::cerr << "ERROR in Util::embeddingName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisGoodRun(UInt_t id) const { if(id<nGRCuts) return mGRCutName[id]; else{ std::cerr << "ERROR in Util::GRCutName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString analysisTofEff(UInt_t id) const { if(id<nTOFEFFCuts) return mTOFEFFName[id]; else{ std::cerr << "ERROR in Util::analysisTofEff(UInt_t id): id out of range" << std::endl; return TString("");} }      
      inline TString dataSetName(UInt_t id) const { if(id<nDataSets) return mDataSetName[id]; else{ std::cerr << "ERROR in Util::dataSetName(UInt_t id): id out of range" << std::endl; return TString("");} }
      inline TString dataTagName(UInt_t id) const { if(id<nDataTag) return mDataTagName[id]; else{ std::cerr << "ERROR in Util::dataTagName(UInt_t id): id out of range" << std::endl; return TString("");} }


      inline Double_t mass(int name) const{ return mParticleMass[name]; }
      inline Double_t c() const{ return mSpeedOfLight; }
      inline Double_t p0() const{ return mBeamMomentum; }
      inline Double_t pi() const{ return mPi; }
      inline Double_t epsiolon() const{ return mEpsilon; }

      inline Double_t rpZPosition(int rpId) const { return mRpZPosition[rpId]; };
      inline int branchPerRp(int rpId) const { return mBranchPerRp[rpId]; };
      inline int oppositeBranch(int br) const { return mOppositeBranch[br]; };
      inline int sidePerRp(int rpId) const { return mSidePerRp[rpId]; };
      inline int stationOrderPerRp(int rpId) const { return mStationOrderPerRp[rpId]; };
      inline int rpPerBranchStationOrder(int br, int st) const { return mRpPerBranchStationOrder[br][st]; }
      inline int branchPerBranchConfiguration(int brConf, int side) const { return mBranchPerBranchConfiguration[brConf][side]; }
      inline int planeToCoor(int plane) const { return plane%2 ? X : Y;} 
      inline int mapBranchConfiguration(int estOri, int wstOri) const { return mBranchConfigMap[estOri][wstOri];}
      
      Double_t binomialCoeff(UInt_t, UInt_t) const;
      Double_t bkgdFraction(const TH1*, Double_t=0.1, Double_t=0.16, Double_t=0.24, TF1* =nullptr) const;
      Double_t integratePtMiss(const TH1*, Double_t=0.1) const;
      vector<float> getBinsVectorF(const TAxis *) const;
      vector<double> getBinsVectorD(const TAxis *) const;
      TH1D* bkgdHistogram(const TH2*, Double_t=0.1, Double_t=0.16, Double_t=0.24, Int_t=0, vector<TF1*> * = nullptr) const;
      TH1D* backgroundFromSameSignTemplate( TH2D* const, Double_t, TH1D* const ) const;
      void subtractBackground(TH1*, const TH1*) const;

      TVector3 fitLine(const vector<TVector3> trackPoints, double positionOfFit);

   private:
      // Labels, names etc. (defined as TString to gain higher functionality than const char*, e.g. defined "+" operator)
      TString* mSideName;
      TString* mCoordinateName;
      TString* mTriggerName;
      TString* mArmName;
      TString* mBranchName;
      TString* mBranchesConfigurationName;
      TString* mRpName;
      TString* mRpConfigName;
      TString* mPlaneName;
      TString* mStationName;
      TString* mParticleName;
      TString* mNeutralParticleName;
      TString* mTpcTrackTypeName;
      TString* mBunchCrossingTypeName;
      TString* mChargeSum2TrksName;
      TString* mSignName;
      TString* mEfficiencyName;
      TString* mCutName;
      TString* mV0CutName;
      TString* mJPSICutName;
      TString* mJPSI2CutName;
      TString* mDataSetName;
      TString* mDataTagName;
      TString* mTOFEFFName;
      TString* mGRCutName;
      TString* mEmbeddingName;
          
      Double_t mParticleMass[nParticles]; // GeV/c^2
      const Double_t mSpeedOfLight; // m/s
      const Double_t mBeamMomentum; // GeV/c
      const Double_t mPi;
      const Double_t mEpsilon;
      
      Double_t mRpZPosition[nRomanPots];
      BRANCH_ID mBranchPerRp[nRomanPots];
      BRANCH_ID mOppositeBranch[nBranches];
      SIDE mSidePerRp[nRomanPots];
      STATION_ORDER mStationOrderPerRp[nRomanPots];
      RP_ID mRpPerBranchStationOrder[nBranches][nStationPerSide];
      BRANCH_ID mBranchPerBranchConfiguration[nBranchesConfigurations][nSides];   
      BRANCHES_CONFIGURATION_ID mBranchConfigMap[nRpOrientations][nRpOrientations]; 
};

   void SumDistance2(int &, double *, double & sum, double * par, int);
   void line(double z, double *p, double &x, double &y);
   double distance2(double x,double y,double z, double *p);

#endif
