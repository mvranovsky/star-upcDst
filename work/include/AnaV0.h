#ifndef AnaV0_h
#define AnaV0_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "StUPCV0.h"
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>
#include <memory>



using namespace std;
using namespace UTIL;

class AnaV0 : public Ana{
   public:
      AnaV0(TFile *outFile);
      ~AnaV0(){};
      
      void Make() override;
      void Init() override;

   private:
      UInt_t mRunNumber;

      TVector3 mRPpTBalance;
      
      // control plots
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

      // all control histograms
      TH2F *hEtaPhi, *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH1D *hPt, *hEta, *hEtaCut, *hDcaZ, *hDcaXY, *hNfitHits, *hPosZ, *hNhitsDEdx, *hHasPrimVtx, *hNVertices, *hTotQ, *hNumberRPTracks, *hGlobalTracks;
      TH1D *hPosZCut, *hVtxDiff, *hVtxDiffAfter, *hNPairV0, *hSameTrackPair, *hBothFlags, *hNTracksTpc, *hNTracksTof;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;
      Util* mUtil;
      //topology cuts
      TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
      TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;
      TH2F *hDecayLPointingA, *hDecayLPointingACut, *hArmenterosPodolanski, *hInvMassEta;
      TH2F *hEtaVtxZ, *hEtaVtxZCut;
      TH1D *hMultipleGoodV0;

      TH1D* hInvMassTof1, *hInvMassTof2, *hInvMassTof1AfterPicking1, *hInvMassTof2AfterPicking1;



      //TH1D *hMSquared[3][3][4];

      // RP ADC and TAC
      //TH1D *hRpAdc[2*nRomanPots];
      //TH1D *hRpAdcInWindow[2*nRomanPots];
      //TH1D *hRpTac[2*nRomanPots];
      //void removeDuplicates(vector<pair<int, int>>& vec, TVector3 vtx, double * beamline, float const bField); 
      //int findValuePosition(map<int, int> myMap, int value);
      vector<pair<int, int>> filterPairs(vector<pair<int, int>>& pairs);
      bool shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b);
      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillTopologyCutsBefore(const StUPCV0& V0);
      void fillTopologyCutsAfter(const StUPCV0& V0);
      int hasGoodTPCnSigma(const StUPCTrack *trk);
      void fillFinalPlots();
      void fillPrimVtxInfo(const StUPCTrack* trk1, const StUPCTrack *trk2);
      pair<int,int> resizePairs(vector<pair<int,int>> hadronPairV0);
      void fillNSigmaPlots(const StUPCTrack *trk);
      void fillEtaVtxPlots(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillBeamlineInfo();




      bool is2pions(const StUPCTrack *trk1, const StUPCTrack *trk2);
      bool isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2);
      void fillBeamlineInfo(const StUPCEvent *mUpcEvt);
      //TLorentzVector get4Momentum(const StUPCTrack *trk1, const StUPCTrack *trk);


      void SaveMissingMomenta(TVector3 missP);
      bool oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool lambda(const StUPCTrack *trk1, const StUPCTrack *trk2);


      double bField;
      double beamline[4]; 

      /*
      void CalculatePID();
      Int_t PIDStartegy(int strategy);
      bool IsPairOf(int type);

      void CalculateTOFEff(unsigned int tagID);
      void FillMSquared();
      */



};

#endif 
