#ifndef AnaV0Mult_h
#define AnaV0Mult_h

// include headers
#include "Util.h"
#include "UpcDstLibreries.h"
#include "RecTree.h"
#include "Ana.h"
#include "StUPCV0.h"
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>


using namespace std;
using namespace UTIL;

class AnaV0Mult : public Ana{
   public:
      AnaV0Mult(TFile *outFile);
      ~AnaV0Mult(){};
      
      void Make() override;
      void Init() override;

   private:
      UInt_t mRunNumber;

      TVector3 mRPpTBalance;
      
      // control plots
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut;
      TH2F* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 


      TH2F *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH1D *hPt, *hDcaZ, *hDcaXY, *hNfitHits, *hNhitsDEdx, *hTOFTracks, *hNVertices, *hTotQ, *hNumberRPTracks, *hGlobalTracks;
      TH1D *hVtxDiff, *hNPairV0, *hSameTrackPair, *hBothFlags;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;

      Util* mUtil;
      //topology cuts
      TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
      TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;
      TH2F *hDecayLPointingA, *hDecayLPointingACut, *hArmenterosPodolanski, *hInvMassEta;


      //TH1D *hMSquared[3][3][4];

      // RP ADC and TAC
      //TH1D *hRpAdc[2*nRomanPots];
      //TH1D *hRpAdcInWindow[2*nRomanPots];
      //TH1D *hRpTac[2*nRomanPots];
      //void removeDuplicates(vector<pair<int, int>>& vec, TVector3 vtx, double * beamline, float const bField); 
      //int findValuePosition(map<int, int> myMap, int value);
      vector<pair<int, int>> filterPairs(vector<pair<int, int>>& pairs);
      bool shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b);
      void saveNSigmaCorr(const StUPCTrack *trk);
      void fillGoodTrackCuts(const StUPCTrack* trk);
      void fillTopologyCutsBefore(const StUPCV0& V0);
      void fillTopologyCutsAfter(const StUPCV0& V0);
      int hasGoodTPCnSigma(const StUPCTrack *trk);
      void fillFinalPlots();
      void fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);


      bool is2pions(const StUPCTrack *trk1, const StUPCTrack *trk2);
      bool isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2);

      void SaveMissingMomenta(TVector3 missP);
      bool oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool lambda(const StUPCTrack *trk1, const StUPCTrack *trk2);


};

#endif 
