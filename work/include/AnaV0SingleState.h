#ifndef AnaV0SingleState_h
#define AnaV0SingleState_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "StUPCSelectV0Modified.h"
#include "StUPCV0.h"
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>


using namespace std;
using namespace UTIL;

class AnaV0SingleState : public Ana{
   public:
      AnaV0SingleState(TFile *outFile);
      ~AnaV0SingleState(){};
      
      void Make() override;
      void Init() override;

   private:
      UInt_t mRunNumber;

      TVector3 mRPpTBalance;
      
      // control plots
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];


      TH2F *hEtaPhi, *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH1D *hPt, *hEta, *hDcaZ, *hDcaXY, *hNfitHits, *hPosZ, *hNhitsDEdx, *hTOFTracks, *hNVertices, *hTotQ, *hNumberRPTracks, *hGlobalTracks;
      TH1D *hPosZCut, *hVtxDiff, *hNPairV0, *hSameTrackPair, *hBothFlags;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;

      StUPCSelectV0Modified *mSelectV0; // selector for tracks from V0 candidates
      Util* mUtil;
      //topology cuts
      TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
      TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;
      TH2F *hDecayLPointingA, *hDecayLPointingACut, *hArmenterosPodolanski, *hInvMassEta;

      TH1D *hEta1Tag, *hEta1Probe, *hpT1Tag, *hpT1Probe, *hPhi1Tag, *hPhi1Probe;
      TH2F *hEtaPhi1Tag, *hEtaPhi1Probe;
      TH1D *hEta2Tag, *hEta2Probe, *hpT2Tag, *hpT2Probe, *hPhi2Tag, *hPhi2Probe;
      TH2F *hEtaPhi2Tag, *hEtaPhi2Probe;

      TH2F *hPtPhi1, *hPtVz1, *hVzPhi1, *hInvMassPhi1, *hInvMassVz1, *hInvMassPt1; 
      TH2F *hPtPhi2, *hPtVz2, *hVzPhi2, *hInvMassPhi2, *hInvMassVz2, *hInvMassPt2; 
      TH2F *hEtaPhiProbeYesToF, *hEtaPhiProbeNoToF;
      TH2F *hEtaPhiProbeYesToFLarge, *hEtaPhiProbeNoToFLarge;

      // plots around eta = 0
      TH1D *hNSigmaPiProbe1, *hNSigmaPProbe1, *hNSigmaKProbe1, *hNSigmaEProbe1; 
      TH1D *hNSigmaPiProbe2, *hNSigmaPProbe2, *hNSigmaKProbe2, *hNSigmaEProbe2; 
      TH1D *hNSigmaPiTag1, *hNSigmaPTag1, *hNSigmaKTag1, *hNSigmaETag1; 
      TH1D *hNSigmaPiTag2, *hNSigmaPTag2, *hNSigmaKTag2, *hNSigmaETag2; 



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

      bool is2pions(const StUPCTrack *trk1, const StUPCTrack *trk2);
      bool isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2);
      void fillBeamlineInfo(const StUPCEvent *mUpcEvt);
      //TLorentzVector get4Momentum(const StUPCTrack *trk1, const StUPCTrack *trk);


      void SaveMissingMomenta(TVector3 missP);
      bool oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool lambda(const StUPCTrack *trk1, const StUPCTrack *trk2);
      void fillTag1(const StUPCTrack* trk, Double_t pT, Double_t Vz, Double_t invMass);
      void fillProbe1(const StUPCTrack* trk, Double_t pT, Double_t Vz, Double_t invMass);
      void fillTag2(const StUPCTrack* trk, Double_t pT, Double_t Vz, Double_t invMass);
      void fillProbe2(const StUPCTrack* trk, Double_t pT, Double_t Vz, Double_t invMass);


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
