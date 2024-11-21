#ifndef TofEffMult_h
#define TofEffMult_h

#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class TofEffMult : public Ana{
   public:
      TofEffMult(TFile *outfile);
      ~TofEffMult(); 

      void Init() override;
      void Make() override;

   private:
   
      int tpcCounter;

      double pTMiss, invMass, deltaPhiProtons, pTState, phiState, thetaState;
      double probePt, probeEta, probePhi;
      double tagPt, tagEta, tagPhi;
      bool probeTofHit, probeTof;

      StUPCV0* getV0(int itrk1,const StUPCTrack* trk1,int itrk2, const StUPCTrack* trk2);
      void fillTopologyCutsBefore(const StUPCV0& V0);
      void fillTopologyCutsAfter(const StUPCV0& V0);
      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillNSigmaPlots(const StUPCTrack *trk);
      void fillEtaVtxPlots(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillBeamlineInfo();


      double bField;
      double beamline[4]; 


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

      TH1D *hInvMassTof1, *hInvMassTof2;



      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

};

#endif
