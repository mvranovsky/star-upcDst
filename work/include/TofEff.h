#ifndef TofEff_h
#define TofEff_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibreries.h"


using namespace std;
using namespace UTIL;

class TofEff : public Ana{
   public:
      TofEff(TFile *outfile);
      ~TofEff(); 

      void Init() override;
      void Make() override;

   private:
   
      void CalculateTOFEff(unsigned int tagID);
      StUPCV0* getV0(int itrk1,const StUPCTrack* trk1,int itrk2, const StUPCTrack* trk2);
      void fillTopologyCutsBefore(const StUPCV0& V0);
      void fillTopologyCutsAfter(const StUPCV0& V0);
      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillNSigmaPlots(const StUPCTrack *trk);
      void fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);



      // all control histograms
      TH2F *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH1D *hPt, *hDcaZ, *hDcaXY, *hNfitHits, *hNhitsDEdx, *hHasPrimVtx, *hNVertices, *hTotQ, *hNumberRPTracks, *hGlobalTracks;
      TH1D *hVtxDiff, *hVtxDiffAfter, *hNPairV0, *hSameTrackPair, *hBothFlags, *hNTracksTpc, *hNTracksTof;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;
      Util* mUtil;
      //topology cuts
      TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength;
      TH1D *hDcaDaughtersCut, *hDcaBeamlineCut, *hPointingAngleCut, *hDecayLengthCut;
      TH2F *hDecayLPointingA, *hDecayLPointingACut, *hArmenterosPodolanski, *hInvMassEta;
      TH1D *hMultipleGoodV0;

      TH1D *hInvMassTof1, *hInvMassTof2;

      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut;
      TH2F* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 



      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

};

#endif
