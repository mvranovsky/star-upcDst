#ifndef AnaJPSI_h
#define AnaJPSI_h

#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class AnaJPSI : public Ana{
   public:
      AnaJPSI(TFile *outfile);
      ~AnaJPSI(); 

      void Init() override;
      void Make() override;

   private:
   
      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillTrackQualityCutsAfter(const StUPCTrack* trk);
      void fillNSigmaPlots(const StUPCTrack *trk);
      bool goodQualityTrack(const StUPCTrack *trk);
      bool sameVertex(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2);
      void fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);


      // all control histograms
      TH2F *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH2D *hInvMassEta;
      TH1D *hDcaZ, *hDcaZCut, *hDcaXY, *hDcaXYCut, *hNfitHits, *hNfitHitsCut, *hNhitsDEdx, *hNhitsDEdxCut, *hNVertices, *hTotQ;
      TH1D *hSameTrackPair, *hNTracksTpc, *hNTracksTof, *hNTracksRP, *hNTracksBEMC;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;
      TH1D *hPt, *hPtCut;
      TH1D *hInvMassJPsi, *hInvMassJPsiBcg, *hTrackQualityFlow;
      TH1D *hEtaDifference;
      Util* mUtil;

      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut;
      TH2F* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 



      vector<int> tracksBEMC;

      int tpcCounter;
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

};

#endif
