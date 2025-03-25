#ifndef EmbeddingJPsi_h
#define EmbeddingJPsi_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibreries.h"
#include "RunDef.h"


using namespace std;
using namespace UTIL;

class EmbeddingJPsi : public Ana{
   public:
      EmbeddingJPsi(TFile *outfile);
      ~EmbeddingJPsi(); 

      void Init() override;
      void Make() override;

   private:
   
      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillTrackQualityCutsAfter(const StUPCTrack* trk);
      void fillNSigmaPlots(const StUPCTrack *trk);
      bool goodQualityTrack(const StUPCTrack *trk);
      bool sameVertex(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2);
      // all control histograms
      TH2F *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH1D *hBranchRP;
      TH1D *hDcaZ, *hDcaZCut, *hDcaXY, *hDcaXYCut, *hNfitHits, *hNfitHitsCut, *hNhitsDEdx, *hNhitsDEdxCut, *hNVertices, *hTotQ;
      TH1D *hSameTrackPair, *hNTracksTpc, *hNTracksTof, *hNTracksRP, *hNTracksBEMC, *hNClustersBEMC;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;
      TH1D *hPIDChiee, *hPIDChipp, *hPIDChipipi, *hPIDChikk;
      TH2F *hPIDChiep, *hPIDChiek, *hPIDChiepi, *hPIDChipip;
      TH1D *hPt, *hPtCut;
      TH1D *hInvMassJPsi, *hInvMassJPsiBcg, *hTrackQualityFlow;
      Util* mUtil;

      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut, *hEtaBemc, *hEtaBemcCut;
      TH2F* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 


      vector<int> tracksBEMC;

      int tpcCounter;
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

};

#endif
