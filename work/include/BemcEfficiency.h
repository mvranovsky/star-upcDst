#ifndef BemcEfficiency_h
#define BemcEfficiency_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibraries.h"
#include "StPicoHelix.h"


using namespace std;
using namespace UTIL;

class BemcEfficiency : public Ana{
   public:
      BemcEfficiency(TFile *outfile);
      ~BemcEfficiency(); 

      void Init() override;
      void Make() override;

   private:
   
      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillTrackQualityCutsAfter(const StUPCTrack* trk);
      void fillNSigmaPlots(const StUPCTrack *trk);
      bool goodQualityTrack(const StUPCTrack *trk);
      bool sameVertex(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool chiSquarePID(const StUPCTrack *trk1, const StUPCTrack *trk2);
      bool exactly1RPTrack(int &side);
      bool fiducialVolume(const StUPCRpsTrack *trackRP, int side);
      void fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      bool isProjectedToBemc(const StUPCTrack *trk);
      void loadCuts();

      // all control histograms
      TH2F *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH2F *hNSigmaEE1, *hNSigmaEE2, *hNSigmaPP1, *hNSigmaPP2, *hNSigmaKK1, *hNSigmaKK2, *hNSigmaPiPi1, *hNSigmaPiPi2;
      TH2D *hInvMassEta, *hInvMassBemcEta;
      TH1D *hBranchRP;
      TH1D *hDcaZ, *hDcaZCut, *hDcaXY, *hDcaXYCut, *hNfitHits, *hNfitHitsCut, *hNhitsDEdx, *hNhitsDEdxCut, *hNVertices, *hTotQ;
      TH1D *hSameTrackPair, *hNTracksTpc, *hNTracksTOF, *hNTracksRP;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;
      TH1D *hPIDChiee, *hPIDChipp, *hPIDChipipi, *hPIDChikk;
      TH1D *hPt, *hPtCut;
      TH1D *hTrackQualityFlow;
      Util* mUtil;
      TH1D* hVtxZByFillNum;

      TH1D* hDeltaDipAngle, *hDeltaDipAngleCut, *hInvMass, *hInvMassCut;

      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut, *hEtaBemc, *hEtaBemcCut;
      TH2D* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 

      // virtual photon plots
      TH1D* hPtMissing, *hPtMissingCut, *hPhotonMomX, *hPhotonMomY, *hPhotonMomXBcg, *hPhotonMomYBcg, *hPtMissingBcg, *hPtMissingBcgCut;


      vector<int> tracksTOF;

      int tpcCounter;
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];


      // cuts
      double VERTEXZRANGE,MAXDCAZ ,MINDCAZ, MAXDCAXY, MAXETA, MINPIDCHIPP, MINPIDCHIPIPI, MINPIDCHIKK, MAXPIDCHIEE;
      int MINNHITSFIT, MINNHITSDEDX;
      double MAXDELTADIPANGLE, MAXINVMASS;

};

#endif
