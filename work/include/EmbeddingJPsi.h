#ifndef EmbeddingJPsi_h
#define EmbeddingJPsi_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibraries.h"
#include "RunDef.h"
#include "StUPCV0.h"
#include "TParticle.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

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
      void fillBemcInfo(const StUPCTrack *trk);
      void fillBemcInfoAll();
      void trueMCPeak();
      void loadCuts();
      bool isBemcHit(const StUPCTrack *trk); // MC simulator of BEMC efficiency
      double bemcEfficiency(double pT); //function which returns bemc efficiency based on fit to real data
      void runControlOfBemcEfficiency();
      // all control histograms
      TH2F *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPPicorr, *hNSigmaKPcorr, *hNSigmaKPicorr;
      TH1D *hBranchRP;
      TH1D *hDcaZ, *hDcaZCut, *hDcaXY, *hDcaXYCut, *hNfitHits, *hNfitHitsCut, *hNhitsDEdx, *hNhitsDEdxCut, *hNVertices, *hTotQ;
      TH1D *hSameTrackPair, *hNTracksTpc, *hNTracksTof, *hNTracksBEMC, *hNClustersBEMC;
      TH1D *hNSigmaPi, *hNSigmaP, *hNSigmaK, *hDEdxSignal;
      TH1D *hPIDChiee, *hPIDChipp, *hPIDChipipi, *hPIDChikk;
      TH1D *hPt, *hPtCut;
      TH1D *hInvMassJPsi, *hTrackQualityFlow;
      Util* mUtil;
      TH2F *hNSigmaEE1, *hNSigmaEE2, *hNSigmaPP1, *hNSigmaPP2, *hNSigmaKK1, *hNSigmaKK2, *hNSigmaPiPi1, *hNSigmaPiPi2;

      TH1D *hInvMassJPsiMC;
      TH1D *hPtPair;

      
      
      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut, *hEtaBemc, *hEtaBemcCut;
      TH2F* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 

      TH1D *hEClusters, *hBemcTrackPhi, *hBemcClusterPhi, *hEClustersAll, *hBemcClusterPhiAll;
      TH2D *hEClusterTrack, *hBemcEtaPhiTrack, *hBemcEtaPhiCluster, *hBemcEtaPhiClusterAll;
      TH1D *hClusterMatched;
      
      Double_t beamline[4];
      Double_t bField;
      
      TGraph *gBEControl; // graph for control of bemc efficiency function
      
      vector<int> tracksBEMC;

      int tpcCounter;
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

      // cuts
      double MAXETA, MINPIDCHIPP, MINPIDCHIPIPI, MINPIDCHIKK, MAXPIDCHIEE;
      int MINNHITSFIT, MINNHITSDEDX;


};

#endif
