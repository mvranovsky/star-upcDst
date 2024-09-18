#ifndef RpMCAna_h
#define RpMCAna_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "ElasticAna.h"

using namespace std;
using namespace UTIL;

class RpMCAna : public Ana {
   public:
      RpMCAna(TFile *outFile);
      ~RpMCAna();
      
      void Make() override;
      void Init() override;

      void SetRpPosition(TVector3 (&corr)[nRomanPots], TVector3 (&offsets)[nRomanPots]);
      void SetMCInfo(double (&mc_vrtx)[nCoordinates], double (&mc_p)[nCoordinates][nSides]);
      bool IsTrackInRp(int side);
      void ProjectToRP(double& x, double& y, int RP, int side);
      bool AreTracksInRp();

   private:

      void runEmbedding(bool embedding);
      void runMCEff(int set);

      TVector3 mRPpTBalance;
      
      EmbedMaker* embedMaker;
      StRPEvent *mEmbedEvt;
      ElasticAna *mElAna[DATA];

      double *mVertex, (*mMomentum)[nSides];
      TVector3 *mCorrection;
      TVector3 *mOffSet;

      // control plots
      enum { kAll = 1, kInside, kTrigger, kOneTrack, kDist, kMax};
      TH1D *hMCEffFlow[DATA]; // DATA = nSets { MC, MCZB}
      static const TString mMCEffFlowCutsName[kMax]; 

      // RP efficiency
      TH1D *hMCtoRecoDiff[DATA]; // DATA = nSets { MC, MCZB}
      TH2D *hEffPxPy[3][nBranches][DATA];
      TH2D *hEffXY[3][nRomanPots][DATA];
      TH2D *hEffXYOffSub[3][nRomanPots][DATA];
      TH2D *hNTpPerRP[DATA], *hNClustersPerPlane[DATA];
};

#endif
