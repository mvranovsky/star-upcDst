#ifndef ElasticAna_h
#define ElasticAna_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class ElasticAna : public Ana{
   public:
      ElasticAna(TFile *outfile);
      ~ElasticAna(); 

      void Init() override;
      void Make() override;

      TVector3 CalculateOffsetCorrForRP(unsigned int Rp, unsigned int iter);

      void InitRPMCInfo(); 
      void SetRPMCInfo(double *mc_vrtx, double (*mc_p)[nSides]);

   private:
      // Control plots
      TH1D *hElasticAnaFlow; 
      enum ANA_CUT { ALL = 1, TRIG, TWORPTRKS, INFID, COLLINEAR, nCuts };
      static const TString mCutName[nCuts];
      enum { kAll = 1, kTtrig, kETPattern, kFourPoints, kTRange,kColinear, kFiducial, kDCut, kElastic, kMax};
      static const TString mElAnaCutsName[kMax];

      TH1D *hTheta[2][2];
      TH1D *hDeltaTheta[2];
      TH2D *hThetaCorr;
      TH1D *hDCutR;
      TH2D *hDCut;

      TH2D *hVertexExtraction[nArms][nCoordinates -1];
      TH1D *hVertexZExtraction[nArms], *hDeltaProjection[nCoordinates-1][nBranches], *hDeltaFitProjection[nRomanPots][nCoordinates-1];
      TH1D *hClusterLength[nRomanPots], *hClusterEnergy[nRomanPots], *hClusterPerPlain[nRomanPots][nPlanes], *hNClusterPerRP[nRomanPots];
      TH1D *hT;

      // pX, pY and x,y RP accaptance
      TH2D *hEffPxPy[nBranches];
      TH2D *hEffXY[nRomanPots];

      // RP ADC and TAC
      TH1D *hRpAdc[2*nRomanPots];
      TH1D *hRpAdcInWindow[2*nRomanPots];
      TH1D *hRpTac[2*nRomanPots];

      // alignment calculation
      vector<double> alignment[nRomanPots][nXYCoordinates];
      TH1D *hAligEvents[nRomanPots][nAligIteration], *hAligXCorr[nRomanPots][nAligIteration], *hAligYCorr[nRomanPots][nAligIteration]; 

      UInt_t mRunNumber;
      
      typedef struct track_t{
        float posX[4];   // point position at RP 
        float posY[4];   // point position at RP 
        float posZ[4];   // point position at RP
        unsigned int nPoints; 
        float thX;
        float thY;
        float phi;
        float t;
        float X0;
        float Y0;  // X,Y track position at z=0
      } TRACK_T;

      bool IsInGeoWindow( double phi){ return ( abs(abs( phi )- PHI_CENT ) < PHI_WIDTH); }
      track_t makeTrack(bitset<8> rpsBits, int copyRP);
      bool IsElasticEvent();
      void FillClusterInfo();
      void runMCValidation();
      void runVertexExtraction();
      void closerTest(const vector<TVector3> trackPoints, int arm);
      void FillAccaptancePlots();
      void runAlignment();
};

#endif
