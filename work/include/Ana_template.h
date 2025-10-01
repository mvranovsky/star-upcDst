#ifndef Ana_template_h
#define Ana_template_h

// include headers
#include "UpcDstLibraries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class Ana_template : public Ana{
   public:
      Ana_template(TFile *outFile);
      ~Ana_template(){};
      
      void Make() override;
      void Init() override;

   private:
      UInt_t mRunNumber;

      TVector3 mRPpTBalance;
      
      // control plots
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];

      TH1D *hMSquared[3][3][4];

      // RP ADC and TAC
      TH1D *hRpAdc[2*nRomanPots];
      TH1D *hRpAdcInWindow[2*nRomanPots];
      TH1D *hRpTac[2*nRomanPots];

      void CalculatePID();
      Int_t PIDStartegy(int strategy);
      bool IsPairOf(int type);

      void CalculateTOFEff(unsigned int tagID);
      void SaveMissingMomenta(TVector3 missP);
      void FillMSquared();
      
      //TOF eff study
      TTree *tofEffTree;
      double pTMiss, invMass, deltaPhiProtons, pTState, phiState;
      double probePt, probeEta, probePhi;
      double tagPt, tagEta, tagPhi;
      bool probeTofHit, probeTof;
      unsigned int nProbes;
};

#endif 