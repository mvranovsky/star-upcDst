#ifndef AnaBP_h
#define AnaBP_h

// include headers
#include "UpcDstLibraries.h"
#include "Util.h"
#include "RecTree.h"
#include "Ana.h"

using namespace std;
using namespace UTIL;

class AnaBP : public Ana{
   public:
      AnaBP(TFile *outFile);
      ~AnaBP(){};
      
      void Make() override;
      void Init() override;

   private:
      UInt_t mRunNumber;

      TVector3 mRPpTBalance;
      
      // control plots
      TH1D *hTriggerBits;
      TH1D *hPIDStats[2];
      Util *mUtil;


      TH2F *hEtaPhi, *hRPcorr[2], *hRPcorrWest[2], *hRPcorrEast[2], *hNSigmaPiPcorr, *hNSigmaPiKcorr, *hNSigmaPiecorr, *hNSigmaPKcorr, *hNSigmaPecorr, *hNSigmaKecorr, *hNSigmaPUPCrr, *hNSigmaKPcorr, *hNSigmaKUPCrr;
      TH1D *hPt, *hEta, *hDcaZ, *hDcaXY, *hNfitHits, *hPosZ, *hNhitsDEdx, *hTOFTracks, *hNVertices, *hTotQ, *hNumberRPTracks, *hGlobalTracks;

      //TH1D *hMSquared[3][3][4];

      // RP ADC and TAC
      //TH1D *hRpAdc[2*nRomanPots];
      //TH1D *hRpAdcInWindow[2*nRomanPots];
      //TH1D *hRpTac[2*nRomanPots];

      void saveNSigmaCorr(const StUPCTrack* trk);
      void fillGoodTrackCuts(const StUPCTrack* trk);
      int hasGoodTPCnSigma(const StUPCTrack *trk);

      bool is2pions(const StUPCTrack *trk1, const StUPCTrack *trk2);
      bool isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2);

      void SaveMissingMomenta(TVector3 missP);
      bool oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2);
      bool lambda(const StUPCTrack *trk1, const StUPCTrack *trk2);
      /*
      void CalculatePID();
      Int_t PIDStartegy(int strategy);
      bool IsPairOf(int type);

      void CalculateTOFEff(unsigned int tagID);
      void FillMSquared();
      */

};

#endif 