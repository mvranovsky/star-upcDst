#ifndef AnaGoodRun_h
#define AnaGoodRun_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibreries.h"



using namespace std;
using namespace UTIL;

class AnaGoodRun : public Ana{
   public:
      AnaGoodRun(TFile *outfile);
      ~AnaGoodRun(); 

      void Init() override;
      void Make() override;

   private:
   
      bool goodQualityTrack(const StUPCTrack *trk);
      bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]);


      vector<int> tracksTPC;
      vector<int> tracksBEMC;

      TH1D *hTriggerBits;

      map<unsigned int, TVector3> corrections[nRomanPots];
      map<unsigned int, TVector3> offsets[nRomanPots];

};

#endif
