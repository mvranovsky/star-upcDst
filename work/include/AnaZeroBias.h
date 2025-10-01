#ifndef AnaZeroBias_h
#define AnaZeroBias_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibraries.h"
#include "StUPCTrack.h" 
#include "TLorentzVector.h"


using namespace std;
using namespace UTIL;

class AnaZeroBias : public Ana {
   public:
      AnaZeroBias(TFile *outfile);
      ~AnaZeroBias(); 

      void Init() override;
      void Make() override;

      void fillRunTree();
      
      bool goodQualityTrack(const StUPCTrack *trk);
      bool isJPsiTrigger() const { return JPsiTrigger; }
      double getInstantLuminosity(int RUNNUMBER) const {return mInstLumiPerRun.count(RUNNUMBER) ? mInstLumiPerRun.at(RUNNUMBER) : 0; }
      double getLuminosity(int RUNNUMBER) const {return mLumiPerRun.count(RUNNUMBER) ? mLumiPerRun.at(RUNNUMBER) : 0; }

   private:
   
      bool goodPID(const StUPCTrack *trk1, const StUPCTrack *trk2);
      void addZBEventAll() { nEventsZBAll += 1; }

      TH1D *hNTracksTof, *hNTracksBemc, *hNVertices;

      vector<int> tracksTPC;
      vector<int> tracksBEMC;

      TH1D *hTriggerBits;
      TH2F *hTpcEtaPhi, *hBemcEtaPhi;
      TH1D *hTpcEta, *hBemcEta;
      TH1D *hTpcPhi, *hBemcPhi;

      map<unsigned int, TVector3> corrections[nRomanPots];
      map<unsigned int, TVector3> offsets[nRomanPots];

      map<int,double> mInstLumiPerRun;
      map<int,int> mNEventsLumiFile;
      map<int,double> mLumiPerRun;
      int nEventsZBAll;

      int nEventsAll, nEventsJPsi;  // all events in the run, events that have JPsi trigger
      bool JPsiTrigger;

};

#endif
