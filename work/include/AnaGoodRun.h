#ifndef AnaGoodRun_h
#define AnaGoodRun_h

#include "Util.h"
#include "RecTree.h"
#include "Ana.h"
#include "UpcDstLibraries.h"
#include "StUPCTrack.h" 
#include "TLorentzVector.h"


using namespace std;
using namespace UTIL;

class AnaGoodRun : public Ana {
   public:
      AnaGoodRun(TFile *outfile);
      ~AnaGoodRun(); 

      void Init() override;
      void Make() override;

      void fillRunTree();
      bool areRPsCloseEnough(int mRunNumber);
      
      bool goodQualityTrack(const StUPCTrack *trk);
      bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]);
      Double_t getAverageBemcTracks() const { return (1.0*nTracksBEMC/nEvents); }
      Double_t getAverageBemcClusters() const { return (1.0*nClustersBEMC/nEvents); }
      Double_t getAverageTpcTracks() const {  return (1.0*nTracksTPC/nEvents); }
      Double_t getAverageTofTracks() const { return (1.0*nTracksTOF/nEvents); }
      Double_t getAverageVertices() const { return (1.0*nVertices/nEvents); }
      Double_t getAverageTpcEta() const { return (tpcEtaSum/nTracksTPC); }
      Double_t getAverageBemcEta() const { return (bemcEtaSum/nTracksBEMC); }
      Double_t getAverageTpcPhi() const { return (tpcPhiSum/nTracksTPC); }
      Double_t getAverageBemcPhi() const { return (bemcPhiSum/nTracksBEMC); }
      int getJPsiTriggerEvents() const { return nEventsJPsi; }
      int getNEventsZBVetoAll() const { return nEventsZBVetoAll; }
      int getNEventsZBVetoPassed() const { return nEventsZBVetoPassed; }
      int getNEventsAll() const { return nEventsAll; }
      int getNEventsLumiFile(int runnum) const {return mNEventsLumiFile.count(runnum) ? mNEventsLumiFile.at(runnum) : 0; }
      bool isJPsiTrigger() const { return JPsiTrigger; }
      double getInstantLuminosity(int RUNNUMBER) const {return mInstLumiPerRun.count(RUNNUMBER) ? mInstLumiPerRun.at(RUNNUMBER) : 0; }
      double getLuminosity(int RUNNUMBER) const {return mLumiPerRun.count(RUNNUMBER) ? mLumiPerRun.at(RUNNUMBER) : 0; }

   private:
   
      void addAverageBemcTracks(Double_t var) { nTracksBEMC += var; }
      void addAverageBemcClusters(Double_t var) { nClustersBEMC += var; }  
      void addAverageTpcTracks(Double_t var) { nTracksTPC += var; }
      void addAverageTOFTracks(Double_t var) { nTracksTOF += var; }
      void addAverageVertices(Double_t var) { nVertices += var; }
      void addTpcEtaValue(Double_t var)  { tpcEtaSum += var; }
      void addBemcEtaValue(Double_t var) { bemcEtaSum += var; }
      void addTpcPhiValue(Double_t var)  { tpcPhiSum += var; }
      void addBemcPhiValue(Double_t var) { bemcPhiSum += var; }
      void addJPsiTriggerEvent() { nEventsJPsi += 1;}
      void addEvent() { nEvents += 1; }
      void readLumiFile(bool isZB);
      void triggerVetoEfficiency();
      void addZBEventAll() { nEventsZBVetoAll += 1; }



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

      int nTracksBEMC, nClustersBEMC, nTracksTPC, nTracksTOF, nEvents, nVertices;
      Double_t tpcEtaSum, bemcEtaSum;
      Double_t tpcPhiSum, bemcPhiSum;
      int nEventsZBVetoAll, nEventsZBVetoPassed;  //all events for zerobias trigger,  ZB events that pass condition for JPsi 
      int nEventsAll, nEventsJPsi;  // all events in the run, events that have JPsi trigger
      bool JPsiTrigger;

};

#endif
