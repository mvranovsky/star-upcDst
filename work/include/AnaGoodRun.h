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
      
      bool goodQualityTrack(const StUPCTrack *trk);
      bool LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]);
      Double_t getAverageBemcTracks() const { return (1.0*nTracksBEMC/nEvents); }
      Double_t getAverageBemcClusters() const { return (1.0*nClustersBEMC/nEvents); }
      Double_t getAverageTpcTracks() const {  return (1.0*nTracksTPC/nEvents); }
      Double_t getAverageTOFTracks() const { return (1.0*nTracksTOF/nEvents); }
      Double_t getAverageVertices() const { return (1.0*nVertices/nEvents); }
      Double_t getAverageTpcEta() const { return (tpcEtaSum/nTracksTPC); }
      Double_t getAverageBemcEta() const { return (bemcEtaSum/nTracksBEMC); }
      Double_t getAverageTpcPhi() const { return (tpcPhiSum/nTracksTPC); }
      Double_t getAverageBemcPhi() const { return (bemcPhiSum/nTracksBEMC); }
      int getJPsiTriggerEvents() const { return nEventsJPsi; }
      int getNEventsAll() const { return nEventsAll; }
      int getNEventsPassed() const { return nEventsPassed; }
      bool isJPsiTrigger() const { return JPsiTrigger; }
      int getLuminosity(int RUNNUMBER) {return mInstLumiPerRun[RUNNUMBER];}

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
      void readLumiFile();

      vector<int> tracksTPC;
      vector<int> tracksBEMC;

      TH1D *hTriggerBits;
      TH2F *hTpcEtaPhi, *hBemcEtaPhi;
      TH1D *hTpcEta, *hBemcEta;
      TH1D *hTpcPhi, *hBemcPhi;

      map<unsigned int, TVector3> corrections[nRomanPots];
      map<unsigned int, TVector3> offsets[nRomanPots];

      map<int,int> mInstLumiPerRun;

      Long64_t nTracksBEMC, nClustersBEMC, nTracksTPC, nTracksTOF, nEvents, nVertices;
      Double_t tpcEtaSum, bemcEtaSum;
      Double_t tpcPhiSum, bemcPhiSum;
      int nEventsAll, nEventsPassed, nEventsJPsi;
      bool JPsiTrigger;


};

#endif
