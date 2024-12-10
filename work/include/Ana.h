#ifndef Ana_h
#define Ana_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"
#include "StUPCV0.h"

using namespace std;
using namespace UTIL;

class Ana{
   public:
      Ana(TFile *outFile);
      virtual ~Ana();

      bool RPInFidRange(double x, double y) const;
      bool IsInRpRange(double x, double y, int rpId, TVector3 offSet) const;


      //inline bool RPInFidRange(double x, double y) const { return (abs(y) < fpPyMax[EU] && abs(y) > fpPyMin[EU] && x > fpPxMin[EU] && (x + fpPxCenter[EU])*(x + fpPxCenter[EU]) + y*y < fpPRadius[EU]) ? true : false;}
      //inline bool RPInElFidRange(double x, double y) const { return (abs(y) < fpElPyMax && abs(y) > fpElPyMin && x > fpElPxMin && x < fpElPxMax) ? true : false;}
      //inline bool IsInRpRange(double x, double y, int rpId, TVector3 offSet) const { return (abs(y) < abs(offSet[Y]) || abs(y) > abs(offSet[Y]) + 0.03 || x < -0.02 || x > 0.02) ? false : true; };

      inline bool IsGoodTrack(const StUPCTrack *trk) const { return ( trk->getNhitsFit() > minNHitsFit && trk->getNhitsDEdx() > minNHitsDEdx && trk->getPt() > minPt);}
      inline bool IsGoodTofTrack(const StUPCTrack *trk) const {return (trk->getTofTime() > 0 && trk->getTofPathLength() > 0);}
      
      inline double vertexEtaRange(int rSide, int idx ) const{ return mRecTree->IsVertexSet(idx) ? etaVertexSlope*mRecTree->getVertexZInCm(idx) + (2*rSide - 1)*etaVertexShift : (2*rSide - 1)*maxEta; } // 0 for min, 1 for max
      inline bool IsGoodEtaTrack(const StUPCTrack *trk, int idx) const { return (abs(trk->getEta()) < maxEta && trk->getEta() > vertexEtaRange(0,idx) && trk->getEta() < vertexEtaRange(1, idx) );}
      int hasGoodTPCnSigma(const StUPCTrack *trk); 
      bool CheckTriggers(const vector<int> *triggerArray, StUPCEvent *mUpcEvt, TH1D *hTriggerBits) const;




      void AnaRpTracks(StRPEvent *event);
      virtual void Make(){cout<<"Hi my name is make"<<endl;};
      virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
      void SetEvent(StUPCEvent *upcEvt, StRPEvent *rpEvt, StRPEvent *mcEvt);

      inline void SetAnaName(TString name){ anaName = name; }
      inline void SetTriggers(const vector<int> *trigg){ trigger = trigg; }

      void fillTrackQualityCuts(const StUPCTrack* trk);
      void fillNSigmaPlots(const StUPCTrack *trk);
      void fillEtaVtxPlots(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ);
      void SaveEventInfo(const StUPCEvent *upcEvt);
      void SaveRPinfo(const StUPCRpsTrack *trackRP, unsigned int iSide);
      //void SaveTrackInfo(const StUPCTrack *trk, unsigned int iTrack);
      void SaveTrackInfo(const StUPCTrack *trk, TLorentzVector hadron ,unsigned int iTrack);
      void SaveStateInfo(TLorentzVector state,int totQ, unsigned int iState);
      void SaveVertexInfo(const StUPCV0* V0, unsigned int iVtx);
      void SaveZdcInfo(const StUPCEvent *upcEvt);
      void SaveBbcInfo(const StUPCEvent *upcEvt);
      void SaveTriggerInfo(const StUPCEvent *upcEvt, const StRPEvent *rpEvt);
      void saveRpTrigBit(const StRPEvent *rpEvt);
      bool IsRpTrigBit(const StRPEvent *rpEvt, unsigned int iRp);
      void fillBeamlineInfo();

      void resetInfo();
   
   protected:
      TFile *mOutFile;
      RecTree *mRecTree; 

      StUPCEvent *mUpcEvt;
      StRPEvent *mRpEvt, *mMcEvt;

      Util* mUtil;

      vector<unsigned int> mRpTrackIdVec_perSide[nSides];
      vector<unsigned int> mRpTrackIdVec_perBranch[nBranches];
      vector<unsigned int> mTrackPointIdVec[nRomanPots];
      vector<unsigned int> mClusterIdVec[nRomanPots*nPlanes];

      // Control plots
      TH1D *hAnalysisFlow; 

      TH1D* hEta,*hEtaCut, *hPosZ, *hPosZCut;
      TH2F* hEtaPhi, *hEtaPhiCut, *hEtaVtxZ, *hEtaVtxZCut; 

      TString anaName;
      const vector<int> *trigger;

      vector<int> hadronID, tagID;


      double bField;
      double beamline[4]; 


};

#endif
