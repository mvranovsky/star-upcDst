#ifndef Ana_h
#define Ana_h

// include headers
#include "UpcDstLibreries.h"
#include "Util.h"
#include "RecTree.h"

using namespace std;
using namespace UTIL;

class Ana{
   public:
      Ana(TFile *outFile);
      virtual ~Ana();

      bool RPInFidRange(double x, double y) const;
      bool IsInRpRange(double x, double y, int rpId, TVector3 offSet) const;
      bool IsGoodTrack(const StUPCTrack *trk) const;
      bool IsGoodTofTrack(const StUPCTrack *trk) const;
      bool CheckTriggers(const vector<int> *triggerArray, StUPCEvent *mUpcEvt, TH1D *hTriggerBits) const;

      void AnaRpTracks(StRPEvent *event);

      virtual void Make(){cout<<"Hi my name is make"<<endl;};
      virtual void Init(){cout<<"Hi I should not be there but I am"<<endl;};
      void SetEvent(StUPCEvent *upcEvt, StRPEvent *rpEvt, StRPEvent *mcEvt);

      inline void SetAnaName(TString name){ anaName = name; }
      inline void SetTriggers(const vector<int> *trigg){ trigger = trigg; }

      void SaveEventInfo(const StUPCEvent *upcEvt);
      void SaveRPinfo(const StUPCRpsTrack *trackRP, unsigned int iSide);
      void SaveTrackInfo(const StUPCTrack *trk, unsigned int iTrack);
      void SaveTrackInfo(const StUPCTrack *trk, TLorentzVector hadron ,unsigned int iTrack);
      void SaveStateInfo(TLorentzVector state,int totQ, unsigned int iState);
      void SaveVertexInfo(const TVector3 vrtx, Double_t dcaDaughters, Double_t dcaBeamline, Double_t pointingAngle, Double_t decayLength, Double_t vertexDiff, unsigned int iVtx);
      void SaveZdcInfo(const StUPCEvent *upcEvt);
      void SaveBbcInfo(const StUPCEvent *upcEvt);
      void SaveTriggerInfo(const StUPCEvent *upcEvt, const StRPEvent *rpEvt);
      void saveRpTrigBit(const StRPEvent *rpEvt);
      bool IsRpTrigBit(const StRPEvent *rpEvt, unsigned int iRp);

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

      TString anaName;
      const vector<int> *trigger;

};

#endif
