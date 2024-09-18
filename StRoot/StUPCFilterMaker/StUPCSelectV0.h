
#ifndef StUPCSelectV0_h
#define StUPCSelectV0_h

#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"

class StUPCSelectV0 {

 public:

  int selectTracks(std::vector<StUPCTrack>& tracks, std::vector<UChar_t>& sel, StUPCEvent *upcEvent, TVector3 const & vertex, TList*& mHistList);
  bool loadHists(TList*& mHistList);


  TH1D *hDcaDaughters, *hDcaBeamline, *hPointingAngle, *hDecayLength, *hEta, *hEtaAfter;
  TH2D *hDecayLPointingA, *hDcaBeamDaughters;

};

#endif

