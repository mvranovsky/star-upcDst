
#ifndef StUPCSelectV0_h
#define StUPCSelectV0_h

#include "TH1D.h"
#include "TH2D.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCV0.h"
#include "TVector3.h"
#include <vector> 

using namespace std;


class StUPCSelectV0 {

 public:

  int selectTracks(std::vector<StUPCTrack>& tracks, std::vector<UChar_t>& sel, StUPCEvent *upcEvent, TVector3 const & vertex, vector<TH1D*>& hists1D, vector<TH2D*>& hists2D);


};

#endif

