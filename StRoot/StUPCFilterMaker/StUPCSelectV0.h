
#ifndef StUPCSelectV0_h
#define StUPCSelectV0_h

class StUPCSelectV0 {

 public:

  int selectTracks(std::vector<StUPCTrack>& tracks, std::vector<UChar_t>& sel, StUPCEvent *upcEvent, TVector3 const & vertex);

};

#endif

