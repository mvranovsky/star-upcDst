
#ifndef StUPCSelectV0_h
#define StUPCSelectV0_h

class StUPCSelectV0 {

 public:

  int selectTracks(std::vector<StUPCTrack>& tracks, std::vector<UChar_t>& sel, StUPCEvent *upcEvent, TVector3 const & vertex, vector<TH1D*>& hists1D, vector<TH2D*>& hists2D);

 private:
  void fillHists(StUPCV0& K0);  


};

#endif

