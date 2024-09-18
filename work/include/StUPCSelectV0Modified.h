
#ifndef StUPCSelectV0Modified_h
#define StUPCSelectV0Modified_h

class StUPCSelectV0Modified {

 public:

  std::vector<int> selectTracks(std::vector<int> tracks, StUPCEvent *upcEvent, int nVertices);
  std::vector<int> selectTracksK0(std::vector<int> tracks, StUPCEvent *upcEvent, TVector3 const & vertex);
  std::vector<int> selectTracksLambda(std::vector<int> pionTracks, std::vector<int> protonTracks, StUPCEvent *upcEvent, TVector3 const & vertex);

 //private:

    //Util* mUtil;

};

#endif

