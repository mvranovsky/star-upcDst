#ifndef EmbedMaker_h
#define EmbedMaker_h

//_____________________________________________________________________________
//    Class to do embedding of forward-scattered protons (from ppSim)
//    Author: Truhlar Tomas
//_____________________________________________________________________________

#include "UpcDstLibreries.h"

using namespace std;

class EmbedMaker{

 public:
  enum {ErrorCode = -9999,
	kMAXSEC = 2 ,   /// 2 sides
	kMAXCHAIN = 4 , /// 4 chains/planes
	kMAXSVX = 6 ,   
	kMAXSEQ = 8 ,   /// 8 sequencers/roman pots
	kMAXSTRIP = 128 } ;

  EmbedMaker();
  virtual ~EmbedMaker();

  void MakeTracks(float blue_beamenergy, float yellow_beamenergy);
  void MakePMTs();
  
  void Init(TFile *outfile);

  void setMCEvent(StRPEvent* event){ MCEvent = event;}
  void setZBEvent(StRPEvent* event){ ZBEvent = event;}
  void setRPEvent(StRPEvent* event){ rpEvent = event;}
  void setVertex(double X, double Y, double Z){ mXYZ_IP[0] = X; mXYZ_IP[1] = Y; mXYZ_IP[2] = Z; }
  void setEmbedding(bool var){ mDoEmbedding = var;}
  void setAfterburner(bool var){ mDoAfterburner = var;}
  void setOffsetCorrections(TVector3 corr[kMAXSEQ]);

  StRPEvent *getMCEvent() const {return MCEvent;}
  StRPEvent *getZBEvent() const {return ZBEvent;}
  StRPEvent *getRPEvent() const {return rpEvent;}
  void getVertex(double &X, double &Y, double &Z){ X = mXYZ_IP[0]; Y = mXYZ_IP[1]; Z = mXYZ_IP[2]; }
  bool getEmbedding(){ return mDoEmbedding;}
  bool getAfterburner(){ return mDoAfterburner;}


  private:

  enum RP_STATIONS_PER_BRANCH { kRP1, kRP2, kStationsPerBranch };
  enum SILICON_PLANES_PER_COORDINATE { kFirst, kSecond };
  enum COORDINATES { kX, kY, kCoordinates, kZ = kCoordinates };
  enum BeamDirection { kEast, kWest};

  bool mDoEmbedding;
  bool mDoAfterburner;

  // correction for each RP: E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D
  TVector3 correction[kMAXSEQ];

  // additional constants describing the Roman Pot Phase II* system
  static const int kBranches = 4; // number of branches in RP system (EU, ED, WU and WD)
  static const int kPlanesPerCoordinate = kMAXCHAIN/2; // number of Si planes for one spatial coordinate (assumes kMAXCHAIN is even!)

  // reconstruction parameters
  static const int kMaxClusterLength = 5;
  static const int kMaxNumberOfClusterPerPlane = 5;
  static const int kMaxPedestalTAC = 100;
  double kMaxPitchesToMatch;
  static const double kEmin[kMAXSEQ][kMaxClusterLength];

  // handy array
  static const unsigned int kPlanes[kCoordinates][kPlanesPerCoordinate];
  static const int kRpInBranch[kBranches][kStationsPerBranch];
  static const double kPitch[kCoordinates];
  static const double STRIP_PITCH[ kMAXSEQ * kMAXCHAIN ];
  static const double RP_POS_Z[ kMAXSEQ * kMAXCHAIN ];
  static const short orientations2[ kMAXCHAIN * kMAXSEQ ];

  // structures
  struct StRpsHit { // structure representing reconstructed coordinate (one) of track-point
    double mPositionXY;
    double mPositionZ;
    int mClusterId[kPlanesPerCoordinate];
    int mPlanesUsed;
    bool mGolden;
  };


  TH2D *eneSimHits;
  TH2D *numSimHits;
  TH2D *posSimHits;
  TH2D *hmatchDist;
  TH2D *numClusters;
  TH2D *numClustersUsed;
  TH2D *posClusters;
  TH1D *lenClusters;
  TH1D *eneClusters;
  TH2D *numPoints;
  TH2D *posPoints;
  TH2D *hHitSizeMap[kMAXSEQ];
  TH1D *hGoodTrackPoints[kMAXSEQ];

  
  // methods
  void formTracks(const vector<unsigned int> * , const float, const float ) const;
  void formTrackPoints(vector<unsigned int> * ) const;
  vector<EmbedMaker::StRpsHit> formHits(const unsigned int, const int) const;
  void preselectClusters(const unsigned int, const int coordinate, vector<double>*, vector<unsigned int>*, vector<unsigned int>*, vector<unsigned int>*) const;
  Int_t classifyClustersCase(vector<double>*) const;
  Bool_t matchClusters(const int, const int, const vector<double>*, std::vector<unsigned int>*) const;
  Bool_t areMatched(const int, const double, const double, double* = nullptr) const;
  Double_t timeFromTAC(const int, const int, const int, const int) const;

  StRPEvent *MCEvent ; /// pointer for fetching data from MC input ( = ppSim output)
  StRPEvent *ZBEvent ; /// pointer for fetching StRPEvent from ZB data
  StRPEvent *rpEvent ; /// pointer to StRPEvent where embedding sample is stored

  typedef pair<Int_t, Double_t> HitChannel ; /// HitChannel : a pair in which first -> Position ; second -> Energy

  static Bool_t hitcompare (HitChannel A, HitChannel B) { return (A.first<B.first); }

  Double_t  mPedave[kMAXSEQ][kMAXCHAIN][kMAXSVX][kMAXSTRIP] ;
  Double_t  mPedrms[kMAXSEQ][kMAXCHAIN][kMAXSVX][kMAXSTRIP] ;

  unsigned char mRpStatus[kMAXSEQ] ;

  UChar_t mSiliconBunch ;

  // K. Yip (2015-10-22) : Adding variables for Accelerator and Skew parameters in the respective database ;
  double mSkew_param[kMAXSEQ][2][4] ; // 4 parameters for each PMT and there are 2 PMT's for each of the 8 RP
  double mXYZ_IP[3]; /* collision coordinates at the IP; 0=X, 1=Y, 2=Z */
  double mThetaXY_tilt[2]; /* tilt angles of the beam at collision; 0=X, 1=Y */
  double mDistanceFromIPtoDX[2]; /* distance from the IP to the DX magnet in the East and West; 0=E, 1=W */
  double mLDX[2];     /* length of DX in the East and West; 0=E, 1=W */
  double mBendingAngle[2];     /* DX bending angles in the East and West; 0=E, 1=W */
  double mConversion_TAC_time ; /* converting the TAC tick to time (second) */

};

#endif

