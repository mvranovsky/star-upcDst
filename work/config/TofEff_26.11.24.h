#ifndef RUNDEF_H
#define RUNDEF_H

#include <vector>
#include <bitset>
#include "TString.h"


const bool runAnaBP = false;
const bool runAnaV0 = false;
const bool runAnaV0Mult = false;
const bool runMCAna = false;
const bool runTofEff = true;
const bool runTofEffMult = false;

const unsigned int nAligIteration = 4;
const bool DEBUG = false;
const bool V0Control = false; 
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 
// 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL 
// 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const std::vector<int> triggerID = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
const std::vector<int> CEPtriggers = { 570705 }; //, 590705 //{ 570701, 570705, 570711, 590701, 590705, 590708};


const int nPlanesUsed = 3;
const bool AFTERBURNER = true;

//const TString nameOfTree[] = { TString("recTreeBP"), TString("recTreeV0"), TString("recTreeV0Control")};
//const bitset<16> treeBits[] = { bitset<16>(string("0000000000011111")), bitset<16>(string("0000011100010001")), bitset<16>(string("0000011100010001"))};
//postupne odzadu: eventinfo, vertexinfo(1), stateinfo(1), centralhadronsinfo(2), rptracksinfo(2), zdc,bbc info, trigger efficiency info, TRUEMC info, vertexinfo(multiple), stateinfo(multiple), centralhadronsinfo(multiple)

const TString nameOfAnaBPTree = "recTreeBP";
const std::bitset<16> AnaBPTreeBits = std::bitset<16>(std::string("0000011100010001"));

const TString nameOfAnaV0Tree = "recTreeAnaV0";
const std::bitset<16> AnaV0TreeBits = std::bitset<16>(std::string("0000000000001111")); //RPs info turned off

const TString nameOfAnaV0MultTree = "recTreeAnaV0Mult";
const std::bitset<16> AnaV0MultTreeBits = std::bitset<16>(std::string("0000011100010001")); 

const TString nameOfTofEffTree = "recTreeTofEff";
const std::bitset<16> TofEffTreeBits = std::bitset<16>(std::string("0000000000001111")); //RPs info turned off

const TString nameOfTofEffMultTree = "recTreeTofEffMult";
const std::bitset<16> TofEffMultTreeBits = std::bitset<16>(std::string("0000011100010001")); 


const TString YAxisDescription = "counts";
//const TString pathToOfLumiFile = "/star/u/truhlar/star-upcDst/work/lists/luminosityForZB.list";
const TString nameOfOffSetFile = "OffSetsRun17.list";
const TString nameOfOffSetCorrectionFile = "OffSetsCorrectionsRun17.list";


// TPC good track quality cuts
const int minNHitsFit = 25;
const int minNHitsDEdx = 15;
const double minPt = 0.2;
const double maxDcaXY = 1.5;
const double minDcaZ = -1.0;
const double maxDcaZ = 1.0;
const bool TOF2Tracks = false;  // condition whether both tracks are to be matched to ToF or not
const bool usePrimVtx = false;  // condition whether to use or not to use primary vertex in the analysis

//topology cuts for V0 selection
const double vertexDiffMax = 5.;
const double maxDcaDaughters = 1.5;
const double maxDcaBeamLine = 1.5;
const double minPointingAngle = 0.97;
const double maxDecayLengthHypo = 3.;

// specific conditions for special eta-vtxZ cut
const double maxEta = 0.9;
const double minEta = 0.;
const double vertexRangeForEVz = 100;
const double etaVertexSlope = -1/250.0;
const double etaVertexShift = 0.9;

// vertex range in z-coordinate
const double vertexRange = 60.0; // cm
//const double exclusivityCut = 0.7;
//const double exclusivityCut = 0.1;

//const double sqrtOfTwo = 1.41421356237;

// Plot setting
const double textSize = 0.04;
const double labelSize = 0.05;
const double labelTextSize = 0.06;
const int fontStyle = 42;
const TString yAxisTitle = "counts";

const int mainColor = 4;
const int mainMarker = 20;

const int bckgColor = 2;


#endif //ifndef RUNDEF_H
