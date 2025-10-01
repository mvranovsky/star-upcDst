#ifndef RUNDEF_H
#define RUNDEF_H

#include <vector>
#include <bitset>
#include "TString.h"


const bool runAnaBP = false;
const bool runAnaV0 = false;
const bool runAnaV0Mult = false;
const bool runTofEff = false;   
const bool runTofEffMult = false;
const bool runAnaJPsi = false;
const bool runAnaJPSI = false;
const bool runAnaGoodRun = true;
const bool runEmbeddingJPsi = false;
const bool runSysStudy = false;
const bool runSysStudyEmbedding = false;


const unsigned int nAligIteration = 4;
const bool runMCAna = false;
const bool DEBUG = false;
const bool V0Control = false;
const bool analysisWithRPs = false;
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 
// 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL 
// 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const std::vector<int> triggerID = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
const std::vector<int> CEPtriggers = { 570705 }; //, 590705 //{ 570701, 570705, 570711, 590701, 590705, 590708};
           
const std::vector<int> JPSItriggers = { 570209, 570219, 570229 };

const std::vector<int> ZBtriggers = {570704};
const int nPlanesUsed = 3;
const bool AFTERBURNER = true;

//const TString nameOfTree[] = { TString("recTreeBP"), TString("recTreeV0"), TString("recTreeV0Control")};
//const bitset<16> treeBits[] = { bitset<16>(string("0000000000011111")), bitset<16>(string("0000011100010001")), bitset<16>(string("0000011100010001"))};
//postupne odzadu: eventinfo, vertexinfo(1), stateinfo(1), centralhadronsinfo(2), rptracksinfo(2), zdc,bbc info, trigger efficiency info, TRUEMC info, vertexinfo(multiple), stateinfo(multiple), centralhadronsinfo(multiple)

const TString nameOfAnaBPTree = "recTreeBP";
const std::bitset<16> AnaBPTreeBits = std::bitset<16>(std::string("0000011100010001"));
const TString nameOfAnaBPDir = "AnaBPPlots";

const TString nameOfAnaV0Tree = "recTreeAnaV0";
const std::bitset<16> AnaV0TreeBits = std::bitset<16>(std::string("0000000000001111")); //RPs info turned off
const TString nameOfAnaV0Dir = "AnaV0Plots";

const TString nameOfAnaV0MultTree = "recTreeAnaV0Mult";
const std::bitset<16> AnaV0MultTreeBits = std::bitset<16>(std::string("0000011100010001")); 
const TString nameOfAnaV0MultDir = "AnaV0MultPlots";

const TString nameOfTofEffTree = "recTreeTofEff";
const std::bitset<16> TofEffTreeBits = std::bitset<16>(std::string("0000000000001111")); //RPs info turned off
const TString nameOfTofEffDir = "TofEffPlots";

const TString nameOfTofEffMultTree = "recTreeTofEffMult";
const std::bitset<16> TofEffMultTreeBits = std::bitset<16>(std::string("0000011100010001")); 
const TString nameOfTofEffMultDir = "TofEffMultPlots";

const TString nameOfAnaJPsiTree = "recTreeAnaJPsi";
const std::bitset<16> AnaJPsiTreeBits = std::bitset<16>(std::string("0011000000011111")); 
const TString nameOfAnaJPsiDir = "AnaJPsiPlots";

const TString nameOfAnaJPSITree = "recTreeAnaJPSI";
const std::bitset<16> AnaJPSITreeBits = std::bitset<16>(std::string("0011010000010111")); 
const TString nameOfAnaJPSIDir = "AnaJPSIPlots";

const TString nameOfEmbeddingJPsiTree = "recTreeEmbeddingJPsi";
const std::bitset<16> EmbeddingJPsiTreeBits = std::bitset<16>(std::string("0011000000001101")); 
const TString nameOfEmbeddingJPsiDir = "EmbeddingJPsiPlots";

const TString nameOfAnaGoodRunTree = "recAnaGoodRun";
const std::bitset<16> AnaGoodRunTreeBits = std::bitset<16>(std::string("0100000000000000")); 
const TString nameOfAnaGoodRunDir = "AnaGoodRunPlots";

const TString nameOfZBTree = "recTreeZB";
const std::bitset<16> ZBTreeBits = std::bitset<16>(std::string("0000100000001001")); 


const TString offsetFilePath = "/star/u/mvranovsk/star-upcDst/work/OffSetsRun17.list";
const TString offsetCorrectionsFilePath = "/star/u/mvranovsk/star-upcDst/work/OffSetsCorrectionsRun17.list";

const TString YAxisDescription = "counts";

const bool TOF2Tracks = false;  // condition whether both tracks are to be matched to ToF or not
const bool usePrimVtx = false;  // condition whether to use or not to use primary vertex in the analysis

// TPC good track quality cuts
const int minNHitsFit = 15;
const int minNHitsDEdx = 15;
const double minPt = 0.2;
const double maxDcaXY = 1.5;
const double minDcaZ = -1.0;
const double maxDcaZ = 1.0;

// loose conditions for Systematic study
const int minNHitsFitLoose = 12;
const int minNHitsDEdxLoose = 12;
const double maxDcaXYLoose = 1.8;
const double minDcaZLoose = -1.2;
const double maxDcaZLoose = 1.2;
const double maxEtaLoose = 1.0;
const double vertexRangeLoose = 120;
const double maxPidChiEELoose = 11;
const double minPidChiPPLoose = 8;
const double minPidChiPiPiLoose = 8;
const double minPidChiKKLoose = 8;

//topology cuts for V0 selection
const double vertexDiffMax = 5.;
const double maxDcaDaughters = 1.5;
const double maxDcaBeamLine = 1.5;
const double minPointingAngle = 0.97;
const double maxDecayLengthHypo = 3.;

// specific conditions for special eta-vtxZ cut
const double maxEta = 0.9;
const double vertexRange = 100.0; // cm
const double etaVertexSlope = -1/250.0;
const double etaVertexShift = 0.9;

// vertex range in z-coordinate
//const double exclusivityCut = 0.7;
//const double exclusivityCut = 0.1;

// cuts for JPsi analysis
const double maxPidChiEE = 9;
const double minPidChiPP = 10;
const double minPidChiPiPi = 10;
const double minPidChiKK = 10;


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
