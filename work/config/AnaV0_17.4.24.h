#ifndef RUNDEF_H
#define RUNDEF_H

//const bool runMAINANA = false;
const bool runAnaBP = false;
const bool runAnaV0 = true;
const bool runAnaV0Control = false;
//const bool runRPMCANA = false;

const bool runStudy[] = { runAnaBP, runAnaV0, runAnaV0Control};

const unsigned int nAligIteration = 4;
const bool DEBUG = false;
const bool V0Control = false;
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 
// 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL 
// 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const vector<int> triggerID = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
const vector<int> CEPtriggers = { 570705, 590705 }; //{ 570701, 570705, 570711, 590701, 590705, 590708};
//const int CEPtriggers[] = { 570701, 590711 }; //{ 570701, 570705, 570711, 590701, 590705, 590708};
//const int CEPtriggers[] = { 570701, 590711 }; //{ 570701, 570705, 570711, 590701, 590705, 590708};


const int nPlanesUsed = 3;
const bool AFTERBURNER = true;

const TString nameOfTree[] = { TString("recTreeBP"), TString("recTreeV0"), TString("recTreeV0Control")};
const bitset<16> treeBits[] = { bitset<16>(string("0000000000011111")), bitset<16>(string("0000011100010001")), bitset<16>(string("0000011100010001"))};
//postupne odzadu: eventinfo, vertexinfo(1), stateinfo(1), centralhadronsinfo(2), rptracksinfo(2), zdc,bbc info, trigger efficiency info, TRUEMC info, vertexinfo(multiple), stateinfo(multiple), centralhadronsinfo(multiple)
const TString nameOfAnaBPTree = "recTreeBP";
//const bitset<16> AnaBPTreeBits = bitset<16>(string("0000000000011111"));
const bitset<16> AnaBPTreeBits = bitset<16>(string("0000011100010001"));
const TString nameOfAnaV0Tree = "recTreeV0";
const bitset<16> AnaV0TreeBits = bitset<16>(string("0000011100010001")); 
const TString nameOfAnaV0ControlTree = "recTreeV0Control";
const bitset<16> AnaV0ControlTreeBits = bitset<16>(string("0000011100010001")); 

const TString YAxisDescription = "counts";
//const TString pathToOfLumiFile = "/star/u/truhlar/star-upcDst/work/lists/luminosityForZB.list";
const TString nameOfOffSetFile = "OffSetsRun17.list";
const TString nameOfOffSetCorrectionFile = "OffSetsCorrectionsRun17.list";

// RP MC 
//const string pathToMC = "/star/data01/pwg_tasks/upc02/run17MC/";
//const int nMcEventsPerZbEvent = 10;
//const double maxMcToRecoDist = 0.08; // max distance between MC generated and reconstructed is 80 MeV in px, py space
    
// Elastic event conditions:
//const unsigned int nSigma = 3;
//const double sigma = 0.00013; // rad, for Main period
//const double PHI_CENT  = 90.0*0.0174533; // geometrical window
//const double PHI_WIDTH = 12.0*0.0174533; // geometrical window
//const double DCUT = 3*0.0007; // in meters

// TPC good track quality cuts
const int minNHitsFit = 25;
const int minNHitsDEdx = 15;
const double minPt = 0.2;
const double maxDcaXY = 1.5;
const double minDcaZ = -1.0;
const double maxDcaZ = 1.0;
const double maxEta = 0.7;

//topology cuts for V0 selection
const double maxDcaDaughters = 1.5;
const double maxDcaBeamLine = 1.5;
const double minPointingAngle = 0.995;
const double maxDecayLengthHypo = 3;

// vertex range in z-coordinate
const double vertexRange = 80.0;
// Exclusivity cut
//const double exclusivityCut = 0.7;
const double exclusivityCut = 0.1;

const double sqrtOfTwo = 1.41421356237;

// Plot setting
const double textSize = 0.04;
const double labelSize = 0.05;
const double labelTextSize = 0.06;
const int fontStyle = 42;
const TString yAxisTitle = "Normalized counts";

const int mainColor = 4;
const int mainMarker = 20;

const int bckgColor = 2;

#endif //ifndef RUNDEF_H
