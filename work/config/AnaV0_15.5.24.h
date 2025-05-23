#ifndef RUNDEF_H
#define RUNDEF_H

#include <vector>
#include <bitset>
#include "TString.h"


//const bool runMAINANA = false;
const bool runAnaBP = false;
const bool runAnaV0 = false;
const bool runAnaV0Control = false;
const bool runAnaV0SingleState = true;
//const bool runRPMCANA = false;

const bool runStudy[] = { runAnaBP, runAnaV0, runAnaV0Control};

const unsigned int nAligIteration = 4;
const bool DEBUG = false;
const bool V0Control = false; // used in runAnaV0, runAnaSingleState: true=kPrimary, false=kV0
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 
// 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL 
// 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
const std::vector<int> triggerID = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
const std::vector<int> CEPtriggers = { 570705, 590705 }; //{ 570701, 570705, 570711, 590701, 590705, 590708};
//const int CEPtriggers[] = { 570701, 590711 }; //{ 570701, 570705, 570711, 590701, 590705, 590708};
//const int CEPtriggers[] = { 570701, 590711 }; //{ 570701, 570705, 570711, 590701, 590705, 590708};


const int nPlanesUsed = 3;
const bool AFTERBURNER = true;

//const TString nameOfTree[] = { TString("recTreeBP"), TString("recTreeV0"), TString("recTreeV0Control")};
//const bitset<16> treeBits[] = { bitset<16>(string("0000000000011111")), bitset<16>(string("0000011100010001")), bitset<16>(string("0000011100010001"))};
//postupne odzadu: eventinfo, vertexinfo(1), stateinfo(1), centralhadronsinfo(2), rptracksinfo(2), zdc,bbc info, trigger efficiency info, TRUEMC info, vertexinfo(multiple), stateinfo(multiple), centralhadronsinfo(multiple)
const TString nameOfAnaBPTree = "recTreeBP";
//const bitset<16> AnaBPTreeBits = bitset<16>(string("0000000000011111"));
const std::bitset<16> AnaBPTreeBits = std::bitset<16>(std::string("0000011100010001"));
const TString nameOfAnaV0Tree = "recTreeV0";
const std::bitset<16> AnaV0TreeBits = std::bitset<16>(std::string("0000011100010001")); 
const TString nameOfAnaV0ControlTree = "recTreeV0Control";
const std::bitset<16> AnaV0ControlTreeBits = std::bitset<16>(std::string("0000011100010001")); 
const TString nameOfAnaV0SingleStateTree = "recTreeV0SingleState";
const std::bitset<16> AnaV0SingleStateTreeBits = std::bitset<16>(std::string("0000000000001111")); //RPs info turned off



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
const double minPt = 0.15;
const double maxDcaXY = 1.5;
const double minDcaZ = -1.0;
const double maxDcaZ = 1.0;
const double maxEta = 0.9;
const bool TOF2Tracks = false;

//topology cuts for V0 selection
const double vertexDiffMax = 3.;
const double maxDcaDaughters = 1.5;
const double maxDcaBeamLine = 1.5;
const double minPointingAngle = 0.95;
const double maxDecayLengthHypo = 3.;

// vertex range in z-coordinate
const double vertexRange = 100.0;
// Exclusivity cut
//const double exclusivityCut = 0.7;
const double exclusivityCut = 0.1;

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


const std::vector<std::pair<TString,bool>> plotsV0 = {
    {"invMassK0s", true},
    {"invMassLambda", true},
    {"hNSigmaPiPcorr", true},
    {"hNSigmaPiKcorr",true},
    {"hNSigmaPKcorr",true},
    {"hNPairV0", true},
    {"hPosZ", true},
    {"hPosZCut", true},
    {"hVtxDiff", true},
    {"hDcaDaughters", true},
    {"hDcaDaughtersCut", true},
    {"hDcaBeamline", true},
    {"hDcaBeamlineCut", true},
    {"hPointingAngle", true},
    {"hPointingAngleCut", true},
    {"hDecayLength", true},
    {"hDecayLengthCut", true},
    {"hEta",true},
    {"hPt",true},
    {"hNfitHits",true},
    {"hNhitsDEdx", true},
    {"hAnalysisFlow", true},
    {"hNSigmaPi", true},
    {"hNSigmaP", true},
    {"hNSigmaK", true},
    {"hDecayLPointingA", true},
    {"hDecayLPointingACut", true},
    {"hArmenterosPodolanski", true},
    {"hInvMassEta", true}
};

const std::vector<std::pair<TString,bool>> plotsV0Control = {
    {"invMassK0s", true}, 
    {"invMassLambda", true},
    {"hNSigmaPiPcorr", true},
    {"hNSigmaPiKcorr",true},
    {"hNSigmaPKcorr",true},
    {"hNPairV0", true},
    {"hPosZ", true},
    {"hPosZCut", true},
    {"hVtxDiff", true},
    {"hDcaDaughters", true},
    {"hDcaDaughtersCut", true},
    {"hDcaBeamline", true},
    {"hDcaBeamlineCut", true},
    {"hPointingAngle", true},
    {"hPointingAngleCut", true},
    {"hDecayLength", true},
    {"hDecayLengthCut", true},
    {"hEta",true},
    {"hPt",true},
    {"hNfitHits",true},
    {"hNhitsDEdx", true},
    {"hAnalysisFlow", true},
    {"hNSigmaPi", true},
    {"hNSigmaP", true},
    {"hNSigmaK", true},
    {"hDecayLPointingA", true},
    {"hDecayLPointingACut", true},
    {"hArmenterosPodolanski", true},
    {"hInvMassEta", true}
};

const std::vector<std::pair<TString,bool>> plotsBP = {
    {"invMassK0s", true},
    {"invMassLambda", true},
    {"hNSigmaPiPcorr", true},
    {"hNSigmaPiKcorr",true},
    {"hNSigmaPKcorr",true},
    {"hNPairV0", true},
    {"hPosZ", true},
    {"hPosZCut", true},
    {"hVtxDiff", true},
    {"hDcaDaughters", true},
    {"hDcaDaughtersCut", true},
    {"hDcaBeamline", true},
    {"hDcaBeamlineCut", true},
    {"hPointingAngle", true},
    {"hPointingAngleCut", true},
    {"hDecayLength", true},
    {"hDecayLengthCut", true},
    {"hEta",true},
    {"hPt",true},
    {"hNfitHits",true},
    {"hNhitsDEdx", true},
    {"hAnalysisFlow", true},
    {"hNSigmaPi", true},
    {"hNSigmaP", true},
    {"hNSigmaK", true},
    {"hDecayLPointingA", true},
    {"hDecayLPointingACut", true},
    {"hArmenterosPodolanski", true},
    {"hInvMassEta", true}
};



const TString cutDescriptions[23] = {
    TString("Number of pairs per event"),
    TString("Vertex_{Z}"),
    TString("|Vertex_{Z}| < 80 cm"), 
    TString("Vertex difference"), 
    TString("|DCA_{daughters}|"), 
    TString("|DCA_{daughters}| < ") + std::to_string(maxDcaDaughters) + TString(" cm"), 
    TString("|DCA_{beamline}|") + std::to_string(maxDcaBeamLine) + TString(" cm"), 
    TString("|DCA_{beamline}| < ") + std::to_string(minPointingAngle), 
    TString("cos(pointingAngle)") + std::to_string(maxDecayLengthHypo) + TString(" cm"), 
    TString("cos(pointingAngle) > "), 
    TString("|Decay Length|"), 
    TString("|Decay Length| < "), 
    TString("#eta distribution"), 
    TString("pT > 0.2 GeV"), 
    TString("N^{fit}_{hits} > 25"), 
    TString("N^{dE/dx}_{hits} > 15"), 
    TString(""), 
    TString(""), 
    TString(""), 
    TString(""),
    TString(""), 
    TString(""),
    TString("")
};


#endif //ifndef RUNDEF_H
