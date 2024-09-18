#ifndef PlotUtil_h
#define PlotUtil_h

// include headers
#include "Libreries.h"
#include "RecTree.h"
using namespace std;

#include "Util.h"
using namespace UTIL;
// include definitions of plot style
#include "../config/__CONFIG__"

// include plot macros
#include "PlotProbOfRetainEvent.h"
#include "PlotRPMCPlots.h"
#include "PlotTrigEff.h"
#include "PlotVertexStudy.h"
#include "PlotTofQA.h"
#include "PlotEmbedQA.h"
#include "PlotMainAna.h"

TFile *inFile, *embFile, *outFile;
Util *mUtil;

TTree *mTree[nStudies];
RecTree *mRecTree[nStudies];

TCanvas *canvas;
TLegend *legend;
TPaveText *text;

// general stuff
bool runPlots(TString inputFile, TString embedFile);
bool Init(TString inputFile, TString embedFile);
bool ConnectInput( TString inputFile, TFile **file); 
TFile *CreateOutputFile(const string& out);
void CreateCanvas(TCanvas **canvas, TString canvasName);
void SetHistStyle(TH1* hist, Int_t color, Int_t markStyle);
void SetTH2Style(TH2* hist);
void CreateLegend(TLegend **legend);
void DrawSTARInternal(double xl = 0.75, double yl = 0.89, double xr = 0.88, double yr = 0.93);
void CreateText(double xl, double yl, double xr, double yr);
void SetGPad(double xl = 0.09, double yl = 0.16, double xr = 0.09, double yr = 0.015); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
void DrawFiducial();
void SetLineStyle(TLine* line);

#endif
