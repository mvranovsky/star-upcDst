//
//root4star [0] .L bichselG10.C 
//root4star [1] .x bichselG10.C("/gpfs01/star/pwg/truhlar/Run17_P20ic/mainAna0304/merged/StRP_production_0000.root") 
//
#if !defined(__CINT__)
// code that should be seen ONLY by the compiler
#else
#if !defined(__CINT__) || defined(__MAKECINT__)
// code that should be seen by the compiler AND rootcint
#else
// code that should always be seen
#endif
#endif
//________________________________________________________________________________
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include <stdio.h>
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TClassTable.h"
#include "StBichsel/Bichsel.h"
#include "StBichsel/StdEdxModel.h"
#include "TLegend.h"
#include "TROOT.h"
#else
class Bichsel;
#endif
Bichsel *m_Bichsel = 0;
const Int_t NMasses = 10;
const Double_t Masses[NMasses] = {0.13956995,
               0.493677,
               0.93827231,
               1.87561339,
               0.51099907e-3,
               0.1056584,
               2.80925,
               2.80923, //GEANT3
               3.727417, //GEANT3
               0.13956995,
};
const Int_t   Index[NMasses] = { 2,    3,   4,   5,   0,    1,  6,    7,       8,    -2};
//const Int_t   Index[NMasses] = { 4,    3,   2,   0,   5,    1,  6,    7,       8,    -2};
const Int_t   Colors[NMasses] = { 1,    2,   4,   7,   6,    3,  9,    30,       8,    -2};
const Char_t *Names[NMasses] = {"Pion","Kaon","Proton","Deuteron","Electron","#mu","t","He3","#alpha","2#pi"};
//const Char_t *Names[NMasses] = {"Proton", "Kaon","Pion","Electron", "Deuteron","#mu","t","He3","#alpha","2#pi"};
const Int_t NF = 5;  //         0       1    2     3     4      5   6     7
const Char_t *FNames[8] = {"Girrf","Sirrf","Bz","B70","B60","B70M","dNdx","BzM"};
const Int_t Nlog2dx = 3;
const Double_t log2dx[Nlog2dx] = {0,1,2};
//________________________________________________________________________________
Double_t bichselZ(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*charge*charge*TMath::Exp(m_Bichsel->GetMostProbableZ(TMath::Log10(poverm),dx2)));//TMath::Exp(7.81779499999999961e-01));
   //return charge*charge*TMath::Log10(m_Bichsel->GetI70(TMath::Log10(poverm),1.));
}
//________________________________________________________________________________
Double_t bichselZM(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*charge*charge*TMath::Exp(m_Bichsel->GetMostProbableZM(TMath::Log10(poverm),dx2)));//TMath::Exp(7.81779499999999961e-01));
   //return charge*charge*TMath::Log10(m_Bichsel->GetI70(TMath::Log10(poverm),1.));
}
//________________________________________________________________________________
Double_t bichsel70(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   // return  TMath::Log10(scale*charge*charge*m_Bichsel->GetI70M(TMath::Log10(poverm),dx2));//TMath::Exp(7.81779499999999961e-01));
   return charge*charge*TMath::Log10(m_Bichsel->GetI70(TMath::Log10(poverm),1.));
}
//________________________________________________________________________________
Double_t bichsel70M(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*charge*charge*m_Bichsel->GetI70M(TMath::Log10(poverm),dx2));//TMath::Exp(7.81779499999999961e-01));
}
//________________________________________________________________________________
Double_t dNdx(Double_t *x,Double_t *par) {
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t scale = 1;
   Double_t mass = par[0];
   if (mass < 0) {mass = - mass; scale = 2;}
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   Double_t dx2 = 1;
   if (par[1] > 1.0) {
      charge = 2;
      poverm *= charge;
      dx2 = TMath::Log2(5.);
   }
   return  TMath::Log10(scale*StdEdxModel::instance()->dNdx(poverm,charge));//TMath::Exp(7.81779499999999961e-01));
}
#ifndef __CINT__
//________________________________________________________________________________
Double_t aleph70P(Double_t *x,Double_t *par) {
   /* 
       W.Blum, L. Rolandi "Particle Detection with Drift Chambers", page 246, eq. (9.5)
       F_g(v) = p[0]/beta**p[3]*(p[1] - beta**p[3] - log(p[2] + (beta*gamma)**-p[4]);
       F_g(v) = p[0]*(1/beta**p[3]*(p[1] - log(p[2] + 1/(beta*gamma)**p[4])) - 1) 
   */
   Double_t bg = x[0];
   Double_t b2inv = 1. + 1./(bg*bg);
   Double_t beta  = 1./TMath::Sqrt(b2inv);
   Double_t dEdx = par[0]*(-1 + TMath::Power(beta,-par[3])*(par[1] - TMath::Log(TMath::Max(1e-10,par[2] + TMath::Power(bg,-par[4])))));
   return dEdx;
};
//________________________________________________________________________________
Double_t aleph70(Double_t *x,Double_t *par) {
   static const Double_t dEdxMIP = 2.39761562607903311;
   static Double_t MIPBetaGamma = 4.;
#if 0
   struct Par_t {Int_t h, N; Double_t xmin, xmax, pars[10];};
   const Par_t Par[9] = {
      /* name          h n+1   xmin   xmax      pars[n+1] */
      /* electron */{  0,  4,   3.0,   6.0,{     0.14105,    -0.09078,     0.01901,    -0.00128,           0,           0,           0,           0,           0,           0}},
      /*     muon */{  1,  2,   0.0,   4.5,{    -0.00689,     0.00256,           0,           0,           0,           0,           0,           0,           0,           0}},
      /*     pion */{  2, -1,   0.0,   4.5,{     0.00000,           0,           0,           0,           0,           0,           0,           0,           0,           0}},
      /*     kaon */{  3,  9,  -0.1,   3.7,{     0.00869,     0.21918,    -0.88919,     1.30023,    -0.97075,     0.41298,    -0.10214,     0.01379,    -0.00079,           0}},
      /*   proton */{  4,  9,  -0.6,   3.3,{     0.03052,    -0.02423,    -0.05636,    -0.11585,     0.41292,    -0.38956,     0.16837,    -0.03494,     0.00283,           0}},
      /* deuteron */{  5,  8,  -1.0,   2.9,{     0.03523,    -0.10625,     0.04182,     0.07820,    -0.03816,    -0.02735,     0.01940,    -0.00304,           0,           0}},
      /*   triton */{  6, 10,  -1.0,   2.8,{     0.03092,    -0.07846,     0.01668,     0.00331,     0.12771,    -0.05480,    -0.13897,     0.13928,    -0.04715,     0.00555}},
      /*      He3 */{  7,  9,  -0.8,   2.9,{     0.09158,    -0.07816,     0.07039,     0.00578,    -0.16160,     0.26547,    -0.18728,     0.05988,    -0.00710,           0}},
      /*    alpha */{  8, 10,  -0.8,   2.9,{     0.09366,    -0.08276,     0.06191,     0.02631,    -0.17044,     0.30536,    -0.26867,     0.11847,    -0.02505,     0.00201}}
   };
   //  const Double_t ppar[7]     = { 0.0857988,   9.46138, 0.000206655,     2.12222,       0.974,    -1, 0.13957}; /* pion */
   static Double_t ppar[7]     = { 0.0857988,   9.46138, 0.000206655,     2.12222,       0.974,    -1, 0.13957}; /* pion */
#else
   //  static Double_t ppar[7]     = { 0.0762,  10.632, 0.134e-4, 1.863,  1.948,    -1, -1}; /* Aleph parameters from Alice TPC TDR */
   //  static Double_t ppar[7] = { 0.08942,     8.91971,     0.00024,     2.27383,     1.54174,    -1.00000, 0}; /* pion */
   static Double_t ppar[7]     = {0.12337,     6.61371,     0.00201,     2.27381,     1.54174,    -1.00000, 0 }; /* g All */
#endif
   static Double_t Norm = dEdxMIP/aleph70P(&MIPBetaGamma,ppar);
   Int_t hyp = (Int_t ) par[0];
   Int_t h = Index[hyp];
   Double_t ScaleL10 = 0;
   if (h < 0) {
      h = -h;
      ScaleL10 = TMath::Log10(2.);
   }
   Double_t pove   = TMath::Power(10.,x[0]);
   Double_t mass = Masses[hyp];
   Double_t poverm = pove/mass; 
   Double_t charge = 1.;
   if (h > 6) {
      charge = 2;
      poverm *= charge;
   }
   Double_t bg = poverm;
   /* 
       W.Blum, L. Rolandi "Particle Detection with Drift Chambers", page 246, eq. (9.5)
       F_g(v) = p[0]/beta**p[3]*(p[1] - beta**p[3] - log(p[2] + (beta*gamma)**-p[4]));
       F_g(v) = p[0]*(1/beta**p[3]*(p[1] - log(p[2] + 1/(beta*gamma)**p[4])) - 1) 
   */
   Double_t dEdxL10 =  TMath::Log10(Norm*aleph70P(&bg,ppar));
#if 0
   if (Par[h].N > 0) {
      TString fName(Form("pol%i",Par[h].N-1));
      TF1 *f = (TF1 *) gROOT->GetListOfFunctions()->FindObject(fName);
      if (! f) {
         f = new TF1(fName,fName,0,1);
      }
      f->SetParameters(&Par[h].pars[0]);
      Double_t bgL10 = TMath::Log10(bg);
      dEdxL10 += f->Eval(bgL10);
   }
#endif
   return 2*TMath::Log10(charge) + dEdxL10 + ScaleL10;
}
#endif /* __CINT__ */

void bichselG10(TString input) {  
   if (gClassTable->GetID("StBichsel") < 0 || !m_Bichsel)
   {
      gSystem->Load("libTable");
      gSystem->Load("St_base");
      gSystem->Load("StarClassLibrary");
      gSystem->Load("StBichsel");
      m_Bichsel = Bichsel::Instance();
   }
   const Char_t *type="Bz";
   //TString inputFileLocation = "/gpfs01/star/pwg/truhlar/Run17_P20ic/" + input + "/28FC38D8C6BCFC2677ACA1B26D01EC8D_985.root";
   TString inputFileLocation = "/gpfs01/star/pwg/truhlar/Run17_P20ic/" + input + "/merged/StRP_production_0000.root";
   TFile* data = TFile::Open(inputFileLocation, "read");
   if (!data)
   {
      cout<<"Error: cannot open "<<inputFileLocation<<endl;
      return;
   }

   TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
   tree->Draw("TMath::Log10(dEdxInKevCm0):TMath::Log10(momentumInGeV0)>>hDEdx(200,-1,1,100,0.1,2)"); 
   tree->Draw("TMath::Log10(dEdxInKevCm1):TMath::Log10(momentumInGeV1)>>+hDEdx"); 
   TH2D* hdEdx = (TH2D*)gPad->GetPrimitive("hDEdx");


   TCanvas *cCanvas2D = new TCanvas("hDEdx","hDEdx",800,700);
   gPad->SetMargin(0.12,0.12,0.12,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(0);
   gStyle->SetOptDate(0);
   gStyle->SetLineWidth(2);      //axis line
   gStyle->SetFrameLineWidth(2); //frame line
   gPad->SetTickx();
   gPad->SetTicky(); 
   gStyle->SetOptStat("");
   cCanvas2D->SetGridx(0);
   cCanvas2D->SetGridy(0);
   hdEdx->GetXaxis()->SetTitleFont(42);
   hdEdx->GetYaxis()->SetTitleFont(42);
   hdEdx->GetZaxis()->SetTitleFont(42);
   hdEdx->GetXaxis()->SetLabelFont(42);
   hdEdx->GetYaxis()->SetLabelFont(42);
   hdEdx->GetZaxis()->SetLabelFont(42);
   hdEdx->GetYaxis()->SetRangeUser(0.1,1.7);
   hdEdx->GetXaxis()->SetRangeUser(-0.8,1.0);
   hdEdx->GetXaxis()->SetLabelSize(0.05);
   hdEdx->GetYaxis()->SetLabelSize(0.05);
   hdEdx->GetZaxis()->SetLabelSize(0.05);
   hdEdx->GetXaxis()->SetTitleSize(0.05);
   hdEdx->GetYaxis()->SetTitleSize(0.05);
   hdEdx->GetZaxis()->SetTitleSize(0.05);
   hdEdx->GetXaxis()->SetTitleOffset(1.07);
   hdEdx->GetYaxis()->SetTitleOffset(1.0);
   hdEdx->GetZaxis()->SetTitleOffset(0.50);
   hdEdx->SetStats(0); 
   hdEdx->SetTitle(" ; log_{10} p [GeV] ;log_{10} dE/dx [keV/cm] ");
   hdEdx->Draw("colz");
   cCanvas2D->SetLogz(1);
   
   TPaveText *textPub = new TPaveText(0.15,0.88,0.7,0.95,"brNDC");
   textPub -> SetTextSize(0.05);
   textPub -> SetFillColor(0);
   textPub -> SetTextFont(42);
   textPub -> SetTextAlign(12);
   textPub -> AddText("p + p #rightarrow p + h^{+}h^{-} + p     #sqrt{s} = 510 GeV");
   textPub -> Draw("same");

   TString Type(type);
   TLegend *leg = new TLegend(0.65,0.55,0.8,0.8);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetTextFont(42);
   Double_t xmax = 0.9;
   //  for (int h = 0; h < NMasses; h++) { // Masses
   for (int h = 0; h < NF; h++) 
   { // Masses
      Int_t f = 3;
      if      (Type.Contains("BzM",TString::kIgnoreCase))  f = 7;
      else if (Type.Contains("Bz",TString::kIgnoreCase))   f = 2;
      else if (Type.Contains("I70M",TString::kIgnoreCase)) f = 5;
      else if (Type.Contains("I70",TString::kIgnoreCase))  f = 3;
      else if (Type.Contains("I60",TString::kIgnoreCase))  f = 4;
      else if (Type.Contains("N",TString::kIgnoreCase))    f = 6;

      Int_t dx = 1;
      Char_t *FunName = Form("%s%s%i",FNames[f],Names[h],(int)log2dx[dx]);
      cout << "Make " << FunName << endl;
      Double_t xmin = -0.9;
      //    if (h == 0 || h >= 5) xmin = -0.75;
      if (h == 1) xmin = -0.80;
      if (h == 2) xmin = -0.60;
      if (h == 3) xmin = -0.30;
      TF1 *func = 0;
      if      (f == 3) func = new TF1(FunName,bichsel70,xmin, xmax,2);
      else if (f == 2) func = new TF1(FunName,bichselZ ,xmin, xmax,2);
      else if (f == 5) func = new TF1(FunName,bichsel70M ,xmin, xmax,2);
      else if (f == 6) func = new TF1(FunName,dNdx ,xmin, xmax,2);
      else if (f == 7) func = new TF1(FunName,bichselZM,xmin, xmax,2);
      else { return;}

      if (h == 9) func->SetParameter(0,-Masses[h]);
      else       func->SetParameter(0,Masses[h]);

      func->SetParameter(1,1.);

      if (h >= 7 && h < 9) func->SetParameter(1,2.);

      Int_t color = Colors[h];
      #if 1
         func->SetLineColor(color);
         func->SetMarkerColor(color);
         func->SetLineStyle(1);
         func->SetLineWidth(4);
      #endif
      func->Draw("same");
      leg->AddEntry(func,Names[h]);
      #if !defined( __CINT__) && defined(__Aleph__)
         TF1 *fA = new TF1(Form("Aleph%s",Names[h]),aleph70,xmin,xmax, 1);
         fA->SetParameter(0,h);
         fA->SetLineColor(color);
         fA->SetMarkerColor(color);
         fA->SetLineStyle(4);
         fA->Draw("same");
         leg->AddEntry(fA,Names[h]);
      #endif
   }
   leg->Draw("same");

   cCanvas2D->Update();
   cCanvas2D->SaveAs("dEdx.png");
}
