#include "include/Libreries.h"
#include "include/Util.h"

const double textSize = 0.04;
const double labelSize = 0.05;
const int fontStyle = 42;
TString yAxisTitle = "Normalized counts";

using namespace std;

void SetHistStyle(TH1* hist, Int_t color, Int_t markStyle);
void createCanvas(TCanvas **canvas, TString histName);
void createLegend(TLegend **legend);

//_____________________________________________________________________________
void finalize(const string inFileString) {
   //load to prevent error messages about missing dictionary class
   //gSystem->Load("/star/u/truhlar/star-upcDst/build/libstar-upc.so");
   cout<<"Ahoj"<<endl;
   TFile *infile = TFile::Open(inFileString.c_str(), "UPDATE");
   if(!infile){
      cout<<"ERROR: file "<< inFileString<<" not find..."<<endl;
      return;
   }
   Util *mUtil = new Util();
   cout<<"Util created"<<endl;
/*
   TH1D *hVertexZExtraction[nArms], *hDeltaProjection[nCoordinates-1][nBranches], *hDeltaFitProjection[nRomanPots][nCoordinates-1];
   TH1I *hClusterLength[nRomanPots], *hClusterEnergy[nRomanPots], *hClusterPerPlain[nRomanPots][nPlanes], *hNClusterPerRP[nRomanPots];
*/
   vector<TString> setName = { TString("MC"), TString("MCZB"), TString("DATA") };
   vector<TString> anaName = { TString("PureMc"), TString("MC+ZB"), TString("ElasticAna") };
   vector<int> colorSet = { 4, 2, 1};
   vector<int> markerStyle = { 26, 22, 20};
   TH1D *hist[anaName.size()];

   vector<TString> histNames = { TString("hThetaSumX"), TString("hThetaSumY"), TString("hDeltaTheta"), TString("hDCutR") }; //, TString("")
   vector<TString> histPerRPNames = { TString("hClusterLength"), TString("hClusterEnergy"), TString("hNClusterPerRP") }; //, TString("")

   TCanvas *canvas;
   TLegend *legend;
   cout<<"Draw 1D histos"<<endl;
   // Draw 1D histos
   for (unsigned int iHist = 0; iHist < histNames.size(); ++iHist)
   {
      createCanvas(&canvas,histNames[iHist]);
      canvas->SetLogy();
      createLegend(&legend);
      double yMax = 0;
      for (unsigned int iSet = 0; iSet < anaName.size(); ++iSet)
      {
         hist[iSet] = (TH1D*)infile->Get( anaName[iSet] + "/" + histNames[iHist])->Clone("hist" + anaName[iSet]);
         if(!hist[iSet]){
            cout<<"Error hist is not loaded"<<endl; return;
         }
         hist[iSet]->Scale(1.0/double(hist[iSet]->Integral()));
         SetHistStyle(hist[iSet], colorSet[iSet], markerStyle[iSet]);
         if(hist[iSet]->GetMaximum() > yMax)
            yMax = hist[iSet]->GetMaximum();
         //hist[iSet]->SetMaxiumum(0, yMax*1.2); 
         if( iSet == 0)
            hist[iSet]->Draw("hist E");
         else
            hist[iSet]->Draw("same hist E");
      
         legend->AddEntry(hist[iSet], setName[iSet],"ple");
      }

      legend->Draw("same");
      cout<<"Write: "<<histNames[iHist]<<endl;
      canvas->Update();
      canvas->Write(histNames[iHist]);
   }

   cout<<"Debug"<<endl;
   /*
  // Draw 1D histos[nRomanPots]
   int RPordered[] = { 2, 0, 5, 7, 3, 1, 4, 6};
   int rp;
   for (unsigned int iHist = 0; iHist < histPerRPNames.size(); ++iHist)
   {
      createCanvas(&canvas,histPerRPNames[iHist]);
      canvas->Divide(4,2);
      for (int iRp = 0; iRp < nRomanPots; ++iRp)
      {
         canvas->cd(iRp+1);
         canvas->SetLogy();
         createLegend(&legend);
         double yMax = 0;
         rp = RPordered[iRp];
         for (unsigned int iSet = 0; iSet < anaName.size(); ++iSet)
         {
            cout<<"Getting: "<<anaName[iSet] + "/" + histNames[iHist] + "_" + mUtil->rpName(rp)<<endl;
            hist[iSet] = (TH1D*)infile->Get( anaName[iSet] + "/" + histNames[iHist] + "_" + mUtil->rpName(rp))->Clone("hist" + anaName[iSet]);
            if(!hist[iSet]){
               cout<<"Error hist is not loaded"<<endl; return;
            }
            hist[iSet]->Scale(1.0/double(hist[iSet]->Integral()));
            SetHistStyle(hist[iSet], colorSet[iSet], markerStyle[iSet]);
            if(hist[iSet]->GetMaximum() > yMax)
               yMax = hist[iSet]->GetMaximum();
            //hist[iSet]->SetMaxiumum(0, yMax*1.2); 
            if( iSet == 0)
               hist[iSet]->Draw("hist E");
            else
               hist[iSet]->Draw("same hist E");
         
            legend->AddEntry(hist[iSet], setName[iSet],"ple");
         }

         legend->Draw("same");
      }
      canvas->Update();
      canvas->Write(histPerRPNames[iHist]);
   }
   */
}

void createCanvas(TCanvas **canvas, TString histName)
{
   *canvas = new TCanvas(histName,histName,800,700);
   gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetTickx();
   gPad->SetTicky(); 
   gPad->SetLogy(0);
}

void SetHistStyle(TH1* hist, Int_t color, Int_t markStyle)
{
   hist->GetYaxis()->SetTitle(yAxisTitle);
   hist->SetStats(false);
   hist->GetXaxis()->SetTitleFont(fontStyle);
   hist->GetXaxis()->SetTitleFont(fontStyle);
   hist->GetXaxis()->SetLabelFont(fontStyle);
   hist->GetYaxis()->SetLabelFont(fontStyle);
   hist->GetXaxis()->SetLabelSize(labelSize);
   hist->GetYaxis()->SetLabelSize(labelSize);
   hist->GetXaxis()->SetTitleSize(labelSize);
   hist->GetYaxis()->SetTitleSize(labelSize);
   hist->GetXaxis()->SetTitleOffset(0.9);
   hist->GetYaxis()->SetTitleOffset(1.3);
   //hist->GetYaxis()->SetRangeUser(0, 0.4); 
   hist->SetLineColor(color);
   hist->SetLineStyle(1);
   hist->SetLineWidth(1);  
   hist->SetMarkerSize(1);
   hist->SetMarkerColor(color);
   hist->SetMarkerStyle(markStyle);
}//SetGraphStyle

void createLegend(TLegend **legend)
{
   *legend = new TLegend(0.74, 0.74, 0.97, 0.89);
   (*legend)->SetFillStyle(0);
   (*legend)->SetBorderSize(0);
   (*legend)->SetTextSize(textSize);
   (*legend)->SetTextFont(42);
   (*legend)->SetMargin(0.1);   
}

