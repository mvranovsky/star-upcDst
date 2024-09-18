#include "../include/PlotUtil.h"

bool runPlots(TString inputFile, TString embedFile)
{
   if( !Init(inputFile, embedFile) ){
      cout<<"PlotUtil::runPlots() could not init input files..."<<endl;
      return false;
   }

   if( runMAINANA )
      runMainAnaPlots();

   if( runRPMCANA )
      runRPMCPlots();

   if( runFULLZB ){
      runProbOfRetainEvent();
      runDsmEffStudy( mRecTree[kFULLZB], nameOfFullZBSet);
   }

   if( runTRIGEFF ){
      runDsmEffStudy( mRecTree[kTRIGEFF], nameOfTrigEffSet);
      runTofTrigStudy();
      runRpTrigStudy();
   }

   if( runTOFQA )
      runTofQA();

   if( runVERTEXSTUDY )
      runVertexStudy();

   if( runEMBEDQA && runMAINANA )
      runEmbedQA();

   outFile->Close(); 
   if( runEMBEDQA )
      embFile->Close(); 
   inFile->Close();      

   return true;
}


bool ConnectInput( TString inputFile, TFile **file) 
{
   cout << "Input from root file: "<< inputFile << endl;
   *file = TFile::Open(inputFile, "read");
   if(!*file)
   {
      cout<< "PlotUtil::ConnectInput() Couldn't open "<< inputFile <<endl;
      return false;
   } 

   return true;
}//ConnectInput


//_____________________________________________________________________________
TFile *CreateOutputFile(const string& out) {

   TFile *outputFile = TFile::Open(out.c_str(), "recreate");
   if(!outputFile) 
      return 0x0;

   return outputFile;
}//CreateOutputFile


void CreateCanvas(TCanvas **canvas, TString canvasName)
{
   *canvas = new TCanvas(canvasName,canvasName,800,700);
   gPad->SetMargin(0.13,0.03,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetTickx();
   gPad->SetTicky(); 
   //gPad->SetLogy(0);

   // setup the drawing style
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetNumberContours(99);
   gStyle->SetPalette(55);
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
}//SetHistStyle

void SetTH2Style(TH2* hist)
{
   hist->SetStats(false);
   hist->GetXaxis()->SetTitleFont(fontStyle);
   hist->GetYaxis()->SetTitleFont(fontStyle);
   hist->GetZaxis()->SetTitleFont(fontStyle);
   hist->GetXaxis()->SetLabelFont(fontStyle);
   hist->GetYaxis()->SetLabelFont(fontStyle);
   hist->GetZaxis()->SetLabelFont(fontStyle);
   hist->GetXaxis()->SetLabelSize(textSize);
   hist->GetYaxis()->SetLabelSize(textSize);
   hist->GetZaxis()->SetLabelSize(textSize);
   hist->GetXaxis()->SetTitleSize(labelSize);
   hist->GetYaxis()->SetTitleSize(labelSize);
   hist->GetZaxis()->SetTitleSize(labelSize);
   hist->GetXaxis()->SetTitleOffset(0.7);
   hist->GetYaxis()->SetTitleOffset(0.7);
   hist->GetZaxis()->SetTitleOffset(0.9);
}//SetTH2Style

void CreateLegend(TLegend **legend)
{
   *legend = new TLegend(0.74, 0.74, 0.97, 0.89);
   (*legend)->SetFillStyle(0);
   (*legend)->SetBorderSize(0);
   (*legend)->SetTextSize(textSize);
   (*legend)->SetTextFont(42);
   (*legend)->SetMargin(0.1);   
}

void DrawSTARInternal(double xl, double yl, double xr, double yr)
{
   TPaveText *textSTAR;
   textSTAR = new TPaveText(xl, yl, xr, yr,"brNDC");
   textSTAR -> SetTextSize(textSize);
   textSTAR -> SetFillColor(0);
   textSTAR -> SetTextFont(62);
   textSTAR->AddText("STAR Internal");
   textSTAR -> Draw("same");
}

void CreateText(double xl, double yl, double xr, double yr)
{
   text = new TPaveText(xl, yl, xr, yr,"brNDC");
   text -> SetTextSize(textSize);
   text -> SetFillColor(0);
   text -> SetTextFont(fontStyle);
   text -> SetTextAlign(32);   
}

bool Init(TString inputFile, TString embedFile){
   //connect input file
   if(!ConnectInput(inputFile, &inFile))
   {
      cout << "PlotUtil::Init() Could not open file: "<< inputFile <<endl; 
      return false;
   }
   
   if( runEMBEDQA ){
      if(!ConnectInput(embedFile, &embFile))
      {
         cout << "PlotUtil::Init() Could not open file: "<< embedFile <<endl; 
         return false;
      }
      mTree[kEMBEDQA] = dynamic_cast<TTree*>( embFile->Get( nameOfTree[kEMBEDQA] ) );
      if (!mTree[kEMBEDQA])
      {
         cout<<"Error: cannot open "<<nameOfTree[kEMBEDQA]<<endl;
         return false;
      }
      mRecTree[kEMBEDQA] = new RecTree(mTree[kEMBEDQA], treeBits[kEMBEDQA] );      
   }

   //open output file
   outFile = CreateOutputFile("FinalPlots.root"); 
   if(!outFile) 
   {
      cout << "Can not open output file." << endl; 
      return false;
   }

   mUtil = new Util();

   if( runFULLZB )
      if(!readLumiFile()) 
         return false;

   for (int iStudy = 0; iStudy < nStudies; ++iStudy)
   {
      if( iStudy == kEMBEDQA || iStudy == kALIGNMENT )
         continue;
      if( !runStudy[iStudy] || nameOfTree[iStudy]=="")
         continue;

      mTree[iStudy] = dynamic_cast<TTree*>( inFile->Get( nameOfTree[iStudy] ) );
      if (!mTree[iStudy])
      {
         cout<<"Error: cannot open "<<nameOfTree[iStudy]<<endl;
         return false;
      }
      mRecTree[iStudy] = new RecTree(mTree[iStudy], treeBits[iStudy] );
   }   


   return true;
}

 
void SetGPad(double xl, double yl, double xr, double yr)
{
   gPad->SetMargin(xl,yl,xr,yr); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
   gPad->SetTickx();
   gPad->SetTicky(); 
   gPad->SetLogy(0);
   gStyle->SetOptStat("");   
}

void DrawFiducial()
{
   const Int_t n = 100;
   Double_t x[n], y[n];
   Double_t tmp;
   for(int i = 0; i < n; ++i){
      x[i] = 0.175 + (0.27*i)/n;
      //tmp = (x[i] -1.163)*(x[i]-1.163) - 0.464*0.464;
      tmp = (x[i] +0.6)*(x[i]+0.6) - 1.25;
      y[i] = -sqrt(abs(tmp));
   }
   TGraph* gr = new TGraph(n,x,y);
   gr->SetLineWidth(4);
   gr->Draw("same");


   TLine *left02 = new TLine(-0.27,-0.4,0.445,-0.4);
   SetLineStyle(left02);
   left02->Draw("same");


   TLine *left01 = new TLine(-0.27,-0.8,-0.27,-0.4);
   SetLineStyle(left01);
   left01->Draw("same");

   left01 = new TLine(-0.27,-0.8,0.185,-0.8);
   SetLineStyle(left01);
   left01->Draw("same");          
   // UP
   left02 = new TLine(-0.27,0.4,0.445,0.4);
   SetLineStyle(left02);
   left02->Draw("same");


   left01 = new TLine(-0.27,0.4,-0.27,0.8);
   SetLineStyle(left01);
   left01->Draw("same");

   left01 = new TLine(-0.27,0.8,0.185,0.8);
   SetLineStyle(left01);
   left01->Draw("same");

   for(int i = 0; i < n; ++i){
      x[i] = 0.175 + (0.27*i)/n;
      //tmp = (x[i] -1.31)*(x[i]-1.31) - 0.725*0.725;
      tmp = (x[i] +0.6)*(x[i]+0.6) - 1.25;
      y[i] = sqrt(abs(tmp));
   }
   gr = new TGraph(n,x,y);
   gr->SetLineWidth(4);
   gr->Draw("same"); 
}

void SetLineStyle(TLine* line)
{
   line->SetLineStyle(1);
   line->SetLineColor(1);
   line->SetLineWidth(4);
}//SetLineStyle