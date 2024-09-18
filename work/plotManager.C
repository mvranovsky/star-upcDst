#include "include/Libreries.h"
#include "include/PlotUtil.h"

using namespace std;

void plotManager( TString inputFile = "/gpfs01/star/pwg/mvranovsk/Run17_P20ic/AnaV0_23.3.24/merged/StRP_production_0000.root", 
   TString embedFile ="/gpfs01/star/pwg/truhlar/Embedding/embed2510/merged/StRP_production_0000.root"){    

   vector<TString> sourceFiles;  
   sourceFiles.push_back(TString("src/Util.cxx"));
   sourceFiles.push_back(TString("src/PlotUtil.cxx"));
   sourceFiles.push_back(TString("src/RecTree.cxx"));
   //sourceFiles.push_back(TString("src/PlotProbOfRetainEvent.cxx"));
   //sourceFiles.push_back(TString("src/PlotRPMCPlots.cxx")); 
   //sourceFiles.push_back(TString("src/PlotTrigEff.cxx")); 
   //sourceFiles.push_back(TString("src/PlotEmbedQA.cxx")); 
   //sourceFiles.push_back(TString("src/PlotVertexStudy.cxx")); 
   //sourceFiles.push_back(TString("src/PlotTofQA.cxx"));  
   sourceFiles.push_back(TString("src/PlotMainAna.cxx"));
   sourceFiles.push_back(TString("src/PlotAnaV0.cxx"));
   sourceFiles.push_back(TString("src/PlotAnaBP.cxx"));
   cout<<"Loading libraries..."<<endl;
   for(unsigned int i=0; i<sourceFiles.size(); ++i)  
      gROOT->ProcessLine(".L "+sourceFiles[i]+"+g");
   cout<<"Libraries loaded!"<<endl;

   cout<<"Create plots..."<<endl;
   if( runPlots(inputFile, embedFile) ){
      cout<<"PlotManager finished succesfully!"<<endl;
      cout<<"Check FinalPlots.root for the plots!"<<endl;
   }else{
      cout<<"There is a problem. Solve it! And try it again..."<<endl; 
   }
}
