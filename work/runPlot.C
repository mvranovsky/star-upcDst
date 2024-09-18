#include "include/PlotUtil.h"

using namespace std;

void runPlots( TString inputFile = "/gpfs01/star/pwg/truhlar/Run17_P20ic/ana3010/merged/StRP_production_0000.root", 
   TString embedFile ="/gpfs01/star/pwg/truhlar/Embedding/embed2510/merged/StRP_production_0000.root"){    

   cout<<"Create plots..."<<endl;
   if( PlotUtil::runPlots(inputFile, embedFile) ){
      cout<<"PlotManager finished succesfully!"<<endl;
      cout<<"Check FinalPlots.root for the plots!"<<endl;
   }else{
      cout<<"There is problem. Solve it! And try it again..."<<endl; 
   }
}
