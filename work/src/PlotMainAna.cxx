#include "../include/PlotMainAna.h"

void runMainAnaPlots()
{
   vector<TString> stratName = { TString("noM2"), TString("noPt"), TString("StrictDEdx") };
   vector<TString> setName = { TString("BeforeTotCharge"), TString("AfterTotCharge"), TString("Exclusive") };
   vector<TString> partName = { TString("All pairs"), TString("#pi^{+}#pi^{-} (dE/dx)"),
                           TString("K^{+}K^{-} (dE/dx)"), TString("p#bar{p} (dE/dx)") };

   vector<int> colorSet = { 4, 209, 2, 1};

   TH1D *hist[ partName.size() ];
   for (unsigned int pid = 0; pid < stratName.size(); ++pid)
      for (unsigned int stat = 0; stat < setName.size(); ++stat)
      {
         CreateCanvas(&canvas, stratName[pid] + "_" + setName[stat]);
         canvas->SetLogy();
         CreateLegend(&legend);
         
         // Load and plot hist
         for (unsigned int part = 0; part < partName.size(); ++part)
         {
            hist[part] = (TH1D*)inFile->Get( Form("PID/hMSquared_%i_%i_%i",pid, stat, part) );
            if(!hist[part]){
               cout<<"Error: cannot loaded in PlotMainAna::runMainAnaPlots()"<<endl; 
               continue;
            } 
            SetHistStyle(hist[part], colorSet[part], 20);
            hist[part]->SetTitle(";m^{2}_{TOF} [GeV^{2}];Number of events");
            if( part == 0)
               hist[part]->Draw("");
            else
               hist[part]->Draw("same");

            legend->AddEntry(hist[part], partName[part],"l");
         }
         legend->Draw("same");

         TLine *mLine = new TLine(0.15,0,0.15,hist[0]->GetMaximum()/2);
         mLine->SetLineStyle(10);
         mLine->SetLineColor(2);
         mLine->SetLineWidth(4);
         mLine->Draw("same");

         mLine = new TLine(0.6,0,0.6,hist[0]->GetMaximum()/2);
         mLine->SetLineStyle(10);
         mLine->SetLineColor(1);
         mLine->SetLineWidth(4);
         mLine->Draw("same");

         canvas->Update();
         canvas->Write();
      }
}


