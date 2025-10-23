// saveResults.C
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TString.h"

void SaveResults(const char* dirName) {
    // Open the input ROOT file
    TFile* inFile = TFile::Open("AnalysisOutput.root", "READ");
    if (!inFile || inFile->IsZombie()) {
        printf("Error: Cannot open AnalysisOutput.root\n");
        return;
    }

    // Get the histogram
    TH1* hist = (TH1*)inFile->Get("hInvMassSpectrum");
    if (!hist) {
        printf("Error: Histogram hInvMassSpectrum not found\n");
    }

    // Get the graph
    TGraph* graph = (TGraph*)inFile->Get("gBEControl");
    if (!graph) {
        printf("Error: Graph gBEControl not found\n");
    }


    // get histograms for bemc efficiency study
    TDirectory* dir = (TDirectory*)inFile->GetDirectory("EmbeddingJPsiPlots");
    TGraphAsymmErrors* gEffEmbed;
    TGraphAsymmErrors* gEffMc;
    if(dir){

        TH1* hBemcPtAllEmbed = (TH1*)dir->Get("hBemcPtAllEmbed");
        TH1* hBemcPtHitEmbed = (TH1*)dir->Get("hBemcPtHitEmbed");
        TH1* hBemcPtAllMc = (TH1*)dir->Get("hBemcPtAllMc");
        TH1* hBemcPtHitMc = (TH1*)dir->Get("hBemcPtHitMc");

        if(hBemcPtAllEmbed && hBemcPtHitEmbed){

            TEfficiency* effEmbed = new TEfficiency(*hBemcPtHitEmbed, *hBemcPtAllEmbed);
            gEffEmbed = (TGraphAsymmErrors*)effEmbed->CreateGraph();
            gEffEmbed->SetName("gEffEmbed");
        }else {
            printf("Error: One or more BEMC histograms not found\n");
        }

        if(hBemcPtAllMc && hBemcPtHitMc){
            TEfficiency* effMc = new TEfficiency(*hBemcPtHitMc, *hBemcPtAllMc);
            gEffMc = (TGraphAsymmErrors*)effMc->CreateGraph();
            gEffMc->SetName("gEffMc");
        } else {
            printf("Error: One or more BEMC histograms not found\n");
        }
    } else {
        printf("Error: Directory EmbeddingJPsiPlots not found\n");

    }




    // Open the output file (update mode so we don’t overwrite)
    TFile* outFile = TFile::Open("BemcStudy.root", "UPDATE");
    if (!outFile || outFile->IsZombie()) {
        printf("Error: Cannot open BemcStudy.root\n");
        inFile->Close();
        return;
    }

    // Create directory
    outFile->cd();
    outFile->mkdir(dirName);
    outFile->cd(dirName);

    // Write objects
    if (hist) hist->Write();
    if (graph) graph->Write();
    if(gEffEmbed) gEffEmbed->Write();
    if(gEffMc) gEffMc->Write();

    // Clean up
    outFile->Close();
    inFile->Close();
}
