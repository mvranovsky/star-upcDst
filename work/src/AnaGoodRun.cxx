#include "AnaGoodRun.h"

//_____________________________________________________________________________
AnaGoodRun::AnaGoodRun(TFile *outFile): Ana(outFile){}

AnaGoodRun::~AnaGoodRun(){
   //if(mUtil) delete mUtil;
}


void AnaGoodRun::Make(){

   tracksBEMC.clear();
   tracksTPC.clear();


   //trigger
   hAnalysisFlow->Fill(GRALL);

   
   if(CheckTriggers(&ZBtriggers,mUpcEvt,hTriggerBits)){

      nEventsAll++;
      
      
      UShort_t dsm0 = mUpcEvt->getLastDSM0(); //trigger bit Tof Mult
      Bool_t tof_mult = dsm0 & (1 << 5);      //true pokud bit 5 byl 1, false naopak.
      UShort_t dsm1 = mUpcEvt->getLastDSM1(); //trigger bit BBC E
      Bool_t bbcE = dsm1 & (1 << 1);      //true pokud bit 1 byl 1, false naopak.
      UShort_t dsm2 = mUpcEvt->getLastDSM1(); //trigger bit BBC W
      Bool_t bbcW = dsm2 & (1 << 2);      //true pokud bit 2 byl 1, false naopak.
      
      if((tof_mult == false && bbcW == false && bbcE == false))   nEventsPassed++;
      
   }
   
   if(CheckTriggers(&JPSItriggers, mUpcEvt, hTriggerBits)){
      JPsiTrigger = true;
      addJPsiTriggerEvent();
   }

   
   // track loop
   int nTofTracks = 0;
   vector<Double_t> tpcEta, tpcPhi, tpcNHitsFit, tpcNHitsDEdx, tpcNSigmaP, tpcNSigmaPi, tpcNSigmaK, tpcNSigmaE, bemcEta, bemcPhi;
   for (int iTrk = 0; iTrk < mUpcEvt->getNumberOfTracks(); ++iTrk){
      const StUPCTrack *trk = mUpcEvt->getTrack(iTrk);

      if(!trk->getFlag(StUPCTrack::kPrimary)) // original star-upcDst
         continue;

      if(!goodQualityTrack(trk))
         continue;

      
      tracksTPC.push_back(iTrk);

      hAnalysisFlow->Fill(GRGOODTRACKTPC);

      tpcEta.push_back(trk->getEta());
      tpcPhi.push_back(trk->getPhi());
      tpcNHitsFit.push_back(trk->getNhitsFit());
      tpcNHitsDEdx.push_back(trk->getNhitsDEdx());
      tpcNSigmaP.push_back(trk->getNSigmasTPCProton());
      tpcNSigmaPi.push_back(trk->getNSigmasTPCPion());
      tpcNSigmaK.push_back(trk->getNSigmasTPCKaon());
      tpcNSigmaE.push_back(trk->getNSigmasTPCElectron());


      if( !trk->getFlag(StUPCTrack::kTof) )
         nTofTracks++;
      if( !trk->getFlag(StUPCTrack::kBemc) )
         continue;

      hAnalysisFlow->Fill(GRGOODTRACKBEMC);
      tracksBEMC.push_back(iTrk);
      bemcEta.push_back(trk->getBemcEta());
      bemcPhi.push_back(trk->getBemcPhi());
   }

   if(tracksTPC.size() == 0){
      //cout << "Found " << tracksBEMC.size() << " tracks in BEMC. Leaving event." << endl;
      return;
   }

   addAverageBemcTracks(tracksBEMC.size());
   addAverageTpcTracks(tracksTPC.size());
   addAverageTOFTracks(nTofTracks);
   addAverageBemcClusters(mUpcEvt->getNumberOfClusters());
   addAverageVertices(mUpcEvt->getNumberOfVertices());
   addEvent();


   for(auto &trk : tracksTPC){
      const StUPCTrack *track = mUpcEvt->getTrack(trk);
      if(!track)  continue;

      addTpcEtaValue(track->getEta());
      addTpcPhiValue(track->getPhi());
      hTpcEta->Fill(track->getEta());
      hTpcPhi->Fill(track->getPhi());
      hTpcEtaPhi->Fill(track->getEta(), track->getPhi());
   }

   for(auto &trk : tracksBEMC){
      const StUPCTrack *track = mUpcEvt->getTrack(trk);
      if(!track)  continue;

      addBemcEtaValue(track->getBemcEta());
      addBemcPhiValue(track->getBemcPhi());
      hBemcEta->Fill(track->getBemcEta());
      hBemcPhi->Fill(track->getBemcPhi());
      hBemcEtaPhi->Fill(track->getBemcEta(), track->getBemcPhi());
   }

   
   if(DEBUG){
      cout << "Finished AnaGoodRun::Make()" << endl;
   }
}



void AnaGoodRun::Init(){

   mUtil = new Util();
   
   
   if( DEBUG )
     cout<<"AnaGoodRun::Init() called"<<endl;
     
     mOutFile->cd();
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nGRCuts-1, 1, nGRCuts);
   for(int tb=1; tb<nGRCuts; ++tb) {
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisGoodRun(tb));
   }
   
   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   
   //mRecTree = new RecTree(nameOfAnaGoodRunTree, AnaGoodRunTreeBits, false); 

   hTpcEtaPhi = new TH2F("hTpcEtaPhi", "hTpcEtaPhi; #eta; #phi", 100, -1.0, 1.0, 100, -3.14, 3.14);
   hBemcEtaPhi = new TH2F("hBemcEtaPhi", "hBemcEtaPhi; #eta; #phi", 100, -1.0, 1.0, 100, -3.14, 3.14);
   hTpcEta = new TH1D("hTpcEta", "hTpcEta; #eta; counts", 100, -1.5, 1.5);
   hBemcEta = new TH1D("hBemcEta", "hBemcEta; #eta; counts", 100, -1.5, 1.5);
   hTpcPhi = new TH1D("hTpcPhi", "hTpcPhi; #phi; counts", 100, -3.14, 3.14);
   hBemcPhi = new TH1D("hBemcPhi", "hBemcPhi; #phi; counts", 100, -3.14, 3.14);


   // load offset file
   if(!LoadOffsetFile(offsetFilePath, offsets )){
      cout << "Could not load file with offsets!" << endl;
      return;
   }

   // load offset correction file 
   if(!LoadOffsetFile(offsetCorrectionsFilePath, corrections ) ){
      cout << "Could not load file with offset corrections!" << endl;
      return;
   }

   readLumiFile();
   
   nTracksBEMC = 0;
   nClustersBEMC = 0;
   nTracksTPC = 0;
   nTracksTOF = 0;
   nEvents = 0;
   nVertices = 0;
   tpcEtaSum = 0.0;
   bemcEtaSum = 0.0;
   tpcPhiSum = 0.0;
   bemcPhiSum = 0.0;
   nEventsAll = 0;
   nEventsPassed = 0;
   JPsiTrigger = false;
   cout << "Finished AnaGoodRun::Init()" << endl;

}

bool AnaGoodRun::goodQualityTrack(const StUPCTrack *trk){

   if( !(trk->getDcaZ() > minDcaZ && trk->getDcaZ() < maxDcaZ) )  //DCA z
      return false;

   if( !(trk->getDcaXY() < maxDcaXY) ) //DCA xy
      return false;

   if( !(trk->getNhitsFit() > minNHitsFit) )  //NhitsFit
      return false;

   if( !(trk->getNhitsDEdx() > minNHitsDEdx) ) //NhitsdEdx
      return false;

   return true;

}

void AnaGoodRun::readLumiFile(){
    ifstream lumiFile;

    lumiFile.open( "/star/u/truhlar/star-upcDst/work/lists/luminosityForZB.list" );
    if (!lumiFile.is_open() )
    {
       cerr << "\nERROR in PlotManager::readLumiFile(): Problems with opening a file: " << endl;
       return;
    }

    int runNumber, time0, time1, fillNumber;
    double lumi, prescale, livetime;
    unsigned int timeDiff;
    string line;

    while ( getline(lumiFile, line) )
    {
       stringstream ss(line);
       ss >> runNumber >> time0 >> time1 >> fillNumber >> lumi >> prescale >> livetime;
       timeDiff = time1-time0;
       double instantinousLumi = prescale*lumi*1000000/double(livetime*timeDiff);
       if( instantinousLumi < 60 || instantinousLumi > 160 ){
          continue;
       }
       mInstLumiPerRun[runNumber] = instantinousLumi;

    }

    lumiFile.close();

    return;
}  




bool AnaGoodRun::LoadOffsetFile(TString fileName, map<unsigned int, TVector3> (&offsets)[nRomanPots]){

   string line;
   ifstream file( fileName );
   if (!file.is_open() )
   {
      cout << "\n ERROR in LoadOffsets(): Problems with opening a file: "<< fileName<< endl;
      return false;
   }   
   int RUNNUMBER;
   double var[2*nRomanPots]; // placeholder for offsets or offsets correction in x and y
   while ( getline(file, line) )
   {
      stringstream ss(line);
      ss >> RUNNUMBER;
      for (int i = 0; i < 2*nRomanPots; ++i)
         ss >> var[i]; 

      for (int iRP = 0; iRP < nRomanPots; ++iRP)
         offsets[iRP][RUNNUMBER].SetXYZ(var[2*iRP], var[2*iRP + 1], 0.0); 
   }

   file.close();
   return true;
}
