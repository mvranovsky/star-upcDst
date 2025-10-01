#include "AnaGoodRun.h"

//_____________________________________________________________________________
AnaGoodRun::AnaGoodRun(TFile *outFile): Ana(outFile){}

AnaGoodRun::~AnaGoodRun(){
   //if(mUtil) delete mUtil;
}


void AnaGoodRun::Make(){


   nEventsAll++;
   if(CheckTriggers(&JPSItriggers, mUpcEvt, hTriggerBits)){
      JPsiTrigger = true;
      addJPsiTriggerEvent();
   }

   if(CheckTriggers(&ZBtriggers, mUpcEvt, hTriggerBits)){
      addZBEventAll();
      triggerVetoEfficiency();
   }

   tracksBEMC.clear();
   tracksTPC.clear();

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


}



void AnaGoodRun::Init(){

   mUtil = new Util();
   
   
   if( DEBUG )
   cout<<"AnaGoodRun::Init() called"<<endl;
   
   mOutFile->cd();
   mRecTree = new RecTree(nameOfAnaGoodRunTree, AnaGoodRunTreeBits, false);
   mOutFile->mkdir(nameOfAnaGoodRunDir);
   mOutFile->cd(nameOfAnaGoodRunDir);
   
   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   
   
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

   readLumiFile(true);
   readLumiFile(false);
   
   nTracksBEMC = 0;
   nClustersBEMC = 0;
   nTracksTPC = 0;
   nTracksTOF = 0;
   nEvents = 0;
   tpcEtaSum = 0.0;
   bemcEtaSum = 0.0;
   tpcPhiSum = 0.0;
   bemcPhiSum = 0.0;
   nEventsAll = 0;
   nEventsZBVetoAll = 0;
   nEventsZBVetoPassed = 0;
   JPsiTrigger = false;
   cout << "Finished AnaGoodRun::Init()" << endl;

}



void AnaGoodRun::fillRunTree(){

   cout << "Filling runnumber " << mUpcEvt->getRunNumber() << endl;
   if( nEvents <= 0 || !mUpcEvt ){
      // convert inputPath to int
      // inputPath looks like /path/to/file/RUNNUMBER.root
      // get the last part of the path, which is the file name
      // get the run number from the file name


      mRecTree->setRunNumber(0);
      mRecTree->setAtLeast1JPsiTrigger(0);
      mRecTree->setRPsClose(0);
      mRecTree->setNEventsAll(0);
      mRecTree->setNEventsZBVetoAll(0);
      mRecTree->setNEventsZBVetoPassed(0);
      mRecTree->setNEventsLumiFile(0);
      mRecTree->setNEventsJPsi(0);
      mRecTree->setLuminosity(0.0);
      mRecTree->setLuminosityError(0.0);
      mRecTree->setInstLumi(0.0);
      mRecTree->setNTracksBemc(0.0);
      mRecTree->setNClustersBemc(0.0);
      mRecTree->setNTracksTpc(0.0);
      mRecTree->setNTracksTof(0.0);
      mRecTree->setNVertices(0.0);
      mRecTree->setTpcEtaAverage(0.0);
      mRecTree->setBemcEtaAverage(0.0);
      mRecTree->setTpcPhiAverage(0.0);
      mRecTree->setBemcPhiAverage(0.0);
   }else{

      if(isJPsiTrigger()){
         mRecTree->setAtLeast1JPsiTrigger(1);
      }else{
         mRecTree->setAtLeast1JPsiTrigger(0);
      }

      if(areRPsCloseEnough(mUpcEvt->getRunNumber())){
         mRecTree->setRPsClose(1);
      }else{
         mRecTree->setRPsClose(0);
      }

      mRecTree->setRunNumber(mUpcEvt->getRunNumber());
      mRecTree->setNEventsAll(getNEventsAll());
      mRecTree->setNEventsZBVetoAll( getNEventsZBVetoAll() );
      mRecTree->setNEventsZBVetoPassed(getNEventsZBVetoPassed());
      mRecTree->setNEventsJPsi(getJPsiTriggerEvents());
      mRecTree->setNEventsLumiFile(getNEventsLumiFile(mUpcEvt->getRunNumber()));
      mRecTree->setLuminosity(getLuminosity(mUpcEvt->getRunNumber()));
      mRecTree->setInstLumi(getInstantLuminosity(mUpcEvt->getRunNumber()));
      mRecTree->setLuminosityError(0.0); // not used yet
      mRecTree->setTpcEtaAverage(getAverageTpcEta());
      mRecTree->setBemcEtaAverage(getAverageBemcEta());
      mRecTree->setTpcPhiAverage(getAverageTpcPhi());
      mRecTree->setBemcPhiAverage(getAverageBemcPhi());
      //average values
      mRecTree->setAverageTracksBemc(getAverageBemcTracks());
      mRecTree->setAverageClustersBemc(getAverageBemcClusters());
      mRecTree->setAverageTracksTpc(getAverageTpcTracks());
      mRecTree->setAverageTracksTof(getAverageTofTracks());
      mRecTree->setAverageVertices(getAverageVertices());

      cout << "average bemc tracks: " << getAverageBemcTracks() << endl;
      cout << "average tpc tracks: " << getAverageTpcTracks() << endl;
      cout << "average tof tracks: " << getAverageTofTracks() << endl;
      cout << "veto: " << getNEventsZBVetoPassed() << " / " << getNEventsZBVetoAll() << endl;

      /*
      cout << "Luminosity for this run: " << getLuminosity(mUpcEvt->getRunNumber()) << endl;
      cout << "Instant Luminosity: " << getInstantLuminosity(mUpcEvt->getRunNumber()) << endl;
      cout << "Veto efficiency: " << getVetoEfficiency() << endl;
      cout << "Topology efficiency: " << getTopologyEfficiency() << endl;
      */
   }

   mRecTree->FillRecTree();

}



void AnaGoodRun::triggerVetoEfficiency(){
   //
   
   UShort_t dsm0 = mUpcEvt->getLastDSM0(); //trigger bit Tof Mult
   Bool_t tof_mult = dsm0 & (1 << 5);      //true pokud bit 5 byl 1, false naopak.
   UShort_t dsm1 = mUpcEvt->getLastDSM1(); //trigger bit BBC E
   Bool_t bbcE = dsm1 & (1 << 1);      //true pokud bit 1 byl 1, false naopak.
   UShort_t dsm2 = mUpcEvt->getLastDSM1(); //trigger bit BBC W
   Bool_t bbcW = dsm2 & (1 << 2);      //true pokud bit 2 byl 1, false naopak.

   if((tof_mult == false && bbcW == false && bbcE == false))   
      nEventsZBVetoPassed++;

}



bool AnaGoodRun::areRPsCloseEnough(int mRunNumber)
{

   // check if mRunNumber is in mOffset
   if(offsets[0].find(mRunNumber) == offsets[0].end()){
      //cout << "RunNumber: " << mRunNumber << " not found in offsets file!" << endl;
      return false;
   }

   for(int iRP = 0; iRP < nRomanPots; iRP++){

      double offsetY = offsets[iRP][mRunNumber].Y() + corrections[iRP][mRunNumber].Y();
      if( offsetY > 0.04){
         //cout << "Offset Y is too big for this run: " << offsetY << endl;
         return false;
      }
   }
   return true;
}


void AnaGoodRun::readLumiFile( bool isZB ){
   ifstream lumiFile;
   if(isZB){
      lumiFile.open( "/star/u/mvranovsk/star-upcDst/work/lists/luminosityForZB.list" );
   }else{
      lumiFile.open( "/star/u/mvranovsk/star-upcDst/work/lists/lum_perrun_JPsi_HTTP.txt" );
   }

   if (!lumiFile.is_open() ){
      cerr << "\nERROR in AnaGoodRun::readLumiFile(): Problems with opening a file: " << endl;
      return;
   }


   //18149030(runnumber) 12908182(time0) 12909102(time1) 20931(fillnumber) 4.48836088583608e-08(lumi) 1689030.12500000(prescale) 0.706462022811793(livetime) ZDC-trgonly(zdctrg) 0.996090861818113(bs[0]) 1.58758396379823e-08(bs[1]) 5.47782395356853e-09(bs[2]) 3594(nEvents) 438.630043124484(bs[3]) 1271.2384122(bs[4]) 1.83613529620265e-08(bs[5]) 

   int runNumber, time0, time1, fillNumber;
   double lumi, prescale, livetime;
   unsigned int timeDiff;
   string line, zdctrg;
   Double_t bs[6];
   int nEventsLumiFile;

   while ( getline(lumiFile, line) )
   {
      stringstream ss(line);
      ss >> runNumber >> time0 >> time1 >> fillNumber >> lumi >> prescale >> livetime >> zdctrg >> bs[0] >> bs[1] >> bs[2] >> nEventsLumiFile >> bs[3] >> bs[4] >> bs[5];   
      timeDiff = time1-time0;
      double instantinousLumi = prescale*lumi*1000000/double(livetime*timeDiff);
      if( instantinousLumi < 60 || instantinousLumi > 160 ){
         continue;
      }
      // i am interested in instLumi for ProbRetainEvent, in lumi for actual luminosity measurement
      if( isZB ){
         mInstLumiPerRun[runNumber] = instantinousLumi;
      }else{
         mLumiPerRun[runNumber] = lumi;
         mNEventsLumiFile[runNumber] = nEventsLumiFile;
      }

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
