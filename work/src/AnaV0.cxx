#include "AnaV0.h"

AnaV0::AnaV0(TFile *outFile): Ana(outFile){}

void AnaV0::Make()
{

   vector<int> tracksV0;
   for (int itrk = 0; itrk < mUpcEvt->getNumberOfTracks(); ++itrk){
       const StUPCTrack* trk = mUpcEvt->getTrack(itrk);
       hAnalysisFlow->Fill(V0ALL);
       if(!CheckTriggers(&triggerID, mUpcEvt, hTriggerBits))
          return;
       hAnalysisFlow->Fill(V0TRIG);
       
       if(TOF2Tracks && !trk->getFlag( StUPCTrack::kTof ) )
        continue;

       //conditions for track quality
       fillGoodTrackCuts(trk);
       if (trk->getNhitsFit() < minNHitsFit)
           continue;
       if (trk->getNhitsDEdx() < minNHitsDEdx)
           continue;
       if (trk->getPt() < minPt)
           continue; 

       //condition for pseudorapidity and filling histos
       Double_t eta = trk->getEta();
       Double_t phi = trk->getPhi();
       hEtaPhi->Fill(eta, phi);
       hEta->Fill(eta);
       if(abs(eta) > maxEta) continue;
       hAnalysisFlow->Fill(V0ETA);


       saveNSigmaCorr(trk);
       //particle identification: proton or pion
       if(!(hasGoodTPCnSigma(trk) == PROTON || hasGoodTPCnSigma(trk) == PION))
           continue;
       hAnalysisFlow->Fill(V0PID);

       if (trk->getFlag( StUPCTrack::kPrimary ) && trk->getFlag( StUPCTrack::kV0 )){
           cout << "We got a track with kPrimary and kV0 flag." << endl;
           hBothFlags->Fill(1);
       }

       if( trk->getFlag( StUPCTrack::kV0 ) && !V0Control){
           tracksV0.push_back(itrk);
       } else if(trk->getFlag( StUPCTrack::kPrimary ) && V0Control ){
           tracksV0.push_back(itrk);
       } else {
           continue; 
       }

       hAnalysisFlow->Fill(V0FLAG);
   }//global tracks loop

   mRecTree->setNGoodTpcTrks( tracksV0.size() );
   SaveEventInfo(mUpcEvt);
   SaveTriggerInfo(mUpcEvt, mRpEvt);

   //info about the beamline and mag field
   double bField = mUpcEvt->getMagneticField();
   double beamline[4]; 
   beamline[0] = mUpcEvt->getBeamXPosition();
   beamline[2] = mUpcEvt->getBeamXSlope();
   beamline[1] = mUpcEvt->getBeamYPosition();
   beamline[3] = mUpcEvt->getBeamYSlope();


   TLorentzVector state, hadron1, hadron2;
   TVector3 vertex, vertex0, primaryVertex;
   vertex0 = {0,0,0};
   int HADRON1, HADRON2, totalCharge;
   int vertexIdTrk1, vertexIdTrk2;
   Bool_t tofMatchTrk1, tofMatchTrk2;
   Double_t vertexDiff;
   vector<pair<int, int>> hadronPairV0;
   //first double loop is over all the V0 candidates
   //outer tracks loop
   for (unsigned int itrk = 0; itrk < tracksV0.size(); ++itrk){
       const StUPCTrack* trk1 = mUpcEvt->getTrack(tracksV0[itrk]);
       tofMatchTrk1 = trk1->getFlag( StUPCTrack::kTof );
       vertexIdTrk1 = trk1->getVertexId();
       HADRON1 = hasGoodTPCnSigma(trk1);

       //inner loop of tracks
       for (unsigned int jtrk = itrk+1; jtrk < tracksV0.size(); ++jtrk){

           const StUPCTrack* trk2 = mUpcEvt->getTrack(tracksV0[jtrk]);
           tofMatchTrk2 = trk2->getFlag( StUPCTrack::kTof );
           vertexIdTrk2 = trk2->getVertexId();
           HADRON2 = hasGoodTPCnSigma(trk2);

           // checking for at least 1 tof match
           if(!tofMatchTrk1 && !tofMatchTrk2)
              continue;
           // checking for wrong identification
           if (HADRON1 == PROTON && HADRON2 == PROTON)
               continue; 

            //if both tracks originate from a upcDst vertex, calculate the difference between prodvertexhypo and upcDst vtx
            const StUPCVertex* vtx = mUpcEvt->getVertex(vertexIdTrk1);
            if(vertexIdTrk1 != vertexIdTrk2)
                continue;
            if(vtx){
                primaryVertex = {vtx->getPosX(),vtx->getPosY(), vtx->getPosZ()};
            } else continue;
            //hAnalysisFlow->Fill(V0SAMEVTX);


           //define V0 pair
           StUPCV0 K0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), tracksV0[itrk], tracksV0[jtrk], primaryVertex, beamline, bField, false);
           vertex = K0.prodVertexHypo();

           //same cuts as were in StUPCSelectV0 except K0L1, K0L2 -> K0
           fillTopologyCutsBefore(K0);
           if ( !(K0.dcaDaughters() < maxDcaDaughters && K0.DCABeamLine() < maxDcaBeamLine && K0.pointingAngleHypo()> minPointingAngle && K0.decayLengthHypo()< maxDecayLengthHypo)  )
              continue;
           fillTopologyCutsAfter(K0);
           hAnalysisFlow->Fill(V0PAIR);



            // condition for difference between 
            vertexDiff = abs(vtx->getPosZ() - vertex.Z());
            hVtxDiff->Fill(vertexDiff);
            if(vertexDiff > vertexDiffMax)
                continue;

           //cut for z position of vertex- using prodvertexhypo
           hPosZ->Fill(vertex.Z());
           if (abs(vertex.Z()) > vertexRange )
               continue;
           hPosZCut->Fill(vertex.Z()); 
           hAnalysisFlow->Fill(V0ZVERTEX);
           hadronPairV0.push_back(make_pair(tracksV0[itrk],tracksV0[jtrk]));
       }
   }

   if (hadronPairV0.size() != 0)
       hadronPairV0 = filterPairs(hadronPairV0);

   //tree can only hold max 5 states per event
   hNPairV0->Fill(hadronPairV0.size());

   if (hadronPairV0.size() > 5){
       //cout<< "too many pairs. Number of pairs: " << hadronPairV0.size() << endl; 
       hadronPairV0.resize(5);
   }
   int j = 0;
   for (int i = 0; i < hadronPairV0.size(); ++i){
        //checking for pair made of one track
       if (hadronPairV0[i].first == hadronPairV0[i].second){
          hSameTrackPair->Fill(1);
          cout << "We got a pair with same track for both. " << endl;     
       }

      const StUPCTrack* trk1 = mUpcEvt->getTrack(hadronPairV0[i].first);
      const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronPairV0[i].second);

      HADRON1 = hasGoodTPCnSigma(trk1);
      HADRON2 = hasGoodTPCnSigma(trk2);

      hTotQ->Fill( trk1->getCharge() + trk2->getCharge() );
      if (oppositePair(trk1, trk2))
         hAnalysisFlow->Fill(V0OPPOSITE);

      trk1->getLorentzVector(hadron1, mUtil->mass(HADRON1));
      trk2->getLorentzVector(hadron2, mUtil->mass(HADRON2));
      state = hadron1 + hadron2;

      totalCharge = trk1->getCharge() + trk2->getCharge();
           
      if (HADRON1 == HADRON2){
         mRecTree->setPairID(K0S,j);
         if(totalCharge == 0)
           hAnalysisFlow->Fill(V0PIPI);
      } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && lambda(trk1, trk2) ){ 
         mRecTree->setPairID(LAMBDA,j);
         if(totalCharge == 0)
           hAnalysisFlow->Fill(V0PPI);
      } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && !lambda(trk1, trk2) ){
         mRecTree->setPairID(LAMBDABAR,j);
         if(totalCharge == 0)
           hAnalysisFlow->Fill(V0PIPBAR); 
      } else continue;

      StUPCV0 V0(trk1, trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronPairV0[i].first,hadronPairV0[i].second, vertex0, beamline, bField, false);  
      vertex = V0.prodVertexHypo();
      
      // define vertex difference, if possible
      vertexIdTrk1 = trk1->getVertexId();
      vertexIdTrk2 = trk2->getVertexId();
      const StUPCVertex* vtx = mUpcEvt->getVertex(vertexIdTrk1);
      if(vertexIdTrk1 == vertexIdTrk2 && vtx){
        vertexDiff = vertex.Z() - vtx->getPosZ();
      } else{
        vertexDiff = -999;
      }


      // saving info for tracks and vertices
      SaveVertexInfo(vertex, V0.dcaDaughters(), V0.DCABeamLine(), V0.pointingAngleHypo(), V0.decayLengthHypo(), vertexDiff , j);
      SaveTrackInfo(trk1,j*2);
      SaveTrackInfo(trk2,j*2 +1);
      
        // save info about tracks-> the first one always has match with ToF 
      if( trk1->getFlag( StUPCTrack::kTof ) ){
          //trk1 = tag
          SaveTrackInfo(trk1,2*j);
          fillTag(trk1);
          //trk2 = probe
          SaveTrackInfo(trk2,2*j+1);
          fillProbe(trk2);
      } else{
          //trk2 = tag
          SaveTrackInfo(trk2,2*j);  
          fillTag(trk2);
          //trk1 = probe
          SaveTrackInfo(trk1,2*j+1);
          fillProbe(trk1);  
      }
      // save info about state
      SaveStateInfo(state,totalCharge, j);
      ++j;
   }
   //cerr << "Got through the whole thing" << endl;
   if(mRecTree->getInvMass(0) == -9999){
    resetInfo();
    return;
   }

   // save event info
   mRecTree->FillRecTree();

   // create an Armenteros-Podolanski plot
   fillFinalPlots();



   // end with clearing all the variables for rec tree 
   resetInfo();
}//end of make()

void AnaV0::Init()
{
    //initialize class for selecting V0 candidates 
    //mSelectV0 = new StUPCSelectV0Modified();

    //init database to read beam-line parameters
    //mDbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    //mDbMk->Init();

   mUtil = new Util();


   if( DEBUG )
      cout<<"AnaV0::Init() called"<<endl;

   mOutFile->cd();
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nV0SelectionCuts-1, 1, nV0SelectionCuts);
   for(int tb=1; tb<nV0SelectionCuts; ++tb) {
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisV0SelectionName(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   mRecTree = new RecTree(nameOfAnaV0Tree, AnaV0TreeBits, false); 

   mOutFile->mkdir("PID")->cd();

   hNSigmaPiPcorr = new TH2F("hNSigmaPiPcorr","n_{#sigma} pions against protons", 100, -25, 25, 100, -25, 25);
   hNSigmaPiPcorr->GetXaxis()->SetTitle("n#sigma_{#pi}");
   hNSigmaPiPcorr->GetYaxis()->SetTitle("n#sigma_{p}");

   hNSigmaPiKcorr = new TH2F("hNSigmaPiKcorr","n_{#sigma} pions against kaons", 100, -25, 25, 100, -25, 25);
   hNSigmaPiKcorr->GetXaxis()->SetTitle("n#sigma_{#pi}");
   hNSigmaPiKcorr->GetYaxis()->SetTitle("n#sigma_{K}");

   hNSigmaPiecorr = new TH2F("hNSigmaPiecorr","n_{#sigma} pions against electrons", 100, -25, 25, 100, -25, 25);
   hNSigmaPiecorr->GetXaxis()->SetTitle("n#sigma_{#pi}");
   hNSigmaPiecorr->GetYaxis()->SetTitle("n#sigma_{e}");

   hNSigmaPPicorr = new TH2F("hNSigmaPPicorr","n_{#sigma} protons against pions", 100, -25, 25, 100, -25, 25);
   hNSigmaPPicorr->GetXaxis()->SetTitle("n#sigma_{p}");
   hNSigmaPPicorr->GetYaxis()->SetTitle("n#sigma_{#pi}");

   hNSigmaPKcorr = new TH2F("hNSigmaPKcorr","n_{#sigma} protons against kaons", 100, -25, 25, 100, -25, 25);
   hNSigmaPKcorr->GetXaxis()->SetTitle("n#sigma_{p}");
   hNSigmaPKcorr->GetYaxis()->SetTitle("n#sigma_{K}");

   hNSigmaPecorr = new TH2F("hNSigmaPecorr","n_{#sigma} protons against electrons", 100, -25, 25, 100, -25, 25);
   hNSigmaPecorr->GetXaxis()->SetTitle("n#sigma_{p}");
   hNSigmaPecorr->GetYaxis()->SetTitle("n#sigma_{e}");

   hNSigmaKPicorr = new TH2F("hNSigmaKPicorr","n_{#sigma} kaons against pions", 100, -25, 25, 100, -25, 25);
   hNSigmaKPicorr->GetXaxis()->SetTitle("n#sigma_{K}");
   hNSigmaKPicorr->GetYaxis()->SetTitle("n#sigma_{#pi}");

   hNSigmaKPcorr = new TH2F("hNSigmaKPcorr","n_{#sigma} kaons against protons", 100, -25, 25, 100, -25, 25);
   hNSigmaKPcorr->GetXaxis()->SetTitle("n#sigma_{K}");
   hNSigmaKPcorr->GetYaxis()->SetTitle("n#sigma_{p}");

   hNSigmaKecorr = new TH2F("hNSigmaKecorr","n_{#sigma} kaons against electrons", 100, -25, 25, 100, -25, 25);
   hNSigmaKecorr->GetXaxis()->SetTitle("n#sigma_{K}");
   hNSigmaKecorr->GetYaxis()->SetTitle("n#sigma_{e}");

   hNSigmaPi = new TH1D("hNSigmaPi", "n_{#sigma} pions; n_{#sigma} #pi [-]; counts", 200, -10, 10 );

   hNSigmaP = new TH1D("hNSigmaP", "n_{#sigma} protons; n_{#sigma} p [-]; counts", 200, -10, 10 );

   hNSigmaK = new TH1D("hNSigmaK", "n_{#sigma} kaons; n_{#sigma} K [-]; counts", 200, -10, 10 );

   hDEdxSignal = new TH1D("hDEdxSignal", "dE/dx signal; ln#frac{dE}{dx} [MeV/cm];counts", 200, 0, 10 );

   mOutFile->cd();
   mOutFile->mkdir("RPinfo")->cd();

   //initialize 2D graphs of before and after fiducial RP cut
   hRPcorr[0] = new TH2F("hRPcorr","p_{y} vs p_{x} of protons in RP", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorr[0]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorr[0]->GetYaxis()->SetTitle("p_{y} [GeV]");
   hRPcorr[1] = new TH2F("hRPcorr_fid","p_{y} vs p_{x} of proton in RP fiducial region", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorr[1]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorr[1]->GetYaxis()->SetTitle("p_{y} [GeV]");

   hRPcorrWest[0] = new TH2F("hRPcorrWest","p_{y} vs p_{x} of protons in RP (west)", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorrWest[0]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrWest[0]->GetYaxis()->SetTitle("p_{y} [GeV]");
   hRPcorrWest[1] = new TH2F("hRPcorrWest_fid","p_{y} vs p_{x} of proton in RP fiducial region (west)", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorrWest[1]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrWest[1]->GetYaxis()->SetTitle("p_{y} [GeV]");

   hRPcorrEast[0] = new TH2F("hRPcorrEast","p_{y} vs p_{x} of protons in RP (east)", 200, -0.7, 0.8, 150, -1, 1);
   hRPcorrEast[0]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrEast[0]->GetYaxis()->SetTitle("p_{y} [GeV]");
   hRPcorrEast[1] = new TH2F("hRPcorrEast_fid","p_{y} vs p_{x} of proton in RP fiducial region (east)",  200, -0.7, 0.8, 150, -1, 1);
   hRPcorrEast[1]->GetXaxis()->SetTitle("p_{x} [GeV]");
   hRPcorrEast[1]->GetYaxis()->SetTitle("p_{y} [GeV]");


   /*
   hPIDStats[0] = new TH1D("PIDStatsTotCh", "PID stats before tot. charge cut", 4, -1.5, 2.5);
   hPIDStats[1] = new TH1D("PIDStatsExlsv", "PID stats before pT exclusive cut", 4, -1.5, 2.5);
   for (int pid = 0; pid < 3; ++pid)
      for (int stat = 0; stat < 3; ++stat)
         for (int part = 0; part < 4; ++part)
            hMSquared[pid][stat][part] = new TH1D( Form("hMSquared_%i_%i_%i",pid, stat, part), Form("hMSquared_%i_%i_%i",pid, stat, part), 200, -0.5, 1.5);
   mOutFile->cd();
    */
   
   //co su tieto grafy adc?
   /*mOutFile->mkdir("CPT2noBBCL")->cd();
   for (int iRp = 0; iRp < 2*nRomanPots; ++iRp)
   {
      hRpAdc[iRp]= new TH1D( mUtil->rpName(iRp/2) + Form("_%i_ADC",iRp%2), "ADC", 100, 0, 600);
      hRpAdcInWindow[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_ADCinTAC",iRp%2), "ADC in TAC window", 100, 0, 600);
      hRpTac[iRp]= new TH1D(mUtil->rpName(iRp/2) + Form("_%i_TAC",iRp%2), "TAC", 100, 0, 2000);
   }
   */
   mOutFile->mkdir("cuts")->cd();

   hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks",100,-2,2,100,-4,4);
   hEtaPhi->GetXaxis()->SetTitle("#eta");
   hEtaPhi->GetYaxis()->SetTitle("#varphi");
   hEtaPhi->SetTitle("Pseudorapidity and azimuthal angle distribution");

   hNumberRPTracks = new TH1D("NumberRPTracks", "Number of Tracks in RPs per event", 400, 0, 40);
   hNumberRPTracks->GetXaxis()->SetTitle("Number of tracks in RPs");
   hNumberRPTracks->GetYaxis()->SetTitle(YAxisDescription);
   hNumberRPTracks->SetTitle("Number of tracks in RPs per event");

   hGlobalTracks = new TH1D("NumberGlobalTracks", "Number of Tracks in RPs per event;n^{global}_tracks [-]", 50, 0, 50);
   hGlobalTracks->GetYaxis()->SetTitle(YAxisDescription);

   hPt = new TH1D("hPt", "Transverse momentum of hadrons", 30, 0, 3);
   hPt->SetTitle("Distribution of p_{T}");
   hPt->GetXaxis()->SetTitle("p_{T} [GeV]");
   hPt->GetYaxis()->SetTitle(YAxisDescription);

   hEta = new TH1D("hEta", "Pseudorapidity", 60, -2, 2);
   hEta->SetTitle("Distribution of pseudorapidity ");
   hEta->GetXaxis()->SetTitle("#eta");
   hEta->GetYaxis()->SetTitle(YAxisDescription);

   hDcaZ =  new TH1D("hDcaZ", "DcaZ", 100, -3,3); 
   hDcaZ->SetTitle("Distribution of DCA_{z} ");
   hDcaZ->GetXaxis()->SetTitle("DCA_{z} [cm]");
   hDcaZ->GetYaxis()->SetTitle(YAxisDescription);

   hDcaXY =  new TH1D("hDcaXY", "DcaXY", 60, -0.5,5.5); 
   hDcaXY->SetTitle("Distribution of DCA_{xy} ");
   hDcaXY->GetXaxis()->SetTitle("DCA_{xy} [cm]");
   hDcaXY->GetYaxis()->SetTitle(YAxisDescription);
 
   hNfitHits = new TH1D("hNfitHits", "NfitHits", 45, 5, 50);
   hNfitHits->SetTitle("Distribution of number of TPC hits");
   hNfitHits->GetXaxis()->SetTitle("hits");
   hNfitHits->GetYaxis()->SetTitle(YAxisDescription);
 
   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx", 45, 5, 50);
   hNhitsDEdx->SetTitle("Distribution of number of hits dE/dx ");
   hNhitsDEdx->GetXaxis()->SetTitle("Number of hits dE/dx");
   hNhitsDEdx->GetYaxis()->SetTitle(YAxisDescription);
   
   hVtxDiff =  new TH1D("hVtxDiff", "Difference in Vtx_{z} between prodVertexHypo and UPCVertex_{z}", 50, 0, 5); 
   hVtxDiff->SetTitle("Difference in determination of z_{vertex}");
   hVtxDiff->GetXaxis()->SetTitle("Z_{vertex} [cm]");
   hVtxDiff->GetYaxis()->SetTitle(YAxisDescription);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}", 60, -150, 150); 
   hPosZ->SetTitle("Distribution of position of z_{vertex}");
   hPosZ->GetXaxis()->SetTitle("Z_{vertex} [cm]");
   hPosZ->GetYaxis()->SetTitle(YAxisDescription);

   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}", 60, -150, 150); 
   hPosZCut->SetTitle("Distribution of position of z_{vertex}");
   hPosZCut->GetXaxis()->SetTitle("Z_{vertex} [cm]");
   hPosZCut->GetYaxis()->SetTitle(YAxisDescription);

   hTOFTracks = new TH1D("hTOFtracks", "Number of tracks in TOF before cut", 100, 0, 10);
   hTOFTracks->SetTitle("Number of tracks in TOF");
   hTOFTracks->GetXaxis()->SetTitle("Number of tracks");
   hTOFTracks->GetYaxis()->SetTitle(YAxisDescription);
 
   hNVertices = new TH1D("hNVertices", "Number of vertices before cut", 100, 0, 10);
   hNVertices->SetTitle("Number of vertices");
   hNVertices->GetXaxis()->SetTitle("Number of vertices");
   hNVertices->GetYaxis()->SetTitle(YAxisDescription);

   hNPairV0 = new TH1D("hNPairV0", "Number of pairs after cuts", 30, 0, 30);
   hNPairV0->SetTitle("Number of pairs");
   hNPairV0->GetXaxis()->SetTitle("Number of pairs");
   hNPairV0->GetYaxis()->SetTitle(YAxisDescription);

   hTotQ = new TH1D("TotQ", "Total charge of pair", 1000, -10, 10);
   hTotQ->SetTitle("Total charge of pair");
   hTotQ->GetXaxis()->SetTitle("Charge of pair");
   hTotQ->GetYaxis()->SetTitle(YAxisDescription);

   hSameTrackPair = new TH1D("hSameTrackPair", "pair with one track", 2, 0, 2);
   hSameTrackPair->SetTitle("pair with same track");
   hSameTrackPair->GetXaxis()->SetTitle("sameTrackPair [-]");
   hSameTrackPair->GetYaxis()->SetTitle(YAxisDescription);

   hBothFlags = new TH1D("hBothFlags", "track with kPrimary and kV0", 2, 0, 2);
   hBothFlags->SetTitle("track with kPrimary and kV0");
   hBothFlags->GetXaxis()->SetTitle("bothFlags [-]");
   hBothFlags->GetYaxis()->SetTitle(YAxisDescription);

   mOutFile->mkdir("topologyCuts")->cd();

   hDecayLength = new TH1D("hDecayLength", "decaylength", 100, 0., 10.);
   hDecayLength->SetTitle("DecayLength");
   hDecayLength->GetXaxis()->SetTitle("decayLength [cm]");
   hDecayLength->GetYaxis()->SetTitle(YAxisDescription);

   hDecayLengthCut = new TH1D("hDecayLengthCut", "decaylength", 100, 0., 10);
   hDecayLengthCut->SetTitle("DecayLength");
   hDecayLengthCut->GetXaxis()->SetTitle("decayLength [cm]");
   hDecayLengthCut->GetYaxis()->SetTitle(YAxisDescription);

   hPointingAngle = new TH1D("hPointingAngle", "pointing angle", 100, -1, 1);
   hPointingAngle->SetTitle("Pointing angle");
   hPointingAngle->GetXaxis()->SetTitle("cos(pointingAngle) [-]");
   hPointingAngle->GetYaxis()->SetTitle(YAxisDescription);

   hPointingAngleCut = new TH1D("hPointingAngleCut", "pointing angle", 100, -1, 1);
   hPointingAngleCut->SetTitle("Pointing Angle");
   hPointingAngleCut->GetXaxis()->SetTitle("cos(pointingangle) [-]");
   hPointingAngleCut->GetYaxis()->SetTitle(YAxisDescription);

   hDcaDaughters = new TH1D("hDcaDaughters", "dca between pairs", 100, 0., 10);
   hDcaDaughters->SetTitle("DCA between pair");
   hDcaDaughters->GetXaxis()->SetTitle("DCA [cm]");
   hDcaDaughters->GetYaxis()->SetTitle(YAxisDescription);

   hDcaDaughtersCut = new TH1D("hDcaDaughtersCut", "dca between pairs", 100, 0., 10);
   hDcaDaughtersCut->SetTitle("DCA between pair");
   hDcaDaughtersCut->GetXaxis()->SetTitle("DCA [cm]");
   hDcaDaughtersCut->GetYaxis()->SetTitle(YAxisDescription);

   hDcaBeamline = new TH1D("hDcaBeamline", "dca between pair vertex and beamline", 100, 0., 10);
   hDcaBeamline->SetTitle("DCA between pair vertex and beamline");
   hDcaBeamline->GetXaxis()->SetTitle("DCA [cm]");
   hDcaBeamline->GetYaxis()->SetTitle(YAxisDescription);

   hDcaBeamlineCut = new TH1D("hDcaBeamlineCut", "dca between pairs vertex and beamline", 100, 0., 10);
   hDcaBeamlineCut->SetTitle("DCA between pair");
   hDcaBeamlineCut->GetXaxis()->SetTitle("DCA [cm]");
   hDcaBeamlineCut->GetYaxis()->SetTitle(YAxisDescription);

   hDecayLPointingA = new TH2F("hDecayLPointingA", "Decay length and pointing angle corr plot;  Decay length [cm]; Pointing angle [-]", 100,0.,10.,100,-1.5, 1.5);
   hDecayLPointingA->SetTitle("Correlation plot of decay length and pointing angle");

   hDecayLPointingACut = new TH2F("hDecayLPointingACut", "Decay length and pointing angle corr plot;  Decay length [cm]; Pointing angle [-]", 100,0.,10.,100,-1.5, 1.5);
   hDecayLPointingACut->SetTitle("Correlation plot of decay length and pointing angle");

   hArmenterosPodolanski = new TH2F("hArmenterosPodolanski", "Armenteros- Podolanski plot; #frac{(p_{L}^{+} -p_{L}^{-})}{(p_{L}^{+} + p_{L}^{-})}; p_{T} [GeV/c]", 200, -3.,3., 200, 0.,2. );
   hArmenterosPodolanski->SetTitle("Armenteros Podolanski plot");

   hInvMassEta = new TH2F("hInvMassEta", "hInvMassEta; m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; #eta [-]", 100, 0.,1., 200, -1.,1. );
   hInvMassEta->SetTitle("correlation plot of invMass and eta");

    hEtaTag = new TH1D("hEtaTag", "hEtaTag; #eta_{tag} [-]; counts", 20, -1,1);
    hEtaProbe = new TH1D("hEtaProbe", "hEtaProbe; #eta_{probe} [-]; counts", 20, -1,1);

    hpTTag = new TH1D("hpTTag", "hpTTag; p_{T}^{tag} [GeV/c]; counts", 15, 0,1.5);
    hpTProbe = new TH1D("hpTProbe", "hpTProbe; p_{T}^{probe} [GeV/c]; counts", 15, 0,1.5);

    hPhiTag = new TH1D("hPhiTag", "hPhiTag; #phi_{tag} [GeV/c]; counts", 12, -3.15, 3.15);
    hPhiProbe = new TH1D("hPhiProbe", "hPhiProbe; #phi_{probe} [GeV/c]; counts", 12, -3.15, 3.15);

    hEtaPhiTag = new TH2F("hEtaPhiTag", "hEtaPhiTag; #eta_{tag} [-]; #phi_{tag} [-]", 20,-1,1,12,-3.15, 3.15);
    hEtaPhiProbe = new TH2F("hEtaPhiProbe", "hEtaPhiProbe; #eta_{probe} [-]; #phi_{probe} [-]", 20,-1,1,12,-3.15, 3.15);


}


void AnaV0::fillTag(const StUPCTrack* trk) {

    hEtaTag->Fill( trk->getEta() );
    hPhiTag->Fill( trk->getPhi() );
    hpTTag->Fill( trk->getPt() );
    hEtaPhiTag->Fill(trk->getEta() , trk->getPhi() );

}

void AnaV0::fillProbe(const StUPCTrack* trk) {
    
    hEtaProbe->Fill( trk->getEta() );
    hPhiProbe->Fill( trk->getPhi() );
    hpTProbe->Fill( trk->getPt() );
    hEtaPhiProbe->Fill(trk->getEta() , trk->getPhi() );

}


void AnaV0::fillFinalPlots() {

    for (int iPair = 0; iPair < nStates; ++iPair){
        if(mRecTree->getPairID(iPair) == 0){
            hInvMassEta->Fill(mRecTree->getInvMass(iPair), mRecTree->getEta(iPair));
        }

        if( !(mRecTree->getPairID(iPair) == 0 && mRecTree->getInvMass(iPair) > 0.45 && mRecTree->getInvMass(iPair) < 0.55) && !( (mRecTree->getPairID(iPair) == 1 || mRecTree->getPairID(iPair) ==2) && mRecTree->getInvMass(iPair) > 1.1 && mRecTree->getInvMass(iPair) <1.2 ))
            continue;
        Double_t pLPlus, pLMinus, alpha;
        if (mRecTree->getCharge(2*iPair) == 1) {
            pLPlus = mRecTree->getPzInGev(2*iPair);
            pLMinus = mRecTree->getPzInGev(2*iPair + 1);
        } else{
            pLMinus = mRecTree->getPzInGev(2*iPair);
            pLPlus = mRecTree->getPzInGev(2*iPair + 1);            
        }

        alpha = (pLPlus - pLMinus)/(pLPlus + pLMinus);

        hArmenterosPodolanski->Fill(alpha, mRecTree->getPt(iPair));

    }//forloop
}

void AnaV0::fillTopologyCutsBefore(const StUPCV0& V0){

    hDcaDaughters->Fill(V0.dcaDaughters());
    hDcaBeamline->Fill(V0.DCABeamLine());
    hPointingAngle->Fill(V0.pointingAngleHypo());
    hDecayLength->Fill(V0.decayLengthHypo());
    hDecayLPointingA->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());

}


void AnaV0::fillTopologyCutsAfter(const StUPCV0& V0){

    hDcaDaughtersCut->Fill(V0.dcaDaughters());
    hDcaBeamlineCut->Fill(V0.DCABeamLine());
    hPointingAngleCut->Fill(V0.pointingAngleHypo());
    hDecayLengthCut->Fill(V0.decayLengthHypo());
    hDecayLPointingACut->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());
}

bool AnaV0::shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b) {
    TVector3 vtx0 = {0,0,0};
    double bField = mUpcEvt->getMagneticField();
    double beamline[4]; 
    beamline[0] = mUpcEvt->getBeamXPosition();
    beamline[2] = mUpcEvt->getBeamXSlope();
    beamline[1] = mUpcEvt->getBeamYPosition();
    beamline[3] = mUpcEvt->getBeamYSlope();

    const StUPCTrack* currentTrk1 = mUpcEvt->getTrack(a.first);
    const StUPCTrack* currentTrk2 = mUpcEvt->getTrack(a.second);
    
    const StUPCTrack* maxTrk1 = mUpcEvt->getTrack(b.first);
    const StUPCTrack* maxTrk2 = mUpcEvt->getTrack(b.second);


    StUPCV0 currentK0(currentTrk1,currentTrk2, mUtil->mass(hasGoodTPCnSigma(currentTrk1)), mUtil->mass(hasGoodTPCnSigma(currentTrk2)),a.first,a.second, vtx0, beamline, bField, false);
    StUPCV0 maxK0(maxTrk1,maxTrk2, mUtil->mass(hasGoodTPCnSigma(maxTrk1)), mUtil->mass(hasGoodTPCnSigma(maxTrk2)),b.first,b.second, vtx0, beamline, bField, false);


    if (currentK0.dcaDaughters() < maxK0.dcaDaughters()){
        return true;
    }else return false;

}

vector<pair<int, int>> AnaV0::filterPairs(vector<pair<int, int>>& pairs) {
    set<pair<int, int>> result; // Use a set to automatically handle duplicates

    // Normalize pairs: Ensure first element is always less than the second
    for (auto& p : pairs) {
        if (p.first > p.second) {
            swap(p.first, p.second);
        }
    }

    // Sort pairs to bring related pairs closer together
    sort(pairs.begin(), pairs.end());

    for (size_t i = 0; i < pairs.size(); ++i) {
        bool keep = true;
        // Compare with the next pair if it exists
        if (i + 1 < pairs.size()) {
            // Check if the pairs share a value
            if (pairs[i].first == pairs[i + 1].first || pairs[i].second == pairs[i + 1].second ||
                pairs[i].first == pairs[i + 1].second || pairs[i].second == pairs[i + 1].first) {
                // Apply the boolean condition to decide which pair to keep
                if (!shouldKeepPair(pairs[i], pairs[i + 1])) {
                    keep = false; // Don't keep the current pair
                }
            }
        }
        if (keep) {
            result.insert(pairs[i]);
        }
    }

    return vector<pair<int, int>>(result.begin(), result.end());
}

bool AnaV0::lambda(const StUPCTrack *trk1,const StUPCTrack *trk2){
    if (hasGoodTPCnSigma(trk1) == PROTON && trk1->getCharge() > 0)
        return true;
    else if(hasGoodTPCnSigma(trk2) == PROTON && trk2->getCharge() > 0)
        return true;
    else 
        return false;
}



bool AnaV0::oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2){
    Double_t charge;
    charge = trk1->getCharge() + trk2->getCharge();
    if(charge == 0)
        return true;
    else
        return false;
}

bool AnaV0::is2pions(const StUPCTrack *trk1,const StUPCTrack *trk2){

    if((hasGoodTPCnSigma(trk1) == PION) && (hasGoodTPCnSigma(trk2) == PION)){
        return true;
    }
    else
        return false;

}
bool AnaV0::isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2){

   if( hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PROTON ){
       return true;
   }else if(hasGoodTPCnSigma(trk1) == PROTON && hasGoodTPCnSigma(trk2) == PION){
       return true;     
   }
   else
       return false;
}

int AnaV0::hasGoodTPCnSigma(const StUPCTrack *trk){ //zmiernit podmienku na proton, skusit prebehnut aj na K0, Lambda 
  //check for good nSigma of hadron in TPC, 10 means unidentified
    if((trk->getNSigmasTPCProton() < 3) &&  (trk->getNSigmasTPCKaon() > 3) && (trk->getNSigmasTPCPion() > 3))
        return PROTON;
    else if((trk->getNSigmasTPCProton() > 3) && (trk->getNSigmasTPCKaon() < 3) && (trk->getNSigmasTPCPion()>3))
        return KAON;
    else if(trk->getNSigmasTPCPion() < 3)
        return PION;
    else
        return 10;
}

void AnaV0::fillGoodTrackCuts(const StUPCTrack* trk){
    hDcaZ->Fill(trk->getDcaZ());
    hDcaXY->Fill(trk->getDcaXY());
    hNfitHits->Fill(trk->getNhitsFit());
    hNhitsDEdx->Fill(trk->getNhitsDEdx());
    hPt->Fill(trk->getPt());
}


void AnaV0::saveNSigmaCorr(const StUPCTrack *trk){
    Double_t nSigmaPion, nSigmaProton, nSigmaKaon, nSigmaElectron;
    nSigmaPion = trk->getNSigmasTPCPion();
    nSigmaProton = trk->getNSigmasTPCProton();
    nSigmaKaon = trk->getNSigmasTPCKaon();
    nSigmaElectron = trk->getNSigmasTPCElectron();
    hNSigmaPi->Fill(nSigmaPion);
    hNSigmaP->Fill(nSigmaProton);
    hNSigmaK->Fill(nSigmaKaon);
    hDEdxSignal->Fill( trk->getDEdxSignal() );
    hNSigmaPiPcorr->Fill(nSigmaPion, nSigmaProton);
    hNSigmaPiKcorr->Fill(nSigmaPion, nSigmaKaon);
    hNSigmaPiecorr->Fill(nSigmaPion, nSigmaElectron);
    hNSigmaPKcorr->Fill(nSigmaProton, nSigmaKaon);
    hNSigmaPecorr->Fill(nSigmaProton, nSigmaElectron);
    hNSigmaPPicorr->Fill(nSigmaProton, nSigmaPion);
    hNSigmaKecorr->Fill(nSigmaKaon, nSigmaElectron);
    hNSigmaKPcorr->Fill(nSigmaKaon, nSigmaProton);
    hNSigmaKPicorr->Fill(nSigmaKaon, nSigmaPion);
}

/*
void AnaV0::CalculateTOFEff(unsigned int tagID)
{
   // eta cut and PID not implemented
   unsigned int vertexID = mUpcEvt->getVertex(0)->getId();
   vector<int> hadronId;

   TLorentzVector tag, probe, state;
   const StUPCTrack* tagTrk = mUpcEvt->getTrack(tagID);
   int tagCharge = static_cast<int>( tagTrk->getCharge() );
   tagTrk->getLorentzVector(tag, mUtil->mass(PION));
   if( abs(tagTrk->getNSigmasTPC(StUPCTrack::kPion)) > 3)
      return;
   if( abs(tagTrk->getEta()) > 0.7 )
      return;

   for(unsigned int trackID = 0; trackID < mUpcEvt->getNumberOfTracks(); ++trackID)
   {
      if( trackID == tagID)
         continue;

      const StUPCTrack* trk = mUpcEvt->getTrack(trackID);
      if( !trk->getFlag(StUPCTrack::kPrimary) || trk->getVertexId() != vertexID || !IsGoodTrack(trk)) 
         continue;

      if( abs(trk->getNSigmasTPC(StUPCTrack::kPion)) > 3)
         continue;

      if( static_cast<int>( trk->getCharge() ) == -tagCharge  );
         hadronId.push_back(trackID);
   } 

   double missPt, minMissPt;
   unsigned int probeID;
   minMissPt = 999;
   nProbes = hadronId.size();
   for (unsigned int id = 0; id < hadronId.size(); ++id)
   {
      const StUPCTrack* trk = mUpcEvt->getTrack(hadronId[id]);

      trk->getLorentzVector(probe, mUtil->mass(PION));
      state = tag + probe;

      missPt = (state.Vect() + mRPpTBalance).Pt();
      if( missPt < minMissPt){
         minMissPt = missPt;
         probeID = hadronId[id];
      }
   }
   
   if( minMissPt > 2 )
      return;

   const StUPCTrack* probeTrk = mUpcEvt->getTrack(probeID);
   probeTrk->getLorentzVector(probe, mUtil->mass(PION));
   state = tag + probe;

   deltaPhiProtons = abs(mRecTree->getPhiRp(East) - mRecTree->getPhiRp(West)); 
   pTState = state.Pt();
   phiState = abs(state.Phi());
   double deltaPhiPions = abs(tag.Phi() - probe.Phi());
   double piHalf = mUtil->pi() / 2;

   
   if( deltaPhiProtons > piHalf){ // is elastic
      if(pTState > 0.7)
         return;
      if( deltaPhiPions < 2.36)
         return;
   }else{ // is inelastic
      if( pTState < 0.7)
         return;
      if( phiState < 1.14 || phiState > 2.0)
         return;
      if( deltaPhiPions > 2.36)
         return;
   }

   // Fill sample total
   pTMiss = minMissPt;
   invMass = state.M();

   tagPt = tagTrk->getPt();
   tagEta = tagTrk->getEta();
   tagPhi = tagTrk->getPhi();

   probePt = probeTrk->getPt();
   probeEta = probeTrk->getEta();
   probePhi = probeTrk->getPhi();
   probeTof = probeTrk->getFlag(StUPCTrack::kTof);
   probeTofHit = IsGoodTofTrack(probeTrk);

   tofEffTree->Fill();
}

void AnaV0::CalculatePID()
{

   Int_t pairID = -1;
   double chiPair[nParticles];
   for(int iPart = 0; iPart < nParticles; ++iPart)
      chiPair[iPart] = mRecTree->getNSigmaTPC(0, iPart)*mRecTree->getNSigmaTPC(0, iPart) + 
                     mRecTree->getNSigmaTPC(1, iPart)*mRecTree->getNSigmaTPC(1, iPart);

   double deltaTOF = mRecTree->getTofTimeInNs(1) - mRecTree->getTofTimeInNs(0);
   double speedOfLight2 = mUtil->c()*mUtil->c();
   double speedOfLight4 = speedOfLight2*speedOfLight2; 
   double length1Squared = mRecTree->getTofLengthInCm(0)*mRecTree->getTofLengthInCm(0)/(100*100); // convert TOFlength from cm to m
   double length2Squared = mRecTree->getTofLengthInCm(1)*mRecTree->getTofLengthInCm(1)/(100*100); // convert TOFlength from cm to m
   double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
   double deltaTime4 = deltaTime2*deltaTime2;
   double oneOverMomentum1sq = 1/(mRecTree->getMomentumInGeV(0)*mRecTree->getMomentumInGeV(0));
   double oneOverMomentum2sq = 1/(mRecTree->getMomentumInGeV(1)*mRecTree->getMomentumInGeV(1));
   double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
   double bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
   double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
   double mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);
   double pT[] = { mRecTree->getPtInGev(0), mRecTree->getPtInGev(1)};

   if(chiPair[PION] > 9 && chiPair[KAON] > 9 && chiPair[PROTON] < 9 && mSquared > 0.6) // it is... proton!
   {
      if(pT[0] > 0.4 && pT[1] > 0.4 && (pT[0] < 1.1 || pT[1] < 1.1) )
         pairID = PROTON;
   }
   else if(chiPair[PION] > 9 && chiPair[KAON] < 9 && chiPair[PROTON] > 9 && mSquared > 0.15) // it is... kaon!
   {
      if(pT[0] > 0.3 && pT[1] > 0.3 && (pT[0] < 0.7 || pT[1] < 0.7) )
         pairID = KAON;
   }
   else if( chiPair[PION] < 12) // it is... pion!
   {
      pairID = PION;
   }
   
   mRecTree->setPairID(pairID);
   mRecTree->setMSquared(mSquared);
}

Int_t AnaV0::PIDStartegy(int strategy)
{
   // 0 - the main PID stragy without mSquared
   // 1 - the main PID stragy without mSquared and without pT cuts
   // else return -1
   if( strategy < 0 || strategy > 1)
      return -1;

   Int_t pairID = -1;
   double chiPair[nParticles];
   for(int iPart = 0; iPart < nParticles; ++iPart)
      chiPair[iPart] = mRecTree->getNSigmaTPC(0, iPart)*mRecTree->getNSigmaTPC(0, iPart) + 
                     mRecTree->getNSigmaTPC(1, iPart)*mRecTree->getNSigmaTPC(1, iPart);
   double pT[] = { mRecTree->getPtInGev(0), mRecTree->getPtInGev(1)};

   if(chiPair[PION] > 9 && chiPair[KAON] > 9 && chiPair[PROTON] < 9) // it is... proton!
   {
      if(pT[0] > 0.4 && pT[1] > 0.4 && (pT[0] < 1.1 || pT[1] < 1.1) && strategy == 0)
         pairID = PROTON;
      if( strategy == 1)
         pairID = PROTON;
   }
   else if(chiPair[PION] > 9 && chiPair[KAON] < 9 && chiPair[PROTON] > 9) // it is... kaon!
   {
      if(pT[0] > 0.3 && pT[1] > 0.3 && (pT[0] < 0.7 || pT[1] < 0.7) && strategy == 0 )
         pairID = KAON;
      if( strategy == 1)
         pairID = KAON;
   }
   else if( chiPair[PION] < 12) // it is... pion!
   {
      pairID = PION;
   }
   
   return pairID;
}

bool AnaV0::IsPairOf(int type)
{
   // 0 - PION
   // 1 - KAON
   // 2 - PROTON
   if( type < 0 || type > 2)
      return false;

   return (mRecTree->getNSigmaTPC(0, type) < 2 && mRecTree->getNSigmaTPC(1, type) < 2);

}
*/
void AnaV0::SaveMissingMomenta(TVector3 missP)
{
   mRecTree->setPtMissing( missP.Pt() );
   mRecTree->setPxMissing( missP.X() );
   mRecTree->setPyMissing( missP.Y() );   
}

/*
void AnaV0::FillMSquared()
{
   int totCharge = mRecTree->getCharge(0) + mRecTree->getCharge(1);
   double pTMissing = mRecTree->getPtMissing();
   double mSquared = mRecTree->getMSquared();
   
   // hMSquared[strategy][state][particle]
   // particle: All Pion Kaon Proton
   // state: before tot charge | after tot charge | exclusive
   // strategy: standard PID | standard PID w/O pT | selective dEdx 

   // fill hMSquared for "before tot charge cut"
   for (int iState = 0; iState < 3; ++iState)
   {
      // iState == 0 -> fill hMSquared for "before tot charge cut"
      // iState == 1 -> fill hMSquared for "after tot charge cut"
      // iState == 2 -> fill hMSquared for "after exclusive cut"

      if( iState == 1 && totCharge != 0 )
         continue;

      if( iState == 2 && pTMissing > exclusivityCut )
         continue;

      for (int iStrategy = 0; iStrategy < 3; ++iStrategy) // for each strategy
         hMSquared[iStrategy][iState][0]->Fill(mSquared);

      int partTypeStrategyOne = PIDStartegy(0) + 1;
      int partTypeStrategyTwo = PIDStartegy(1) + 1;

      if( partTypeStrategyOne > 0 )
         hMSquared[0][iState][partTypeStrategyOne]->Fill(mSquared);

      if( partTypeStrategyTwo > 0 )
         hMSquared[1][iState][partTypeStrategyTwo]->Fill(mSquared);

      for(int iPart = 0; iPart < nParticles; ++iPart)
         if( IsPairOf(iPart))
            hMSquared[2][iState][iPart+1]->Fill(mSquared);
   }
}
*/