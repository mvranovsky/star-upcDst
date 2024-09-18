#include "AnaV0Control.h"

AnaV0Control::AnaV0Control(TFile *outFile): Ana(outFile){}

void AnaV0Control::Make()
{
   hAnalysisFlow->Fill(V0ALL);
   if(!CheckTriggers(&triggerID, mUpcEvt, hTriggerBits))
      return;
   hAnalysisFlow->Fill(V0TRIG);

   vector<int> tracksPrimary;
   for (int itrk = 0; itrk < mUpcEvt->getNumberOfTracks(); ++itrk){
       const StUPCTrack* trk = mUpcEvt->getTrack(itrk);
       
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
       hEta->Fill(eta);
       hEtaPhi->Fill(eta, phi);
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

       //mozno zmenit pre tracksV0 : else if -> if?
       if( trk->getFlag( StUPCTrack::kPrimary ) )
           tracksPrimary.push_back(itrk);
       hAnalysisFlow->Fill(V0FLAG);
   }//global tracks loop

   mRecTree->setNGoodTpcTrks( tracksPrimary.size() );
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
   TVector3 vertex, vertex0;
   vertex0 = {0,0,0};
   int HADRON1, HADRON2, totalCharge;
   int vertexIdTrk1, vertexIdTrk2;
   Bool_t tofMatchTrk1, tofMatchTrk2;
   Double_t vertexDiff;
   int pairId = 0;
   vector<pair<int, int>> hadronPairPrimary;

   
   //loop over all tracks with primary flag
   for (unsigned int itrk = 0; itrk < tracksPrimary.size(); itrk++){
       const StUPCTrack* trk1 = mUpcEvt->getTrack(tracksPrimary[itrk]);
       tofMatchTrk1 = trk1->getFlag( StUPCTrack::kTof );
       HADRON1 = hasGoodTPCnSigma(trk1);
       vertexIdTrk1 = trk1->getVertexId();

       for (unsigned int jtrk = itrk; jtrk < tracksPrimary.size(); ++jtrk){
           if (itrk == jtrk)
               continue; 
           const StUPCTrack* trk2 = mUpcEvt->getTrack(tracksPrimary[jtrk]);
           tofMatchTrk2 = trk2->getFlag( StUPCTrack::kTof );
           HADRON2 = hasGoodTPCnSigma(trk2);
           vertexIdTrk2 = trk2->getVertexId();

           //check for at least one tof match
           if (!tofMatchTrk1 && !tofMatchTrk2)
               continue;
           //checking for wrong identification
           if(!(HADRON1 == PION || HADRON1 == PROTON))
               continue;
           if(!(HADRON2 == PION || HADRON2 == PROTON))
               continue;
           if (HADRON1 == PROTON && HADRON2 == PROTON)
               continue; 

        

           StUPCV0 K0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), tracksPrimary[itrk],tracksPrimary[jtrk], vertex0, beamline, bField, false);
           vertex = K0.prodVertexHypo();
           

           //same cuts as were in StUPCSelectV0 except K0L1, K0L2 -> K0
           fillTopologyCutsBefore(K0);
           if ( !(K0.dcaDaughters() < maxDcaDaughters && K0.DCABeamLine() < maxDcaBeamLine && (K0.pointingAngleHypo()> minPointingAngle || K0.decayLengthHypo()< maxDecayLengthHypo) ) )
              continue;
           fillTopologyCutsAfter(K0);
           hAnalysisFlow->Fill(V0PAIR);

           
           //save the difference between vertex positions
           const StUPCVertex* vtx = mUpcEvt->getVertex(vertexIdTrk1);
           if (vtx && vertexIdTrk1 == vertexIdTrk2){
               vertexDiff = vtx->getPosZ() - vertex.Z();
               hVtxDiff->Fill(vertexDiff);
           }



           //cut for hypothetical vertex z position to be < 80 cm from the middle of the detector
           //used vertex from UPCVertex possible change?
           hPosZ->Fill(vertex.Z());
           if (abs(vertex.Z()) > vertexRange )
               continue;
           hPosZCut->Fill(vertex.Z()); 
           hAnalysisFlow->Fill(V0ZVERTEX);
           hadronPairPrimary.push_back(make_pair(tracksPrimary[itrk],tracksPrimary[jtrk]));
           pairId++;
       }//end of inner tracks loop
   }//end of outer tracks loop

   //checking for pair duplicates and not empty vector
   if(hadronPairPrimary.size() != 0)
       hadronPairPrimary = filterPairs(hadronPairPrimary);
   
   //tree can only hold max 5 states per event
   hNPairs->Fill(hadronPairPrimary.size());
   if (hadronPairPrimary.size() > 5){
       //cerr<< "too many pairs. Number of pairs: " << hadronPairPrimary.size() << endl; 
       hadronPairPrimary.resize(5);
   }
   for (int i = 0; i < hadronPairPrimary.size(); ++i){
        //checking for pair made of one track
       if (hadronPairPrimary[i].first == hadronPairPrimary[i].second){
          hSameTrackPair->Fill(1);
          cout << "We got a pair with the same track for both. " << endl;     
       }

      const StUPCTrack* trk1 = mUpcEvt->getTrack(hadronPairPrimary[i].first);
      const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronPairPrimary[i].second);

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
         mRecTree->setPairID(K0S,i);
         if(totalCharge == 0)
           hAnalysisFlow->Fill(V0PIPI);
      } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && lambda(trk1, trk2) ){ 
         mRecTree->setPairID(LAMBDA,i);
         if(totalCharge == 0)
           hAnalysisFlow->Fill(V0PPI);
      } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && !lambda(trk1, trk2) ){
         mRecTree->setPairID(LAMBDABAR,i);
         if(totalCharge == 0)
           hAnalysisFlow->Fill(V0PIPBAR); 
      } else continue;

      StUPCV0 V0(trk1, trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronPairPrimary[i].first,hadronPairPrimary[i].second, vertex0, beamline, bField, false);  
      vertex = V0.prodVertexHypo();
      //saving info for tracks and vertices
      SaveVertexInfo(vertex, i);
      SaveTrackInfo(trk1,i*2);
      SaveTrackInfo(trk2,i*2 +1);
      SaveStateInfo(state,totalCharge, i);

   }
   //cerr << "Got through the whole thing" << endl;

   // save event info
   mRecTree->FillRecTree();

}//end of make()

void AnaV0Control::Init()
{
    //initialize class for selecting V0 candidates 
    //mSelectV0 = new StUPCSelectV0Modified();

    //init database to read beam-line parameters
    //mDbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    //mDbMk->Init();

   mUtil = new Util();


   if( DEBUG )
      cout<<"AnaV0Control::Init() called"<<endl;

   mOutFile->cd();
   hAnalysisFlow = new TH1D("AnalysisFlow", "CutsFlow", nV0SelectionCuts-1, 1, nV0SelectionCuts);
   for(int tb=1; tb<nV0SelectionCuts; ++tb) {
      hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisV0SelectionName(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
      TString label; label.Form("%d",triggerID[tb]);
      hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   mRecTree = new RecTree(nameOfAnaV0ControlTree, AnaV0ControlTreeBits, false); 

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
 
   hNfitHits = new TH1D("hNfitHits", "NfitHits", 500, -0.5, 49.5);
   hNfitHits->SetTitle("Distribution of number of TPC hits");
   hNfitHits->GetXaxis()->SetTitle("hits");
   hNfitHits->GetYaxis()->SetTitle(YAxisDescription);
 
   hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx", 500, -0.5, 49.5);
   hNhitsDEdx->SetTitle("Distribution of number of hits dE/dx ");
   hNhitsDEdx->GetXaxis()->SetTitle("Number of hits dE/dx");
   hNhitsDEdx->GetYaxis()->SetTitle(YAxisDescription);
   
   hVtxDiff =  new TH1D("hVtxDiff", "Difference in Vtx_{z} between prodVertexHypo and UPCVertex_{z}", 40, -20, 20); 
   hVtxDiff->SetTitle("Difference in determination of z_{vertex}");
   hVtxDiff->GetXaxis()->SetTitle("Z_{vertex} [cm]");
   hVtxDiff->GetYaxis()->SetTitle(YAxisDescription);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}", 80, -200, 200); 
   hPosZ->SetTitle("Distribution of position of z_{vertex}");
   hPosZ->GetXaxis()->SetTitle("Z_{vertex} [cm]");
   hPosZ->GetYaxis()->SetTitle(YAxisDescription);

   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}", 80, -200, 200); 
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

   hNPairs = new TH1D("hNPairs", "Number of pairs after cuts", 30, 0, 30);
   hNPairs->SetTitle("Number of pairs");
   hNPairs->GetXaxis()->SetTitle("Number of pairs");
   hNPairs->GetYaxis()->SetTitle(YAxisDescription);


   hTotQ = new TH1D("TotQ", "Total charge of pair", 1000, -10, 10);
   hTotQ->SetTitle("Total charge of pair");
   hTotQ->GetXaxis()->SetTitle("Charge of pair");
   hTotQ->GetYaxis()->SetTitle(YAxisDescription);

   hSameTrackPair = new TH1D("hSameTrackPair", "pair with one track", 2, 0, 2);
   hSameTrackPair->SetTitle("pair with same track");
   hSameTrackPair->GetXaxis()->SetTitle("sameTrackPair [-]");
   hSameTrackPair->GetYaxis()->SetTitle(YAxisDescription);

   hBothFlags = new TH1D("hBothFlags", "pair with V0 and Primary flag", 2, 0, 2);
   hBothFlags->SetTitle("pair 2 flags");
   hBothFlags->GetXaxis()->SetTitle("2 flags [-]");
   hBothFlags->GetYaxis()->SetTitle(YAxisDescription);

   mOutFile->mkdir("topologyCuts")->cd();

   hDecayLength = new TH1D("hDecayLength", "decaylength", 200, -10, 10);
   hDecayLength->SetTitle("DecayLength");
   hDecayLength->GetXaxis()->SetTitle("decayLength [cm]");
   hDecayLength->GetYaxis()->SetTitle(YAxisDescription);

   hDecayLengthCut = new TH1D("hDecayLengthCut", "decaylength", 200, -10, 10);
   hDecayLengthCut->SetTitle("DecayLength");
   hDecayLengthCut->GetXaxis()->SetTitle("decayLength [cm]");
   hDecayLengthCut->GetYaxis()->SetTitle(YAxisDescription);

   hPointingAngle = new TH1D("hPointingAngle", "pointing angle", 100, -2, 2);
   hPointingAngle->SetTitle("Pointing angle");
   hPointingAngle->GetXaxis()->SetTitle("cos(pointingAngle) [-]");
   hPointingAngle->GetYaxis()->SetTitle(YAxisDescription);

   hPointingAngleCut = new TH1D("hPointingAngleCut", "pointing angle", 100, -2, 2);
   hPointingAngleCut->SetTitle("Pointing Angle");
   hPointingAngleCut->GetXaxis()->SetTitle("cos(pointingangle) [-]");
   hPointingAngleCut->GetYaxis()->SetTitle(YAxisDescription);

   hDcaDaughters = new TH1D("hDcaDaughters", "dca between pairs", 200, -10, 10);
   hDcaDaughters->SetTitle("DCA between pair");
   hDcaDaughters->GetXaxis()->SetTitle("DCA [cm]");
   hDcaDaughters->GetYaxis()->SetTitle(YAxisDescription);

   hDcaDaughtersCut = new TH1D("hDcaDaughtersCut", "dca between pairs", 200, -10, 10);
   hDcaDaughtersCut->SetTitle("DCA between pair");
   hDcaDaughtersCut->GetXaxis()->SetTitle("DCA [cm]");
   hDcaDaughtersCut->GetYaxis()->SetTitle(YAxisDescription);

   hDcaBeamline = new TH1D("hDcaBeamline", "dca between pair vertex and beamline", 200, -10, 10);
   hDcaBeamline->SetTitle("DCA between pair vertex and beamline");
   hDcaBeamline->GetXaxis()->SetTitle("DCA [cm]");
   hDcaBeamline->GetYaxis()->SetTitle(YAxisDescription);

   hDcaBeamlineCut = new TH1D("hDcaBeamlineCut", "dca between pairs vertex and beamline", 200, -10, 10);
   hDcaBeamlineCut->SetTitle("DCA between pair");
   hDcaBeamlineCut->GetXaxis()->SetTitle("DCA [cm]");
   hDcaBeamlineCut->GetYaxis()->SetTitle(YAxisDescription);


}

void AnaV0Control::fillTopologyCutsBefore(const StUPCV0& V0){

    hDcaDaughters->Fill(V0.dcaDaughters());
    hDcaBeamline->Fill(V0.DCABeamLine());
    hPointingAngle->Fill(V0.pointingAngleHypo());
    hDecayLength->Fill(V0.decayLengthHypo());

}


void AnaV0Control::fillTopologyCutsAfter(const StUPCV0& V0){

    hDcaDaughtersCut->Fill(V0.dcaDaughters());
    hDcaBeamlineCut->Fill(V0.DCABeamLine());
    hPointingAngleCut->Fill(V0.pointingAngleHypo());
    hDecayLengthCut->Fill(V0.decayLengthHypo());

}

bool AnaV0Control::shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b) {
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

vector<pair<int, int>> AnaV0Control::filterPairs(vector<pair<int, int>>& pairs) {
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
/*

int main() {
    vector<pair<int, int>> pairs = {{1, 3}, {2, 4}, {3, 4}, {4, 2}, {5, 6}, {6, 5}};
    auto filteredPairs = filterPairs(pairs);

    for (auto& pair : filteredPairs) {
        cout << "{" << pair.first << ", " << pair.second << "} ";
    }

    return 0;
}




void AnaV0::compareAndStorePair(map<int, int>& pairMap, const pair<int, int>& pair) {
    auto it = pairMap.find(pair.first);
    if (it != pairMap.end()) {
        if (it->second < pair.second) {
            it->second = pair.second;
        }
    } else {
        pairMap[pair.first] = pair.second;
    }
}

void AnaV0::removeDuplicatePairs(vector<pair<int, int>>& pairs) {
    map<int, int> pairMap;
    for (auto pair : pairs) {
        compareAndStorePair(pairMap, pair);
    }
    pairs.clear();
    for (auto it = pairMap.begin(); it != pairMap.end(); ++it) {
        pairs.push_back(std::make_pair(it->first, it->second));
    }
}




void AnaV0::removeDuplicates(vector<pair<int, int>>& vec, TVector3 vtx, double * beamline, float const bField) {
    unordered_map<int, int> maxValueMap; // Map to store the maximum value encountered for each key
    // Iterate over the vector
    vector<pair<int,int>> originalVec = vec;
    int i = 0;
    for (auto it = vec.begin(); it != vec.end(); ) {
        // Check if the current pair has a duplicate value
        if (maxValueMap.find(it->first) != maxValueMap.end() || maxValueMap.find(it->second) != maxValueMap.end()) {
            // If a duplicate value is found, compare the dcaDaughters and keep the pair with smaller dca


            const StUPCTrack* currentTrk1 = mUpcEvt->getTrack(it->first);
            const StUPCTrack* currentTrk2 = mUpcEvt->getTrack(it->second);

            pair<int,int> maxValue;
            if (maxValueMap.find(it->first) != maxValueMap.end() && maxValueMap.find(it->second) != maxValueMap.end()){
                maxValue.first = maxValueMap.find(it->first);
                maxValue.second = maxValueMap.find(it->second);
            }else if(maxValueMap.find(it->first) != maxValueMap.end()) {
                maxValue.first = maxValueMap.find(it->first);
                maxValue.second = originalVec[i].second;
            } else if (maxValueMap.find(it->second) != maxValueMap.end()){
                maxValue.second = maxValueMap.find(it->second);
                maxValue.first = originalVec[i].first;
            }


            const StUPCTrack* maxTrk1 = mUpcEvt->getTrack(maxValue.first);
            const StUPCTrack* maxTrk2 = mUpcEvt->getTrack(maxValue.second);
            cerr << "we good make"<< endl;

            StUPCV0 currentK0(currentTrk1,currentTrk2, mUtil->mass(hasGoodTPCnSigma(currentTrk1)), mUtil->mass(hasGoodTPCnSigma(currentTrk2)),1,2, vtx, beamline, bField, false);
            StUPCV0 maxK0(maxTrk1,maxTrk2, mUtil->mass(hasGoodTPCnSigma(maxTrk1)), mUtil->mass(hasGoodTPCnSigma(maxTrk2)),1,2, vtx, beamline, bField, false);

            if (currentK0.dcaDaughters() < maxK0.dcaDaughters()) {
                maxValueMap.erase(maxValue.first);
                maxValueMap.erase(maxValue.second);
                maxValueMap[it->first] = it->first;
                maxValueMap[it->second] = it->second;
            }

            // Remove the current pair
            it = vec.erase(it);
        }else {
            // If no duplicate value is found, update the values in the map
            maxValueMap[it->first] = it->first;
            maxValueMap[it->second] = it->second;
            ++it;
        }
    }
    ++i;
}

int AnaV0::findValuePosition(map<int, int> myMap, int value) {
    for (int i = 0; i < myMap.size(); ++i) {
        if (myMap[i].first == value || myMap[i].second == value)
            return i;
    }
    // If the value is not found, return -1 or some other sentinel value
    return -1;

}

*/
bool AnaV0Control::lambda(const StUPCTrack *trk1,const StUPCTrack *trk2){
    if (hasGoodTPCnSigma(trk1) == PROTON && trk1->getCharge() > 0)
        return true;
    else if(hasGoodTPCnSigma(trk2) == PROTON && trk2->getCharge() > 0)
        return true;
    else 
        return false;
}



bool AnaV0Control::oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2){
    Double_t charge;
    charge = trk1->getCharge() + trk2->getCharge();
    if(charge == 0)
        return true;
    else
        return false;
}



bool AnaV0Control::is2pions(const StUPCTrack *trk1,const StUPCTrack *trk2){

    if((hasGoodTPCnSigma(trk1) == PION) && (hasGoodTPCnSigma(trk2) == PION)){
        return true;
    }
    else
        return false;

}
bool AnaV0Control::isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2){

   if( hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PROTON ){
       return true;
   }else if(hasGoodTPCnSigma(trk1) == PROTON && hasGoodTPCnSigma(trk2) == PION){
       return true;     
   }
   else
       return false;

}


int AnaV0Control::hasGoodTPCnSigma(const StUPCTrack *trk){ //zmiernit podmienku na proton, skusit prebehnut aj na K0, Lambda 
  //check for good nSigma of hadron in TPC, 10 means unidentified
    if((trk->getNSigmasTPCProton() < 3) && (trk->getNSigmasTPCKaon() > 3) && (trk->getNSigmasTPCPion()>3))
        return PROTON;
    else if((trk->getNSigmasTPCProton() > 3) && (trk->getNSigmasTPCKaon() < 3) && (trk->getNSigmasTPCPion()>3))
        return KAON;
    else if(trk->getNSigmasTPCPion() < 3)
        return PION;
    else
        return 10;

}



void AnaV0Control::fillGoodTrackCuts(const StUPCTrack* trk){
    hDcaZ->Fill(trk->getDcaZ());
    hDcaXY->Fill(trk->getDcaXY());
    hNfitHits->Fill(trk->getNhitsFit());
    hNhitsDEdx->Fill(trk->getNhitsDEdx());
}


void AnaV0Control::saveNSigmaCorr(const StUPCTrack *trk){
    Double_t nSigmaPion, nSigmaProton, nSigmaKaon, nSigmaElectron;
    nSigmaPion = trk->getNSigmasTPCPion();
    nSigmaProton = trk->getNSigmasTPCProton();
    nSigmaKaon = trk->getNSigmasTPCKaon();
    nSigmaElectron = trk->getNSigmasTPCElectron();
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
void AnaV0Control::SaveMissingMomenta(TVector3 missP)
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