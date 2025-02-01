#include "AnaV0Mult.h"

AnaV0Mult::AnaV0Mult(TFile *outFile): Ana(outFile){}

void AnaV0Mult::Make()
{
    // neaktualna verzia, 
   for (int itrk = 0; itrk < mUpcEvt->getNumberOfTracks(); ++itrk){
       const StUPCTrack* trk = mUpcEvt->getTrack(itrk);
       hAnalysisFlow->Fill(V0ALL);
       if(!CheckTriggers(&triggerID, mUpcEvt, hTriggerBits))
          return;
       hAnalysisFlow->Fill(V0TRIG);
       
       if(TOF2Tracks && !(trk->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk) ) )
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
       //hAnalysisFlow->Fill(V0ETA);


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
           hadronID.push_back(itrk);
        } else if(trk->getFlag( StUPCTrack::kPrimary ) && V0Control ){
           hadronID.push_back(itrk);
        } else {
           continue; 
        }

        hAnalysisFlow->Fill(V0FLAG);

        if(trk->getFlag(StUPCTrack::kTof) && IsGoodTofTrack(trk) )
            tagID.push_back(itrk);

   }//global tracks loop

   mRecTree->setNGoodTpcTrks( hadronID.size() );
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
   for (unsigned int itrk = 0; itrk < hadronID.size(); ++itrk){
       const StUPCTrack* trk1 = mUpcEvt->getTrack(hadronID[itrk]);
       tofMatchTrk1 = (trk1->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk1) );
       vertexIdTrk1 = trk1->getVertexId();
       HADRON1 = hasGoodTPCnSigma(trk1);

       //inner loop of tracks
       for (unsigned int jtrk = itrk+1; jtrk < hadronID.size(); ++jtrk){

           const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronID[jtrk]);
           tofMatchTrk2 = (trk2->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk2) );
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
           StUPCV0 K0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronID[itrk], hadronID[jtrk], primaryVertex, beamline, bField, false);
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
           //AnalysisFlow->Fill(V0ZVERTEX);
           hadronPairV0.push_back(make_pair(hadronID[itrk],hadronID[jtrk]));
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
      //SaveVertexInfo(vertex, V0.dcaDaughters(), V0.DCABeamLine(), V0.pointingAngleHypo(), V0.decayLengthHypo(), vertexDiff , j);
      SaveTrackInfo(trk1,hadron1,j*2);
      SaveTrackInfo(trk2,hadron2,j*2 +1);
      
        // save info about tracks-> the first one always has match with ToF 
      if( trk1->getFlag( StUPCTrack::kTof ) ){
          //trk1 = tag
          SaveTrackInfo(trk1,hadron1,2*j);
          //trk2 = probe
          SaveTrackInfo(trk2,hadron2,2*j+1);
      } else{
          //trk2 = tag
          SaveTrackInfo(trk2,hadron2,2*j);  
          //trk1 = probe
          SaveTrackInfo(trk1,hadron1,2*j+1);
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

   if( DEBUG ){
    cout << "Finished AnaV0Mult::Make()" << endl;
   }

   // end with clearing all the variables for rec tree 
   resetInfo();


}//end of make()

void AnaV0Mult::Init()
{
    //initialize class for selecting V0 candidates 
    //mSelectV0 = new StUPCSelectV0Modified();

    //init database to read beam-line parameters
    //mDbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    //mDbMk->Init();

   mUtil = new Util();


   if( DEBUG )
      cout<<"AnaV0Mult::Init() called"<<endl;

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
   mRecTree = new RecTree(nameOfAnaV0MultTree, AnaV0MultTreeBits, false); 

   hEta = new TH1D("hEta", "Pseudorapidity; #eta; counts", 60, -2, 2);
   hEtaCut = new TH1D("hEtaCut", "Pseudorapidity; #eta [-]; counts", 60, -2, 2);

   hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks; #eta; #varphi",100,-2,2,100,-4,4);
   hEtaPhiCut = new TH2F("hEtaPhiCut","Phi vs eta of TOF tracks; #eta ; #varphi",100,-2,2,100,-4,4);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 
   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 

   hEtaVtxZ = new TH2F("hEtaVtxZ", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);
   hEtaVtxZCut = new TH2F("hEtaVtxZCut", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);



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


}


void AnaV0Mult::fillFinalPlots() {

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


void AnaV0Mult::fillEtaVtxPlotsBefore(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ){

   Double_t eta2 = trk2->getEta();
   Double_t phi2 = trk2->getPhi();

   Double_t eta1 = trk1->getEta();
   Double_t phi1 = trk1->getPhi();

   hPosZ->Fill(posZ);

   hEtaVtxZ->Fill(eta1 , posZ );
   hEtaVtxZ->Fill(eta2 , posZ );

   hEtaPhi->Fill(eta1, phi1);
   hEtaPhi->Fill(eta2, phi2);

   hEta->Fill(eta1);
   hEta->Fill(eta2);

}

void AnaV0Mult::fillEtaVtxPlotsAfter(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ){

   Double_t eta2 = trk2->getEta();
   Double_t phi2 = trk2->getPhi();

   Double_t eta1 = trk1->getEta();
   Double_t phi1 = trk1->getPhi();

   hPosZCut->Fill(posZ);

   hEtaVtxZCut->Fill(eta1 , posZ );
   hEtaVtxZCut->Fill(eta2 , posZ );

   hEtaPhiCut->Fill(eta1, phi1);
   hEtaPhiCut->Fill(eta2, phi2);

   hEtaCut->Fill(eta1);
   hEtaCut->Fill(eta2);
}


void AnaV0Mult::fillTopologyCutsBefore(const StUPCV0& V0){

    hDcaDaughters->Fill(V0.dcaDaughters());
    hDcaBeamline->Fill(V0.DCABeamLine());
    hPointingAngle->Fill(V0.pointingAngleHypo());
    hDecayLength->Fill(V0.decayLengthHypo());
    hDecayLPointingA->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());

}


void AnaV0Mult::fillTopologyCutsAfter(const StUPCV0& V0){

    hDcaDaughtersCut->Fill(V0.dcaDaughters());
    hDcaBeamlineCut->Fill(V0.DCABeamLine());
    hPointingAngleCut->Fill(V0.pointingAngleHypo());
    hDecayLengthCut->Fill(V0.decayLengthHypo());
    hDecayLPointingACut->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());
}

bool AnaV0Mult::shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b) {
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

vector<pair<int, int>> AnaV0Mult::filterPairs(vector<pair<int, int>>& pairs) {
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

bool AnaV0Mult::lambda(const StUPCTrack *trk1,const StUPCTrack *trk2){
    if (hasGoodTPCnSigma(trk1) == PROTON && trk1->getCharge() > 0)
        return true;
    else if(hasGoodTPCnSigma(trk2) == PROTON && trk2->getCharge() > 0)
        return true;
    else 
        return false;
}



bool AnaV0Mult::oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2){
    Double_t charge;
    charge = trk1->getCharge() + trk2->getCharge();
    if(charge == 0)
        return true;
    else
        return false;
}

bool AnaV0Mult::is2pions(const StUPCTrack *trk1,const StUPCTrack *trk2){

    if((hasGoodTPCnSigma(trk1) == PION) && (hasGoodTPCnSigma(trk2) == PION)){
        return true;
    }
    else
        return false;

}
bool AnaV0Mult::isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2){

   if( hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PROTON ){
       return true;
   }else if(hasGoodTPCnSigma(trk1) == PROTON && hasGoodTPCnSigma(trk2) == PION){
       return true;     
   }
   else
       return false;
}

int AnaV0Mult::hasGoodTPCnSigma(const StUPCTrack *trk){ //zmiernit podmienku na proton, skusit prebehnut aj na K0, Lambda 
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

void AnaV0Mult::fillGoodTrackCuts(const StUPCTrack* trk){
    hDcaZ->Fill(trk->getDcaZ());
    hDcaXY->Fill(trk->getDcaXY());
    hNfitHits->Fill(trk->getNhitsFit());
    hNhitsDEdx->Fill(trk->getNhitsDEdx());
    hPt->Fill(trk->getPt());
}


void AnaV0Mult::saveNSigmaCorr(const StUPCTrack *trk){
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

void AnaV0Mult::SaveMissingMomenta(TVector3 missP)
{
   mRecTree->setPtMissing( missP.Pt() );
   mRecTree->setPxMissing( missP.X() );
   mRecTree->setPyMissing( missP.Y() );   
}
