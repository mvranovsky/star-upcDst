#include "AnaV0SingleState.h"

AnaV0SingleState::AnaV0SingleState(TFile *outFile): Ana(outFile){}

void AnaV0SingleState::Make()
{
    
    vector<int> tracksV0;
    for (int itrk = 0; itrk < mUpcEvt->getNumberOfTracks(); ++itrk){


        hAnalysisFlow->Fill(V0ALL);
        if(!runMCAna && !CheckTriggers(&triggerID, mUpcEvt, hTriggerBits))
           return;
        hAnalysisFlow->Fill(V0TRIG);
        const StUPCTrack* trk = mUpcEvt->getTrack(itrk);
        
        if( trk->getFlag( StUPCTrack::kV0 ) && V0Control)
            continue;
        if( trk->getFlag( StUPCTrack::kPrimary ) && !V0Control)
            continue;
        if (trk->getFlag( StUPCTrack::kPrimary ) && trk->getFlag( StUPCTrack::kV0 )){
            cout << "We got a track with kPrimary and kV0 flag." << endl;
            hBothFlags->Fill(1);
            continue;
        }

        hAnalysisFlow->Fill(V0FLAG);

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
        if(abs(eta) > maxEta) 
            continue;
        if(abs(eta) < minEta) // avoid the peak in the middle
            continue;
        hAnalysisFlow->Fill(V0ETA);

        saveNSigmaCorr(trk);
        //particle identification: proton or pion
        if(!(hasGoodTPCnSigma(trk) == PROTON || hasGoodTPCnSigma(trk) == PION))
            continue;
        hAnalysisFlow->Fill(V0PID);

        tracksV0.push_back(itrk);

    }//global tracks loop
    mRecTree->setNGoodTpcTrks( tracksV0.size() );
    SaveEventInfo(mUpcEvt);
    
    if(!runMCAna){
        SaveTriggerInfo(mUpcEvt, mRpEvt);
    }

    //info about the beamline and mag field
    fillBeamlineInfo(mUpcEvt);


    TVector3 vertex, vertex0, primaryVertex;
    vertex0 = {0,0,0};
    int HADRON1, HADRON2, totalCharge;
    int vertexIdTrk1, vertexIdTrk2;
    Bool_t tofMatchTrk1, tofMatchTrk2;
    vector<pair<int, int>> hadronPairV0;

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
            // interested only in pi-pi or pi-p
            if (HADRON1 == PROTON && HADRON2 == PROTON)
                continue; 

            
            //if possibly both tracks originate from a upcDst vertex, calculate the difference between prodvertexhypo and upcDst vtx
            const StUPCVertex* vtx = mUpcEvt->getVertex(vertexIdTrk1);
            // same vertex
            double vertexDiff;
            bool primVtx;
            if(vtx && vertexIdTrk1 == vertexIdTrk2){            
                primaryVertex = {vtx->getPosX(),vtx->getPosY(), vtx->getPosZ()};
                StUPCV0 K0prim(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), tracksV0[itrk], tracksV0[jtrk], primaryVertex, beamline, bField, false);
                vertexDiff = K0prim.prodVertexHypo().Z();
                primVtx = true;

            } else{
                primVtx = false;


            }
            
            //define V0 pair
            StUPCV0 K0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), tracksV0[itrk], tracksV0[jtrk], vertex0, beamline, bField, false);
            vertex = K0.prodVertexHypo();


            //same cuts as were in StUPCSelectV0 except K0L1, K0L2 -> K0
            fillTopologyCutsBefore(K0);
            //cout << "beamline: " << beamline[0] << beamline[1] << beamline[2] << beamline[3] << " dca daughters: " << K0.dcaDaughters() << ", dca beamline: " << K0.DCABeamLine() << ", pointing angle: " << K0.pointingAngleHypo() << ", decay length: " << K0.decayLengthHypo() << endl;
            if ( !(K0.dcaDaughters() < maxDcaDaughters && K0.DCABeamLine() < maxDcaBeamLine && ( K0.pointingAngleHypo()> minPointingAngle || K0.decayLengthHypo()< maxDecayLengthHypo) ) )
                continue;
            fillTopologyCutsAfter(K0);
            
            
            // condition for difference between vtx and hypo vertex
            if(primVtx){
                vertexDiff = abs(vertexDiff - vertex.Z() );
                hVtxDiff->Fill(vertexDiff);
                hTOFTracks->Fill(1);
            }else{
                hTOFTracks->Fill(-1);
            }

            hAnalysisFlow->Fill(V0PAIR);

            //cut for z position of vertex- using prodvertexhypo
            hPosZ->Fill(vertex.Z());
            if (abs(vertex.Z()) > vertexRange )
                continue;

            hPosZCut->Fill(vertex.Z()); 
            hAnalysisFlow->Fill(V0ZVERTEX);
            hadronPairV0.push_back(make_pair(tracksV0[itrk],tracksV0[jtrk]));
        }//inner track loop
    }// outer track loop
    hNPairV0->Fill(hadronPairV0.size());

    if(hadronPairV0.size() == 0){
        return;
    }else if (hadronPairV0.size() > 1){
        // first filter if 1 track is not in 2 pairs
        hadronPairV0 = filterPairs(hadronPairV0);
        // if still more than 1 pair, then resize and keep only the first
        // tree can hold only 1 V0 per event
        if(hadronPairV0.size() > 1){
            hadronPairV0.resize(1);
        }
    } 

    //checking for pair made of one track
    if (hadronPairV0[0].first == hadronPairV0[0].second){
        hSameTrackPair->Fill(1);
        cout << "We got a pair with same track for both. " << endl;     
        return;
    }

    const StUPCTrack* trk1 = mUpcEvt->getTrack(hadronPairV0[0].first);
    const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronPairV0[0].second);

    HADRON1 = hasGoodTPCnSigma(trk1);
    HADRON2 = hasGoodTPCnSigma(trk2);

    hTotQ->Fill( trk1->getCharge() + trk2->getCharge() );
    if (oppositePair(trk1, trk2))
        hAnalysisFlow->Fill(V0OPPOSITE);


    totalCharge = trk1->getCharge() + trk2->getCharge(); 
    if (HADRON1 == HADRON2){
        mRecTree->setPairID(K0S,0);
        if(totalCharge == 0)
            hAnalysisFlow->Fill(V0PIPI);
    } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && lambda(trk1, trk2) ){ 
        mRecTree->setPairID(LAMBDA,0);
        if(totalCharge == 0)
            hAnalysisFlow->Fill(V0PPI);
    } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && !lambda(trk1, trk2) ){
        mRecTree->setPairID(LAMBDABAR,0);
        if(totalCharge == 0)
            hAnalysisFlow->Fill(V0PIPBAR); 
    } else return;


    StUPCV0 V0(trk1, trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronPairV0[0].first,hadronPairV0[0].second, vertex0, beamline, bField, false);  
    vertex = V0.prodVertexHypo();

    // saving info about vertex
    SaveVertexInfo(vertex, V0.dcaDaughters(), V0.DCABeamLine(), V0.pointingAngleHypo(), V0.decayLengthHypo(), -999 , 0);
    
    TLorentzVector state, hadron1, hadron2;
    state = V0.lorentzVector();
    hadron1 = V0.lorentzVectorPart1();
    hadron2 = V0.lorentzVectorPart2();


    // save info about tracks-> the first one always has match with ToF 
    if(trk1->getFlag( StUPCTrack::kTof ) && trk2->getFlag( StUPCTrack::kTof ) ){
        SaveTrackInfo(trk1,hadron1,1);
        SaveTrackInfo(trk2,hadron2,0);

    }else if( trk1->getFlag( StUPCTrack::kTof ) ){
        //trk1 = tag
        SaveTrackInfo(trk1, hadron1,0);
        //trk2 = probe
        SaveTrackInfo(trk2,hadron2,1);


    } else if(trk2->getFlag( StUPCTrack::kTof )){
        //trk2 = tag
        SaveTrackInfo(trk2, hadron2,0);  
        //trk1 = probe
        SaveTrackInfo(trk1, hadron1,1);

    }else return;

    // save info about state
    SaveStateInfo(state,totalCharge, 0);
    
    // checking for a problem 
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

void AnaV0SingleState::Init(){

    mUtil = new Util();


    if( DEBUG )
        cout<<"AnaV0SingleState::Init() called"<<endl;

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
    mRecTree = new RecTree(nameOfAnaV0SingleStateTree, AnaV0SingleStateTreeBits, false); 


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
    hNfitHits->GetXaxis()->SetTitle("N^{fit}_{hits} [-]");
    hNfitHits->GetYaxis()->SetTitle(YAxisDescription);
     
    hNhitsDEdx = new TH1D("hNhitsDEdx", "NhitsDEdx", 45, 5, 50);
    hNhitsDEdx->SetTitle("Distribution of number of hits dE/dx ");
    hNhitsDEdx->GetXaxis()->SetTitle("N^{dEdx}_{hits} [-]");
    hNhitsDEdx->GetYaxis()->SetTitle(YAxisDescription);
       
    hVtxDiff =  new TH1D("hVtxDiff", "Difference in Vtx_{z} between prodVertexHypo and UPCVertex_{z}", 60, 0, 30); 
    hVtxDiff->SetTitle("Difference in determination of z_{vertex}");
    hVtxDiff->GetXaxis()->SetTitle("Z_{vertexDiff} [cm]");
    hVtxDiff->GetYaxis()->SetTitle(YAxisDescription);

    hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}", 60, -150, 150); 
    hPosZ->SetTitle("Distribution of position of z_{vertex}");
    hPosZ->GetXaxis()->SetTitle("Vertex_{Z} [cm]");
    hPosZ->GetYaxis()->SetTitle(YAxisDescription);

    hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}", 60, -150, 150); 
    hPosZCut->SetTitle("Distribution of position of z_{vertex}");
    hPosZCut->GetXaxis()->SetTitle("Vertex_{Z} [cm]");
    hPosZCut->GetYaxis()->SetTitle(YAxisDescription);

    hTOFTracks = new TH1D("hGlobalVtx", "Number of tracks in TOF before cut", 4, -2, 2);
    hTOFTracks->SetTitle("Number of pairs w. global tracks");
    hTOFTracks->GetXaxis()->SetTitle("Number of global vtxs");
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

    hInvMassEta = new TH2F("hInvMassEta", "hInvMassEta; m_{#pi^{+} #pi^{-}} [GeV/c^{2}]; #eta [-]", 100, 0.2,1.2, 200, -1.,1. );
    hInvMassEta->SetTitle("correlation plot of invMass and eta");

    /*
    hEta1Tag = new TH1D("hEta1Tag", "hEtaTag; #eta_{tag} [-]; counts", 20, -1,1);
    hEta1Probe = new TH1D("hEta1Probe", "hEtaProbe; #eta_{probe} [-]; counts", 20, -1,1);

    hpT1Tag = new TH1D("hpT1Tag", "hpTTag; p_{T}^{tag} [GeV/c]; counts", 15, 0,1.5);
    hpT1Probe = new TH1D("hpT1Probe", "hpTProbe; p_{T}^{probe} [GeV/c]; counts", 15, 0,1.5);

    hPhi1Tag = new TH1D("hPhi1Tag", "hPhiTag; #phi_{tag} [GeV/c]; counts", 12, -3.15, 3.15);
    hPhi1Probe = new TH1D("hPhi1Probe", "hPhiProbe; #phi_{probe} [GeV/c]; counts", 12, -3.15, 3.15);

    hEtaPhi1Tag = new TH2F("hEtaPhi1Tag", "hEtaPhiTag; #eta_{tag} [-]; #phi_{tag} [-]", 20,-1,1,12,-3.15, 3.15);
    hEtaPhi1Probe = new TH2F("hEtaPhi1Probe", "hEtaPhiProbe; #eta_{probe} [-]; #phi_{probe} [-]", 20,-1,1,12,-3.15, 3.15);

    hEta2Tag = new TH1D("hEta2Tag", "hEtaTag; #eta_{tag} [-]; counts", 20, -1,1);
    hEta2Probe = new TH1D("hEta2Probe", "hEtaProbe; #eta_{probe} [-]; counts", 20, -1,1);

    hpT2Tag = new TH1D("hpT2Tag", "hpTTag; p_{T}^{tag} [GeV/c]; counts", 15, 0,1.5);
    hpT2Probe = new TH1D("hpT2Probe", "hpTProbe; p_{T}^{probe} [GeV/c]; counts", 15, 0,1.5);

    hPhi2Tag = new TH1D("hPhi2Tag", "hPhiTag; #phi_{tag} [GeV/c]; counts", 12, -3.15, 3.15);
    hPhi2Probe = new TH1D("hPhi2Probe", "hPhiProbe; #phi_{probe} [GeV/c]; counts", 12, -3.15, 3.15);

    hEtaPhi2Tag = new TH2F("hEtaPhi2Tag", "hEtaPhiTag; #eta_{tag} [-]; #phi_{tag} [-]", 20,-1,1,12,-3.15, 3.15);
    hEtaPhi2Probe = new TH2F("hEtaPhi2Probe", "hEtaPhiProbe; #eta_{probe} [-]; #phi_{probe} [-]", 20,-1,1,12,-3.15, 3.15);
    */
    //plots for nSigma around eta == 0

    hNSigmaPiProbe1 = new TH1D("hNSigmaPiProbe1", "hNSigmaPiProbe1; n#sigma [-];counts", 20,-10,10);  
    hNSigmaPiTag1 = new TH1D("hNSigmaPiTag1", "hNSigmaPiTag1; n#sigma [-];counts", 20,-10,10);  

    hNSigmaPProbe1 = new TH1D("hNSigmaPProbe1", "hNSigmaPProbe1; n#sigma [-];counts", 20,-10,10);  
    hNSigmaPTag1 = new TH1D("hNSigmaPTag1", "hNSigmaPTag1; n#sigma [-];counts", 20,-10,10);  

    hNSigmaKProbe1 = new TH1D("hNSigmaKProbe1", "hNSigmaKProbe1; n#sigma [-];counts", 20,-10,10);  
    hNSigmaKTag1 = new TH1D("hNSigmaKTag1", "hNSigmaKTag1; n#sigma [-];counts", 20,-10,10);  

    hNSigmaEProbe1 = new TH1D("hNSigmaEProbe1", "hNSigmaEProbe1; n#sigma [-];counts", 20,-10,10);  
    hNSigmaETag1 = new TH1D("hNSigmaETag1", "hNSigmaETag1; n#sigma [-];counts", 20,-10,10);  

    hNSigmaPiProbe2 = new TH1D("hNSigmaPiProbe2", "hNSigmaPiProbe2; n#sigma [-];counts", 20,-10,10);  
    hNSigmaPiTag2 = new TH1D("hNSigmaPiTag2", "hNSigmaPiTag2; n#sigma [-];counts", 20,-10,10);  

    hNSigmaPProbe2 = new TH1D("hNSigmaPProbe2", "hNSigmaPProbe2; n#sigma [-];counts", 20,-10,10);  
    hNSigmaPTag2 = new TH1D("hNSigmaPTag2", "hNSigmaPTag2; n#sigma [-];counts", 20,-10,10);  

    hNSigmaKProbe2 = new TH1D("hNSigmaKProbe2", "hNSigmaKProbe2; n#sigma [-];counts", 20,-10,10);  
    hNSigmaKTag2 = new TH1D("hNSigmaKTag2", "hNSigmaKTag2; n#sigma [-];counts", 20,-10,10);  

    hNSigmaEProbe2 = new TH1D("hNSigmaEProbe2", "hNSigmaEProbe2; n#sigma [-];counts", 20,-10,10);  
    hNSigmaETag2 = new TH1D("hNSigmaETag2", "hNSigmaETag2; n#sigma [-];counts", 20,-10,10);  


    hPtPhi1 = new TH2F("hPtPhiTagAround0", "hPtPhi; p_{T} [GeV/c]; #phi [-]", 180,0.2,2,120,-3.15, 3.15);
    hPtVz1 = new TH2F("hPtVzTagAround0", "hPtVz; p_{T} [GeV/c]; V_{Z} [cm]", 180,0.2,2,200,-100, 100);
    hVzPhi1 = new TH2F("hVzPhiTagAround0", "hVzPhi;  V_{Z} [cm]; #phi [-]", 200,-100, 100, 120, -3.15, 3.15);
    hInvMassPhi1 = new TH2F("hInvMassPhiTagAround0", "hInvMassPhi;  m_{#pi #pi} [GeV/c^{2}]; #phi [-]", 180,0.2,2, 120, -3.15, 3.15);
    hInvMassVz1 = new TH2F("hInvMassVzTagAround0", "hInvMassVz;  m_{#pi #pi} [GeV/c^{2}]; V_{Z} [cm]", 180,0.2,2, 200, -100, 100);
    hInvMassPt1 = new TH2F("hInvMassPtTagAround0", "hInvMassPt;  m_{#pi #pi} [GeV/c^{2}]; p_{T} [GeV/c]", 180,0.2,2, 180, 0.2, 2.);

    hPtPhi2 = new TH2F("hPtPhiProbeAround0", "hPtPhi; p_{T} [GeV/c]; #phi [-]", 180,0.2,2,120,-3.15, 3.15);
    hPtVz2 = new TH2F("hPtVzProbeAround0", "hPtVz; p_{T} [GeV/c]; V_{Z} [cm]", 180,0.2,2,200,-100, 100);
    hVzPhi2 = new TH2F("hVzPhiProbeAround0", "hVzPhi;  V_{Z} [cm]; #phi [-]", 200,-100, 100, 120, -3.15, 3.15);
    hInvMassPhi2 = new TH2F("hInvMassPhiProbeAround0", "hInvMassPhi;  m_{#pi #pi} [GeV/c^{2}]; #phi [-]", 180,0.2,2, 120, -3.15, 3.15);
    hInvMassVz2 = new TH2F("hInvMassVzProbeAround0", "hInvMassVz;  m_{#pi #pi} [GeV/c^{2}]; V_{Z} [cm]", 180,0.2,2, 200, -100, 100);
    hInvMassPt2 = new TH2F("hInvMassPtProbeAround0", "hInvMassPt;  m_{#pi #pi} [GeV/c^{2}]; p_{T} [GeV/c]", 180,0.2,2, 180, 0.2, 2.);
      
    hEtaPhiProbeYesToF = new TH2F("hEtaPhiProbeYesToFAround0", "hEtaPhiProbeYesToF; #eta_{Probe}^{w ToF} [-]; #phi_{Probe}^{w ToF} [-]", 40, -0.2, 0.2, 120, -3.15, 3.15);
    hEtaPhiProbeNoToF = new TH2F("hEtaPhiProbeNoToFAround0", "hEtaPhiProbeNoToF; #eta_{Probe}^{w/o ToF} [-]; #phi_{Probe}^{w/o ToF} [-]", 40, -0.2, 0.2, 120, -3.15, 3.15);

    hEtaPhiProbeYesToFLarge = new TH2F("hEtaPhiProbeYesToF", "hEtaPhiProbeYesToF; #eta_{Probe}^{w ToF} [-]; #phi_{Probe}^{w ToF} [-]", 180, -0.9, 0.9, 120, -3.15, 3.15);
    hEtaPhiProbeNoToFLarge = new TH2F("hEtaPhiProbeNoToF", "hEtaPhiProbeNoToF; #eta_{Probe}^{w/o ToF} [-]; #phi_{Probe}^{w/o ToF} [-]", 180, -0.9, 0.9, 120, -3.15, 3.15);


}

void AnaV0SingleState::fillTag1(const StUPCTrack* trk, Double_t pT, Double_t Vz, Double_t invMass) {

    if(abs(trk->getEta() ) < 0.2){
        hPtPhi1->Fill(pT , trk->getPhi() );
        hPtVz1->Fill(pT, Vz);
        hVzPhi1->Fill(Vz, trk->getPhi() );
        hInvMassPhi1->Fill( invMass,trk->getPhi() );
        hInvMassVz1->Fill( invMass, Vz);
        hInvMassPt1->Fill( invMass, pT);

    }

}

void AnaV0SingleState::fillProbe1(const StUPCTrack* trk, Double_t pT, Double_t Vz, Double_t invMass) {

    if(abs(trk->getEta() ) < 0.2){
        hPtPhi2->Fill(pT , trk->getPhi() );
        hPtVz2->Fill(pT, Vz);
        hVzPhi2->Fill(Vz, trk->getPhi() );
        hInvMassPhi2->Fill( invMass,trk->getPhi() );
        hInvMassVz2->Fill( invMass, Vz);
        hInvMassPt2->Fill( invMass, pT);
    }
}
/*
void AnaV0SingleState::fillTag2(const StUPCTrack* trk) {


    if(abs(trk->getEta()) <0.2 ){
        hEta2Tag->Fill( trk->getEta() );
        hPhi2Tag->Fill( trk->getPhi() );
        hpT2Tag->Fill( trk->getPt() );
        hEtaPhi2Tag->Fill(trk->getEta() , trk->getPhi() );
        hNSigmaPiTag2->Fill(trk->getNSigmasTPCPion());
        hNSigmaPTag2->Fill(trk->getNSigmasTPCProton());
        hNSigmaKTag2->Fill(trk->getNSigmasTPCKaon());
        hNSigmaETag2->Fill(trk->getNSigmasTPCElectron());
    }

}

void AnaV0SingleState::fillProbe2(const StUPCTrack* trk) {
    

    if(abs(trk->getEta()) <0.2 ){
        hEta2Probe->Fill( trk->getEta() );
        hPhi2Probe->Fill( trk->getPhi() );
        hpT2Probe->Fill( trk->getPt() );
        hEtaPhi2Probe->Fill(trk->getEta() , trk->getPhi() );
        hNSigmaPiProbe2->Fill(trk->getNSigmasTPCPion());
        hNSigmaPProbe2->Fill(trk->getNSigmasTPCProton());
        hNSigmaKProbe2->Fill(trk->getNSigmasTPCKaon());
        hNSigmaEProbe2->Fill(trk->getNSigmasTPCElectron());
    }

}
*/
void AnaV0SingleState::fillBeamlineInfo(const StUPCEvent *mUpcEvt){

    if(!runMCAna){
        bField = mUpcEvt->getMagneticField();
        beamline[0] = mUpcEvt->getBeamXPosition();
        beamline[2] = mUpcEvt->getBeamXSlope();
        beamline[1] = mUpcEvt->getBeamYPosition();
        beamline[3] = mUpcEvt->getBeamYSlope();
    }else {
        bField = -4.991; // guess based on real runs
        beamline[0] = 0;
        beamline[2] = 0;
        beamline[1] = 0;
        beamline[3] = 0;
    }

    return;
}



void AnaV0SingleState::fillFinalPlots() {

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



void AnaV0SingleState::fillTopologyCutsBefore(const StUPCV0& V0){

    hDcaDaughters->Fill(V0.dcaDaughters());
    hDcaBeamline->Fill(V0.DCABeamLine());
    hPointingAngle->Fill(V0.pointingAngleHypo());
    hDecayLength->Fill(V0.decayLengthHypo());
    hDecayLPointingA->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());

}


void AnaV0SingleState::fillTopologyCutsAfter(const StUPCV0& V0){

    hDcaDaughtersCut->Fill(V0.dcaDaughters());
    hDcaBeamlineCut->Fill(V0.DCABeamLine());
    hPointingAngleCut->Fill(V0.pointingAngleHypo());
    hDecayLengthCut->Fill(V0.decayLengthHypo());
    hDecayLPointingACut->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());
}

bool AnaV0SingleState::shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b) {
    TVector3 vtx0 = {0,0,0};

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

vector<pair<int, int>> AnaV0SingleState::filterPairs(vector<pair<int, int>>& pairs) {
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

bool AnaV0SingleState::lambda(const StUPCTrack *trk1,const StUPCTrack *trk2){
    if (hasGoodTPCnSigma(trk1) == PROTON && trk1->getCharge() > 0)
        return true;
    else if(hasGoodTPCnSigma(trk2) == PROTON && trk2->getCharge() > 0)
        return true;
    else 
        return false;
}


bool AnaV0SingleState::oppositePair(const StUPCTrack *trk1,const StUPCTrack *trk2){
    Double_t charge;
    charge = trk1->getCharge() + trk2->getCharge();
    if(charge == 0)
        return true;
    else
        return false;
}

bool AnaV0SingleState::is2pions(const StUPCTrack *trk1,const StUPCTrack *trk2){

    if((hasGoodTPCnSigma(trk1) == PION) && (hasGoodTPCnSigma(trk2) == PION)){
        return true;
    }
    else
        return false;

}
bool AnaV0SingleState::isProtonPion(const StUPCTrack *trk1,const StUPCTrack *trk2){

   if( hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PROTON ){
       return true;
   }else if(hasGoodTPCnSigma(trk1) == PROTON && hasGoodTPCnSigma(trk2) == PION){
       return true;     
   }
   else
       return false;
}

int AnaV0SingleState::hasGoodTPCnSigma(const StUPCTrack *trk){ //zmiernit podmienku na proton, skusit prebehnut aj na K0, Lambda 
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

void AnaV0SingleState::fillGoodTrackCuts(const StUPCTrack* trk){
    //hDcaZ->Fill(trk->getDcaZ());
    //hDcaXY->Fill(trk->getDcaXY());
    hNfitHits->Fill(trk->getNhitsFit());
    hNhitsDEdx->Fill(trk->getNhitsDEdx());
    hPt->Fill(trk->getPt());
}


void AnaV0SingleState::saveNSigmaCorr(const StUPCTrack *trk){
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
void AnaV0SingleState::CalculateTOFEff(unsigned int tagID)
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

void AnaV0SingleState::CalculatePID()
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

Int_t AnaV0SingleState::PIDStartegy(int strategy)
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

bool AnaV0SingleState::IsPairOf(int type)
{
   // 0 - PION
   // 1 - KAON
   // 2 - PROTON
   if( type < 0 || type > 2)
      return false;

   return (mRecTree->getNSigmaTPC(0, type) < 2 && mRecTree->getNSigmaTPC(1, type) < 2);

}
*/
void AnaV0SingleState::SaveMissingMomenta(TVector3 missP)
{
   mRecTree->setPtMissing( missP.Pt() );
   mRecTree->setPxMissing( missP.X() );
   mRecTree->setPyMissing( missP.Y() );   
}

/*
void AnaV0SingleState::FillMSquared()
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
