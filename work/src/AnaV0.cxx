#include "AnaV0.h"

AnaV0::AnaV0(TFile *outFile): Ana(outFile){}

void AnaV0::Make()
{
    hAnalysisFlow->Fill(V0ALL);

    if(!runMCAna && !CheckTriggers(&CEPtriggers, mUpcEvt, hTriggerBits))
       return;
    hAnalysisFlow->Fill(V0TRIG);


    for (int itrk = 0; itrk < mUpcEvt->getNumberOfTracks(); ++itrk){

        const StUPCTrack* trk = mUpcEvt->getTrack(itrk);

        
        if( trk->getFlag( StUPCTrack::kV0 ) && V0Control)
            continue;
        if( trk->getFlag( StUPCTrack::kPrimary ) && !V0Control)
            continue;

        hAnalysisFlow->Fill(V0FLAG);

        if(TOF2Tracks && !(trk->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk) ) )
            continue;


        //conditions for track quality
        fillTrackQualityCuts(trk);
        if (trk->getNhitsFit() < minNHitsFit)
            continue;
        if (trk->getNhitsDEdx() < minNHitsDEdx)
            continue;
        if (trk->getPt() < minPt)
            continue; 


        fillNSigmaPlots(trk);
        //particle identification: proton or pion
        if(!(hasGoodTPCnSigma(trk) == PROTON || hasGoodTPCnSigma(trk) == PION))
            continue;

        if(trk->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk) )
            tagID.push_back(itrk);

        hAnalysisFlow->Fill(V0PID);

        hadronID.push_back(itrk);

    }//global tracks loop
    hNTracksTof->Fill( tagID.size() );
    hNTracksTpc->Fill( hadronID.size() );


    SaveEventInfo(mUpcEvt);
    
    if(!runMCAna){
        SaveTriggerInfo(mUpcEvt, mRpEvt);
    }

    //info about the beamline and mag field
    fillBeamlineInfo(mUpcEvt);


    TVector3 vertex, primaryVertex;
    TVector3 vertex0 = {0,0,0};
    int HADRON1, HADRON2, totalCharge;
    int vertexIdTrk1, vertexIdTrk2;
    Bool_t tofMatchTrk1, tofMatchTrk2;
    vector<pair<int, int>> hadronPairV0;
    StUPCV0 *V0;

    // outer loop of tracks
    for (unsigned int itrk = 0; itrk < hadronID.size(); ++itrk){

        const StUPCTrack* trk1 = mUpcEvt->getTrack(hadronID[itrk]);
        tofMatchTrk1 = trk1->getFlag( StUPCTrack::kTof );
        vertexIdTrk1 = trk1->getVertexId();
        HADRON1 = hasGoodTPCnSigma(trk1);

        //inner loop of tracks
        for (unsigned int jtrk = itrk+1; jtrk < hadronID.size(); ++jtrk){

            const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronID[jtrk]);
            tofMatchTrk2 = trk2->getFlag( StUPCTrack::kTof );
            vertexIdTrk2 = trk2->getVertexId();
            HADRON2 = hasGoodTPCnSigma(trk2);

            // checking for at least 1 tof match
            if(!(tofMatchTrk1 && IsGoodTofTrack(trk1)) && !(tofMatchTrk2 && IsGoodTofTrack(trk2)) )
                continue;
            // interested only in pi-pi or pi-p/p - pi
            if (HADRON1 == PROTON && HADRON2 == PROTON)
                continue; 


            // find out if both tracks belong to the same primary vtx
            const StUPCVertex* vtx = trk1->getVertex();
            double vertexDiff;
            if(usePrimVtx && !(vtx && vertexIdTrk1 == vertexIdTrk2) ){
                continue;
            } else if(usePrimVtx && (vtx && vertexIdTrk1 == vertexIdTrk2) ){
                primaryVertex = vtx->getPosVtx();
                V0 = new StUPCV0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronID[itrk], hadronID[jtrk], primaryVertex, beamline, bField, false);
                StUPCV0 V0alt(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronID[itrk], hadronID[jtrk],vertex0, beamline, bField, false);
                hVtxDiff->Fill( abs(V0->prodVertexHypo().Z() - V0alt.prodVertexHypo().Z() ) );
                vertexDiff = abs(V0->prodVertexHypo().Z() - V0alt.prodVertexHypo().Z() );
            }else if(!usePrimVtx && (vtx && vertexIdTrk1 == vertexIdTrk2) ){
                primaryVertex = vtx->getPosVtx();
                V0 = new StUPCV0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronID[itrk], hadronID[jtrk], vertex0, beamline, bField, false);
                StUPCV0 V0alt(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronID[itrk], hadronID[jtrk], primaryVertex, beamline, bField, false);
                hVtxDiff->Fill( abs(V0->prodVertexHypo().Z() - V0alt.prodVertexHypo().Z() ) );
                vertexDiff = abs(V0->prodVertexHypo().Z() - V0alt.prodVertexHypo().Z() );
            }else{
                V0 = new StUPCV0(trk1,trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronID[itrk], hadronID[jtrk], vertex0, beamline, bField, false);
                vertexDiff = -1;
            }


            //same cuts as were in StUPCSelectV0 except K0L1, K0L2 -> K0
            fillTopologyCutsBefore(*V0);
            if ( !(V0->dcaDaughters() < maxDcaDaughters && V0->DCABeamLine() < maxDcaBeamLine && ( V0->pointingAngleHypo()> minPointingAngle || V0->decayLengthHypo()< maxDecayLengthHypo) ) )
                continue;
            fillTopologyCutsAfter(*V0);
            
            hVtxDiffAfter->Fill(vertexDiff);
            TVector3 vertex = V0->prodVertexHypo();
            hAnalysisFlow->Fill(V0PAIR);


            hadronPairV0.push_back(make_pair(hadronID[itrk],hadronID[jtrk]));
        }//inner track loop
    }// outer track loop


    //special eta-vtxZ cut with filling of control plots
    vector<pair<int, int>> hadronPairV0After;

    for (unsigned int i = 0; i < hadronPairV0.size(); ++i){
        const StUPCTrack* trk1 = mUpcEvt->getTrack(hadronPairV0[0].first);
        const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronPairV0[0].second);


        Double_t eta2 = trk2->getEta();
        Double_t phi2 = trk2->getPhi();

        Double_t eta1 = trk1->getEta();
        Double_t phi1 = trk1->getPhi();

        hPosZ->Fill(vertex.Z());
        hEtaVtxZ->Fill(eta1 , vertex.Z() );
        hEtaPhi->Fill(eta1, phi1);
        hEta->Fill(eta1);

        hEtaVtxZ->Fill(eta2 , vertex.Z() );
        hEtaPhi->Fill(eta2, phi2);
        hEta->Fill(eta2);

        if(abs(trk1->getEta() )> maxEta || abs(trk2->getEta() ) > maxEta || abs(vertex.Z()) > vertexRange )
            continue;

        hPosZCut->Fill( vertex.Z() ); 
        hEtaCut->Fill( eta1 );
        hEtaCut->Fill( eta2 );
        hEtaVtxZCut->Fill( eta1, vertex.Z() );
        hEtaVtxZCut->Fill( eta2, vertex.Z() );

        hAnalysisFlow->Fill( V0ETAVTXZ );

        totalCharge = trk1->getCharge() + trk2->getCharge(); 
        hTotQ->Fill( totalCharge );

        if(totalCharge == 0 && hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PION){
            hadronPairV0After.push_back(hadronPairV0[i]);  
            if(trk2->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk2) && trk1->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk1) ){
               hInvMassTof2->Fill( V0->lorentzVector().M() );
               hInvMassTof2->Fill( V0->lorentzVector().M() );
               hInvMassTof1->Fill( V0->lorentzVector().M() );
               hInvMassTof1->Fill( V0->lorentzVector().M() );  
            }
            else{
               hInvMassTof1->Fill( V0->lorentzVector().M() );
            }
        }


    }

    hNPairV0->Fill(hadronPairV0After.size());

    pair<int,int> finalPair;

    if(hadronPairV0After.size() >= 1){
        finalPair = resizePairs( hadronPairV0After );
    }else{ 
        return;
    }


    const StUPCTrack* trk1 = mUpcEvt->getTrack(finalPair.first);
    const StUPCTrack* trk2 = mUpcEvt->getTrack(finalPair.second);

    HADRON1 = hasGoodTPCnSigma(trk1);
    HADRON2 = hasGoodTPCnSigma(trk2);



    if (HADRON1 == HADRON2){
        mRecTree->setPairID(K0S,0);
        if(totalCharge == 0){
            hAnalysisFlow->Fill(V0OPPOSITE);
            hAnalysisFlow->Fill(V0PIPI);
        }
    } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && lambda(trk1, trk2) ){ 
        mRecTree->setPairID(LAMBDA,0);
        if(totalCharge == 0){
            hAnalysisFlow->Fill(V0OPPOSITE);
            hAnalysisFlow->Fill(V0PPI);
        }
    } else if((HADRON1 == PROTON || HADRON2 == PROTON ) && !lambda(trk1, trk2) ){
        mRecTree->setPairID(LAMBDABAR,0);
        if(totalCharge == 0){
            hAnalysisFlow->Fill(V0OPPOSITE);
            hAnalysisFlow->Fill(V0PIPBAR); 
        }
    } else return;


    if(usePrimVtx){
        V0 = new StUPCV0(trk1, trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronPairV0[0].first,hadronPairV0[0].second, trk1->getVertex()->getPosVtx(), beamline, bField, false);  
    }else{
        V0 = new StUPCV0(trk1, trk2, mUtil->mass(HADRON1), mUtil->mass(HADRON2), hadronPairV0[0].first,hadronPairV0[0].second, vertex0, beamline, bField, false);  
    }
    vertex = V0->prodVertexHypo();

    // saving info about vertex
    SaveVertexInfo(V0, 0);
    //fillPrimVtxInfo(trk1, trk2);

    TLorentzVector state, hadron1, hadron2;
    state = V0->lorentzVector();
    hadron1 = V0->lorentzVectorPart1();
    hadron2 = V0->lorentzVectorPart2();


    if(trk1->getCharge() + trk2->getCharge() == 0 && hasGoodTPCnSigma(trk1) == PION && hasGoodTPCnSigma(trk2) == PION){
        if(trk2->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk2) && trk1->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk1) ){
           hInvMassTof2AfterPicking1->Fill( V0->lorentzVector().M() );
           hInvMassTof2AfterPicking1->Fill( V0->lorentzVector().M() );
           hInvMassTof1AfterPicking1->Fill( V0->lorentzVector().M() );
           hInvMassTof1AfterPicking1->Fill( V0->lorentzVector().M() );      
        }
        else{
           hInvMassTof1AfterPicking1->Fill( V0->lorentzVector().M() );
        }
    }


    // save info about tracks-> the first one always has match with ToF 
    if(trk1->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk1) && trk2->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk2) ){
        SaveTrackInfo(trk1,hadron1,1);
        SaveTrackInfo(trk2,hadron2,0);

    }else if( trk1->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk1)){
        //trk1 = tag
        SaveTrackInfo(trk1, hadron1,0);
        //trk2 = probe
        SaveTrackInfo(trk2,hadron2,1);


    } else if(trk2->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk2)){
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

    if(DEBUG){
        cout << "Finished AnaV0::Make()" << endl;
    }
}//end of make()

void AnaV0::Init(){

    mUtil = new Util();


    if( DEBUG )
        cout<<"AnaV0::Init() called"<<endl;

    mOutFile->cd();
    mOutFile->mkdir(nameOfAnaV0Dir);
    mOutFile->cd(nameOfAnaV0Dir);
    
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


    hNSigmaPiPcorr = new TH2F("hNSigmaPiPcorr","n_{#sigma} pions against protons;n#sigma_{#pi};n#sigma_{p}", 100, -25, 25, 100, -25, 25);

    hNSigmaPiKcorr = new TH2F("hNSigmaPiKcorr","n_{#sigma} pions against kaons;n#sigma_{#pi};n#sigma_{K}", 100, -25, 25, 100, -25, 25);

    hNSigmaPiecorr = new TH2F("hNSigmaPiecorr","n_{#sigma} pions against electrons;n#sigma_{#pi};n#sigma_{e}", 100, -25, 25, 100, -25, 25);

    hNSigmaPPicorr = new TH2F("hNSigmaPPicorr","n_{#sigma} protons against pions;n#sigma_{p};n#sigma_{#pi}", 100, -25, 25, 100, -25, 25);

    hNSigmaPKcorr = new TH2F("hNSigmaPKcorr","n_{#sigma} protons against kaons;n#sigma_{p};n#sigma_{K}", 100, -25, 25, 100, -25, 25);

    hNSigmaPecorr = new TH2F("hNSigmaPecorr","n_{#sigma} protons against electrons;n#sigma_{p};n#sigma_{e}", 100, -25, 25, 100, -25, 25);

    hNSigmaKPicorr = new TH2F("hNSigmaKPicorr","n_{#sigma} kaons against pions;n#sigma_{K};n#sigma_{#pi}", 100, -25, 25, 100, -25, 25);

    hNSigmaKPcorr = new TH2F("hNSigmaKPcorr","n_{#sigma} kaons against protons;n#sigma_{K};n#sigma_{p}", 100, -25, 25, 100, -25, 25);

    hNSigmaKecorr = new TH2F("hNSigmaKecorr","n_{#sigma} kaons against electrons;n#sigma_{K};n#sigma_{e}", 100, -25, 25, 100, -25, 25);

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

    hEtaCut = new TH1D("hEtaCut", "Pseudorapidity; #eta [-]; counts", 60, -2, 2);

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
       
    hVtxDiff =  new TH1D("hVtxDiff", "Difference in Vtx_{z} between prodVertexHypo and UPCVertex_{z}", 12, -2, 10); 
    hVtxDiff->SetTitle("Difference in determination of z_{vertex}");
    hVtxDiff->GetXaxis()->SetTitle("Z_{vertexDiff} [cm]");
    hVtxDiff->GetYaxis()->SetTitle(YAxisDescription);

    hVtxDiffAfter =  new TH1D("hVtxDiffAfter", "Difference in Vtx_{z} between prodVertexHypo and UPCVertex_{z}", 12, -2, 10);

    hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}", 60, -150, 150); 
    hPosZ->SetTitle("Distribution of position of z_{vertex}");
    hPosZ->GetXaxis()->SetTitle("Vertex_{Z} [cm]");
    hPosZ->GetYaxis()->SetTitle(YAxisDescription);

    hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}", 60, -150, 150); 
    hPosZCut->SetTitle("Distribution of position of z_{vertex}");
    hPosZCut->GetXaxis()->SetTitle("Vertex_{Z} [cm]");
    hPosZCut->GetYaxis()->SetTitle(YAxisDescription);

    hHasPrimVtx = new TH1D("hHasPrimVtx", "boolean histogram", 3, -1.5, 1.5);
    hHasPrimVtx->SetTitle("Number of pairs w. global tracks");
    hHasPrimVtx->GetXaxis()->SetTitle("1=hasPrimVtx, -1=noPrimVtx");
    hHasPrimVtx->GetYaxis()->SetTitle(YAxisDescription);
     
    hNVertices = new TH1D("hNVertices", "Number of vertices before cut", 100, 0, 10);
    hNVertices->SetTitle("Number of vertices");
    hNVertices->GetXaxis()->SetTitle("Number of vertices");
    hNVertices->GetYaxis()->SetTitle(YAxisDescription);

    hNPairV0 = new TH1D("hNPairV0", "Number of pairs after cuts", 5, -0.5, 4.5);
    hNPairV0->SetTitle("Number of V0 pairs");
    hNPairV0->GetXaxis()->SetTitle("Number of V0 pairs");
    hNPairV0->GetYaxis()->SetTitle(YAxisDescription);

    hTotQ = new TH1D("hTotQ", "Total charge of pair", 3, -1.5, 1.5);
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

    hEtaVtxZ = new TH2F("hEtaVtxZ", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);

    hEtaVtxZCut = new TH2F("hEtaVtxZCut", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);

    hNTracksTof = new TH1D("hNTracksTof", "hNTracksTof; Number of ToF tracks [-]; counts", 11, -0.5, 10.5);

    hNTracksTpc = new TH1D("hNTracksTpc", "hNTracksTpc; Number of TPC tracks [-]; counts", 26, -0.5, 25.5);

    hInvMassTof1 = new TH1D("hInvMassTof1", "hInvMassTof1; m_{#pi^{+} #pi^{-} }; counts",90 , 0.45, 0.54);
    hInvMassTof2 = new TH1D("hInvMassTof2", "hInvMassTof2; m_{#pi^{+} #pi^{-} }; counts",90 , 0.45, 0.54);

    hInvMassTof1AfterPicking1 = new TH1D("hInvMassTof1AfterPicking1", "hInvMassTof1; m_{#pi^{+} #pi^{-} }; counts",90 , 0.45, 0.54);
    hInvMassTof2AfterPicking1 = new TH1D("hInvMassTof2AfterPicking1", "hInvMassTof2; m_{#pi^{+} #pi^{-} }; counts",90 , 0.45, 0.54);

}

void AnaV0::fillTrackQualityCuts(const StUPCTrack* trk){

   hDcaZ->Fill( trk->getDcaZ() );
   hDcaXY->Fill( trk->getDcaXY() );
   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
}



void AnaV0::fillNSigmaPlots(const StUPCTrack *trk){
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
void AnaV0::fillEtaVtxPlots(const StUPCTrack *trk1, const StUPCTrack *trk2, double posZ){

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

void AnaV0::fillBeamlineInfo(const StUPCEvent *mUpcEvt){

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


pair<int,int> AnaV0::resizePairs( vector<pair<int,int>> hadronPairV0 ){


    // first filter if 1 track is not in 2 pairs
    if(hadronPairV0.size() > 1){
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
        return make_pair(-1, -1) ;
    }

    return hadronPairV0[0]; 
}

bool AnaV0::shouldKeepPair(const pair<int, int>& a, const pair<int, int>& b) {
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




void AnaV0::SaveMissingMomenta(TVector3 missP)
{
   mRecTree->setPtMissing( missP.Pt() );
   mRecTree->setPxMissing( missP.X() );
   mRecTree->setPyMissing( missP.Y() );   
}

