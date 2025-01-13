#include "TofEffMult.h"
#include "RunDef.h"

//_____________________________________________________________________________
TofEffMult::TofEffMult(TFile *outFile): Ana(outFile){}

TofEffMult::~TofEffMult(){
   //if(mUtil) delete mUtil;
}


void TofEffMult::Make(){

   hAnalysisFlow->Fill(TOFALL);

   if(!runMCAna && !CheckTriggers(&CEPtriggers, mUpcEvt, hTriggerBits))
      return;

   hAnalysisFlow->Fill(TOFTRIG);

   for (int itrk = 0; itrk < mUpcEvt->getNumberOfTracks(); ++itrk){

      const StUPCTrack* trk = mUpcEvt->getTrack(itrk);

      // conditions on which data we are using- old or new
      if( trk->getFlag( StUPCTrack::kV0 ) && V0Control)
         continue;
      if( trk->getFlag( StUPCTrack::kPrimary ) && !V0Control)
         continue;


      fillTrackQualityCuts(trk);
      if(!IsGoodTrack(trk))
         continue;

      //particle identification: pion
      fillNSigmaPlots(trk);
      if( hasGoodTPCnSigma(trk) != PION )
         continue;

      hadronID.push_back(itrk);

      // ToF hit condition
      if( !( (trk->getFlag( StUPCTrack::kTof) && IsGoodTofTrack(trk) ) || ( trk->getFlag(StUPCTrack::kTof) && runMCAna ) ) )
         continue;

      tagID.push_back(itrk);


   }//global tracks loop
   hNTracksTof->Fill( tagID.size() );
   hNTracksTpc->Fill( hadronID.size() );

   //cout << "hadronID: " << hadronID.size() << endl;
   //cout << "tagID: " << tagID.size() << endl;


   if(tagID.size() == 0 || hadronID.size() == 0 ){
      resetInfo();
      return;
   }

   hAnalysisFlow->Fill(TOFTRACKQUALITY);


   mRecTree->setNGoodTpcTrks( hadronID.size() );
   SaveEventInfo(mUpcEvt);

   if(!runMCAna){
     SaveTriggerInfo(mUpcEvt, mRpEvt);
   }
   //info about the beamline and mag field
   fillBeamlineInfo();

   int idx = 0;
   for (int iTrack = 0; iTrack < (unsigned)tagID.size(); ++iTrack){   //outer track loop of tags
      const StUPCTrack* trk1 = mUpcEvt->getTrack(tagID[iTrack]);
      

      //inner track loop of probes
      for (int jTrack = 0; jTrack < (unsigned)hadronID.size(); ++jTrack){

         if(idx > 4)
            break;

         if(tagID[iTrack] == hadronID[jTrack])
            continue;


         const StUPCTrack* trk2 = mUpcEvt->getTrack(hadronID[jTrack]);


         StUPCV0* V0 = getV0(tagID[iTrack] ,trk1, hadronID[jTrack], trk2);

         if( V0 == nullptr || !V0->isInitialized() )
            continue;

         fillEtaVtxPlotsBefore(trk1, trk2, V0->prodVertexHypo().Z() );
         
         if(abs(V0->prodVertexHypo().Z()) > vertexRangeForEVz  )
            continue;

         SaveVertexInfo(V0, idx);

         //if(abs(V0->prodVertexHypo().Z()) > vertexRange )
         //   continue;


         // special eta-vertexZ cut
         if( !( IsGoodEtaTrack(trk1, idx) && IsGoodEtaTrack(trk2, idx ) ) )
            continue;

         hAnalysisFlow->Fill( TOFETAVTXZ );

         fillEtaVtxPlotsAfter(trk1, trk2, V0->prodVertexHypo().Z() );

         //same cuts as were in StUPCSelectV0 except K0L1, K0L2 -> K0
         fillTopologyCutsBefore(*V0);
         if ( !(V0->dcaDaughters() < maxDcaDaughters && V0->DCABeamLine() < maxDcaBeamLine && ( V0->pointingAngleHypo()> minPointingAngle || V0->decayLengthHypo()< maxDecayLengthHypo) ) )
             continue;
         fillTopologyCutsAfter(*V0);

         hAnalysisFlow->Fill(TOFPAIR);


         hTotQ->Fill( trk1->getCharge() + trk2->getCharge() );
         
         if(trk1->getCharge() + trk2->getCharge() == 0){
            if(trk2->getFlag( StUPCTrack::kTof ) && (IsGoodTofTrack(trk2) || runMCAna ) ){
               hInvMassTof2->Fill( V0->lorentzVector().M() );
               hInvMassTof1->Fill( V0->lorentzVector().M() );      
            }
            else{
               hInvMassTof1->Fill( V0->lorentzVector().M() );
            }
         } else continue;

         hAnalysisFlow->Fill(TOFOPPOSITE);

         TLorentzVector hadron1, hadron2, state;
         trk1->getLorentzVector(hadron1, mUtil->mass(PION));
         trk2->getLorentzVector(hadron2, mUtil->mass(PION));
         
         state = V0->lorentzVector();
         if(state.M() < 0.45 || state.M() > 0.55)
            continue;


         SaveTrackInfo(trk1,hadron1 ,2*idx);
         SaveTrackInfo(trk2,hadron2 ,2*idx + 1);

         SaveStateInfo(V0->lorentzVector(), 0 , idx);
         idx += 1;

      }//inner loop
      if(idx > 4)
         break;
   }//outer loop

   hV0perEvent->Fill(idx);

   // checking for 0 states
   if(idx == 0){
      resetInfo();
      return;
   }

    // save event info
    mRecTree->FillRecTree();

    // end with clearing all the variables for rec tree 
    resetInfo();


   if(DEBUG){
      cout << "Finished TofEffMult::Make()" << endl;
    }
}


StUPCV0* TofEffMult::getV0(int itrk1,const StUPCTrack* trk1,int itrk2, const StUPCTrack* trk2){
 

   // find out if both tracks belong to the same primary vtx
   const StUPCVertex* vtx = trk1->getVertex();
   StUPCV0 *V0;


   if(usePrimVtx && !(vtx && trk1->getVertexId() == trk2->getVertexId() ) ){
      return nullptr;
   } 
   else if(usePrimVtx && (vtx && trk1->getVertexId() == trk2->getVertexId() ) ){
      V0 = new StUPCV0(trk1,trk2, mUtil->mass(PION), mUtil->mass(PION), itrk1, itrk2, vtx->getPosVtx() , beamline, bField, false);
      StUPCV0 V0alt(trk1,trk2, mUtil->mass(PION), mUtil->mass(PION),itrk1, itrk2,{0,0,0}, beamline, bField, false);
      hVtxDiff->Fill( abs(V0->prodVertexHypo().Z() - V0alt.prodVertexHypo().Z() ) );
   }
   else if(!usePrimVtx && (vtx && trk1->getVertexId() == trk2->getVertexId() ) ){
      V0 = new StUPCV0(trk1,trk2, mUtil->mass(PION), mUtil->mass(PION), itrk1, itrk2,{0,0,0} , beamline, bField, false);
      StUPCV0 V0alt(trk1,trk2, mUtil->mass(PION), mUtil->mass(PION),itrk1, itrk2,vtx->getPosVtx(), beamline, bField, false);
      hVtxDiff->Fill( abs(V0->prodVertexHypo().Z() - V0alt.prodVertexHypo().Z() ) );
   }else{
      V0 = new StUPCV0(trk1,trk2, mUtil->mass(PION), mUtil->mass(PION),itrk1, itrk2, {0,0,0}, beamline, bField, false);
   }

   return V0;


}

void TofEffMult::Init(){

   mUtil = new Util();


   if( DEBUG )
     cout<<"TofEffMult::Init() called"<<endl;

   mOutFile->cd();
   hAnalysisFlow = new TH1D("hAnalysisFlow", "CutsFlow", nTOFEFFCuts-1, 1, nTOFEFFCuts);
   for(int tb=1; tb<nTOFEFFCuts; ++tb) {
     hAnalysisFlow->GetXaxis()->SetBinLabel(tb, mUtil->analysisTofEff(tb));
   }

   hTriggerBits = new TH1D("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
   for(int tb=0; tb<nTriggers; ++tb){
     TString label; label.Form("%d",triggerID[tb]);
     hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
   }
   mRecTree = new RecTree(nameOfTofEffMultTree, TofEffMultTreeBits, false); 


   hEta = new TH1D("hEta", "Pseudorapidity; #eta; counts", 60, -2, 2);
   hEtaCut = new TH1D("hEtaCut", "Pseudorapidity; #eta [-]; counts", 60, -2, 2);

   hEtaPhi = new TH2F("hEtaPhi","Phi vs eta of TOF tracks; #eta; #varphi",100,-2,2,100,-4,4);
   hEtaPhiCut = new TH2F("hEtaPhiCut","Phi vs eta of TOF tracks; #eta ; #varphi",100,-2,2,100,-4,4);

   hPosZ =  new TH1D("hPosZ", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 
   hPosZCut =  new TH1D("hPosZCut", "Position of z_{vertex}; Vertex_{Z} [cm]; counts", 60, -150, 150); 

   hEtaVtxZ = new TH2F("hEtaVtxZ", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);
   hEtaVtxZCut = new TH2F("hEtaVtxZCut", "hEtaVtxZ; #eta [-]; V_{Z} [cm]", 40, -1, 1,40 ,-100, 100);


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

   hNTracksTof = new TH1D("hNTracksTof", "hNTracksTof; Number of ToF tracks [-]; counts", 11, -0.5, 10.5);

   hNTracksTpc = new TH1D("hNTracksTpc", "hNTracksTpc; Number of TPC tracks [-]; counts", 26, -0.5, 25.5);

   hMultipleGoodV0 = new TH1D("hMultipleGoodV0", "hMultipleGoodV0; 1 = good V0 not accounted; counts",90 , 0.45, 0.54);
   hInvMassTof1 = new TH1D("hInvMassTof1", "hInvMassTof1; m_{#pi^{+} #pi^{-} }; counts",90 , 0.45, 0.54);
   hInvMassTof2 = new TH1D("hInvMassTof2", "hInvMassTof2; m_{#pi^{+} #pi^{-} }; counts",90 , 0.45, 0.54);

   hV0perEvent = new TH1D("hV0perEvent", "hV0perEvent; V0 per event [-]; counts", 6,-0.5, 5.5);


}

void TofEffMult::fillTopologyCutsBefore(const StUPCV0& V0){

    hDcaDaughters->Fill(V0.dcaDaughters());
    hDcaBeamline->Fill(V0.DCABeamLine());
    hPointingAngle->Fill(V0.pointingAngleHypo());
    hDecayLength->Fill(V0.decayLengthHypo());
    hDecayLPointingA->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());

}


void TofEffMult::fillTopologyCutsAfter(const StUPCV0& V0){

    hDcaDaughtersCut->Fill(V0.dcaDaughters());
    hDcaBeamlineCut->Fill(V0.DCABeamLine());
    hPointingAngleCut->Fill(V0.pointingAngleHypo());
    hDecayLengthCut->Fill(V0.decayLengthHypo());
    hDecayLPointingACut->Fill(V0.decayLengthHypo(), V0.pointingAngleHypo());
}

void TofEffMult::fillTrackQualityCuts(const StUPCTrack* trk){

   hNfitHits->Fill(trk->getNhitsFit() );
   hNhitsDEdx->Fill(trk->getNhitsDEdx() );
   hPt->Fill(trk->getPt());
}


void TofEffMult::fillNSigmaPlots(const StUPCTrack *trk){
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