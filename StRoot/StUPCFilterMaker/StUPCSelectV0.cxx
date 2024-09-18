
//C++
#include <vector>
#include <iostream>

//local classes
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCSelectV0.h"
#include "StUPCV0.h"


using namespace std;

//_____________________________________________________________________________
int StUPCSelectV0::selectTracks(vector<StUPCTrack>& tracks, std::vector<UChar_t>& sel, StUPCEvent *upcEvent, TVector3 const & vertex, TList*& mHistList) {

  //run over global tracks and store all global tracks of interest.
  double massPion = 0.13957061;
  //double massKaon =  0.497611;
  double massProton = 0.93827;  

  loadHists(mHistList);

  int nfound = 0;
  double bField = upcEvent->getMagneticField();
  double beamline[4]; 
  beamline[0] = upcEvent->getBeamXPosition();
  beamline[2] = upcEvent->getBeamXSlope();
  beamline[1] = upcEvent->getBeamYPosition();
  beamline[3] = upcEvent->getBeamYSlope();

  //outer tracks loop
  for(size_t itrk=0; itrk<tracks.size(); itrk++) {

    const StUPCTrack& track1 = tracks[itrk];

    if( track1.getNhits() < 15 ) continue;
    if( track1.getPt() < 0.15 ) continue;
    hEta->Fill( track1.getEta() );
    if( abs(track1.getEta()) > 1.1 ) continue;

    Bool_t matchTof1 =  track1.getFlag( StUPCTrack::kTof );

    //inner tracks loop
    for(size_t jtrk=itrk+1; jtrk<tracks.size(); jtrk++) {

      const StUPCTrack& track2 = tracks[jtrk];

      if( track2.getPt() < 0.15 ) continue;
      if( track2.getEta() > 1.1 ) continue;
      

      Bool_t matchTof2 =  track2.getFlag( StUPCTrack::kTof );

      if ( matchTof1 == kFALSE && matchTof2 == kFALSE ) 
        continue;
      if ( track1.getCharge()+track2.getCharge() != 0 ) 
        continue;
      if ( track2.getNhits()<15 ) 
        continue;

      //cout << "pair" << endl;

      StUPCV0 K0L1(&track1,&track2, massPion, massPion, itrk, jtrk, vertex, beamline, bField, true);
      TVector3 vertex0(0,0,0);
      StUPCV0 K0L2(&track1,&track2, massPion, massPion, itrk, jtrk, vertex0, beamline, bField, true);

      if ( K0L1.dcaDaughters() > 100. &&  K0L2.dcaDaughters() > 100.)
        continue; 

      StUPCV0 K0L( K0L1.dcaDaughters() > K0L2.dcaDaughters() ? &K0L2 : &K0L1);

      bool V0Candidate = false;
      StUPCV0 K0(&track1,&track2, massPion, massPion, itrk, jtrk, vertex, beamline, bField, false);

      hDcaDaughters->Fill(K0.dcaDaughters());
      hDcaBeamline->Fill(K0.DCABeamLine());
      hPointingAngle->Fill(K0.pointingAngleHypo());
      hDecayLength->Fill(K0.decayLengthHypo());
      hDecayLPointingA->Fill(K0.decayLengthHypo(),K0.pointingAngleHypo() );
      hDcaBeamDaughters->Fill(K0.DCABeamLine(), K0.dcaDaughters());

      // the same cut for K0 and Lambda
      if ( !(K0.dcaDaughters() < 3.0 && K0.DCABeamLine() < 2.5 && (K0.pointingAngleHypo()>0.9 || K0.decayLengthHypo()<3.0) ) )
        continue;

      hEtaAfter->Fill(track1.getEta());
      hEtaAfter->Fill(track2.getEta());

      if ( abs(K0.m()-0.495) < 0.05 ) // it's a K0
        V0Candidate = true;

      // check if pair is Lambda/LambdaBar
      StUPCV0 L01(&track1,&track2, massPion, massProton, itrk, jtrk, vertex, beamline, bField, false);
      if ( abs(L01.m()-1.115) < 0.02 ) // it's a Lambda/LambdaBar
        V0Candidate = true;

      // check if pair is LambdaBar/Lambda
      StUPCV0 L02(&track1,&track2, massProton, massPion, itrk, jtrk, vertex, beamline, bField, false);
      if ( abs(L02.m()-1.115) < 0.02 ) // it's a LambdaBar/Lambda
        V0Candidate = true;

      if ( !V0Candidate ) // it is not a V0 canidate, discard the pair
        continue;

      //counter for found V0 candidates
      nfound++;

      //mark the tracks as used
      sel[itrk] |= (1 << StUPCTrack::kV0);
      sel[jtrk] |= (1 << StUPCTrack::kV0);
    }//inner tracks loop
  }//outer tracks loop

  return nfound;

}//selectTracks
//__________________________________________________________________________________________________________
bool StUPCSelectV0::loadHists(TList *& mHistList){

  hDcaDaughters = (TH1D*)mHistList->FindObject("hDcaDaughters");
  hDcaBeamline = (TH1D*)mHistList->FindObject("hDcaBeamline");
  hPointingAngle = (TH1D*)mHistList->FindObject("hPointingAngle");
  hDecayLength = (TH1D*)mHistList->FindObject("hDecayLength");
  hEta = (TH1D*)mHistList->FindObject("hEta");
  hEtaAfter = (TH1D*)mHistList->FindObject("hEtaAfter");


  if(!hDcaBeamline || !hDcaDaughters || !hDecayLength || !hPointingAngle || !hEta || !hEtaAfter){
    cout << "Could not load one of histograms from mHistList to StUPCSelectV0. Leaving..." << endl;
    return false;
  }

  hDecayLPointingA = (TH2D*)mHistList->FindObject("hDecayLPointingA");
  hDcaBeamDaughters = (TH2D*)mHistList->FindObject("hDcaBeamDaughters");

  if(!hDecayLPointingA || !hDcaBeamDaughters){
    cout << "Could not load one of 2D histograms from mHistList to StUPCSelectV0. Leaving..." << endl;
    return false;
  }
  return true;

}
