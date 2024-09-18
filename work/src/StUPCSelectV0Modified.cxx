
//C++
#include <vector>
#include <iostream>

//local classes
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCVertex.h"
#include "StUPCSelectV0Modified.h"
#include "StUPCV0.h"
#include "Util.h"

using namespace std;

//_____________________________________________________________________________
vector<int> StUPCSelectV0Modified::selectTracks(vector<int> tracks, StUPCEvent *upcEvent, int nVertices) {

  //run over global tracks and store all global tracks of interest.
  double massPion = 0.13957061;
  //double massKaon =  0.497611;
  double massProton = 0.93827;  
  vector<int> output, currentV0;
  double bField = upcEvent->getMagneticField();
  double beamline[4]; 
  beamline[0] = upcEvent->getBeamXPosition();
  beamline[2] = upcEvent->getBeamXSlope();
  beamline[1] = upcEvent->getBeamYPosition();
  beamline[3] = upcEvent->getBeamYSlope();

  TVector3 vertex;
  //loop of all vertices, because of difference in defining vertices between PicoDst and UpcDst 
  for (int iVtx = 0; iVtx < nVertices; ++iVtx){ //double counting V0
    vertex.SetXYZ(upcEvent->getVertex(iVtx)->getPosX(), upcEvent->getVertex(iVtx)->getPosY(), upcEvent->getVertex(iVtx)->getPosZ() );
    //outer tracks loop 
    for(size_t itrk=0; itrk<tracks.size(); itrk++) {

      const StUPCTrack* track1 = upcEvent->getTrack(tracks[itrk]);

      Bool_t matchTof1 =  track1->getFlag( StUPCTrack::kTof );

      //inner tracks loop 
      for(size_t jtrk=itrk+1; jtrk<tracks.size(); jtrk++) {

        const StUPCTrack* track2 = upcEvent->getTrack(tracks[jtrk]);

        Bool_t matchTof2 =  track2->getFlag( StUPCTrack::kTof );

        //condition for at least one tof match
        if ( matchTof1 == kFALSE && matchTof2 == kFALSE ) 
          continue;
        //removed condition for opposite sign->is added at the end of make()
        //removed conditions for track quality, because they are redundant

        //cout << "pair" << endl;
        //track indices itrk, jtrk are redundant in StUPCV0
        StUPCV0 K0L1(track1,track2, massPion, massPion, itrk, jtrk, vertex, beamline, bField, true);
        TVector3 vertex0(0,0,0);
        StUPCV0 K0L2(track1,track2, massPion, massPion, itrk, jtrk, vertex0, beamline, bField, true);

        if ( K0L1.dcaDaughters() > 100. &&  K0L2.dcaDaughters() > 100.)
          continue; 

        StUPCV0 K0L( K0L1.dcaDaughters() > K0L2.dcaDaughters() ? &K0L2 : &K0L1);

        StUPCV0 K0(track1,track2, massPion, massPion, itrk, jtrk, vertex, beamline, bField, false);
        //pridat StUPCV0 na lambdy

        // the same cut for K0 and Lambda
        if ( !(K0.dcaDaughters() < 3.0 && K0.DCABeamLine() < 2.5 && (K0.pointingAngleHypo()>0.9 || K0.decayLengthHypo()<3.0) ) )
          continue;


        //counter for found V0 candidates
        //currentV0 = {tracks[itrk], tracks[jtrk]};
        output.push_back(tracks[itrk]);  //save the pair of tracks to output
        output.push_back(tracks[jtrk]);

        //do outputu pridat vertexid
      }//inner tracks loop
    }//outer tracks loop
  }//vertex loop

  return output;

}//selectTracks

/*

vector<int> StUPCSelectV0Modified::selectTracksK0(vector<int> tracks, StUPCEvent *upcEvent, TVector3 const & vertex) {

  //run over global tracks and store all global tracks of interest.
  //double massPion = 0.13957061;
  //double massKaon =  0.497611;
  //double massProton = 0.93827;  
  vector<int> output, currentV0;
  double bField = upcEvent->getMagneticField();
  double beamline[4]; 
  beamline[0] = upcEvent->getBeamXPosition();
  beamline[2] = upcEvent->getBeamXSlope();
  beamline[1] = upcEvent->getBeamYPosition();
  beamline[3] = upcEvent->getBeamYSlope();

  //outer tracks loop of pions
  for(size_t itrk=0; itrk<tracks.size(); itrk++) {

    const StUPCTrack* track1 = upcEvent->getTrack(tracks[itrk]);

    Bool_t matchTof1 =  track1->getFlag( StUPCTrack::kTof );

    //inner tracks loop of pions
    for(size_t jtrk=itrk+1; jtrk<tracks.size(); jtrk++) {

      const StUPCTrack* track2 = upcEvent->getTrack(tracks[jtrk]);

      Bool_t matchTof2 =  track2->getFlag( StUPCTrack::kTof );

      //condition for at least one tof match
      if ( matchTof1 == kFALSE && matchTof2 == kFALSE ) 
        continue;
      //removed condition for opposite sign->is added at the end of make()
      //removed conditions for track quality, because they are redundant

      //cout << "pair" << endl;
      //track indices itrk, jtrk are redundant in StUPCV0
      StUPCV0 K0L1(track1,track2, mass(PION), mass(PION), itrk, jtrk, vertex, beamline, bField, true);
      TVector3 vertex0(0,0,0);
      StUPCV0 K0L2(track1,track2, mass(PION), mass(PION), itrk, jtrk, vertex0, beamline, bField, true);

      if ( K0L1.dcaDaughters() > 100. &&  K0L2.dcaDaughters() > 100.)
        continue; 

      StUPCV0 K0L( K0L1.dcaDaughters() > K0L2.dcaDaughters() ? &K0L2 : &K0L1);

      StUPCV0 K0(track1,track2, mass(PION), mass(PION), itrk, jtrk, vertex, beamline, bField, false);

      // the same cut for K0 and Lambda
      if ( !(K0.dcaDaughters() < 3.0 && K0.DCABeamLine() < 2.5 && (K0.pointingAngleHypo()>0.9 || K0.decayLengthHypo()<3.0) ) )
        continue;


      //counter for found V0 candidates
      currentV0 = {tracks[itrk], tracks[jtrk]};
      output.push_back(currentV0);  //save the pair of tracks to output
      
    }//inner tracks loop
  }//outer tracks loop

  return output;

}//selectTracks

vector<int> StUPCSelectV0Modified::selectTracksLambda(vector<int> pionTracks, vector<int> protonTracks, StUPCEvent *upcEvent, TVector3 const & vertex) {

  //run over global tracks and store all global tracks of interest.
  //double massPion = 0.13957061;
  //double massKaon =  0.497611;
  //double massProton = 0.93827;  
  vector<int> currentV0, output;
  double bField = upcEvent->getMagneticField();
  double beamline[4]; 
  beamline[0] = upcEvent->getBeamXPosition();
  beamline[2] = upcEvent->getBeamXSlope();
  beamline[1] = upcEvent->getBeamYPosition();
  beamline[3] = upcEvent->getBeamYSlope();

  //outer tracks loop
  for(size_t itrk=0; itrk<pionTracks.size(); itrk++) {

    const StUPCTrack* track1 = upcEvent->getTrack(pionTracks[itrk]);
    Bool_t matchTof1 =  track1->getFlag( StUPCTrack::kTof );

    //inner tracks loop
    for(size_t jtrk=0; jtrk<protonTracks.size(); jtrk++) {

      const StUPCTrack* track2 = upcEvent->getTrack(protonTracks[jtrk]);
      Bool_t matchTof2 =  track2.getFlag( StUPCTrack::kTof );

      if ( matchTof1 == kFALSE && matchTof2 == kFALSE ) 
        continue;

      //cout << "pair" << endl;

      //naco tu je KOL1 a KOL2 a neskor KOL? 
      StUPCV0 K0L1(track1,track2, mass(PION), mass(PROTON), itrk, jtrk, vertex, beamline, bField, true);
      TVector3 vertex0(0,0,0);
      StUPCV0 K0L2(track1,track2, mass(PION), mass(PROTON), itrk, jtrk, vertex0, beamline, bField, true);

      if ( K0L1.dcaDaughters() > 100. &&  K0L2.dcaDaughters() > 100.)
        continue; 

      StUPCV0 K0L( K0L1.dcaDaughters() > K0L2.dcaDaughters() ? &K0L2 : &K0L1);

      StUPCV0 K0(track1,track2, mass(PION), mass(PROTON), itrk, jtrk, vertex, beamline, bField, false);

      // the same cut for K0 and Lambda
      if ( !(K0.dcaDaughters() < 3.0 && K0.DCABeamLine() < 2.5 && (K0.pointingAngleHypo()>0.9 || K0.decayLengthHypo()<3.0) ) )
        continue;


      //counter for found V0 candidates
      currentV0 = {pionTracks[itrk], protonTracks[jtrk]};
      output.push_back(currentV0);

    }//inner tracks loop
  }//outer tracks loop

  return output;

}//selectTracks
*/