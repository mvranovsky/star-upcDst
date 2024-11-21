#include "Ana.h"

Ana::Ana(TFile *outFile){
   mOutFile = outFile;
   mUtil = new Util();
}

Ana::~Ana(){
   if(mRecTree) delete mRecTree;
   if(mUtil) delete mUtil;
}

void Ana::SetEvent(StUPCEvent *upcEvt, StRPEvent *rpEvt, StRPEvent *mcEvt){
   mUpcEvt = upcEvt;
   mRpEvt = rpEvt;
   mMcEvt = mcEvt;
}


// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP
//const int Util::triggerID[nTriggers] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
//           570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};



bool Ana::RPInFidRange(double x, double y) const{
   return (abs(y) < 0.8 && abs(y) > 0.4 && x > -0.27 && (x + 0.6)*(x + 0.6) + y*y < 1.25) ? true : false;
}

bool Ana::IsInRpRange(double x, double y, int rpId, TVector3 offSet) const{
   return (abs(y) < abs(offSet[Y]) || abs(y) > abs(offSet[Y]) + 0.03 || x < -0.02 || x > 0.02) ? false : true;
}

bool Ana::CheckTriggers(const vector<int> *triggerArray, StUPCEvent *mUpcEvt, TH1D *hTriggerBits) const{
   if( triggerArray->empty() )
      return true;
   bool triggered = false;
   for(int var = 0; var < nTriggers; ++var)
   {
      if(mUpcEvt->isTrigger(triggerID[var]))
      {  //Checked if it is CPT trigger
         if(hTriggerBits)
            hTriggerBits->Fill(var);
         for (unsigned int i = 0; i < triggerArray->size(); ++i)
            if(triggerID[var] == (*triggerArray)[i])
               triggered=true;
      }
   }
   return triggered;
}


int Ana::hasGoodTPCnSigma(const StUPCTrack *trk){ //zmiernit podmienku na proton, skusit prebehnut aj na K0, Lambda 
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




void Ana::AnaRpTracks(StRPEvent *event)
{
   for (unsigned int i = 0; i < nSides; ++i)
      mRpTrackIdVec_perSide[i].clear();
   for (unsigned int i = 0; i < nBranches; ++i)
      mRpTrackIdVec_perBranch[i].clear();
   for (unsigned int i = 0; i < nRomanPots; ++i)
      mTrackPointIdVec[i].clear();
   for (unsigned int i = 0; i < nRomanPots*nPlanes; ++i)
      mClusterIdVec[i].clear();

   // Loop over all tracks reconstructed in Roman Pots    
   for(unsigned int k = 0; k < event->getNumberOfTracks(); ++k)
   {
      // Get pointer to k-th track in Roman Pot data collection
      StUPCRpsTrack *trk = event->getTrack(k);
      trk->setEvent(event);
      // Get ID of a side in which this k-th track was reconstructed
      int side = trk->branch()<2 ? E : W;

      if( (trk->getTrackPoint(0) ? trk->getTrackPoint(0)->planesUsed()>=nPlanesUsed : false) &&
         (trk->getTrackPoint(1) ? trk->getTrackPoint(1)->planesUsed()>=nPlanesUsed : false))
      {
         mRpTrackIdVec_perSide[side].push_back( k );
         mRpTrackIdVec_perBranch[trk->branch()].push_back( k );
      }
   }

   // Loop over all tracks points reconstructed in Roman Pots
   for(unsigned int k = 0; k < event->getNumberOfTrackPoints(); ++k)
   {
      StUPCRpsTrackPoint *trkPoint = event->getTrackPoint(k);
      if( trkPoint->planesUsed() < nPlanesUsed) // Skip "bad" track points
         continue;
      mTrackPointIdVec[trkPoint->rpId()].push_back( k );
   }

   // Loop over all clusters reconstructed in Roman Pots
   for(unsigned int k = 0; k < event->getNumberOfClusters(); ++k)
   {
      StUPCRpsCluster *cluster = event->getCluster(k);
      mClusterIdVec[cluster->romanPotId()*nPlanes + cluster->planeId()].push_back( k );
   }
}

void Ana::SaveEventInfo(const StUPCEvent *upcEvt)
{
   mRecTree->setTofMult( upcEvt->getTOFMultiplicity());
   mRecTree->setEventNumber( upcEvt->getEventNumber());
   mRecTree->setFillNumber( upcEvt->getFillNumber());
   mRecTree->setBunchCrossId( upcEvt->getBunchCrossId());
   mRecTree->setBunchCrossId7bit( upcEvt->getBunchCrossId7bit());
   mRecTree->setRunNumber( upcEvt->getRunNumber() );
   mRecTree->setNVertecies( upcEvt->getNumberOfVertices() );
}

void Ana::SaveTrackInfo(const StUPCTrack *trk, unsigned int iTrack) // this function uses momentum with trajectory being refitted to primary vertex
{
   mRecTree->setDEdxInKevCm( trk->getDEdxSignal()*1000000 , iTrack); // convert to KeV/cm
   
   //StPicoPhysicalHelix helix(trk->getCurvature() , trk->getDipAngle(), trk->getPhase(), trk->getOrigin(), trk->getCharge() );

   TVector3 momentum;
   trk->getMomentum(momentum);
   mRecTree->setMomentumInGev(  momentum.Mag() , iTrack );
   mRecTree->setPtInGev(  trk->getPt() , iTrack );
   mRecTree->setPxInGev(  momentum.X() , iTrack );
   mRecTree->setPyInGev(  momentum.Y() , iTrack );
   mRecTree->setPzInGev(  momentum.Z() , iTrack );
   mRecTree->setTofTimeInNs(  trk->getTofTime() , iTrack );
   mRecTree->setTofLengthInCm(  trk->getTofPathLength() , iTrack );
   mRecTree->setCharge(  trk->getCharge() , iTrack );
   mRecTree->setDcaXYInCm(  trk->getDcaXY() , iTrack );
   mRecTree->setDcaZInCm(  trk->getDcaZ() , iTrack );
   mRecTree->setNHitsFit(  trk->getNhitsFit() , iTrack );
   mRecTree->setNHitsDEdx(  trk->getNhitsDEdx() , iTrack );
   mRecTree->setEta(  trk->getEta() , iTrack );
   mRecTree->setPhiAngle(  trk->getPhi() , iTrack );
   mRecTree->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kPion) , iTrack, PION );
   mRecTree->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kKaon) , iTrack, KAON );
   mRecTree->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kProton) , iTrack, PROTON );
   if( trk->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk) ){
      mRecTree->setTofHit(1, iTrack);
   } else {
      mRecTree->setTofHit(-1, iTrack);
   }
}


void Ana::SaveTrackInfo(const StUPCTrack *trk, TLorentzVector hadron ,unsigned int iTrack) // this function uses momentum from the picoPhysical helix, not refitted with primary vertex
{
   mRecTree->setDEdxInKevCm( trk->getDEdxSignal()*1000000 , iTrack); // convert to KeV/cm
   
   //StPicoPhysicalHelix helix(trk->getCurvature() , trk->getDipAngle(), trk->getPhase(), trk->getOrigin(), trk->getCharge() );

   TVector3 momentum = hadron.Vect();
   //trk->getMomentum(momentum);
   mRecTree->setMomentumInGev(  momentum.Mag() , iTrack );
   mRecTree->setPtInGev(  sqrt( pow(momentum.X(),2) + pow(momentum.Y(),2) ) , iTrack );
   mRecTree->setPxInGev(  momentum.X() , iTrack );
   mRecTree->setPyInGev(  momentum.Y() , iTrack );
   mRecTree->setPzInGev(  momentum.Z() , iTrack );
   mRecTree->setTofTimeInNs(  trk->getTofTime() , iTrack );
   mRecTree->setTofLengthInCm(  trk->getTofPathLength() , iTrack );
   mRecTree->setCharge(  trk->getCharge() , iTrack );
   mRecTree->setDcaXYInCm(  trk->getDcaXY() , iTrack );
   mRecTree->setDcaZInCm(  trk->getDcaZ() , iTrack );
   mRecTree->setNHitsFit(  trk->getNhitsFit() , iTrack );
   mRecTree->setNHitsDEdx(  trk->getNhitsDEdx() , iTrack );
   mRecTree->setEta(  trk->getEta() , iTrack );
   mRecTree->setPhiAngle(  trk->getPhi() , iTrack );
   mRecTree->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kPion) , iTrack, PION );
   mRecTree->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kKaon) , iTrack, KAON );
   mRecTree->setNSigmaTPC( trk->getNSigmasTPC(StUPCTrack::kProton) , iTrack, PROTON );
   if( trk->getFlag( StUPCTrack::kTof ) && IsGoodTofTrack(trk) ){
      mRecTree->setTofHit(1, iTrack);
   } else {
      mRecTree->setTofHit(-1, iTrack);
   }
}




void Ana::SaveRPinfo(const StUPCRpsTrack *trackRP, unsigned int iSide)
{
   const StUPCRpsTrackPoint* trackPoint = trackRP->getTrackPoint(RP1);

   mRecTree->setRpX( trackPoint->x() , iSide );
   mRecTree->setRpZ( trackPoint->y() , iSide );
   mRecTree->setRpY( trackPoint->z() , iSide );
   mRecTree->setThetaRp( trackRP->thetaRp(2) , iSide );
   mRecTree->setPhiRp( trackRP->phiRp() , iSide );
   mRecTree->setTimeRp( trackRP->time() , iSide );
   mRecTree->setT( trackRP->t(mUtil->p0()) , iSide );
   mRecTree->setPRp( trackRP->p() , iSide );
   mRecTree->setPtRp( trackRP->pt() , iSide );
   mRecTree->setEtaRp( trackRP->eta() , iSide );
   mRecTree->setXi( trackRP->xi(mUtil->p0()) , iSide );
   mRecTree->setPx( trackRP->pVec().X() , iSide );
   mRecTree->setPy( trackRP->pVec().Y() , iSide );
   mRecTree->setPz( trackRP->pVec().Z() , iSide );

/*
   // subtract offset
   int rnNumber, day, run, rpID;
   rpID = trackPoint->rpId();
   rnNumber = upcEvt->getRunNumber() - 18000000;
   day = rnNumber / 1000 - 55;
   run = rnNumber % 1000;
   rpYinner[iSide] = rpY[iSide] - mOffSet[rpID][day][run]; 
*/
}




void Ana::SaveStateInfo(TLorentzVector state, int totQ, unsigned int iState){

   mRecTree->setInvMass( state.M(), iState );
   mRecTree->setTheta( state.Theta(), iState );
   mRecTree->setPhi( state.Phi(), iState );
   mRecTree->setP( state.P(), iState );
   mRecTree->setPt( state.Pt(), iState );
   mRecTree->setRap( state.Rapidity(), iState);
   mRecTree->setTotQ( totQ , iState);
}

void Ana::SaveVertexInfo(const StUPCV0* V0, unsigned int iVtx)
{
   mRecTree->setVertexZInCm( V0->prodVertexHypo().Z(), iVtx );
   mRecTree->setVertexYInCm( V0->prodVertexHypo().Y(), iVtx );
   mRecTree->setVertexXInCm( V0->prodVertexHypo().X(), iVtx );

   mRecTree->setDcaDaughters( V0->dcaDaughters(), iVtx );
   mRecTree->setDcaBeamline( V0->DCABeamLine(), iVtx );
   mRecTree->setPointingAngle( V0->pointingAngleHypo(), iVtx );
   mRecTree->setDecayLength( V0->decayLengthHypo(), iVtx );
   mRecTree->setVertexDiff( -999, iVtx );
}


void Ana::SaveZdcInfo(const StUPCEvent *upcEvt)
{

   for (unsigned int iPmt = 0; iPmt < 3; ++iPmt)
   {
      mRecTree->setZdcAdcEastPmt( upcEvt->getZDCEastADC(iPmt+1), iPmt);
      mRecTree->setZdcAdcWestPmt( upcEvt->getZDCWestADC(iPmt+1), iPmt);
   }
   mRecTree->setZdcTdcEast( upcEvt->getZDCEastTDC() );
   mRecTree->setZdcTdcWest( upcEvt->getZDCWestTDC() );
   mRecTree->setZdcTimeDiff( upcEvt->getZDCTimeDiff() );
   mRecTree->setZdcVertexZ( upcEvt->getZdcVertexZ() );
   mRecTree->setZdcEastRate( upcEvt->getZDCEastRate() );
   mRecTree->setZdcWestRate( upcEvt->getZDCWestRate() );
   mRecTree->setZdcUnAttEast( upcEvt->getZDCUnAttEast() );
   mRecTree->setZdcUnAttWest( upcEvt->getZDCUnAttWest() );

}

void Ana::SaveBbcInfo(const StUPCEvent *upcEvt)
{
   mRecTree->setBbcSmallEast( upcEvt->getBBCSmallEast());
   mRecTree->setBbcSmallWest( upcEvt->getBBCSmallWest());
   mRecTree->setBbcLargeEast( upcEvt->getBBCLargeEast());
   mRecTree->setBbcLargeWest( upcEvt->getBBCLargeWest());
}


void Ana::SaveTriggerInfo(const StUPCEvent *upcEvt, const StRPEvent *rpEvt)
{
   SaveBbcInfo(upcEvt);
   SaveZdcInfo(upcEvt);

   mRecTree->setBbcTrigBit( upcEvt->getBBCSmallEast() > 20 || upcEvt->getBBCSmallWest() > 20 || upcEvt->getBBCLargeEast() > 50 || upcEvt->getBBCLargeWest() > 50);

   mRecTree->setZdcETrigBit(false);
   mRecTree->setZdcWTrigBit(false);
   if( (upcEvt->getZDCEastTDC() > 500 && upcEvt->getZDCEastTDC() < 2700) && (upcEvt->getZDCEastADC(1) > 25 || upcEvt->getZDCEastADC(2) > 25 || upcEvt->getZDCEastADC(3) > 25))
   {
      mRecTree->setZdcTrigBit(true);      
      mRecTree->setZdcETrigBit(true);
   }   
   if( (upcEvt->getZDCWestTDC() > 500 && upcEvt->getZDCWestTDC() < 2700) && (upcEvt->getZDCWestADC(1) > 25 || upcEvt->getZDCWestADC(2) > 25 || upcEvt->getZDCWestADC(3) > 25))
   {
      mRecTree->setZdcTrigBit(true);
      mRecTree->setZdcWTrigBit(true);
   }    


   mRecTree->setTofTrigBit( !( upcEvt->getTOFMultiplicity() < 2 || upcEvt->getTOFMultiplicity() > 10) );
   saveRpTrigBit(rpEvt);

   bool RpEt = upcEvt->getLastDSM0() & (1 << 2);
   bool TofMult0 = upcEvt->getLastDSM0() & (1 << 4);
   bool TofMult2 = upcEvt->getLastDSM0() & (1 << 6);
   bool RpIt = upcEvt->getLastDSM0() & (1 << 12);
   bool BbcE = upcEvt->getLastDSM1() & (1 << 1);
   bool BbcW = upcEvt->getLastDSM1() & (1 << 2);
   bool BbcLE = upcEvt->getLastDSM1() & (1 << 3);
   bool BbcLW = upcEvt->getLastDSM1() & (1 << 4);
   bool ZdcE = upcEvt->getLastDSM1() & (1 << 7);
   bool ZdcW = upcEvt->getLastDSM1() & (1 << 10);

   mRecTree->setRpItDsmBit(RpIt);
   mRecTree->setRpEtDsmBit(RpEt);
   mRecTree->setTofDsmABit( TofMult0 );
   mRecTree->setTofDsmBBit( TofMult2 );
   mRecTree->setTofDsmBit( ( TofMult0 && !TofMult2 ) );
   mRecTree->setRpDsmBit( ( (RpEt && !RpIt) || (!RpEt && RpIt) ) );
   mRecTree->setBbcSmallEDsmBit(BbcE);
   mRecTree->setBbcSmallWDsmBit(BbcW);
   mRecTree->setBbcLargeEDsmBit(BbcLE);
   mRecTree->setBbcLargeWDsmBit(BbcLW);
   mRecTree->setBbcDsmBit( ( BbcE || BbcW || BbcLE || BbcLW ) );
   mRecTree->setZdcEDsmBit(ZdcE);
   mRecTree->setZdcWDsmBit(ZdcW);
   mRecTree->setZdcDsmBit( ( ZdcE || ZdcW ) );

}

void Ana::saveRpTrigBit(const StRPEvent *rpEvt)
{
   for (unsigned int iRp = 0; iRp < nRomanPots; ++iRp)
   {
      mRecTree->setRpTrigBits( IsRpTrigBit(rpEvt, iRp), iRp);
   }

   bool EU = mRecTree->getRpTrigBits(E1U) || mRecTree->getRpTrigBits(E2U);
   bool ED = mRecTree->getRpTrigBits(E1D) || mRecTree->getRpTrigBits(E2D);
   bool WU = mRecTree->getRpTrigBits(W1U) || mRecTree->getRpTrigBits(W2U);
   bool WD = mRecTree->getRpTrigBits(W1D) || mRecTree->getRpTrigBits(W2D);

   bool RP_ET = (EU && WD) || (ED && WU);
   bool RP_IT = (EU && WU) || (ED && WD);

   mRecTree->setRpItTrigBit(RP_IT);
   mRecTree->setRpEtTrigBit(RP_ET);
   mRecTree->setRpTrigBit((RP_ET && !RP_IT) || (!RP_ET && RP_IT));

}

bool Ana::IsRpTrigBit(const StRPEvent *rpEvt, unsigned int iRp)
{
   for (unsigned int iPmt = 0; iPmt < 2; ++iPmt)
      if( rpEvt->tac(iRp, iPmt) > 200 && rpEvt->tac(iRp, iPmt) < 1750 && rpEvt->adc(iRp, iPmt) > 30)
         return true;

   return false;   
}
/*
void Ana::CreateTofClusters()
{
   vector<TofCluster> tofClusterVec;
   int nTofClusters = 0;

   //CLUSTERING BASED ON ETA & PHI OF TOF ELEMENT
   for(unsigned int k=0; k< mUpcEvt->getNumberOfHits(); ++k)
   {
      const StUPCTofHit* hit = &mUpcEvt->getHit(k);
      int tray_A = hit->getTray();
      int module_A = hit->getModule();
      int cell_A = hit->getCell();  

      // skipping vpd hits
      if( (module_A==0 && tray_A==121) || (module_A==0 && tray_A==122) ) 
         continue;

      double phi_A = mTofPhi[tray_A-1][module_A-1][cell_A-1];
      double eta_A = mTofEta[tray_A-1][module_A-1][cell_A-1];
      double x_A = mTofX[tray_A-1][module_A-1][cell_A-1];
      double y_A = mTofY[tray_A-1][module_A-1][cell_A-1];
      double z_A = mTofZ[tray_A-1][module_A-1][cell_A-1];

      TVector2 vec_A( cos(phi_A), sin(phi_A) );

      bool attachToExistingCluster = false;
      int existingClusterIndex = -1;

      for(unsigned int i=0; i<tofClustersVec.size(); ++i)
      {
         for(unsigned int j=0; j<tofClustersVec[i].mTofHitIndexVec.size(); ++j)
         {
            const StUPCTofHit* hit_B = &mUpcEvt->getHit( tofClustersVec[i].mTofHitIndexVec[j] );
            int tray_B = hit_B->getTray();
            int module_B = hit_B->getModule();
            int cell_B = hit_B->getCell();  

            // skipping vpd hits
            if( (module_B==0 && tray_B==121) || (module_B==0 && tray_B==122) ) 
               continue;

            double phi_B = mTofPhi[tray_B-1][module_B-1][cell_B-1];
            double eta_B = mTofEta[tray_B-1][module_B-1][cell_B-1];

            TVector2 vec_B( cos(phi_B), sin(phi_B) );

            double distanceInPhi = vec_A.DeltaPhi( vec_B );
            double distanceInEta = fabs( eta_A - eta_B );

            if( sqrt(distanceInPhi*distanceInPhi + distanceInEta*distanceInEta) < 0.1)
            {
               attachToExistingCluster = true;
               existingClusterIndex = i;
            }
         }
      }

      if( !attachToExistingCluster )
      {
         tofClustersVec.push_back( TofCluster() );
         existingClusterIndex = nTofClusters;
         ++nTofClusters;
      }

      tofClustersVec[existingClusterIndex].mClusterSize++;
      tofClustersVec[existingClusterIndex].mTofHitIndexVec.push_back( k );
      tofClustersVec[existingClusterIndex].mHitXVec.push_back( x_A );
      tofClustersVec[existingClusterIndex].mHitYVec.push_back( y_A );
      tofClustersVec[existingClusterIndex].mHitZVec.push_back( z_A );
   }
}
*/

void Ana::resetInfo() {

   //tracks
   for (int iTrack = 0; iTrack < nHadrons; ++iTrack) {
      mRecTree->setDEdxInKevCm( -9999 , iTrack); // convert to KeV/cm
      mRecTree->setMomentumInGev(  -9999 , iTrack );
      mRecTree->setPtInGev(  -9999, iTrack );
      mRecTree->setPxInGev(  -9999 , iTrack );
      mRecTree->setPyInGev(  -9999 , iTrack );
      mRecTree->setPzInGev(  -9999 , iTrack );
      mRecTree->setTofTimeInNs(  -9999 , iTrack );
      mRecTree->setTofLengthInCm(  -9999 , iTrack );
      mRecTree->setCharge(  -9999 , iTrack );
      mRecTree->setDcaXYInCm(  -9999 , iTrack );
      mRecTree->setDcaZInCm(  -9999 , iTrack );
      mRecTree->setNHitsFit(  -9999 , iTrack );
      mRecTree->setNHitsDEdx(  -9999 , iTrack );
      mRecTree->setEta(  -9999 , iTrack );
      mRecTree->setPhiAngle(  -9999 , iTrack );
      mRecTree->setNSigmaTPC( -9999 , iTrack, PION );
      mRecTree->setNSigmaTPC( -9999 , iTrack, KAON );
      mRecTree->setNSigmaTPC( -9999 , iTrack, PROTON );
   }

   //states
   for (int iState = 0; iState < nStates; ++iState){
      mRecTree->setInvMass( -9999, iState );
      mRecTree->setTheta( -9999, iState );
      mRecTree->setPhi( -9999, iState );
      mRecTree->setP( -9999, iState );
      mRecTree->setPt( -9999, iState );
      mRecTree->setRap( -9999, iState);
      mRecTree->setTotQ( -9999 , iState);
   }

   //vertices
   for (int iVtx = 0; iVtx < nStates; ++iVtx){
      mRecTree->setVertexZInCm( -9999, iVtx );
      mRecTree->setVertexYInCm( -9999, iVtx );
      mRecTree->setVertexXInCm( -9999, iVtx );
      mRecTree->setDcaDaughters( -9999, iVtx );
      mRecTree->setDcaBeamline( -9999, iVtx );
      mRecTree->setPointingAngle( -9999, iVtx );
      mRecTree->setDecayLength( -9999, iVtx );
      mRecTree->setVertexDiff( -9999, iVtx ); 
   }


}
