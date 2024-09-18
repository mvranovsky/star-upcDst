#include "UpcDstLibreries.h"
#include "RunDef.h"
#include "RpMCAna.h"

RpMCAna::RpMCAna(TFile *outFile): Ana(outFile){
   embedMaker = new EmbedMaker();
   embedMaker->Init(outFile);
   embedMaker->setAfterburner(AFTERBURNER);
}

RpMCAna::~RpMCAna(){

   for (int iSet = 0; iSet < DATA; ++iSet)
      if(mElAna[iSet]) delete mElAna[iSet];

   if(embedMaker) delete embedMaker;
}

void RpMCAna::Make()
{
   // check RP trigger
   for (int iSet = 0; iSet < DATA; ++iSet)
   {
      runEmbedding(iSet);
      mElAna[iSet]->SetEvent(mUpcEvt, mEmbedEvt, mMcEvt);
      mElAna[iSet]->SetRPMCInfo(mVertex,mMomentum);
      mElAna[iSet]->Make();
      mElAna[iSet]->saveRpTrigBit(mEmbedEvt);
      runMCEff(iSet);
      if(mEmbedEvt) delete mEmbedEvt;
   }
}

void RpMCAna::Init()
{
   if( DEBUG )
      cout<<"RpMCAna::Init() called"<<endl;
   for (int iSet = 0; iSet < DATA; ++iSet)
   {
      mElAna[iSet] = new ElasticAna(mOutFile);
      mElAna[iSet]->SetAnaName(mUtil->dataSetName(iSet));
      mElAna[iSet]->SetTriggers(&noTriggers);  
      mElAna[iSet]->Init();
      mElAna[iSet]->InitRPMCInfo();

      // RP efficiency
      mOutFile->mkdir("RPMC_"+mUtil->dataSetName(iSet))->cd();
      hMCEffFlow[iSet] = new TH1D("AnaFlow_" + anaName, "CutsFlow", kMax-1, 1, kMax);
      for(int tb=1; tb<kMax; ++tb) 
         hMCEffFlow[iSet]->GetXaxis()->SetBinLabel(tb, mMCEffFlowCutsName[tb-1]);

      hMCtoRecoDiff[iSet] = new TH1D("hMCtoRecoDiff_"+mUtil->dataSetName(iSet),"hMCtoRecoDiff",1000,0,1);
      for (int iBr = 0; iBr < nBranches; ++iBr){
         hEffPxPy[0][iBr][iSet] = new TH2D("hEffPxPy_"+mUtil->branchName(iBr) + "_Total_" + mUtil->dataSetName(iSet),"Sample total",100,-1,1,100,-1,1);
         hEffPxPy[1][iBr][iSet] = new TH2D("hEffPxPy_"+mUtil->branchName(iBr) + "_Passed_" + mUtil->dataSetName(iSet),"Sample passed",100,-1,1,100,-1,1);
         hEffPxPy[2][iBr][iSet] = new TH2D("hEffPxPy_"+mUtil->branchName(iBr) + "_PassedDist_" + mUtil->dataSetName(iSet),"Sample passed",100,-1,1,100,-1,1);         
      }
      for (int iRp = 0; iRp < nRomanPots; ++iRp){
         hEffXY[0][iRp][iSet] = new TH2D("hEffXY_"+mUtil->rpName(iRp) + "_Total_" + mUtil->dataSetName(iSet),"Sample total",30,-0.05,0.05,45,-0.075,0.075);
         hEffXY[1][iRp][iSet] = new TH2D("hEffXY_"+mUtil->rpName(iRp) + "_Passed_" + mUtil->dataSetName(iSet),"Sample passed",30,-0.05,0.05,45,-0.075,0.075);
         hEffXY[2][iRp][iSet] = new TH2D("hEffXY_"+mUtil->rpName(iRp) + "_PassedDist_" + mUtil->dataSetName(iSet),"Sample passed",30,-0.05,0.05,45,-0.075,0.075);
         hEffXYOffSub[0][iRp][iSet] = new TH2D("hEffXYOffSub_"+mUtil->rpName(iRp) + "_Total_" + mUtil->dataSetName(iSet),"Sample total",30,-0.05,0.05,45,-0.075,0.075);
         hEffXYOffSub[1][iRp][iSet] = new TH2D("hEffXYOffSub_"+mUtil->rpName(iRp) + "_Passed_" + mUtil->dataSetName(iSet),"Sample passed",30,-0.05,0.05,45,-0.075,0.075);
         hEffXYOffSub[2][iRp][iSet] = new TH2D("hEffXYOffSub_"+mUtil->rpName(iRp) + "_PassedDist_" + mUtil->dataSetName(iSet),"Sample passed",30,-0.05,0.05,45,-0.075,0.075);
      }
      hNTpPerRP[iSet] = new TH2D("hNTpPerRP_"+mUtil->dataSetName(iSet),"Number of TPs per RP",8,-0.5,7.5,30,-0.5,29.5);
      hNClustersPerPlane[iSet] = new TH2D("hNClustersPerPlane_"+mUtil->dataSetName(iSet),"Number of clusters per plane",32,-0.5,31.5,30,-0.5,29.5);
      for(int iRp=0; iRp<nRomanPots; ++iRp)
      {
         hNTpPerRP[iSet]->GetXaxis()->SetBinLabel(iRp+1, mUtil->rpName(iRp));
         for(int iPl=0; iPl<nPlanes; ++iPl) 
            hNClustersPerPlane[iSet]->GetXaxis()->SetBinLabel(iRp*nPlanes + iPl + 1, mUtil->rpName(iRp)+"_"+mUtil->planeName(iPl));
      }

      mOutFile->cd();
   }
}

void RpMCAna::SetRpPosition(TVector3 (&corr)[nRomanPots], TVector3 (&offsets)[nRomanPots])
{
   mCorrection = corr;
   mOffSet = offsets;
   embedMaker->setOffsetCorrections(corr);
}

void RpMCAna::SetMCInfo(double (&mc_vrtx)[nCoordinates], double (&mc_p)[nCoordinates][nSides])
{
   mVertex = mc_vrtx;
   mMomentum = mc_p;
}

bool RpMCAna::IsTrackInRp(int side)
{
   if( abs(mVertex[Z]) > vertexRange/100 )
      return false;


   int branch = (mMomentum[Z][side] < 0) ? EU : WU; 
   if( mMomentum[Y][side] < 0)
      branch++;

   for (int iRp = 0; iRp < nStationPerSide; ++iRp)
   {
      int RP = mUtil->rpPerBranchStationOrder(branch, iRp); 
      double x, y;
      ProjectToRP(x,y,RP,side);

      //if(!IsInRpRange(x,y,RP,mOffSet[RP]))
      //   return false;
   }

  return true;
}

void RpMCAna::ProjectToRP(double& x, double& y, int RP, int side)
{
   x = mVertex[X]/1000 + (mUtil->rpZPosition(RP) - mVertex[Z])*(mMomentum[X][side]/mMomentum[Z][side]);
   y = mVertex[Y]/1000 + (mUtil->rpZPosition(RP) - mVertex[Z])*(mMomentum[Y][side]/mMomentum[Z][side]);
}


bool RpMCAna::AreTracksInRp()
{
  return IsTrackInRp(East) && IsTrackInRp(West);
}

void RpMCAna::runEmbedding(bool embedding)
{
   embedMaker->setMCEvent(mMcEvt);
   embedMaker->setZBEvent(mRpEvt);
   mEmbedEvt = new StRPEvent(*mRpEvt);
   mEmbedEvt->clearEvent();
   embedMaker->setRPEvent(mEmbedEvt);
   embedMaker->setEmbedding(embedding);
   embedMaker->setVertex( mVertex[X], mVertex[Y], mVertex[Z]);
   embedMaker->MakeTracks(mUtil->p0(),mUtil->p0());
}

void RpMCAna::runMCEff(int set)
{
   int branch;
   AnaRpTracks(mEmbedEvt);
   // plot efficiency as function of momenta
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      hMCEffFlow[set]->Fill(kAll);
      if(!IsTrackInRp(iSide))
         continue;

      hMCEffFlow[set]->Fill(kInside);
      branch = (mMomentum[Z][iSide] < 0) ? EU : WU; 
      if( mMomentum[Y][iSide] < 0)
         branch++;

      //Is there a trigger signal?
      //if( !(IsRpTrigBit(mEmbedEvt, mUtil->rpPerBranchStationOrder( branch, RP1) ) || 
      //   IsRpTrigBit(mEmbedEvt, mUtil->rpPerBranchStationOrder( branch, RP2) )))
      //   continue;
      hMCEffFlow[set]->Fill(kTrigger);
      // fill sample total 
      hEffPxPy[0][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);

      if( mRpTrackIdVec_perBranch[branch].size()!=1 )
         continue;
      hMCEffFlow[set]->Fill(kOneTrack);      
      StUPCRpsTrack* trk = mEmbedEvt->getTrack(mRpTrackIdVec_perBranch[branch][0]);
      if(!trk) continue;

      double MCtoRecoDiff = sqrt( pow(trk->pVec().X() - mMomentum[X][iSide], 2) + pow(trk->pVec().Y() - mMomentum[Y][iSide], 2));
      // plot difference between MC and reconstructed hit
      if(RPInFidRange(mMomentum[X][iSide], mMomentum[Y][iSide]))
         hMCtoRecoDiff[set]->Fill( MCtoRecoDiff);

      hEffPxPy[1][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);
      if(MCtoRecoDiff > maxMcToRecoDist)
         continue;

      hMCEffFlow[set]->Fill(kDist);
      hEffPxPy[2][branch][set]->Fill(mMomentum[X][iSide], mMomentum[Y][iSide]);
   }
   // plot efficiency as function of RP coordinates (x,y)
   for (int iSide = 0; iSide < nSides; ++iSide)
   {
      branch = (mMomentum[Z][iSide] < 0) ? EU : WU; 
      if( mMomentum[Y][iSide] < 0)
         branch++;
      // fill sample total 
      double xMC[nStationPerSide], yMC[nStationPerSide];
      for (int iRp = 0; iRp < nStationPerSide; ++iRp)
      {
         int RP = mUtil->rpPerBranchStationOrder(branch, iRp); 
         ProjectToRP(xMC[iRp],yMC[iRp],RP,iSide);
         hEffXY[0][RP][set]->Fill(xMC[iRp],yMC[iRp]);
         hEffXYOffSub[0][RP][set]->Fill(xMC[iRp]-mOffSet[RP][X],yMC[iRp]-mOffSet[RP][Y]);
         if( IsInRpRange( xMC[iRp], yMC[iRp], RP, mOffSet[RP]))
         { // plot # hits in plane (x vs y) and #TPs in RP in general
            hNTpPerRP[set]->Fill(RP,mTrackPointIdVec[RP].size());
            for (int iPlane = 0; iPlane < nPlanes; ++iPlane)
               hNClustersPerPlane[set]->Fill(RP*nPlanes + iPlane, mClusterIdVec[RP*nPlanes + iPlane].size() );
         }

         if( mTrackPointIdVec[RP].size() !=1 )
            continue;

         StUPCRpsTrackPoint *trkPoint = mEmbedEvt->getTrackPoint(mTrackPointIdVec[RP][0]);
         if(!trkPoint) continue;
         double x = trkPoint->x();
         double y = trkPoint->y();
         hEffXY[1][RP][set]->Fill(xMC[iRp],yMC[iRp]);
         hEffXYOffSub[1][RP][set]->Fill(xMC[iRp]-mOffSet[RP][X],yMC[iRp]-mOffSet[RP][Y]);

         double MCtoRecoDiff = sqrt( pow(x - xMC[iRp], 2) + pow(y - yMC[iRp], 2));
         if(MCtoRecoDiff > 0.01)
            continue;
         hEffXY[2][RP][set]->Fill(xMC[iRp],yMC[iRp]);
         hEffXYOffSub[2][RP][set]->Fill(xMC[iRp]-mOffSet[RP][X],yMC[iRp]-mOffSet[RP][Y]);
      }      
   }

}

const TString RpMCAna::mMCEffFlowCutsName[kMax] = { TString("All"), TString("Inside RP"), TString("Trigger"),
               TString("1 track"), TString("Distance")};
