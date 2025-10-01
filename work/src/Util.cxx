#include "../include/Util.h"

Util::Util(): mSpeedOfLight(299792458), mBeamMomentum(254.867), mPi(3.14159265359), mEpsilon(1e-7){
   mSideName = new TString[nSides];
   mSideName[E] = TString("East");
   mSideName[W] = TString("West");
   
   mCoordinateName = new TString[nCoordinates];
   mCoordinateName[X] = TString("x");
   mCoordinateName[Y] = TString("y");
   mCoordinateName[Z] = TString("z");
   
   mTriggerName = new TString[nTriggers];
   mTriggerName[UPC_v1] = TString("UPC"); // 570702 RP_UPC
   mTriggerName[UPC_v2] = TString("UPC"); // 570712 RP_UPC
   mTriggerName[SDT] = TString("SDT"); // 570703 RP_SDT
   mTriggerName[ET_v1] = TString("ET"); // 570709 RP_ET
   mTriggerName[ET_v2] = TString("ET"); // 570719 RP_ET 
   mTriggerName[CPT2_v1] = TString("CPT2"); // 570701 RP_CPT2
   mTriggerName[CPT2_v2] = TString("CPT2"); // 570711 RP_CPT2
   mTriggerName[CPT2noBBCL] = TString("CPT2noBBCL"); // 570705 RP_CPT2noBBCL 
   mTriggerName[Zerobias] = TString("Zerobias"); // 570704 RP_Zerobias
   mTriggerName[JPsiHTTP_v1] = TString("JPsi*HTTP"); // 570209 JPsi*HTTP 
   mTriggerName[JPsiHTTP_v2] = TString("JPsi*HTTP"); // 570219 JPsi*HTTP
   mTriggerName[JPsiHTTP_v3] = TString("JPsi*HTTP"); // 570229 JPsi*HTTP
   mTriggerName[SDT_RHICf] = TString("SDT"); // 590703 RP_SDT
   mTriggerName[ET_RHICf] = TString("ET"); // 590709 RP_ET
   mTriggerName[CPT2_RHICf] = TString("CPT2"); // 590701 RP_CPT2
   mTriggerName[CPT2noBBCL_RHICf_v1] = TString("CPT2noBBCL"); // 590705 RP_CPT2noBBCL
   mTriggerName[CPT2noBBCL_RHICf_v2] = TString("CPT2noBBCL"); // 590708 RP_CPTnoBBCL
      
   mArmName = new TString[nArms];
   mArmName[EU_WD] = TString("EU-WD");
   mArmName[ED_WU] = TString("ED-WU");
   
   mBranchName = new TString[nBranches];
   mBranchName[EU] = TString("EU");
   mBranchName[ED] = TString("ED");
   mBranchName[WU] = TString("WU");
   mBranchName[WD] = TString("WD");
   
   mBranchesConfigurationName = new TString[nBranchesConfigurations];
   mBranchesConfigurationName[CONF_EU_WU] = TString("EU-WU");
   mBranchesConfigurationName[CONF_ED_WD] = TString("ED-WD");
   mBranchesConfigurationName[CONF_EU_WD] = TString("EU-WD");
   mBranchesConfigurationName[CONF_ED_WU] = TString("ED-WU");
      
   mRpName = new TString[nRomanPots];
   mRpName[E1U] = TString("E1U");
   mRpName[E1D] = TString("E1D");
   mRpName[E2U] = TString("E2U");
   mRpName[E2D] = TString("E2D");
   mRpName[W1U] = TString("W1U");
   mRpName[W1D] = TString("W1D");
   mRpName[W2U] = TString("W2U");
   mRpName[W2D] = TString("W2D");

   mRpConfigName = new TString[nRpConfigurations];
   mRpConfigName[IT] = TString("IT");
   mRpConfigName[ET] = TString("ET");

   mPlaneName = new TString[nPlanes];
   mPlaneName[A] = TString("A");
   mPlaneName[B] = TString("B");
   mPlaneName[C] = TString("C");
   mPlaneName[D] = TString("D");
   
   mStationName = new TString[nStations];
   mStationName[E1] = TString("E1");
   mStationName[E2] = TString("E2");
   mStationName[W1] = TString("W1");
   mStationName[W2] = TString("W2");
   
   mParticleName = new TString[nParticles];
   mParticleName[ELECTRON] = TString("electron");
//   mParticleName[MUON] = TString("muon");
   mParticleName[PION] = TString("pion");
   mParticleName[KAON] = TString("kaon");
   mParticleName[PROTON] = TString("proton");

   mNeutralParticleName = new TString[nNeutralParticles];
   mNeutralParticleName[K0S] = TString("K^{0}_{S}");
   mNeutralParticleName[LAMBDA] = TString("#Lambda");
   mNeutralParticleName[LAMBDABAR] = TString("#bar{#Lambda}");
   
   mTpcTrackTypeName = new TString[nTpcTrkTypes];
   mTpcTrackTypeName[GLO] = TString("global");
   mTpcTrackTypeName[PRI] = TString("primary");
   mTpcTrackTypeName[TOF] = TString("TofMatched");
   mTpcTrackTypeName[QUA] = TString("goodQuality");

   mBunchCrossingTypeName = new TString[nBnchXngsTypes];
   mBunchCrossingTypeName[CB] = TString("collidingBunches");
   mBunchCrossingTypeName[AG] = TString("abortGaps");
   
   mChargeSum2TrksName = new TString[nCharges2Trks];
   mChargeSum2TrksName[OPPO] = TString("oppositeSign");
   mChargeSum2TrksName[SAME] = TString("sameSign");
   
   mSignName = new TString[nSigns];
   mSignName[PLUS] = TString("Plus");
   mSignName[MINUS] = TString("Minus");
   
   mEfficiencyName = new TString[nEffCorrections];
   mEfficiencyName[RPACC] = TString("RpAcc");
   mEfficiencyName[TPCRECOEFF] = TString("TpcRecoEff");
   mEfficiencyName[TOFMATCHEFF] = TString("TofMatchEff");

   mCutName = new TString[nAnalysisCuts];
   mCutName[ALL] = TString("All");
   mCutName[TRIG] = TString("CPT");
   mCutName[TWORPTRKS] = TString("2 RP tracks");
   mCutName[INFID] = TString("RP fiducial");
   mCutName[ONETOFVX] = TString("1 TOF vertex");
   mCutName[ZVERTEX] = TString("|V_{z}| < 80 cm");
   mCutName[TWOTOFTRKS] = TString("2 TOF tracks");
   mCutName[ETA] = TString("|#eta| < 0.7");
   mCutName[OPPOSITE] = TString("TotQ = 0");
   //mCutName[TPCRPVX_MATCHED] = TString("TPC-RP vrtx match");
   mCutName[PIPI] = TString("#pi^{+}#pi^{-}");
   mCutName[PPI] = TString("p#pi^{-}");
   mCutName[PIPBAR] = TString("#pi^{+}#bar{p}");

   mV0CutName = new TString[nV0SelectionCuts];
   mV0CutName[V0ALL] = TString("All");
   mV0CutName[V0TRIG] = TString("CPT");
   //mV0CutName[V0ETA] = TString("|#eta| < 0.7");
   mV0CutName[V0PID] = TString("PID");
   mV0CutName[V0FLAG] = TString("V0 flag");
   mV0CutName[V0PAIR] = TString("V0 pair selection");
   mV0CutName[V0ETAVTXZ] = TString("eta + vertex range");
   mV0CutName[V0OPPOSITE] = TString("Unlike-sign");
   mV0CutName[V0PIPI] = TString("#pi^{+}#pi^{-}");
   mV0CutName[V0PPI] = TString("p#pi^{-}");
   mV0CutName[V0PIPBAR] = TString("#pi^{+}#bar{p}");
   
   mTOFEFFName = new TString[nTOFEFFCuts];
   mTOFEFFName[TOFALL] = TString("All");
   mTOFEFFName[TOFTRIG] = TString("CPT");
   mTOFEFFName[TOFTRACKQUALITY] = TString("Good track quality");
   mTOFEFFName[TOFETAVTXZ] = TString("#eta - Vtx_{Z}");
   mTOFEFFName[TOFPAIR] = TString("topology");
   mTOFEFFName[TOFOPPOSITE] = TString("Q_{tot} = 0");

   mJPSICutName = new TString[nJPSISelectionCuts];
   mJPSICutName[JPSIALL] = TString("All");
   mJPSICutName[JPSITRIG] = TString("J/#Psi trigger");
   mJPSICutName[JPSI1VTX] = TString("1 Vertex");
   mJPSICutName[JPSIBEMC] = TString("2 BEMC tracks");
   mJPSICutName[JPSIVTXZETA] = TString("V_{Z} - #eta cut");
   mJPSICutName[JPSIBACKTOBACK] = TString("Back-to-back");
   mJPSICutName[JPSIPID] = TString("PID");
   mJPSICutName[JPSI1RP] = TString("1 RP track");
   mJPSICutName[JPSIRPFIDCUT] = TString("Fiducial RP cut");
   mJPSICutName[JPSIQTOT] = TString("Q_{tot} = 0");

   mJPSI2CutName = new TString[nJPSI2SelectionCuts];
   mJPSI2CutName[JPSI2ALL] = TString("All");
   mJPSI2CutName[JPSI2TRIG] = TString("J/#Psi trigger");
   mJPSI2CutName[JPSI21VTX] = TString("1 vertex");
   mJPSI2CutName[JPSI24TRACKS] = TString("4 good tracks");
   mJPSI2CutName[JPSI2VTXZETA] = TString("V_{Z} - #eta cut");
   mJPSI2CutName[JPSI2PID] = TString("PID");
   mJPSI2CutName[JPSI2QTOT] = TString("Q_{tot} = 0");


   mEmbeddingName = new TString[nEmbeddingCuts];
   mEmbeddingName[EMBEDDINGALL] = TString("All");
   mEmbeddingName[EMBEDDING2BEMC] = TString("2 BEMC tracks");
   mEmbeddingName[EMBEDDINGBACKTOBACK] = TString("Back-to-back");
   mEmbeddingName[EMBEDDINGETA] = TString("|#eta| < 0.9");
   mEmbeddingName[EMBEDDINGPID] = TString("PID");
   mEmbeddingName[EMBEDDINGQTOT] = TString("Q_{tot} = 0");

   mZBCutName = new TString[nZBCuts];
   mZBCutName[ZBALL] = TString("All");
   mZBCutName[ZBTRIGGER] = TString("RP_zerobias Trigger");
   mZBCutName[ZBBEMCTRACKS] = TString("Good BEMC tracks");
   mZBCutName[ZBPID] = TString("PID");
   mZBCutName[ZBBACKTOBACK] = TString("Back To Back"); 
   mZBCutName[ZBVTXZETA] = TString("V_{Z} - #eta");
   mZBCutName[ZBTRIGGERCONDITION] = TString("Trigger condition");


   mBECutName = new TString[nBECuts];
   mBECutName[BEALL] = TString("All");
   mBECutName[BETRIG] = TString("trigger");
   mBECutName[BE1VTX] = TString("1 Vertex");
   mBECutName[BE2TOF] = TString("2 TOF tracks");
   mBECutName[BEPROJECTION] = TString("Projected to BEMC");
   mBECutName[BEDELTADIPANGLE] = TString("#delta -dip-angle");
   mBECutName[BEETAVTXZ] = TString("#eta - V_{Z}");
   mBECutName[BEPID] = TString("PID");
   mBECutName[BEINVMASS] = TString("m_{ee} < 0.1 GeV/c^{2}");
   mBECutName[BEQTOT] = TString("Q_{tot} = 0");


   mDataSetName = new TString[nDataSets];
   mDataSetName[MC] = TString("PureMc");
   mDataSetName[MCZB] = TString("MCZB");
   mDataSetName[DATA] = TString("Data");

   mDataTagName = new TString[nDataTag];
   mDataTagName[TRUEMC] = TString("TrueMc");
   mDataTagName[RECO] = TString("Reco"); 
   
   mParticleMass[PION] = 0.13957018;
   mParticleMass[KAON] = 0.493667;
   mParticleMass[PROTON] = 0.93827208;
   mParticleMass[ELECTRON] = 0.000510998;
   
   mBranchPerRp[E1U] = EU;
   mBranchPerRp[E2U] = EU;
   mBranchPerRp[E1D] = ED;
   mBranchPerRp[E2D] = ED;
   mBranchPerRp[W1U] = WU;
   mBranchPerRp[W2U] = WU;
   mBranchPerRp[W1D] = WD;
   mBranchPerRp[W2D] = WD;

   mRpZPosition[E1U] = -15.78;
   mRpZPosition[E2U] = -17.58;
   mRpZPosition[E1D] = -15.78;
   mRpZPosition[E2D] = -17.58;
   mRpZPosition[W1U] = 15.78;
   mRpZPosition[W2U] = 17.58;
   mRpZPosition[W1D] = 15.78;
   mRpZPosition[W2D] = 17.58;
   
   mOppositeBranch[EU] = WD;
   mOppositeBranch[ED] = WU;
   mOppositeBranch[WU] = ED;
   mOppositeBranch[WD] = EU;
   
   mSidePerRp[E1U] = E;
   mSidePerRp[E2U] = E;
   mSidePerRp[E1D] = E;
   mSidePerRp[E2D] = E;
   mSidePerRp[W1U] = W;
   mSidePerRp[W2U] = W;
   mSidePerRp[W1D] = W;
   mSidePerRp[W2D] = W;
   
   mStationOrderPerRp[E1U] = RP1;
   mStationOrderPerRp[E2U] = RP2;
   mStationOrderPerRp[E1D] = RP1;
   mStationOrderPerRp[E2D] = RP2;
   mStationOrderPerRp[W1U] = RP1;
   mStationOrderPerRp[W2U] = RP2;
   mStationOrderPerRp[W1D] = RP1;
   mStationOrderPerRp[W2D] = RP2;
   
   mRpPerBranchStationOrder[EU][RP1] = E1U;
   mRpPerBranchStationOrder[EU][RP2] = E2U;
   mRpPerBranchStationOrder[ED][RP1] = E1D;
   mRpPerBranchStationOrder[ED][RP2] = E2D;
   mRpPerBranchStationOrder[WU][RP1] = W1U;
   mRpPerBranchStationOrder[WU][RP2] = W2U;
   mRpPerBranchStationOrder[WD][RP1] = W1D;
   mRpPerBranchStationOrder[WD][RP2] = W2D;

   mBranchPerBranchConfiguration[CONF_EU_WU][E] = EU;
   mBranchPerBranchConfiguration[CONF_EU_WU][W] = WU;
   mBranchPerBranchConfiguration[CONF_ED_WD][E] = ED;
   mBranchPerBranchConfiguration[CONF_ED_WD][W] = WD;
   mBranchPerBranchConfiguration[CONF_EU_WD][E] = EU;
   mBranchPerBranchConfiguration[CONF_EU_WD][W] = WD;
   mBranchPerBranchConfiguration[CONF_ED_WU][E] = ED;
   mBranchPerBranchConfiguration[CONF_ED_WU][W] = WU;  

   mBranchConfigMap[Up][Up] =  CONF_EU_WU;
   mBranchConfigMap[Down][Down] =  CONF_ED_WD;
   mBranchConfigMap[Up][Down] =  CONF_EU_WD;
   mBranchConfigMap[Down][Up] =  CONF_ED_WU;
}

Util::~Util(){
   if(mSideName) delete [] mSideName;
   if(mCoordinateName) delete [] mCoordinateName;
   if(mTriggerName) delete [] mTriggerName;
   if(mArmName) delete [] mArmName;
   if(mBranchName) delete [] mBranchName;
   if(mBranchesConfigurationName) delete [] mBranchesConfigurationName;
   if(mRpName) delete [] mRpName;
   if(mRpConfigName) delete [] mRpConfigName;
   if(mPlaneName) delete [] mPlaneName;
   if(mStationName) delete [] mStationName;
   if(mParticleName) delete [] mParticleName;
   if(mTpcTrackTypeName) delete [] mTpcTrackTypeName;
   if(mBunchCrossingTypeName) delete [] mBunchCrossingTypeName;
   if(mChargeSum2TrksName) delete [] mChargeSum2TrksName;
   if(mSignName) delete [] mSignName;
   if(mEfficiencyName) delete [] mEfficiencyName;
   if(mCutName) delete [] mCutName;
}


Double_t Util::bkgdFraction(const TH1* hPtMiss, Double_t ptMissCut, Double_t fitLimitMin, Double_t fitLimitMax, TF1 *funcReturn) const{
   static int funcId; ++funcId; 
   TString idStr; idStr.Form("func%d", funcId);
   TF1* func = new TF1(idStr, "[0]*x+[1]*x*x", 0, fitLimitMax);
   func->SetNpx(1e3);
   func->SetParameter(0, 1);
   func->SetParameter(1, -1);
   double bkgdFrac = 1e9;
   if(!hPtMiss){
      std::cerr << "ERROR in getBkgdFraction()" << std::endl;
      return bkgdFrac;
   }

   const int firstBinInFitRange = hPtMiss->GetXaxis()->FindBin(fitLimitMin);
   const int lastBinInFitRange = hPtMiss->GetXaxis()->FindBin(fitLimitMax)-1;
   const int nBins = lastBinInFitRange - firstBinInFitRange + 1;
   TH1D* hTemp = new TH1D();
   hPtMiss->Copy( *hTemp );
   hTemp->SetName("VeryTemp");
   const int nRebin = static_cast<int>( nBins/4 ); // so that there are 4 bins to fit
   if( nRebin>1 )
      hTemp->Rebin(nRebin);
   hTemp->Fit( func, "MQNI", "", fitLimitMin, fitLimitMax);
   if( fabs(func->GetParameter(0)-1) < 0.001 && fabs(func->GetParameter(1)+1) < 0.001 ){
      bkgdFrac = 0;
   } else{
      double integralOfFunc = func->Integral(0, ptMissCut);
      if( integralOfFunc < 0 ){
         func->FixParameter(1, 0);
         hTemp->Fit( func, "MQNI", "", fitLimitMin, fitLimitMax);
         integralOfFunc = func->Integral(0, ptMissCut);
      }
      double nBkgdEvents = integralOfFunc / hPtMiss->GetBinWidth(1);
      double nAllEvents = integratePtMiss(hPtMiss, ptMissCut);
      bkgdFrac = nBkgdEvents/nAllEvents;
   }
   if( nRebin>1 ){
      func->SetParameter(0, func->GetParameter(0)/nRebin );
      func->SetParameter(1, func->GetParameter(1)/nRebin );
   }

   delete hTemp;
   
   if(funcReturn) 
      *funcReturn = *func;
   else 
      delete func;
   
   return bkgdFrac;
}


Double_t Util::integratePtMiss(const TH1* hPtMiss, Double_t ptMissCut) const{
   int ptMissCutBin = hPtMiss->FindFirstBinAbove(ptMissCut);
   double integral = hPtMiss->Integral(1, ptMissCutBin-1);
   integral += hPtMiss->GetBinContent(ptMissCutBin) * ( ptMissCut - hPtMiss->GetXaxis()->GetBinLowEdge(ptMissCutBin) ) / hPtMiss->GetXaxis()->GetBinWidth(ptMissCutBin);
   return integral;
}



TH1D* Util::bkgdHistogram(const TH2* hMissPtVsX, Double_t ptMissCut, Double_t fitLimitMin, Double_t fitLimitMax, Int_t mode, vector<TF1*> * funcVec) const{
   std::vector<float> binsVector = getBinsVectorF( hMissPtVsX->GetXaxis() );
   TH1D *hBkgd = new TH1D( TString("Background_")+hMissPtVsX->GetName(), TString("Background_")+hMissPtVsX->GetTitle(), binsVector.size()-1, &(binsVector[0]) );
   //BEGIN             
   TH1D *missingPtHist;
   //previous version:
   missingPtHist = hMissPtVsX->ProjectionY("tmpNameProjectionY");
   //current version (added)
   double backgroundIntegral_FullPtMiss = missingPtHist->Integral( missingPtHist->FindBin( fitLimitMin ), missingPtHist->FindBin( fitLimitMax ) );
   missingPtHist->Reset("ICESM");
   for(int i=0; i<=(hMissPtVsX->GetNbinsX()+1); ++i){
      TH1D* singleBinProjection = hMissPtVsX->ProjectionY("tmpNameProjectionY_singleBin", i, i);
      double singnalRegionIntegral = integratePtMiss( singleBinProjection, ptMissCut );
      if( singnalRegionIntegral>0 )
         missingPtHist->Add( singleBinProjection );
   }
   double backgroundIntegral_SlicedPtMiss = missingPtHist->Integral( missingPtHist->FindBin( fitLimitMin ), missingPtHist->FindBin( fitLimitMax ) );
   //END               
   const double integralLimit_MIN = missingPtHist->GetXaxis()->GetBinLowEdge( missingPtHist->GetXaxis()->FindBin( fitLimitMin ) );
   const double integralLimit_MAX = missingPtHist->GetXaxis()->GetBinUpEdge( missingPtHist->GetXaxis()->FindBin( fitLimitMax ) );
   if( mode==0 ){ // default mode
      TF1 *pTMissExtrapolationFunc = new TF1();
      //double bkgdFrac = bkgdFraction( missingPtHist, ptMissCut, fitLimitMin, fitLimitMax, pTMissExtrapolationFunc);
      const double integralRatio_signalRegion_to_bkgdFreeRegion = pTMissExtrapolationFunc->Integral(0, ptMissCut) / pTMissExtrapolationFunc->Integral(integralLimit_MIN, integralLimit_MAX);
      //     double nBkgdEvents = bkgdFrac*integratePtMiss( missingPtHist, ptMissCut );
      for(int i=0; i<=(hBkgd->GetNbinsX()+1); ++i){
         TH1D* hMissPtProj = hMissPtVsX->ProjectionY("tmpNameProjectionY_singleBin", i, i);
         double nBkgdInPtMissBin = (backgroundIntegral_FullPtMiss / backgroundIntegral_SlicedPtMiss) * integralRatio_signalRegion_to_bkgdFreeRegion * hMissPtProj->Integral( hMissPtProj->FindBin( fitLimitMin ), hMissPtProj->FindBin( fitLimitMax ) );
         double singnalRegionIntegral = integratePtMiss( hMissPtProj, ptMissCut );
         hBkgd->SetBinContent(i, nBkgdInPtMissBin > singnalRegionIntegral ? singnalRegionIntegral : nBkgdInPtMissBin );
         hBkgd->SetBinError(i, 0.0001 * hBkgd->GetBinContent(i) );
//       hBkgd->SetBinError(i, sqrt(hBkgd->GetBinContent(i)) ); //ALERT
      }
      if(funcVec) funcVec->push_back( pTMissExtrapolationFunc );
   } else{
      for(int i=0; i<=(hBkgd->GetNbinsX()+1); ++i){
         TF1 *pTMissExtrapolationFunc = new TF1();
         TH1D* hMissPtProj = hMissPtVsX->ProjectionY("tmpNameProjectionY_singleBin", i, i);
         double bkgdFrac = bkgdFraction( hMissPtProj, ptMissCut, fitLimitMin, fitLimitMax, pTMissExtrapolationFunc);
         double nBkgdEvents = bkgdFrac*integratePtMiss( hMissPtProj, ptMissCut );
         hBkgd->SetBinContent(i, nBkgdEvents);
         hBkgd->SetBinError(i, 0);
         if(funcVec) funcVec->push_back( pTMissExtrapolationFunc );
      }
   }
   return hBkgd;
}



TH1D* Util::backgroundFromSameSignTemplate( TH2D* const hPtMissVsX, Double_t ptMissCut, TH1D* const hPtMissSameSign) const{
   std::vector<float> binsVector = getBinsVectorF( hPtMissVsX->GetXaxis() );
   TH1D *hBkgd = new TH1D( TString("BackgroundFromSameSignTemplate_")+hPtMissVsX->GetName(), TString("BackgroundFromSameSignTemplate_")+hPtMissVsX->GetTitle(), binsVector.size()-1, &(binsVector[0]) );
   const double minFitLimit = 0.16;
   const double maxFitLimit = 0.28;
   const double integralRatio_signalRegion_to_bkgdFreeRegion = integratePtMiss( hPtMissSameSign, ptMissCut ) / hPtMissSameSign->Integral(hPtMissSameSign->FindBin(minFitLimit), hPtMissSameSign->FindBin(maxFitLimit) );
   for(int i=0; i<=(hBkgd->GetNbinsX()+1); ++i){
      TH1D* hMissPtProj = hPtMissVsX->ProjectionY("tmpTmpNameProjectionY_singleBin", i, i);
      hBkgd->SetBinContent(i, integralRatio_signalRegion_to_bkgdFreeRegion * hMissPtProj->Integral(hMissPtProj->FindBin(minFitLimit), hMissPtProj->FindBin(maxFitLimit) ));
      hBkgd->SetBinError(i, 0.0001 * hBkgd->GetBinContent(i) );
   }
   return hBkgd;
}



Double_t Util::binomialCoeff(UInt_t n, UInt_t k) const{
   if(k > n)
      { cerr << "ERROR in Util::binomialCoeff(UInt_t n, UInt_t k):  k>n !!!!" << endl; return 1; }
   double result = 1.0;
   for(unsigned int i=1; i<=k; ++i) 
      result *= static_cast<double>(n-i+1)/static_cast<double>(i);
   return result;
}


void Util::subtractBackground(TH1* hSignalPlusBkgd, const TH1* hBkgd) const{
   if( hSignalPlusBkgd->GetNbinsX() != hBkgd->GetNbinsX() ){
      std::cerr << "ERROR in subtractBackground(): number of bins of signal and background histogram is different. Return" << std::endl;
      return;
   }
   for(int i=0; i<=(hSignalPlusBkgd->GetNbinsX()+1); ++i){
      const double signalPlusBkgdContent = hSignalPlusBkgd->GetBinContent(i);
      if( signalPlusBkgdContent<1e-9 ) continue;
      const double signalPlusBkgdError = hSignalPlusBkgd->GetBinError(i);
      const double backgroundContent = hBkgd->GetBinContent(i);
      const double bkgdFraction = backgroundContent / signalPlusBkgdContent;
      //const double bkgdFractionError = sqrt(bkgdFraction*(1.-bkgdFraction) / signalPlusBkgdContent ); //ALERT binomial error approximation
      
      const double pureSignalContent = signalPlusBkgdContent - backgroundContent;
      if( pureSignalContent<0 ) std::cerr << "WARNING in subtractBackground(): Background larger than signal" << std::endl;
      hSignalPlusBkgd->SetBinContent(i, pureSignalContent>0 ? pureSignalContent : 1e-9 ); //ALERT
//     hSignalPlusBkgd->SetBinError(i, pureSignalContent>0 ? (signalPlusBkgdError*pureSignalContent/signalPlusBkgdContent) : signalPlusBkgdError ); // ALERT     
      double signalError = sqrt( (1.+bkgdFraction)*(1.+bkgdFraction)*signalPlusBkgdError*signalPlusBkgdError   /*+   signalPlusBkgdContent*signalPlusBkgdContent*bkgdFractionError*bkgdFractionError*/ );
      hSignalPlusBkgd->SetBinError(i, signalError );
   }
}      

std::vector<float> Util::getBinsVectorF(const TAxis *axis) const{
   std::vector<float> binsVector;
   for( int i=1; i<=axis->GetNbins(); ++i)
      binsVector.push_back( axis->GetBinLowEdge( i ) );
   binsVector.push_back( axis->GetBinUpEdge( axis->GetNbins() ) );
   return binsVector;
}


std::vector<double> Util::getBinsVectorD(const TAxis *axis) const{
   std::vector<double> binsVector;
   for( int i=1; i<=axis->GetNbins(); ++i)
      binsVector.push_back( axis->GetBinLowEdge( i ) );
   binsVector.push_back( axis->GetBinUpEdge( axis->GetNbins() ) );
   return binsVector;
}

// fit line through n-track points and return position of the fit in z = positionOfFit
TVector3 Util::fitLine(const vector<TVector3> trackPoints, double positionOfFit)
{
   TGraph2D * gr = new TGraph2D();
   //gr->SetPoint(0,x,y,z)
   for (unsigned int iTp = 0; iTp < trackPoints.size(); ++iTp)
      gr->SetPoint(iTp,trackPoints[iTp].X(),trackPoints[iTp].Y(),trackPoints[iTp].Z());
   
   TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
   Double_t arglist[10];
   arglist[0] = -1; // to silence output
   arglist[1] = 0.000001; // tolerance 
   min->ExecuteCommand("SET PRINT",arglist,1); // Silence output
   arglist[0] = 1000; // number of function calls 
   min->SetObjectFit(gr);
   min->SetFCN(SumDistance2);

   double pStart[4] = {1,1,1,1};
   //(Int_t ipar, const char *parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh)
   min->SetParameter(0,"x0",pStart[0],0.01,0,0);
   min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
   min->SetParameter(2,"y0",pStart[2],0.01,0,0);
   min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

   min->ExecuteCommand("MIGRAD",arglist,2);

   // get fit parameters
   double parFit[4];
   for (int i = 0; i < 4; ++i) 
      parFit[i] = min->GetParameter(i);

   TVector3 fitResult;
   fitResult.SetZ(positionOfFit);
   line(fitResult[2],parFit,fitResult[0], fitResult[1]);
   return fitResult;
}

// function to be minimized 
void SumDistance2(int &, double *, double & sum, double * par, int) { 
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) { 
      double d = distance2(x[i],y[i],z[i],par); 
      sum += d;
   }
}

// define the parameteric line equation 
void line(double z, double *p, double &x, double &y) { 
   x = p[0] + p[1]*z; 
   y = p[2] + p[3]*z;
} 

// calculate distance line-point 
double distance2(double x,double y,double z, double *p) { 
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
   ROOT::Math::XYZVector xp(x,y,z); 
   ROOT::Math::XYZVector x0(p[0], p[2], 0. ); 
   ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); 
   ROOT::Math::XYZVector u = (x1-x0).Unit(); 
   double d2 = ((xp-x0).Cross(u)) .Mag2(); 
   return d2; 
}

