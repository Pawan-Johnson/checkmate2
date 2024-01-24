#include "SimpleAnalysisFramework/AnalysisClass.h"

DefineAnalysis(EwkThreeLeptonOnshell2018)

void EwkThreeLeptonOnshell2018::Init() {

  // add regions
  std::cout << "Setting Regions" << std::endl;

  // presel
  addRegions({"presel"});

  // SRs
  addRegions({"DFOSSR1", "DFOSSR2"});

  addRegions({"SFOSOFFSSR1", "SFOSOFFSSR2", "SFOSOFFSSR3", "SFOSOFFSSR4", "SFOSOFFSSR5",
  "SFOSOFFSSR6", "SFOSOFFSSR7", "SFOSOFFSSR8", "SFOSOFFSSR9", "SFOSOFFSSR10",
  "SFOSOFFSSR11", "SFOSOFFSSR12", "SFOSOFFSSR13", "SFOSOFFSSR14", "SFOSOFFSSR15",
  "SFOSOFFSSR16", "SFOSOFFSSR17", "SFOSOFFSSR18", "SFOSOFFSSR19"});

  addRegions({"SFOSONSSR1", "SFOSONSSR2", "SFOSONSSR3", "SFOSONSSR4", "SFOSONSSR5",
  "SFOSONSSR6", "SFOSONSSR7", "SFOSONSSR8", "SFOSONSSR9", "SFOSONSSR10",
  "SFOSONSSR11", "SFOSONSSR12", "SFOSONSSR13", "SFOSONSSR14", "SFOSONSSR15",
  "SFOSONSSR16", "SFOSONSSR17", "SFOSONSSR18", "SFOSONSSR19", "SFOSONSSR20"});

  // CR/VRs
  addRegions({"VR1", "VR2", "VR3"});
  addRegions({"CR1", "CR2", "CR3"});
  addHistogram("cutflow",20,0,20);
  addHistogram("cutflow_DFOS",20,0,20);

  // inclusive SR/VR
  addRegions({"DFOSSR"});
  addRegions({"DFOSSR_loose"});

  addRegions({"SFOSONSSR"});
  addRegions({"SFOSOFFSSR"});
  addRegions({"SFOSOFFSSR_lowMll"});
  addRegions({"SFOSOFFSSR_highMll"});

  addRegions({"SFOSONSSR_j0"});
  addRegions({"SFOSOFFSSR_j0"});
  addRegions({"SFOSOFFSSR_lowMll_j0"});
  addRegions({"SFOSOFFSSR_highMll_j0"});

  addRegions({"SFOSONSSR_j1inc"});
  addRegions({"SFOSOFFSSR_j1inc"});
  addRegions({"SFOSOFFSSR_lowMll_j1inc"});
  addRegions({"SFOSOFFSSR_highMll_j1inc"});

  addRegions({"VR"});
  addRegions({"VR_j0"});
  addRegions({"VR_j1inc"});
}

void EwkThreeLeptonOnshell2018::ProcessEvent(AnalysisEvent *event) {
  //====================================================================================================
  // get baseline (+combi) objects
  // cfr. STConfig: https://gitlab.cern.ch/atlas-phys-susy-wg/AnalysisSUSYToolsConfigurations/-/blob/master/ANA-SUSY-2019-09_Onshell.conf
  auto combiElectrons     = event->getElectrons(4.5, 2.47, ELooseBLLH);
  auto baselineElectrons  = event->getElectrons(10., 2.47, ELooseBLLH | EZ05mm);
  auto combiMuons         = event->getMuons(3.0, 2.7, MuMedium | MuNotCosmic | MuQoPSignificance);
  auto baselineMuons      = event->getMuons(10., 2.5, MuMedium | MuZ05mm | MuQoPSignificance | MuNotCosmic);
  auto baselineJets       = event->getJets(20., 4.5);
  auto metVec             = event->getMET();
  double met              = metVec.Et();
  double metSig           = event->getMETSignificance();

  auto combiLeptons = combiElectrons + combiMuons;
  int nCombiLeptons = combiLeptons.size();

  auto baselineLeptons = baselineElectrons + baselineMuons;
  int nBaselineLeptons = baselineLeptons.size();

  //====================================================================================================
  // overlap removal
  auto radiusCalcLepton = [] (const AnalysisObject& lepton, const AnalysisObject&) { return std::min(0.4, 0.04 + 10./lepton.Pt()); };
  // 1. e-mu: shared ID track
  baselineElectrons = overlapRemoval(baselineElectrons, baselineMuons, 0.01);
  // 2a. j-e: deltaR < 0.2
  baselineJets      = overlapRemoval(baselineJets, baselineElectrons, 0.2);
  // 2b. j-mu: deltaR < 0.2 but only LessThanThreeTrack
  baselineJets      = overlapRemoval(baselineJets, baselineMuons, 0.2, LessThan3Tracks);
  // 3. e/mu-j: (62.5:dR=0.2, 27.8: dR=0.4) variable radius: no leptons removed below/above 28/60
  baselineElectrons = overlapRemoval(baselineElectrons, baselineJets, radiusCalcLepton);
  baselineMuons     = overlapRemoval(baselineMuons, baselineJets, radiusCalcLepton);

  //====================================================================================================

  baselineLeptons  = baselineElectrons + baselineMuons;
  nBaselineLeptons = baselineLeptons.size();

  //====================================================================================================
  // signal objects
  auto signalElectrons = filterObjects(baselineElectrons, 10, 2.47, EMediumLH | ED0Sigma5 | EZ05mm | EIsoFCTight);
  auto signalMuons     = filterObjects(baselineMuons, 10, 2.5, MuD0Sigma3 | MuZ05mm | MuIsoFCTightFR);
  auto signalJets      = filterObjects(baselineJets, 20, 2.8, JVT120Jet);
  auto bjets           = filterObjects(signalJets, 20, 2.5, BTag85MV2c10);
  auto signalLeptons   = signalElectrons + signalMuons;

  //====================================================================================================

  int nLeptons   = signalLeptons.size();
  int nElectrons = signalElectrons.size();
  int nMuons     = signalMuons.size();
  int nJets      = signalJets.size();
  int nBjets     = bjets.size();

  // pT order the objects
  sortObjectsByPt(signalLeptons);
  sortObjectsByPt(signalJets);

  //====================================================================================================
  // define preselection cuts

  int cutflowcounter = 0;
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);

  if ( nCombiLeptons != 3 )             return;
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);

  if ( nBaselineLeptons != 3 )          return;
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);

  if ( nLeptons != 3 )                  return
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);
  
  if ( (signalLeptons[0].Pt()<25) )     return;
  if ( (signalLeptons[1].Pt()<20) )     return;
  if ( (signalLeptons[2].Pt()<10) )     return;
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);

  if ( nBjets>0 )                       return;
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);
  
  if ( met<50 )                         return;
  fill("cutflow", cutflowcounter++);
  fill("cutflow_DFOS", cutflowcounter++);
  
  accept("presel");

  //====================================================================================================

  // calculate variables for SR definitions
  int hasSFOS = 0.;
  double mll3lZ = 9999999.;
  double mll3ldZ = 9999999.;
  double mTWZ = 0.;
  float Zmass = 91.188;

  for ( signed int i = 0; i < nLeptons; i++ ) {
    for ( signed int j = 0; j < nLeptons; j++ ) {
      for ( signed int k = 0; k < nLeptons; k++ ) {
        if ( i == j || j == k || i == k ) {
          continue;
        }
        if ( (signalLeptons[i].charge() != signalLeptons[j].charge()) &&
             (signalLeptons[i].type() == signalLeptons[j].type()) ) {
          hasSFOS = 1;
          float tmp_mll = ( signalLeptons[i] + signalLeptons[j] ).M();
          if ( fabs(tmp_mll - Zmass) < mll3ldZ ) {
            mll3ldZ = fabs(tmp_mll - Zmass);
            mll3lZ = tmp_mll;
            mTWZ = calcMT(signalLeptons[k], metVec);
          }
        }
      }
    }
  }

  double HT = sumObjectsPt(signalJets);
  double HTlep = sumObjectsPt(signalLeptons);
  double mlll = (signalLeptons[0]+signalLeptons[1]+signalLeptons[2]).M();

  bool passDFOS = false;
  double dRnear = 0;
  for (signed int i = 0; i < nLeptons; i++) {
    for (signed int j = 0; j < nLeptons; j++) {
      if ((signalLeptons[i].charge() == signalLeptons[j].charge()) &&
          (signalLeptons[i].type() == signalLeptons[j].type()) && (i != j)) {
        for (signed int k = 0; k < nLeptons; k++) {
          if (k != i && k != j &&
              (signalLeptons[k].charge() != signalLeptons[j].charge()) &&
              (signalLeptons[k].type() != signalLeptons[j].type())) {
            passDFOS = true;

            if (fabs(signalLeptons[k].DeltaPhi(signalLeptons[i])) <= 
                fabs(signalLeptons[k].DeltaPhi(signalLeptons[j]))) {
              dRnear = signalLeptons[i].DeltaR(signalLeptons[k]);
            } else {
              dRnear = signalLeptons[j].DeltaR(signalLeptons[k]);
            }
          }
        }
      }
    }
  }

  //====================================================================================================
  // Regions

  bool SFOSOFFSSR1(false), SFOSOFFSSR2(false),  SFOSOFFSSR3(false), SFOSOFFSSR4(false), SFOSOFFSSR5(false),
    SFOSOFFSSR6(false), SFOSOFFSSR7(false), SFOSOFFSSR8(false), SFOSOFFSSR9(false), SFOSOFFSSR10(false),
    SFOSOFFSSR11(false), SFOSOFFSSR12(false), SFOSOFFSSR13(false), SFOSOFFSSR14(false), SFOSOFFSSR15(false),
    SFOSOFFSSR16(false), SFOSOFFSSR17(false), SFOSOFFSSR18(false), SFOSOFFSSR19(false);

  bool SFOSONSSR1(false), SFOSONSSR2(false), SFOSONSSR3(false), SFOSONSSR4(false), SFOSONSSR5(false),
    SFOSONSSR6(false), SFOSONSSR7(false), SFOSONSSR8(false), SFOSONSSR9(false), SFOSONSSR10(false),
    SFOSONSSR11(false), SFOSONSSR12(false), SFOSONSSR13(false), SFOSONSSR14(false), SFOSONSSR15(false),
    SFOSONSSR16(false), SFOSONSSR17(false), SFOSONSSR18(false), SFOSONSSR19(false), SFOSONSSR20(false);

  bool DFOSSR1(false), DFOSSR2(false);

  bool VR1(false), VR2(false), VR3(false);

  bool presel_DFOS = true;
  // DFOS SRs
  if   ( passDFOS == !1 ) presel_DFOS = false;
  if (presel_DFOS) fill("cutflow_DFOS", cutflowcounter++);
  if (  metSig<8 ) presel_DFOS = false;
  if (presel_DFOS) fill("cutflow_DFOS", cutflowcounter++);

  if ( passDFOS == 1 && signalLeptons[2].Pt()>15 && metSig>8 ) {
    if ( nJets==0 && dRnear<1.2                                            ) { DFOSSR1 = true; accept("DFOSSR1"); }
    if ( (nJets==1 || nJets==2) && dRnear<1.0 && signalLeptons[2].Pt()>20  ) { DFOSSR2 = true; accept("DFOSSR2"); }
  }

  bool presel_SFOS = true;

  if   ( hasSFOS == !1 ) presel_SFOS = false;
  if (presel_SFOS) fill("cutflow", cutflowcounter++);
  if   (  fabs(mlll-Zmass)<15 ) presel_SFOS = false;
  if (presel_SFOS) fill("cutflow", cutflowcounter++);

  /// Nominal SFOS SRs and VRs
  if (hasSFOS == 1 && fabs(mlll-Zmass)>15 ){

    // Offshell
    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ <100              && met >50   && met <100  ) { SFOSOFFSSR1 = true; accept("SFOSOFFSSR1"); }
    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ <100              && met >100  && met <150  ) { SFOSOFFSSR2 = true; accept("SFOSOFFSSR2"); }
    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ <100              && met >150               ) { SFOSOFFSSR3 = true; accept("SFOSOFFSSR3"); }

    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ >100 && mTWZ <160 && met >50   && met <100  ) { SFOSOFFSSR4 = true; accept("SFOSOFFSSR4"); }
    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ >100 && mTWZ <160 && met >100               ) { SFOSOFFSSR5 = true; accept("SFOSOFFSSR5"); }
    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ >160              && met >50   && met <100  ) { SFOSOFFSSR6 = true; accept("SFOSOFFSSR6"); }
    if ( mll3lZ>12   &&   mll3lZ<75    && nJets==0  && mTWZ >160              && met >100               ) { SFOSOFFSSR7 = true; accept("SFOSOFFSSR7"); }

    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ <50               && met >50 && met <100    ) { SFOSOFFSSR8 = true; accept("SFOSOFFSSR8"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ >50  && mTWZ <100 && met >50 && met <100    ) { SFOSOFFSSR9 = true; accept("SFOSOFFSSR9"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ <100              && met >100 && met <150   ) { SFOSOFFSSR10 = true; accept("SFOSOFFSSR10"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ <100              && met >150               ) { SFOSOFFSSR11 = true; accept("SFOSOFFSSR11"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ >100 && mTWZ <160 && met >50 && met <100    ) { SFOSOFFSSR12 = true; accept("SFOSOFFSSR12"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ >100 && mTWZ <160 && met >100 && met <150   ) { SFOSOFFSSR13 = true; accept("SFOSOFFSSR13"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ >100 && mTWZ <160 && met >150               ) { SFOSOFFSSR14 = true; accept("SFOSOFFSSR14"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ >160              && met >50 && met <150    ) { SFOSOFFSSR15 = true; accept("SFOSOFFSSR15"); }
    if ( mll3lZ>12 && mll3lZ<75    && nJets>0  && HT <200 && mTWZ >160              && met >150               ) { SFOSOFFSSR16 = true; accept("SFOSOFFSSR16"); }

    if ( mll3lZ>105                && nJets==0 && mTWZ >100              && met >50 && met <100    ) { SFOSOFFSSR17 = true; accept("SFOSOFFSSR17"); }
    if ( mll3lZ>105                && nJets==0 && mTWZ >100              && met >100 && met <200   ) { SFOSOFFSSR18 = true; accept("SFOSOFFSSR18"); }
    if ( mll3lZ>105                && nJets==0 && mTWZ >100              && met >200               ) { SFOSOFFSSR19 = true; accept("SFOSOFFSSR19"); }

    // Onshell
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >100 && mTWZ <160 && met >50  && met <100   ) { SFOSONSSR1 = true; accept("SFOSONSSR1"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >100 && mTWZ <160 && met >100 && met <150   ) { SFOSONSSR2 = true; accept("SFOSONSSR2"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >100 && mTWZ <160 && met >150 && met <200   ) { SFOSONSSR3 = true; accept("SFOSONSSR3"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >100 && mTWZ <160 && met >200               ) { SFOSONSSR4 = true; accept("SFOSONSSR4"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >160              && met >50 && met <150    ) { SFOSONSSR5 = true; accept("SFOSONSSR5"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >160              && met >150 && met <200   ) { SFOSONSSR6 = true; accept("SFOSONSSR6"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >160              && met >200 && met <350   ) { SFOSONSSR7 = true; accept("SFOSONSSR7"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets==0 &&  mTWZ >160              && met >350               ) { SFOSONSSR8 = true; accept("SFOSONSSR8"); }

    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >100 && mTWZ <160 && met >100 && met <150   ) { SFOSONSSR9 = true; accept("SFOSONSSR9"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >100 && mTWZ <160 && met >150 && met <250   ) { SFOSONSSR10 = true; accept("SFOSONSSR10"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >100 && mTWZ <160 && met >250 && met <300   ) { SFOSONSSR11 = true; accept("SFOSONSSR11"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >100 && mTWZ <160 && met >300               ) { SFOSONSSR12 = true; accept("SFOSONSSR12"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >160              && met >50 && met <150    ) { SFOSONSSR13 = true; accept("SFOSONSSR13"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >160              && met >150 && met <250   ) { SFOSONSSR14 = true; accept("SFOSONSSR14"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >160              && met >250 && met <400   ) { SFOSONSSR15 = true; accept("SFOSONSSR15"); }
    if ( mll3lZ>75 && mll3lZ<105   && nJets>0  && HT <200 &&  mTWZ >160              && met >400               ) { SFOSONSSR16 = true; accept("SFOSONSSR16"); }

    if ( mll3lZ>75 && mll3lZ<105  && nJets>0  && HT >200 &&  HTlep<350  &&  mTWZ >100    && met >150 && met <200   ) { SFOSONSSR17 = true; accept("SFOSONSSR17"); }
    if ( mll3lZ>75 && mll3lZ<105  && nJets>0  && HT >200 &&  HTlep<350  &&  mTWZ >100    && met >200 && met <300   ) { SFOSONSSR18 = true; accept("SFOSONSSR18"); }
    if ( mll3lZ>75 && mll3lZ<105  && nJets>0  && HT >200 &&  HTlep<350  &&  mTWZ >100    && met >300 && met <400   ) { SFOSONSSR19 = true; accept("SFOSONSSR19"); }
    if ( mll3lZ>75 && mll3lZ<105  && nJets>0  && HT >200 &&  HTlep<350  &&  mTWZ >100    && met >400               ) { SFOSONSSR20 = true; accept("SFOSONSSR20"); }

  }

  if ( hasSFOS ==1 && fabs(mlll-Zmass)>15 && met>100 && mll3lZ>75 && mll3lZ<105 && mTWZ<100 && mTWZ>20 ) {
    if ( nJets==0                   ) { VR1 = true; accept("VR1"); }
    if ( nJets>0   && HT<200        ) { VR2 = true; accept("VR2"); }
    if ( nJets>0   && HT>200        ) { VR3 = true; accept("VR3"); }
  }

  if ( hasSFOS ==1 && fabs(mlll-Zmass)>15 && met<100 && mll3lZ>75 && mll3lZ<105 && mTWZ<100 && mTWZ>20 ) {
    if ( nJets==0                   ) { accept("CR1"); }
    if ( nJets>0   && HT<200        ) { accept("CR2"); }
    if ( nJets>0   && HT>200        ) { accept("CR3"); }
  }

  //====================================================================================================
  // inclusive SR/VR

  bool SFOSOFFSSR(false), SFOSONSSR(false), SFOSOFFSSR_lowMll(false), SFOSOFFSSR_highMll(false);

  if ( SFOSOFFSSR1 || SFOSOFFSSR2 || SFOSOFFSSR3 || SFOSOFFSSR4 || SFOSOFFSSR5 ||
       SFOSOFFSSR6 || SFOSOFFSSR7 || SFOSOFFSSR8 || SFOSOFFSSR9 || SFOSOFFSSR10 ||
       SFOSOFFSSR11 || SFOSOFFSSR12 || SFOSOFFSSR13 || SFOSOFFSSR14 || SFOSOFFSSR15 ||
       SFOSOFFSSR16 || SFOSOFFSSR17 || SFOSOFFSSR18 || SFOSOFFSSR19 ) { SFOSOFFSSR=true; accept("SFOSOFFSSR"); }

  if ( SFOSONSSR1 || SFOSONSSR2 || SFOSONSSR3 || SFOSONSSR4 || SFOSONSSR5 ||
       SFOSONSSR6 || SFOSONSSR7 || SFOSONSSR8 || SFOSONSSR9 || SFOSONSSR10 ||
       SFOSONSSR11 || SFOSONSSR12 || SFOSONSSR13 || SFOSONSSR14 || SFOSONSSR15 ||
       SFOSONSSR16 || SFOSONSSR17 || SFOSONSSR18 || SFOSONSSR19 || SFOSONSSR20 ) { SFOSONSSR=true; accept("SFOSONSSR"); }

  if ( SFOSOFFSSR1 || SFOSOFFSSR2 || SFOSOFFSSR3 || SFOSOFFSSR4 || SFOSOFFSSR5 ||
       SFOSOFFSSR6 || SFOSOFFSSR7 || SFOSOFFSSR8 || SFOSOFFSSR9 || SFOSOFFSSR10 ||
       SFOSOFFSSR11 || SFOSOFFSSR12 || SFOSOFFSSR13 || SFOSOFFSSR14 || SFOSOFFSSR15 ||
       SFOSOFFSSR16 ) { SFOSOFFSSR_lowMll=true; accept("SFOSOFFSSR_lowMll"); }

  if ( SFOSOFFSSR17 || SFOSOFFSSR18 || SFOSOFFSSR19 ) { SFOSOFFSSR_highMll=true; accept("SFOSOFFSSR_highMll"); }

  if ( SFOSOFFSSR && nJets==0 ) accept("SFOSOFFSSR_j0");
  if ( SFOSONSSR && nJets==0 ) accept("SFOSONSSR_j0");
  if ( SFOSOFFSSR_lowMll && nJets==0 ) accept("SFOSOFFSSR_lowMll_j0");
  if ( SFOSOFFSSR_highMll && nJets==0 ) accept("SFOSOFFSSR_highMll_j0");

  if ( SFOSOFFSSR && nJets>0 ) accept("SFOSOFFSSR_j1inc");
  if ( SFOSONSSR && nJets>0 ) accept("SFOSONSSR_j1inc");
  if ( SFOSOFFSSR_lowMll && nJets>0 ) accept("SFOSOFFSSR_lowMll_j1inc");
  if ( SFOSOFFSSR_highMll && nJets>0 ) accept("SFOSOFFSSR_highMll_j1inc");

  if ( DFOSSR1 || DFOSSR2 ) accept("DFOSSR");
  if ( passDFOS == 1 && signalLeptons[2].Pt()>15 ) accept("DFOSSR_loose");

  if ( VR1 || VR2 || VR3 )  accept("VR");

  //====================================================================================================
  // Fill optional ntuple
  ntupVar("hasSFOS",hasSFOS);
  ntupVar("passDFOS",passDFOS);
  ntupVar("met",met);
  ntupVar("mll3lZ",mll3lZ);
  ntupVar("mlll",mlll);
  ntupVar("mTWZ",mTWZ);
  ntupVar("nJets",nJets);
  ntupVar("nBjets",nBjets);
  ntupVar("nLeptons",nLeptons);
  ntupVar("nElectrons",nElectrons);
  ntupVar("nMuons",nMuons);
  ntupVar("HT",HT);
  ntupVar("HTlep",HTlep);
  ntupVar("dRnear",dRnear);
  ntupVar("metSig",metSig);
  ntupVar("susyProcess",event->getSUSYChannel());
  ntupVar("mcDSID",event->getMCNumber());
  ntupVar("mcWeights",event->getMCWeights());

  return;
}
