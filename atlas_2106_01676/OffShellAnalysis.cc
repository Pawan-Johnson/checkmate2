#include "SimpleAnalysisFramework/AnalysisClass.h"

DefineAnalysis(EwkThreeLeptonOffshell2018)

void EwkThreeLeptonOffshell2018::Init() {
   // add regions
   addRegions({"3L","Presel","Presel_offshell"});
   addRegions({"Presellow_0j", "Presellow_nj", "Preselhigh_0j", "Preselhigh_nj"});
   addRegions({"SRlow_0j", "SRlow_nj", "SRhigh_0j", "SRhigh_nj"});
   addRegions({"CRWZ_0j", "CRWZ_nj", "VRWZ_0j", "VRWZ_nj", "VRWZ_nj_lowmll", "VRttbar"});
   addRegions({"SRlow_0jb", "SRlow_0jc", "SRlow_0jd", "SRlow_0je", "SRlow_0jf1", "SRlow_0jf2", "SRlow_0jg1", "SRlow_0jg2"});
   addRegions({"SRlow_njb", "SRlow_njc", "SRlow_njd", "SRlow_nje", "SRlow_njf1", "SRlow_njf2", "SRlow_njg1", "SRlow_njg2"});
   addRegions({"SRhigh_0jb", "SRhigh_0jc", "SRhigh_0jd", "SRhigh_0je", "SRhigh_0jf1", "SRhigh_0jf2", "SRhigh_0jg1", "SRhigh_0jg2"});
   addRegions({"SRhigh_nja", "SRhigh_njb", "SRhigh_njc", "SRhigh_njd", "SRhigh_nje", "SRhigh_njf", "SRhigh_njg"});
   //
   // Book 1/2D histograms
   addHistogram("cutflow",20,0,20);
}

void EwkThreeLeptonOffshell2018::ProcessEvent(AnalysisEvent *event) {
   //====================================================================================================
   // get baseline objects
   // cfr. STconfig: https://gitlab.cern.ch/atlas-phys-susy-wg/AnalysisSUSYToolsConfigurations/-/blob/master/ANA-SUSY-2019-09_Offshell.conf
   auto baselineElectrons  = event->getElectrons(4.5, 2.47, ELooseBLLH);
   auto baselineMuons      = event->getMuons(3.0, 2.5, MuMedium | MuNotCosmic | MuQoPSignificance);
   auto baselineJets       = event->getJets(20., 4.5);
   auto metVec             = event->getMET();
   double met              = metVec.Et();
   double metSig           = event->getMETSignificance();

   //====================================================================================================
   auto baselineleptons = baselineElectrons + baselineMuons;
   int nbaselineleptons = baselineleptons.size();

   // start cutflow counter
   int cutflowcounter = 0;
   fill("cutflow", cutflowcounter++);
   bool DEBUG=false;

   // 3 baseline leptons _before_ overlap removal
   if (nbaselineleptons != 3) return;
   fill("cutflow", cutflowcounter++);

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
   baselineleptons = baselineElectrons + baselineMuons;
   nbaselineleptons = baselineleptons.size();

   // 3 baseline leptons _after_ overlap removal
   if (nbaselineleptons != 3) return;
   fill("cutflow", cutflowcounter++);

   //====================================================================================================
   // signal objects
   auto electrons = filterObjects(baselineElectrons, 4.5, 2.47, EMediumLH | ED0Sigma5 | EZ05mm | EIsoGradient); // e: Gradient
   auto muons     = filterObjects(baselineMuons, 3.0, 2.5, MuD0Sigma3 | MuZ05mm | MuIsoFCLoose);                // m: FCLoose
   auto jets20    = filterObjects(baselineJets, 20, 2.8, JVT120Jet);
   auto jets      = filterObjects(baselineJets, 30, 2.8, JVT120Jet);
   auto bjets     = filterObjects(jets20, 20., 2.5, BTag85MV2c10);                                              // b: MV2c10, fc85
   auto leptons   = electrons + muons;

   //====================================================================================================
   int nleptons   = leptons.size();
   int nelectrons = electrons.size();
   int nmuons     = muons.size();
   int njets      = jets.size(); // use nJets_30, not 20
   int nbjets     = bjets.size(); // uses jets20

   // pT order the objects
   sortObjectsByPt(leptons);
   sortObjectsByPt(jets);

   if (DEBUG) {
      std::cout << "nleptons: " << nleptons << std::endl;
      for (auto lep : leptons) { std::cout << "lep: (type:" << lep.type() << "; id: (" << lep.id() << ") --> "; lep.Print(); }
      std::cout << "njets: " << njets << std::endl;
      for (auto jet : jets) { std::cout << "jet: "; jet.Print(); }
      std::cout << "met: "; metVec.Print();
      std::cout << "metSig: " << metSig << std::endl;
   }

   //====================================================================================================
   // Baseline cuts
   //
   // exactly 3 baseline and signal leptons (_after_ overlap removal)
   if (nleptons != 3) return;
   fill("cutflow", cutflowcounter++);
   accept("3L");

   //====================================================================================================
   //do mll assignment & mindR_OSSF + mindR_3L
   double mZ  = 91.1876;
   double mll = -999.;
   double mdiff = 1e6;
   double minmll = 1e6;
   double maxmll = -999;
   double mindR_OSSF = 1e6;
   double mindR_3L = 1e6;
   int nOSSF = 0;
   int i0, i1, i2;
   int iZ1(-1), iZ2(-1), iW(-1); // closest to mZ order
   int jZ1(-1), jZ2(-1), jW(-1); // minmll order

   // lepton permutations Z1, Z2, W
   std::string index="012";
   do {
      i0 = index[0]-'0'; i1 = index[1]-'0'; i2 = index[2]-'0';  // char2int

      // only Z1,Z2 with Z1>Z2
      if (i0<i1) continue;

      double idR = leptons[i0].DeltaR(leptons[i1]);
      // mindR_3L
      if (idR < mindR_3L) mindR_3L = idR;

      // OSSF pairs only
      if ( (leptons[i0].charge() == leptons[i1].charge()) || (leptons[i0].type() != leptons[i1].type()) ) continue;
      nOSSF++;

      double imll = (leptons[i0]+leptons[i1]).M();
      double imdiff = fabs(imll-mZ);

      // closest to mZ
      if (imdiff < mdiff) {
         mdiff = imdiff;
         mll = imll;
         iZ1 = i0; iZ2 = i1; iW = i2;
      }
      // minmll
      if (imll < minmll) {
         minmll = imll;
         jZ1 = i0; jZ2 = i1; jW = i2;
      }
      // maxmll
      if (imll > maxmll) maxmll = imll;
      // mindR_OSSF
      if (idR < mindR_OSSF) mindR_OSSF = idR;

   } while (std::next_permutation(index.begin(), index.end()));

   if (DEBUG) {
      std::cout << "mll: " << mll << std::endl;
      std::cout << "minmll: " << minmll << std::endl;
      std::cout << "maxmll: " << maxmll << std::endl;
   }

   // at least one OSSF pair
   if (nOSSF==0) return;
   fill("cutflow", cutflowcounter++);

   //====================================================================================================
   // at this point everything should be assigned!
   if (iZ1==-1 || iZ2==-1 || iW==-1 || jZ1==-1 || jZ2==-1 || jW==-1) {
      std::cout << "Z(ll)W(lv) lepton assignment failed!" << std::endl;
      for (int ilep=0; ilep<(int)leptons.size(); ilep++) {
         for (int jlep=ilep+1; jlep<3; jlep++) {
            printf("ilep %4d jlep %4d | ich %4d jch %4d | itype %4d jtype %4d | mll %6.2f\n", ilep, jlep, leptons[ilep].charge(), leptons[jlep].charge(), leptons[ilep].type(), leptons[jlep].type(), (leptons[ilep]+leptons[jlep]).M());
         }
      }
      return;
   }
   fill("cutflow", cutflowcounter++);

   //====================================================================================================
   // More baseline cuts
   // resonance vetoes in low met regions
   if ( met<200. && ((minmll>3.0 && minmll<3.2) || (minmll>9.&& minmll<12.)) ) return;
   fill("cutflow", cutflowcounter++);
   // MET (high met regions) or pT (other)
   if ( ! ( met>200. || leptons[2].Pt()>10. ) ) return;
   fill("cutflow", cutflowcounter++);
   // dR cleaning
   if (mindR_3L < 0.4) return;
   fill("cutflow", cutflowcounter++);
   // PLV on 3rd lepton
   if (leptons[2].type()==ELECTRON && !leptons[2].pass(EIsoPLVTight)) return;
   if (leptons[2].type()==MUON && !leptons[2].pass(MuIsoPLVTight)) return;
   fill("cutflow", cutflowcounter++);

   double m3l(-1), mT_mZ(-1), mT_minmll(-1), mT2_100_minmll(-1), dRWMET(-1), pT3L(-1), pT3LoverMET(-1), mW_pZB(-1);
   int lfaketype(-1);
   bool anticonv(false), efake(0), mfake(0);

   //====================================================================================================
   // define additional kinematic variables
   m3l = (leptons[0]+leptons[1]+leptons[2]).M();
   mT_mZ = calcMT(leptons[iW], metVec);
   mT_minmll = calcMT(leptons[jW], metVec);
   mT2_100_minmll = calcAMT2(leptons[jZ1]+leptons[jZ2], leptons[jW], metVec, 100, 100);
   dRWMET = leptons[jW].DeltaR(metVec);
   pT3L = (leptons[0]+leptons[1]+leptons[2]).Pt();
   pT3LoverMET = pT3L / met;
   mW_pZB = 0.;
   TLorentzVector wnu_pzb;
   wnu_pzb.SetXYZM(metVec.Px(), metVec.Py(), (leptons[jZ1]+leptons[jZ2]).Pz() - leptons[jW].Pz(), 0.); // nu_pz assuming W & Z are balanced in pZ.
   mW_pZB = (leptons[jW] + wnu_pzb).M();
   efake = (leptons[2].type()==ELECTRON)?1:0;
   mfake = (leptons[2].type()==MUON)?1:0;
   lfaketype = (efake)?1:2;
   if ( (efake==true && fabs(m3l-mZ)>20. && mindR_OSSF<2.4 && mindR_OSSF>0.6) || (mfake==true) ) anticonv = true;
   
   //====================================================================================================
   // Baseline regions (b-veto added only after CRs/VRs)
   // Preselected: 3L, >=1 OSSF, SUSY2/16 pT/MET sel, dR cleaning (0.4), PLV iso 3rd lepton
   fill("cutflow", cutflowcounter++); 
   accept("Presel");
   // Preselected: 3L, >=1 OSSF, SUSY2/16 pT/MET sel, dR cleaning (0.4), PLV iso 3rd lepton && offshell
   if (mll < 75. && maxmll < 75.) {
      fill("cutflow", cutflowcounter++); 
      accept("Presel_offshell");
   }

   //====================================================================================================
   // Control & validation regions
   // CRWZ
   if      (mll>81. && mll<101. && njets==0 && met<50. && leptons[2].Pt()>10. && mT_mZ>50. && nbjets==0) accept("CRWZ_0j");
   else if (mll>81. && mll<101. && njets>0  && met<50. && leptons[2].Pt()>10. && mT_mZ>50. && nbjets==0) accept("CRWZ_nj");

   if (mll < 75. && maxmll < 75.) {
      // VRWZ
      if (metSig>1.5 && nbjets==0) {
         if      (minmll>12. && minmll<75. && njets==0 && met<50. && leptons[2].Pt()>10. && mT_minmll>60. && mT_minmll<90. && anticonv && mW_pZB>75. && dRWMET>2.6) accept("VRWZ_0j");
         else if (minmll>12. && minmll<75. && njets>0  && met<80. && leptons[2].Pt()>10. && mT_minmll>60. && mT_minmll<90. && anticonv && pT3LoverMET>0.5         ) accept("VRWZ_nj");
         else if (minmll>1.  && minmll<9.  && njets>0  && met>80. &&                        mT_minmll>30.                              && pT3LoverMET>0.3         ) accept("VRWZ_nj_lowmll");
      }

      // VRtt
      if (minmll>1. && minmll<75. && met>50. && nbjets>0) accept("VRttbar");
   }

   //====================================================================================================
   // b-jet veto
   if (nbjets>0) return;
   fill("cutflow", cutflowcounter++); 

   //====================================================================================================
   // Baseline regions
   bool presellow_0j(false), presellow_nj(false), preselhigh_0j(false), preselhigh_nj(false);

   bool minmllgt12 = (minmll>12.?1:0);

   if ( njets==0 && maxmll < 75.) { // 0j
      if ( met<50. && metSig>1.5 && anticonv && minmllgt12 && minmll<75. ) { presellow_0j = true; accept("Presellow_0j"); }
      else if ( met>50. && metSig>3.0 && minmllgt12 && minmll<75. )        { preselhigh_0j = true; accept("Preselhigh_0j"); }
   }
   else if (maxmll < 75.) {         // nj
      bool softlep(true);
      for (int ilep=0; ilep<nleptons; ilep++) { if ( (leptons[ilep].type()==ELECTRON && leptons[ilep].Pt()<4.5) || (leptons[ilep].type()==MUON && leptons[ilep].Pt()<3.0) ) { softlep = false; break; } }
      //
      if ( met<200. && metSig>3.0 && anticonv && minmllgt12 && minmll<75. ) { presellow_nj = true; accept("Presellow_nj"); }
      else if ( met>200. && metSig>3.0 && softlep && minmll<75. )           { preselhigh_nj = true; accept("Preselhigh_nj"); }
   }

   if (!(presellow_0j||presellow_nj||preselhigh_0j||preselhigh_nj)) return;
   fill("cutflow", cutflowcounter++); 

   //====================================================================================================
   // Signal regions
   //
   bool SRlow_0jb(false), SRlow_0jc(false), SRlow_0jd(false), SRlow_0je(false), SRlow_0jf1(false), SRlow_0jf2(false), SRlow_0jg1(false), SRlow_0jg2(false);
   bool SRlow_njb(false), SRlow_njc(false), SRlow_njd(false), SRlow_nje(false), SRlow_njf1(false), SRlow_njf2(false), SRlow_njg1(false), SRlow_njg2(false);
   bool                    SRhigh_0jb(false), SRhigh_0jc(false), SRhigh_0jd(false), SRhigh_0je(false), SRhigh_0jf1(false), SRhigh_0jf2(false), SRhigh_0jg1(false), SRhigh_0jg2(false);
   bool SRhigh_nja(false), SRhigh_njb(false), SRhigh_njc(false), SRhigh_njd(false), SRhigh_nje(false), SRhigh_njf(false),                      SRhigh_njg(false);
   
   // SRlow_0j
   if (presellow_0j) {
      if (maxmll<60. && minmll<40. && leptons[2].Pt()>10. && mT_minmll<60. && pT3LoverMET<1.3) {
         if      (minmll>12. && minmll<15. && mT_minmll<50. && pT3LoverMET<1.1 && mT2_100_minmll<115. && mindR_3L<1.6) { SRlow_0jb = true; accept("SRlow_0jb"); }
         else if (minmll>15. && minmll<20. && mT_minmll<50. && pT3LoverMET<1.1 && mT2_100_minmll<120. && mindR_3L<1.6) { SRlow_0jc = true; accept("SRlow_0jc"); }
         else if (minmll>20. && minmll<30. && mT_minmll<50. && pT3LoverMET<1.1 && mT2_100_minmll<130. && mindR_3L<1.6) { SRlow_0jd = true; accept("SRlow_0jd"); }
         else if (minmll>30.)                                                                                          { SRlow_0je = true; accept("SRlow_0je"); }
      }
      else if (minmll>40. && minmll<75. && leptons[2].Pt()>15. && (mT_minmll<60.||mT_minmll>90.) && pT3LoverMET<1.4 && m3l>100) {
         if (minmll<60.) {
            if (mT_minmll<60.)      { SRlow_0jf1 = true; accept("SRlow_0jf1"); }
            else if (mT_minmll>90.) { SRlow_0jf2 = true; accept("SRlow_0jf2"); }
         }
         else if (minmll>60.) {
            if (mT_minmll<60.)      { SRlow_0jg1 = true; accept("SRlow_0jg1"); }
            else if (mT_minmll>90.) { SRlow_0jg2 = true; accept("SRlow_0jg2"); }
         }
      }
   }
   // SRlow_nj
   else if (presellow_nj) {
      if (maxmll<60. && minmll<40. && leptons[2].Pt()>10. && mT_minmll<60. && pT3LoverMET<1.0) {
         if      (minmll>12. && minmll<15. && mT_minmll<50. && mT2_100_minmll<115. && mindR_3L<1.6) { SRlow_njb = true; accept("SRlow_njb"); }
         else if (minmll>15. && minmll<20. && mT_minmll<50. && mT2_100_minmll<120. && mindR_3L<1.6) { SRlow_njc = true; accept("SRlow_njc"); }
         else if (minmll>20. && minmll<30. && mT_minmll<50. && mT2_100_minmll<130. && mindR_3L<1.6) { SRlow_njd = true; accept("SRlow_njd"); }
         else if (minmll>30.)                                                                       { SRlow_nje = true; accept("SRlow_nje"); }
      }
      else if (minmll>40. && leptons[2].Pt()>15. && (mT_minmll<60.||mT_minmll>90.) && pT3LoverMET<1.2) {
         if (minmll<60.) {
            if (mT_minmll<60.)      { SRlow_njf1 = true; accept("SRlow_njf1"); }
            else if (mT_minmll>90.) { SRlow_njf2 = true; accept("SRlow_njf2"); }
         }
         else if (minmll>60. && minmll<75.) {
            if (mT_minmll<60.)      { SRlow_njg1 = true; accept("SRlow_njg1"); }
            else if (mT_minmll>90.) { SRlow_njg2 = true; accept("SRlow_njg2"); }
         }
      }
   }
   // SRhigh_0j
   else if (preselhigh_0j && leptons[0].Pt()>25. && leptons[1].Pt()>15. && leptons[2].Pt()>10.) {
      if      (minmll>12. && minmll<15. && mT2_100_minmll<115. && mT_minmll<50.) { SRhigh_0jb  = true; accept("SRhigh_0jb"); }
      else if (minmll>15. && minmll<20. && mT2_100_minmll<120. && mT_minmll<50.) { SRhigh_0jc  = true; accept("SRhigh_0jc"); }
      else if (minmll>20. && minmll<30. && mT2_100_minmll<130. && mT_minmll<60.) { SRhigh_0jd  = true; accept("SRhigh_0jd"); }
      else if (minmll>30. && minmll<40. && mT2_100_minmll<140. && mT_minmll<60.) { SRhigh_0je  = true; accept("SRhigh_0je"); }
      else if (minmll>40. && minmll<60. && mT2_100_minmll<160. && mT_minmll<70.) { SRhigh_0jf1 = true; accept("SRhigh_0jf1"); }
      else if (minmll>40. && minmll<60. && mT2_100_minmll<160. && mT_minmll>90.) { SRhigh_0jf2 = true; accept("SRhigh_0jf2"); }
      else if (minmll>60. && minmll<75. && mT2_100_minmll<175. && mT_minmll<70.) { SRhigh_0jg1 = true; accept("SRhigh_0jg1"); }
      else if (minmll>60. && minmll<75. && mT2_100_minmll<175. && mT_minmll>90.) { SRhigh_0jg2 = true; accept("SRhigh_0jg2"); }
   }
   // SRhigh_nj
   else if (preselhigh_nj && pT3LoverMET<1.0) {
      if      (minmll>1.  && minmll<12. && mT2_100_minmll<112. && pT3LoverMET<0.2) { SRhigh_nja = true; accept("SRhigh_nja"); }
      else if (minmll>12. && minmll<15. && mT2_100_minmll<115. && pT3LoverMET<0.2) { SRhigh_njb = true; accept("SRhigh_njb"); }
      else if (minmll>15. && minmll<20. && mT2_100_minmll<120. && pT3LoverMET<0.3) { SRhigh_njc = true; accept("SRhigh_njc"); }
      else if (minmll>20. && minmll<30. && mT2_100_minmll<130. && pT3LoverMET<0.3) { SRhigh_njd = true; accept("SRhigh_njd"); }
      else if (minmll>30. && minmll<40. && mT2_100_minmll<140. && pT3LoverMET<0.3) { SRhigh_nje = true; accept("SRhigh_nje"); }
      else if (minmll>40. && minmll<60. && mT2_100_minmll<160. && pT3LoverMET<1.0) { SRhigh_njf = true; accept("SRhigh_njf"); }
      else if (minmll>60. && minmll<75. && mT2_100_minmll<175. && pT3LoverMET<1.0) { SRhigh_njg = true; accept("SRhigh_njg"); }
   }

   // grouped
   if (            SRlow_0jb ||SRlow_0jc ||SRlow_0jd ||SRlow_0je ||SRlow_0jf1 ||SRlow_0jf2 ||SRlow_0jg1 ||SRlow_0jg2)  accept("SRlow_0j");
   if (            SRlow_njb ||SRlow_njc ||SRlow_njd ||SRlow_nje ||SRlow_njf1 ||SRlow_njf2 ||SRlow_njg1 ||SRlow_njg2)  accept("SRlow_nj");
   if (            SRhigh_0jb||SRhigh_0jc||SRhigh_0jd||SRhigh_0je||SRhigh_0jf1||SRhigh_0jf2||SRhigh_0jg1||SRhigh_0jg2) accept("SRhigh_0j");
   if (SRhigh_nja||SRhigh_njb||SRhigh_njc||SRhigh_njd||SRhigh_nje||SRhigh_njf              ||SRhigh_njg              ) accept("SRhigh_nj");

   //====================================================================================================
   // Fill optional ntuple
   //
   ntupVar("leptons",      leptons);
   ntupVar("jets",         jets);
   ntupVar("met",          met);
   ntupVar("metSig",       metSig);
   //
   ntupVar("nleptons",     nleptons);
   ntupVar("nelectrons",   nelectrons);
   ntupVar("nmuons",       nmuons);
   ntupVar("njets",        njets);
   ntupVar("nbjets",       nbjets);
   ntupVar("nOSSF",        nOSSF);
   //
   std::vector<int> idx_mZ     = {iZ1, iZ2, iW};
   std::vector<int> idx_minmll = {jZ1, jZ2, jW};
   ntupVar("idx_mZ",       idx_mZ);
   ntupVar("idx_minmll",   idx_minmll);
   //
   ntupVar("mll",          mll);
   ntupVar("minmll",       minmll);
   ntupVar("maxmll",       maxmll);
   ntupVar("mT_mZ",        mT_mZ);
   ntupVar("mT_minmll",    mT_minmll);
   ntupVar("mT2_100_minmll", mT2_100_minmll);
   ntupVar("mindR_OSSF",   mindR_OSSF);
   ntupVar("mindR_3L",     mindR_3L);
   ntupVar("m3l",          m3l );
   ntupVar("pT3L",         pT3L);
   ntupVar("pT3LoverMET",  pT3LoverMET);
   ntupVar("mW_pZB",       mW_pZB);
   //
   std::vector<int> isPLVTight(3,0);
   for (int ilep=0; ilep<nleptons; ilep++) {
      if (leptons[ilep].type()==ELECTRON) isPLVTight[ilep] = leptons[ilep].pass(EIsoPLVTight);
      else if (leptons[ilep].type()==MUON) isPLVTight[ilep] = leptons[ilep].pass(MuIsoPLVTight);
      else std::cout << "lepton type isn't 0 or 1" << std::endl;
   }
   ntupVar("isPLVTight",     isPLVTight);
   ntupVar("l3type",         lfaketype);
   //
   ntupVar("SusyProcess",  event->getSUSYChannel());
   ntupVar("DSID",         event->getMCNumber());
   ntupVar("mcweight",     event->getMCWeights());

   //====================================================================================================
   return;
}

