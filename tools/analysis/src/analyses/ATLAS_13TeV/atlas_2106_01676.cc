#include "atlas_2106_01676.h"
// AUTHOR: K. Rolbiecki
//  EMAIL: krolb@fuw.edu.pl
void Atlas_2106_01676::initialize()
{
  setAnalysisName("atlas_2106_01676");
  setInformation(""
                 "# electroweakinos, 3 leptons, WZ, WH, on+off-shell\n"
                 "");
  setLuminosity(139.0 * units::INVFB);
  bookSignalRegions("SR-WZ-01;SR-WZ-02;SR-WZ-03;SR-WZ-04;SR-WZ-05;SR-WZ-06;SR-WZ-07;SR-WZ-08;SR-WZ-09;SR-WZ-10;SR-WZ-11;SR-WZ-12;SR-WZ-13;SR-WZ-14;SR-WZ-15;SR-WZ-16;SR-WZ-17;SR-WZ-18;SR-WZ-19;SR-WZ-20;SR-Wh-SF-01;SR-Wh-SF-02;SR-Wh-SF-03;SR-Wh-SF-04;SR-Wh-SF-05;SR-Wh-SF-06;SR-Wh-SF-07;SR-Wh-SF-08;SR-Wh-SF-09;SR-Wh-SF-10;SR-Wh-SF-11;SR-Wh-SF-12;SR-Wh-SF-13;SR-Wh-SF-14;SR-Wh-SF-15;SR-Wh-SF-16;SR-Wh-SF-17;SR-Wh-SF-18;SR-Wh-SF-19;SR-Wh-DF-1;SR-Wh-DF-2;SR-offWZ-low-0jb;SR-offWZ-low-0jc;SR-offWZ-low-0jd;SR-offWZ-low-0je;SR-offWZ-low-0jf1;SR-offWZ-low-0jf2;SR-offWZ-low-0jg1;SR-offWZ-low-0jg2;SR-offWZ-low-njb;SR-offWZ-low-njc;SR-offWZ-low-njd;SR-offWZ-low-nje;SR-offWZ-low-njf1;SR-offWZ-low-njf2;SR-offWZ-low-njg1;SR-offWZ-low-njg2;SR-offWZ-high-0jb;SR-offWZ-high-0jc;SR-offWZ-high-0jd;SR-offWZ-high-0je;SR-offWZ-high-0jf1;SR-offWZ-high-0jf2;SR-offWZ-high-0jg1;SR-offWZ-high-0jg2;SR-offWZ-high-nja;SR-offWZ-high-njb;SR-offWZ-high-njc;SR-offWZ-high-njd;SR-offWZ-high-nje;SR-offWZ-high-njf;SR-offWZ-high-njg");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}

void Atlas_2106_01676::analyze()
{

  // ===============================================================================================================
  // Baseline Selection
  baselineElectrons = filterPhaseSpace(electronsLoose, 4.5, -2.47, 2.47);
  baselineMuons = filterPhaseSpace(muonsCombined, 3., -2.5, 2.5);
  // TODO: How to avoid puleup xinustion z0sint < 0.5
  auto baselineJets = filterPhaseSpace(jets, 30, -4.5, 4.5);
  // TODO: Photon with tight identification criteria. Less important photon conversion low probab.

  // ===============================================================================================================
  // Trigger Selection. ! Where should this go?
  std::sort(baselineElectrons.begin(), baselineElectrons.end(), sortByPTEl);
  std::sort(baselineMuons.begin(), baselineMuons.end(), sortByPTMu);
  double met = missingET->P4().Perp();
  trigger = false;
  if (baselineElectrons.size() > 1 and baselineElectrons[1]->PT > 18.)
    trigger = true;
  else if (baselineMuons.size() and baselineMuons[0]->PT > 27.3)
    trigger = true;
  else if (baselineMuons.size() > 1 and baselineMuons[1]->PT > 14.7)
    trigger = true;
  else if (baselineMuons.size() > 2 and baselineMuons[2]->PT > 6.5)
    trigger = true;
  else if (met > 200.)
    trigger = true; //! Check this
  // TODO: Intercation vertex with mini two tracks with pt > 0.5 GeV
  countCutflowEvent("00_all");

  // ==============================================================================================================
  // Ovelap Removal
  // TODO: Remove electrons and muons sharing ID Track
  baselineElectrons = overlapRemoval(baselineElectrons, baselineMuons, 0.1);
  baselineJets = overlapRemoval(baselineJets, baselineElectrons, 0.2, "y");
  baselineJets = overlapRemovalJetMuonAndTracks(baselineJets, baselineMuons, 0.4, 3);
  baselineElectrons = overlapRemoval(baselineElectrons, baselineJets, 0.4, "y");
  baselineMuons = overlapRemoval(baselineMuons, baselineJets, 0.4, "y");

  //! This was in the code. I am not sure where it is mentioned in the paper
  baselineElectrons = Isolate_leptons_with_inverse_track_isolation_cone(baselineElectrons, tracks, towers, 0.3, 10., 0.2, 0.15, 0.2, true);
  baselineMuons = Isolate_leptons_with_inverse_track_isolation_cone(baselineMuons, tracks, towers, 0.3, 10., 0.2, 0.06, 0.06, false);

  // ================================================================================================================
  // Missing ET
  missingET->addMuons(baselineMuons); //! Is it computed with the new baseline electrons and jets defined previously? Test both
  pTmiss = missingET->P4();
  double met = missingET->P4().Perp();
  // auto ETmisssig = missingpT / ( sigma_L*sqrt(1 - rho_LT*rho_LT)) //! Find def of sigmaL and rhoLT

  // ================================================================================================================
  // Signal Objects
  auto signalJets = filterPhaseSpace(baselineJets, 30, -2.8, 2.8); //! Loose Quality Criteria? JVT120?
  //! bjets selection
  auto signalElectrons = baselineElectrons; //! transverse impact parameter?
  auto signalMuons = baselineMuons;         //! transverse impact parameter?

  // TODO: FNP coorections to leptons

  // =================================================================================================================
  // Scale Eff between simulation snad data
  MCCorrections();

  countCutflowEvent("00_all");

  // =================================================================================================================
  // Preselction Varaibles

  // Leptons
  // ! Feature: Over load + operator to do this
  for (int i = 0; i < signalElectrons.size(); i++)
  {
    FinalStateObject *lep = newFinalStateObject(signalElectrons[i]);
    signalLeptons.push_back(lep);
  }
  for (int i = 0; i < signalMuons.size(); i++)
  {
    FinalStateObject *lep = newFinalStateObject(signalMuons[i]);
    signalLeptons.push_back(lep);
  }
  for (int i = 0; i < baselineElectrons.size(); i++)
  {
    FinalStateObject *lep = newFinalStateObject(baselineElectrons[i]);
    signalLeptons.push_back(lep);
  }
  for (int i = 0; i < baselineMuons.size(); i++)
  {
    FinalStateObject *lep = newFinalStateObject(baselineMuons[i]);
    signalLeptons.push_back(lep);
  }
  std::sort(signalLeptons.begin(), signalLeptons.end(), sortByPT);
  std::sort(baselineLeptons, begin(), baselineLeptons.end(), sortByPT);

  // b veto
  bveto = false;
  for (int i = 0; i < jets_off.size(); i++)
    if (fabs(jets_off[i]->Eta) < 2.5 && checkBTag(jets_off[i]))
    {
      bveto = true;
      break;
    }

  // if (!trigger) return;
  // countCutflowEvent("03_trigger");

  // SFOS and mll
  SFOS = false;
  if (signalElectrons.size() == 2 and signalElectrons[0]->Charge * signalElectrons[1]->Charge < 0)
  {
    SFOS = true;
    mll = (signalElectrons[0]->P4() + signalElectrons[1]->P4()).M();
  }
  else if (signalMuons.size() == 2 and signalMuons[0]->Charge * signalMuons[1]->Charge < 0)
  {  
    SFOS = true;
    mll = (signalMuons[0]->P4() + signalMuons[1]->P4()).M();
  }

  // HT
  HT = 0.;
  for (int i = 0; i < jets_off.size(); i++)
    HT += jets_off[i]->PT;

  // mt
  double mt = 0.;
  if (SFOS)
    mt = mT(signalLeptons[2]->P4(), pTmiss, 0.);

  // =================================================================================================================
  // Signal Regions
  if (preselection_onWZ(true))
  {
    if (jets_off.size() == 0)
    {
      countCutflowEvent("onWZ_09_0j");
      if (mt > 100. and mt < 160.)
      {
        countCutflowEvent("onWZ_10_mt<160");
        if (met < 100.)
          countSignalEvent("SR-WZ-01");
        else if (met < 150.)
          countSignalEvent("SR-WZ-02");
        else if (met < 200.)
          countSignalEvent("SR-WZ-03");
        else
          countSignalEvent("SR-WZ-04");
      }
      else if (mt > 160.)
      {
        countCutflowEvent("onWZ_11_mt>160");
        if (met < 150.)
          countSignalEvent("SR-WZ-05");
        else if (met < 200.)
          countSignalEvent("SR-WZ-06");
        else if (met < 350.)
          countSignalEvent("SR-WZ-07");
        else
          countSignalEvent("SR-WZ-08");
      }
    }
    else if (jets_off.size() > 0 and HT < 200.)
    {
      countCutflowEvent("onWZ_19_njHTlow");
      if (mt > 100. and mt < 160.)
      {
        countCutflowEvent("onWZ_20_mt<160");
        if (met > 100 and met < 150.)
          countSignalEvent("SR-WZ-09");
        else if (met > 150. and met < 250)
          countSignalEvent("SR-WZ-10");
        else if (met > 250. and met < 300)
          countSignalEvent("SR-WZ-11");
        else if (met > 300.)
          countSignalEvent("SR-WZ-12");
      }
      else if (mt > 160.)
      {
        countCutflowEvent("onWZ_21_mt>160");
        if (met < 150.)
          countSignalEvent("SR-WZ-13");
        else if (met < 250.)
          countSignalEvent("SR-WZ-14");
        else if (met < 400.)
          countSignalEvent("SR-WZ-15");
        else
          countSignalEvent("SR-WZ-16");
      }
    }
    else if (jets_off.size() > 0 and HT < 200.)
    {
      countCutflowEvent("onWZ_29_njHThigh");
      if (leptons[0]->PT + leptons[1]->PT + leptons[2]->PT < 350.)
      {
        countCutflowEvent("onWZ_30_njHTlep");
        if (mt > 100.)
        {
          countCutflowEvent("onWZ_31_mt>100");
          if (met > 150. and met < 200.)
            countSignalEvent("SR-WZ-17");
          else if (met > 200. and met < 300.)
            countSignalEvent("SR-WZ-18");
          else if (met > 300. and met < 400.)
            countSignalEvent("SR-WZ-19");
          else if (met > 400.)
            countSignalEvent("SR-WZ-20");
        }
      }
    }
  }

  return;
}

void Atlas_2106_01676::finalize()
{
  // Whatever should be done after the run goes here
}

bool Atlas_2106_01676::sortByPTEl(Electron *i, Electron *j) { return (i->PT > j->PT); }
bool Atlas_2106_01676::sortByPTMu(Muon *i, Muon *j) { return (i->PT > j->PT); }

// ! Validate this. How exactly are tracks and towers associated?
int Atlas_2106_01676::nTracksJet(Jet *jet, std::vector<Track *> tracks)
{
  int nTracks = 0;
  for (std::vector<Track *>::iterator it = tracks.begin(); it != tracks.end(); it++)
    for (int part = 0; part < jet->Particles.GetEntries(); part++)
      if (jet->Particles.At(part) == (*it)->Particle && (*it)->PT > 0.5)
        nTracks++;

  return nTracks;
}

// ! Validate this.
std::vector<Jet *> Atlas_2106_01676::overlapRemovalJetMuonAndTracks(std::vector<Jet *> cand_jets, std::vector<Muon *> cand_muons, double deltaR, int nTracks_max)
{

  std::vector<Jet *> passed;
  for (std::vector<Jet *>::iterator jet = cand_jets.begin(); jet != cand_jets.end(); jet++)
  {
    bool iso = true;
    for (std::vector<Muon *>::iterator mu = cand_muons.begin(); mu != cand_muons.end(); mu++)
      if ((*jet)->P4().DeltaR((*mu)->P4()) < deltaR and nTracksJet(*jet, tracks) > nTracks_max)
        iso = false;
    if (iso)
      passed.push_back(*jet);
  }

  return passed;
}

void Atlas_2106_01676::MCCorrections()
{
  vector<vector<double>> pt_edges;
  vector<double> effs;
  /*
  Doubts
  * Why loop over true particles instead of baseline?
  */
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (abs(true_particles[t]->PID == 11) and true_particles[t]->Status == 1)
    { // Is final state electron.
      double pt = true_particles[t]->PT;

      // Baseline selection
      if pt
        < 4.5 continue;
      if true_particles
        [i] -> Eta > 2.47 continue;

      double eff = rand() / RAND_MAX + 1.;

      for (int i = 0; i < effs.size(); i++)
      {
        if (pt > pt_edges[i][0] and pt < pt_edges[i][1])
        {
          if (eff < effs[i])
          {
            // TODO: Select This Electon
          }
        }
      }
    }
    else if (abs(true_particles[t]->PID == 13) and true_particles[t]->Status == 1)
    {
      // Is a final state muon
      double pt = true_particles[t]->PT;

      // Baseline Conditions
      if pt
        < 3. continue;
      if true_particles
        [i] -> Eta > 2.5 continue;

      double eff = rand() / RAND_MAX + 1.;

      for (int i = 0; i < effs.size(); i++)
      {
        if (pt > pt_edges[i][0] and pt < pt_edges[i][1])
        {
          if (eff < effs[i])
          {
            // TODO: Select This Muon
          }
        }
      }
    }
  }
}

bool Atlas_2106_01676::preselection_onWZ(bool cutflow)
{

  if (baselineElectrons.size() + baselineMuons.size() != 3)
    return false;
  if (signalElectrons.size() + signalMuons.size() != 3)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_01_3leptons");

  //! Trigger Condition to be implemented here

  if (signalLeptons[0]->PT < 25. or signalLeptons[1]->PT < 20. or signalLeptons[2]->PT < 10.)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_02_pTlep");

  if (pTmiss.Perp() < 50.)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_03_ETmiss");

  if (bveto)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_04_bveto");

  if (mll < 12.)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_05_resveto");

  if (!SFOS)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_06_SFOS");

  if (mllmax < 75. or mllmax > 105.)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_07_onZ");
  
  if (fabs((leptons[0]->P4() + leptons[1]->P4() + leptons[2]->P4()).M() - 91.2) < 15.)
    return false;
  if (cutflow)
    countCutflowEvent("onWZ_08_Zveto");

  
}


{
  // Setup mt2 evaluation object with top quarks. Modfified from AnalysisBase
  mt2_bisect::mt2 mt2_event;
  TLorentzVector zeroVector = TLorentzVector(0., 0., 0., 0.);
  // If no invis is given, use missingET. Note that pmiss[0] is irrelvant, which is why we set it to -1.
  double pmiss[] = {-1, missingET->P4().Px(), missingET->P4().Py()};
  if (invis != zeroVector)
  {
    pmiss[0] = -1;
    pmiss[1] = invis.Px();
    pmiss[2] = invis.Py();
  }

  // Construct arrays that mt2_bisect needs as input and start evaluation
  double p1[3] = {0.0, vis1.Px(), vis1.Py()};
  double p2[3] = {0.0, vis2.Px(), vis2.Py()};
  mt2_event.set_momenta(p1, p2, pmiss);
  mt2_event.set_mn(m_inv);
  return mt2_event.get_mt2();
}