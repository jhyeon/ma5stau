#include "SampleAnalyzer/User/Analyzer/ATLAS_SUSY_2018_04.h"
#include <iostream>
#include <ctime>
using namespace MA5;
using namespace std;

// Overlap Removal
template<typename T1, typename T2> std::vector<const T1*>
  Removal(std::vector<const T1*> &v1, std::vector<const T2*> &v2,
  const double &drmin)
{
  // Determining which objects should be removed
  std::vector<bool> mask(v1.size(),false);
  for (unsigned int j=0;j<v1.size();j++)
    for (unsigned int i=0;i<v2.size();i++)
      if (v2[i]->dr(v1[j]) < drmin)
      {
        mask[j]=true;
        break;
      }

  // Building the cleaned container
  std::vector<const T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
}
// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool ATLAS_SUSY_2018_04::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  INFO << "    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
  INFO << "    <>    Analysis: ATLAS SUSY 2018 04                                  <>" << endmsg;
  INFO << "    <>    SUSY, stau, hadronic decaying ditau + MET @ 13 TeV, 139 fb^-1 <>" << endmsg;
  INFO << "    <>    arXiv:1911.06660                                              <>" << endmsg;
  INFO << "    <>    Recasted by : Jongwon Lim, Chih-Ting Lu,                      <>" << endmsg;
  INFO << "    <>                  Jae-hyeon Park, Jiwon Park                      <>" << endmsg;
  INFO << "    <>    Contact     : jongwon.lim@cern.ch, timluyu@gmail.com,         <>" << endmsg; 
  INFO << "    <>                  jhpark@kias.re.kr, jiwon.park@cern.ch           <>" << endmsg;
  INFO << "    <>    Based on MadAnalysis 5 v1.8 and above                         <>" << endmsg;
  INFO << "    <>    For more information, see                                     <>" << endmsg;
  INFO << "    <>    http://madanalysis.irmp.ucl.ac.be/wiki/PublicAnalysisDatabase <>" << endmsg;
  INFO << "    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;

  // ========================= //
  // ===== Signal region ===== //
  // ========================= //

  Manager()->AddRegionSelection("SRlow");
  Manager()->AddRegionSelection("SRhigh");

  std::string SRlow[]     = {"SRlow"};
  std::string SRhigh[]    = {"SRhigh"};
  std::string SRlowhigh[] = {"SRlow","SRhigh"};

  // ====================== //
  // ===== Selections ===== //
  // ====================== //

  // Baseline cut 
  Manager()->AddCut("Baseline cut",            SRlowhigh);

  // Trigger and offline cuts
  Manager()->AddCut("asymmetric di-$\\tau$ trigger",     SRlow);
  Manager()->AddCut("di-$\\tau +E^{miss}_{T}$ trigger",  SRhigh);

  // Common selections
  Manager()->AddCut("2 medium $\\tau$ (OS) and 3rd medium $\\tau$ veto",   SRlowhigh);
  Manager()->AddCut("b-jet veto",              SRlowhigh);
  Manager()->AddCut("light lepton veto",       SRlowhigh);
  Manager()->AddCut("Z/H veto",                SRlowhigh);

  // SR Low Mass selections
  Manager()->AddCut("$75 < E^{miss}_{T} < 150$ GeV",  SRlow);
  Manager()->AddCut("2 tight $\\tau$ (OS)",           SRlow);
  Manager()->AddCut("$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad],low", SRlow);
  Manager()->AddCut("$\\Delta R(\\tau_{1},\\tau_{2})<3.2$,low",            SRlow);
  Manager()->AddCut("$m_{T2}>70$ GeV,low",                                 SRlow);

  // SR High Mass selections
  Manager()->AddCut("$\\geq 1$ tight $\\tau$",           SRhigh);
  Manager()->AddCut("$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad],high", SRhigh);
  Manager()->AddCut("$\\Delta R(\\tau_{1},\\tau_{2})<3.2$,high",            SRhigh);
  Manager()->AddCut("$m_{T2}>70$ GeV,high",                                 SRhigh);

  // ====================== //
  // ===== Histograms ===== //
  // ====================== //

  Manager()->AddHisto("SRlow_mT2", 5,70.0,120., "SRlow");

  Manager()->AddHisto("SRhigh_mT2", 5,70.0,220.,  "SRhigh");

  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void ATLAS_SUSY_2018_04::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files){}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool ATLAS_SUSY_2018_04::Execute(SampleFormat& sample, const EventFormat& event)
{
  // Event weight
  double myWeight;
  
  if      ( Configuration().IsNoEventWeight() ) myWeight=1.;
  else if ( event.mc()->weight()!=0. ) myWeight=event.mc()->weight();
  else return false;
  myWeight=1.;
  Manager()->InitializeForNewEvent(myWeight);

  // Security for empty events
  if ( event.rec()==0 ) return true;
  event_num++;

  DEBUG << "============== Event " << event_num << endmsg;



  // ================================ //
  // ===== Event reconstruction ===== //
  // ================================ //

  //// Jets ////
  std::vector <const RecJetFormat*> Jets;
  unsigned int nb = 0;
  double HT = 0.;
  for( unsigned int ij=0; ij<event.rec()->jets().size(); ij++ ){
    const RecJetFormat *Jet = &(event.rec()->jets()[ij]);
    double eta = Jet->abseta();
    double pt = Jet->pt();

    if( pt>20. && eta<2.8 ) { 
      Jets.push_back(Jet);

      if( Jet->btag() && pt>20.0 && eta<2.5 ){
        nb++;
        HT += Jet->pt();
      }
    }
  }


  //// Electrons ////
  std::vector<const RecLeptonFormat*> SignalElectrons;
  for( unsigned int ie=0; ie<event.rec()->electrons().size(); ie++ ){
    const RecLeptonFormat *Lep = &(event.rec()->electrons()[ie]);

    // Kinematics
    double eta = Lep->abseta();
    double pt = Lep->pt();

    // Isolation
    double iso_dR = std::min(10./pt,0.2);
    double iso_tracks = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(), iso_dR, 0., IsolationEFlow::TRACK_COMPONENT);
    double iso_all    = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(), 0.2, 0., IsolationEFlow::ALL_COMPONENTS);
    bool iso = (iso_tracks<0.15 && iso_all<0.20);
    if( pt>200. ) iso = (iso_all<std::max(0.015, 3.5/pt));

    // Signal electrons
    if( eta < 2.47 && pt > 17. && iso) SignalElectrons.push_back(Lep);
  }


  //// Muons ////
  std::vector<const RecLeptonFormat*> SignalMuons;
  for( unsigned int im=0; im<event.rec()->muons().size(); im++ ){
    const RecLeptonFormat *Lep = &(event.rec()->muons()[im]);

    // Kinematics
    double eta = Lep->abseta();
    double pt = Lep->pt();

    // Isolation
    double iso_dR = std::min(10./pt,0.3);
    double iso_tracks = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(), iso_dR, 0., IsolationEFlow::TRACK_COMPONENT);
    double iso_all    = PHYSICS->Isol->eflow->relIsolation(Lep,event.rec(), 0.2, 0., IsolationEFlow::ALL_COMPONENTS);
    bool iso = (iso_tracks<0.15 && iso_all<0.30);

    // Signal leptons
    if( eta < 2.7 && pt > 14. && iso ) SignalMuons.push_back(Lep);
  }


  //// Taus ////
  std::vector<const RecTauFormat*> SignalTaus;
  for( unsigned int it=0; it<event.rec()->taus().size(); it++ ){
    const RecTauFormat * CurrentTau = &(event.rec()->taus()[it]);

    if( CurrentTau->pt() > 20. && abs(CurrentTau->eta()) < 2.5
    //if( CurrentTau->pt() * 0.9 > 20. && abs(CurrentTau->eta()) < 2.5
        && (abs(CurrentTau->eta()) < 1.37 || abs(CurrentTau->eta()) > 1.52)){
      SignalTaus.push_back(CurrentTau);
      //if(SignalTaus.size() == 0 && CurrentTau->pt() > 50) SignalTaus.push_back(CurrentTau);
      //else if(SignalTaus.size() >= 1 && CurrentTau->pt() > 40) SignalTaus.push_back(CurrentTau);
    }
  }


  // Electron-electron overlap removal 
  SignalElectrons = Removal(SignalElectrons, SignalElectrons, 0.05);

  // Electron-tau overlap removal
  SignalTaus = Removal(SignalTaus, SignalElectrons, 0.2);
  // Muon-tau overlap removal
  SignalTaus = Removal(SignalTaus, SignalMuons, 0.2);

  // Muon-electron overlap removal
  SignalElectrons = Removal(SignalElectrons, SignalMuons, 0.01);

  // Electron-jet overlap removal
  //Jets = Removal(Jets, SignalElectrons, 0.2);
  Jets = PHYSICS->Isol->JetCleaning(Jets, SignalElectrons, 0.2);

  // Jet-electron overlap removal
  //Jets = PHYSICS->Isol->JetCleaning(Jets, SignalElectrons, 0.4);
  SignalElectrons = Removal(SignalElectrons, Jets, 0.2);

  // Muon-jet overlap removal
  //Jets = Removal(Jets, SignalMuons, 0.2);
  Jets = PHYSICS->Isol->JetCleaning(Jets, SignalMuons, 0.2);
  // Jet-muon overlap removal
  //Jets = PHYSICS->Isol->JetCleaning(Jets, SignalMuons, 0.4);
  SignalMuons = Removal(SignalMuons, Jets, 0.4);

  // Tau-jet overlap removal
  Jets = Removal(Jets, SignalTaus, 0.2);

  //// MET ////
  MALorentzVector pTmiss = event.rec()->MET().momentum();
  pTmiss.SetPz(0.); // Remove eta component

  //MALorentzVector pTmiss(0,0,0,0);
  //for( unsigned int i=0; i<event.rec()->photons().size(); i++ ){
  //  MALorentzVector temp = event.rec()->photons()[i].momentum();
  //  pTmiss = pTmiss + temp;
  //}

  //for( unsigned int i=0; i<event.rec()->muons().size(); i++ ){
  //  MALorentzVector temp = event.rec()->muons()[i].momentum();
  //  pTmiss = pTmiss + temp;
  //}
  //for( unsigned int i=0; i<event.rec()->jets().size(); i++ ){
  //  MALorentzVector temp = event.rec()->jets()[i].momentum();
  //  pTmiss = pTmiss + temp;
  //}

  //for( unsigned int i=0; i<Jets.size(); i++ ){
  //  MALorentzVector temp = Jets.at(i)->momentum();
  //  pTmiss = pTmiss + temp;
  //}
  //for( unsigned int i=0; i<SignalMuons.size(); i++ ){
  //  MALorentzVector temp = SignalMuons.at(i)->momentum();
  //  pTmiss = pTmiss + temp;
  //}
  //for( unsigned int i=0; i<SignalElectrons.size(); i++ ){
  //  MALorentzVector temp = SignalElectrons.at(i)->momentum();
  //  pTmiss = pTmiss + temp;
  //}
  double MET = pTmiss.Pt();
  //pTmiss = -1 * pTmiss;
  //pTmiss.SetPz(0.);
  //cout << MET << " " << event.rec()->MET().pt() << endl;
  //cout << MET << ", ";
  //cout << event.rec()->MET().pt() << ", ";

  DEBUG << "    * Event reconstructed properly..." << endmsg;


  // ============================ //
  // ===== Event selection ===== //
  // =========================== //

  //srand((unsigned int)time(NULL));
  srand(event_num);
  // Select '2018' events according to the lumi weight
  // 139 /fb = 36.2 + 44.3 + 58.5
  bool is2018 = (rand() % 100)/100. < 0.421;

  //// Baseline cut ////
  if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt()>50 && SignalTaus[1]->pt()>40,"Baseline cut") ) return true;
  //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && 0.9* SignalTaus[0]->pt()>50 && 0.9* SignalTaus[1]->pt()>40,"Baseline cut") ) return true;

  double eff=myWeight*0.8;
  Manager()->SetCurrentEventWeight(eff);

  //// SRlow cut-1 : asymmetric di-tau trigger. ////
  //bool effl = double(rand())/RAND_MAX < 0.8;
  if( is2018 ){
    if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt()>95 && SignalTaus[1]->pt()>75,"asymmetric di-$\\tau$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt()>95 && SignalTaus[1]->pt()>75 && effl,"asymmetric di-$\\tau$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && 0.9* SignalTaus[0]->pt()>95 && 0.9* SignalTaus[1]->pt()>75 && effl,"asymmetric di-$\\tau$ trigger") )
      return true;
  }
  else{
    if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt()>95 && SignalTaus[1]->pt()>60,"asymmetric di-$\\tau$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt()>95 && SignalTaus[1]->pt()>60 && effl,"asymmetric di-$\\tau$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && 0.9* SignalTaus[0]->pt()>95 && 0.9* SignalTaus[1]->pt()>60 && effl,"asymmetric di-$\\tau$ trigger") )
      return true;
  }


  //// SRhigh cut-1 : di-tau +mET trigger. ////
  //bool effh = double(rand())/RAND_MAX < 0.8;
  if( is2018 ){
    if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt() > 75 && SignalTaus[1]->pt() > 40 && MET > 150., "di-$\\tau +E^{miss}_{T}$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt() > 75 && SignalTaus[1]->pt() > 40 && effh && MET > 150., "di-$\\tau +E^{miss}_{T}$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && 0.9* SignalTaus[0]->pt() > 75 && 0.9* SignalTaus[1]->pt() > 40 && effh && MET > 150., "di-$\\tau +E^{miss}_{T}$ trigger") )
      return true;
  }
  else{
    if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt() > 50 && SignalTaus[1]->pt() > 40 && MET > 150., "di-$\\tau +E^{miss}_{T}$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && SignalTaus[0]->pt() > 50 && SignalTaus[1]->pt() > 40 && effh && MET > 150., "di-$\\tau +E^{miss}_{T}$ trigger") )
    //if( !Manager()->ApplyCut(SignalTaus.size() > 1 && 0.9* SignalTaus[0]->pt() > 50 && 0.9* SignalTaus[1]->pt() > 40 && effh && MET > 150., "di-$\\tau +E^{miss}_{T}$ trigger") )
      return true;
  }

  //// Common cut-1 : 2 medium taus (OS) and 3rd medium tau veto. ////
  int tau_charge_sum=0;
  for( unsigned int i=0; i<SignalTaus.size(); i++ ){
    if( SignalTaus[i]->charge() > 0 ) tau_charge_sum++;
    if( SignalTaus[i]->charge() < 0 ) tau_charge_sum--;
  }
  if( !Manager()->ApplyCut(SignalTaus.size() == 2 && tau_charge_sum == 0, "2 medium $\\tau$ (OS) and 3rd medium $\\tau$ veto") ) return true;


  //// Common cut-2 : b-jet veto. ////
  if( !Manager()->ApplyCut(nb == 0,"b-jet veto") ) return true;


  //// Common cut-3 : light lepton veto. ////  
  if( !Manager()->ApplyCut(SignalElectrons.size() == 0 && SignalMuons.size() == 0, "light lepton veto") ) return true;


  //// Common cut-4 : Z/H veto. //// 
  double mtata=0;
  mtata=(SignalTaus[1]->momentum()+SignalTaus[0]->momentum()).M();
  if( !Manager()->ApplyCut(mtata > 120.0, "Z/H veto") ) return true;

  
  //// SRlow cut-2 : 75 < ETmiss < 150 GeV. ////
  if( !Manager()->ApplyCut((MET > 75. && MET < 150.), "$75 < E^{miss}_{T} < 150$ GeV") ) return true;


  double tight_low=1565./2228.0;
  Manager()->SetCurrentEventWeight(eff*tight_low);


  //// SRlow cut-3 : 2 tight taus (OS). //// 
  // 2228 and 1565 are Nraw before and after 2-tight-tau cut in SR-lowMass at
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-04/tabaux_01.pdf
  //bool tight1 = double(rand())/RAND_MAX < sqrt(1565./2228.0);
  //bool tight2 = double(rand())/RAND_MAX < sqrt(1565./2228.0);
  // bool two_tight_ratio = (rand() % 100)/100. < 0.7;

  //if( !Manager()->ApplyCut(tight1 && tight2, "2 tight $\\tau$ (OS)") ) return true;
  //if( !Manager()->ApplyCut(two_tight_ratio, "2 tight $\\tau$ (OS)") ) return true;
  if( !Manager()->ApplyCut(true, "2 tight $\\tau$ (OS)") ) return true;


  // Histograms
  //Manager()->FillHisto("SRlow_mT2",mt2_low);


  //// SRlow cut-4 : |dphi(ta1,ta2)|>0.8 [rad]. ////
  double DeltaPhiTau = 999999.9;
  DeltaPhiTau = SignalTaus[0]->dphi_0_pi(SignalTaus[1]);
  if( !Manager()->ApplyCut((fabs(DeltaPhiTau)) > 0.8, "$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad],low") ) return true;


  //// SRlow cut-5 : dR(ta1,ta2)<3.2. ////
  if( !Manager()->ApplyCut(((SignalTaus[0]->momentum()).DeltaR(SignalTaus[1]->momentum()) < 3.2), "$\\Delta R(\\tau_{1},\\tau_{2})<3.2$,low"))
    return true;


  //// SRlow cut-6 : mT2>70 GeV. ////
  double mt2_high = PHYSICS->Transverse->MT2(SignalTaus[0],SignalTaus[1],event.rec()->MET(),0.);
  if( !Manager()->ApplyCut(mt2_high > 70, "$m_{T2}>70$ GeV,low") ) return true;

  //eff of 2 tight tau = 1565./2228.0 = p, from 120 GeV SRlow
  //eff of 1 tight tau = sqrt(p)
  //eff of at least 1 tau = p^2 + 2*(1-p)*p = 0.973791
  Manager()->SetCurrentEventWeight(eff*0.973791);


  //// SRhigh cut-3 : >= 1 tight tau. ////
  //tight1 = double(rand())/RAND_MAX < sqrt(1565./2228.0);
  //tight2 = double(rand())/RAND_MAX < sqrt(1565./2228.0);
 
  //if(!Manager()->ApplyCut(SignalTaus.size() == 2,"$\\geq 1$ tight $\\tau$")) return true;
  //if( !Manager()->ApplyCut(tight1 or tight2,"$\\geq 1$ tight $\\tau$") ) return true;
  if( !Manager()->ApplyCut(true, "$\\geq 1$ tight $\\tau$") ) return true;


  //// SRhigh cut-4 : |dphi(ta1,ta2)|>0.8 [rad]. ////
  DeltaPhiTau = SignalTaus[0]->dphi_0_pi(SignalTaus[1]);
  if( !Manager()->ApplyCut((fabs(DeltaPhiTau)) > 0.8, "$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad],high") ) return true;


  //// SRhigh cut-5 : dR(ta1,ta2)<3.2. ////
  if( !Manager()->ApplyCut(((SignalTaus[0]->momentum()).DeltaR(SignalTaus[1]->momentum()) < 3.2), "$\\Delta R(\\tau_{1},\\tau_{2})<3.2$,high")) 
    return true;


  //// SRhigh cut-6 : mT2>70 GeV. ////
  //MALorentzVector tau1 = SignalTaus[0]->momentum();
  //MALorentzVector tau2 = SignalTaus[1]->momentum();
  //MALorentzVector tau1 = 0.9* SignalTaus[0]->momentum();
  //MALorentzVector tau2 = 0.9* SignalTaus[1]->momentum();
  //double mt2_high = PHYSICS->Transverse->MT2(&tau1, &tau2, pTmiss,0.);
  //double mt2_high = PHYSICS->Transverse->MT2(SignalTaus[0],SignalTaus[1],pTmiss,0.);
  if( !Manager()->ApplyCut(mt2_high > 70, "$m_{T2}>70$ GeV,high") ) return true;

  // Histograms
//  Manager()->FillHisto("SRhigh_mT2",mt2_high);


  return true;
}

