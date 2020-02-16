#include "SampleAnalyzer/User/Analyzer/ATLAS_SUSY_2018_04.h"
#include "SampleAnalyzer/User/Analyzer/lester_mt2_bisect.h"
using namespace MA5;
using namespace std;

// Overlap Removal
template<typename T1, typename T2> std::vector<const T1*>
  Removal(std::vector<const T1*> &v1, std::vector<const T2*> &v2,
  const double &drmin)
{
  // Determining with objects should be removed
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
  INFO << "    <>    Recasted by : Chih-Ting Lu                                    <>" << endmsg;
  INFO << "    <>    Contact     : timluyu@gmail.com                               <>" << endmsg;
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

  // Common selections
  Manager()->AddCut("2 medium $\\tau$ (OS)",   SRlowhigh);
  Manager()->AddCut("3rd medium $\\tau$ veto", SRlowhigh);
  Manager()->AddCut("b-jet veto",              SRlowhigh);
  Manager()->AddCut("light lepton veto",       SRlowhigh);
  Manager()->AddCut("Z/H veto",                SRlowhigh);

  // SR low Mass selections
  Manager()->AddCut("asymmetric di-$\\tau$ trigger",  SRlow);
  Manager()->AddCut("$75 < E^{miss}_{T} < 150$ GeV",  SRlow);
  Manager()->AddCut("2 tight $\\tau$ (OS)",           SRlow);
  Manager()->AddCut("$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad]", SRlow);
  Manager()->AddCut("$\\Delta R(\\tau_{1},\\tau_{2})<3.2$",            SRlow);
  Manager()->AddCut("$m_{T2}>70$ GeV",                                 SRlow);

  // SR high Mass selections
  Manager()->AddCut("di-$\\tau +E^{miss}_{T}$ trigger",  SRhigh);
  Manager()->AddCut("$E^{miss}_{T} > 150$ GeV",          SRhigh);
  Manager()->AddCut("$\\geq 1$ tight $\\tau$",           SRhigh);
  Manager()->AddCut("$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad]", SRhigh);
  Manager()->AddCut("$\\Delta R(\\tau_{1},\\tau_{2})<3.2$",            SRhigh);
  Manager()->AddCut("$m_{T2}>70$ GeV",                                 SRhigh);

  // ====================== //
  // ===== Histograms ===== //
  // ====================== //

  Manager()->AddHisto("SRlow_MET", 10,0.0,150., "SRlow");
  Manager()->AddHisto("SRlow_mT2", 5,30.0,70.,  "SRlow");

  Manager()->AddHisto("SRhigh_MET", 10,50.0,100., "SRhigh");
  Manager()->AddHisto("SRhigh_mT2", 2.5,50.0,70., "SRhigh");

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
  Manager()->InitializeForNewEvent(myWeight);

  // Security for empty events
  if ( event.rec()==0 ) return true;
  event_num++;

//  std::vector<const RecJetFormat*>    BaseJets,  SignalBJets,   SignalnonBJets, SignalJets;
//  std::vector<const RecLeptonFormat*> BaseMuons, BaseElectrons, SignalMuons,    SignalElectrons;
//  std::vector<const RecTauFormat*>    BaseTaus,  SignalTaus;
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

    // Signal leptons
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
    if( eta < 2.47 && pt > 14. && iso ) SignalMuons.push_back(Lep);
  }


  //// Taus ////
  std::vector<const RecTauFormat*> SignalTaus;
  for( unsigned int it=0; it<event.rec()->taus().size(); it++ ){
    const RecTauFormat * CurrentTau = &(event.rec()->taus()[it]);

    if( CurrentTau->pt() > 10. &&
        (abs(CurrentTau->eta()) < 2.5) )
      SignalTaus.push_back(CurrentTau);
  }



  // Electron-tau overlap removal
  SignalElectrons = Removal(SignalElectrons, SignalTaus, 0.2);
  // Muon-tau overlap removal
  SignalMuons = Removal(SignalMuons, SignalTaus, 0.2);
  // Jet-electron overlap removal
  Jets = PHYSICS->Isol->JetCleaning(Jets, SignalElectrons, 0.4);
  // Jet-muon overlap removal
  Jets = PHYSICS->Isol->JetCleaning(Jets, SignalMuons, 0.4);

  // Tau-jet overlap removal
  SignalTaus = Removal(SignalTaus, Jets, 0.2);


  //// MET ////
  MALorentzVector pTmiss = event.rec()->MET().momentum();
  pTmiss.SetPz(0.); // Remove eta component
  double MET = pTmiss.Pt();

  DEBUG << "    * Event reconstructed properly..." << endmsg;


  // ============================ //
  // ===== Event selection ===== //
  // =========================== //

  //// Common cut-1 : 2 medium taus (OS). ////
  int tau_charge_sum=0;
  for( unsigned int i=0; i<SignalTaus.size(); i++ ){
    if( SignalTaus[i]->charge() > 0 ) tau_charge_sum++;
    if( SignalTaus[i]->charge() < 0 ) tau_charge_sum--;
  }
  for( unsigned int ii=0; ii<SignalTaus.size(); ii++ ){
  if( !Manager()->ApplyCut((SignalTaus.size() > 1) && tau_charge_sum==0 && SignalTaus[ii]->pt() > 20. && (abs(SignalTaus[ii]->eta()) < 1.37 || (abs(SignalTaus[ii]->eta()) > 1.52 && abs(SignalTaus[ii]->eta()) < 2.5)), "2 medium $\\tau$ (OS)") ) return true; }


  //// Common cut-2 : 3rd medium tau veto. ////
  if( !Manager()->ApplyCut(SignalTaus.size() < 3, "3rd medium $\\tau$ veto") ) return true;


  //// Common cut-3 : b-jet veto. ////
  if( !Manager()->ApplyCut(nb < 1,"b-jet veto") ) return true;


  //// Common cut-4 : light lepton veto. ////  
  if( !Manager()->ApplyCut(SignalElectrons.size() < 1 && SignalMuons.size() < 1, "light lepton veto") ) return true;

  //// Common cut-5 : Z/H veto. //// 
  double mtata=0;
  if( SignalTaus.size()==2 ) mtata=(SignalTaus[1]->momentum()+SignalTaus[0]->momentum()).M();
  if( !Manager()->ApplyCut(mtata > 120.0, "Z/H veto") ) return true;

  //// SRlow cut-1 : asymmetric di-tau trigger. ////

  
  //// SRLow cut-2 : 75 < ETmiss < 150 GeV. ////
  if( !Manager()->ApplyCut((MET > 75. && MET < 150.), "$75 < E^{miss}_{T} < 150$ GeV") ) return true;


  //// SRLow cut-3 : 2 tight tau (OS). //// 
  if( !Manager()->ApplyCut(SignalTaus.size()==2 && SignalTaus[0]->pt()>95 && SignalTaus[1]->pt()>75,"2 tight $\\tau$ (OS)") )
    return true;

  //// SRLow cut-4 : |dphi(ta1,ta2)|>0.8 [rad]. ////
  double DeltaPhiTau = 999999.9;
  if( SignalTaus.size()==2 ) DeltaPhiTau = SignalTaus[0]->dphi_0_pi(SignalTaus[1]);
  if( !Manager()->ApplyCut((fabs(DeltaPhiTau)) > 0.8, "$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad]") ) return true;


  //// SRLow cut-5 : dR(ta1,ta2)<3.2. ////
  if( !Manager()->ApplyCut(((SignalTaus[0]->momentum()).DeltaR(SignalTaus[1]->momentum()) < 3.2), "$\\Delta R(\\tau_{1},\\tau_{2})<3.2$") ) 
    return true;

  //// SRLow cut-6 : mT2>70 GeV. ////LorentzVector tau_low1 = SignalTaus[0]->momentum();
  MALorentzVector tau_low1 = SignalTaus[0]->momentum();
  MALorentzVector tau_low2 = SignalTaus[1]->momentum();
  double mt2_low = asymm_mt2_lester_bisect::get_mT2(tau_low1.M(), tau_low1.Px(), tau_low1.Py(),
                                                    tau_low2.M(), tau_low2.Px(), tau_low2.Py(),
                                                    pTmiss.Px(), pTmiss.Py(), 1., 1.);
  if( !Manager()->ApplyCut(mt2_low > 70, "$m_{T2}>70$ GeV") ) return true;



  //// SRhigh cut-1 : di-tau +mET trigger. ////


  //// SRhigh cut-2 : ETmiss > 150 GeV. ////
  if( !Manager()->ApplyCut((MET > 150.),"E^{miss}_{T} > 150$ GeV") ) return true;

  //// SRhigh cut-3 : >= 1 tight tau. //// 
  if( !Manager()->ApplyCut(SignalTaus.size() == 2 && SignalTaus[0]->pt() > 75 && SignalTaus[1]->pt() > 40, "\\geq 1$ tight $\\tau") )
    return true;

  //// SRhigh cut-4 : |dphi(ta1,ta2)|>0.8 [rad]. ////
  if( SignalTaus.size() > 1 ) DeltaPhiTau = SignalTaus[0]->dphi_0_pi(SignalTaus[1]);
  if( !Manager()->ApplyCut((fabs(DeltaPhiTau)) > 0.8, "$|\\Delta\\phi(\\tau_{1},\\tau_{2})|>0.8$ [rad]") ) return true;

  //// SRhigh cut-5 : dR(ta1,ta2)<3.2. ////
  if( !Manager()->ApplyCut(((SignalTaus[0]->momentum()).DeltaR(SignalTaus[1]->momentum()) < 3.2), "$\\Delta R(\\tau_{1},\\tau_{2})<3.2$"))  
    return true;

  //// SRhigh cut-6 : mT2>70 GeV. ////
  MALorentzVector tau_high1 = SignalTaus[0]->momentum();
  MALorentzVector tau_high2 = SignalTaus[1]->momentum();
  double mt2_high = asymm_mt2_lester_bisect::get_mT2(tau_high1.M(), tau_high1.Px(), tau_high1.Py(),
                                                     tau_high2.M(), tau_high2.Px(), tau_high2.Py(),
                                                     pTmiss.Px(), pTmiss.Py(), 1., 1.);
  if( !Manager()->ApplyCut(mt2_high > 70, "$m_{T2}>70$ GeV") ) return true;


  return true;
}

