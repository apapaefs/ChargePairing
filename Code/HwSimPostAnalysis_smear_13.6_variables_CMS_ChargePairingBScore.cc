#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <numeric>
//ROOT include files
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <chrono>
#include <random>
#include <TH1F.h>
#include <TLorentzVector.h>

//Fastjet headers
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <fastjet/tools/JHTopTagger.hh>
#include <fastjet/Selector.hh>

//Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>

//custom headers
#include "TopHist.h"
#include "complex_d.h"

//btagscore
#include "btag_score_lib.h"

using namespace std;
using namespace fastjet;
using namespace pdf;

//----------------------------------------------------------------------
// Some four-vector operators
//----------------------------------------------------------------------
double dot(fastjet::PseudoJet p1, fastjet::PseudoJet p2);
double deltaR(fastjet::PseudoJet p1, fastjet::PseudoJet p2);

/* jet to lepton mistag */
double Pjet_to_lepton(double pt);

/* jet to photon mistag */
double Pjet_to_photon(double pt);

//----------------------------------------------------------------------// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
// command line parameters
//----------------------------------------------------------------------
char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);

//----------------------------------------------------------------------
// Analysis functions
//----------------------------------------------------------------------

// smearing of jets, leptons and photons. 
fastjet::PseudoJet smear_jet(fastjet::PseudoJet jet_in);
fastjet::PseudoJet smear_lepton(fastjet::PseudoJet lepton_in, int lepton_id);
fastjet::PseudoJet smear_photon(fastjet::PseudoJet photon_in);

// acceptance efficiency for leptons, photons, jets
bool lepton_efficiency_accept(fastjet::PseudoJet lepton_in, int lepton_id);
bool photon_efficiency_accept(fastjet::PseudoJet photon_in);
bool jet_efficiency_accept(fastjet::PseudoJet jet_in);
bool btag_hadrons(fastjet::PseudoJet jet);


//analysis functions
double analyze_event(fastjet::PseudoJet photon1, fastjet::PseudoJet photon2, fastjet::PseudoJet cjet, fastjet::PseudoJet bjet, fastjet::PseudoJet lepton, fastjet::PseudoJet etmiss, double evweight_i);
double Pb_to_b(double pt);
double Pb_to_c(double pt);
double Pc_to_b(double pt);
double Pjet_to_b(double pt);
// jet to lepton mistag prob
double Pjet_to_photon(double pt);
double btag_weight(fastjet::PseudoJet jet, bool btag, bool ctag);
double atag_weight(fastjet::PseudoJet jet, bool btag, bool ctag);
double ctag_weight(fastjet::PseudoJet jet, bool btag, bool ctag);
std::pair<fastjet::PseudoJet, fastjet::PseudoJet> get_Wvectors(fastjet::PseudoJet plepton, fastjet::PseudoJet pmiss);
std::pair<fastjet::PseudoJet, fastjet::PseudoJet> get_Wvectors_GTX(fastjet::PseudoJet plepton, fastjet::PseudoJet pmiss);
std::pair<complex_ld, complex_ld> pnuz_fromw(fastjet::PseudoJet plepton, fastjet::PseudoJet pmiss);
std::pair<complex_d, complex_d> quadsolve(complex_ld a, complex_ld b, complex_ld c);


void sort_by_btag_desc(std::vector<fastjet::PseudoJet>& jets);

//get all pairings:
void generatePairings(int* items, int itemcount, int start);

//check if integer is in vector of integers
bool is_in(int a, vector<int> intvec);

// ATLAS-style smearer 
fastjet::PseudoJet smearer_j_ATLAS_Benj(double E, double pt, double pz, double phi, double eta);

// CMS-style smearer
fastjet::PseudoJet smearer_j_CMS_Benj(double E, double pt, double pz, double phi, double eta);

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// IDs of B-hadrons used by the btag_hadrons function
int bhadronid[105] = {5122, -5122, 15122, -15122, 5124, -5124, 5334, -5334, 5114, -5114, 5214, -5214, 5224, -5224, 5112, -5112, 5212, -5212, 5222, -5222, 15322, -15322, 15312, -15312, 15324, -15324, 15314, -15314, 5314, -5314, 5324, -5324, 5132, -5132, 5232, -5232, 5312, -5312, 5322, -5322, 551, 10555, 100551, 200551, 553, 557, 555, 100555, 200555, 20523, -20523, 20513, -20513, 20543, -20543, 20533, -20533, 511, 521, -511, -521, 531, -531, 541, -541, 513, 523, -513, -523, 533, -533, 543, -543, 10513, 10523, -10513, -10523, 10533, -10533, 10543, -10543, 10511, 10521, -10511, -10521, 10531, -10531, 10541, -10541, 20513, 20523, -20513, -20523, 20533, -20533, 20543, -20543, 515, 525, -515, -525, 535, -535, 545, -545};

//particle masses
double mw = 80.4;
double mtop = 173.0;
double mhiggs1 = 120.0; //changed to correspond to shifted distributions due to neutrinos
double mhiggs2 = 115.0; //changed to correspond to shifted distributions due to neutrinos
double mhiggs3 = 110.0; //changed to correspond to shifted distributions due to neutrinos


/*
 * CREATE ROOT CHAIN TO READ IN THE FILES
 */
TChain t("Data");


/* 
 * DECLARE RANDOM NUMBERS
 */ 
TRandom3 rnd;
TRandom3 rndint;

/***** 
 ***** SWITCHES FOR SMEARING/EFFICIENCIES
 *****/
bool donotsmear_jets = 1;
bool donotsmear_leptons = 1;
bool donot_apply_efficiency = 1;
bool donotsmear_photons = 1;

/* 
 * cut counters
 */ 
double passed_lepton_cuts(0.);   


/* 
 * CUTS DEFINED HERE IN GEV
*/
 
/* JET CUTS */

double cut_pt_jet(20.0); //pt cut for jets
double cut_eta_jet(4.0); //pseudo-rapidity cut for jets

/* B-jet CUTS */
double cut_pt_bjet(20.0); //pt cut for b-jets
double cut_eta_bjet(2.5); //pseudo-rapidity cut for b-jets
double cut_dRbbmin(0.3); //minimum delta R between b-jets
double cut_pt_bjet1(80.0); //pt cut for b-jets
double cut_pt_bjet2(50.0); //pt cut for b-jets
double cut_pt_bjet3(40.0); //pt cut for b-jets
double cut_pt_bjet4(40.0); //pt cut for b-jets


/* RECO Higgs CUTS */ 
double cut_pt_higgs1(100); //minimum pt for reco higgs 1
double cut_pt_higgs2(100); //minimum pt for reco higgs 2
double cut_pt_higgs3(80); //minimum pt for reco higgs 3
double cut_chisq_min(50.); //maximum of the minimum (best) chi-squared
double cut_DeltaM_min(20); //maximum DeltaM_min
double cut_DeltaM_med(30); //maximum DeltaM_med
double cut_DeltaM_max(40); //maximum DeltaM_max
double cut_dR_higgses(3.5); //maximum delta R beteween reco Higgses
double cut_dR_hbbreco(3.5); //maximum delta R between b-jets in a reco Higgs
double cut_m6b(350.); //invariant mass of 6 b-jets cut

// The b-jet charge identification probability, i.e. the probability to CORRECTLY identify the charge of a b-jet:
double Pb_charge(1.0);
// whether to consider or ignore the b-jet charge
int chargetagging(1); 

bool perfect_tagging = true;

std::vector<std::vector<int>> pairs_of_six;

/* 
 * Class for PseudoJet btagscore
 */
class MyInfo: public PseudoJet::UserInfoBase {
public:
  MyInfo(double bscore, int trueb) : _btagscore(bscore), _trueb(trueb) {}
  double btagscore() const {return _btagscore;} //the b-tag score
  int trueb() const { return _trueb; } //whether it is a true b-jet or not
protected:
  double _btagscore;
  int _trueb;
};
  
int main(int argc, char *argv[]) {

  //take command line options
  char* output;
  char* infile = "";
  if(argv[1]) { infile = argv[1]; } else { cout << "Use: ./HwSimAnalysis [input] [options]" << endl; exit(1); }

  //set the variables and addresses to be read from root file
  //total number of particles in an event
  int numparticles;
  //total number of jets in an event
  int numJets;
  //total number of b-jets in an event
  int numbJets;
  //total number of photons in an event
  int numPhotons;
  //total number of Electrons in an event
  int numElectrons;
  //total number of Positrons in an event
  int numPositrons;
  //total number of Muons in an event
  int numMuons;
  //total number of Muons in an event
  int numantiMuons;
  /** particle information in the order: 
   * 4 momenta (E,x,y,z), id, other info
   **/
  double objects[8][1000];
  /* the event weight */ 
  double evweight;
  /* The vector of jets */
  double theJets[5][100];
  /* The vector of jets */
  double thebJets[5][100];
  /* The vector of photons */
  double thePhotons[4][100];
  /* The vector of muons */
  double theMuons[4][100];
  /* The vector of electrons */
  double theElectrons[4][100];
   /* The vector of anti-muons */
  double theantiMuons[4][100];
  /* The vector of positrons */
  double thePositrons[4][100];
  /* The missing energy four-vector */
  double theETmiss[4];
  /* c-Tagging containers corresponding to the jets */
  double cTag[100];


  /* 
   * SET THE ROOT BRANCH ADDRESSES
   */
  //comment out numparticles and objects if not saving them
  //t.SetBranchAddress("numparticles",&numparticles);
  //t.SetBranchAddress("objects",&objects);
  
  t.SetBranchAddress("evweight",&evweight);
  
  t.SetBranchAddress("theJets", &theJets);
  t.SetBranchAddress("numJets", &numJets);
  t.SetBranchAddress("cTag", &cTag);

  t.SetBranchAddress("thebJets", &thebJets);
  t.SetBranchAddress("numbJets", &numbJets);
        
  t.SetBranchAddress("thePhotons", &thePhotons);
  t.SetBranchAddress("numPhotons", &numPhotons);
        
  t.SetBranchAddress("theMuons", &theMuons);
  t.SetBranchAddress("numMuons", &numMuons);
        
  t.SetBranchAddress("theantiMuons", &theantiMuons);
  t.SetBranchAddress("numantiMuons", &numantiMuons);
  
  t.SetBranchAddress("thePositrons", &thePositrons);
  t.SetBranchAddress("numPositrons", &numPositrons);
        
  t.SetBranchAddress("theElectrons", &theElectrons);
  t.SetBranchAddress("numElectrons", &numElectrons);
        
  t.SetBranchAddress("theETmiss", &theETmiss);


  /* Set up random number
   * generator
   */ 
  rnd.SetSeed(14101983);


  /* Add up all the input 
   * files to the chain
   */ 
  string stringin = "";
  ifstream inputlist;
  if (std::string(infile).find(".input") != std::string::npos) {
    inputlist.open(infile);
    if(!inputlist) {  cerr << "Error: Failed to open input file " << infile << endl; exit(1); }
    while(inputlist) { 
      inputlist >> stringin; 
      if(stringin!="") { t.Add(stringin.c_str()); 
        cout << "Adding " << stringin.c_str() << endl;
      }
      stringin = "";
    }
    inputlist.close();
  } else if (std::string(infile).find(".root") != std::string::npos) {
    cout << "Adding " << infile << endl;
    t.Add(infile);
  }

  /* Get Number of events
   * and print
   */
  int EventNumber(int(t.GetEntries()));
  cout << "Total number of events in " << infile << " : " << EventNumber << endl;

  // warn if the perfect tagging flag is on
  if(perfect_tagging) cout << "WARNING: Perfect tagging of b-jets and c-jets and no mis-tagging of them enabled" << endl;

  /* 
   * -b: USED TO REANALYZE PREVIOUSLY PASSED EVENTS ONLY, DEFAULT IS ALL EVENTS 
   */
  
  //whether the analysis performed is level-2 or level-3 
  bool basic = true;
  if(cmdOptionExists(argv, argv+argc, "-b")) {
    cout << "Looking for .evp2 file, running over all events" << endl;
    basic = false;
  }


  /* 
   * -t: ADD AN EXTENSIN TAG TO YOUR OUTPUT FILES
   */
  string tag;
  tag = "";
  if(cmdOptionExists(argv, argv+argc, "-t")) {
    tag = getCmdOption(argv, argv + argc, "-t");
    tag = "-" + tag;
    cout << "Adding tag: " << tag << endl;
  }

  /* 
   * -n: RUN FROM START OF FILE UP TO A GIVEN NUMBER OF EVENTS
   */
  char * switch_maxevents;
  char * switch_minevents;
  int maxevents(0), minevents(0);
  if(cmdOptionExists(argv, argv+argc, "-n")) {
    switch_maxevents = getCmdOption(argv, argv + argc, "-n");  
    maxevents=(atoi(switch_maxevents));	       
    if(maxevents > EventNumber) { maxevents = EventNumber; } 
    cout << "Analyzing up to " << maxevents << endl;
    if(maxevents < 1 || maxevents > 1E10) { cout << "Error: maxevents must be in the range [1,1E10]" << endl; exit(1); } 
  }
  
  /* 
   * -nmax: RUN FROM START OF FILE UP TO A GIVEN NUMBER OF EVENTS, TO BE USED IN CONJUNCTION WITH -nmin
   */
  //maximum number of events to analyze
  if(cmdOptionExists(argv, argv+argc, "-nmax") && !cmdOptionExists(argv, argv+argc, "-n")) {
   switch_maxevents = getCmdOption(argv, argv + argc, "-nmax");  
    maxevents=(atoi(switch_maxevents));	       
    if(maxevents > EventNumber) { maxevents = EventNumber; } 
    cout << "Analyzing up to " << maxevents << endl;
    if(maxevents < 1 || maxevents > 1E10) { cout << "Error: maxevents must be in the range [1,1E10]" << endl; exit(1); } 
  } 
  if(!cmdOptionExists(argv, argv+argc, "-nmax") && !cmdOptionExists(argv, argv+argc, "-n")) { maxevents = EventNumber; }

  /* 
   * -nmin: RUN FROM POINT nmin OF FILE UP TO A GIVEN NUMBER OF EVENTS SPECIFIED BY -nmax
   */
  //starting number of events to analyse
  if(cmdOptionExists(argv, argv+argc, "-nmin")) {
    switch_minevents = getCmdOption(argv, argv + argc, "-nmin");  
    minevents=(atoi(switch_minevents));	       
    if(minevents > maxevents) { minevents = 0; }
    cout << "Analyzing from " << minevents << endl;
    if(minevents < 1 || minevents > 1E10) { cout << "Error: minevents must be in the range [1,1E10]" << endl; exit(1); } 
  }

  /* 
   * -pbc: CHANGE THE CHARGE TAGGING PROBABILITY
   */
  char * switch_pbc;
  if(cmdOptionExists(argv, argv+argc, "-pbc")) {
    switch_pbc = getCmdOption(argv, argv + argc, "-pbc");
    Pb_charge=(atof(switch_pbc));	       
    if(atof(switch_pbc) > 1.0 || atof(switch_pbc) < 0.5) { cout << "Charge identification probability has to lie in [0.5, 1.0]" << endl; exit(1); }
    cout << "B-jet charge identification probability: " << Pb_charge << endl;

  }

  /* 
   * -chtag: SWITCH CHARGE TAGGING ON AND OFF COMPLETELY: 1 is on, 0 is off
   */
  char * switch_chtag;
  if(cmdOptionExists(argv, argv+argc, "-chtag")) {
    switch_chtag = getCmdOption(argv, argv + argc, "-chtag");
    chargetagging=(atoi(switch_chtag));	       
    if(atoi(switch_chtag) != 1 && atoi(switch_chtag) != 0) { cout << "Choose whether charge tagging is on (1) or off (0)" << endl; exit(1); }
    cout << "Charge pairing switch: " << chargetagging << endl;
  }

  /* 
   * CREATE THE OUTPUT FILE STRINGS 
   */
  string outnew = "";
  outnew = std::string(infile);
  string replacement = tag + ".top";
  boost::replace_all(outnew, ".root", replacement);
  boost::replace_all(outnew, ".input", replacement);        
  char* output2 = new char[outnew.length() + 1];
  //  cout << outnew.c_str() << endl;
  strcpy (output2, outnew.c_str());
  output = output2;
  
  char* output_dat;
  string outnew2 = "";
  outnew2 = std::string(infile);
  replacement = tag + ".smearCMS.dat";
  boost::replace_all(outnew2, ".root", replacement);
  boost::replace_all(outnew2, ".input", replacement);        
  char* output3 = new char[outnew2.length() + 1];
  strcpy (output3, outnew2.c_str());
  output_dat = output3;
  ofstream outdat(output_dat, ios::out);

  //load events that have passed the second stage of analysis
  //if basic = false;
  string ineventpass;
  ifstream inevt;
  string inevt_curr;
  int passed_event[20000];  
  int npassed_previous(0);
  if(basic == false) { 
    ineventpass = std::string(infile);
    replacement = tag + ".evp";

    boost::replace_all(ineventpass,".input", replacement);
    boost::replace_all(ineventpass,".root", replacement);
    inevt.open(ineventpass.c_str());
    if(!inevt) { cerr << "Error: Cannot open "<< ineventpass.c_str() << endl; exit(1); } 
    for(int ii = 0; ii < 1000; ii++) { passed_event[ii] = -1; }
    while(inevt) { 
      inevt >> inevt_curr;
      // cout << inevt_curr.c_str() << endl;
      passed_event[npassed_previous] = atoi(inevt_curr.c_str());
      npassed_previous++;
    }
  }
  //for(int pp = 0; pp < npassed_previous; pp++) { coust << passed_event[pp] << endl; }
 
  string outeventpass = ""; 
  ofstream outevp;

  if(basic == false) { 
    outeventpass = std::string(infile);
    replacement = tag + ".evp2";
    boost::replace_all(outeventpass,".root", replacement);
    boost::replace_all(outeventpass,".input", replacement);
    boost::replace_all(outeventpass,".top", replacement);
    outevp.open(outeventpass.c_str());
  } else if(basic == true) {
    outeventpass = std::string(infile);
    replacement = tag + ".evp";
    boost::replace_all(outeventpass,".root", replacement);
    boost::replace_all(outeventpass,".input", replacement);
    boost::replace_all(outeventpass,".top", replacement);
    outevp.open(outeventpass.c_str());
  }

  /*
   * PREPARES THE OUTPUT ARRAY FOR *_var.root: USED FOR FURTHER ANALYSIS
   */
  std::cout << "Preparing Root Tree for event variables" << endl;
  TTree* Data2;
  TFile* dat2;
  string fnameroot = std::string(infile);
  replacement = tag + "_var.smearCMS.root";
  boost::replace_all(fnameroot,".root", replacement);
  boost::replace_all(fnameroot,".input", replacement);
  dat2 = new TFile(fnameroot.c_str(), "RECREATE");
  Data2 = new TTree ("Data2", "Data Tree");
  //variables to fill in the .root file
  double variables[21]; 
  double weight;
  Data2->Branch("variables", &variables, "variables[21]/D");
  Data2->Branch("weight", &weight);

 
  /* 
   * COUNTERS FOR NUMBER OF EVENTS THAT PASS CUTS
   */
  double pass_6b(0); //passed reconstruction of 6 b-jets
  double pass_ptb(0); //passed minimum pt of all 6 b-jets
  double pass_drbb(0); //passed dR > 0.3 for all pairings
  double pass_pthiggses(0); //passed reco higgses pT cuts
  double pass_chisq(0); //passed the chi-squared minimum cut
  double pass_DeltaM(0); //passed delta M cut
  double pass_dRhiggses(0); //passed reco higgses dR cuts
  double pass_dRbbhiggses(0); //passed reco higgses dR(b,b) cuts
  double pass_m6b(0); //passed m6b cut
  
  double passcuts(0); //passed all cuts
  double eventcount(0); //counting the events (no weights)
  double total_event_in(0); //counting the events in (no weights)
  double total_weight_in(0);//the total weight before analysis


  /* 
   * PARAMETERS AND SWITCHES
   */
 
  /* 
   * HISTOGRAMS DEFINED HERE 
   */
  TopHist h_dummy(10,output,"dummy histo", 0,1);
  TopHist h_pT_jets(60,output,"pT of all jets",0, 300);
  TopHist h_pT_leptons(60,output,"pT of leptons",0, 300);
  TopHist h_pT_b(60,output,"pT of reco b jets",0, 300);
  TopHist h_pT_b1(60,output,"pT of reco b jet 1",0, 300);
  TopHist h_pT_b2(60,output,"pT of reco b jet 2",0, 300);
  TopHist h_pT_b3(60,output,"pT of reco b jet 3",0, 300);
  TopHist h_pT_b4(60,output,"pT of reco b jet 4",0, 300);
  TopHist h_pT_b5(60,output,"pT of reco b jet 5",0, 300);
  TopHist h_pT_b6(60,output,"pT of reco b jet 6",0, 300);
  TopHist h_chisq(60,output,"chisq min",0, 150);

  TopHist h_DeltaM_min(60,output,"Delta M min",0, 150);
  TopHist h_DeltaM_med(60,output,"Delta M med",0, 150);
  TopHist h_DeltaM_max(60,output,"Delta M max",0, 150);
  TopHist h_pT_h1(60,output,"pT of Higgs 1",0, 300);
  TopHist h_pT_h2(60,output,"pT of Higgs 2",0, 300);
  TopHist h_pT_h3(60,output,"pT of Higgs 3",0, 300);
  TopHist h_pT_dRhh(60,output,"delta R between Higgs bosons",0, 3.14153);
  TopHist h_m6b(100,output,"6b invariant mass",0, 1000);

  TopHist h_numbjets(12,output,"number of bJets",-2, 10);

  
  TopHist h_chargesum_initial(240,output,"6 highest-pT b-Jet INITIAL charge sum (absolute value)",-2, 10);
  TopHist h_chargesumabs_initial(240,output,"6 highest-pT b-Jet INITIAL charge_sum_abs (sum of absolute values)",-2, 10);
  TopHist h_chargesum_final(240,output,"6 highest-pT b-Jet FINAL charge sum (absolute value)",-2, 10);
  TopHist h_chargesumabs_final(240,output,"6 highest-pT b-Jet FINAL charge_sum_abs (sum of absolute values)",-2, 10);


  bool RESET_WEIGHTS_TO_UNITY = false;
  if(RESET_WEIGHTS_TO_UNITY) {
    cout << "WARNING: RESETTING ALL WEIGHTS TO = 1" << endl;
  }

  
  int listind[6] = {0,1,2,3,4,5};
  generatePairings(listind, 6, 0);

  for(int pp = 0; pp < pairs_of_six.size(); pp++) {
    for(int ps = 0; ps < pairs_of_six[pp].size(); ps++) { 
      cout << pairs_of_six[pp][ps] << ",";
    }
    cout << endl;
  }





  /*
   * Construct the maps for the btag score 
   */ 

  std::string filePath_NonB = "../btagscore/JetBtagDeepFlavB_NonB_Distributions.root";

  std::string filePath_B ="../btagscore/JetBtagDeepFlavB_B_Distributions.root";

  std::map<
    std::pair<int, int>,
    std::pair<std::vector<double>, std::vector<double>>
    >   JetBtagDeepFlavB_B_Distributions = read_btagscore_histograms(filePath_B);

  std::map<
    std::pair<int, int>,
    std::pair<std::vector<double>, std::vector<double>>
    >   JetBtagDeepFlavB_NonB_Distributions = read_btagscore_histograms(filePath_NonB);
  
  /*
   *
   * LOOP OVER EVENTS
   * AND
   * PERFORM ANALYSIS
   *
   */

  bool perform_analysis_on_event = false;
  for(int ii = minevents; ii < maxevents; ii++) {

    bool passed_all_cuts = true;
    
    /* IF LEVEL 3 ANALYSIS THEN
     * CHECK IF EVENT IS IN .evp FILE
     */ 
    perform_analysis_on_event = false;
    if(basic == false) { 
       for(int pp = 0; pp < npassed_previous; pp++) { if(ii == passed_event[pp]) { perform_analysis_on_event = true; } }
    }
    if(!perform_analysis_on_event && basic == false) { continue; }

    /* GRAB EVENT ENTRY
     * FROM ROOT FILE
     * AND PRINT EVENT NUMBER
     */
    t.GetEntry(ii);
    if(RESET_WEIGHTS_TO_UNITY) { 
      evweight = 1.0;
    }

    if(ii%10 == 0) { cout << "Event number: " << ii << "\r" << flush; }

    /*
     * PUSH BACK JETS, LEPTONS & PHOTONS INTO PSEUDOJETS
     */
    
    vector<fastjet::PseudoJet> bJets_unsort, Jets_unsort, Electrons_unsort, Positrons_unsort, Muons_unsort, AntiMuons_unsort, Photons_unsort, Leptons_unsort;

    fastjet::PseudoJet bjetcan, jetcan, photoncan;
    for(int jj = 0; jj < numJets; jj++) {
      jetcan = fastjet::PseudoJet(theJets[1][jj], theJets[2][jj], theJets[3][jj], theJets[0][jj]);
      //jetcan.set_user_index(cTag[jj]);
      
      fastjet::PseudoJet jetcan_smeared = smearer_j_CMS_Benj(jetcan.e(), jetcan.perp(), jetcan.pz(),
							     jetcan.phi(), jetcan.eta());
      //get the btagscore for the non-b-jet
      double btagscore = 0;
      //check if within tagging region
      if(fabs(jetcan_smeared.eta()) < 2.5 && jetcan_smeared.perp() > 20 && jetcan_smeared.perp() < 1000) { btagscore = jet_btagscore(jetcan_smeared.perp(), fabs(jetcan_smeared.eta()), JetBtagDeepFlavB_NonB_Distributions); }
      //add it to the MyInfo Class for the PseudoJet:
      jetcan_smeared.set_user_info(new MyInfo(btagscore, 0));
      jetcan_smeared.set_user_index( (rnd.Rndm() < 0.5) ? -1 : 1 );  //random charge for non-bjets
      //jetcan_smeared.set_user_index(0); //charge zero for non-bjets
      //cout << "jetcan_smeared.user_info<MyInfo>().btagscore()= " << jetcan_smeared.user_info<MyInfo>().btagscore() << endl;
      if(jetcan.perp() > cut_pt_jet && fabs(jetcan.eta()) < cut_eta_jet && jet_efficiency_accept(jetcan)) { 
	Jets_unsort.push_back(jetcan_smeared);
      }
    }
    for(int jj = 0; jj < numbJets; jj++) {
      bjetcan = fastjet::PseudoJet(thebJets[1][jj], thebJets[2][jj], thebJets[3][jj], thebJets[0][jj]);
      fastjet::PseudoJet bjetcan_smeared = smearer_j_CMS_Benj(bjetcan.e(), bjetcan.perp(), bjetcan.pz(),
							      bjetcan.phi(), bjetcan.eta());
      //cout << "thebJets[4][jj] = " << thebJets[4][jj] << endl;
      //get the btagscore for the non-b-jet
      double btagscore = 0;
      //check if within tagging region
      if(fabs(bjetcan_smeared.eta()) < 2.5 && bjetcan_smeared.perp() > 20 && bjetcan_smeared.perp() < 1000) { btagscore = jet_btagscore(bjetcan_smeared.perp(), fabs(bjetcan_smeared.eta()), JetBtagDeepFlavB_B_Distributions); }
      //add it to the MyInfo Class for the PseudoJet:
      bjetcan_smeared.set_user_info(new MyInfo(btagscore, 1));
      bjetcan_smeared.set_user_index(thebJets[4][jj]);
      //cout << "bjetcan_smeared.user_info<MyInfo>().btagscore()= " << bjetcan_smeared.user_info<MyInfo>().btagscore() << endl;
      if(bjetcan.perp() > cut_pt_bjet && fabs(bjetcan.eta()) < cut_eta_bjet && jet_efficiency_accept(bjetcan)) { 
	bJets_unsort.push_back(bjetcan_smeared);
	//Set the charge
	bJets_unsort[bJets_unsort.size()-1].set_user_index(thebJets[4][jj]);
	//cout << "bJets_unsort[bJets_unsort.size()-1].user_index() = " << bJets_unsort[bJets_unsort.size()-1].user_index() << endl;
      }
    }
  
    for(int jj = 0; jj < numPhotons; jj++) {
      photoncan = fastjet::PseudoJet(thePhotons[1][jj], thePhotons[2][jj], thePhotons[3][jj], thePhotons[0][jj]);
	 Photons_unsort.push_back(smear_photon(photoncan));
     }
     for(int jj = 0; jj < numElectrons; jj++) {
       fastjet::PseudoJet Electron = fastjet::PseudoJet(theElectrons[1][jj], theElectrons[2][jj], theElectrons[3][jj], theElectrons[0][jj]);
       if(lepton_efficiency_accept(Electron, 11)) {
	 Electron.set_user_index(11);
	 Electrons_unsort.push_back(smear_lepton(Electron, 11));
	 Leptons_unsort.push_back(smear_lepton(Electron, 11));
       }
     }
     for(int jj = 0; jj < numPositrons; jj++) {
       fastjet::PseudoJet Positron = fastjet::PseudoJet(thePositrons[1][jj], thePositrons[2][jj], thePositrons[3][jj], thePositrons[0][jj]);
       if(lepton_efficiency_accept(Positron, 11)) {
	 Positron.set_user_index(-11);
	 Positrons_unsort.push_back(smear_lepton(Positron, -11));
	 Leptons_unsort.push_back(smear_lepton(Positron, -11));
       }
     }
     for(int jj = 0; jj < numMuons; jj++) {
       fastjet::PseudoJet Muon = fastjet::PseudoJet(theMuons[1][jj], theMuons[2][jj], theMuons[3][jj], theMuons[0][jj]);
       if(lepton_efficiency_accept(Muon, 13)) {
	 Muon.set_user_index(13);
	 Muons_unsort.push_back(smear_lepton(Muon, 13));
	 Leptons_unsort.push_back(smear_lepton(Muon, 13));
       }
     }
     for(int jj = 0; jj < numantiMuons; jj++) {
       fastjet::PseudoJet antiMuon = fastjet::PseudoJet(theantiMuons[1][jj], theantiMuons[2][jj], theantiMuons[3][jj], theantiMuons[0][jj]);
       if(lepton_efficiency_accept(antiMuon, 13)) { 
	 antiMuon.set_user_index(-13);
	 AntiMuons_unsort.push_back(smear_lepton(antiMuon, -13));
	 Leptons_unsort.push_back(smear_lepton(antiMuon, -13));
       }
     }


     /*
      * SORT RECONSTRUCTED OBJETS BY PT
      */
     vector<fastjet::PseudoJet> bJets, Jets, Electrons, Positrons, Muons, AntiMuons, Photons, Leptons, cJets, LightJets;

     //create a std::vector with all jets
     std::vector<fastjet::PseudoJet> allJets;
     allJets.reserve(Jets_unsort.size() + bJets_unsort.size());
     allJets.insert(allJets.end(), Jets_unsort.begin(), Jets_unsort.end());
     allJets.insert(allJets.end(), bJets_unsort.begin(), bJets_unsort.end());
     sort_by_btag_desc(allJets);

     //test: loop through allJets and print their btagscore and whether they are true b-jets
     /*for(int aj = 0; aj < allJets.size(); aj++) {
       cout << aj << "\t" << allJets[aj].user_info<MyInfo>().btagscore() << "\t" <<  allJets[aj].user_info<MyInfo>().trueb() << endl;
       }*/

     Jets = sorted_by_pt(Jets_unsort);
     bJets = sorted_by_pt(bJets_unsort);
     Photons = sorted_by_pt(Photons_unsort);
     numJets = Jets.size();
     numbJets = bJets.size();
     numPhotons = Photons.size();
     
     Electrons = sorted_by_pt(Electrons_unsort);
     Positrons = sorted_by_pt(Positrons_unsort);
     Muons = sorted_by_pt(Muons_unsort);
     AntiMuons = sorted_by_pt(AntiMuons_unsort);
     Leptons = sorted_by_pt(Leptons_unsort);

     numElectrons = Electrons.size();
     numPositrons = Positrons.size();
     numMuons = Muons.size();
     numantiMuons = AntiMuons.size();

     h_numbjets.thfill(numbJets);
     
     int numLeptons = numElectrons + numPositrons + numMuons + numantiMuons;


     //fill in the input weight:
     total_weight_in += evweight;
     total_event_in++;
     
     fastjet::PseudoJet ETmiss = fastjet::PseudoJet(theETmiss[1], theETmiss[2], theETmiss[3], theETmiss[0]);

     /*
      * FILL IN THE HISTOGRAMS
      */
     for(int jj = 0; jj < numJets; jj++) {
       h_pT_jets.thfill(Jets[jj].perp(), evweight);
     }

     for(int jj = 0; jj < numLeptons; jj++) {
       h_pT_leptons.thfill(Leptons[jj].perp(), evweight);
     }

     /* 
      * CUTS START HERE
      */

     //Fill in the number of b-jets and number of photons until we reach two of each
     fastjet::PseudoJet bJet1, bJet2, bJet3, bJet4, bJet5, bJet6;
     
     //select the six highest-pT b-jets
     //if(numbJets >= 6 && passed_all_cuts) { //change this to allJets
     if(allJets.size() >= 6 && passed_all_cuts) {
       /*bJet1 = bJets[0];
       bJet2 = bJets[1];
       bJet3 = bJets[2];
       bJet4 = bJets[3];
       bJet5 = bJets[4];
       bJet6 = bJets[5];*/
       bJet1 = allJets[0];
       bJet2 = allJets[1];
       bJet3 = allJets[2];
       bJet4 = allJets[3];
       bJet5 = allJets[4];
       bJet6 = allJets[5];
       /*cout << "bJets:" << endl;
       cout << bJet1 << endl;
       cout << bJet2 << endl;
       cout << bJet3 << endl;
       cout << bJet4 << endl;
       cout << bJet5 << endl;
       cout << bJet6 << endl;*/

       //evweight *= btag_weight(bJets[0], 1, 0) * btag_weight(bJets[1], 1, 0) * btag_weight(bJets[2], 1, 0) * btag_weight(bJets[3], 1, 0) *btag_weight(bJets[4], 1, 0) * btag_weight(bJets[5], 1, 0);
     } else continue;
     pass_6b+=evweight;

      //check that the distance between any two b-jets is larger than DeltaR = 0.3
     bool dRbbfail = false;
     for(int b1=0; b1 < 6; b1++) {
       for(int b2=0; b2 < 6; b2++) {
	 if(b1!=b2) { 
	   //double dRbb = deltaR(bJets[b1],bJets[b2]);//change to allJets
	   double dRbb = deltaR(allJets[b1],allJets[b2]);
	   if(dRbb < cut_dRbbmin) dRbbfail = true;
	 }
       }
     }
     if(dRbbfail == false && passed_all_cuts) {
       pass_drbb+=evweight;
     }
     else {
       passed_all_cuts = false;
     }


     //check that the pT of the b-jets is larger than cut_pt_bjet_1/2/3/4 and 5/6 larger than cut_pt_bjet
     if(bJet1.perp() > cut_pt_bjet1 && bJet2.perp() > cut_pt_bjet2 && bJet3.perp() > cut_pt_bjet3 && bJet4.perp() > cut_pt_bjet4 && bJet5.perp() > cut_pt_bjet && bJet6.perp() > cut_pt_bjet && passed_all_cuts) {
       pass_ptb+=evweight;
     }
     else {
       passed_all_cuts = false;
     }

     /* At this point:
      * we are only considering events with at least 6 SINGLY b-tagged jets 
      * Note that: there will be events with 4/2 SINGLY b-tagged jets and 1/2 DOUBLY-CHARGED B-JETS, RESPECTIVELY
      * count each combination: all 6 b-jets identified with the same charge, then 5-1, 4-2, 3-3
      * THEN, FOR THE CASE OF THE 6 SINGLY b-tagged jets:
      * - 6-0: if all are identified as the same charge, perform the same analysis as before
      * - 5-1/4-2/3-3: allow combinations only of the one with opposite charge only in finding the best combination
      */
     int charge_sum(0);
     int charge_sum_abs(0);

     //get the initial charge sum:
     for(int bc=0; bc < 6; bc++) {
       charge_sum += allJets[bc].user_index(); 
       charge_sum_abs += abs(allJets[bc].user_index());
     }
     charge_sum = abs(charge_sum);
     h_chargesum_initial.thfill(charge_sum,evweight);
     h_chargesumabs_initial.thfill(charge_sum_abs,evweight);

     //if(charge_sum_abs != 6 || charge_sum != 0) { continue; }
     //if(abs(allJets[0].user_index()) != 1 || abs(allJets[1].user_index()) !=1 || abs(allJets[2].user_index()) != 1 || abs(allJets[3].user_index()) != 1 || abs(allJets[4].user_index()) != 1 || abs(allJets[5].user_index()) != 1) { continue; }
     
     /* RANDOMLY CHANGE THE CHARGE OF A B-JET WITH A PROBABILITY 1-Pb_charge */
     for(int bc=0; bc < 6; bc++) {
       if(rnd.Rndm() > Pb_charge) {
	 allJets[bc].set_user_index( -1*allJets[bc].user_index() ); 
       }
     }
     // recalculate charge sum after flips
     charge_sum = 0;
     charge_sum_abs = 0;
     for(int bc=0; bc < 6; bc++) {
       charge_sum += allJets[bc].user_index();
       charge_sum_abs += fabs(allJets[bc].user_index());
     }
     //cout << bJets[0].user_index() << "\t" << bJets[1].user_index() << "\t" << bJets[2].user_index() << "\t" << bJets[3].user_index() << "\t" << bJets[4].user_index() << "\t" << bJets[5].user_index() << endl << endl;
     // we don't care if positive or negative, so just take the absolute value
     charge_sum = abs(charge_sum);
     h_chargesum_final.thfill(charge_sum,evweight);
     h_chargesumabs_final.thfill(charge_sum_abs,evweight);


     //loop over the 15 possible pairings of the 6 b-jets and calculate invariant mass for each "higgs boson candiate"
     double chisq_min(1E99);
     int mincombo(-1);
     std::vector<double> mbb; //reco higgs masses for optimal combo
     std::vector<fastjet::PseudoJet> ph; //reco higgses for optimal combo
     for(int pp = 0; pp < pairs_of_six.size(); pp++) {
       fastjet::PseudoJet bb1;
       fastjet::PseudoJet bb2;
       fastjet::PseudoJet bb3;
       if(chargetagging) {
	 /* check whether the current charge combination should be considered or not 
	    in the minimization of the chisq variable. skip if not compatible with charge
	    if all the jets are identified to have the same charge or we have the 5/1 case, 
	    we just associate them randomly as before.
	    The real difference comes in when we have 4/2 or 3/3 and then we want to avoid having combos
	    of same-charge b-jets. 
	 */
	 //first calculate the charges of each combination: 
	 int charge_combo[3]= {allJets[pairs_of_six[pp][0]].user_index() + allJets[pairs_of_six[pp][1]].user_index(), allJets[pairs_of_six[pp][2]].user_index() + allJets[pairs_of_six[pp][3]].user_index(), allJets[pairs_of_six[pp][4]].user_index() + allJets[pairs_of_six[pp][5]].user_index()};
	 /* then ensure that the sum of the absolute values is minimized: 
	    for +++--- (3/3) this should be zero
	    for 4/2 this should be 2
	    for 5/1 this should be 4 (this is always the case, so the loop doesn't do anything here)
	 */
	 //cout << "abs(charge_combo[0]) + abs(charge_combo[1]) + abs(charge_combo[2]) = " <<  abs(charge_combo[0]) + abs(charge_combo[1]) + abs(charge_combo[2]) << " charge_sum = " << charge_sum << endl;
	 if( abs(charge_combo[0]) + abs(charge_combo[1]) + abs(charge_combo[2]) != charge_sum) continue;
       }
       bb1 = allJets[pairs_of_six[pp][0]] + allJets[pairs_of_six[pp][1]]; 
       bb2 = allJets[pairs_of_six[pp][2]] + allJets[pairs_of_six[pp][3]];
       bb3 = allJets[pairs_of_six[pp][4]] + allJets[pairs_of_six[pp][5]];
       /*cout << bJets[pairs_of_six[pp][0]] << endl;
       cout << bJets[pairs_of_six[pp][1]] << endl;
       cout << bb1 << endl;*/
       double mbb1 = bb1.m();
       double mbb2 = bb2.m();
       double mbb3 = bb3.m();
       //for(int ps = 0; ps < pairs_of_six[pp].size(); ps++) cout << pairs_of_six[pp][ps] << ",";
       double chisq_combo = sqrt(pow ((mbb1-mhiggs1), 2)+ pow( (mbb2-mhiggs2), 2) + pow((mbb3-mhiggs3),2)); 
       //cout << "\tchisq = " << chisq_combo << endl;
       if(chisq_combo < chisq_min) { mbb.clear(); chisq_min = chisq_combo; mincombo = pp; mbb.push_back(fabs((mbb1-mhiggs1))); mbb.push_back(fabs((mbb2-mhiggs2))); mbb.push_back(fabs((mbb3-mhiggs3))); ph.push_back(bb1); ph.push_back(bb2); ph.push_back(bb3); }
     }
     //no combo found
     if(mincombo == -1) continue;

     //cout << "min. chi-sq = " << chisq_min << " for combo " << mincombo << endl;
     std::vector<double> DeltaM_unsort; //store the differences of the optimal combination with the Higgs mass in ascending order (DeltaM_min, DeltaM_med, DeltaM_max)
     for(int m = 0; m < 3; m++) { DeltaM_unsort.push_back(mbb[m]); /*cout << "mbb = " << mbb[m] << endl;*/ }

     //DeltaM_index contains the indices of DeltaM_unsort in assending order
     //so: DeltaM_unsort[DeltaM_index[0]] is DeltaM_min, 1 is DeltaM_med and 2 is DeltaM_max
     std::vector<int> DeltaM_index(DeltaM_unsort.size());
     std::size_t n(0);
     std::generate(std::begin(DeltaM_index), std::end(DeltaM_index), [&]{ return n++; });
     std::sort(  std::begin(DeltaM_index), std::end(DeltaM_index), [&](int i1, int i2) { return DeltaM_unsort[i1] < DeltaM_unsort[i2]; } );
     ph = sorted_by_pt(ph);


     //impose further cuts:

     //chi-sq cut:
     if(chisq_min < cut_chisq_min && passed_all_cuts) {
       pass_chisq+=evweight;
     }
     else {
       passed_all_cuts = false;
     }
       
     //DeltaM cuts:
     /*cout << "DeltaM:" << endl;
       cout << DeltaM_unsort[DeltaM_index[0]] << "\t" << DeltaM_unsort[DeltaM_index[1]] << "\t" << DeltaM_unsort[DeltaM_index[2]] << endl;*/
     if(DeltaM_unsort[DeltaM_index[0]] < cut_DeltaM_min && DeltaM_unsort[DeltaM_index[1]] < cut_DeltaM_med && DeltaM_unsort[DeltaM_index[2]] < cut_DeltaM_max && passed_all_cuts) {
       pass_DeltaM+=evweight;
     }
     else {
       passed_all_cuts = false;
     }
     
     //pT of reconstructed Higgs bosons:
     if(ph[0].perp() > cut_pt_higgs1 && ph[1].perp() > cut_pt_higgs2 && ph[2].perp() > cut_pt_higgs3 && passed_all_cuts) {
       pass_pthiggses+=evweight;
     }
     else {
       passed_all_cuts = false;
     }

     //delta R(bb) between b's in reco higgses:
     if(deltaR(allJets[pairs_of_six[mincombo][0]], allJets[pairs_of_six[mincombo][1]]) < cut_dR_hbbreco && deltaR(allJets[pairs_of_six[mincombo][2]], allJets[pairs_of_six[mincombo][3]]) < cut_dR_hbbreco && deltaR(allJets[pairs_of_six[mincombo][4]], allJets[pairs_of_six[mincombo][5]]) < cut_dR_hbbreco && passed_all_cuts) { pass_dRbbhiggses+=evweight;
     }
     else {
       passed_all_cuts = false;
     }

     //delta R between Higgs bosons:
     if(deltaR(ph[0], ph[1]) < cut_dR_higgses && deltaR(ph[0], ph[2]) < cut_dR_higgses && deltaR(ph[1], ph[2]) < cut_dR_higgses && passed_all_cuts) {
       pass_dRhiggses+=evweight;
     }
     else {
       passed_all_cuts = false;
     }

     double m6b = (allJets[0]+allJets[1]+allJets[2]+allJets[3]+allJets[4]+allJets[5]).m();
     if(m6b > cut_m6b && passed_all_cuts) {
       pass_m6b+=evweight;
     }
     else {
       passed_all_cuts = false;
     }
     
     /*
     * DOES THE EVENT PASS ALL THE CUTS?
     * IF SO INCREMENT THE WEIGHT
     */
     if(passed_all_cuts) { 
       passcuts+=evweight;
       eventcount++;
     }

     /* 
      * calculate variables for the _var.root file and plot:
      */

     
     
     /*
      * Fill in the _var.root file for further analysis.
      */
     variables[0] = evweight;
     variables[1] = bJet1.perp();
     variables[2] = bJet2.perp();
     variables[3] = bJet3.perp();
     variables[4] = bJet4.perp();
     variables[5] = bJet5.perp();
     variables[6] = bJet6.perp();
     variables[7] = m6b;
     variables[8] = chisq_min;
     variables[9] = DeltaM_unsort[DeltaM_index[0]];
     variables[10] = DeltaM_unsort[DeltaM_index[1]];
     variables[11] = DeltaM_unsort[DeltaM_index[2]];
     variables[12] = ph[0].perp();
     variables[13] = ph[1].perp();
     variables[14] = ph[2].perp();
     variables[15] = deltaR(ph[0], ph[1]);
     variables[16] = deltaR(ph[0], ph[2]);
     variables[17] = deltaR(ph[1], ph[2]);
     variables[18] = deltaR(allJets[pairs_of_six[mincombo][0]], allJets[pairs_of_six[mincombo][1]]);
     variables[19] = deltaR(allJets[pairs_of_six[mincombo][2]], allJets[pairs_of_six[mincombo][3]]);
     variables[20] = deltaR(allJets[pairs_of_six[mincombo][4]], allJets[pairs_of_six[mincombo][5]]);
       

     weight = evweight;
     
     Data2->Fill();

     
     /* fill in Histograms: */
     h_pT_b.thfill(allJets[0].perp());
     h_pT_b.thfill(allJets[1].perp());
     h_pT_b.thfill(allJets[2].perp());
     h_pT_b.thfill(allJets[3].perp());
     h_pT_b.thfill(allJets[4].perp());
     h_pT_b.thfill(allJets[5].perp());
     h_pT_b1.thfill(allJets[0].perp());
     h_pT_b2.thfill(allJets[1].perp());
     h_pT_b3.thfill(allJets[2].perp());
     h_pT_b4.thfill(allJets[3].perp());
     h_pT_b5.thfill(allJets[4].perp());
     h_pT_b6.thfill(allJets[5].perp());
     h_chisq.thfill(chisq_min);
     h_DeltaM_min.thfill(DeltaM_unsort[DeltaM_index[0]]);
     h_DeltaM_med.thfill(DeltaM_unsort[DeltaM_index[1]]);
     h_DeltaM_max.thfill(DeltaM_unsort[DeltaM_index[2]]);
     h_pT_h1.thfill(ph[0].perp());
     h_pT_h2.thfill(ph[1].perp());
     h_pT_h3.thfill(ph[2].perp());
     h_pT_dRhh.thfill(deltaR(ph[0], ph[1]));
     h_pT_dRhh.thfill(deltaR(ph[0], ph[2]));
     h_pT_dRhh.thfill(deltaR(ph[1], ph[2]));
     h_m6b.thfill(m6b);
     
     /* IF EVENT HAS PASSED CUTS
      * PRINT TO .evp or .evp2 FILE 
      * INCREMENT AND CONTINUE
      */
     outevp << ii << endl;
     
  } /* LOOP OVER EVENTS ENDS HERE 
     * ENDS HERE
     */
  Data2->GetCurrentFile();
  Data2->Write();
  dat2->Close();
  cout << "A root tree has been written to the file: " << fnameroot << endl;
		
   /* OUTPUT HISTOGRAMS
		 * HERE AND 
		 * FINISH
   */
  h_dummy.plot(output,1,0);
		// "inclusive plot"
  h_pT_jets.add(output,1,0);
  h_pT_leptons.add(output,1,0);
  h_pT_b.add(output,1,0);
  h_pT_b1.add(output,1,0);
  h_pT_b2.add(output,1,0);
  h_pT_b3.add(output,1,0);
  h_pT_b4.add(output,1,0);
  h_pT_b5.add(output,1,0);
  h_pT_b6.add(output,1,0);
  h_chisq.add(output,1,0);
  h_DeltaM_min.add(output,1,0);
  h_DeltaM_med.add(output,1,0);
  h_DeltaM_max.add(output,1,0);
  h_pT_h1.add(output,1,0);
  h_pT_h2.add(output,1,0);
  h_pT_h3.add(output,1,0);
  h_pT_dRhh.add(output,1,0);
  h_m6b.add(output,1,0);
  h_chargesumabs_initial.add(output,1,0);
  h_chargesum_initial.add(output,1,0);
  h_chargesumabs_final.add(output,1,0);
  h_chargesum_final.add(output,1,0);
  h_numbjets.add(output,1,0);
  
  cout << "------------------" << endl;
  cout << "total weight in =\t\t\t\t\t\t" <<  total_weight_in << endl;
  cout << "total MC events in =\t\t\t\t\t\t" << total_event_in << endl;
  cout << "------------------" << endl;
  cout << "cuts/counters:" << endl;
  cout << "6bs:\t\t\t\t\t\t\t\t" <<  pass_6b << endl;
  cout << "6bs with dR(b,b) > " << cut_dRbbmin << "\t\t\t\t\t\t" << pass_drbb << endl;
  cout << "6bs with pT > [" << cut_pt_bjet1 << ", " << cut_pt_bjet2 << ", " << cut_pt_bjet3 << ", " << cut_pt_bjet4 << ", " << cut_pt_bjet << ", " << cut_pt_bjet <<  "]\t\t\t\t" << pass_ptb << endl;
  cout << "chisq minimum < " << cut_chisq_min << "\t\t\t\t\t\t" << pass_chisq << endl;
  cout << "DeltaM(min,med,max) < [" << cut_DeltaM_min << ", " << cut_DeltaM_med << ", " << cut_DeltaM_max << "]\t\t\t\t" << pass_DeltaM << endl;
  cout << "Three reco Higgses with pT > [" << cut_pt_higgs1 << ", " << cut_pt_higgs2 << ", " << cut_pt_higgs3 << "]\t\t\t" << pass_pthiggses << endl;
  cout << "DeltaR(b,b) in reco Higgses < " << cut_dR_hbbreco << "\t\t\t\t" << pass_dRbbhiggses << endl;
  cout << "dR between reco Higgses < " << cut_dR_higgses << "\t\t\t\t\t" << pass_dRhiggses << endl;
  cout << "m6b > " << cut_m6b << "\t\t\t\t\t\t\t" << pass_m6b << endl;
  cout << "------------------" << endl;
  cout << "total weight out =\t\t\t\t\t\t" <<  passcuts << endl;
  cout << "actual MC events = \t\t\t\t\t\t" << eventcount << endl;
  cout << "efficiency =\t\t\t\t\t\t\t" <<  passcuts/total_weight_in << endl;
  cout << "------------------" << endl;
  outdat << passcuts/total_weight_in << endl;
  return 0;
}


//check if an integer is in the vector<int>:
bool is_in(int a, vector<int> intvec) {
  bool is_it_in = 0;
  for(int iv = 0; iv < intvec.size(); iv++) {
    if(intvec[iv] == a) {
      is_it_in = 1;
    }
  }
  return is_it_in;
}


double analyze_event(fastjet::PseudoJet photon1, fastjet::PseudoJet photon2, fastjet::PseudoJet cjet, fastjet::PseudoJet bjet, fastjet::PseudoJet lepton, fastjet::PseudoJet etmiss, double evweight_i) {
  

  //construct all relevant variables:
  

  //passed_mcevents++;


  return evweight_i;
}


double Pb_to_b(double pt) {
  if(perfect_tagging) { return 1.0; }
  return 1.00;//0.75;
}

double Pb_to_c(double pt) {
  if(perfect_tagging) { return 0.; } 
  return 0.125; 
}

double Pc_to_c(double pt) {
  if(perfect_tagging) { return 1.0; }
  return 0.2;
}

double Pc_to_b(double pt) {
  if(perfect_tagging) { return 0.0; }
  return 0.1;
}

double Pjet_to_b(double pt) {
  if(perfect_tagging) { return 0.0; }
  return 0.01;
}

double Pjet_to_c(double pt) {
  if(perfect_tagging) { return 0.0; }
  return 0.005;
}

// jet to lepton mistag prob
double Pjet_to_photon(double pt) {
  if(perfect_tagging) { return 0.0; }
  double pval(0.);
  double alpha = 0.01;
  double beta = 1/30.0;

  pval = alpha*exp(-beta * pt);

  return pval;
}

double ctag_weight(fastjet::PseudoJet jet, bool btag, bool ctag) {
  double weight = 0;
  if(btag) { 
    weight = Pb_to_c(jet.perp());
    return weight; 
  }
  else if(ctag) { 
    weight = Pc_to_c(jet.perp());
    return weight;
  }
  else { // light jet
    weight = Pjet_to_c(jet.perp());
    return weight; 
  }
  return weight; 
}


double btag_weight(fastjet::PseudoJet jet, bool btag, bool ctag) {
  double weight = 0;
  if(btag) { 
    weight = Pb_to_b(jet.perp());
    return weight; 
  }
  else if(ctag) { 
    weight = Pc_to_b(jet.perp());
    return weight;
  }
  else { // light jet
    weight = Pjet_to_b(jet.perp());
    return weight; 
  }
  return weight; 
}

double atag_weight(fastjet::PseudoJet jet, bool btag, bool ctag) {
  double weight = 0;
  if(btag) { 
    weight = 0.;
  }
  else if(ctag) { 
    weight = 0.;
  }
  else {
    weight = Pjet_to_photon(jet.perp());
  }
  return weight; 
}


double dot(fastjet::PseudoJet p1, fastjet::PseudoJet p2) {
  return (p1.e() * p2.e() - p1.px() * p2.px() - p1.py() * p2.py() - p1.pz() * p2.pz() );
}


double deltaR(fastjet::PseudoJet p1, fastjet::PseudoJet p2) { 
  double dphi_tmp; 

  dphi_tmp = p2.phi() - p1.phi();
  if(dphi_tmp > M_PI) 
    dphi_tmp = 2 * M_PI - dphi_tmp;
  else if( dphi_tmp < - M_PI)  
    dphi_tmp = 2 * M_PI + dphi_tmp;
  
  //  return sqrt(sqr(p1.eta() - p2.eta()) + sqr(dphi_tmp));
  return sqrt(sqr(p1.rap() - p2.rap()) + sqr(dphi_tmp));
}

//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  ostr << "e, pt, y, phi =" 
       << " " << setw(10) <<  jet.e()  
       << " " << setw(10) << jet.perp() 
       << " " << setw(6) <<  jet.rap()  
       << " " << setw(6) <<  jet.phi()  
       << ", mass = " << setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;
}
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


bool btag_hadrons(fastjet::PseudoJet jet) {
  bool btagged(false);
  /* search constintuents of jets for b-mesons */
  for(int cc = 0; cc < jet.constituents().size(); cc++) { 
    for(int bb = 0; bb < 105; bb++) { 
      if(jet.constituents()[cc].user_index() == bhadronid[bb]) { 
	btagged = true;
	//	cout << "Jet B-tagged!" << endl;
	//	cout << jet << endl;
      }
    }
  }
  return btagged;
}

std::pair<complex_ld, complex_ld> pnuz_fromw(fastjet::PseudoJet plepton, fastjet::PseudoJet pmiss) {
  
  double l = sqrt(sqr(plepton.px()) + sqr(plepton.py()) + sqr(plepton.pz()));
  double A = pow(mw,2)/2. + plepton.px()*pmiss.px() + plepton.py()*pmiss.py();
  double zl = plepton.pz();
  double a = pow(l,2) - pow(zl,2);
  double b = - 2 * A * zl;
  double c = - pow(A,2) + pow(l,2) * pow(pmiss.perp(),2);
  
  std::pair<complex_ld, complex_ld> result = quadsolve(a,b,c);
  return result;
}

std::pair<fastjet::PseudoJet, fastjet::PseudoJet> get_Wvectors(fastjet::PseudoJet plepton, fastjet::PseudoJet pmiss) { 
  std::pair<fastjet::PseudoJet, fastjet::PseudoJet> wvecs;

  //get the two solutions for the pnuz
  std::pair<complex_ld, complex_ld> pnuz_fromw_res = pnuz_fromw(plepton, pmiss);

  //reconstruct the two ws
  double Ew1 = sqrt( pow((plepton.px()+pmiss.px()),2) + pow(plepton.py()+pmiss.py(),2) + pow(plepton.pz() + real(pnuz_fromw_res.first),2) + pow(mw,2) );
  double Ew2 = sqrt( pow((plepton.px()+pmiss.px()),2) + pow(plepton.py()+pmiss.py(),2) + pow(plepton.pz() + real(pnuz_fromw_res.second),2) + pow(mw,2) );
    
  fastjet::PseudoJet wvec1 = fastjet::PseudoJet( plepton.px()+pmiss.px(), plepton.py()+pmiss.py(), plepton.pz() + real(pnuz_fromw_res.first), Ew1 );
  fastjet::PseudoJet wvec2 = fastjet::PseudoJet( plepton.px()+pmiss.px(), plepton.py()+pmiss.py(), plepton.pz() + real(pnuz_fromw_res.second), Ew2 );
  wvecs.first = wvec1;
  wvecs.second = wvec2;
  return wvecs; 
}


std::pair<fastjet::PseudoJet, fastjet::PseudoJet> get_Wvectors_GTX(fastjet::PseudoJet plepton, fastjet::PseudoJet pmiss) { 
  std::pair<fastjet::PseudoJet, fastjet::PseudoJet> wvecs;

  //get the two solutions for the pnuz
  std::pair<complex_ld, complex_ld> pnuz_fromw_res = pnuz_fromw(plepton, pmiss);

  fastjet::PseudoJet pnu1 = fastjet::PseudoJet( pmiss.px(), pmiss.py(), real(pnuz_fromw_res.first), sqrt( sqr(pmiss.px()) + sqr(pmiss.py()) + sqr(real(pnuz_fromw_res.first)) ) );
  fastjet::PseudoJet pnu2 = fastjet::PseudoJet( pmiss.px(), pmiss.py(), real(pnuz_fromw_res.second), sqrt( sqr(pmiss.px()) + sqr(pmiss.py()) + sqr(real(pnuz_fromw_res.second)) ) );

  fastjet::PseudoJet wvec1 = pnu1 + plepton;
  fastjet::PseudoJet wvec2 = pnu2 + plepton;
    
  wvecs.first = wvec1;
  wvecs.second = wvec2;
  return wvecs; 
}
	



fastjet::PseudoJet smear_jet(fastjet::PseudoJet jet_in) {
  if(donotsmear_jets) { return jet_in; }

  fastjet::PseudoJet smeared; 
  double smearing = 20, smeared_pt(0);

  double pt = jet_in.perp();
  double eta = fabs(jet_in.eta());
  double sigma(0);
  
  //double a, b, S, C;
  //if(eta < 0.8) { a = 3.2; b = 0.07; S = 0.74; C = 0.05; }
  //if(eta > 0.8 && eta < 1.2) { a = 3.0; b = 0.07; S = 0.81; C = 0.05; }
  //if(eta > 1.2 && eta < 2.8) { a = 3.3; b = 0.08; S = 0.54; C = 0.05; }
  //if(eta > 2.8 /*&& eta < 3.6*/) { a = 2.8; b = 0.11; S = 0.83; C = 0.05; }

  //double mu_pileup = 40;
  //double N = a + b * mu_pileup;

  //sigma = pt * sqrt( sqr(N)/sqr(pt) + sqr(S) / pt + sqr(C) );*/

  //MRM:
  sigma = 1.0 * sqrt(pt);
  
  smeared_pt = fabs(rnd.Gaus(0,sigma));
  double theta = rnd.Rndm()*M_PI;
  double phi = rnd.Rndm()*2.*M_PI;

  
  double deltaE = - jet_in.e() + sqrt( sqr(jet_in.e()) + sqr(smeared_pt) + 2 * (smeared_pt*sin(theta)*cos(phi)*jet_in.px() + smeared_pt*sin(theta)*sin(phi)*jet_in.py() + smeared_pt*cos(theta)*jet_in.pz()));

  fastjet::PseudoJet smearing_vector(smeared_pt*sin(theta)*cos(phi),smeared_pt*sin(theta)*sin(phi), smeared_pt*cos(theta), deltaE);
  
  smeared = jet_in + smearing_vector;  
  
  return smeared;
}

fastjet::PseudoJet smear_photon(fastjet::PseudoJet photon_in) {
  if(donotsmear_photons) { return photon_in; }

  fastjet::PseudoJet smeared;
  double smeared_pt = 0;
  //double smear_frac = 0.1E-2;
  //double smear_sampling = 0.15;
  //FROM MRM:
  double smear_frac = 0.17E-2;
  double smear_sampling = 0.20;
  double sigma(smear_sampling * sqrt(photon_in.perp()) + smear_frac*photon_in.perp());

  smeared_pt = fabs(rnd.Gaus(0,sigma));
  double theta = rnd.Rndm()*M_PI;
  double phi = rnd.Rndm()*2.*M_PI;

  //cout << smeared_pt*sin(theta)*cos(phi) <<  " " << smeared_pt*sin(theta)*sin(phi) << " " << smeared_pt*cos(theta) << " " << smeared_pt << endl;
  //cout << smeared_pt*sin(theta)*cos(phi) << " " << smeared_pt*sin(theta)*sin(phi) << "  " <<  smeared_pt*cos(theta) << endl;

  
  double deltaE = - photon_in.e() + sqrt( sqr(photon_in.e()) + sqr(smeared_pt) + 2 * (smeared_pt*sin(theta)*cos(phi)*photon_in.px() + smeared_pt*sin(theta)*sin(phi)*photon_in.py() + smeared_pt*cos(theta)*photon_in.pz()));
  
  fastjet::PseudoJet smearing_vector(smeared_pt*sin(theta)*cos(phi),smeared_pt*sin(theta)*sin(phi), smeared_pt*cos(theta), deltaE);
  smeared = photon_in + smearing_vector;
  //cout << "smeared mass = " << smeared.m() << endl;
  return smeared;

}



fastjet::PseudoJet smear_lepton(fastjet::PseudoJet lepton_in, int lepton_id) {

  if(donotsmear_leptons) { return lepton_in; }
   
    
  fastjet::PseudoJet smeared;
  double smeared_pt = 0;
  double smearing = 20.;

  double pt = lepton_in.perp();
  double lepton_energy = lepton_in.e();
  double eta = fabs(lepton_in.eta());
  double sigma(0);

  //see ATL-PHYS-PUB-2013-009
  if(lepton_id == 13) {
    double sigma_id = 0;
    double sigma_ms = 0;
    double sigma_cb = 0;
    double a1, a2, b0, b1, b2;
    
    if(eta < 0.18) { a1 = 0.01061; a2 = 0.000157; }
    if(eta > 0.18 && eta < 0.36) { a1 = 0.01084; a2 = 0.000153; }
    if(eta > 0.36 && eta < 0.54) { a1 = 0.01124; a2 = 0.000150; }
    if(eta > 0.54 && eta < 0.72) { a1 = 0.01173; a2 = 0.000149; }
    if(eta > 0.72 && eta < 0.90) { a1 = 0.01269; a2 = 0.000148; }
    if(eta > 0.90 && eta < 1.08) { a1 = 0.01406; a2 = 0.000161; }
    if(eta > 1.08 && eta < 1.26) { a1 = 0.01623; a2 = 0.000192; }
    if(eta > 1.26 && eta < 1.44) { a1 = 0.01755; a2 = 0.000199; } 
    if(eta > 1.44 && eta < 1.62) { a1 = 0.01997; a2 = 0.000232; } 
    if(eta > 1.62 && eta < 1.80) { a1 = 0.02453; a2 = 0.000261; }
    if(eta > 1.80 && eta < 1.98) { a1 = 0.03121; a2 = 0.000297; }
    if(eta > 1.98 && eta < 2.16) { a1 = 0.03858; a2 = 0.000375; }
    if(eta > 2.16 && eta < 2.34) { a1 = 0.05273; a2 = 0.000465; }
    if(eta > 2.34 && eta < 2.52) { a1 = 0.05329; a2 = 0.000642; }
    if(eta > 2.52 /*&& eta < 2.70*/) { a1 = 0.05683; a2 = 0.000746; }

    if(eta < 1.05) { b1 = 0.02676; b2 = 0.00012; }
    if(eta > 1.05) { b1 = 0.03880; b2 = 0.00016; }

    sigma_id = pt * sqrt( a1 + sqr(a2 * pt) );
    sigma_ms = pt * sqrt( sqr(b0/pt) + sqr(b1) + sqr(b2*pt) );
    sigma = (sigma_id * sigma_ms)/sqrt( sqr(sigma_id) + sqr(sigma_ms) ); //sigma_cb

  }


  if(lepton_id == 11) {
    double sigma = 0;
    if(eta < 1.4) { sigma = sqrt( sqr(0.3) + sqr(0.10 * sqrt(lepton_energy)) + sqr( 0.010 * lepton_energy ) ); }
    if(eta > 1.4 /* && eta < 2.47 */) { sigma = sqrt( sqr(0.3) + sqr(0.15 * sqrt(lepton_energy)) + sqr( 0.015 * lepton_energy ) ); }
  }

  smeared_pt = fabs(rnd.Gaus(0,sigma));
  double theta = rnd.Rndm()*M_PI;
  double phi = rnd.Rndm()*2.*M_PI;

  double deltaE = - lepton_in.e() +  sqrt( sqr(lepton_in.e()) + sqr(smeared_pt) + 2 * (smeared_pt*sin(theta)*cos(phi)*lepton_in.px() + smeared_pt*sin(theta)*sin(phi)*lepton_in.py() + smeared_pt*cos(theta)*lepton_in.pz()));

  fastjet::PseudoJet smearing_vector(smeared_pt*sin(theta)*cos(phi),smeared_pt*sin(theta)*sin(phi), smeared_pt*cos(theta), deltaE);
  
  smeared = lepton_in + smearing_vector;  
  
  //smeared.reset(smeared.px(), smeared.py(), smeared.pz(), eprime);

  smeared.set_user_index(lepton_id);
  
  return smeared;
}

bool lepton_efficiency_accept(fastjet::PseudoJet lepton_in, int lepton_id) {
  bool accepted(1);
  if(donot_apply_efficiency) { return accepted; }

  double pt = lepton_in.perp();
  double eta = fabs(lepton_in.eta());
  
  double epsilon = 0;
  if(lepton_id == 11) {
    epsilon = 0.85 - 0.191 * exp(1 - pt/20);
  }
  if(lepton_id == 13) {
    if(eta<0.1) { epsilon = 0.54; }
    if(eta>0.1) { epsilon = 0.97; } 
  }
  double random_num = rnd.Rndm();
  //  cout << lepton_id << " " << pt << " " << eta << " " << random_num << " " << epsilon << endl;
  if(random_num > epsilon) { accepted = 0; }
  return accepted;
}
bool photon_efficiency_accept(fastjet::PseudoJet photon_in) {
  bool accepted(1);
  if(donot_apply_efficiency) { return accepted; }

  double pt = photon_in.perp();
  double eta = fabs(photon_in.eta());
  
  double epsilon = 0;
  
  epsilon = 0.76 - 1.98 * exp(-pt/16.1);
 
  double random_num = rnd.Rndm();
  if(random_num > epsilon) { accepted = 0; }
  return accepted;
}
bool jet_efficiency_accept(fastjet::PseudoJet jet_in) {
    bool accepted(1);
    if(donot_apply_efficiency) { return accepted; }

    double pt = jet_in.perp();
    double epsilon = 0;

    epsilon = 0.75 + (0.95 - 0.75) * pt / (50. - 20.);
    if(epsilon < 0) { epsilon = 0; }
    if(epsilon > 1.0) { epsilon = 1.0; }

    
    double random_num = rnd.Rndm();
    if(random_num > epsilon) { accepted = 0; }
    return accepted;
    
}

                     
                     
std::pair<complex_d, complex_d> quadsolve(complex_ld a, complex_ld b, complex_ld c) {
  std::pair<complex_d, complex_d> res;
  complex_ld D(0.);
  complex_ld ac = a*c;
  complex_ld fourac = complex_ld(4 * real(ac), 4 * imag(ac));
  D =  sqrt(b*b - fourac);
  complex_ld twoa = complex_ld( 2 * real(a), 2 * imag(a));
  res.first = (-b + D)/twoa;
  res.second = (-b - D)/twoa;

  return res;
}

// start is the current position in the list, advancing by 2 each time
// pass 0 as start when calling at the top level
void generatePairings(int* items, int itemcount, int start)
{
  vector<int> items_complete;
    if(itemcount & 1)
        return; // must be an even number of items
    // is this a complete pairing?
    if(start == itemcount)
    {
        // output pairings:
        int i;
	//items_complete.resize(itemcount);
        for(i = 0; i<itemcount; i+=2)
        {
	  printf("[%d, %d] ", items[i], items[i+1]);
	  items_complete.push_back(items[i]);
	  items_complete.push_back(items[i+1]);
        }
	//global!
	pairs_of_six.push_back(items_complete);
        //printf("\n");
        return;
    }

    // for the next pair, choose the first element in the list for the
    // first item in the pair (meaning we don't have to do anything 
    // but leave it in place), and each of the remaining elements for
    // the second item:
    int j;
    for(j = start+1; j<itemcount; j++)
    {
        // swap start+1 and j:
        int temp = items[start+1];
        items[start+1] = items[j];
        items[j] = temp;

        // recurse:
        generatePairings(items, itemcount, start+2);

        // swap them back:
        temp = items[start+1];
        items[start+1] = items[j];
        items[j] = temp;
    }

}


fastjet::PseudoJet smearer_j_CMS_Benj(double E, double pt, double pz, double phi, double eta){
       double New_E=0;
       double New_pt=0;
       if(abs(eta)<=3.0){
          New_E=sqrt( pow(E*0.05, 2)  + E*pow(1.5, 2) );
       }
       if( 3.0<abs(eta) && abs(eta)<=5.0){
          New_E=sqrt( pow(E*0.130, 2)  + E*pow(2.7, 2) );
       }
       if(abs(eta)<=0.5 &&  0.1 <pt){
         New_pt=sqrt(pow(0.06,2) + pow(pt*0.0013, 2));
       }
       if(0.5 < abs(eta) and abs(eta) <= 1.5 and 0.1 < pt){
          New_pt=sqrt(pow(0.10,2) + pow(pt*0.0017, 2));
       }
       if(1.5 < abs(eta) && abs(eta)<=2.5 &&  0.1<pt){
         New_pt=sqrt(pow(0.25,2) + pow(pt*0.0031, 2));
       }

       // construct a trivial random generator engine from a time-based seed:
       unsigned seed_E = std::chrono::system_clock::now().time_since_epoch().count();
       std::default_random_engine generator_E (seed_E);
       std::normal_distribution<double> distribution_E (0.0, New_E);

       unsigned seed_pt = std::chrono::system_clock::now().time_since_epoch().count();
       std::default_random_engine generator_pt (seed_pt);
       std::normal_distribution<double> distribution_pt (0.0, New_pt);

       //Use Benj idea based on Root function
       TRandom *distribution_E_Benj  = new TRandom();
       double shift_E=distribution_E_Benj->Gaus(0.0, New_E);

       //Use Benj idea based on Root function
       TRandom *distribution_pt_Benj  = new TRandom();
       double shift_pt=distribution_pt_Benj->Gaus(0.0, New_pt);


       TLorentzVector Momentum;
       Momentum.SetPtEtaPhiE( (E +shift_E)/cosh(eta), eta, phi,  E +shift_E);

//       fastjet::PseudoJet jet_final(pz, (pt + distribution_pt(generator_pt))*sin(phi), (pt + distribution_pt(generator_pt))*cos(phi), E + distribution_E(generator_E));

//       fastjet::PseudoJet jet_final((pt + shift_pt)*cos(phi), (pt + shift_pt)*sin(phi), pz, E + shift_E);

       fastjet::PseudoJet jet_final(Momentum.Px(), Momentum.Py(), Momentum.Pz(), Momentum.E());

       return jet_final;
}


fastjet::PseudoJet smearer_j_ATLAS_Benj(double E, double pt, double pz, double phi, double eta){
       double New_E=0;
       double New_pt=0;
//First modify the energy
       if(abs(eta)<=1.7){
         New_E=sqrt(pow(E,2)*pow(0.0302,2) + E*pow(0.5205,2)+pow(1.59,2));
       }
       if(1.7<abs(eta) && abs(eta)<=3.2){
         New_E=sqrt(pow(E,2)*pow(0.05,2) + E*pow(0.706,2));
       }
       if(3.2<abs(eta) && abs(eta)<=4.9){
         New_E=sqrt(pow(E,2)*pow(0.0942,2) + E);
       }
//Second modify the pt
       if(abs(eta)<=0.5 and 0.1<pt){
          New_pt=sqrt(pow(0.06, 2) + pow(pt, 2)*pow(0.0013,2));
       }
       if(abs(eta)<=1.5 and 0.5 < abs(eta) and 0.1<pt){
         New_pt=sqrt(pow(0.10, 2) + pow(pt, 2)*pow(0.0017,2));
       }
       if(abs(eta)<=2.5 and 1.5 < abs(eta) and 0.1<pt){
         New_pt=sqrt(pow(0.25, 2) + pow(pt, 2)*pow(0.0031,2));
       }

       // construct a trivial random generator engine from a time-based seed:
       unsigned seed_E = std::chrono::system_clock::now().time_since_epoch().count();
       std::default_random_engine generator_E (seed_E);
       std::normal_distribution<double> distribution_E (0.0, New_E);

       //Use Benj idea based on Root function
       TRandom *distribution_E_Benj  = new TRandom();
       double shift_E=distribution_E_Benj->Gaus(0.0, New_E);


       unsigned seed_pt = std::chrono::system_clock::now().time_since_epoch().count();
       std::default_random_engine generator_pt (seed_pt);
       std::normal_distribution<double> distribution_pt (0.0, New_pt);

       //Use Benj idea based on Root function
       TRandom *distribution_pt_Benj  = new TRandom();
       double shift_pt=distribution_pt_Benj->Gaus(0.0, New_pt);

       TLorentzVector Momentum;
       Momentum.SetPtEtaPhiE( (E +shift_E)/cosh(eta), eta, phi,  E +shift_E);

//       fastjet::PseudoJet jet_final(pz, (pt + distribution_pt(generator_pt))*sin(phi), (pt + distribution_pt(generator_pt))*cos(phi), E + distribution_E(generator_E));

//       fastjet::PseudoJet jet_final((pt + shift_pt)*cos(phi), (pt + shift_pt)*sin(phi), pz, E + shift_E);

       fastjet::PseudoJet jet_final(Momentum.Px(), Momentum.Py(), Momentum.Pz(), Momentum.E());


       return jet_final;
}

void sort_by_btag_desc(std::vector<fastjet::PseudoJet>& jets) {
  std::sort(jets.begin(), jets.end(),
            [](fastjet::PseudoJet const& a, fastjet::PseudoJet const& b) {
              return a.user_info<MyInfo>().btagscore() > b.user_info<MyInfo>().btagscore();
            });
}
