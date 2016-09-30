#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TString.h"

#include "StPicoDstMaker/StPicoEvent.h"

#include "StPicoBesNetParticleHists.h"
#include "StPicoBesNetParticleCuts.h"

#define TRACK_THN 0
#define EVENT_THN 0

class StPicoPrescales;
ClassImp(StPicoBesNetParticleHists)

// _________________________________________________________
StPicoBesNetParticleHists::StPicoBesNetParticleHists() : TNamed("StPicoBesNetParticleHists", "StPicoBesNetParticleHists"),
  mAnalysisIdx(0),
  mNCentralityBins(8),
  mEventList(NULL), mQAList(NULL), mQARunList(NULL),
  mHnEvent(NULL), mHnTrack(NULL),
//  mPrescales(NULL), 
  mNRuns(0) {
}

// _________________________________________________________
StPicoBesNetParticleHists::StPicoBesNetParticleHists(const Char_t* name) : TNamed(name, name),
  mAnalysisIdx(0),
  mNCentralityBins(8),
  mEventList(NULL), mQAList(NULL), mQARunList(NULL), 
  mHnEvent(NULL), mHnTrack(NULL),
  //  mPrescales(NULL), 
  mNRuns(0) {
}

/*
 * ---------------------------------------------------------------------------------
 *                    Static Const Members - public
 * ---------------------------------------------------------------------------------
 */

// -- Energy / Analysis / RefMult names
const Int_t   StPicoBesNetParticleHists::kNEnergies          = 8;
const Char_t* StPicoBesNetParticleHists::kNnergies[]         = {  "7",   "11",   "14",   "19",   "27",   "39",   "62", "200"};
const Char_t* StPicoBesNetParticleHists::kExactEnergies[]    = {"7.7", "11.5", "14.5", "19.6", "27.0", "39.0", "62.4", "200"};

const Int_t   StPicoBesNetParticleHists::kNNames             = kMaxAnalysisType;
const Char_t* StPicoBesNetParticleHists::kName[]             = {"NetCharge", "NetProton", "NetKaon"};
const Char_t* StPicoBesNetParticleHists::kNameShort[]        = {"charge", "proton", "kaon"};

const Int_t   StPicoBesNetParticleHists::kNRefMult           = kMaxAnalysisType;
const Char_t* StPicoBesNetParticleHists::kNameRefMult[]      = {"RefMult2", "RefMult3", "RefMult4"};
const Char_t* StPicoBesNetParticleHists::kNameRefMultShort[] = {"refmult2", "refmult3", "refmult4"};

const Int_t   StPicoBesNetParticleHists::kNParticles         = 2;
const Char_t* StPicoBesNetParticleHists::kParticleName[kMaxAnalysisType][2]   = {{"neg", "pos"}, {"pbar", "p" }, {"Km",  "Kp" }};
const Char_t* StPicoBesNetParticleHists::kParticleTitle[kMaxAnalysisType][2]  = {{"neg.", "pos."}, {"p#bar", "p"},{"K^-", "K^+"}};
// XX fix for P and K

// -- Histogram sets : nSets / names / titles
const Int_t   StPicoBesNetParticleHists::kNQASets       = 2;
const Char_t* StPicoBesNetParticleHists::kQANames[]     = { "before", "after"};
const Char_t* StPicoBesNetParticleHists::kQATitles[]    = { " (before cuts)", " (after cuts)"};
const Int_t   StPicoBesNetParticleHists::kNQARunSets    = 4;
const Char_t* StPicoBesNetParticleHists::kQARunNames[]  = { "after_good", "after_bad", 
							    "before_good", "before_bad"};
const Char_t* StPicoBesNetParticleHists::kQARunTitles[] = { " (after cuts - good runs)", " (after cuts - bad runs)", 
							    " (before cuts - good runs)", " (before cuts - bad runs)"};

const Int_t   StPicoBesNetParticleHists::kNMultSets     = 5;
const Char_t* StPicoBesNetParticleHists::kMultNames[]   = {"_base", 
							   "_before_TOF_rejected",  "_before_TOF", 
							   "_after_TOF_rejected", "_after_TOF"};
const Char_t* StPicoBesNetParticleHists::kMultTitles[]  = {" (before)",
							   " (before TOF cuts rejected)", " (before TOF cuts)", 
							   " (after TOF cuts rejected)", " (after TOF cuts)"};

// -- Histogram/THnSparse binnings
const Int_t    StPicoBesNetParticleHists::kNHnEvent     = 11;
const Int_t    StPicoBesNetParticleHists::kBinHnEvent[] = {  10,   41,   41,    401,   21,   601,   601,   3001,    2,    201,   21};
const Double_t StPicoBesNetParticleHists::kMinHnEvent[] = {-0.5, -2.0, -3.0, -100.0,  0.0,   0.0,   0.0,    0.0, -0.5, -100.0,  0.0};
const Double_t StPicoBesNetParticleHists::kMaxHnEvent[] = { 9.5,  2.0,  1.0,  100.0,  2.0, 600.0, 600.0, 3000.0,  1.5,  100.0, 20.0};

const Int_t    StPicoBesNetParticleHists::kNHnTrack     = 12;
const Int_t    StPicoBesNetParticleHists::kBinHnTrack[] = {  10,  34,   21,    3,  50,   51,   51,   51, 101,    2,    2, 81};
const Double_t StPicoBesNetParticleHists::kMinHnTrack[] = {-0.5, 0.1, -1.0, -1.5,   0,  0.0,  0.0,  0.0, 0.0, -0.5, -0.5, -4};
const Double_t StPicoBesNetParticleHists::kMaxHnTrack[] = { 9.5, 3.0,  1.0,  1.5,   5, 50.0, 50.0, 50.0, 1.0,  1.5,  1.5,  4};


// -- Centrality Names

const Int_t   StPicoBesNetParticleHists::kNCentralityBins = 11;
//                                                                    0       1         2         3         4         5         6         7         8         9         10
const Char_t* StPicoBesNetParticleHists::kCentralityTitles[]    = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};
const Char_t* StPicoBesNetParticleHists::kCentralityMaxNames[] = {   "5",    "10",     "20",     "30",     "40",     "50",    "60",      "70",     "80",     "90",     "100"}; 

/*
 * ---------------------------------------------------------------------------------
 *
 * ---------------------------------------------------------------------------------
 */

// _________________________________________________________
StPicoBesNetParticleHists::~StPicoBesNetParticleHists() {
  // -- Destructor

  //    Note that histograms are owned by mOutFile. They will be destructed 
  //    when the file is closed.

}

// _________________________________________________________
void StPicoBesNetParticleHists::init(TList* outList, const Char_t* title,
				     Int_t analysisIdx, Int_t nCentralityBins) {
  // -- init method to set up internal lists /hists

  // path to lists of triggers prescales
  // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
  // const char * prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";
  // mPrescales = new StPicoPrescales(prescalesFilesDirectoryName); // fix dir name
  // mNRuns = mPrescales->numberOfRuns();
   
  mAnalysisIdx = analysisIdx;
  mNCentralityBins = nCentralityBins; 

  // -----------------------------------------------------------------------
  // -- Define Lists
  // -----------------------------------------------------------------------

  // -- event list
  outList->Add(new TList);
  mEventList = static_cast<TList*>(outList->Last());
  mEventList->SetOwner(kTRUE);
  mEventList->SetName("BesEventHists");

  // -- event list
  outList->Add(new TList);
  mAnaList = static_cast<TList*>(outList->Last());
  mAnaList->SetOwner(kTRUE);
  mAnaList->SetName("BesAnaHists");

  // -- QA list
  outList->Add(new TList);
  mQAList = static_cast<TList*>(outList->Last());
  mQAList->SetOwner(kTRUE);
  mQAList->SetName("BesQAHists");

  // // -- Run QA list
  // outList->Add(new TList);
  // mQARunList = static_cast<TList*>(outList->Last());
  // mQARunList->SetOwner(kTRUE);
  // mQARunList->SetName("BESQARunHists");

  // -----------------------------------------------------------------------
  // -- Define Histograms
  // -----------------------------------------------------------------------
  InitializeEventStats();
  InitializeMultiplicityStats();
  InitializeEventQAHists();
  InitializeTrackQAHists();

  // InitializeRunByRunEventQAHists();
  // InitializeRunByRunTrackQAHists();

  AddHistSetCent("Dist",       title);
  AddHistSetCent("Dist_lower", title);
  AddHistSetCent("Dist_upper", title);

  // ------------------------------------------------------------------
  // -- Get event container
  // ------------------------------------------------------------------
#if EVENT_THN
  mEventList->Add(new THnSparseD("hnEvent", "cent:vx:vy:vz:shiftedVtx:nRefMultX:nRefMultXCorr:nTracks:isRejected:vzVPD:deltaVz", 
				 StPicoBesNetParticleHists::kNHnEvent, 
				 StPicoBesNetParticleHists::kBinHnEvent, 
				 StPicoBesNetParticleHists::kMinHnEvent, StPicoBesNetParticleHists::kMaxHnEvent));  

  mHnEvent = static_cast<THnSparseD*>(mEventList->Last());
  mHnEvent->Sumw2(); 
  mHnEvent->GetAxis(0)->SetTitle("centrality");
  mHnEvent->GetAxis(1)->SetTitle("#it{v}_{x} (cm)");
  mHnEvent->GetAxis(2)->SetTitle("#it{v}_{y} (cm)");
  mHnEvent->GetAxis(3)->SetTitle("#it{v}_{z} (cm)");
  mHnEvent->GetAxis(4)->SetTitle("shifted vtx (cm)");
  mHnEvent->GetAxis(5)->SetTitle("refMultX");
  mHnEvent->GetAxis(6)->SetTitle("refMultXCorr");
  mHnEvent->GetAxis(7)->SetTitle("nTracks");
  mHnEvent->GetAxis(8)->SetTitle("isRejected");
  mHnEvent->GetAxis(9)->SetTitle("#it{v}_{z}^{vpd} (cm)");
  mHnEvent->GetAxis(10)->SetTitle("#Delta#it{v} (cm)");
#endif

  // ------------------------------------------------------------------
  // -- Get tracks container
  // ------------------------------------------------------------------
#if TRACK_THN
  mEventList->Add(new THnSparseD("hnTrack", "cent:pt:eta:sign:dca:nHitsDedx:nHitsFit:nFitPoss:ratio:isInRefMult:isTrackAccepted;nSigmaProton",
				 StPicoBesNetParticleHists::kNHnTrack, 
				 StPicoBesNetParticleHists::kBinHnTrack, 
				 StPicoBesNetParticleHists::kMinHnTrack, StPicoBesNetParticleHists::kMaxHnTrack));  

  mHnTrack = static_cast<THnSparseD*>(mEventList->Last());
  mHnTrack->Sumw2(); 
  mHnTrack->GetAxis(0)->SetTitle("centrality");
  mHnTrack->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  mHnTrack->GetAxis(2)->SetTitle("#eta");
  mHnTrack->GetAxis(3)->SetTitle("sign");
  mHnTrack->GetAxis(4)->SetTitle("dca");
  mHnTrack->GetAxis(5)->SetTitle("nHitsDedx");
  mHnTrack->GetAxis(6)->SetTitle("nHitsFit");
  mHnTrack->GetAxis(7)->SetTitle("nFitPoss");
  mHnTrack->GetAxis(8)->SetTitle("nHitsFit/nFitPoss");
  mHnTrack->GetAxis(9)->SetTitle("isInRefMult");
  mHnTrack->GetAxis(10)->SetTitle("isTrackAccepted");
  mHnTrack->GetAxis(11)->SetTitle("nSigmaProton");
#endif

}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillHnEvent(Double_t *aEvent) {
  // -- Fill Track THnSparse
  
#if EVENT_THN
  mHnEvent->Fill(aEvent);
#endif
}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillHnTrack(Double_t *aTrack) {
  // -- Fill Track THnSparse
  
#if TRACK_THN
  mHnTrack->Fill(aTrack);
#endif
}

//________________________________________________________________________
void StPicoBesNetParticleHists::InitializeEventStats() {
  // -- Initialize event statistics histograms
  
  mEventList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", 
			   StPicoBesNetParticleCuts::kNEventStat, -0.5, StPicoBesNetParticleCuts::kNEventStat-0.5));
  TH1F *hEventStat0 = static_cast<TH1F*>(mEventList->Last());

  mEventList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", 
			   StPicoBesNetParticleCuts::kNEventStat, -0.5, StPicoBesNetParticleCuts::kNEventStat-0.5));
  TH1F *hEventStat1 = static_cast<TH1F*>(mEventList->Last());

  for (Int_t ii=0; ii < StPicoBesNetParticleCuts::kNEventStat -1; ii++) {
    hEventStat0->GetXaxis()->SetBinLabel(ii+1,  StPicoBesNetParticleCuts::kEventStatNames[ii]);
    hEventStat1->GetXaxis()->SetBinLabel(ii+1,  StPicoBesNetParticleCuts::kEventStatNames[ii]);
  }
}
 
//________________________________________________________________________
Bool_t StPicoBesNetParticleHists::FillEventStats(Int_t *aEventCuts) {
  // -- Fill event / centrality statistics 
  //    Also checks if event is rejected
  // 
  //   aEventCuts is an array with one entry for every cut
  //    0 - event passes cut
  //    1 - event doesn't pass cut

  Bool_t isRejected = kFALSE;
  
  // -- Fill event statistics
  for (Int_t idx = 0; idx < StPicoBesNetParticleCuts::kNEventStat; ++idx) {
    if (aEventCuts[idx])
      isRejected = kTRUE;
    else
      (static_cast<TH1F*>(mEventList->FindObject("hEventStat0")))->Fill(idx);
  }
  
  for (Int_t idx = 0; idx < StPicoBesNetParticleCuts::kNEventStat; ++idx) {
    if (aEventCuts[idx])
      break;
    (static_cast<TH1F*>(mEventList->FindObject("hEventStat1")))->Fill(idx);
  }
  
  return isRejected;
}

//________________________________________________________________________
void StPicoBesNetParticleHists::AddHistSetCent(const Char_t *name, const Char_t *title)  {
   // -- Add histogram sets for particle and anti-particle
   //    dependence : centrality

   TString sName(name);
   TString sTitle(title);

   // -- Add List
   mAnaList->Add(new TList);
   TList *list = static_cast<TList*>(mAnaList->Last());
   list->SetName(Form("f%s", name));
   list->SetOwner(kTRUE);

   // -- Create Titles
   TString sNetTitle(Form("N_{%s} - N_{%s}", 
			  StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][1], 
			  StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][0]));
   TString sSumTitle(Form("N_{%s} + N_{%s}", 
			  StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][1],
			  StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][0]));
   
   // -----------------------------------------------------------------------------------------------
   
   // -- Add Particle / Anti-Particle Distributions
   for (Int_t idxPart = 0; idxPart < StPicoBesNetParticleHists::kNParticles; ++idxPart) {
     list->Add(new TH2D(Form("h%s%s", name, StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][idxPart]), 
			Form("N_{%s} : %s;Centrality;N_{%s}", 
			     StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][idxPart], sTitle.Data(), 
			     StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][idxPart]),
			mNCentralityBins, -0.5, mNCentralityBins-0.5, 601, -0.5, 600.49));
   } // for (Int_t idxPart = 0; idxPart < 2; ++idxPart) {
 
   // -- Add Particle vs Anti-Particle Distribution
   list->Add(new TH2D(Form("h%s%s%s", name, 
			   StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][0], 
			   StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][1]), 
		      Form("N_{%s} vs N_{%s} : %s;N_{%s};N_{%s}", 
			   StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][0], 
			   StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][1], sTitle.Data(), 
			   StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][0], 
			   StPicoBesNetParticleHists::kParticleTitle[mAnalysisIdx][1]),
		      601, -0.5, 600.49, 601, -0.5, 600.49)); 

   // -- Add NetParticle Distributions
   list->Add(new TH2D(Form("h%s%s", name, StPicoBesNetParticleHists::kName[mAnalysisIdx]), 
		      Form("%s : %s;Centrality;%s", sNetTitle.Data(), sTitle.Data(), sNetTitle.Data()), 
		      mNCentralityBins, -0.5, mNCentralityBins-0.5, 301, -300.5, 300.49));

   // -- Add NetParticle vs SumParticle
   list->Add(new TH2D(Form("h%s%sOverSum", name, StPicoBesNetParticleHists::kName[mAnalysisIdx]), 
		      Form("(%s)/(%s) : %s;Centrality;(%s)/(%s)", sNetTitle.Data(), sSumTitle.Data(), 
			   sTitle.Data(), sNetTitle.Data(), sSumTitle.Data()), 
		      mNCentralityBins, -0.5, mNCentralityBins-0.5, 41, -2.5, 2.49));
   return;
}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillHistSetCent(const Char_t *name, Int_t* np, Int_t cent)  {
  // -- Fill histogram sets for particle and anti-particle
  //    dependence : centrality 
  
  // -- Get List
   TList *list = static_cast<TList*>(mAnaList->FindObject(Form("f%s",name)));
  
  // -- Get Centrality Bin
  Float_t centralityBin = Float_t(cent);

  // -----------------------------------------------------------------------------------------------

  Int_t sumNp   = np[1]+np[0];  // particle + anti-particles
  Int_t deltaNp = np[1]-np[0];  // particle - anti-particles

  Double_t deltaNpOverSumNp = (sumNp == 0.) ? 0. : deltaNp/Double_t(sumNp);

  // -- Fill Particle / Anti-Particle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%s%s", name, 
					    StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][0]))))->Fill(centralityBin, np[0]);
  (static_cast<TH2D*>(list->FindObject(Form("h%s%s", name, 
					    StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][1]))))->Fill(centralityBin, np[1]);

  (static_cast<TH2D*>(list->FindObject(Form("h%s%s%s", name, 
					    StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][0], 
					    StPicoBesNetParticleHists::kParticleName[mAnalysisIdx][1]))))->Fill(np[0], np[1]);

  // -- Fill NetParticle Distributions
  (static_cast<TH2D*>(list->FindObject(Form("h%sNet%s", name, 
					    StPicoBesNetParticleHists::kName[mAnalysisIdx]))))->Fill(centralityBin, deltaNp);

  // -- Fill NetParticle vs SumParticle
  (static_cast<TH2D*>(list->FindObject(Form("h%s%sOverSum", name, 
					    StPicoBesNetParticleHists::kName[mAnalysisIdx]))))->Fill(centralityBin, deltaNpOverSumNp);

  // -----------------------------------------------------------------------------------------------

  return;
}

//________________________________________________________________________
void StPicoBesNetParticleHists::InitializeMultiplicityStats() {
  // -- Initialize multiplicity statistics histograms

  for (Int_t ii = 0; ii < StPicoBesNetParticleHists::kNMultSets; ++ii) {

    mEventList->Add(new TList);
    TList *list = static_cast<TList*>(mEventList->Last());
    list->SetName(Form("f%s_multHists%s", StPicoBesNetParticleHists::kName[mAnalysisIdx], 
		       StPicoBesNetParticleHists::kMultNames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH1F(Form("hCentralityStat%s", StPicoBesNetParticleHists::kMultNames[ii]),               
		       Form("Centrality statistics%s;Centrality Bins;Events",
			    StPicoBesNetParticleHists::kMultTitles[ii]),
		       mNCentralityBins, -0.5, mNCentralityBins-0.5));      

    TH1F* hCentralityStat = static_cast<TH1F*>(list->Last());    
    for (Int_t jj = 0; jj < mNCentralityBins; jj++)
      hCentralityStat->GetXaxis()->SetBinLabel(jj+1, StPicoBesNetParticleHists::kCentralityTitles[jj]);
  
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    list->Add(new TH1F(Form("hRefMultStat%s", StPicoBesNetParticleHists::kMultNames[ii]),                    
		       Form("RefMult Statistics%s;RefMult;Events", StPicoBesNetParticleHists::kMultTitles[ii]),
		       601, 0., 600.));

    list->Add(new TH1F(Form("h%sStat%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultNames[ii]), 
		       Form("%s Statistics%s;%s;Events", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultTitles[ii], 
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),  
		       601, 0., 600.));

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH2F(Form("h%s_nGlobalTracks%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultNames[ii]),         
		       Form("%s vs nGlobalTracks%s;%s;nGlobalTracks", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultTitles[ii], 
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
		       601, 0., 600., 2501, 0., 2500.));
    
    list->Add(new TH2F(Form("h%s_nPrimaryTracks%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultNames[ii]),        
		       Form("%s vs nPrimaryTracks%s;%s;nPrimaryTracks", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultTitles[ii], 
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
		       601, 0., 600., 1001, 0., 1000.));
    
    list->Add(new TH2F(Form("h%s_nTOFMatch%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultNames[ii]),   
		       Form("%s vs nTOFMatch%s;%s;nTOFMatch", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultTitles[ii], 
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
		       601, 0., 600., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hRefMult_nGlobalTracks%s", StPicoBesNetParticleHists::kMultNames[ii]),          
		       Form("RefMult vs nGlobalTracks%s;RefMult;nGlobalTracks", StPicoBesNetParticleHists::kMultTitles[ii]),
		       601, 0., 600., 2501, 0., 2500.));

    list->Add(new TH2F(Form("hRefMult_nPrimaryTracks%s", StPicoBesNetParticleHists::kMultNames[ii]),         
		       Form("RefMult vs nPrimaryTracks%s;RefMult;nPrimaryTracks", StPicoBesNetParticleHists::kMultTitles[ii]),
		       601, 0., 600., 1001, 0., 1000.));

    list->Add(new TH2F(Form("hRefMult_nTOFMatch%s", StPicoBesNetParticleHists::kMultNames[ii]),              
		       Form("RefMult vs nTOFMatch%s;RefMult;nTOFMatch", StPicoBesNetParticleHists::kMultTitles[ii]),
		       601, 0., 600., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("h%s_%sCorr%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], StPicoBesNetParticleHists::kMultNames[ii]),          
		       Form("%s vs %sCorr%s;%s;%sCorr", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], StPicoBesNetParticleHists::kMultTitles[ii],
			    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
		       601, 0., 600., 601, 0., 600.));
    
    list->Add(new TH2F(Form("h%s_RefMult%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultNames[ii]),               
		       Form("%s vs RefMult%s;%s;RefMult", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
			    StPicoBesNetParticleHists::kMultTitles[ii], StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
		       601, 0., 600., 601, 0., 600.));
    
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    list->Add(new TH2F(Form("hnPrimaryTracks_nGlobalTracks%s", StPicoBesNetParticleHists::kMultNames[ii]),   
		       Form("nPrimaryTracks vs nGlobalTracks%s;nPrimaryTracks;nGlobalTracks", StPicoBesNetParticleHists::kMultTitles[ii]),
		       1001, 0., 1000., 2501, 0., 2500.));

    list->Add(new TH2F(Form("hnPrimaryTracks_nTOFMatch%s", StPicoBesNetParticleHists::kMultNames[ii]),       
		       Form("nPrimaryTracks vs nTOFMatch%s; nPrimaryTracks;nTOFMatch", StPicoBesNetParticleHists::kMultTitles[ii]),
		       1001, 0., 1000., 601, 0., 600.));    
    
    list->Add(new TH2F(Form("hnGlobalTracks_nTOFMatch%s", StPicoBesNetParticleHists::kMultNames[ii]),        
		       Form("nGlobalTracks vs nTOFMatch%s;nGlobalTracks;nTOFMatch", StPicoBesNetParticleHists::kMultTitles[ii]),
		       2501, 0., 2500., 601, 0., 600.));
      
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Add this only for the last set 

    if (ii == StPicoBesNetParticleHists::kNMultSets-1) {
      list->Add(new TProfile(Form("h%sMean", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
			     Form("%s Mean vs Centrality;Centrality Bins;<%s>", 
				  StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
				  StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
			     mNCentralityBins, -0.5, mNCentralityBins-0.5));      
      TH1F* hRefMultXMean = static_cast<TH1F*>(list->Last());
      
      list->Add(new TProfile(Form("h%sCorrMean", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),        
			     Form("%sCorr Mean vs Centrality;Centrality Bins;<%sCorr>", 
				  StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
				  StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]),
			     mNCentralityBins, -0.5, mNCentralityBins-0.5));       
      TH1F* hRefMultXCorrMean = static_cast<TH1F*>(list->Last());
      
      // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      
      for (Int_t jj = 0; jj < mNCentralityBins; jj++) {
	hRefMultXMean->GetXaxis()->SetBinLabel(jj+1, StPicoBesNetParticleHists::kCentralityTitles[jj]);
	hRefMultXCorrMean->GetXaxis()->SetBinLabel(jj+1, StPicoBesNetParticleHists::kCentralityTitles[jj]);
      }
    } // if (ii == StPicoBesNetParticleHists::kNMultSets-1) {
  }    
}
 
//________________________________________________________________________
void StPicoBesNetParticleHists::FillMultiplicityStats(Double_t *aMult, Int_t mode) {
  // -- Fill multiplicity statistics of accepted events

  TList* list = static_cast<TList*>(mEventList->FindObject(Form("f%s_multHists%s", 
								StPicoBesNetParticleHists::kName[mAnalysisIdx], 
								StPicoBesNetParticleHists::kMultNames[mode])));

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  (static_cast<TH1F*>(list->FindObject(Form("hCentralityStat%s",         
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[0]);

  (static_cast<TH1F*>(list->FindObject(Form("hRefMultStat%s",            
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[1]);  
  (static_cast<TH1F*>(list->FindObject(Form("h%sStat%s",                 
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[2]);

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  (static_cast<TH2F*>(list->FindObject(Form("h%s_nGlobalTracks%s",       
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[2], aMult[4]);
  (static_cast<TH2F*>(list->FindObject(Form("h%s_nPrimaryTracks%s",      
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[2], aMult[5]);
  (static_cast<TH2F*>(list->FindObject(Form("h%s_nTOFMatch%s",         
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[2], aMult[6]);
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult_nGlobalTracks%s",  
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[1], aMult[4]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult_nPrimaryTracks%s", 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[1], aMult[5]);
  (static_cast<TH2F*>(list->FindObject(Form("hRefMult_nTOFMatch%s",      
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[1], aMult[6]);
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  
  (static_cast<TH2F*>(list->FindObject(Form("h%s_%sCorr%s",  
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx],
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[2], aMult[3]);
  (static_cast<TH2F*>(list->FindObject(Form("h%s_RefMult%s", 
					    StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[2], aMult[1]);
  (static_cast<TH2F*>(list->FindObject(Form("hnPrimaryTracks_nGlobalTracks%s", 
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[5], aMult[4]);
  
  (static_cast<TH2F*>(list->FindObject(Form("hnPrimaryTracks_nTOFMatch%s",     
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[5], aMult[6]);
  (static_cast<TH2F*>(list->FindObject(Form("hnGlobalTracks_nTOFMatch%s",     
					    StPicoBesNetParticleHists::kMultNames[mode]))))->Fill(aMult[4], aMult[6]);
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  if (mode == StPicoBesNetParticleHists::kNMultSets-1) {
    (static_cast<TProfile*>(list->FindObject(Form("h%sMean", 
						  StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]))))->Fill(aMult[0], aMult[2]);
    (static_cast<TProfile*>(list->FindObject(Form("h%sCorrMean", 
						  StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]))))->Fill(aMult[0], aMult[3]);
  }
}

//________________________________________________________________________
void StPicoBesNetParticleHists::InitializeEventQAHists() {
  // -- Initialize event QA histograms

  for (Int_t ii = 0; ii < StPicoBesNetParticleHists::kNQASets; ++ii) {
    mQAList->Add(new TList);
    TList *list = static_cast<TList*>(mQAList->Last());
    list->SetName(Form("f%s_eventHists_%s", 
		       StPicoBesNetParticleHists::kName[mAnalysisIdx], 
		       StPicoBesNetParticleHists::kQANames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    list->Add(new TH1F(Form("vx_%s", StPicoBesNetParticleHists::kQANames[ii]),        
		       Form("#it{v}_{x}%s;#it{v}_{x} (cm);Events", StPicoBesNetParticleHists::kQATitles[ii]), 
		       StPicoBesNetParticleHists::kBinHnEvent[1], 
		       StPicoBesNetParticleHists::kMinHnEvent[1], StPicoBesNetParticleHists::kMaxHnEvent[1]));
    
    list->Add(new TH1F(Form("vy_%s", StPicoBesNetParticleHists::kQANames[ii]),        
		       Form("#it{v}_{y}%s;#it{v}_{y} (cm);Events", StPicoBesNetParticleHists::kQATitles[ii]),
		       StPicoBesNetParticleHists::kBinHnEvent[2], 
		       StPicoBesNetParticleHists::kMinHnEvent[2], StPicoBesNetParticleHists::kMaxHnEvent[2]));
    
    list->Add(new TH1F(Form("vz_%s", StPicoBesNetParticleHists::kQANames[ii]),        
		       Form("#it{v}_{z}%s;#it{v}_{z} (cm);Events", StPicoBesNetParticleHists::kQATitles[ii]), 
		       StPicoBesNetParticleHists::kBinHnEvent[3], 
		       StPicoBesNetParticleHists::kMinHnEvent[3], StPicoBesNetParticleHists::kMaxHnEvent[3]));
    
    list->Add(new TH1F(Form("shiftedVr_%s", StPicoBesNetParticleHists::kQANames[ii]), 
		       Form("#it{v}_{r}^{shifted}%s;#it{v}_{r}^{shifted} (cm);Events", StPicoBesNetParticleHists::kQATitles[ii]), 
		       StPicoBesNetParticleHists::kBinHnEvent[4], 
		       StPicoBesNetParticleHists::kMinHnEvent[4], StPicoBesNetParticleHists::kMaxHnEvent[4]));
    
    list->Add(new TH2F(Form("vxvy_%s", StPicoBesNetParticleHists::kQANames[ii]),      
		       Form("#it{v}_{x} vs #it{v}_{y}%s;#it{v}_{x} (cm); #it{v}_{y} (cm)", StPicoBesNetParticleHists::kQATitles[ii]), 
		       StPicoBesNetParticleHists::kBinHnEvent[1], 
		       StPicoBesNetParticleHists::kMinHnEvent[1], StPicoBesNetParticleHists::kMaxHnEvent[1], 
		       StPicoBesNetParticleHists::kBinHnEvent[2], 
		       StPicoBesNetParticleHists::kMinHnEvent[2], StPicoBesNetParticleHists::kMaxHnEvent[2]));
    
    list->Add(new TH1F(Form("vzVpd_%s", StPicoBesNetParticleHists::kQANames[ii]),     
		       Form("#it{v}_{z}^{vpd}%s;#it{v}_{z}^{vpd} (cm);Events", StPicoBesNetParticleHists::kQATitles[ii]),
		       StPicoBesNetParticleHists::kBinHnEvent[9], 
		       StPicoBesNetParticleHists::kMinHnEvent[9], StPicoBesNetParticleHists::kMaxHnEvent[9]));
    
    list->Add(new TH1F(Form("deltaVz_%s", StPicoBesNetParticleHists::kQANames[ii]),   
		       Form("#Delta#it{v}_{z}%s;#Delta#it{v}_{z} (cm);Events", StPicoBesNetParticleHists::kQATitles[ii]),
		       StPicoBesNetParticleHists::kBinHnEvent[10], 
		       StPicoBesNetParticleHists::kMinHnEvent[10], StPicoBesNetParticleHists::kMaxHnEvent[10]));
  }
}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillEventQAHists(Double_t *aEvent, Int_t mode) {
  // -- Fill event QA histograms      

  TList* list = static_cast<TList*>(mQAList->FindObject(Form("f%s_eventHists_%s", 
							     StPicoBesNetParticleHists::kName[mAnalysisIdx], 
							     StPicoBesNetParticleHists::kQANames[mode])));

  (static_cast<TH1F*>(list->FindObject(Form("vx_%s",        StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[1]);
  (static_cast<TH1F*>(list->FindObject(Form("vy_%s",        StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[2]);
  (static_cast<TH1F*>(list->FindObject(Form("vz_%s",        StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[3]);
  (static_cast<TH1F*>(list->FindObject(Form("shiftedVr_%s", StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[3]);
  (static_cast<TH2F*>(list->FindObject(Form("vxvy_%s",      StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[1], aEvent[2]);
  (static_cast<TH1F*>(list->FindObject(Form("vzVpd_%s",     StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[9]);
  (static_cast<TH1F*>(list->FindObject(Form("deltaVz_%s",   StPicoBesNetParticleHists::kQANames[mode]))))->Fill(aEvent[10]);
}

//________________________________________________________________________
void StPicoBesNetParticleHists::InitializeTrackQAHists() {
  // -- Initialize track QA histograms

  for (Int_t ii = 0; ii < StPicoBesNetParticleHists::kNQASets; ++ii) {
    mQAList->Add(new TList);
    TList *list = static_cast<TList*>(mQAList->Last());
    list->SetName(Form("f%s_trackHists_%s", 
		       StPicoBesNetParticleHists::kName[mAnalysisIdx], 
		       StPicoBesNetParticleHists::kQANames[ii]));
    list->SetOwner(kTRUE);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    for (Int_t idxSign = 0 ; idxSign < 2; ++idxSign) { 
      list->Add(new TH1F(Form("pt_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),                
			 Form("#it{p}_{T}%s;#it{p}_{T} (GeV/#it{c});Tracks",
			      StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[1], 
			 StPicoBesNetParticleHists::kMinHnTrack[1], StPicoBesNetParticleHists::kMaxHnTrack[1]));
      
      list->Add(new TH1F(Form("eta_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),               
			 Form("#eta%s;#eta;Tracks", 
			      StPicoBesNetParticleHists::kQATitles[ii]), 
			 StPicoBesNetParticleHists::kBinHnTrack[2], 
			 StPicoBesNetParticleHists::kMinHnTrack[2], StPicoBesNetParticleHists::kMaxHnTrack[2]));
      
      list->Add(new TH1F(Form("dca_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),               
			 Form("DCA%s;DCA (cm);Tracks", 
			      StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[4], 
			 StPicoBesNetParticleHists::kMinHnTrack[4], StPicoBesNetParticleHists::kMaxHnTrack[4]));
      
      list->Add(new TH1F(Form("nHitsDedx_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),         
			 Form("nHitsDedx%s;nHitsDedx;Tracks", 
			      StPicoBesNetParticleHists::kQATitles[ii]), 
			 StPicoBesNetParticleHists::kBinHnTrack[5], 
			 StPicoBesNetParticleHists::kMinHnTrack[5], StPicoBesNetParticleHists::kMaxHnTrack[5]));

      list->Add(new TH1F(Form("nHitsFit_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),          
			 Form("nHitsFit%s;nHitsFit;Tracks", 
			      StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[6], 
			 StPicoBesNetParticleHists::kMinHnTrack[6], StPicoBesNetParticleHists::kMaxHnTrack[6]));
      
      list->Add(new TH1F(Form("nHitsFit_nFitPoss_%s_%d",
			      StPicoBesNetParticleHists::kQANames[ii], idxSign), 
			 Form("nHitsFit/nFitPoss%s;nHitsFit/nFitPoss;Tracks", StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[8],
			 StPicoBesNetParticleHists::kMinHnTrack[8], StPicoBesNetParticleHists::kMaxHnTrack[8]));
      
      list->Add(new TH1F(Form("nSigmaP_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),           
			 Form("nSigmaProton%s;n#sigma_{proton};Tracks",
			      StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[11], 
			 StPicoBesNetParticleHists::kMinHnTrack[11], StPicoBesNetParticleHists::kMaxHnTrack[11]));
      
      list->Add(new TH2F(Form("dca_pt_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),            
			 Form("DCA vs #it{p}_{T}%s;DCA (cm);#it{p}_{T} (GeV/#it{c})", StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[4], 
			 StPicoBesNetParticleHists::kMinHnTrack[4], StPicoBesNetParticleHists::kMaxHnTrack[4], 
			 StPicoBesNetParticleHists::kBinHnTrack[1], 
			 StPicoBesNetParticleHists::kMinHnTrack[1], StPicoBesNetParticleHists::kMaxHnTrack[1]));
      
      list->Add(new TH2F(Form("nHitsDedx_pt_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),      
			 Form("nHitsDedx vs #it{p}_{T}%s;nHitsDedx;#it{p}_{T} (GeV/#it{c})", StPicoBesNetParticleHists::kQATitles[ii]), 
			 StPicoBesNetParticleHists::kBinHnTrack[5], 
			 StPicoBesNetParticleHists::kMinHnTrack[5], StPicoBesNetParticleHists::kMaxHnTrack[5], 
			 StPicoBesNetParticleHists::kBinHnTrack[1], 
			 StPicoBesNetParticleHists::kMinHnTrack[1], StPicoBesNetParticleHists::kMaxHnTrack[1]));
      
      list->Add(new TH2F(Form("nHitsFit_pt_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),       
			 Form("nHitsFit vs #it{p}_{T}%s;nHitsFit;#it{p}_{T} (GeV/#it{c})", StPicoBesNetParticleHists::kQATitles[ii]),  
			 StPicoBesNetParticleHists::kBinHnTrack[6], 
			 StPicoBesNetParticleHists::kMinHnTrack[6], StPicoBesNetParticleHists::kMaxHnTrack[6], 
			 StPicoBesNetParticleHists::kBinHnTrack[1], 
			 StPicoBesNetParticleHists::kMinHnTrack[1], StPicoBesNetParticleHists::kMaxHnTrack[1]));

      list->Add(new TH2F(Form("nSigmaP_pt_%s_%d", StPicoBesNetParticleHists::kQANames[ii], idxSign),        
			 Form("nSigmaProton vs #it{p}_{T}%s;n#sigma_{proton};#it{p}_{T} (GeV/#it{c})", StPicoBesNetParticleHists::kQATitles[ii]),
			 StPicoBesNetParticleHists::kBinHnTrack[11], 
			 StPicoBesNetParticleHists::kMinHnTrack[11], StPicoBesNetParticleHists::kMaxHnTrack[11], 
			 StPicoBesNetParticleHists::kBinHnTrack[1],  
			 StPicoBesNetParticleHists::kMinHnTrack[1],  StPicoBesNetParticleHists::kMaxHnTrack[1]));
    }
  }
}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillTrackQAHists(Double_t *aTrack, Int_t mode) {
  // -- Fill track QA histograms

  Int_t idxSign = (aTrack[3] < 0) ? 0 : 1;

  TList* list = static_cast<TList*>(mQAList->FindObject(Form("f%s_trackHists_%s", 
							     StPicoBesNetParticleHists::kName[mAnalysisIdx], 
							     StPicoBesNetParticleHists::kQANames[mode])));
  (static_cast<TH1F*>(list->FindObject(Form("pt_%s_%d",                
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[1]);
  (static_cast<TH1F*>(list->FindObject(Form("pcm_%s_%d",               
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[12]);
  (static_cast<TH1F*>(list->FindObject(Form("eta_%s_%d",               
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[2]);
  (static_cast<TH1F*>(list->FindObject(Form("dca_%s_%d",               
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[4]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsDedx_%s_%d",         
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[5]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsFit_%s_%d",         
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[6]);
  (static_cast<TH1F*>(list->FindObject(Form("nHitsFit_nFitPoss_%s_%d",
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[8]);
  (static_cast<TH1F*>(list->FindObject(Form("nSigmaP_%s_%d",        
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[11]);
  
  (static_cast<TH2F*>(list->FindObject(Form("dca_pt_%s_%d",       
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[4], aTrack[1]);
  (static_cast<TH2F*>(list->FindObject(Form("nHitsDedx_pt_%s_%d",
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[5], aTrack[1]);
  (static_cast<TH2F*>(list->FindObject(Form("nHitsFit_pt_%s_%d", 
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[6], aTrack[1]);
  (static_cast<TH2F*>(list->FindObject(Form("nSigmaP_pt_%s_%d",  
					    StPicoBesNetParticleHists::kQANames[mode], idxSign))))->Fill(aTrack[11], aTrack[1]);
}

#if 0 
//________________________________________________________________________
void StPicoBesNetParticleHists::InitializeRunByRunEventQAHists() {
  // -- Initialize run-by-run event distributions
  

  // << nRunsAll


  for (Int_t ii = 0; ii < StPicoBesNetParticleHists::kNQARunSets; ++ii) {
    mQARunList->Add(new TList);
    TList *list = static_cast<TList*>(mQARunList->Last());
    list->SetName(Form("f%s_runByRunEventHists_%s", 
		       StPicoBesNetParticleHists::kName[mAnalysisIdx], 
		       StPicoBesNetParticleHists::kQARunNames[ii]));
    list->SetOwner(kTRUE);

    // -- Over all runs
    list->Add(new TProfile(Form("pAllRefMult_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<RefMult>%s;All;<RefMult>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   1, -0.5, 0.5));

    list->Add(new TProfile(Form("pAll%s_%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<%s>%s;All;<%s>", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunTitles[ii], StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]), 
			   1, -0.5, 0.5));
    
    list->Add(new TProfile(Form("pAll%s_%s", StPicoBesNetParticleHists::kName[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<%s>%s;All;<%s>", StPicoBesNetParticleHists::kName[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunTitles[ii], StPicoBesNetParticleHists::kName[mAnalysisIdx]), 
			   1, -0.5, 0.5));
    
    list->Add(new TProfile(Form("pAllneg_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<neg>%s;All;<neg>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   1, -0.5, 0.5));
    
    list->Add(new TProfile(Form("pAllpos_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<pos>%s;All;<pos>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   1, -0.5, 0.5));

    // -- Run by run
    list->Add(new TProfile(Form("pRefMult_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<RefMult>%s;Run;<RefMult>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));

    list->Add(new TProfile(Form("p%s_%s", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx],
				StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<%s>%s;Run;<%s>", StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunTitles[ii], StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx]), 
			   nRunsAll, -0.5, nRunsAll-0.5));
    
    list->Add(new TProfile(Form("p%s_%s", StPicoBesNetParticleHists::kName[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<%s>%s;Run;<%s>", StPicoBesNetParticleHists::kName[mAnalysisIdx], 
				StPicoBesNetParticleHists::kQARunTitles[ii], StPicoBesNetParticleHists::kName[mAnalysisIdx]), 
			   nRunsAll, -0.5, nRunsAll-0.5));
    
    list->Add(new TProfile(Form("pneg_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<neg>%s;Run;<neg>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));

    list->Add(new TProfile(Form("ppos_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<pos>%s;Run;<pos>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));
  }
}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillRunByRunEventQAHists(Double_t *aEvent, Int_t mode, Int_t runIdx) {
  // -- Fill run-by-run event distributions
  
  TList* list = static_cast<TList*>(mQARunList->FindObject(Form("f%s_runByRunEventHists_%s", 
								StPicoBesNetParticleHists::kName[mAnalysisIdx], 
								StPicoBesNetParticleHists::kQARunNames[mode])));
  
  (static_cast<TProfile*>(list->FindObject(Form("pAllRefMult_%s", 
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(0., aEvent[0]);
  (static_cast<TProfile*>(list->FindObject(Form("pAll%s_%s",      
						StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(0., aEvent[1]);
  (static_cast<TProfile*>(list->FindObject(Form("pAll%s_%s",      
						StPicoBesNetParticleHists::kName[mAnalysisIdx], 
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(0., aEvent[2]);
  (static_cast<TProfile*>(list->FindObject(Form("pAllneg_%s",     
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(0., aEvent[3]);
  (static_cast<TProfile*>(list->FindObject(Form("pAllpos_%s",     
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(0., aEvent[4]);
  

  (static_cast<TProfile*>(list->FindObject(Form("pRefMult_%s",   
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(runIdx, aEvent[0]);
  (static_cast<TProfile*>(list->FindObject(Form("p%s_%s",       
						StPicoBesNetParticleHists::kNameRefMult[mAnalysisIdx], 
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(runIdx, aEvent[1]);
  (static_cast<TProfile*>(list->FindObject(Form("p%s_%s",       
						StPicoBesNetParticleHists::kName[mAnalysisIdx], 
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(runIdx, aEvent[2]);
  (static_cast<TProfile*>(list->FindObject(Form("pneg_%s",      
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(runIdx, aEvent[3]);
  (static_cast<TProfile*>(list->FindObject(Form("ppos_%s",      
						StPicoBesNetParticleHists::kQARunNames[mode]))))->Fill(runIdx, aEvent[4]);
}

//________________________________________________________________________
void StPicoBesNetParticleHists::InitializeRunByRunTrackQAHists() {
  // -- Initialize run-by-run event distributions
  
  for (Int_t ii = 0; ii < StPicoBesNetParticleHists::kNQARunSets; ++ii) {

    mQARunList->Add(new TList);
    TList *list = static_cast<TList*>(mQARunList->Last());
    list->SetName(Form("f%s_runByRunTrackHists_%s", 
		       StPicoBesNetParticleHists::kName[mAnalysisIdx], 
		       StPicoBesNetParticleHists::kQARunNames[ii]));
    list->SetOwner(kTRUE);
    
    list->Add(new TProfile(Form("pPt_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<p_T>%s;Run;<p_T>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));

    list->Add(new TProfile(Form("pPt_neg_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<p_T^{neg}>%s;Run;<p_T^{neg}>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));
    
    list->Add(new TProfile(Form("pPt_pos_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<p_T^{pos}>%s;Run;<p_T^{pos}>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));

    list->Add(new TProfile(Form("pDca_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<dca>%s;Run;<dca>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));
    
    list->Add(new TProfile(Form("pNSigmaProton_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
			   Form("<nSigmaProton>%s;Run;<nSigmaProton>", StPicoBesNetParticleHists::kQARunTitles[ii]), 
			   nRunsAll, -0.5, nRunsAll-0.5));

    
    list->Add(new TH2D(Form("hPt_neg_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
		       Form("p_T^{neg}%s;Run;p_T^{neg}", StPicoBesNetParticleHists::kQARunTitles[ii]), 
		       nRunsAll, -0.5, nRunsAll-0.5,
		       StPicoBesNetParticleHists::kBinHnTrack[1], 
		       StPicoBesNetParticleHists::kMinHnTrack[1], StPicoBesNetParticleHists::kMaxHnTrack[1]));
    
    list->Add(new TH2D(Form("hPt_pos_%s", StPicoBesNetParticleHists::kQARunNames[ii]), 
		       Form("p_T^{pos}%s;Run;p_T^{pos}", StPicoBesNetParticleHists::kQARunTitles[ii]), 
		       nRunsAll, -0.5, nRunsAll-0.5,
		       StPicoBesNetParticleHists::kBinHnTrack[1],
		       StPicoBesNetParticleHists::kMinHnTrack[1], StPicoBesNetParticleHists::kMaxHnTrack[1]));
  }
}

//________________________________________________________________________
void StPicoBesNetParticleHists::FillRunByRunTrackQAHists(Double_t *aTrack, Int_t isBadRun, Int_t mode, Int_t runIdx) {
  // -- Fill run-by-run track distributions

  Int_t idx = mode*2 + isBadRun;

  Int_t idxSign = (aTrack[3] < 0) ? 0 : 1;  

  TList* list = static_cast<TList*>(mQARunList->FindObject(Form("f%s_runByRunTrackHists_%s", 
								StPicoBesNetParticleHists::kName[mAnalysisIdx], 
								StPicoBesNetParticleHists::kQARunNames[idx])));

  (static_cast<TProfile*>(list->FindObject(Form("pPt_%s",    
						StPicoBesNetParticleHists::kQARunNames[idx]))))->Fill(runIdx, aTrack[1]);
  (static_cast<TProfile*>(list->FindObject(Form("pDca_%s",   
						StPicoBesNetParticleHists::kQARunNames[idx]))))->Fill(runIdx, aTrack[4]);  
  
  (static_cast<TProfile*>(list->FindObject(Form("pPt_%s_%s", 
						StPicoBesNetParticleHists::kParticleName[idxSign], 
						StPicoBesNetParticleHists::kQARunNames[idx]))))->Fill(runIdx, aTrack[1]);
  (static_cast<TH2D*>(list->FindObject(Form("hPt_%s_%s",     
					    StPicoBesNetParticleHists::kParticleName[idxSign], 
					    StPicoBesNetParticleHists::kQARunNames[idx]))))->Fill(runIdx, aTrack[1]);
  
}
#endif
