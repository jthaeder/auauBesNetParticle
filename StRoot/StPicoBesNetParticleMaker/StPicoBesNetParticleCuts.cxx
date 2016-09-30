#include <limits>

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEvent.h"

#include "StPicoBesNetParticleCuts.h"
#include "StRefMultCorr/StRefMultCorr.h"

ClassImp(StPicoBesNetParticleCuts)

// _________________________________________________________
StPicoBesNetParticleCuts::StPicoBesNetParticleCuts() : StPicoCutsBase("BesCutsBase"), 
  mRefmultCorr(NULL),
  mEventStatMax(10),
  mAnalysisIdx(0),
  mVxShift(0.), mVyShift(0.),
  mVrMax(std::numeric_limits<float>::max()),
  mNCentralityBinsMax(std::numeric_limits<float>::max()),
  mNTOFMatchMin(0),
  mUseCutRefMultVsNTOFMatch(kFALSE),
  mCutRefMultVsNTOFMatchA(1.),
  mCutRefMultVsNTOFMatchB(0.),
  mPtMidPoint(1.),
  mNHitsDedxMin(1.) {
  // -- default constructor

  for (Int_t idx = 0; idx < kPicoPIDMax; ++idx) {
    mEtaRange[idx][0] = std::numeric_limits<float>::lowest();
    mEtaRange[idx][1] = std::numeric_limits<float>::max();
    mYRange[idx][0] = std::numeric_limits<float>::lowest();
    mYRange[idx][1] = std::numeric_limits<float>::max();
   }
}

// _________________________________________________________
StPicoBesNetParticleCuts::StPicoBesNetParticleCuts(const Char_t *name) : StPicoCutsBase(name),
  mRefmultCorr(NULL),
  mEventStatMax(10),
  mAnalysisIdx(0),
  mVxShift(0.), mVyShift(0.),
  mVrMax(std::numeric_limits<float>::max()),
  mNCentralityBinsMax(std::numeric_limits<float>::max()),
  mNTOFMatchMin(0),
  mUseCutRefMultVsNTOFMatch(kFALSE),
  mCutRefMultVsNTOFMatchA(1.),
  mCutRefMultVsNTOFMatchB(0.),
  mPtMidPoint(1.),
  mNHitsDedxMin(1.) {
  // -- constructor

  for (Int_t idx = 0; idx < kPicoPIDMax; ++idx) {
    mEtaRange[idx][0] = std::numeric_limits<float>::lowest();
    mEtaRange[idx][1] = std::numeric_limits<float>::max();
    mYRange[idx][0] = std::numeric_limits<float>::lowest();
    mYRange[idx][1] = std::numeric_limits<float>::max();
   }
}

const Int_t   StPicoBesNetParticleCuts::kNEventStat    = 9;
const Char_t* StPicoBesNetParticleCuts::kEventStatNames[] =  {"all", "bad run", 
							      "trigger", "#it{v}_{z} < #it{v}_{z}^{max}", 
							      "#it{v}_{z}-#it{v}_{z}^{vpd} < 3cm", 
							      "shifted vtx <1", 
							      "centrality", 
							      "nTOFMatch>2", "nTOFMatch cut"};


// _________________________________________________________
StPicoBesNetParticleCuts::~StPicoBesNetParticleCuts() { 
  // destructor
  
}

// =======================================================================


// _________________________________________________________
void StPicoBesNetParticleCuts::init(StRefMultCorr* refmultCorr) {
  // -- initial cuts class

  // -- max number of event cuts
  mEventStatMax = kNEventStat;
				
  // -- set pointer to refMultCorr
  mRefmultCorr = refmultCorr;
}

// _________________________________________________________
const Char_t* StPicoBesNetParticleCuts::getCutsTitle() {
  // retrun a string of cuts for this setting 
  
  Int_t particleIdx = kPion;
  if (mAnalysisIdx == kNetProton)  
    particleIdx = kProton;
  else if (mAnalysisIdx == kNetKaon)  
    particleIdx = kKaon;

  TString sTitle("");
  sTitle += (mAnalysisIdx == kNetCharge) ? Form("%.2f<#eta<%.2f", mEtaRange[particleIdx][0], mEtaRange[particleIdx][1]) :
    Form("%.2f<y<%.2f", mYRange[particleIdx][0], mYRange[particleIdx][1]);

  sTitle += Form(" #it{p}_{T} [%.1f,%.1f]", mPtRange[particleIdx][0], mPtRange[particleIdx][1]);   

  return sTitle.Data();
}

// _________________________________________________________
Bool_t StPicoBesNetParticleCuts::isGoodBesEvent(StPicoDst const * const picoDst, Int_t *aEventCuts) {
  // -- method to check if good BES  event
  //    calls is  StPicoCutsBase::isGoodEvent
  
  // -- get picoDst event
  StPicoEvent* picoEvent = mPicoDst->event();

  // -- Check cuts implenented in base class
  //    -- 0 - before event cuts
  //    not used -- 1 - is bad run
  //    -- 2 - No Trigger fired
  //    -- 3 - Vertex z outside cut window
  //    -- 4 - Vertex z - vertex_z(vpd) outside cut window
  isGoodEvent(picoDst, aEventCuts); 
  
  // -- bad runs / trigger 
  // ------------------------------------------------------------------
  
  // -- 1 - bad run
  Int_t cutIdx = 1;
  if (mRefmultCorr->isBadRun(picoEvent->runId())) 
    aEventCuts[cutIdx] = 1;
  
  // -- 5 - Vr Cuts - shifted vertex cut
  cutIdx = 5;
  Double_t shiftedVtx = TMath::Sqrt((picoEvent->primaryVertex().x()-mVxShift)*
				    (picoEvent->primaryVertex().x()-mVxShift) + 
				    (picoEvent->primaryVertex().y()-mVyShift)*
				    (picoEvent->primaryVertex().y()-mVyShift));
  if (shiftedVtx >= mVrMax) 
    aEventCuts[cutIdx] = 1;
    
  // -- 6 - cut for centrality range
  ++cutIdx;
  Int_t centralityBin = 8 - mRefmultCorr->getCentralityBin9();

  if (centralityBin >= mNCentralityBinsMax || centralityBin == -1)
    aEventCuts[cutIdx] = 1;
    
  // -- 7 - cut for nTOFMatch >2
  ++cutIdx;
  if (picoEvent->nBTOFMatch() <= mNTOFMatchMin)
    aEventCuts[cutIdx] = 1;
  
  // -- 8 - cut for nTOFMatch vs nRefMult
  ++cutIdx;
  if (mUseCutRefMultVsNTOFMatch && 
      mCutRefMultVsNTOFMatchA*picoEvent->refMult() - mCutRefMultVsNTOFMatchB > picoEvent->nBTOFMatch() ) 
    aEventCuts[cutIdx] = 1;

  // ------------------------------------------------------------------

  // -- check if event is rejected
  Bool_t isGoodEvent = true;
  for (Int_t ii = 0; ii < mEventStatMax; ++ii)
    if  (aEventCuts[ii])
      isGoodEvent = false;
        
  return isGoodEvent;
}

// _________________________________________________________
Bool_t StPicoBesNetParticleCuts::isGoodBesTrack(StPicoTrack const * const trk, Double_t *aTrack) {
  // -- method to check if good BES track
  //    calls is  StPicoCutsBase::isGoodTrack

  // -- Set which particle to use 
  // ------------------------------------------------------------------
  Int_t particleIdx = kPion;
  if (mAnalysisIdx == kNetProton)  
    particleIdx = kProton;
  else if (mAnalysisIdx == kNetKaon)  
    particleIdx = kKaon;
  
  // -- Get Variables
  // ------------------------------------------------------------------
  aTrack[1] = trk->gPt();

  // -- Eta for Net-Charge
  if (mAnalysisIdx == 0)  
    aTrack[2] = trk->pMom().pseudoRapidity();

  // -- Rapidity for identified particles
  else {
    Float_t ecm = TMath::Sqrt(trk->pMom().mag2() + mHypotheticalMass2[particleIdx]);
    aTrack[2] = 0.5 * TMath::Log( (ecm +trk->pMom().z()) / (ecm-trk->pMom().z()) );
  }

  aTrack[3] = trk->charge();

  StPhysicalHelixD helix = trk->dcaGeometry().helix();
  helix.moveOrigin(helix.pathLength(mPrimVtx));
  
  aTrack[4] = (mPrimVtx - helix.origin()).mag();   // dca
  
  aTrack[5] = Double_t(trk->nHitsDedx());
  aTrack[6] = Double_t(TMath::Abs(trk->nHitsFit()));
  aTrack[7] = Double_t(trk->nHitsMax());
  aTrack[8] = Double_t((1+aTrack[6])/(1+aTrack[7]));   // ratio 

  // is accepted  aTrack[9] = 
  
  if (mAnalysisIdx == 1)  
    aTrack[10] = TMath::Abs(trk->nSigmaProton());
  else if (mAnalysisIdx == 2)
    aTrack[10] = TMath::Abs(trk->nSigmaKaon());
	 
  
  // -- is track accepted flag - kinematics
  // ------------------------------------------------------------------
  Bool_t isTrackAcceptedKin = (aTrack[1] >= mPtRange[particleIdx][0] 
			       && aTrack[1] < mPtRange[particleIdx][1]) ? kTRUE : kFALSE;


  if (isTrackAcceptedKin && mAnalysisIdx == 0) {
    isTrackAcceptedKin = (aTrack[2] >= mEtaRange[particleIdx][0]
			  && aTrack[2] < mEtaRange[particleIdx][1]) ? kTRUE : kFALSE;
  }
  else if (isTrackAcceptedKin && mAnalysisIdx > 0) {
    isTrackAcceptedKin = (aTrack[2] >= mYRange[particleIdx][0] 
			  && aTrack[2] < mYRange[particleIdx][1]) ? kTRUE : kFALSE;

  }
   
  // -- is track accepted flag - clusters/dca
  // ------------------------------------------------------------------
  Bool_t isTrackAcceptedCut = (aTrack[6]    > mNHitsFitMin
			       && aTrack[4] < mPrimaryDCAtoVtxMax
			       && aTrack[5] > mNHitsDedxMin
			       && aTrack[8] > mNHitsFitnHitsMax) ? kTRUE : kFALSE;
            
  // -- is track accepted flag - PID
  // ------------------------------------------------------------------
  Bool_t isTrackAcceptedPid = kTRUE;

  if (mAnalysisIdx > 0) {
    // -- PID for net-proton and net-kaon
    Bool_t isTrackAcceptedPidTPC = isTPCHadron(trk, particleIdx); 

    Bool_t isTrackAcceptedPidTOF = (aTrack[1] > mPtMidPoint) ? isTOFHadron(trk, getTofBeta(trk), particleIdx) : kTRUE;

    /* FIX ME
    Double_t mSquareRange[3][2]          = { {-1, -1}, { 0.6, 1.2}, { 0.15, 0.4} };
    Float_t beta    = pico->Tracks_mBTofBeta[idxTrack]/20000.;
    beta = (beta <= 0) ? -999 : beta;
    
    Float_t mSquare = (beta == -999. || beta <= 1.e-5) ? -999. : pow(pcm,2)*(pow(1/beta,2)-1);
    */

    
    isTrackAcceptedPid = (isTrackAcceptedPidTPC && isTrackAcceptedPidTOF);
  }

  else {
    // -- is track Spallation proton/anti-proton for net-charge
    Bool_t isTrackSpallationProton = (aTrack[1] >  mPtRange[kProton][0] && aTrack[1] < mPtRange[kProton][1]
				      && TMath::Abs(trk->nSigmaProton()) <  mTPCNSigmaMax[kProton]) ? kTRUE : kFALSE;  
        
    isTrackAcceptedPid = !isTrackSpallationProton;
  }
    
  // -->> is track accepted  - clusters/dca && kinematics && PID
  Bool_t isTrackAccepted = (isTrackAcceptedKin && isTrackAcceptedCut && isTrackAcceptedPid);    
  aTrack[9] = Double_t(isTrackAccepted);

  return isTrackAccepted;
}


