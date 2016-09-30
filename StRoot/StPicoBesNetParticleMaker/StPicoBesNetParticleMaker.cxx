#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TProfile.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"

#include "StRefMultCorr/StRefMultCorr.h"

#include "StPicoBesNetParticleCuts.h"
#include "StPicoBesNetParticleHists.h"
#include "StPicoBesNetParticleMaker.h"

ClassImp(StPicoBesNetParticleMaker)

// _________________________________________________________
StPicoBesNetParticleMaker::StPicoBesNetParticleMaker(char const* name, StPicoDstMaker* picoMaker, 
			       char const* outputBaseFileName) :
  StMaker(name), mPicoDst(NULL), mBesCuts(NULL), mBesHists(NULL),

  mBField(0.), mPrimVtx(NULL),
  mOutList(NULL),

  mOutputFileBaseName(outputBaseFileName), 
  mOutputFileList(NULL),

  mPicoDstMaker(picoMaker), mPicoEvent(NULL), 

  mOrder(8), mNNp(3), mNp(NULL), mFact(NULL),

  mEnergyIdx(-1), mAnalysisIdx(-1), 
  mQaMode(0), 

  mUseModeChargeSeparation(0),

  mNRefMultX(-1), mNRefMultXCorr(-1), mCentralityBin(-1) {
  // -- constructor
}

// _________________________________________________________
StPicoBesNetParticleMaker::~StPicoBesNetParticleMaker() {
   // -- destructor 
  
  if (mBesCuts)
    delete mBesCuts;
  mBesCuts = NULL;

  if (mBesHists)
    delete mBesHists;
  mBesHists = NULL;

  for (Int_t ii = 0 ; ii < mNNp; ++ii) {
    if (mNp[ii])
      delete[] mNp[ii];
    mNp[ii] = NULL;
  }
  if (mNp)
    delete[] mNp;
  mNp = NULL;

  for (int i = 0; i <= mOrder; i++) { 
    for (int j = 0; j <= mOrder; j++) { 
      for (int k = 0; k <= mOrder; k++) { 
	for (int h = 0; h <= mOrder; h++) {   
	  for (int idxCent = 0; idxCent <= 9; idxCent++) {
	    if (mFact[i][j][k][h][idxCent])
	      delete mFact[i][j][k][h][idxCent];
	    mFact[i][j][k][h][idxCent] = NULL;
	  }
	  if (mFact[i][j][k][h])
	    delete mFact[i][j][k][h];
	  mFact[i][j][k][h] = NULL;
	} // for (int h = 0; h <= mOrder; h++) {   
	if (mFact[i][j][k])
	  delete[] mFact[i][j][k];
	mFact[i][j][k] = NULL;
      } // for (int k = 0; k <= mOrder; k++) { 
      if (mFact[i][j])
	delete[] mFact[i][j];
      mFact[i][j] = NULL;
    } // for (int j = 0; j <= mOrder; j++) { 
    if (mFact[i])
      delete[] mFact[i];
    mFact[i] = NULL;
  } // for (int i = 0; i <= mOrder; i++) { 

  if (mFact)
    delete[] mFact;
  mFact = NULL;

}

// _________________________________________________________
Int_t StPicoBesNetParticleMaker::Init() {
  // -- Inhertited from StMaker 
  //    Initializes everything

  // -- Check if analysis was set up properly in macro
  if (mEnergyIdx < 0) {
    LOG_WARN << "Energy Idx not set." << endm;
    return kStWarn;
  }

  if (mAnalysisIdx < 0) {
    LOG_WARN << "Analysis Idx not set." << endm;
    return kStWarn;
  }

  // -- File which holds list of histograms
  mOutputFileList = new TFile(Form("%s.%s.root", mOutputFileBaseName.Data(), GetName()), "RECREATE");
  mOutputFileList->SetCompressionLevel(1);

  // -- Disable automatic adding of objects to file
  bool oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  // -- Add list which holds all histograms  
  mOutList = new TList();
  mOutList->SetName(GetName());
  mOutList->SetOwner(true);

  // -- Get RefMult Corr
  mRefMultCorr = new StRefMultCorr(StPicoBesNetParticleHists::kNameRefMultShort[mAnalysisIdx]);
  
  // -- Check for cut class
  if (!mBesCuts)
    mBesCuts = new StPicoBesNetParticleCuts;

  // -- Initialize cuts class
  mBesCuts->init(mRefMultCorr);

  // -- Initialize histogram class
  mBesHists = new StPicoBesNetParticleHists(Form("besHists_%s",GetName()));
  mBesHists->init(mOutList, mBesCuts->getCutsTitle(),
		  mAnalysisIdx, mBesCuts->getNCentralityBinsMax()); 

  // -- Initialize net particle variables
  initNetParticle();

  TH1::AddDirectory(oldStatus);

  cout << StPicoBesNetParticleHists::kName[mAnalysisIdx] << ": .... " << mBesCuts->getCutsTitle() << endl;

  // -- reset event to be in a defined state
  resetEvent();

  return kStOK;
}

// _________________________________________________________
Int_t StPicoBesNetParticleMaker::initNetParticle() {
  // -- init analysis members

  // -- Counting array
  mNp = new Int_t*[mNNp];
  for (Int_t ii = 0 ; ii < mNNp; ++ii)
    mNp[ii] = new Int_t[2];

  Int_t nCent = 10;    // -- ## nCent to be added

  // --  Factorial moment array of TProfiles
  mFact = new TProfile*****[mOrder+1];

  for (int i = 0; i <= mOrder; i++) { 
    mFact[i] = new TProfile****[mOrder+1];

    for (int j = 0; j <= mOrder; j++) {
      mFact[i][j] = new TProfile***[mOrder+1];

      for (int k = 0; k <= mOrder; k++) { 
	mFact[i][j][k] = new TProfile**[mOrder+1];

	for (int h = 0; h <= mOrder; h++) { 
	  mFact[i][j][k][h] = new TProfile*[nCent];
	}
      }
    }
  }

  // -- Local Factorial Moment array of TProfiles 
  mOutList->Add(new TList);
 
  TList *fijkhList = static_cast<TList*>(mOutList->Last());
  fijkhList->SetName(Form("f%sFijkh", StPicoBesNetParticleHists::kName[mAnalysisIdx]));
  fijkhList->SetOwner(kTRUE);

  for (int idxCent = 0; idxCent <= 9; idxCent++) 
    for (int i = 0; i <= mOrder; i++) 
      for (int j = 0; j <= mOrder; j++) 
	for (int k = 0; k <= mOrder; k++) 
	  for (int h = 0; h <= mOrder; h++) 
	    if ((i+j+k+h) <= mOrder) {
	      fijkhList->Add(new TProfile(Form("f%d%d%d%d_Cent%d", i, j, k, h, idxCent),
					  Form("f%d%d%d%d_Cent%d", i, j, k, h, idxCent), 
					  1001, -0.5, 1000.5));
	      mFact[i][j][k][h][idxCent] = static_cast<TProfile*>(fijkhList->Last());
	    }
  
   return 0;
}

// _________________________________________________________
Int_t StPicoBesNetParticleMaker::Finish() {
  // -- Inhertited from StMaker 

  mOutputFileList->cd();
  mOutList->Write(mOutList->GetName(), TObject::kSingleKey);
  mOutputFileList->Close();

  return kStOK;
}

// _________________________________________________________
void StPicoBesNetParticleMaker::resetEvent() {
  // -- reset event

  // -- Reset N particles/anti-particles
  for (Int_t ii = 0; ii < mNNp; ++ii) 
    for (Int_t jj = 0; jj < 2; ++jj)
      mNp[ii][jj] = 0;

}

// _________________________________________________________
void StPicoBesNetParticleMaker::Clear(Option_t *opt) {
  // -- Inhertited from StMaker 

  resetEvent();
}

// _________________________________________________________
Int_t StPicoBesNetParticleMaker::Make() {
  // -- Inhertited from StMaker 

  if (!mPicoDstMaker) {
    LOG_WARN << " StPicoBesNetParticleMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) {
    LOG_WARN << " StPicoBesNetParticleMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  
  Int_t iReturn = kStOK;

  // -- Setup Event and check if event should be processed
  if (setupEvent()) {

    // -- call track loop for analysis
    iReturn = makeNetParticle();

  } // if (setupEvent()) {
    
  
  return (kStOK && iReturn);
}

// _________________________________________________________
Int_t StPicoBesNetParticleMaker::makeNetParticle() {
  // -- Track loop for the analysis
  //    Fill particle hists and factorial moments

  Int_t nTracks = mPicoDst->numberOfTracks(); 

  // ------------------------------------------------------------------
  // -- Track loop - track multiplicity
  // ------------------------------------------------------------------
  for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
    StPicoTrack* trk = mPicoDst->track(idxTrack);      
    
    Bool_t isTrackAccepted = kTRUE;
    
    Double_t aTrack[13];
    if (!trk || ! mBesCuts->isGoodBesTrack(trk, aTrack)) 
      isTrackAccepted = kFALSE;
 
    // -- fill ThnSparse - tracks
    // ------------------------------------------------------------------
    aTrack[0] = Double_t(mCentralityBin);
 
    // -- Fill Track THn
    mBesHists->FillHnTrack(aTrack);

    // -- Fill Track QA Hists
    mBesHists->FillTrackQAHists(aTrack, 0);
    //    mBesHists->FillRunByRunTrackHists(aTrack, isBadRun, 0, runIdx);

    // -- reject track
    // ------------------------------------------------------------------
    if (!isTrackAccepted)
      continue;
    
    mBesHists->FillTrackQAHists(aTrack, 1);
    //    mBesHists->FillRunByRunTrackHists(aTrack, isBadRun, 1, runIdx);

    // ------------------------------------------------------------------
    // -- Add up for event multiplicity
    // ------------------------------------------------------------------
    //  idxPart = 0 -> anti particle
    //  idxPart = 1 -> particle
    Int_t idxPart    = (aTrack[3] < 0) ? 0 : 1;
    Int_t idxEtaSign = (aTrack[2]  < 0) ? 0 : 1;
    
#if USE_RANDOM_EFF      
    // -- discard a random amount of tracks
    if (gRandom->Rndm() > randomEff[idxPart][centrality])
      continue;
#endif
    
    // -- Apply Charge separation 
    //    -> default: off  => 0  
    //    -> positive particles from postive eta / negative particles from negative eta => 1
    //    -> positive particles from negative eta / negative particles from postive eta => 2
    
    Int_t count = 0;
    if (mUseModeChargeSeparation == 0) 
      ++count;
    else if ( (mUseModeChargeSeparation == 1) && (idxPart == idxEtaSign) ) 
      ++count;
    else if ( (mUseModeChargeSeparation == 2) && (idxPart != idxEtaSign) ) 
      ++count;
    
    // -- in full pt Range
    mNp[0][idxPart] += count;
    
    // -- divide in 2 parts
    if (aTrack[1] < mBesCuts->getPtMidPoint())
      mNp[1][idxPart] += count;
    else
      mNp[2][idxPart] += count;
  } 


  // -- Fill histograms
  // ------------------------------------------------------------------
  mBesHists->FillHistSetCent("Dist",       mNp[0], mCentralityBin);
  mBesHists->FillHistSetCent("Dist_lower", mNp[1], mCentralityBin);
  mBesHists->FillHistSetCent("Dist_upper", mNp[2], mCentralityBin);
  
  // ------------------------------------------------------------------
  
  // Double_t aEventRun[] = { Double_t(mPicoEvent->refMult()), 
  // 			   Double_t(mNRefMultX),
  // 			   Double_t(mNp[0][1] - mNp[0][0]), 
  // 			   Double_t(mNp[0][0]), Double_t (mNp[0][1])};
  // FIX ME LATER  FillRunByRunEventHists(aEventRun, isBadRun, runIdx); 
  
  // -- Filling of factorial moments

    for (int i = 0; i <= mOrder; i++) 
      for (int j = 0; j <= mOrder; j++) 
	for (int k = 0; k <= mOrder; k++) 
	  for (int h = 0; h <= mOrder; h++) 
	    if ((i+j+k+h) <= mOrder) {
	      mFact[i][j][k][h][mCentralityBin+1]->Fill(mNRefMultXCorr, NN(mNp[1][1],i) * NN(mNp[2][1],j) * NN(mNp[1][0],k) * NN(mNp[2][0],h));
	      mFact[i][j][k][h][0]->Fill(               mNRefMultXCorr, NN(mNp[1][1],i) * NN(mNp[2][1],j) * NN(mNp[1][0],k) * NN(mNp[2][0],h));
	    }

    // ------------------------------------------------------------------

    return kStOK;
}

//________________________________________________________________________
Double_t NN(Double_t num, Int_t order) {
  return (order == 0) ? 1 : NN(num,order-1)*(num-order+1);
}

// _________________________________________________________
bool StPicoBesNetParticleMaker::setupEvent() {
  // -- fill members from pico event, check for good eventa and fill event statistics

  // -- Reset event to be in a defined state
  resetEvent();

  mPicoEvent = mPicoDst->event();
  
  mBField = mPicoEvent->bField();
  mPrimVtx = mPicoEvent->primaryVertex();
  
  // -- Get refmultX number of event
  if (mAnalysisIdx == 0)
    mNRefMultX = mPicoEvent->refMult2();
  else if (mAnalysisIdx == 1)
    mNRefMultX = mPicoEvent->refMult3();
  else if (mAnalysisIdx == 2)
    mNRefMultX = mPicoEvent->refMult4();
  else
    mNRefMultX = 0;

  // -- Setup refMultCorr
  mRefMultCorr->init(mPicoEvent->runId());
  mRefMultCorr->initEvent(mNRefMultX, mPrimVtx.z());
  
  // -- Correct refmult vz dependent 
  mNRefMultXCorr = mRefMultCorr->getRefMultCorr();

  // -- get CentralityBin
  mCentralityBin = 8 - mRefMultCorr->getCentralityBin9();

  // ----------------------------------------------------
  // -- EVENT CUTS
  // ----------------------------------------------------

  // -- Array with one entry for every cut
  //    0 - event passes cut
  //    1 - event doesn't pass cut
  Int_t aEventStat[StPicoBesNetParticleCuts::kNEventStat];
  
  // -- Check for good event
  bool isGoodEvent = mBesCuts->isGoodBesEvent(mPicoDst, aEventStat);

  // ----------------------------------------------------
  // -- Fill event QA histograms
  // ----------------------------------------------------

  //  needed of the run wise event cuts 
  //Int_t isBadRun = (aEventStat[1]) ? 1 : 0;

  // -- don't cut on bad runs for qaMode
  if (mQaMode)
    aEventStat[1] = 0;

  // -- Fill event statistics histograms
  mBesHists->FillEventStats(aEventStat);


  Bool_t isRejectedWithoutTOF = (aEventStat[0] || aEventStat[1] || aEventStat[2] || aEventStat[3] || 
				 aEventStat[4] || aEventStat[5] || aEventStat[6]) ? kTRUE : kFALSE;
  

  // -- Create arrays  aEvent - aMult
  // ----------------------------------------------------
  Double_t shiftedVtx = TMath::Sqrt((mPrimVtx.x()-mBesCuts->getVxShift())*
				    (mPrimVtx.x()-mBesCuts->getVxShift()) + 
				    (mPrimVtx.y()-mBesCuts->getVyShift())*
				    (mPrimVtx.y()-mBesCuts->getVyShift()));
  

  Double_t aEvent[11] = {Double_t(mCentralityBin), 
			 mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z(),
			 shiftedVtx, 
			 Double_t(mNRefMultX),
			 Double_t(mNRefMultXCorr),
			 
			 Double_t(mPicoDst->numberOfTracks()),

			 Double_t(!isGoodEvent),
			 mPicoEvent->vzVpd(),
			 mPrimVtx.z()-mPicoEvent->vzVpd()};
      
  Double_t aMult[7]  = {Double_t(mCentralityBin),
			Double_t(mPicoEvent->refMult()), 
			Double_t(mNRefMultX),
			Double_t(mNRefMultXCorr),
			
			Double_t(mPicoEvent->numberOfGlobalTracks()),
			Double_t(0),    // FIXME nPrimaryTracks),
			Double_t(mPicoEvent->nBTOFMatch())};
  
  
  // -- Fill Event THn
  mBesHists->FillHnEvent(aEvent);
    
  // -- Fill Event QA Hists
  mBesHists->FillEventQAHists(aEvent, 0);
  
  // -- Fill multiplicty stats - before event cuts
  mBesHists->FillMultiplicityStats(aMult, 0);
  
  // -- Reject Event - but not TOF cuts
  // ------------------------------------------------------------------
  if (isRejectedWithoutTOF) {
    mBesHists->FillMultiplicityStats(aMult, 1);
    return isGoodEvent;
  }
  
  // -- Fill multiplicty stats - after event cuts (but not TOF cuts)
  mBesHists->FillMultiplicityStats(aMult, 2);
  
  // -- Rejected Event
  // ------------------------------------------------------------------
  if (!isGoodEvent) {
    mBesHists->FillMultiplicityStats(aMult, 3);
    return isGoodEvent; 
  }
  
  // -- Fill multiplicty stats
  mBesHists->FillMultiplicityStats(aMult, 4);
  
  // -- Fill eventHists
  mBesHists->FillEventQAHists(aEvent, 1);
  
  // ------------------------------------------------------------------
  
  return isGoodEvent;
}



