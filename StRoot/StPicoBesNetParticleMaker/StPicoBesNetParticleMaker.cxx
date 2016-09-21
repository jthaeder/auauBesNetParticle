#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoPrescales/StPicoPrescales.h"

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
  mOutputFileList(NULL) 

  mPicoDstMaker(picoMaker), mPicoEvent(NULL), 
  mEventCounter(0), 

  mOrder(8), mNNp(3), mNp(NULL), mFact(NULL)

  mEneryIdx(-1), mAnalysisIdx(-1), 
  mQaMode(0), 

  mUseModeChargeSeparation(0) {
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
  
  

  for (int idxCent = 0; idxCent <= 9; idxCent++) {
    if (mFact[i][j][k][h])
      delete mFact[i][j][k][h];
    mFact[i][j][k][h] = NULL;
    
    for (int h = 0; h <= mOrder; h++) { 
      if (mFact[i][j][k])
	delete[] mFact[i][j][k];
      mFact[i][j][k] = NULL;
      
      for (int k = 0; k <= mOrder; k++) { 
	if (mFact[i][j])
	  delete[] mFact[i][j];
	mFact[i][j] = NULL;
      
	for (int j = 0; j <= mOrder; j++) {
	  if (mFact[i])
	    delete[] mFact[i];
	  mFact[i] = NULL;
	}
      }
    }
  }

  if (mFact)
    delete[] mFact;
  mFact = NULL;

}

// _________________________________________________________
Int_t StPicoBesNetParticleMaker::Init() {
  // -- Inhertited from StMaker 
  //    Initializes everything

  // -- Check for cut class
  if (!mBesCuts)
    mBesCuts = new StPicoBesNetParticleCuts;
  mBesCuts->init();

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

#if 0 ... MOVE SOMEWHERE ELSE
  // -- get general tile for histograms
  TString sTitle("");
  sTitle += (!analysisIdx) ? Form("%.2f<#eta<%.2f", etaAbsRange[analysisIdx][0], etaAbsRange[analysisIdx][1]) : 
    Form("%.2f<y<%.2f", yAbsRange[analysisIdx][0], yAbsRange[analysisIdx][1]); 		   
  sTitle += Form(" #it{p}_{T} [%.1f,%.1f]", ptRange[analysisIdx][0], ptRange[analysisIdx][1]);
#endif

  // -- Get StRefMultCorr
  mRefMultCorr = new StRefMultCorr(StBesHists::kNameRefMultShort[mAnalysisIdx]);

  // -- Initialize histogram class
  mBesHists = new StBesHists(Form("besHists_%s",GetName()));
  mBesHists->init(mOutList, sTitle.Data());     // -- no title yet

  // -- Initialize net particle variables
  initNetParticle();

  TH1::AddDirectory(oldStatus);

  cout << StBesHists::kName[analysisIdx] << ": .... " << sTitle << " useChargeSeparation "<<  useModeChargeSeparation << endl;

  // -- reset event to be in a defined state
  resetEvent();

  return kStOK;
}

// _________________________________________________________
Int_t initNetParticle() {
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
  fijkhList->SetName(Form("f%sFijkh", StBesHists::kName[mAnalysisIdx]));
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
  
   return;
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

  if (setupEvent()) {
    //    UInt_t nTracks = mPicoDst->numberOfTracks();

    // -- call method of daughter class
    iReturn = MakeNetParticle();

    // -- fill basic event histograms - for good events
    ///     mBesHists->fillGoodEventHists(*mPicoEvent);

  } // if (setupEvent()) {
    
  // -- fill basic event histograms - for all events
  /// mBEesists->fillEventHists(*mPicoEvent);
  
  return (kStOK && iReturn);
}

Int_t StPicoBesNetParticleMaker::MakeBES() {
  
  UInt_t nTracks = mPicoDst->numberOfTracks();   // XX member?
    // ------------------------------------------------------------------
    // -- Track loop - track multiplicity
    // ------------------------------------------------------------------
    for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
      StPicoTrack* trk = mPicoDst->track(iTrack);      


      /*
      Float_t pxcm    = pico->Tracks_mPMomentum_mX1[idxTrack];
      Float_t pycm    = pico->Tracks_mPMomentum_mX2[idxTrack];
      Float_t pzcm    = pico->Tracks_mPMomentum_mX3[idxTrack];
      if (pxcm == 0 && pycm == 0 && pzcm == 0) 
	continue;
      
      Float_t pt      = TMath::Sqrt(pxcm*pxcm + pycm*pycm);
      Float_t pcm     = TMath::Sqrt(pt*pt + pzcm*pzcm);
      Float_t ecm     = TMath::Sqrt(pcm*pcm + masses[analysisIdx]*masses[analysisIdx]);
      
      Float_t eta     = 0.5 * TMath::Log((pcm+pzcm)/(pcm-pzcm));
      Float_t DCA     = pico->Tracks_mGDca[idxTrack]/1000.;

      Float_t y       = 0.5 * TMath::Log((ecm+pzcm)/(ecm-pzcm));

      Int_t nHitsDedx = pico->Tracks_mNHitsDedx[idxTrack];
      Int_t nHitsFit  = pico->Tracks_mNHitsFit[idxTrack];
      Int_t nFitPoss  = pico->Tracks_mNHitsMax[idxTrack];
      
      Float_t ratio   = (1+fabs(nHitsFit))/(1+nFitPoss);
      
      Float_t sign    = (nHitsFit > 0) ? +1 : -1;
      
      Float_t nSigma[3];
      nSigma[0]       = pico->Tracks_mNSigmaProton[idxTrack]/100.;
      nSigma[1]       = pico->Tracks_mNSigmaProton[idxTrack]/100.;
      nSigma[2]       = pico->Tracks_mNSigmaKaon[idxTrack]/100.;
           
      Float_t beta    = pico->Tracks_mBTofBeta[idxTrack]/20000.;
      beta = (beta <= 0) ? -999 : beta;
            
      Float_t mSquare = (beta == -999. || beta <= 1.e-5) ? -999. : pow(pcm,2)*(pow(1/beta,2)-1);
      */

      // -- is in RefMult flag
      Bool_t isInRefMult = (TMath::Abs(eta) > etaAbsRangeRefMult[analysisIdx][0] 
			    && TMath::Abs(eta) <= etaAbsRangeRefMult[analysisIdx][1]  
			    && TMath::Abs(nHitsFit) > nHitsFitRefMult[analysisIdx]  
			    && DCA < dcaMaxRefMult[analysisIdx]) ? kTRUE : kFALSE;
      
      // -- is track accepted flag - kinematics
      Bool_t isTrackAcceptedKin = (pt > ptRange[analysisIdx][0] && pt < ptRange[analysisIdx][1]) ? kTRUE : kFALSE;

      if (isTrackAcceptedKin && analysisIdx == 0) 
	isTrackAcceptedKin = (eta > etaAbsRange[analysisIdx][0] && eta < etaAbsRange[analysisIdx][1]) ? kTRUE : kFALSE;
      else if (isTrackAcceptedKin && analysisIdx > 0) 
	isTrackAcceptedKin = (y > yAbsRange[analysisIdx][0] && y < yAbsRange[analysisIdx][1]) ? kTRUE : kFALSE;

      // -- is track accepted flag - clusters/dca
      Bool_t isTrackAcceptedCut = (TMath::Abs(nHitsFit) > nHitsFitMin[analysisIdx] 
				   && DCA < dcaMax[analysisIdx] 
				   && nHitsDedx > nHitsDedxMin[analysisIdx] 
				   && ratio > ratioNHitsFitNFitPossMin[analysisIdx]) ? kTRUE : kFALSE;
      
      // -- is track accepted flag - PID
      Bool_t isTrackAcceptedPid = kTRUE;
      
      if (analysisIdx > 0) {
	// -- PID for net-proton . net-kaon
	Bool_t isTrackAcceptedPidTPC = (TMath::Abs(nSigma[analysisIdx]) < nSigmaMax[analysisIdx]) ? kTRUE : kFALSE;
	Bool_t isTrackAcceptedPidTOF = kTRUE;

	if (pt > ptMidPoint[analysisIdx])
	  isTrackAcceptedPidTOF = (mSquare > mSquareRange[analysisIdx][0] 
				   && mSquare < mSquareRange[analysisIdx][1]) ? kTRUE : kFALSE;
	
	isTrackAcceptedPid = (isTrackAcceptedPidTPC && isTrackAcceptedPidTOF);
      }
      else {
	// -- is track Spallation proton/anti-proton
	Bool_t isTrackSpallationProton = (pt > ptRangeSpallation[0] && pt < ptRangeSpallation[1] 
					  && TMath::Abs(nSigma[0]) < nSigmaProtonMaxSpallation)  ? kTRUE : kFALSE;
	
	isTrackAcceptedPid = isTrackSpallationProton;
      }
   
      // -->> is track accepted  - clusters/dca && kinematics && PID
      Bool_t isTrackAccepted = (isTrackAcceptedKin && isTrackAcceptedCut && isTrackAcceptedPid);

 
      // -- fill ThnSparse - tracks
      // ------------------------------------------------------------------
      Double_t aTrack[13] = {Double_t(centrality), pt, eta, sign, DCA,
			     Double_t(nHitsDedx), Double_t(TMath::Abs(nHitsFit)), Double_t(nFitPoss), ratio, 
			     Double_t(isInRefMult), Double_t(isTrackAccepted), nSigma[analysisIdx], pcm};
 
#if TRACK_THN      
      fHnTrackUnCorr->Fill(aTrack);
#endif

      FillTrackHists(aTrack, 0);
      FillRunByRunTrackHists(aTrack, isBadRun, 0, runIdx);

      // -- reject track
      // ------------------------------------------------------------------
      if (!isTrackAccepted)
	continue;
      
      FillTrackHists(aTrack, 1);
      FillRunByRunTrackHists(aTrack, isBadRun, 1, runIdx);

      // ------------------------------------------------------------------
      // -- Add up for event multiplicity
      // ------------------------------------------------------------------
      //  idxPart = 0 -> anti particle
      //  idxPart = 1 -> particle
      Int_t idxPart    = (sign < 0) ? 0 : 1;
      Int_t idxEtaSign = (eta  < 0) ? 0 : 1;

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
      if (useModeChargeSeparation == 0) 
	++count;
      else if ( (useModeChargeSeparation == 1) && (idxPart == idxEtaSign) ) 
	++count;
      else if ( (useModeChargeSeparation == 2) && (idxPart != idxEtaSign) ) 
	++count;

      // -- in full pt Range
      fNp[0][idxPart] += count;
      
      // -- divide in 2 parts
      if (pt < ptMidPoint[analysisIdx])
	fNp[1][idxPart] += count;
      else
	fNp[2][idxPart] += count;
    } // for (Int_t idxTrack = 0; idxTrack < nTracks; idxTrack++)  {
}


// _________________________________________________________
bool StPicoBesNetParticleMaker::setupEvent() {
  // -- fill members from pico event, check for good eventa and fill event statistics

  // -- Reset event to be in a defined state
  resetEvent();

  mPicoEvent = mPicoDst->event();
  
  mBField = mPicoEvent->bField();
  mPrimVtx = mPicoEvent->primaryVertex();
  
  // -- arry with one entry for every cut
  //    0 - event passes cut
  //    1 - event doesn't pass cut
  int aEventStat[mBesCuts->eventStatMax()];
  
  bool bResult = mBesCuts->isGoodEvent(mPicoDst, aEventStat);

  // -- fill event statistics histograms
  fillEventStats(aEventStat);

  return bResult;
}

// _________________________________________________________
void StPicoBesNetParticleMaker::initializeEventStats() {
  // -- Initialize event statistics histograms
  
  const char *aEventCutNames[] = {"all", "good run", "trigger", "#it{v}_{z}", "#it{v}_{z}-#it{v}^{VPD}_{z}", "accepted"};

  mOutList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", mBesCuts->eventStatMax(), -0.5, mBesCuts->eventStatMax()-0.5));
  TH1F *hEventStat0 = static_cast<TH1F*>(mOutList->Last());

  mOutList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", mBesCuts->eventStatMax(), -0.5, mBesCuts->eventStatMax()-0.5));
  TH1F *hEventStat1 = static_cast<TH1F*>(mOutList->Last());

  for (unsigned int ii = 0; ii < mBesCuts->eventStatMax(); ii++) {
    hEventStat0->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
    hEventStat1->GetXaxis()->SetBinLabel(ii+1, aEventCutNames[ii]);
  }

  //  hEventStat0->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", aCentralityMaxNames[9-1]));
  //  hEventStat1->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", aCentralityMaxNames[9-1]));
}

//________________________________________________________________________
void StPicoBesNetParticleMaker::fillEventStats(int *aEventStat) {
  // -- Fill event statistics 

  TH1F *hEventStat0 = static_cast<TH1F*>(mOutList->FindObject("hEventStat0"));
  TH1F *hEventStat1 = static_cast<TH1F*>(mOutList->FindObject("hEventStat1"));

  for (unsigned int idx = 0; idx < mBesCuts->eventStatMax() ; ++idx) {
    if (!aEventStat[idx])
      hEventStat0->Fill(idx);
  }
  
  for (unsigned int idx = 0; idx < mBesCuts->eventStatMax(); ++idx) {
    if (aEventStat[idx])
      break;
    hEventStat1->Fill(idx);
  }
}

