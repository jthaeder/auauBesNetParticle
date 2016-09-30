#ifndef StPicoBesNetParticleMaker_h
#define StPicoBesNetParticleMaker_h

#include "StChain/StMaker.h"
#include "StarClassLibrary/StLorentzVectorF.hh"

/* **************************************************
 *  Base class for HF analysis
 *
 *  Usage:
 *  - Implement in daughter class
 *     InitBES()
 *     MakeBES()
 *     ClearBES()
 *     FinishBES()
 *  - Do not ovewrite Init, Make, Clear, Finish which are inhertited from StMaker
 *  - Set StBESCuts class via setBESBaseCuts(...) in run macro
 *
 * **************************************************
 *
 *  Initial Authors:
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *          **Jochen Thaeder  (jmthader@lbl.gov)
 *
 *  Contributing Authors
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

class TTree;
class TFile;
class TChain;
class TProfile;

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;

class StRefMultCorr;

class StPicoBesNetParticleCuts;
class StPicoBesNetParticleHists;

class StPicoBesNetParticleMaker : public StMaker 
{
  public:
    StPicoBesNetParticleMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName);

    virtual ~StPicoBesNetParticleMaker();
    
    void setBesCuts(StPicoBesNetParticleCuts* cuts) { mBesCuts = cuts; }
    void setEnergyIdx(Int_t i)                      { mEnergyIdx = i; }
    void setAnalysisIdx(Int_t i)                    { mAnalysisIdx = i; }
    void setQaMode(Int_t i)                         { mQaMode = i; }

    Int_t Init();
    Int_t Make();
    void  Clear(Option_t *opt="");
    Int_t Finish();

  protected:
    // -- protected members ------------------------

    StPicoDst                 *mPicoDst;

    StPicoBesNetParticleCuts  *mBesCuts;
    StPicoBesNetParticleHists *mBesHists;

    float                      mBField;
    StThreeVectorF             mPrimVtx;

    TList                     *mOutList;

  private:

    // -- reset everything on start of new event
    void  resetEvent();
    
    // -- setup event - and do event cuts
    bool  setupEvent();
    
    // -- initialize event stat histograms
    void  initializeEventStats();

    // -- fill event stat histograms
    void  fillEventStats(int *aEventStat);

    // -- initialize net particle arrays
    Int_t initNetParticle();

    // -- Runs the actual analysis
    Int_t makeNetParticle();
    
    // -- Helper method for factorial moments calculation
    Double_t NN(Double_t num, Int_t order);


    // -- private members ------------------------
    TString         mOutputFileBaseName; // base name for output files
                                         //   for histList -> <mOutputFileBaseName>.GetName().root

    TFile*          mOutputFileList;     // ptr to file saving the list of histograms

    StPicoDstMaker* mPicoDstMaker;       // ptr to picoDst maker

    StPicoEvent*    mPicoEvent;          // ptr to picoDstEvent

    StRefMultCorr*  mRefMultCorr;        // ptr to refMultCorr


    // -----------------------------------------------------------------------

    Int_t           mOrder;              // Max order of higher order distributions

    Int_t           mNNp;                // N sets of arrays of particle/anti-particle counts
                                         //   0 is all 
                                         //   1,2 for arbitrary subset

    Int_t          **mNp;                // Array of particle/anti-particle counts

    TProfile   ******mFact;              // Local Factorial Moment array of TProfiles
                                         // for order 8 (+1 as starts with idx = 1) + 10 cent bins

    Int_t            mEnergyIdx;         // energyIdx
                                         // - for different beam energies

    Int_t            mAnalysisIdx;       // 0 - net-charge
                                         // 1 - net-proton
                                         // 2 - net-kaon

    Int_t            mQaMode;            // 0 - NO QA MODE (default)
                                         // 1 - QA MODE
    // =======================================================================
    Int_t            mUseModeChargeSeparation;  // 0 - off
                                                // 1 - positive charge -> positive eta
                                                // 2 - positive charge -> negative eta


    // -----------------------------------------------------------------------

    Int_t           mNRefMultX;          // refMuly<X> of event according to analysis
                                         //  <X> - 2 for Charge
                                         //  <X> - 3 for Proton
                                         //  <X> - 4 for Kaon

    Int_t           mNRefMultXCorr;      // vz corrected refMuly<X> of event according to analysis
                                         //  <X> - 2 for Charge
                                         //  <X> - 3 for Proton
                                         //  <X> - 4 for Kaon

    Int_t           mCentralityBin;      // Centrality Bin

    ClassDef(StPicoBesNetParticleMaker, 0)
};

#endif
