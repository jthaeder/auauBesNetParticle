#ifndef StPicoBesNetParticleCuts__h
#define StPicoBesNetParticleCuts__h

/* **************************************************
 *  Cut class for BES analysis
 *  - Based on PicoCuts class 
 *
 *  Initial Authors:  
 *            Xin Dong        (xdong@lbl.gov)
 *          **Jochen Thaeder  (jmthader@lbl.gov)   
 *
 *  Contributing Authors
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

class StRefMultCorr;

#include "StPicoCutsBase/StPicoCutsBase.h"

class StPicoBesNetParticleCuts : public StPicoCutsBase
{
 public:
  
  StPicoBesNetParticleCuts();
  StPicoBesNetParticleCuts(const Char_t *name);
  ~StPicoBesNetParticleCuts();
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  virtual void init() { initBase(); }   // not needed
  
  void init(StRefMultCorr* refmultCorr);
 
  Bool_t isGoodBesEvent(StPicoDst const * const picoDst, Int_t *aEventCuts); 
  Bool_t isGoodBesTrack(StPicoTrack const * const trk,   Double_t *aTrack); 

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- SETTER for CUTS - EVENT CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  void SetVxShift(Float_t f) {mVxShift = f;}
  void SetVyShift(Float_t f) {mVyShift = f;}
  void SetVrMax(Float_t f) {mVrMax = f;}
  
  void SetNCentralityBinsMax(Float_t f) {mNCentralityBinsMax = f;}
  
  void SetNTOFMatchMin(Int_t i) {mNTOFMatchMin = i;} 

  void SetCutRefMultVsNTOFMatch(Float_t f1, Float_t f2) {
    mUseCutRefMultVsNTOFMatch = kTRUE;
    mCutRefMultVsNTOFMatchA = f1;
    mCutRefMultVsNTOFMatchB = f2;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- SETTER for CUTS - TRACK CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  



  

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- GETTER for single CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  const Char_t* getCutsTitle();

  const Int_t   getNCentralityBinsMax() { return mNCentralityBinsMax; }

  const Float_t getPtMidPoint()         { return mPtMidPoint; }


  const Float_t getVxShift()   { return mVxShift; }
  const Float_t getVyShift()   { return mVyShift; }


  static const Int_t    kNEventStat;
  static const Char_t* kEventStatNames[];

 private:
  
  enum analysisType {kNetCharge, kNetProton, kNetKaon, kMaxAnalysisType};

  StPicoBesNetParticleCuts(StPicoBesNetParticleCuts const &);       
  StPicoBesNetParticleCuts& operator=(StPicoBesNetParticleCuts const &); 

  // ------------------------------------------
  // -- Members
  // ------------------------------------------
  
  StRefMultCorr* mRefmultCorr;       // Ptr to refMultCorr

  Int_t          mEventStatMax;      // max number of cuts  
  Int_t          mAnalysisIdx;       // analysisIdx

  // ------------------------------------------
  // -- BES event cuts
  // ------------------------------------------

  Float_t mVxShift;                  // Shift of Vx 
  Float_t mVyShift;                  // Shift of Vy 
  Float_t mVrMax;                    // Max Vr

  Int_t   mNCentralityBinsMax;       // Max CentralityBins from cuts

  Int_t   mNTOFMatchMin;             // Min N for TOFmatch hits - pile up removal

  Bool_t  mUseCutRefMultVsNTOFMatch; // Enable linear cut
  Float_t mCutRefMultVsNTOFMatchA;   // Parameter for linear cut ;  nTOFMatch vs nRefMult
  Float_t mCutRefMultVsNTOFMatchB;   // mCutRefMultVsNTOFMatchA * refMult() - mCutRefMultVsNTOFMatchB > nBTOFMatch

  // ------------------------------------------
  // -- BES track cuts
  // ------------------------------------------

  
  Float_t mPtMidPoint;                // midPoint between two pid areas: TPC pid and TPC+TOF pid

  Float_t mEtaRange[kPicoPIDMax][2];  // eta Range

  Float_t mYRange[kPicoPIDMax][2];    // rapidity range
 
  Int_t   mNHitsDedxMin;              // min number of hits used for dEdx

  ClassDef(StPicoBesNetParticleCuts,0)

};
#endif
