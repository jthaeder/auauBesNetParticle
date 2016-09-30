#ifndef StPicoBesNetParticleHists__h
#define StPicoBesNetParticleHists__h

/* **************************************************
 *  A class to create and save production QA
 *  histograms.
 *
 * **************************************************
 *
 *  Initial Authors:
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *          **Giacomo Contin  (gcontin@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

#include "TNamed.h"
#include "TList.h"

class TH1F;
class TH2F;
class THnSparseD;

class TFile;
class TString;
class StPicoPrescales;

class StPicoEvent;

class StPicoBesNetParticleHists: public TNamed
{
 public:
  StPicoBesNetParticleHists();
  StPicoBesNetParticleHists(const char* name);

  virtual ~StPicoBesNetParticleHists();

  void init(TList *outList, const Char_t* title, Int_t analysisIdx, Int_t nCentralityBins);

  // -- Fill THnSparse
  void FillHnEvent(Double_t *aEvent);
  void FillHnTrack(Double_t *aTrack);

  // -- Event Statistics
  void InitializeEventStats();
  Bool_t FillEventStats(Int_t *aEventCuts);

  // -- Histogram sets for particle and anti-particle
  void AddHistSetCent(const Char_t *name, const Char_t *title);
  void FillHistSetCent(const Char_t *name, Int_t* np, Int_t cent);

  // -- Multiplicity statistics of accepted events (is that true?)
  void InitializeMultiplicityStats();
  void FillMultiplicityStats(Double_t *aMult, Int_t mode);
  
  // -- Event QA Histograms
  void InitializeEventQAHists();
  void FillEventQAHists(Double_t *aEvent, Int_t mode);
    
  // -- Track QA Histograms
  void InitializeTrackQAHists();
  void FillTrackQAHists(Double_t *aTrack, Int_t mode);

#if 0
  // -- Currently not implemented
  void InitializeRunByRunEventHists();
  void FillRunByRunEventHists(Double_t *aEvent, Int_t mode, Int_t runIdx);

  void InitializeRunByRunTrackHists();
  void FillRunByRunTrackHists(Double_t *aTrack, Int_t isBadRun, Int_t mode, Int_t runIdx);
#endif

  /*
   * ---------------------------------------------------------------------------------
   *                    Static Const Members - public
   * ---------------------------------------------------------------------------------
   */
 
  enum particleCharge {kPOS, kNEG, kNET, kParticleCharge};
  enum analysisType {kNetCharge, kNetProton, kNetKaon, kMaxAnalysisType};

  // -- Energy / Analysis / RefMult names
  static const Int_t   kNEnergies;
  static const Char_t* kNnergies[];
  static const Char_t* kExactEnergies[];
  
  static const Int_t   kNNames;
  static const Char_t* kName[];
  static const Char_t* kNameShort[];
  
  static const Int_t   kNRefMult;
  static const Char_t* kNameRefMult[];
  static const Char_t* kNameRefMultShort[];
  
  static const Int_t   kNParticles; 
  static const Char_t* kParticleName[kMaxAnalysisType][2];
  static const Char_t* kParticleTitle[kMaxAnalysisType][2];

  // -- Histogram sets : nSets / names / titles

  static const Int_t   kNQASets;
  static const Char_t* kQANames[];
  static const Char_t* kQATitles[];
  
  static const Int_t   kNQARunSets;
  static const Char_t* kQARunNames[];
  static const Char_t* kQARunTitles[]; 
  
  static const Int_t   kNMultSets;
  static const Char_t* kMultNames[];
  static const Char_t* kMultTitles[];

  // -- Histogram/THnSparse binnings
  static const Int_t    kNHnEvent;
  static const Int_t    kBinHnEvent[];
  static const Double_t kMinHnEvent[];
  static const Double_t kMaxHnEvent[];
  
  static const Int_t    kNHnTrack;  
  static const Int_t    kBinHnTrack[];
  static const Double_t kMinHnTrack[];
  static const Double_t kMaxHnTrack[];
 
  static const Int_t    kNCentralityBins;
  static const Char_t*  kCentralityTitles[];
  static const Char_t*  kCentralityMaxNames[];


 private:
  Int_t            mAnalysisIdx;
  Int_t            mNCentralityBins;

  TList*           mEventList;
  TList*           mAnaList;
  TList*           mQAList;
  TList*           mQARunList;

  THnSparseD*      mHnEvent;               //  THnSparseD : event
  THnSparseD*      mHnTrack;               //  THnSparseD : probe particles
  Int_t            mNRuns;

  // StPicoPrescales* mPrescales;
  
  ClassDef(StPicoBesNetParticleHists, 1)
};
#endif
