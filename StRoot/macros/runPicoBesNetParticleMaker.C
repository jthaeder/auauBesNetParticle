
/* **************************************************
 *   Run LambdaC Maker in different modes
 * --------------------------------------------------
 * run as :
 *  root -l -b -q StRoot/macros/loadSharedBesLibraries.C StRoot/macros/runPicoBesNetParticleMaker.C++
 *   or
 *  root -l -b -q StRoot/macros/runPicoBesNetParticleMaker.C++
 *
 * -------------------------------------------------- 
 *  - Different modes to use the  class
 *    - StPicoHFMaker::kAnalyze - don't write candidate trees, just fill histograms
 *        inputFile : fileList of PicoDst files or single picoDst file
 *        outputFile: baseName for outfile 
 *    - StPicoHFMaker::kWrite   - write candidate trees
 *        inputFile : path to single picoDist file
 *        outputFile: baseName for outfile 
 *    - StPicoHFMaker::kRead    - read candidate trees and fill histograms
 *        inputFile : fileList of PicoDst files
 *        outputFile: baseName for outfile 
 *
 * --------------------------------------------------
 *  Authors:  Jochen Thaeder  (jmthader@lbl.gov)
 *
 * **************************************************
 */

#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"

#include "StMaker.h"
#include "StChain.h"

#include "StPicoDstMaker/StPicoDstMaker.h"

#include "StPicoBesNetParticleMaker/StPicoBesNetParticleCuts.h"
#include "StPicoBesNetParticleMaker/StPicoBesNetParticleMaker.h"

#include "macros/loadSharedBesLibraries.C"

#include <iostream>
#include <ctime>
#include <cstdio>

using namespace std;

#else
class StChain;
#endif

StChain *chain;

void runPicoBesNetParticleMaker(const Char_t *inputFile="test.list", 
				const Char_t *outputFile="outputBaseName", 
				const unsigned int energyIdx = 2,    /* 14.5 GeV */
				const unsigned int analysisIdx = 0,  /* NetCharge */
				const unsigned int qaMode = 1){      /* QA mode on */
  // -- Check STAR Library. Please set SL_version to the original star library used in the production 
  //    from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
  string SL_version = "SL16d";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFLambdaCMaker.C. Exiting..."<<endl;
      exit(1);
  }
  
  Int_t nEvents = 100000000;

#ifdef __CINT__
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

  chain = new StChain();
  
  cout << "Energy Mode   " << energyIdx   << " -> " << StPicoBesNetParticleHists::kExactEnergies[energyIdx] << endl;
  cout << "Analysis Mode " << analysisIdx << " -> " << StPicoBesNetParticleHists::kName[analysisIdx] << endl;
  cout << "QA Mode       " << qaMode << endl; 

  // -- Check if inputFile is either list or single file
  TString sInputFile(inputFile);
  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }
   
  // --Create classes
  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, sInputFile, "picoDstMaker");

  StPicoBesNetParticleMaker* picoBesNetParticleMaker = StPicoBesNetParticleMaker("picoBesNetParticleMaker",picoDstMaker, outputFile);

  // ---------------------------------------------------
  // -- Set analysis type and energy
  // ---------------------------------------------------
  picoBesNetParticleMaker->setEnergyIdx(energyIdx);
  picoBesNetParticleMaker->setAnalysisIdx(analysisIdx);
  picoBesNetParticleMaker->setQaMode(qaMode);

  // ---------------------------------------------------
  // -- Set cuts for Net-Particle analysis
  // ---------------------------------------------------
  StPicoBesNetParticleCuts* besCuts = new StPicoBesNetParticleCuts(Form("%s_Cuts", StPicoBesNetParticleHists::kName[analysisIdx]));
  picoBesNetParticleMaker->setBesCuts(besCuts);

  // >>>>>>>>>>>>>>>>>>>>>>------------------------------<<<<<<<<<<<<<<<<<<<<<<
  // -- Event Cuts
  // >>>>>>>>>>>>>>>>>>>>>>------------------------------<<<<<<<<<<<<<<<<<<<<<<

  // -- Trigger
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  if (energyIdx == 0 ) { 
    besCuts->addTriggerId(290001);
    besCuts->addTriggerId(290004);
  }
  else if (energyIdx == 1 ) { 
    besCuts->addTriggerId(310004);
    besCuts->addTriggerId(310014);
  }
  else if (energyIdx == 2 ) {
    besCuts->addTriggerId(440015);
    besCuts->addTriggerId(440016);
  }
  else if (energyIdx == 3 ) { 
    besCuts->addTriggerId(340001);
    besCuts->addTriggerId(340011);
    besCuts->addTriggerId(340021);
  }
  else if (energyIdx == 4 ) {
    besCuts->addTriggerId(360001);
  }
  else if (energyIdx == 5 ) { 
    besCuts->addTriggerId(280001);
  }
  else if (energyIdx == 6 ) {
    besCuts->addTriggerId(270001);
    besCuts->addTriggerId(270011);
    besCuts->addTriggerId(270021);
  }
  else if (energyIdx == 7 ) { 
    besCuts->addTriggerId(260001);
    besCuts->addTriggerId(260011);
    besCuts->addTriggerId(260021);
    besCuts->addTriggerId(260031);
  }

  // -- Vertex cuts
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  if (energyIdx == 0)
    besCuts->setCutVzMax(50.);
  else 
    besCuts->setCutVzMax(30.);

  if (energyIdx > 4)
    besCuts->setCutVzVpdVzMax(3.);

  if (energyIdx == 2)
    besCuts->SetVyShift(-0.89);

  besCuts->SetVrMax(1.);
  
  // -- Centrality cuts
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  besCuts->SetNCentralityBinsMax(8);

  // -- Pile up removal cuts
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  besCuts->SetNTOFMatchMin(2);
  besCuts->SetCutRefMultVsNTOFMatch(0.78, 10.2);

  // >>>>>>>>>>>>>>>>>>>>>>------------------------------<<<<<<<<<<<<<<<<<<<<<<
  // -- Track Cuts
  // >>>>>>>>>>>>>>>>>>>>>>------------------------------<<<<<<<<<<<<<<<<<<<<<<

  besCuts->setCutNHitsFitnHitsMax(0.52);
  besCuts->setCutNHitsFitMin(20); 
  besCuts->setCutNHitsDedxMin(5); 

  besCuts->setCutPrimaryDCAtoVtxMax(1.);

  // -- Net-Charge
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  if (analysisIdx == 0) {
    besCuts->setCutPtRange(0.2, 2., StPicoCutsBase::kPion);
    besCuts->setCutPtMidPoint(1.);
    
    besCuts->setCutEtaRange(-0.5, 0.5., StPicoCutsBase::kPion);

    // -- spallation protona
    besCuts->setCutPtRange(0.2, 0.4, StPicoCutsBase::kProton);
    besCuts->setCutTPCNSigma(2, StPicoCutsBase::kProton);
  }
  
  // -- Net-Proton
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  else if (analysisIdx == 1) {
    besCuts->setCutPtRange(0.4, 2., StPicoCutsBase::kProton);
    besCuts->setCutPtMidPoint(0.8.);
    
    besCuts->setCutYRange(-0.5, 0.5., StPicoCutsBase::kProton);

    besCuts->setCutTPCNSigma(3, StPicoCutsBase::kProton);
    besCuts->setCutTOFmSquaredRange(0.6, 1.2, StPicoCutsBase::kProtom);
  }
  // -- Net-Kaon
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  else if (analysisIdx == 2) {
    besCuts->setCutPtRange(0.2, 1.6, StPicoCutsBase::kKaon);
    besCuts->setCutPtMidPoint(0.4);
    
    besCuts->setCutYRange(-0.5, 0.5., StPicoCutsBase::kKaon);

    besCuts->setCutTPCNSigma(2, StPicoCutsBase::kKaon);
    besCuts->setCutTOFmSquaredRange(0.15, 0.4, StPicoCutsBase::kKaon);
  }
    
  // ========================================================================================

  std::clock_t start = std::clock(); // getting starting time 
  chain->Init();
  cout << "chain->Init();" << endl;
  int total = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++) {
    if(i%1000==0)
      cout << "Working on eventNumber " << i << endl;
    
    chain->Clear();
    int iret = chain->Make(i);
    
    if (iret) { cout << "Bad return code!" << iret << endl; break;}
    
    total++;
  }
  
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  double duration = (double) (std::clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
  cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;
  
  delete chain;

  // -- clean up if in read mode
  if (makerMode == StPicoHFMaker::kRead)
    gSystem->Exec(Form("rm -f %s", sInputFile.Data()));
}

