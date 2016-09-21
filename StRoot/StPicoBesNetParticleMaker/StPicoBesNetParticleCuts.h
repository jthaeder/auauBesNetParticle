#ifndef StBESCUTS_H
#define StBESCUTS_H

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

#include "StPicoCutsBase/StPicoCutsBase.h"

class StBESCuts : public StPicoCutsBase
{
 public:
  
  StBESCuts();
  StBESCuts(const Char_t *name);
  ~StBESCuts();
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  virtual void init() { initBase(); }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  //
  bool isClosePair(StHFPair const & pair) const;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- SETTER for CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  
  //  void setCutSecondaryPair(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
  //float cosThetaMin, float massMin, float massMax); 

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- GETTER for single CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  //  const float&    cutSecondaryPairDcaDaughtersMax()       const;

 private:
  
  StBESCuts(StBESCuts const &);       
  StBESCuts& operator=(StBESCuts const &); 

  // ------------------------------------------
  // -- Pair cuts for secondary pair
  // ------------------------------------------
  //  float mSecondaryPairDcaDaughtersMax;

  ClassDef(StBESCuts,1)
};

//inline const float&    StBESCuts::cutSecondaryPairDcaDaughtersMax()       const { return mSecondaryPairDcaDaughtersMax; }

#endif
