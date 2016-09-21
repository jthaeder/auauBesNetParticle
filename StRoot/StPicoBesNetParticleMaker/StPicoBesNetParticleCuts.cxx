#include <limits>

#include "StBESCuts.h"

#include "StHFPair.h"
#include "StHFTriplet.h"

ClassImp(StBESCuts)

// _________________________________________________________
StBESCuts::StBESCuts() : StPicoCutsBase("HFCutsBase"), 
//  mSecondaryPairDcaDaughtersMax(std::numeric_limits<float>::max()), 
{
  // -- default constructor
}

// _________________________________________________________
StBESCuts::StBESCuts(const Char_t *name) : StPicoCutsBase(name), 
					 //					 mSecondaryPairDcaDaughtersMax(std::numeric_limits<float>::max()), 
{
  // -- constructor
}

// _________________________________________________________
StBESCuts::~StBESCuts() { 
  // destructor
  
}

// =======================================================================

#if 0
// _________________________________________________________
bool StBESCuts::isClosePair(StHFPair const & pair) const {
  // -- check for a pair which is close in dca w/o mass constraint,
  //    using secondary vertex cuts
  return ( std::cos(pair.pointingAngle()) > mSecondaryPairCosThetaMin &&
	   pair.decayLength() > mSecondaryPairDecayLengthMin && pair.decayLength() < mSecondaryPairDecayLengthMax &&
	   pair.dcaDaughters() < mSecondaryPairDcaDaughtersMax);
}
#endif
