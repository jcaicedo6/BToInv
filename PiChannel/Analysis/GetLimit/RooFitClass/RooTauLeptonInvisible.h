///////////////////////////////////////////////////////////////////////////////
//
//  RooTauLeptonInvisible class implementation
//  Roofit class for 2 or 3 body decays in the pseudo tau rest frame
//
//  Distribution of x = 2E/m_tau for
//        type  0:  tau->e- + nu_tau + anti-nu_e
//        type  1:  tau->e- + alpha 
//
//  Contributors: Eduard De La Cruz Burelo, CINVESTAV IPN Mexico, July 2019
//
////////////////////////////////////////////////////////////////////////////////


#include "TPseudoRestFrame.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"

class RooTauLeptonInvisible : public RooAbsPdf {

 public:

  RooTauLeptonInvisible(const char *name, const char *title, RooAbsReal &_x);
  RooTauLeptonInvisible(const char *name, const char *title, RooAbsReal &_x, RooAbsReal &_mb);
  RooTauLeptonInvisible(const  RooTauLeptonInvisible& other, const char *name=0);
  virtual TObject* clone(const char *newname) const { return new RooTauLeptonInvisible(*this,newname);}
  inline virtual ~RooTauLeptonInvisible(){}
  
 protected:
  
  TPseudoRestFrame f;
  RooRealProxy mb;
  RooRealProxy x;
  Double_t evaluate() const;
  
 private:
  ClassDef(RooTauLeptonInvisible,1) 
};
