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
#include "RooFit.h"
#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "RooTauLeptonInvisible.h" 
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"
#include "RooRealConstant.h"

using namespace std;
 
//ClassImp(RooTauLeptonInvisible);

////////////////////////////////////////////////////////////////////////////////
 
RooTauLeptonInvisible::RooTauLeptonInvisible(const char *name, const char *title,
					     RooAbsReal& _x) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  mb("mb","Observable",this,(RooRealVar&)RooRealConstant::value(-1))
 {
   f.Init_PseudoRestFrameFuntions(0);
 }

RooTauLeptonInvisible::RooTauLeptonInvisible(const char *name, const char *title,
					     RooAbsReal& _x, RooAbsReal & _mb) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  mb("mb","Observable",this,_mb)
 {
   f.Init_PseudoRestFrameFuntions(1);
 }


////////////////////////////////////////////////////////////////////////////////
 
RooTauLeptonInvisible::RooTauLeptonInvisible(const RooTauLeptonInvisible& other, const char* name) :
  RooAbsPdf(other,name), x("x",this,other.x), mb("mb",this,other.mb), f(other.f)
 {
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 
 Double_t RooTauLeptonInvisible::evaluate() const
 {
   Double_t xval = x;
   Double_t mass_boson = mb;
   Double_t fval = 0;
   if(f._type == 0 ) fval = f.Unnormalized_PRF_lepton_nu_nu_PDF(xval);
   if(f._type == 1 ) fval =  f.Unnormalized_PRF2Body_PDF(xval,mass_boson);
   return fval;
 }
 
