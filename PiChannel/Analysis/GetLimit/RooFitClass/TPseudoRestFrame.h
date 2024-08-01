///////////////////////////////////////////////////////////////////////////////
//
//  TPseudoRestFrame class implementation
//   Distribution of x = 2E/m_tau for 2 or 3 body decays in the
//   pseudo tau rest frame:
//        type  0:  tau->e- + nu_tau + anti-nu_e
//        type  1:  tau->e- + alpha 
//
//  Contributors: Eduard De La Cruz Burelo, CINVESTAV IPN Mexico, July 2019
//
////////////////////////////////////////////////////////////////////////////////
#ifndef TPseudoRestFrame_h
#define TPseudoRestFrame_h

#include <iostream>
#include <cmath>
#include <TROOT.h>
#include <TMath.h>
#include <TH2D.h>

using namespace std;

class TPseudoRestFrame {

 public:
  
  TPseudoRestFrame();
  TPseudoRestFrame(const TPseudoRestFrame& other);
  inline virtual ~TPseudoRestFrame(){};
  
 protected:
  Bool_t _init_xy=false;
  Bool_t _init_uv=false;
  Double_t _Ecm; /// Energy on the CM frame
  Double_t _E_tau; // Tau energy
  Double_t _m_tau; // Mass of the tau
  Double_t _m_e; // Mass of the electron
  Double_t _beta; // boost magnitude
  Double_t _gammav;
  
  Int_t _first_bin_u;
  Int_t _last_bin_u;
  Int_t _first_bin_v;
  Int_t _last_bin_v;
  
  Double_t _umin;
  Double_t _umax;
  Double_t _vmin;
  Double_t _vmax;
  
  Int_t _first_bin_x;
  Int_t _last_bin_x;
  Int_t _first_bin_y;
  Int_t _last_bin_y;

  Double_t _xmin;
  Double_t _xmax;
  Double_t _ymin;
  Double_t _ymax;

  TH2D *_hUV;
  TH2D *_hXY;

 public: 
  Int_t _type;
  void build_Join_Pdf_SM();
  void build_Join_Pdf_UV();
  void Init_PseudoRestFrameFuntions(Int_t type);
  Double_t xrf(Double_t M,Double_t m,Double_t m1) const;
  Double_t Pdf_ZW(Double_t w, Double_t z, Double_t x0) const;
  Double_t Aproximate_Integral(Double_t z,Double_t x0) const;
  Double_t Unnormalized_PRF2Body_PDF(Double_t x,Double_t boson_mass) const;
  Double_t GetIntegralShiftedY(Double_t z) const;
  Double_t Unnormalized_PRF_lepton_nu_nu_PDF(Double_t x) const;
};

#endif
