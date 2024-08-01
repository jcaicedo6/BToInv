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
#define TPseudoRestFrame_cxx
#ifndef TPseudoRestFrame_h
#include "TPseudoRestFrame.h"
#endif

#include "UVHisto_m_0_1M.C"
#include "XY_1M_SM.C"
 
TPseudoRestFrame::TPseudoRestFrame()
{
}

TPseudoRestFrame::TPseudoRestFrame(const TPseudoRestFrame& other) :
  _type(other._type),
  _init_xy(other._init_xy),
  _init_uv(other._init_uv),
  _hXY(other._hXY),
  _hUV(other._hUV),
  _Ecm(other._Ecm), 
  _E_tau(other._E_tau),
  _m_tau(other._m_tau),
  _m_e(other._m_e), 
  _beta(other._beta),
  _gammav(other._gammav),
  _first_bin_u(other._first_bin_u),
  _last_bin_u(other._last_bin_u),
  _first_bin_v(other._first_bin_v),
  _last_bin_v(other._last_bin_v),
  _umin(other._umin),
  _umax(other._umax),
  _vmin(other._vmin),
  _vmax(other._vmax),
  _first_bin_x(other._first_bin_x),
  _last_bin_x(other._last_bin_x),
  _first_bin_y(other._first_bin_y),
  _last_bin_y(other._last_bin_y),
  _xmin(other._xmin),
  _xmax(other._xmax),
  _ymin(other._ymin),
  _ymax(other._ymax)
{
}

void TPseudoRestFrame::build_Join_Pdf_SM()
{
  // cout<<"*********************************"<<endl<<endl;
  // cout<<"  Building SM PDF ... done "<<endl<<endl;
  // cout<<"*********************************"<<endl;

  _hXY = new TH2D("_hXY","X vs Y",500,-0.1,1.1,500,-2.0,20.0); // Histogram to build join pdf for tau->e + anti_nu_e nu_tau
  XY_1M_SM(_hXY);
  
  Double_t sc = _hXY->Integral("width");
  _hXY->Scale(1.0/sc);
  _hXY->Smooth(1,"k5b");

  _first_bin_x = _hXY->GetXaxis()->GetFirst();
  _last_bin_x = _hXY->GetXaxis()->GetLast();
  _first_bin_y = _hXY->GetYaxis()->GetFirst();
  _last_bin_y = _hXY->GetYaxis()->GetLast();

  _xmin = _hXY->GetXaxis()->GetBinLowEdge(_first_bin_x);
  _xmax = _hXY->GetXaxis()->GetBinUpEdge(_last_bin_x);
  _ymin = _hXY->GetYaxis()->GetBinLowEdge(_first_bin_y);
  _ymax = _hXY->GetYaxis()->GetBinUpEdge(_last_bin_y);

  _init_xy=true;
}


void TPseudoRestFrame::build_Join_Pdf_UV()
{
  // cout<<"*********************************"<<endl<<endl;
  // cout<<"  Building Non-SM PDF ... done "<<endl<<endl;
  // cout<<"*********************************"<<endl;

  _hUV = new TH2D("_hUV","U vs V",1000,-1.0,1.0,1000,-2.0,2.0); // Histogram to build join pdf for tau->e + boson
  UVHisto_m_0_1M(_hUV);  
  Double_t sc = _hUV->Integral("width");
  _hUV->Scale(1.0/sc);
  _hUV->Smooth(1,"k5b");
    
  _first_bin_u = _hUV->GetXaxis()->GetFirst();
  _last_bin_u = _hUV->GetXaxis()->GetLast();
  _first_bin_v = _hUV->GetYaxis()->GetFirst();
  _last_bin_v = _hUV->GetYaxis()->GetLast();

  _umin = _hUV->GetXaxis()->GetBinLowEdge(_first_bin_u);
  _umax = _hUV->GetXaxis()->GetBinUpEdge(_last_bin_u);
  _vmin = _hUV->GetYaxis()->GetBinLowEdge(_first_bin_v);
  _vmax = _hUV->GetYaxis()->GetBinUpEdge(_last_bin_v);

  _init_uv=true;
}
  

void TPseudoRestFrame::Init_PseudoRestFrameFuntions(Int_t type)
{
  // cout<<"*********************************"<<endl<<endl;
  // cout<<"  Uploading functions ... done "<<endl<<endl;
  // cout<<"*********************************"<<endl;

  _Ecm = 10.5794;
  _E_tau = _Ecm/2.0;
  _m_tau = 1.776;
  _m_e = 0.00051;
  _beta = TMath::Sqrt(_E_tau*_E_tau - _m_tau*_m_tau)/_E_tau;
  _gammav = 1.0/TMath::Sqrt(1.0 - _beta*_beta);

  _type = type;
  if(_type==0) {build_Join_Pdf_SM();}
  if(_type==1) {build_Join_Pdf_UV();}

  
}

Double_t TPseudoRestFrame::xrf(Double_t M,Double_t m,Double_t m1) const
{
  Double_t x0 = (M*M - m*m + m1*m1)/(M*M);
  return x0;
}

Double_t TPseudoRestFrame::Pdf_ZW(Double_t w, Double_t z,Double_t x0) const
{

  if(isinf(z) || isinf(w)) return 0;
  
  Double_t A = _beta*_gammav*_gammav*x0;
  Double_t factor = A*(1.0 - _beta*w);
  Double_t v = z/factor;
  Double_t Jacobian = fabs(-1.0/factor);
  Double_t hv = 0;
  
  if((w>_umin && w<_umax) && (v>_vmin && v<_vmax))  hv = _hUV->Interpolate(w,v); 
  else hv = 0;
  
  Double_t fxy = Jacobian*hv;

  return fxy;
  
}

/*Double_t TPseudoRestFrame::Aproximate_Integral(Double_t z,Double_t x0) const
{
 //Trapezoidal rule
 Int_t n = 10000;
  Double_t a = _umin;
  Double_t b = _umax;
  Double_t delta = (b - a);
  Double_t h = delta/(n-1); 
  Double_t fval = 0.5*( Pdf_ZW(a,z,x0)+ Pdf_ZW(b,z,x0));
  //Double_t fval = 0;
  Double_t Integral = 0;

  for(Int_t k = 2; k<n; k++)
    {
      fval += Pdf_ZW(a+h*(k-1),z,x0);
      

      //Integral+=B*fval; 
      //w_tmp+=B;
    }
     
  fval*=h;
    
  Double_t IntegralX = fval;
  return IntegralX;
  }*/

/*Double_t TPseudoRestFrame::Aproximate_Integral(Double_t z,Double_t x0) const
{
  //Simpson's Rule
  Int_t n = 5000;
  Double_t a = _umin;
  Double_t b = _umax;
  Double_t delta = (b - a);
  Double_t h = delta/(n-1); 
  Double_t fval = (1/3)*( Pdf_ZW(a,z,x0)+ Pdf_ZW(b,z,x0));
  //Double_t fval = 0;
  Double_t Integral = 0;

  for(Int_t k = 1; k<n; k++)
    {
      if( (k%2)==0 ) fval += 2.0*Pdf_ZW(a+h*(k-1),z,x0);
      if( (k%2)!=0 ) fval += 4.0*Pdf_ZW(a+h*(k-1),z,x0);
      

      //Integral+=B*fval; 
      //w_tmp+=B;
    }
     
  fval*=h;
    
  Double_t IntegralX = fval;
  return IntegralX;
}*/

/*Double_t TPseudoRestFrame::Aproximate_Integral(Double_t z,Double_t x0) const
{
  //Simpson's Rule
  Int_t n = 5000;
  Double_t a = _umin;
  Double_t b = _umax;
  Double_t delta = (b - a);
  Double_t h = delta/(n-1); 
  Double_t fval = 0;//(1/3)*( Pdf_ZW(a,z,x0)+ Pdf_ZW(b,z,x0));
  //Double_t fval = 0;
  Double_t GaussSum  = 0.0;
  Double_t Integral = 0;

  for(Int_t k = 0; k<n; k++)
    {
      //GaussSum += Pdf_ZW(h*(k+0.5-0.5/TMath::Sqrt(3.)),z,x0) +  Pdf_ZW(h*(k+0.5+0.5/TMath::Sqrt(3.)),z,x0);
      GaussSum_5th += Pdf_ZW(h*(k+0.5-0.5/TMath::Sqrt(3.)),z,x0) +  Pdf_ZW(h*(k+0.5+0.5/TMath::Sqrt(3.)),z,x0);
      
      

      //Integral+=B*fval; 
      //w_tmp+=B;
    }
  fval=GaussSum;   
  fval*=h;
    
  Double_t IntegralX = fval;
  return IntegralX;
  }*/

/*Double_t TPseudoRestFrame::Aproximate_Integral(Double_t z,Double_t x0) const
{
  Double_t integral = Pdf_ZW()->Integral(_umin,_umax);
  }*/
 

Double_t TPseudoRestFrame::Aproximate_Integral(Double_t z,Double_t x0) const
{
  Int_t n = 5000;
  Double_t delta = (_umax - _umin);
  Double_t B = delta/(3.0*n); 
  Double_t w_tmp = _umin;
  Double_t fval = 0;
  Double_t Integral = 0;

  for(Int_t i = 0; i<=n;i++)
    {
      fval = Pdf_ZW(w_tmp,z,x0);
      if( i==0 || i==n) fval = 1.0*fval; 
      if( (i%2)==0 ) fval = 2.0*fval;
      if( (i%2)!=0 ) fval = 4.0*fval;

      Integral+=B*fval; 
      w_tmp+=B;
    }
    
  Double_t IntegralX = Integral;
  return IntegralX;
  }

Double_t TPseudoRestFrame::Unnormalized_PRF2Body_PDF(Double_t x,Double_t boson_mass) const
{
  // The Integral function use the variable x = E/2m_tau
  Double_t x0 = xrf(_m_tau,boson_mass,_m_e);
  Double_t fval =  Aproximate_Integral(x - x0,x0); 
  return fval;
}

//**********************************************
//**********************************************
// SM functions
//
//**********************************************

Double_t TPseudoRestFrame::GetIntegralShiftedY(Double_t z) const
{
  //Let's compute the integral for f(x,y)~histogram(x,y)
  //  for a fix y value
  //  g(z) = Int f(x,z - y) dx ~ Sum f(x,z - y) d_x
  // 

  Int_t nbin_x = _last_bin_x;
  TH1D *h_px = new TH1D("h_px","projection in X",nbin_x,_xmin,_xmax);
  Double_t fxy = 0;
  for(Int_t i = _first_bin_x; i<=_last_bin_x;i++)
    {
      Double_t x_eval = h_px->GetBinCenter(i);
      Double_t y_eval = z - x_eval;
      Int_t y_bin = _hXY->GetYaxis()->FindBin(y_eval);
      fxy = _hXY->GetBinContent(i,y_bin);
      h_px->SetBinContent(i,fxy);
    }

  Double_t Norm = 1.0;//_hXY->Integral();  
  Double_t IntegralX = Norm*(h_px->Integral());
    
  delete h_px;
  return z>0 ? IntegralX : 0;
}

Double_t TPseudoRestFrame::Unnormalized_PRF_lepton_nu_nu_PDF(Double_t x) const
{
  // The Integral function use the variable x = E/2m_tau
  return GetIntegralShiftedY(x);
}




