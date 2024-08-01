#include <TROOT.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "DataGen3x1.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <string.h>

using namespace RooFit;

//#include "/home/belle2/johancol/tau3x1/RooFitClass/RooTauLeptonInvisible.h"
//Adding library: make changes accordingly
//#pragma cling load("/home/belle2/johancol/tau3x1/RooFitClass/RooTauLeptonInvisible.so")


//Global variables

Double_t CMS_E = 10.58;
Double_t tau_m = 1.777;
Double_t h_m = 0.13957;
Double_t e_m = 0.000511;

TChain *chSpdf;
TChain *chSMpdf;
TChain *chQQpdf;

TChain *chS ;
TChain *chSM;
TChain *chQQ;

TTree *treeSpdf;
TTree *treeSMpdf;
TTree *treeQQpdf;

TTree *treeS;
TTree *treeSM;
TTree *treeQQ;

DataGen3x1 *tSpdf;
DataGen3x1 *tSMpdf;
DataGen3x1 *tQQpdf;

DataGen3x1 *tS;
DataGen3x1 *tSM;
DataGen3x1 *tQQ;


//Functions
void SetData();
//Double_t getLimit(Double_t lum, Int_t type, Int_t pdf_type);
Double_t getLimit(Double_t lum, Int_t type);
Double_t getLimit2D(Double_t lum, Int_t pdf_type);
TVector3 calculateThrust(const std::vector<TVector3> momenta);
Double_t get_X_PseudoRestFrame(TLorentzVector piAll, TLorentzVector el);
void GetHistogram( DataGen3x1 *t, TH1D *hS, Int_t type, Int_t sample);
void GetHistogram2D( DataGen3x1 *t, TH2D *hS, Int_t sample);

//========================================
//  Histograms definition
//  Q : qqbar bkg
//  B : tau+tau- bkg
//  QB: Q+B
  
Double_t M_min2Low = -5.0;
Double_t M_min2Up = 28.0;

Double_t M_max2Low = -1.0;
Double_t M_max2Up = 27.0;

  
Double_t Xps_Low = 0.1;
Double_t Xps_Up = 1.2;

Double_t q2_Low = -3;
Double_t q2_Up = 25;

  
Double_t M_minLow = -5.0;
Double_t M_minUp = 25.0;

Double_t M_maxLow = -3.0;
Double_t M_maxUp = 25.0;


Int_t nbins = 200;
Int_t nbinsH = 15;
Int_t nbinsL = 15;


// =============================

// For Mmin^2 limit

// For data ------

TH1D *hS_DMmin2 = new TH1D("hS_DMmin2","hS_DMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hS_DstarGammaMmin2 = new TH1D("hS_DstarGammaMmin2","hS_DstarGammaMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hS_DstarPi0Mmin2 = new TH1D("hS_DstarPi0Mmin2","hS_DstarPi0Mmin2",nbins,M_min2Low,M_min2Up);

TH1D *hBkgDMmin2 = new TH1D("hBkgDMmin2","hBkgDMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hBkgDstarGammaMmin2 = new TH1D("hBkgDstarGammaMmin2","hBkgDstarGammaMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hBkgDstarPi0Mmin2 = new TH1D("hBkgDstarPi0Mmin2","hBkgDstarPi0Mmin2",nbins,M_min2Low,M_min2Up);

// For pdf   ------

TH1D *hSD_Mmin2pdf = new TH1D("hSD_Mmin2pdf","hSD_Mmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hSDstarGamma_Mmin2pdf = new TH1D("hSDstarGamma_Mmin2pdf","hSDstarGamma_Mmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hSDstarPi0_Mmin2pdf = new TH1D("hSDstarPi0_Mmin2pdf","hS_DstarPi0Mmin2pdf",nbins,M_min2Low,M_min2Up);

TH1D *hBkgD_Mmin2pdf = new TH1D("hBkgD_Mmin2pdf","hBkgD_Mmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hBkgDstarGamma_Mmin2pdf = new TH1D("hBkgDstarGamma_Mmin2pdf","hBkgDstarGamma_Mmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hBkgDstarPi0_Mmin2pdf = new TH1D("hBkgDstarPi0_Mmin2pdf","hBkgDstarPi0_Mmin2pdf",nbins,M_min2Low,M_min2Up);

TH1D *hB1_Mmin2pdf = new TH1D("hB1_Mmin2pdf","hB1_Mmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hB_Mmin2pdf = new TH1D("hB_Mmin2pdf","hB_Mmin2pdf",nbins,M_min2Low,M_min2Up);

// =============================

// For Mmax^2 limit

// For data ------

TH1D *hSD_Mmax2 = new TH1D("hS_DMmax2","hS_DMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hSDstarGamma_Mmax2 = new TH1D("hS_DstarGammaMmax2","hS_DstarGammaMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hSDstarPi0_Mmax2 = new TH1D("hS_DMTauMmax2","hS_DMTauMmax2",nbins,M_max2Low,M_max2Up);

TH1D *hBkgD_Mmax2 = new TH1D("hBkgD_Mmax2","hBkgD_Mmax2",nbins,M_max2Low,M_max2Up);
TH1D *hBkgDstarGammaMmax2 = new TH1D("hBkgDstarGammaMmax2","hBkgDstarGammaMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hBkgDstarPi0Mmax2 = new TH1D("hBkgDstarPi0Mmax2","hBkgDstarPi0Mmax2",nbins,M_max2Low,M_max2Up);

// For pdf ------

TH1D *hSD_Mmax2pdf = new TH1D("hSD_Mmax2pdf","hSD_Mmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hSDstarGamma_Mmax2pdf = new TH1D("hSDstarGamma_Mmax2pdf","hSDstarGamma_Mmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hSDstarPi0_Mmax2pdf = new TH1D("hSDstarPi0_Mmax2pdf","hSDstarPi0_Mmax2pdf",nbins,M_max2Low,M_max2Up);

TH1D *hBkgD_Mmax2pdf = new TH1D("hBkgD_Mmax2pdf","hBkgD_Mmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hBkgDstarGamma_Mmax2pdf = new TH1D("hBkgDstarGamma_Mmax2pdf","hBkgDstarGamma_Mmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hBkgDstarPi0_Mmax2pdf = new TH1D("hBkgDstarPi0_Mmax2pdf","hBkgDstarPi0_Mmax2pdf",nbins,M_max2Low,M_max2Up);

TH1D *hB1_Mmax2pdf = new TH1D("hB1_Mmax2pdf","hB1_Mmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hB_Mmax2pdf = new TH1D("hB_Mmax2pdf","hB_Mmax2pdf",nbins,M_max2Low,M_max2Up);


// =============================

// For Xps limit

// For data ------

TH1D *hSD_Xps = new TH1D("hS_DXps","hS_DXps",nbins,Xps_Low,Xps_Up);
TH1D *hS_DstarGammaXps= new TH1D("hS_DstarGammaXps","hS_DstarGammaXps",nbins,Xps_Low,Xps_Up);
TH1D *hS_DstarPi0Xps= new TH1D("hS_DstarPi0Xps","hS_DstarPi0Xps",nbins,Xps_Low,Xps_Up);

TH1D *hBkgDXps = new TH1D("hBkgDXps","hBkgDXps",nbins,Xps_Low,Xps_Up);
TH1D *hBkgDstarGammaXps= new TH1D("hBkgDstarGammaXps","hBkgDstarGammaXps",nbins,Xps_Low,Xps_Up);
TH1D *hBkgDstarPi0Xps = new TH1D("hBkgDstarPi0Xps","hBkgDstarPi0Xps",nbins,Xps_Low,Xps_Up);

// For pdf ------

TH1D *hSD_Xpspdf = new TH1D("hSD_Xpspdf","hSD_Xpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hSDstarGamma_Xpspdf = new TH1D("hSDstarGamma_Xpspdf","hSDstarGamma_Xpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hSDstarPi0_Xpspdf = new TH1D("hSDstarPi0_Xpspdf","hSDstarPi0_Xpspdf",nbins,Xps_Low,Xps_Up);

TH1D *hBkgD_Xpspdf = new TH1D("hBkgD_Xpspdf","hBkgD_Xpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hBkgDstarGamma_Xpspdf = new TH1D("hBkgDstarGamma_Xpspdf","hBkgDstarGamma_Xpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hBkgDstarPi0_Xpspdf = new TH1D("hBkgDstarPi0_Xpspdf","hBkgDstarPi0_Xpspdf",nbins,Xps_Low,Xps_Up);

TH1D *hB1_Xpspdf = new TH1D("hB1_Xpspdf","hB1_Xpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hB_Xpspdf = new TH1D("hB_Xpspdf","hB_Xpspdf",nbins,Xps_Low,Xps_Up);


// For q2 limit

// For data ------

TH1D *hSD_q2 = new TH1D("hS_Dq2 ","hS_Dq2 ",nbins,q2_Low,q2_Up);
TH1D *hS_DstarGammaq2= new TH1D("hS_DstarGammaq2","hS_DstarGammaq2",nbins,q2_Low,q2_Up);
TH1D *hS_DstarPi0q2= new TH1D("hS_DstarPi0q2","hS_DstarPi0q2",nbins,q2_Low,q2_Up);

TH1D *hBkgDq2 = new TH1D("hBkgDq2","hBkgDq2",nbins,q2_Low,q2_Up);
TH1D *hBkgDstarGammaq2= new TH1D("hBkgDstarGammaq2","hBkgDstarGammaq2",nbins,q2_Low,q2_Up);
TH1D *hBkgDstarPi0q2 = new TH1D("hBkgDstarPi0q2","hBkgDstarPi0q2",nbins,q2_Low,q2_Up);

// For pdf ------

TH1D *hSD_q2pdf = new TH1D("hSD_q2pdf","hSD_q2pdf",nbins,q2_Low,q2_Up);
TH1D *hSDstarGamma_q2pdf = new TH1D("hSDstarGamma_q2pdf","hSDstarGamma_q2pdf",nbins,q2_Low,q2_Up);
TH1D *hSDstarPi0_q2pdf = new TH1D("hSDstarPi0_q2pdf","hSDstarPi0_q2pdf",nbins,q2_Low,q2_Up);

TH1D *hBkgD_q2pdf = new TH1D("hBkgD_Xpspdf","hBkgD_Xpspdf",nbins,q2_Low,q2_Up);
TH1D *hBkgDstarGamma_q2pdf = new TH1D("hBkgDstarGamma_q2pdf","hBkgDstarGamma_q2pdf",nbins,q2_Low,q2_Up);
TH1D *hBkgDstarPi0_q2pdf = new TH1D("hBkgDstarPi0_q2pdf","hBkgDstarPi0_q2pdf",nbins,q2_Low,q2_Up);

TH1D *hB1_q2pdf = new TH1D("hB1_q2pdf","hB1_q2pdf",nbins,q2_Low,q2_Up);
TH1D *hB_q2pdf = new TH1D("hB_q2pdf","hB_q2pdf",nbins,q2_Low,q2_Up);


// Defining the branchings and number of events
Double_t Br_B_to_l_nu_D0 = 2*0.0235;
Double_t Br_B_to_l_nu_D0star = 2*0.0558;
Double_t Br_D0_to_pi_K = 0.03947;
//Double_t Br_B_to_l_nu_D0star = 2*0.0549;
Double_t Br_D0star_to_D0_gamma = 0.353;
Double_t Br_D0star_to_D0_pi0 = 0.647;
Double_t Br_pi0_to_2gamma = 1;
//Double_t Br_pi0_to_2gamma = 0.9882;

Double_t Br_tagD = Br_B_to_l_nu_D0*Br_D0_to_pi_K;
Double_t Br_tagDstarGamma = Br_B_to_l_nu_D0star*Br_D0star_to_D0_gamma*Br_D0_to_pi_K;
Double_t Br_tagDstarPi0 = Br_B_to_l_nu_D0star*Br_D0star_to_D0_pi0*Br_D0_to_pi_K*Br_pi0_to_2gamma;


Double_t Luminosity = 1.08e6; //Normalized to the SemiLepTag_Dstar_pi0 
Double_t CS_ee_BB = 0.565e6; //in fb
Double_t Ns_generated = 100000;//*(Br_B_to_l_nu_D0);
//Double_t NB_generated = 2*Luminosity*CS_ee_BB*Br_B_to_l_nu_D0*Br_B_to_tau_pinu;




void FindLimitCls_BToInvCombined()
{
  SetData();
  Double_t limit[4];
  //Double_t lum = 361.6; 
  //Double_t lum = 403.9; //(361.6 + 42.3 (off 4S)) 
  //Double_t lum = 10000; 
  //Double_t lum = 30000; 
  Double_t lum = 50000; 
  

  // Depending on the kind of pdf we can estimnate the UL in normalized POI, Number of
  // Events, or directly the branching ratio, the pdfs used are:

  //   PDF1 = Snu*mu*effrel*S_pdf + Snu*pi_taunu_pdf   <- without background
  //   PDF2 = S*S_pdf + Snu*pi_taunu_pdf + B*pi_nunu_pdf  <- Fitting number of events
  //   PDF3 = fac*mu*S_pdf + Snu*pi_taunu_pdf + B*pi_nunu_pdf   <- Branching ratio as POI

  // 1:Dtag, 2:Dstar_Gamma_tag  3:Dstar_Pi0_tag
  //Int_t type_tag =1;
  //Int_t type_tag =2;
  //Int_t type_tag =3;


  // 1:Mmin2, 2:Mmax2,  3:Xps,  4:q2
  //limit[0] = getLimit(lum,1);
  //limit[1] = getLimit(lum,2);
  limit[2] = getLimit(lum,3);
  //limit[3] = getLimit(lum,4);
  //limit[3] = getLimit2D(lum,pdf_type);

  cout<<" **************************** "<<endl<<endl;
  //cout<<" Limit : "<<limit[0]<<endl;
  //cout<<" Limit : "<<limit[1]<<endl;
  cout<<" Limit : "<<limit[2]<<endl;
  //cout<<" Limit : "<<limit[3]<<endl;
  //cout<<" Limit : "<<limit[4]<<endl;

  
  return;
}


void SetData()
{
          //-----------------------------------------------
    //                              PDF
    //----------------------------------------------
    
    TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_m00_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m00_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m05_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m05_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m07_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m07_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m10_Mmin2_pdf.root");
    TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_m10_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m12_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m12_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m14_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m14_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m16_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_m16_Mmax2_pdf.root");

    
    TFile *fSMPinuMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_pinu_Mmin2_pdf.root");
    TFile *fSMPinuMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_pinu_Mmax2_pdf.root");
    TFile *fSMTauMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_Mmin2_pdf.root");
    TFile *fSMTauMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_Mmax2_pdf.root");
    TFile *fBMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_Mmin2_pdf.root");
    TFile *fBMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_Mmax2_pdf.root");

    
    //----------------------------------------------------
    //                         DATA
    //----------------------------------------------------

    
    //TFile *fSMPinuMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_pinu_Mmin2.root");
    //TFile *fSMPinuMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_pinu_Mmax2.root");
    //TFile *fSMTauMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_tau_Mmin2.root");
    //TFile *fSMTauMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_tau_Mmax2.root");
    //TFile *fBMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hB_Mmin2.root");
    //TFile *fBMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hB_Mmax2.root");
    
    
        
     //============================================
  //   This is for 2D pdf (Mmin^2, Mmax^2)
  //============================================
    //-----------------------------------------------
    //                              Data
    //----------------------------------------------
    
  
    //TFile *fSMpinu_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_pinu_2D.root");
    //TFile *fSMtau_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_tau_2D.root");
    //TFile *fB_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hB_2D.root");

    
    //----------------------------------------------------
    //                         PDF
    //----------------------------------------------------

    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m00_pdf.root");
    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m05_pdf.root");
    TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_2D_m07_pdf.root");
    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m10_pdf.root");
    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m12_pdf.root");
    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m14_pdf.root");
    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m16_pdf.root");
    //TFile *fS_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_2D_m16_pdf.root");
  
    TFile *fSMpinu_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_pinu_2D_pdf.root");
    TFile *fSMtau_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_2D_pdf.root");
    TFile *fB_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_2D_pdf.root");

  return; 

}

Double_t getLimit(Double_t lum, Int_t type)
{

  Int_t NsD = 0;
  Int_t NsDstarGamma = 0;
  Int_t NsDstarPi0 = 0;
  //Int_t Npitaunu = 0;
  Int_t Nb = 0;

  RooRealVar *x;

  RooDataHist *dSD;
  RooDataHist *dSDstarGamma;
  RooDataHist *dSDstarPi0;
  RooDataHist *dB;
    
  
 
  RooDataHist *hist_dataB;
  //RooDataHist *hist_data;
    
    

  RooHistPdf *s1_pdf;
  RooHistPdf *s2_pdf;
  RooHistPdf *s3_pdf;
  RooHistPdf *b_pdf;
  
    //------------------------------------------------------------------------
    //                      Mmin and Mmax   PDF
    //------------------------------------------------------------------------
    //******************* D_tag
    //
    TFile *fSDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m00_Mmin2_pdf.root");
    TFile *fSDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m00_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m10_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m10_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m20_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m20_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m30_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m30_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m40_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m40_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m45_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m45_Mmax2_pdf.root");

    //**************** Dstar_gamma_tag
    TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Mmin2_pdf.root");
    TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m45_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m45_Mmax2_pdf.root");
    
    //**************** Dstar_pi0_tag
    TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Mmin2_pdf.root");
    TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m45_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m45_Mmax2_pdf.root");    
    
    // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_D_Mmin2_pdf.root");
    TFile *fBkgStagDMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_D_Mmax2_pdf.root");
    TFile *fBkgStagDstarGammaMmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Mmin2_pdf.root");
    TFile *fBkgStagDstarGammaMmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Mmax2_pdf.root");
    TFile *fBkgStagDstarPi0Mmin2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Mmin2_pdf.root");
    TFile *fBkgStagDstarPi0Mmax2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Mmax2_pdf.root");

    
    //----------------------------------------------------
    //                         DATA
    //----------------------------------------------------

    
    //TFile *fSDMPinuMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_DM_pinu_Mmin2.root");
    //TFile *fSDMPinuMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_DM_pinu_Mmax2.root");
    //TFile *fSDMTauMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_DM_tau_Mmin2.root");
    //TFile *fSDMTauMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_DM_tau_Mmax2.root");
    //TFile *fBkgStagDstarPi0min2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hBkgStag_Dstar_pi0_min2.root");
    //TFile *fBkgStagDstarPi0max2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hBkgStag_Dstar_pi0_max2.root");
    
    //------------------------------------------------------------------------
    //                      Xps  PDF
    //------------------------------------------------------------------------
    //******************* D_tag
    //
    //TFile *fSDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m00_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m10_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m20_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m30_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m40_Xps_pdf.root");
    TFile *fSDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m45_Xps_pdf.root");
    
    //******************* Dstar_gamma_tag
    //
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Xps_pdf.root");
    TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m45_Xps_pdf.root");

    //******************* Dstar_pi0_tag
    //
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Xps_pdf.root");
    TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m45_Xps_pdf.root");

     // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_D_Xps_pdf.root");
    TFile *fBkgStagDstarGammaXpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Xps_pdf.root");
    TFile *fBkgStagDstarPi0Xpspdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Xps_pdf.root");
    
    
    //------------------------------------------------------------------------
    //                      q2  PDF
    //------------------------------------------------------------------------
    //******************* D_tag
    //
    //TFile *fSDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m00_q2_pdf.root");
    //TFile *fSDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m10_q2_pdf.root");
    //TFile *fSDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m20_q2_pdf.root");
    //TFile *fSDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m30_q2_pdf.root");
    //TFile *fSDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m40_q2_pdf.root");
    TFile *fSDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_D_m45_q2_pdf.root");
    
    //******************* Dstar_gamma_tag
    //
    //TFile *fSDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_q2_pdf.root");
    //TFile *fSDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_q2_pdf.root");
    //TFile *fSDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_q2_pdf.root");
    //TFile *fSDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_q2_pdf.root");
    //TFile *fSDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_q2_pdf.root");
    TFile *fSDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m45_q2_pdf.root");

    //******************* Dstar_pi0_tag
    //
    //TFile *fSDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_q2_pdf.root");
    //TFile *fSDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_q2_pdf.root");
    //TFile *fSDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_q2_pdf.root");
    //TFile *fSDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_q2_pdf.root");
    //TFile *fSDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_q2_pdf.root");
    TFile *fSDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m45_q2_pdf.root");

     // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_D_q2_pdf.root");
    TFile *fBkgStagDstarGammaq2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_q2_pdf.root");
    TFile *fBkgStagDstarPi0q2pdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_q2_pdf.root");
    
    
    
  if(type==1)
    {
      // Histograms for signal and bkg D tag
      hSD_Mmin2pdf = (TH1D*)fSDMmin2pdf->Get("Mmin2");
      hBkgD_Mmin2pdf = (TH1D*)fBkgStagDMmin2pdf->Get("Mmin2");
      // Histograms for signal and bkg DstarGamma tag
      hSDstarGamma_Mmin2pdf = (TH1D*)fSDstarGammaMmin2pdf->Get("Mmin2");
      hBkgDstarGamma_Mmin2pdf = (TH1D*)fBkgStagDstarGammaMmin2pdf->Get("Mmin2");
      // Histograms for signal and bkg DstarPi0 tag
      hSDstarPi0_Mmin2pdf = (TH1D*)fSDstarPi0Mmin2pdf->Get("Mmin2");
      hBkgDstarPi0_Mmin2pdf = (TH1D*)fBkgStagDstarPi0Mmin2pdf->Get("Mmin2");

      // Adding bkg histograms
      hB1_Mmin2pdf->Add(hBkgD_Mmin2pdf,hBkgDstarGamma_Mmin2pdf,1,1);
      hB_Mmin2pdf->Add(hB1_Mmin2pdf,hBkgDstarPi0_Mmin2pdf,1,1);

      //Number of events
      NsD = hSD_Mmin2pdf->GetEntries();
      NsDstarGamma = hSDstarGamma_Mmin2pdf->GetEntries();
      NsDstarPi0 = hSDstarPi0_Mmin2pdf->GetEntries();

      Nb = hB_Mmin2pdf->GetEntries();
      //NbD = hBkgD_Mmin2pdf->GetEntries();
      //NbDstarGamma = hBkgDstarGamma_Mmin2pdf->GetEntries();
      //NbDstarPi0 = hBkgDstarPi0_Mmin2pdf->GetEntries();

      
      // For Mmin^2
      x = new RooRealVar("x", "x", M_min2Low, M_min2Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_Mmin2pdf));  // B->Pi+alpha Dtag
      //dBD = new RooDataHist("dB", "dB", *x, Import(*hBkgD_Mmin2pdf));  // BKG

      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_Mmin2pdf));  // B->Pi+alpha DstarGamma tag
      //dBDstarGamma = new RooDataHist("dB", "dB", *x, Import(*hBkgDstarGamma_Mmin2pdf));  // BKG

      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_Mmin2pdf));  // B->Pi+alpha  Dstar Pi0 tag
      //dBDstarPi0 = new RooDataHist("dB", "dB", *x, Import(*hBkgDstarPi0_Mmin2pdf));  // BKG

      dB = new RooDataHist("dB", "dB", *x, Import(*hB_Mmin2pdf));  // BKG
    }
  
  if(type==2)
    { 
      hSD_Mmax2pdf = (TH1D*)fSDMmax2pdf->Get("Mmax2");
      hBkgD_Mmax2pdf = (TH1D*)fBkgStagDMmax2pdf->Get("Mmax2");

      hSDstarGamma_Mmax2pdf = (TH1D*)fSDstarGammaMmax2pdf->Get("Mmax2");
      hBkgDstarGamma_Mmax2pdf = (TH1D*)fBkgStagDstarGammaMmax2pdf->Get("Mmax2");
  
      hSDstarPi0_Mmax2pdf = (TH1D*)fSDstarPi0Mmax2pdf->Get("Mmax2");
      hBkgDstarPi0_Mmax2pdf = (TH1D*)fBkgStagDstarPi0Mmax2pdf->Get("Mmax2");

      // Adding bkg histograms
      hB1_Mmax2pdf->Add(hBkgD_Mmax2pdf,hBkgDstarGamma_Mmax2pdf,1,1);
      hB_Mmax2pdf->Add(hB1_Mmax2pdf,hBkgDstarPi0_Mmax2pdf,1,1);

      //Number of events
      NsD = hSD_Mmax2pdf->GetEntries();
      NsDstarGamma = hSDstarGamma_Mmax2pdf->GetEntries();
      NsDstarPi0 = hSDstarPi0_Mmax2pdf->GetEntries();

      Nb = hB_Mmax2pdf->GetEntries();
      //NbD = hBkgD_Mmax2pdf->GetEntries();
      //NbDstarGamma = hBkgDstarGamma_Mmax2pdf->GetEntries();
      //NbDstarPi0 = hBkgDstarPi0_Mmax2pdf->GetEntries();

      
      // For Mmin^2
      x = new RooRealVar("x", "x", M_max2Low, M_max2Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_Mmax2pdf));  // B->Pi+alpha Dtag
      //dBD = new RooDataHist("dBD", "dBD", *x, Import(*hBkgD_Mmax2pdf));  // BKG

      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_Mmax2pdf));  // B->Pi+alpha DstarGamma tag
      //dBDstarGamma = new RooDataHist("dBDstarGamma", "dBDstarGamma", *x, Import(*hBkgDstarGamma_Mmax2pdf));  // BKG

      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_Mmax2pdf));  // B->Pi+alpha  Dstar Pi0 tag
      //dBDstarPi0 = new RooDataHist("dBDstarPi0", "dBDstarPi0", *x, Import(*hBkgDstarPi0_Mmax2pdf));  // BKG

      dB = new RooDataHist("dB", "dB", *x, Import(*hB_Mmax2pdf));  // BKG
    }

 if(type==3)
    {
      hSD_Xpspdf = (TH1D*)fSDXpspdf->Get("Xps");
      hBkgD_Xpspdf = (TH1D*)fBkgStagDXpspdf->Get("Xps");

      hSDstarGamma_Xpspdf = (TH1D*)fSDstarGammaXpspdf->Get("Xps");
      hBkgDstarGamma_Xpspdf = (TH1D*)fBkgStagDstarGammaXpspdf->Get("Xps");
  
      hSDstarPi0_Xpspdf = (TH1D*)fSDstarPi0Xpspdf->Get("Xps");
      hBkgDstarPi0_Xpspdf = (TH1D*)fBkgStagDstarPi0Xpspdf->Get("Xps");

      // Adding bkg histograms
      hB1_Xpspdf->Add(hBkgD_Xpspdf,hBkgDstarGamma_Xpspdf,1,1);
      hB_Xpspdf->Add(hB1_Xpspdf,hBkgDstarPi0_Xpspdf,1,1);

      //Number of events
      NsD = hSD_Xpspdf->GetEntries();
      NsDstarGamma = hSDstarGamma_Xpspdf->GetEntries();
      NsDstarPi0 = hSDstarPi0_Xpspdf->GetEntries();

      Nb = hB_Xpspdf->GetEntries();
      //NbD = hBkgD_Xpspdf->GetEntries();
      //NbDstarGamma = hBkgDstarGamma_Xpspdf->GetEntries();
      //NbDstarPi0 = hBkgDstarPi0_Xpspdf->GetEntries();

      
      // For Mmin^2
      x = new RooRealVar("x", "x", Xps_Low, Xps_Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_Xpspdf));  // B->Pi+alpha Dtag
      //dBD = new RooDataHist("dBD", "dBD", *x, Import(*hBkgD_Xpspdf));  // BKG

      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_Xpspdf));  // B->Pi+alpha DstarGamma tag
      //dBDstarGamma = new RooDataHist("dBDstarGamma", "dBDstarGamma", *x, Import(*hBkgDstarGamma_Xpspdf));  // BKG

      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_Xpspdf));  // B->Pi+alpha  Dstar Pi0 tag
      //dBDstarPi0 = new RooDataHist("dBDstarPi0", "dBDstarPi0", *x, Import(*hBkgDstarPi0_Mmax2pdf));  // BKG

      dB = new RooDataHist("dB", "dBD", *x, Import(*hB_Xpspdf));  // BKG
    }
    
    if(type==4)
    {
      hSD_q2pdf = (TH1D*)fSDq2pdf->Get("q2");
      hBkgD_q2pdf = (TH1D*)fBkgStagDq2pdf->Get("q2");

      hSDstarGamma_q2pdf = (TH1D*)fSDstarGammaq2pdf->Get("q2");
      hBkgDstarGamma_q2pdf = (TH1D*)fBkgStagDstarGammaq2pdf->Get("q2");
  
      hSDstarPi0_q2pdf = (TH1D*)fSDstarPi0q2pdf->Get("q2");
      hBkgDstarPi0_q2pdf = (TH1D*)fBkgStagDstarPi0q2pdf->Get("q2");

      // Adding bkg histograms
      hB1_q2pdf->Add(hBkgD_q2pdf,hBkgDstarGamma_q2pdf,1,1);
      hB_q2pdf->Add(hB1_q2pdf,hBkgDstarPi0_q2pdf,1,1);

      //Number of events
      NsD = hSD_q2pdf->GetEntries();
      NsDstarGamma = hSDstarGamma_q2pdf->GetEntries();
      NsDstarPi0 = hSDstarPi0_q2pdf->GetEntries();

      Nb = hB_q2pdf->GetEntries();
      //NbD = hBkgD_Xpspdf->GetEntries();
      //NbDstarGamma = hBkgDstarGamma_Xpspdf->GetEntries();
      //NbDstarPi0 = hBkgDstarPi0_Xpspdf->GetEntries();

      
      // For Mmin^2
      x = new RooRealVar("x", "x", q2_Low, q2_Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_q2pdf));  // B->Pi+alpha Dtag
      //dBD = new RooDataHist("dBD", "dBD", *x, Import(*hBkgD_Xpspdf));  // BKG

      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_q2pdf));  // B->Pi+alpha DstarGamma tag
      //dBDstarGamma = new RooDataHist("dBDstarGamma", "dBDstarGamma", *x, Import(*hBkgDstarGamma_Xpspdf));  // BKG

      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_q2pdf));  // B->Pi+alpha  Dstar Pi0 tag
      //dBDstarPi0 = new RooDataHist("dBDstarPi0", "dBDstarPi0", *x, Import(*hBkgDstarPi0_Mmax2pdf));  // BKG

      dB = new RooDataHist("dB", "dBD", *x, Import(*hB_q2pdf));  // BKG
    }
    

  Nb = lum*Nb/Luminosity;
  //NbD = lum*NbD/Luminosity;
  //NbDstarGamma = lum*NbDstarGamma/Luminosity;
  //NbDstarPi0 = lum*NbDstarPi0/Luminosity;
  //Ns = lum*Ns/Luminosity;
  //Nbsm = lum*Nbsm/Luminosity;
  
  s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",*x,*dSD,2);
  s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",*x,*dSDstarGamma,2);
  s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",*x,*dSDstarPi0,2);

  b_pdf = new RooHistPdf("b_pdf","bkg pdf",*x,*dB,2);



  
  TCanvas *c = new TCanvas("c","c",800,600);
  c->Divide(2, 2);

  
  if(type==1){

      //TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c->cd(1);
      // Customize histogram attributes
      hSD_Mmin2pdf->SetLineColor(kBlue);
      hSD_Mmin2pdf->SetFillColor(kBlue);
      hSD_Mmin2pdf->SetFillStyle(3001);  // Diagonal lines
      
      hSDstarGamma_Mmin2pdf->SetLineColor(kRed);
      hSDstarGamma_Mmin2pdf->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
      hSDstarGamma_Mmin2pdf->SetMarkerStyle(21);
      hSDstarGamma_Mmin2pdf->SetMarkerSize(0.7);
      
      hSDstarPi0_Mmin2pdf->SetLineColor(kGreen);

       // Draw histograms
      hSD_Mmin2pdf->Draw("HIST");
      hSDstarGamma_Mmin2pdf->Draw("HISTsame");
      hSDstarPi0_Mmin2pdf->Draw("HISTsame");
      hB_Mmin2pdf->Draw("HISTsame");
      //c1->Draw();
      c->Draw();
  }
    
  if(type==2){
      //TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c->cd(1);
      // Customize histogram attributes
      hSD_Mmax2pdf->SetLineColor(kBlue);
      hSD_Mmax2pdf->SetFillColor(kBlue);
      hSD_Mmax2pdf->SetFillStyle(3001);  // Diagonal lines
      
      hSDstarGamma_Mmin2pdf->SetLineColor(kRed);
      hSDstarGamma_Mmin2pdf->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
      hSDstarGamma_Mmin2pdf->SetMarkerStyle(21);
      hSDstarGamma_Mmin2pdf->SetMarkerSize(0.7);
      
      hSDstarPi0_Mmin2pdf->SetLineColor(kGreen);

       // Draw histograms
      hSD_Mmax2pdf->Draw("HIST");
      hSDstarGamma_Mmax2pdf->Draw("HISTsame");
      hSDstarPi0_Mmax2pdf->Draw("HISTsame");
      hB_Mmax2pdf->Draw("HISTsame");
      //c1->Draw();
      c->Draw();
  }
  
  if(type==3){
      //TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c->cd(1);
      // Customize histogram attributes
      hSD_Xpspdf->SetLineColor(kBlue);
      hSD_Xpspdf->SetFillColor(kBlue);
      hSD_Xpspdf->SetFillStyle(3001);  // Diagonal lines
      
      hSDstarGamma_Xpspdf->SetLineColor(kRed);
      hSDstarGamma_Xpspdf->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
      hSDstarGamma_Xpspdf->SetMarkerStyle(21);
      hSDstarGamma_Xpspdf->SetMarkerSize(0.7);
      
      hSDstarPi0_Xpspdf->SetLineColor(kGreen);

       // Draw histograms
      hSD_Xpspdf->Draw("HIST");
      hSDstarGamma_Xpspdf->Draw("HISTsame");
      hSDstarPi0_Xpspdf->Draw("HISTsame");
      hB_Xpspdf->Draw("HISTsame");
      //c1->Draw();
      c->Draw();
  }
  
  if(type==4){
      //TCanvas *c1 = new TCanvas("c1","c1",800,600);
      c->cd(1);
      // Customize histogram attributes
      hSD_q2pdf->SetLineColor(kBlue);
      hSD_q2pdf->SetFillColor(kBlue);
      hSD_q2pdf->SetFillStyle(3001);  // Diagonal lines
      
      hSDstarGamma_q2pdf->SetLineColor(kRed);
      hSDstarGamma_q2pdf->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
      hSDstarGamma_q2pdf->SetMarkerStyle(21);
      hSDstarGamma_q2pdf->SetMarkerSize(0.7);
      
      hSDstarPi0_q2pdf->SetLineColor(kGreen);

       // Draw histograms
      hSD_q2pdf->Draw("HIST");
      hSDstarGamma_q2pdf->Draw("HISTsame");
      hSDstarPi0_q2pdf->Draw("HISTsame");
      hB_q2pdf->Draw("HISTsame");
      //c1->Draw();
      c->Draw();
  }
  
  // Let's generate the pseudo-data

  hist_dataB = b_pdf->generateBinned(RooArgList(*x),Nb);  

 
  hist_dataB->Print();
  //hist_dataSM-Print();
  //hist_dataB.Print();
    
  RooDataHist* hist_data = hist_dataB;

  //hist_data->add(*hist_dataSM);
  //hist_data->add(*hist_dataB);
  hist_data->Print();  
  ////////////////////////////
    
    
    
    

  //if(type==1) hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x), Import(*hxminQBBOPdf));
  //if(type==2) hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x), Import(*hxmaxQBBOPdf));
  //if(type==3) hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x), Import(*hXpsQBBOPdf));
  // hist_data->Print();
  
  //==============================================

  //==============================================
  //  Some tests
  
  // // This may help to check the data by fitting including some signal
  // RooFitResult *fitmass = sb_pdf.fitTo(*hist_data, Extended(),Save(kTRUE), Strategy(0), NumCPU(4)); //For small values of signal, the Likelihood shos some bias
  // //sb_pdf.chi2FitTo(*hist_data);  //This is to test a chi2 fit which show smallest bias than Likelihood 

  // TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  // RooPlot* framer = x->frame() ;
  // hist_data->plotOn(framer) ;
  // sb_pdf.plotOn(framer) ;
  // //sb_pdf.plotOn(framer,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  // sb_pdf.plotOn(framer,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed)) ;
  // //w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed)) ;
  // framer->Draw() ;
  // cR1->Draw();


  //return 0;
  // //=============================================
  // //  Limit
  // //=============================================
    
    

    
    
  RooWorkspace w("w");

  w.import(*s1_pdf);
  w.import(*s2_pdf);
  w.import(*s3_pdf);
  w.import(*b_pdf);

  w.factory("fac[1]");
  w.factory("f1[1]");
  w.factory("f2[1]");
  w.factory("f3[1]");


  //w.factory("effenunu[1]");

  Double_t eff_signal_D = 1;
  Double_t eff_signal_DstarGamma = 1;
  Double_t eff_signal_DstarPi0 = 1;

  Double_t factor =1;
  Double_t f1 =1;
  Double_t f2 =1;
  Double_t f3 =1;
  //Lets extract the efficiencies
  
    
    if(type==1)
    {
      eff_signal_D = (hSD_Mmin2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hSDstarGamma_Mmin2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hSDstarPi0_Mmin2pdf->GetEntries())/Ns_generated;
     
    }
    if(type==2)
    {
      eff_signal_D = (hSD_Mmax2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hSDstarGamma_Mmax2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hSDstarPi0_Mmax2pdf->GetEntries())/Ns_generated;
    }
    
    if(type==3)
    {
      eff_signal_D = (hSD_Xpspdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hSDstarGamma_Xpspdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hSDstarPi0_Xpspdf->GetEntries())/Ns_generated;
      //eff_enunu = (hXpsBPdf->GetEntries())/Ntau_generated;
    }

    if(type==4)
    {
      eff_signal_D = (hSD_q2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hSDstarGamma_q2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hSDstarPi0_q2pdf->GetEntries())/Ns_generated;
      //eff_enunu = (hXpsBPdf->GetEntries())/Ntau_generated;
    }

  


  //Double_t effrelVal = eff_signal/(eff_taupinu);

  //cout<<" --> Ns : "<<hxminPdf->GetEntries()<<"  Nnu : "<<hxminBPdf->GetEntries()<<"  Nb : "<<hxminQBOPdf->GetEntries()<<endl;
  cout<<" Signal efficiency for D : "<<eff_signal_D<<" , DstarGamma:  "<<eff_signal_DstarGamma<<"  and DstarPi0 "<<eff_signal_DstarPi0<<endl;
  //cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_sm<<endl;
  cout<<" Data : "<<hist_data->numEntries()<<endl;

  cout<<"NsD: "<<NsD<<endl;
  cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
  cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
  cout<<"Nb: "<<Nb<<endl;

    // Creating a new variable in order to estimate the branching ratio as a POI
  cout<<" Br_tagD : "<<Br_tagD<<", Br_tagDstarGamma: "<<Br_tagDstarGamma<<", Br_tagDstarPi0: "<<Br_tagDstarPi0<<endl;

  Double_t B_eff1 = 1;
  Double_t B_eff2 = 1;
  Double_t B_eff3 = 1;
  Double_t lum_CS = 1;

  B_eff1 = Br_tagD*eff_signal_D ;
  B_eff2 = Br_tagDstarGamma*eff_signal_DstarGamma;
  B_eff3 = Br_tagDstarPi0*eff_signal_DstarPi0; 
  lum_CS = lum*CS_ee_BB;

  factor = 2*lum_CS*(B_eff1 + B_eff2 + B_eff3);

  //factor =  2*lum*CS_ee_BB*((Br_tagD*eff_signal_D)+(Br_tagDstarGamma*eff_signal_DstarGamma)+(Br_tagDstarPi0+eff_signal_DstarPi0));

  cout<<" factor: (2*lum*CS_ee_BB*[Br_tagD*eff_signal_D+Br_tagDstarGamma*eff_signal_DstarGamma+Br_tagDstarPi0*eff_signal_DstarPi0) = "<<factor<<endl;

  f1 = 1.0 *NsD/(NsD+NsDstarGamma+NsDstarPi0);
  f2 = 1.0 *NsDstarGamma/(NsD+NsDstarGamma+NsDstarPi0);
  f3 = 1.0 *NsDstarPi0/(NsD+NsDstarGamma+NsDstarPi0);

  cout<<" f1 = "<<f1<<endl;
  cout<<" f2 = "<<f2<<endl;
  cout<<" f3 = "<<f3<<endl;
  //w.var("effrel")->setVal(effrelVal);
  w.var("fac")->setVal(factor);
  w.var("f1")->setVal(f1);
  w.var("f2")->setVal(f2);
  w.var("f3")->setVal(f3);
  //w.var("effenunu")->setVal(eff_enunu);


    // For Mmin 
  if (type == 1){

    if (lum == 403.9){
      w.factory("mu[4.0e-6,-1.0,1.0]");
      //w.factory("Snu[0.0,0.0,10.0]");
      w.factory("B[43,2.0,35.0]");
    


      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(150);
  
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);
      
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
          
      // to m_x = 0  
      w.var("mu")->setMax(4.0e-5);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
      
    }
      
    if (lum == 10000){
      w.factory("mu[0.0,1.0]");
      //w.factory("Snu[0.0,0.0,10.0]");
      w.factory("B[739.0,600.0,900.0]");
    


      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(150);
  
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);
      

    
      w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
          
      // to m_x = 0  
      w.var("mu")->setMax(2.0e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
      
    }
        
    if (lum == 30000){
      
      w.factory("mu[0,1]");
      //w.factory("S[1e3,0,2.5e+03]");
      w.factory("B[1e2,0,2.0e+02]");

  
      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
  
          
      // to m_x = 0  
      w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.5e-2);
    }
      
        
    if (lum == 50000){
      
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[4e-8,-1.0,1.0]");
      //w.factory("S[0.0,0.0,10.0]");
      w.factory("B[5412.0,2000.0,7000.0]");
  

  
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);

      //w.var("mu")->setMin(0);
      //w.var("mu")->setMin(-0.4e-6);
      // to m_x = 1.0
      //w.var("mu")->setMin(-0.9e-6);
      // to m_x = 3.0
     //w.var("mu")->setMin(-0.9e-6);
      // to m_x = 4.0
      w.var("mu")->setMin(-0.5e-6);
  
  
          
          
      w.var("mu")->setMax(0.6e-6);
      // to m_x = 5.0
      //w.var("mu")->setMax(2e-6);
    }
  }
   
    // For Mmax2   
  if (type == 2){
        
    if (lum == 361.6){
      w.factory("mu[0.0,1.0]");
      //w.factory("Snu[0.0,0.0,10.0]");
      w.factory("B[27,2.0,35.0]");
    


      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(150);
  
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);
      
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
          
      // to m_x = 0  
      w.var("mu")->setMax(4.0e-5);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
      
    }

    if (lum == 10000){
      w.factory("mu[0.0,1.0]");
      w.factory("Snu[0.0,0.0,10.0]");
      w.factory("B[739.0,600.0,900.0]");
  

      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(150);
  
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);
      

    
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
      //w.var("mu")->setMax(1.5e-1);
      
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
          
    }
        
    if (lum == 30000){
      
      w.factory("mu[0,1]");
      w.factory("S[1e3,0,2.5e+03]");
      w.factory("B[1e2,0,2.0e+02]");

  
      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
  
          
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(6.5e-2);
    }  
        
    if (lum == 50000){
      
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[0.0,1.0]");
      //w.factory("S[0.0,0.0,10.0]");
      w.factory("B[5412.0,2000.0,7000.0]");
  

  
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      //w.var("mu")->setVal(2.5e-7);
      w.var("mu")->setMin(0);
  
  
          
          
      //w.var("mu")->setMax(2.5e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(5e-6);
    }

  }

    // For Xps   
  if (type == 3){

    if (lum == 403.9){
      w.factory("mu[4.0e-6,-1.0,1.0]");
      //w.factory("Snu[0.0,0.0,10.0]");
      w.factory("B[43,2.0,55.0]");
    


      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(150);
  
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);
      
      //w.var("mu")->setVal(-10.0e-6);
      w.var("mu")->setMin(-10e-6);
          
      // to m_x = 0  
      w.var("mu")->setMax(2e-5);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
      
    }
        
    if (lum == 10000){
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[0.0,1.0]");
     // w.factory("S[0.0,0.0,10.0]");
      w.factory("B[3696.0,3300.0,3900.0]");

  
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("Snu")->setMin(0);
    
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
      //w.var("mu")->setMax(1.5e-1);
      
      // to m_x = 0  
      w.var("mu")->setMax(4.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
          
    }
        
    if (lum == 30000){
      
      w.factory("mu[0,1]");
      //w.factory("S[1e3,0,2.5e+03]");
      w.factory("B[1e2,0,2.0e+02]");

      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
  
          
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(6.5e-2);
    }  
        
    if (lum == 50000){
      
       //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[4e-8, -1.0,1.0]");
      //w.factory("S[0.0,0.0,10.0]");
      w.factory("B[3048.0,2000.0,5000.0]");
  

  
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      //w.var("mu")->setMin(0);
      // to m_x = 0.0    
      w.var("mu")->setMin(-3e-7);

  
  
      // to m_x = 0.0    
      //w.var("mu")->setMax(7e-7);
      // to m_x = 4.0
      //w.var("mu")->setMax(4e-7);
      // to m_x = 5.0
      w.var("mu")->setMax(5e-7);
    }

  }
    
  // For q2   
  if (type == 4){

    if (lum == 403.9){
      w.factory("mu[0.0,1.0]");
      //w.factory("Snu[0.0,0.0,10.0]");
      w.factory("B[66,2.0,100.0]");
    


      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(150);
  
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);
      
      //w.var("mu")->setVal(1.0e-5);
      //w.var("mu")->setMin(0);
          
      w.var("mu")->setMin(0);
          
      // to m_X = 0.0
      //w.var("mu")->setMin(-1e-5);
      // to m_X = 2.0
      //w.var("mu")->setMin(-5e-6);
      // to m_x = 5.0
      //w.var("mu")->setMin(-2.0e-6);
  
  
      // to m_X = 0.0   
     // w.var("mu")->setMax(2.5e-5);  
      // to m_X = 2.0  
      w.var("mu")->setMax(1.3e-5);
      
    }
        
    if (lum == 10000){
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[0.0,1.0]");
     // w.factory("S[0.0,0.0,10.0]");
      w.factory("B[3696.0,3300.0,3900.0]");

  
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("Snu")->setMin(0);
    
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
      //w.var("mu")->setMax(1.5e-1);
      
      // to m_x = 0  
      w.var("mu")->setMax(4.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
          
    }
        
    if (lum == 30000){
      
      w.factory("mu[0,1]");
      //w.factory("S[1e3,0,2.5e+03]");
      w.factory("B[1e2,0,2.0e+02]");

      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      w.var("mu")->setVal(1.0e-5);
      
      //w.var("mu")->setMin(0);
  
          
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(6.5e-2);
    }  
        
    if (lum == 50000){
      
       //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[4e-8, -1.0,1.0]");
      //w.factory("S[0.0,0.0,10.0]");
      w.factory("B[3048.0,2000.0,5000.0]");
  

  
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

          
          
      //w.var("S")->setVal(Ns);
      //w.var("S")->setMin(0);


      //w.var("mu")->setMin(0);
      // to m_x = 0.0    
      //w.var("mu")->setMin(-3e-7);
      w.var("mu")->setMin(-2e-7);
  
  
      // to m_x = 0.0    
      //w.var("mu")->setMax(6.0e-7);
      // to m_x = 4.0
      //w.var("mu")->setMax(8e-7);
      // to m_x = 5.0
      w.var("mu")->setMax(3e-7);
    }

  }


    // Now, 'pdf' is composite PDF:
    // pdf(x) = Ns * (f1 * S1(x) + f2 * S2(x) + (1 - f1 - f2) * S3(x)) + Nb * B(x)



    w.factory("expr::S('fac*mu',fac, mu)") ;
    w.factory("SUM::s_pdf(f1*s1_pdf,f2*s2_pdf,f3*s3_pdf)");
    w.factory("SUM::model(S*s_pdf, B*b_pdf)");
    //w.factory("SUM::model(S*s_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(mu*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;


  
  
  //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf)") ;

  //w.factory("expr::S('B*mu*effrel',B, mu[0.00001,0,0.02],effrel)") ;
  //w.factory("SUM::model(S*s_pdf,B*b_pdf)") ;
   
  w.pdf("model")->fitTo(*hist_data, Extended()) ;
  Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  // w.var("mu")->setMax(100*mu_val);


  c->cd(2);
  //TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
  RooPlot* frame = w.var("x")->frame() ;
  hist_data->plotOn(frame) ;
  w.pdf("model")->plotOn(frame) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
  //w.pdf("model")->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
  frame->Draw() ;
  //cR11->Draw();

  
  
  // Sets up b-only model 
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*w.pdf("model"));
  b_modelNM.SetParametersOfInterest(*w.var("mu"));
  //b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("Stau"),*w.var("B")));
  //b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("B")));

  b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("B")));
  //b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("S"),*w.var("B")));
  b_modelNM.SetObservables(*w.var("x"));  //This is for 1D limit
  w.var("mu")->setVal(0.0);
  //w.var("mu")->setVal(1.5e-5);
  //w.var("mu")->setMin(0);
  //w.var("mu")->setMax(6.0e-7);
  b_modelNM.SetSnapshot(*w.var("mu"));     // sets up b hypothesis as s = 0

  
   
  // create the alternate (s+b) ModelConfig with given value of s
  double s = 1;
  RooStats::ModelConfig sb_modelNM("S+B_modelNM", &w);
  sb_modelNM.SetPdf(*w.pdf("model"));
  sb_modelNM.SetObservables(*w.var("x"));  //This is for 1D limit
  sb_modelNM.SetParametersOfInterest(*w.var("mu"));
  RooRealVar* poi = (RooRealVar*) sb_modelNM.GetParametersOfInterest()->first();
  w.var("mu")->setVal(s);
  sb_modelNM.SetSnapshot(*w.var("mu"));    // sets up sb hypothesis with given s

  b_modelNM.Print();
  sb_modelNM.Print();
   
  // Test statistic \lambda(s) = -log L(s,\hat\hat{b})/L(\hat{s},\hat{b})
  RooStats::ProfileLikelihoodTestStat profll(*sb_modelNM.GetPdf());
  //profll.SetOneSidedDiscovery(true);
 

  // RooStatsWorkbook Documentation
  //Instantiate a Profile Likelihood interval calculator

  RooStats::ProfileLikelihoodCalculator plCalc(*hist_data, sb_modelNM);
  plCalc.SetConfidenceLevel(0.68);
  // Set the test statistic for the calculator
  //plCalc.SetTestStatistic(profll);
  // Perform the profile likelihood calculation
  // Perform the profile likelihood calculation
  RooStats::LikelihoodInterval* interval = plCalc.GetInterval();

  if (interval) {
      double lowerLimitL = interval->LowerLimit(*poi);
      double upperLimitL = interval->UpperLimit(*poi);
      
      std::cout << "Lower limit: " << lowerLimitL << std::endl;
      std::cout << "Upper limit: " << upperLimitL << std::endl;

      cout << "RESULT: " << 100*plCalc.ConfidenceLevel() << "% interval is : ["<< lowerLimitL << ", "<< upperLimitL << "]" <<endl;
  } else {
      std::cerr << "Failed to retrieve interval." << std::endl;
  }
  //RooStats::LikelihoodInterval* interval = plCalc.GetInterval();
  //double lowerLimitL = interval->LowerLimit(*poi);
  //double upperLimitL = interval->UpperLimit(*poi);
  


  // Now compute using asymptotic formula; b is alt, sb is null
  RooStats::AsymptoticCalculator ac(*hist_data, b_modelNM, sb_modelNM);
  ac.SetOneSided(true);     // Should want one sided (true) for limit
  RooStats::AsymptoticCalculator::SetPrintLevel(-1);
  RooStats::HypoTestResult* asympResult = ac.GetHypoTest();

  asympResult->SetPValueIsRightTail(false);         // appears to do nothing?!
  double asymp_pb = 1. - asympResult->AlternatePValue(); 
  asympResult->SetPValueIsRightTail(true);
  double asymp_psb = asympResult->NullPValue();

  cout << "Results based on asymptotic formulae:" << endl;
  cout << "psb = " << asymp_psb << endl;
  cout << "pbr  = " << asymp_pb << endl;
  double asymp_clb = 1. - asymp_pb;
  double asymp_clsb = asymp_psb;
  double asymp_cls = 9999.;
  if ( asymp_clb > 0 ) {
    asymp_cls = asymp_clsb/asymp_clb;
  }
  cout << "cls = " << asymp_cls << endl;
  cout << endl;

  
  // create hypotest inverter passing the desired calculator (hc or ac)
  RooStats::HypoTestInverter calc(ac);
  calc.SetVerbose(false);
  calc.SetConfidenceLevel(0.95);
  bool useCLs = true;
  calc.UseCLs(useCLs);
  if (useCLs) { profll.SetOneSided(true); }

  int npoints = 150;  // number of points to scan
  // min and max for scan (better to choose smaller intervals)
  double poimin = poi->getMin();
  double poimax = poi->getMax();
  cout << "Doing a fixed scan  in interval : " 
       << poimin << " , " << poimax << endl;
  calc.SetFixedScan(npoints, poimin, poimax);
  RooStats::HypoTestInverterResult* r = calc.GetInterval();

  double upperLimit = r->UpperLimit();
  double ulError = r->UpperLimitEstimatedError();
  cout << "The computed upper limit is: " << upperLimit 
       << " +/- " << ulError << endl;
  
  // compute expected limit
  cout << "Expected upper limits using b-only model : " << endl;
  cout << "median limit = " << r->GetExpectedUpperLimit(0) << endl;
  cout << "med. limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << endl;
  cout << "med. limit (+2 sig) " << r->GetExpectedUpperLimit(2) << endl;
  cout << "med. limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
  cout << "med. limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
  cout << endl;


   




  RooStats::LikelihoodIntervalPlot *plot = new RooStats::LikelihoodIntervalPlot(interval);
  plot->SetNPoints(50); // Use this to reduce sampling granularity (trades speed for precis
  //TCanvas *cR2 = new TCanvas("cR2","cR2",800,600);
  c->cd(3);
  //plot->SetRange(-0.2e-6, 0.3e-6);
  //plot->SetRange(-5e-6,2.0e-6);
  plot->Draw("TF1"); gPad->Draw();

  // plot result of the scan 
  RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
  //TCanvas* c2 = new TCanvas("HypoTestInverter Scan"); 
  //TCanvas *cR3 = new TCanvas("cR3","cR3",800,600);
  c->cd(4);
  //cR3->SetLogy(false);
  plot2->Draw("CLb 2CL");
  //cR3->SaveAs("~/tau_lalpha/eChannel/GetLimit/PlotsCLS/SimpleCLsLimitMmax_m16.pdf");

  Double_t limitval = r->GetExpectedUpperLimit(0);
  return limitval;
  //return;

}

/*
Double_t getLimit2D(Double_t lum, Int_t pdf_type)
{

  Int_t Ns = 0;
  Int_t Npitaunu = 0;
  Int_t Nb = 0;
  //Int_t Ntau =0;
  //Int_t Nbsm =0;

  RooRealVar *x;
  RooRealVar *y;

  RooDataHist *dS;
  //RooDataHist *dBSM;
  RooDataHist *dSMpitaunu;
  //RooDataHist *dSMtau;
  RooDataHist *dB;

    
  RooHistPdf *s_pdf;
  //RooHistPdf *sm_pdf;
  RooHistPdf *pi_taunu_pdf;
  //RooHistPdf *sm_tau_pdf;
  //RooHistPdf *bsm_pdf;
  RooHistPdf *b_pdf;
    

  RooDataHist *hist_dataSM;
  RooDataHist *hist_dataB;
  //RooDataHist *hist_data;

  //============================================
  //   This is for 2D pdf (Mmin^2, Mmax^2)
  //============================================
    //-----------------------------------------------
    //                              Data
    //----------------------------------------------
    
  
    //TFile *fSMpinu_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_pinu_2D.root");
    //TFile *fSMtau_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hSM_tau_2D.root");
    //TFile *fB_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hB_2D.root");

    
    //----------------------------------------------------
    //                         PDF
    //----------------------------------------------------

   
    //TFile *fS_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_2D_m00_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_2D_m10_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_2D_m20_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_2D_m30_pdf.root");
    TFile *fS_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_2D_m40_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hS_2D_m45_pdf.root");
    

  
    TFile *fPitaunu_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hpi_taunu_2D_pdf.root");
    //TFile *fpinunu_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_2D_pdf.root");
    TFile *fB_2DPdf = new TFile("~/BToInv/KChannel/Analysis/Files/RootFiles/hpi_nunu_2D_pdf.root");
    
    

    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m00");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m10");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m20");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m30");
    hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m40");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m45");
   
    
    hPitaunu2Dpdf = (TH2D*)fPitaunu_2DPdf->Get("hpi_taunu_2D");
    //hSMTau2Dpdf = (TH2D*)fSMtau_2DPdf->Get("hSM_tau_2D_1");
    hB2Dpdf = (TH2D*)fB_2DPdf->Get("hpi_nunu_2D");
    
      
      
    //hS2D = (TH2D*)fSm00_2D->Get("hS_2D_m00_2");
    //hSMPinu2D = (TH2D*)fSMpinu_2D->Get("hSM_pinu_2D_2");
    //hSMTau2D = (TH2D*)fSMtau_2D->Get("hSM_tau_2D_2");
    //hB2D = (TH2D*)fB_2D->Get("hB_2D_2");
    
    //hBSM2Dpdf->Add(hSMTau2Dpdf,hB2Dpdf,1,1);
    //hBSM2D->Add(hSMTau2D,hB2D,1,1);
    //hSM2D->Add(hSMPinu2D,hSMTau2D,1,1);
    
    Ns = hS2Dpdf->GetEntries();
    Npitaunu = hPitaunu2Dpdf->GetEntries();
    //Ntau = hSMTau2Dpdf->GetEntries();
    Nb = hB2Dpdf->GetEntries();
    //Nbsm = hBSM2Dpdf->GetEntries();
    
  //cout<<" Ratio of Nbsm/Nenunu :" <<Nb/Nenunu<<endl;
  
  //Int_t NevtsGen = lum*(Nb+((7.78571)*Npitaunu))/Luminosity;
  //Normalizing for a luminosity of 3.56249e4/ab (1M events of B->piunu)
  Npitaunu = lum*((0.849421)*Npitaunu)/Luminosity;
  //Npitaunu = lum*Npitaunu/Luminosity;
 
  Nb = lum*Nb/Luminosity;
  Ns = lum*Ns/Luminosity;
  
    
  //cout<<" Ratio of Nb/Nenunu :" <<Nb/Nenunu<<endl;
  //cout<<"  NevtsGen : "<<NevtsGen<<endl;
  cout<<"  Ns : "<<Ns<<endl;
  cout<<"  Npitaunu : "<<Npitaunu<<endl;
  cout<<"  Npinunu : "<<Nb<<endl;

 
    
  x = new RooRealVar("x", "x", M_minLow, M_minUp);
  y = new RooRealVar("y", "y", M_maxLow, M_maxUp);
  x->setBins(nbinsL);
  y->setBins(nbinsH);
  
  dS = new RooDataHist("dS", "dS", RooArgList(*x,*y), Import(*hS2Dpdf));
  dSMpitaunu = new RooDataHist("dSMpitaunu", "dSMpitaunu", RooArgList(*x,*y), Import(*hPitaunu2Dpdf));
  //dSMtau = new RooDataHist("dSMtau", "dSMtau", RooArgList(*x,*y), Import(*hSMTau2Dpdf));
  //dBSM = new RooDataHist("dBSM", "dBSM", RooArgList(*x,*y), Import(*hBSM2Dpdf));
  dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hB2Dpdf));
  //dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hMQBBOPdf));
  
  s_pdf = new RooHistPdf("s_pdf","signal pdf",RooArgList(*x,*y),*dS,2);
  pi_taunu_pdf = new RooHistPdf("pi_taunu_pdf","taupinu 2D pdf",RooArgList(*x,*y),*dSMpitaunu,2);
  //sm_tau_pdf = new RooHistPdf("sm_tau_pdf","tau background pdf",RooArgList(*x,*y),*dSMtau,2);
  //bsm_pdf = new RooHistPdf("bsm_pdf","tau+background",RooArgList(*x,*y),*dBSM,2); 
  b_pdf = new RooHistPdf("b_pdf","pinunu pdf",RooArgList(*x,*y),*dB,2);

  //RooRealVar NsR("NsR", "fraction signal", Ns, 0, 3*Ns);
  //RooRealVar NenunuR("NenunuR", "fraction enunu", Nenunu, 0, 3*NevtsGen);
  //RooRealVar NtauR("NtauR", "fraction tau", Ntau, 0, 3*NevtsGen);
  //RooRealVar NsmR("NsmR", "fraction", Nsm, 0 , 3*NevtsGen);
  //RooRealVar NbR("NbR", "fraction", Nb, 0 , 3*NevtsGen);
  //RooRealVar NbR("NbR", "fraction", NevtsGen, 0 , 5*NevtsGen);
      
  //RooAddPdf sb_pdf("sb_pdf", "Signal+Bkg", RooArgList(*s_pdf, *b_pdf), RooArgList(NsR, NbR));
  //RooAddPdf sb_pdf("sb_pdf", "Signal+Bkg", RooArgList(*s_pdf, *sm_pdf, *b_pdf), RooArgList(NsR, NenunuR, NbR));
  //RooAddPdf smb_pdf("smb_pdf", "Signal+Bkg", RooArgList(*sm_pdf, *b_pdf), RooArgList(NenunuR, NbR));
  //RooAddPdf smb_pdf("smb_pdf", "Signal+tau+Bkg", RooArgList(*sm_pinu_pdf, *sm_tau_pdf,  *b_pdf), RooArgList(NenunuR, NtauR, NbR));
  //RooAddPdf smb_pdf("smb_pdf", "Signal+tau+Bkg", RooArgList(*sm_pdf,  *b_pdf), RooArgList(NenunuR, NsmR));
    
  RooRandom::randomGenerator()->SetSeed(1000);
    
    
    
    ////////////////////////////
  
  if(pdf_type == 1){

    hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x,*y),Npitaunu);
    hist_dataB = b_pdf->generateBinned(RooArgList(*x,*y),1);  


  }
  else{
    hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x,*y),Npitaunu);
    hist_dataB = b_pdf->generateBinned(RooArgList(*x,*y),Nb);  
      
    
  }
  hist_dataSM->Print();
  hist_dataB->Print();

  RooDataHist* hist_data = hist_dataSM;

    //hist_data->add(*hist_dataSM);
  hist_data->add(*hist_dataB);

  hist_data->Print(); 


  // RooDataHist* hist_data;
  // hist_data = smb_pdf.generateBinned(RooArgList(*x,*y),NevtsGen);
  //hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x,*y), Import(*hMQBBOPdf));
  // hist_data->Print();
  
  //==============================================

  //==============================================
  //  Some tests
  
  // // This may help to check the data by fitting including some signal
  // RooFitResult *fitmass = sb_pdf.fitTo(*hist_data, Extended(),Save(kTRUE), Strategy(0), NumCPU(4)); //For small values of signal, the Likelihood shos some bias
  // //sb_pdf.chi2FitTo(*hist_data);  //This is to test a chi2 fit which show smallest bias than Likelihood 

  // TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  // RooPlot* framer = x->frame() ;
  // hist_data->plotOn(framer) ;
  // sb_pdf.plotOn(framer) ;
  // //sb_pdf.plotOn(framer,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  // sb_pdf.plotOn(framer,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed)) ;
  // //w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed)) ;
  // framer->Draw() ;
  // cR1->Draw();


  //return;
  // //=============================================
  // //  Limit
  // //=============================================


  RooWorkspace w("w");

  if(pdf_type == 1){
    w.import(*s_pdf);
    w.import(*pi_taunu_pdf);
    //w.import(*sm_tau_pdf);
    //w.import(*bsm_pdf);
    //w.import(*b_pdf);



    w.factory("effrel[1]");

  }
  if(pdf_type == 2){

    w.import(*s_pdf);
    w.import(*pi_taunu_pdf);
    //w.import(*sm_tau_pdf);
    //w.import(*bsm_pdf);
    w.import(*b_pdf);



    //w.factory("effrel[1]");
    //w.factory("fac[1]");
  }
  else{

    w.import(*s_pdf);
    w.import(*pi_taunu_pdf);
    //w.import(*sm_tau_pdf);
    //w.import(*bsm_pdf);
    w.import(*b_pdf);



    //w.factory("effrel[1]");
    w.factory("fac[1]");
  }

  //Lets extract the efficiencies
  
  Double_t eff_signal = (hS2Dpdf->GetEntries())/Ns_generated;
  Double_t eff_taupinu = (hPitaunu2Dpdf->GetEntries())/NB_generated;
  //Double_t eff_taupinu = (hPitaunu2Dpdf->GetEntries()*(7.78571))/NB_generated;
  //Double_t eff_tau = (hSMTau2Dpdf->GetEntries())/Ntau_generated;

  //Double_t eff_bkg = (hxminQBPdf->GetEntries())/Ntau_generated;
  Double_t effrelVal = eff_signal/(eff_taupinu);

  //cout<<" --> Ns : "<<hS2Dpdf->GetEntries()<<"  Ntaupinu : "<<hPitaunu2Dpdf->GetEntries()<<"  Nb : "<<hB2Dpdf->GetEntries()<<endl;
  cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_taupinu<<endl;
  //cout<<" Data : "<<hist_data->numEntries()<<"  vs  "<<5*(Nb+Nenunu+Ntau)<<endl;

    
  Double_t factor =  2*lum*CS_ee_BB*Br_B_to_l_nu_D0;
  
  cout<<" factor: (2*lum*CS_ee_BB*Br_B_to_l_nu_D0) = "<<factor<<endl;
  


  if(pdf_type == 1){

    w.var("effrel")->setVal(effrelVal);
    if (lum == 10000){
      w.factory("mu[0,1]");
      w.factory("Snu[100,0,270]");
      //w.factory("B[20,0,40]");
  


      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);
 
      w.var("Snu")->setVal(Npitaunu);
      w.var("Snu")->setMin(0);
    

  
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
      //w.var("mu")->setMax(1.5e-1);
      
        // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
        //w.var("mu")->setMax(2.2e-1);
        // to m_x = 4.0
        //w.var("mu")->setMax(1.4e-1);
        // to m_x = 5.0
      w.var("mu")->setMax(6.0e-2);
    
    }
    
    if (lum == 30000){
    
      w.factory("mu[0,1]");
      w.factory("Snu[1e3,0,2.5e+03]");
      //w.factory("B[1e2,0,2.0e+02]");

 
      //w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

        
        
      w.var("Snu")->setVal(Npitaunu);
      w.var("Snu")->setMin(0);


      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
 
        
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(6.5e-2);
    }
    
    if (lum == 50000){
    
        
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[0,1]");
      w.factory("Snu[1e3,0,2.5e+03]");
      //w.factory("B[1e3,0,2.0e+03]");
      //w.factory("Snu[1167]");
      //w.factory("B[149]");

 
      //w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

        
        
      w.var("Snu")->setVal(Npitaunu);
      //w.var("Snu")->setMin(0);


      w.var("mu")->setVal(2.5e-3);
      //w.var("mu")->setMin(0);
 
 
        
        
      //w.var("mu")->setMax(2.5e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(5.0e-2);
    }
  
    //w.factory("expr::S('fac*mu',fac, mu)") ;

    w.factory("expr::S('Snu*mu*effrel',Snu, mu,effrel)") ;
    w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;

  }

  if(pdf_type == 2){

    //w.var("effrel")->setVal(effrelVal);
    if (lum == 10000){
      w.factory("mu[0,1]");
      w.factory("Snu[100,0,270]");
      //w.factory("B[20,0,40]");
  


      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);
 
      w.var("Snu")->setVal(Npitaunu);
      w.var("Snu")->setMin(0);
    

  
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
      //w.var("mu")->setMax(1.5e-1);
      
        // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
        //w.var("mu")->setMax(2.2e-1);
        // to m_x = 4.0
        //w.var("mu")->setMax(1.4e-1);
        // to m_x = 5.0
      w.var("mu")->setMax(6.0e-2);
    
    }
    
    if (lum == 30000){
    
      w.factory("mu[0,1]");
      w.factory("Snu[1e3,0,2.5e+03]");
      //w.factory("B[1e2,0,2.0e+02]");

 
      //w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

        
        
      w.var("Snu")->setVal(Npitaunu);
      w.var("Snu")->setMin(0);


      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
 
        
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(6.5e-2);
    }
    
    if (lum == 50000){
    
        
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[0,1]");
      w.factory("Snu[1e3,0,2.5e+03]");
      w.factory("B[1e3,0,2.0e+03]");
      //w.factory("Snu[1167]");
      //w.factory("B[149]");

 
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

        
        
      w.var("Snu")->setVal(Npitaunu);
      //w.var("Snu")->setMin(0);


      w.var("mu")->setVal(2.5e-7);
      //w.var("mu")->setMin(0);
 
 
        
        
      //w.var("mu")->setMax(2.0e-6);
      // to m_x = 5.0
      w.var("mu")->setMax(8.0e-7);
    }
  

    w.factory("SUM::model(mu*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;;

  }

  if(pdf_type == 3){
    //w.var("effrel")->setVal(effrelVal);
    w.var("fac")->setVal(factor);

    if (lum == 10000){
      w.factory("mu[0,1]");
      w.factory("Snu[100,0,270]");
      //w.factory("B[20,0,40]");
  


      w.var("B")->setVal(Nb);
      w.var("B")->setMin(0);
 
      w.var("Snu")->setVal(Npitaunu);
      w.var("Snu")->setMin(0);
    

  
      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
      //w.var("mu")->setMax(1.5e-1);
      
        // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
        //w.var("mu")->setMax(2.2e-1);
        // to m_x = 4.0
        w.var("mu")->setMax(1.4e-1);
        // to m_x = 5.0
      //w.var("mu")->setMax(6.0e-2);
    
    }
    
    if (lum == 30000){
    
      w.factory("mu[0,1]");
      w.factory("Snu[1e3,0,2.5e+03]");
      //w.factory("B[1e2,0,2.0e+02]");

 
      //w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

        
        
      w.var("Snu")->setVal(Npitaunu);
      w.var("Snu")->setMin(0);


      //w.var("mu")->setVal(1.0e-5);
      w.var("mu")->setMin(0);
 
        
      // to m_x = 0  
      //w.var("mu")->setMax(1.5e-1);
      // to m_x = 1.0
      //w.var("mu")->setMax(2.0e-1);
      // to m_x = 2.0, 3.0
      //w.var("mu")->setMax(2.2e-1);
      // to m_x = 4.0
      //w.var("mu")->setMax(1.4e-1);
      // to m_x = 5.0
      w.var("mu")->setMax(6.5e-2);
    }
    
    if (lum == 50000){
    
        
      //w.factory("mu[-0.5e-6,1]");
      w.factory("mu[0,1]");
      w.factory("Snu[1e3,0,2.5e+03]");
      w.factory("B[1e3,0,2.0e+03]");
      //w.factory("Snu[1167]");
      //w.factory("B[149]");

 
      w.var("B")->setVal(Nb);
      //w.var("B")->setMin(0);

        
        
      w.var("Snu")->setVal(Npitaunu);
      //w.var("Snu")->setMin(0);


      w.var("mu")->setVal(2.5e-7);
      //w.var("mu")->setMin(0);
 
 
        
        
      //w.var("mu")->setMax(2.0e-6);
      // to m_x = 5.0
      w.var("mu")->setMax(8.0e-7);
    }
  
    w.factory("expr::S('fac*mu',fac, mu)") ;
    w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;;

  }


   
  w.pdf("model")->fitTo(*hist_data, Extended()) ;

  //Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  //w.var("mu")->setMax(100*mu_val);
 
  if(pdf_type == 2 && pdf_type == 3){

    TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
    RooPlot* frame = w.var("x")->frame() ;
    hist_data->plotOn(frame) ;
    w.pdf("model")->plotOn(frame) ;
    w.pdf("model")->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
    w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
    w.pdf("model")->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
    frame->Draw() ;
    cR11->Draw();

  

  }
  else{
  
    TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
    RooPlot* frame = w.var("x")->frame() ;
    hist_data->plotOn(frame) ;
    w.pdf("model")->plotOn(frame) ;
    w.pdf("model")->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
    //w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
    w.pdf("model")->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
    frame->Draw() ;
    cR11->Draw();

  }

  // b-only model construction
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*w.pdf("model"));
  b_modelNM.SetParametersOfInterest(*w.var("mu"));

  if(pdf_type == 1){
    b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu")));
  }
  else{
    b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("B")));
  }

  b_modelNM.SetObservables(RooArgSet(*x,*y)); //This is for 2D limit

  w.var("mu")->setVal(0.0);
  b_modelNM.SetSnapshot(*w.var("mu"));     // sets up b hypothesis as s = 0

   
  // create the alternate (s+b) ModelConfig with given value of s
  double s = 1;
  RooStats::ModelConfig sb_modelNM("S+B_modelNM", &w);
  sb_modelNM.SetPdf(*w.pdf("model"));
  sb_modelNM.SetObservables(RooArgSet(*x,*y)); //This is for 2D limit
  sb_modelNM.SetParametersOfInterest(*w.var("mu"));
  RooRealVar* poi = (RooRealVar*) sb_modelNM.GetParametersOfInterest()->first();
  w.var("mu")->setVal(s);
  sb_modelNM.SetSnapshot(*w.var("mu"));    // sets up sb hypothesis with given s

  b_modelNM.Print();
  sb_modelNM.Print();
   
  // Test statistic \lambda(s) = -log L(s,\hat\hat{b})/L(\hat{s},\hat{b})
  RooStats::ProfileLikelihoodTestStat profll(*sb_modelNM.GetPdf());

  // Now compute using asymptotic formula; b is alt, sb is null
  RooStats::AsymptoticCalculator ac(*hist_data, b_modelNM, sb_modelNM);
  ac.SetOneSided(true);     // Should want one sided (true) for limit
  RooStats::AsymptoticCalculator::SetPrintLevel(-1);
  RooStats::HypoTestResult* asympResult = ac.GetHypoTest();

  asympResult->SetPValueIsRightTail(false);         // appears to do nothing?!
  double asymp_pb = 1. - asympResult->AlternatePValue(); 
  asympResult->SetPValueIsRightTail(true);
  double asymp_psb = asympResult->NullPValue();

  cout << "Results based on asymptotic formulae:" << endl;
  cout << "psb = " << asymp_psb << endl;
  cout << "pbr  = " << asymp_pb << endl;
  double asymp_clb = 1. - asymp_pb;
  double asymp_clsb = asymp_psb;
  double asymp_cls = 9999.;
  if ( asymp_clb > 0 ) {
    asymp_cls = asymp_clsb/asymp_clb;
  }
  cout << "cls = " << asymp_cls << endl;
  cout << endl;

  
  // create hypotest inverter passing the desired calculator (hc or ac)
  RooStats::HypoTestInverter calc(ac);
  calc.SetVerbose(false);
  calc.SetConfidenceLevel(0.95);
  bool useCLs = true;
  calc.UseCLs(useCLs);
  if (useCLs) { profll.SetOneSided(true); }

  int npoints = 100;  // number of points to scan
  // min and max for scan (better to choose smaller intervals)
  double poimin = poi->getMin();
  double poimax = poi->getMax();
  cout << "Doing a fixed scan  in interval : " 
       << poimin << " , " << poimax << endl;
  calc.SetFixedScan(npoints, poimin, poimax);
  RooStats::HypoTestInverterResult* r = calc.GetInterval();

  double upperLimit = r->UpperLimit();
  double ulError = r->UpperLimitEstimatedError();
  cout << "The computed upper limit is: " << upperLimit 
       << " +/- " << ulError << endl;
  
  // compute expected limit
  cout << "Expected upper limits using b-only model : " << endl;
  cout << "median limit = " << r->GetExpectedUpperLimit(0) << endl;
  cout << "med. limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << endl;
  cout << "med. limit (+2 sig) " << r->GetExpectedUpperLimit(2) << endl;
  cout << "med. limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
  cout << "med. limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
  cout << endl;

  RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
  TCanvas *cR3 = new TCanvas("cR3","cR3",800,600);
  cR3->SetLogy(false);
  plot2->Draw("CLb 2CL");
  //cR3->SaveAs("~/tau_lalpha/eChannel/GetLimit/PlotsCLS/SimpleCLsLimit2D_m10.pdf"); 

  Double_t limitval = r->GetExpectedUpperLimit(0);
  return limitval;

}*/

  
TVector3 calculateThrust(const std::vector<TVector3> momenta)
{

  decltype(momenta.begin()) p;
  decltype(momenta.begin()) q;
  decltype(momenta.begin()) loopcount;

  const auto begin = momenta.begin();
  const auto end = momenta.end();

  double sump = 0;
  for (p = begin; p != end; p++)
    sump += (*p).Mag();


  TVector3 Axis;

  // Thrust and thrust vectors

  double Thru = 0;
  for (p = begin; p != end; p++) {
    TVector3 rvec(*p);
    if (rvec.z() <= 0.0) rvec = -rvec;

    double s = rvec.Mag();
    if (s != 0.0) rvec *= (1 / s);

    for (loopcount = begin; loopcount != end; loopcount++) {
      TVector3 rprev(rvec);
      rvec = TVector3(); // clear

      for (q = begin; q != end; q++) {
        const TVector3 qvec(*q);
        rvec += (qvec.Dot(rprev) >= 0) ? qvec : - qvec;
      }

      for (q = begin; q != end; q++) {
        const TVector3 qvec(*q);
        if (qvec.Dot(rvec) * qvec.Dot(rprev) < 0) break;
      }

      if (q == end) break;
    }

    double ttmp = 0.0;
    for (q = begin; q != end; q++) {
      const TVector3 qvec = *q;
      ttmp += std::fabs(qvec.Dot(rvec));
    }
    ttmp /= (sump * rvec.Mag());
    rvec *= 1 / rvec.Mag();
    if (ttmp > Thru) {
      Thru = ttmp;
      Axis = rvec;
    }
  }
  Axis *= Thru;
  return Axis;
}

Double_t get_X_PseudoRestFrame(TLorentzVector piAll, TLorentzVector el)
{
  Double_t Ebeam = 10.58;
  Double_t E_tau = Ebeam/2.0;
  Double_t M_tau = 1.7768;
  Double_t P_tau = TMath::Sqrt(E_tau*E_tau - M_tau*M_tau);
  Double_t betav = P_tau/E_tau;
  TVector3 boostPS=betav*(piAll.BoostVector().Unit());
  el.Boost(boostPS);
  Double_t Xval = 2.0*el.E()/M_tau;
  return Xval;
}


void GetHistogram( DataGen3x1 *t, TH1D *hS, Int_t type, Int_t sample)
{
 
  Bool_t IsSignal = false;
  TRandom3 ra;
  ra.SetSeed(1);
  Double_t x0 = (tau_m/CMS_E)*(tau_m/CMS_E);
  
  std::vector<TVector3> momenta;
  TVector3 mom(0,0,0);
  TLorentzVector ptmp(0,0,0,0);
  TLorentzVector pa,pb;

  Int_t nentries = t->fChain->GetEntries();
  cout<<"No. events in the tree: "<<nentries<<endl;

  for(Long64_t jentry=0; jentry!=nentries;++jentry)
    {
      Long64_t ientry = t->LoadTree(jentry);
      if (ientry < 0) break;
      Int_t nb = t->fChain->GetEntry(jentry);

      momenta.clear();
      pa.Clear();
      pb.Clear();
      
      for(Int_t ik = 0;ik<t->nchparts;ik++)
	{
	  mom.Clear();
	  mom.SetXYZ(t->px[ik],t->py[ik],t->pz[ik]);
	  momenta.push_back(mom);
	}      

      TVector3 Thrust = calculateThrust(momenta);

      Int_t nac=0;
      Int_t nbc=0;

      TLorentzVector patmp(0,0,0,0);
      TLorentzVector pbtmp(0,0,0,0);

      
      for(Int_t ik = 0;ik<t->nchparts;ik++)
	{

	  mom.Clear();
	  ptmp.Clear();
	  patmp.Clear();
	  pbtmp.Clear();

	  //Smearing
	  Double_t pT = TMath::Sqrt((t->px[ik])*(t->px[ik]) + (t->py[ik])*(t->py[ik]));
	  Double_t sigma = 0.01*pT/TMath::Sqrt(2);
	  Double_t px = t->px[ik] + ra.Gaus(0,sigma);
	  Double_t py = t->py[ik] + ra.Gaus(0,sigma);
	  Double_t pz = t->pz[ik] + ra.Gaus(0,sigma);
	  Double_t mp = h_m;

	  if(fabs(t->pdgID[ik])==11) mp = e_m;
	  if(fabs(t->pdgID[ik])==211) mp = h_m;
		  
	  Double_t Energy = TMath::Sqrt(px*px + py*py + pz*pz + mp*mp);
	  
	  mom.SetXYZ(t->px[ik],t->py[ik],t->pz[ik]);
	  ptmp.SetPxPyPzE(px,py,pz,Energy);

	  //ptmp.SetPxPyPzE(t->px[ik],t->py[ik],t->pz[ik],Ep);

	  if(mom.Dot(Thrust)>0)
	    {
	      patmp+=ptmp;
	      nac++;
	    }
	  else
	    {
	      pbtmp+=ptmp;
	      nbc++;
	    }
	}

      Bool_t passed = false;

      if(nac==3 && nbc==1)
	{
	  pb = patmp;
	  pa = pbtmp;
	  passed = true;
	}

      if(nac==1 && nbc==3)
	{
	  pb = pbtmp;
	  pa = patmp;
	  passed = true;
	}
	  
      if(!passed) continue;
      
      Double_t Etau = CMS_E/2;
      Double_t E3pi = pa.E();
      Double_t P3pi = pa.Vect().Mag();
      Double_t M3pi = pa.M();
        
      TVector3 n_a = (1.0/CMS_E)*pa.Vect();
      TVector3 n_b = (1.0/CMS_E)*pb.Vect();
      Double_t ab = n_a.Dot(n_b);
  
      Double_t za = pa.E()/CMS_E;
      Double_t zb = pb.E()/CMS_E;
      
      Double_t a2 = n_a.Mag2();
      Double_t b2 = n_b.Mag2();
      
      Double_t wa = (zb*zb - zb - b2 - 2*ab);
      Double_t wb = (za*za - za + a2);
  
      TVector3 H;
      H = wa*n_a + wb*n_b; 

      TVector3 acrossb = n_a.Cross(n_b);
  
      Double_t A1 = b2;
      Double_t A2 = a2;
      Double_t A3 = 2.0*ab;
      Double_t B1 = 2.0*(n_b.Dot(H));
      Double_t B2 = 2.0*(n_a.Dot(H));
      Double_t C1 = 4.0*acrossb.Mag2();
      Double_t D1 = H.Dot(H) -C1*(0.5 - za)*(0.5 - za);
    
      Double_t A = A1;
      Double_t B = -B1 + C1 - 2*A1*x0 - A3*x0;
      Double_t C = A3*x0*x0 + A2*x0*x0 + A1*x0*x0 + B2*x0 + B1*x0 + D1;

      if( B*B - 4*A*C < 0 ) { continue;}
	    
      Double_t Ymax = CMS_E*CMS_E*(-B + TMath::Sqrt(B*B - 4*A*C))/(2*A);
      Double_t Ymin = CMS_E*CMS_E*(-B - TMath::Sqrt(B*B - 4*A*C))/(2*A);  

      Double_t xPseudo = get_X_PseudoRestFrame(pb, pa);

      // Bool_t IsTauTau = false;
      // Bool_t IsPiPi0 = false;
      // Bool_t IsInv = false;
      
      // for(Int_t kk=0;kk<t->ndaughtersTauN;kk++)
      // 	{
      // 	  if(fabs(t->daughterTauN[kk])==15) IsTauTau = true;
      // 	  if(fabs(t->daughterTauN[kk])==111) IsPiPi0 = true;
      // 	  if(fabs(t->daughterTauN[kk])==94144) IsInv = true;		  
      // 	}
      // for(Int_t kk=0;kk<t->ndaughtersTauP;kk++)
      // 	{
      // 	  if(fabs(t->daughterTauP[kk])==15) IsTauTau = true;
      // 	  if(fabs(t->daughterTauP[kk])==111) IsPiPi0 = true;
      // 	  if(fabs(t->daughterTauP[kk])==94144) IsInv = true;		  
      // 	}
	      
      //if(IsTauTau) continue;
      //if(IsPiPi0) continue;
      //if(IsInv) continue;//cout<<"  signal ------* "<<endl;
      //if(IsInv) cout<<"  signal ------* "<<endl;

      if(sample==0)
	{
	  if(type == 1) hS->Fill(Ymax);
	  if(type == 2) hS->Fill(Ymin);
	  if(type == 3) hS->Fill(xPseudo);
	}	  
      
      if(sample==1)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(!Is_enunu) continue;

	  if(type == 1) hS->Fill(Ymax);
	  if(type == 2) hS->Fill(Ymin);
	  if(type == 3) hS->Fill(xPseudo);
	  
	}

      if(sample==2)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(Is_enunu) continue;

	  if(type == 1) hS->Fill(Ymax);
	  if(type == 2) hS->Fill(Ymin);
	  if(type == 3) hS->Fill(xPseudo);
	  
	}

      if(sample==3)
	{
	  if(type == 1) hS->Fill(Ymax);
	  if(type == 2) hS->Fill(Ymin);
	  if(type == 3) hS->Fill(xPseudo);
	}	  
      
    }

  return;
}


void GetHistogram2D( DataGen3x1 *t, TH2D *hS, Int_t sample)
{
 
  Bool_t IsSignal = false;
  TRandom3 ra;
  ra.SetSeed(1);
  Double_t x0 = (tau_m/CMS_E)*(tau_m/CMS_E);
  
  std::vector<TVector3> momenta;
  TVector3 mom(0,0,0);
  TLorentzVector ptmp(0,0,0,0);
  TLorentzVector pa,pb;

  Int_t nentries = t->fChain->GetEntries();
  cout<<"No. events in the tree: "<<nentries<<endl;

  for(Long64_t jentry=0; jentry!=nentries;++jentry)
    {
      Long64_t ientry = t->LoadTree(jentry);
      if (ientry < 0) break;
      Int_t nb = t->fChain->GetEntry(jentry);

      momenta.clear();
      pa.Clear();
      pb.Clear();
      
      for(Int_t ik = 0;ik<t->nchparts;ik++)
	{
	  mom.Clear();
	  mom.SetXYZ(t->px[ik],t->py[ik],t->pz[ik]);
	  momenta.push_back(mom);
	}      

      TVector3 Thrust = calculateThrust(momenta);

      Int_t nac=0;
      Int_t nbc=0;

      TLorentzVector patmp(0,0,0,0);
      TLorentzVector pbtmp(0,0,0,0);

      
      for(Int_t ik = 0;ik<t->nchparts;ik++)
	{

	  mom.Clear();
	  ptmp.Clear();
	  patmp.Clear();
	  pbtmp.Clear();

	  //Smearing
	  Double_t pT = TMath::Sqrt((t->px[ik])*(t->px[ik]) + (t->py[ik])*(t->py[ik]));
	  Double_t sigma = 0.01*pT/TMath::Sqrt(2);
	  Double_t px = t->px[ik] + ra.Gaus(0,sigma);
	  Double_t py = t->py[ik] + ra.Gaus(0,sigma);
	  Double_t pz = t->pz[ik] + ra.Gaus(0,sigma);
	  Double_t mp = h_m;

	  if(fabs(t->pdgID[ik])==11) mp = e_m;
	  if(fabs(t->pdgID[ik])==211) mp = h_m;
		  
	  Double_t Energy = TMath::Sqrt(px*px + py*py + pz*pz + mp*mp);
	  
	  mom.SetXYZ(t->px[ik],t->py[ik],t->pz[ik]);
	  ptmp.SetPxPyPzE(px,py,pz,Energy);

	  //ptmp.SetPxPyPzE(t->px[ik],t->py[ik],t->pz[ik],Ep);

	  if(mom.Dot(Thrust)>0)
	    {
	      patmp+=ptmp;
	      nac++;
	    }
	  else
	    {
	      pbtmp+=ptmp;
	      nbc++;
	    }
	}

      Bool_t passed = false;

      if(nac==3 && nbc==1)
	{
	  pb = patmp;
	  pa = pbtmp;
	  passed = true;
	}

      if(nac==1 && nbc==3)
	{
	  pb = pbtmp;
	  pa = patmp;
	  passed = true;
	}
	  
      if(!passed) continue;
      
      Double_t Etau = CMS_E/2;
      Double_t E3pi = pa.E();
      Double_t P3pi = pa.Vect().Mag();
      Double_t M3pi = pa.M();
        
      TVector3 n_a = (1.0/CMS_E)*pa.Vect();
      TVector3 n_b = (1.0/CMS_E)*pb.Vect();
      Double_t ab = n_a.Dot(n_b);
  
      Double_t za = pa.E()/CMS_E;
      Double_t zb = pb.E()/CMS_E;
      
      Double_t a2 = n_a.Mag2();
      Double_t b2 = n_b.Mag2();
      
      Double_t wa = (zb*zb - zb - b2 - 2*ab);
      Double_t wb = (za*za - za + a2);
  
      TVector3 H;
      H = wa*n_a + wb*n_b; 

      TVector3 acrossb = n_a.Cross(n_b);
  
      Double_t A1 = b2;
      Double_t A2 = a2;
      Double_t A3 = 2.0*ab;
      Double_t B1 = 2.0*(n_b.Dot(H));
      Double_t B2 = 2.0*(n_a.Dot(H));
      Double_t C1 = 4.0*acrossb.Mag2();
      Double_t D1 = H.Dot(H) -C1*(0.5 - za)*(0.5 - za);
    
      Double_t A = A1;
      Double_t B = -B1 + C1 - 2*A1*x0 - A3*x0;
      Double_t C = A3*x0*x0 + A2*x0*x0 + A1*x0*x0 + B2*x0 + B1*x0 + D1;

      if( B*B - 4*A*C < 0 ) { continue;}
	    
      Double_t Ymax = CMS_E*CMS_E*(-B + TMath::Sqrt(B*B - 4*A*C))/(2*A);
      Double_t Ymin = CMS_E*CMS_E*(-B - TMath::Sqrt(B*B - 4*A*C))/(2*A);  

      Double_t xPseudo = get_X_PseudoRestFrame(pb, pa);

      // Bool_t IsTauTau = false;
      // Bool_t IsPiPi0 = false;
      // Bool_t IsInv = false;
      
      // for(Int_t kk=0;kk<t->ndaughtersTauN;kk++)
      // 	{
      // 	  if(fabs(t->daughterTauN[kk])==15) IsTauTau = true;
      // 	  if(fabs(t->daughterTauN[kk])==111) IsPiPi0 = true;
      // 	  if(fabs(t->daughterTauN[kk])==94144) IsInv = true;		  
      // 	}
      // for(Int_t kk=0;kk<t->ndaughtersTauP;kk++)
      // 	{
      // 	  if(fabs(t->daughterTauP[kk])==15) IsTauTau = true;
      // 	  if(fabs(t->daughterTauP[kk])==111) IsPiPi0 = true;
      // 	  if(fabs(t->daughterTauP[kk])==94144) IsInv = true;		  
      // 	}
	      
      //if(IsTauTau) continue;
      //if(IsPiPi0) continue;
      //if(IsInv) continue;//cout<<"  signal ------* "<<endl;
      //if(IsInv) cout<<"  signal ------* "<<endl;

      //Bool_t Is_enunu = false;
      //Let's check the decays
      //if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
      //if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

      if(sample==0)
	{
	  hS->Fill(Ymin,Ymax);
	}	  
      
      if(sample==1)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(!Is_enunu) continue;

	  hS->Fill(Ymin,Ymax);
	  
	}

      if(sample==2)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(Is_enunu) continue;

	  hS->Fill(Ymin,Ymax);
	  
	}

      if(sample==3)
	{
	  hS->Fill(Ymin,Ymax);
	}	  
      

    }

  return;
}
