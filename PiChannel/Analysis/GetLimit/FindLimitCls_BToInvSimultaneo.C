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
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <string.h>

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

//#include "/home/belle2/johancol/tau3x1/RooFitClass/RooTauLeptonInvisible.h"
//Adding library: make changes accordingly
//#pragma cling load("/home/belle2/johancol/tau3x1/RooFitClass/RooTauLeptonInvisible.so")


//Global variables

Double_t CMS_E = 10.58;
Double_t tau_m = 1.777;
Double_t h_m = 0.13957;
Double_t e_m = 0.000511;

TChain *chS_Dpdf;
TChain *chS_DMpdf;
TChain *chQQpdf;

TChain *chS_D ;
TChain *chS_DM;
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
//void getLimit(Double_t lum, Int_t type, Int_t pdf_type);
void getLimit(Double_t lum, Int_t type);
TH1* normalize(TH1* h);
//Double_t getLimit2D(Double_t lum, Int_t pdf_type);
TVector3 calculateThrust(const std::vector<TVector3> momenta);
Double_t get_X_PseudoRestFrame(TLorentzVector piAll, TLorentzVector el);
void GetHistogram( DataGen3x1 *t, TH1D *hS_D, Int_t type, Int_t sample);
void GetHistogram2D( DataGen3x1 *t, TH2D *hS_D, Int_t sample);

//========================================
//  Histograms definition
//  Q : qqbar bkg
//  B : tau+tau- bkg
//  QB: Q+B
  
Double_t M_min2Low = -5.0;
Double_t M_min2Up = 28.0;

Double_t M_max2Low = -1.0;
Double_t M_max2Up = 27.0;

  
Double_t Xps_Low = 0.0;
Double_t Xps_Up = 1.2;

  
Double_t M_minLow = -5.0;
Double_t M_minUp = 27.0;

Double_t M_maxLow = -1.0;
Double_t M_maxUp = 27.0;


Int_t nbins = 100;
Int_t nbinsH = 15;
Int_t nbinsL = 15;


// =============================
// For 2D limit

//  For data ------
TH2D *hS_D2D = new TH2D("hS_D2D","hS_D2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hS_DMPinu2D = new TH2D("hS_DMPinu2D","hS_DMPinu2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
//TH2D *hS_DMTau2D = new TH2D("hS_DMTau2D","hS_DMTau2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hB2D = new TH2D("hB2D","hB2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
//TH2D *hBSM2D = new TH2D("hBSM2D","hBSM2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp); 


//  For Pdf  ------
TH2D *hS_D2Dpdf = new TH2D("hS_D2Dpdf","hS_D2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hhBkgStagD2Dpdf = new TH2D("hhBkgStagD2Dpdf","hhBkgStagD2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
//TH2D *hS_DMTau2Dpdf = new TH2D("hS_DMTau2Dpdf","hS_DMTau2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hB2Dpdf = new TH2D("hB2Dpdf","hB2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp); 
//TH2D *hBSM2Dpdf = new TH2D("hBSM2Dpdf","hBSM2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);  




// =============================

// For Mmin^2 limit

// For data ------

TH1D *hS_DMmin2 = new TH1D("hS_DMmin2","hS_DMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hS_DstarGammaMmin2 = new TH1D("hS_DstarGammaMmin2","hS_DstarGammaMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hS_DstarPi0Mmin2 = new TH1D("hS_DstarPi0Mmin2","hS_DstarPi0Mmin2",nbins,M_min2Low,M_min2Up);

TH1D *hBkgStagDMmin2 = new TH1D("hBkgStagDMmin2","hBkgStagDMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hBkgStagDstarGammaMmin2 = new TH1D("hBkgStagDstarGammaMmin2","hBkgStagDstarGammaMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hBkgStagDstarPi0Mmin2 = new TH1D("hBkgStagDstarPi0Mmin2","hBkgStagDstarPi0Mmin2",nbins,M_min2Low,M_min2Up);

// For pdf   ------
TH1D *hS_DMmin2pdf = new TH1D("hS_DMmin2pdf","hS_DMmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hS_DstarGammaMmin2pdf = new TH1D("hS_DstarGammaMmin2pdf","hS_DstarGammaMmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hS_DstarPi0Mmin2pdf = new TH1D("hS_DstarPi0Mmin2pdf","hS_DstarPi0Mmin2pdf",nbins,M_min2Low,M_min2Up);

TH1D *hBkgStagDMmin2pdf = new TH1D("hBkgStagDMmin2pdf","hBkgStagDMmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hBkgStagDstarGammaMmin2pdf = new TH1D("hBkgStagDstarGammaMmin2pdf","hBkgStagDstarGammaMmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hBkgStagDstarPi0Mmin2pdf = new TH1D("hBkgStagDstarPi0Mmin2pdf","hBkgStagDstarPi0Mmin2pdf",nbins,M_min2Low,M_min2Up);

// =============================

// For Mmax^2 limit

// For data ------

TH1D *hS_DMmax2 = new TH1D("hS_DMmax2","hS_DMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hS_DstarGammaMmax2 = new TH1D("hS_DstarGammaMmax2","hS_DstarGammaMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hS_DstarPi0Mmax2 = new TH1D("hS_DMTauMmax2","hS_DMTauMmax2",nbins,M_max2Low,M_max2Up);

TH1D *hBkgStagDMmax2 = new TH1D("hBkgStagDMmax2","hBkgStagDMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hBkgStagDstarGammaMmax2 = new TH1D("hBkgStagDstarGammaMmax2","hBkgStagDstarGammaMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hBkgStagDstarPi0Mmax2 = new TH1D("hBkgStagDstarPi0Mmax2","hBkgStagDstarPi0Mmax2",nbins,M_max2Low,M_max2Up);

// For pdf ------

TH1D *hS_DMmax2pdf = new TH1D("hS_DMmax2pdf","hS_DMmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hS_DstarGammaMmax2pdf = new TH1D("hS_DstarGammaMmax2pdf","hS_DstarGammaMmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hS_DstarPi0Mmax2pdf = new TH1D("hS_DstarPi0Mmax2pdf","hS_DstarPi0Mmax2pdf",nbins,M_max2Low,M_max2Up);

TH1D *hBkgStagDMmax2pdf = new TH1D("hBkgStagDMmax2pdf","hBkgStagDMmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hBkgStagDstarGammaMmax2pdf = new TH1D("hBkgStagDstarGammaMmax2pdf","hBkgStagDstarGammaMmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hBkgStagDstarPi0Mmax2pdf = new TH1D("hBkgStagDstarPi0Mmax2pdf","hBkgStagDstarPi0Mmax2pdf",nbins,M_max2Low,M_max2Up);


// =============================

// For Xps limit

// For data ------

TH1D *hS_DXps = new TH1D("hS_DXps","hS_DXps",nbins,Xps_Low,Xps_Up);
TH1D *hS_DstarGammaXps= new TH1D("hS_DstarGammaXps","hS_DstarGammaXps",nbins,Xps_Low,Xps_Up);
TH1D *hS_DstarPi0Xps= new TH1D("hS_DstarPi0Xps","hS_DstarPi0Xps",nbins,Xps_Low,Xps_Up);

TH1D *hBkgStagDXps = new TH1D("hBkgStagDXps","hBkgStagDXps",nbins,Xps_Low,Xps_Up);
TH1D *hBkgStagDstarGammaXps= new TH1D("hBkgStagDstarGammaXps","hBkgStagDstarGammaXps",nbins,Xps_Low,Xps_Up);
TH1D *hBkgStagDstarPi0Xps = new TH1D("hBkgStagDstarPi0Xps","hBkgStagDstarPi0Xps",nbins,Xps_Low,Xps_Up);

// For pdf ------

TH1D *hS_DXpspdf = new TH1D("hS_DXpspdf","hS_DXpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hS_DstarGammaXpspdf = new TH1D("hS_DstarGammaXpspdf","hS_DstarGammaXpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hS_DstarPi0Xpspdf = new TH1D("hS_DstarPi0Xpspdf","hS_DstarPi0Xpspdf",nbins,Xps_Low,Xps_Up);

TH1D *hBkgStagDXpspdf = new TH1D("hBkgStagDXpspdf","hBkgStagDXpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hBkgStagDstarGammaXpspdf = new TH1D("hBkgStagDstarGammaXpspdf","hBkgStagDstarGammaXpspdf",nbins,Xps_Low,Xps_Up);
TH1D *hBkgStagDstarPi0Xpspdf = new TH1D("hBkgStagDstarPi0Xpspdf","hBkgStagDstarPi0Xpspdf",nbins,Xps_Low,Xps_Up);


//Double_t Br_tau_to_3pi_X = 0.152;
Double_t Br_B_to_l_nu_D0 = 2*0.0235;
Double_t Br_D0_to_pi_K = 0.0395;
Double_t Br_B_to_l_nu_D0star = 2*0.0549;
Double_t Br_D0star_to_D0_gamma = 0.353;
Double_t Br_D0star_to_D0_pi0 = 0.547;
Double_t Br_pi0_to_2gamma = 0.9882;

Double_t Br_tagD = Br_B_to_l_nu_D0*Br_D0_to_pi_K;
Double_t Br_tagDstarGamma = Br_D0star_to_D0_gamma*Br_D0_to_pi_K;
Double_t Br_tagDstarPi0 = Br_D0star_to_D0_pi0*Br_D0_to_pi_K*Br_pi0_to_2gamma;

Double_t Br_B_to_tau_pinu = 0.000109*0.1091;
Double_t Br_B_to_pinunu = 0.000014;
Double_t Br_B_to_X = 2.5e-6;
Double_t Luminosity = 8.20e5; //Normalized to the SemiLepTag_Dstar_pi0 
Double_t CS_ee_BB = 0.565e6; //in fb
Double_t Ns_generated = 100000;//*(Br_B_to_l_nu_D0);
Double_t Ns_generated_BR = 2*Luminosity*CS_ee_BB*Br_B_to_l_nu_D0*Br_B_to_X;
Double_t NB_generated = 2*Luminosity*CS_ee_BB*Br_B_to_l_nu_D0*Br_B_to_tau_pinu;



void FindLimitCls_BToInvSimultaneo()
{
  SetData();
  Double_t limit[4];
  //Double_t lum = 10000; 
  //Double_t lum = 30000; 
  Double_t lum = 50000; 
  

  // Depending on the kind of pdf we can estimnate the UL in normalized POI, Number of
  // Events, or directly the branching ratio, the pdfs used are:

  //   PDF1 = Snu*mu*effrel*S_pdf + Snu*pi_taunu_pdf   <- without background
  //   PDF2 = S*S_pdf + Snu*pi_taunu_pdf + B*pi_nunu_pdf  <- Fitting number of events
  //   PDF3 = fac*mu*S_pdf + Snu*pi_taunu_pdf + B*pi_nunu_pdf   <- Branching ratio as POI

  //Int_t pdf_type = 1;
  //Int_t pdf_type = 2;
  //Int_t pdf_type = 3;

  getLimit(lum,3);

  // 1:Mmin2, 2:Mmax2  3:Xps
  //limit[0] = getLimit(lum,1,pdf_type);
  ///limit[1] = getLimit(lum,2,pdf_type);
  //limit[2] = getLimit(lum,3,pdf_type);
  //limit[3] = getLimit2D(lum,pdf_type);

  cout<<" **************************** "<<endl<<endl;
  //cout<<" Limit : "<<limit[0]<<endl;
  //cout<<" Limit : "<<limit[1]<<endl;
  //cout<<" Limit : "<<limit[2]<<endl;
  //cout<<" Limit : "<<limit[3]<<endl;

  
  return;
}


void SetData()
{
    //------------------------------------------------------------------------
    //                      Mmin and Mmax   PDF
    //------------------------------------------------------------------------
    //******************* D_tag
    //
    TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m00_Mmin2_pdf.root");
    TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m00_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m10_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m10_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m20_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m20_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m30_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m30_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m40_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m40_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m50_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m50_Mmax2_pdf.root");

    //**************** Dstar_gamma_tag
    TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Mmin2_pdf.root");
    TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m50_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m50_Mmax2_pdf.root");
    
    //**************** Dstar_pi0_tag
    TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Mmin2_pdf.root");
    TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m50_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m50_Mmax2_pdf.root");    
    
    // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_Mmin2_pdf.root");
    TFile *fBkgStagDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_Mmax2_pdf.root");
    TFile *fBkgStagDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Mmin2_pdf.root");
    TFile *fBkgStagDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Mmax2_pdf.root");
    TFile *fBkgStagDstarPi0min2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Mmin2_pdf.root");
    TFile *fBkgStagDstarPi0max2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Mmax2_pdf.root");

    
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
    TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m00_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m10_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m20_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m30_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m40_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m50_Xps_pdf.root");
    
    //******************* Dstar_gamma_tag
    //
    TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m50_Xps_pdf.root");

    //******************* Dstar_pi0_tag
    //
    TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m50_Xps_pdf.root");

     // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_Xps_pdf.root");
    TFile *fBkgStagDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Xps_pdf.root");
    TFile *fBkgStagDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Xps_pdf.root");
    
  return; 

}

//Double_t getLimit(Double_t lum, Int_t type, Int_t pdf_type)
void getLimit(Double_t lum, Int_t type)
{

  Int_t NsD = 0;
  Int_t NsDstarGamma =0;
  Int_t NsDstarPi0 =0;
  Int_t NBkgD = 0;
  Int_t NBkgDstarGamma = 0;
  Int_t NBkgDstarPi0 = 0;
  //Int_t Ntau = 0;
  //Int_t Nbsm = 0;

  RooRealVar *x;

  RooDataHist *dSD;
  RooDataHist *dSDstarGamma;
  RooDataHist *dSDstarPi0;
  RooDataHist *dBkgStagD;
  RooDataHist *dBkgStagDstarGamma;
  RooDataHist *dBkgStagDstarPi0;
    
  
  RooDataHist *hist_dataBD;
  RooDataHist *hist_dataBDstarGamma;
  RooDataHist *hist_dataBDstarPi0;
  RooDataHist *hist_dataSD;
  RooDataHist *hist_dataSDstarGamma;
  RooDataHist *hist_dataSDstarPi0;
  //RooDataHist *hist_data;
    


  RooHistPdf *sD_pdf;
  RooHistPdf *sDstarGamma_pdf;
  RooHistPdf *sDstarPi0_pdf;
  RooHistPdf *BkgD_pdf;
  RooHistPdf *BkgDstarGamma_pdf;
  RooHistPdf *BkgDstarPi0_pdf;
  
    //------------------------------------------------------------------------
    //                      Mmin and Mmax   PDF
    //------------------------------------------------------------------------
    //******************* D_tag
    //
    TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m00_Mmin2_pdf.root");
    TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m00_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m10_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m10_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m20_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m20_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m30_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m30_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m40_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m40_Mmax2_pdf.root");
    //TFile *fSDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m50_Mmin2_pdf.root");
    //TFile *fSDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m50_Mmax2_pdf.root");

    //**************** Dstar_gamma_tag
    TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Mmin2_pdf.root");
    TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Mmax2_pdf.root");
    //TFile *fSDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m50_Mmin2_pdf.root");
    //TFile *fSDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m50_Mmax2_pdf.root");
    
    //**************** Dstar_pi0_tag
    TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Mmin2_pdf.root");
    TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Mmax2_pdf.root");
    //TFile *fSDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m50_Mmin2_pdf.root");
    //TFile *fSDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m50_Mmax2_pdf.root");    
    
    // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_Mmin2_pdf.root");
    TFile *fBkgStagDMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_Mmax2_pdf.root");
    TFile *fBkgStagDstarGammaMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Mmin2_pdf.root");
    TFile *fBkgStagDstarGammaMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Mmax2_pdf.root");
    TFile *fBkgStagDstarPi0Mmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Mmin2_pdf.root");
    TFile *fBkgStagDstarPi0Mmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Mmax2_pdf.root");

    
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
    TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m00_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m10_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m20_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m30_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m40_Xps_pdf.root");
    //TFile *fSDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_m50_Xps_pdf.root");
    
    //******************* Dstar_gamma_tag
    //
    TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m00_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m10_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m20_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m30_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m40_Xps_pdf.root");
    //TFile *fSDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_gamma_m50_Xps_pdf.root");

    //******************* Dstar_pi0_tag
    //
    TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m00_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m10_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m20_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m30_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m40_Xps_pdf.root");
    //TFile *fSDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_Dstar_pi0_m50_Xps_pdf.root");

     // Background files
    //**************** Dstar_pi0_tag
    TFile *fBkgStagDXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_Xps_pdf.root");
    TFile *fBkgStagDstarGammaXpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_Xps_pdf.root");
    TFile *fBkgStagDstarPi0Xpspdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_pi0_Xps_pdf.root");
    
  if(type==1)
    {
      // Signal D_tag, Dstar_Gamma_tag and Dstar_pi0_tag
      hS_DMmin2pdf = (TH1D*)fSDMmin2pdf->Get("Mmin2");
      hS_DstarGammaMmin2pdf = (TH1D*)fSDstarGammaMmin2pdf->Get("Mmin2");
      hS_DstarPi0Mmin2pdf = (TH1D*)fSDstarPi0Mmin2pdf->Get("Mmin2");
      // Backgrounds D_tag, Dstar_Gamma_tag and Dstar_pi0_tag
      hBkgStagDMmin2pdf = (TH1D*)fBkgStagDMmin2pdf->Get("Mmin2");
      hBkgStagDstarGammaMmin2pdf = (TH1D*)fBkgStagDstarGammaMmin2pdf->Get("Mmin2");
      hBkgStagDstarPi0Mmin2pdf = (TH1D*)fBkgStagDstarPi0Mmin2pdf->Get("Mmin2");
  
     
      
      //hS_DMmin2 = (TH1D*)fSDMmin2->Get("Mmin2");
      //hS_DMPinuMmin2 = (TH1D*)fSDMPinuMmin2->Get("Mmin2");
      //hS_DMTauMmin2 = (TH1D*)fSDMTauMmin2->Get("Mmin2");
      //hBMmin2 = (TH1D*)fBkgStagDstarPi0min2->Get("Mmin2");
      
      
      
      // Number of events

      NsD = hS_DMmin2pdf->GetEntries();
      NsDstarGamma = hS_DstarGammaMmin2pdf->GetEntries();
      NsDstarPi0 = hS_DstarPi0Mmin2pdf->GetEntries();

      NBkgD = hBkgStagDMmin2pdf->GetEntries();
      NBkgDstarGamma = hBkgStagDstarGammaMmin2pdf->GetEntries();
      NBkgDstarPi0 = hBkgStagDstarPi0Mmin2pdf->GetEntries();

      
      //Int_t c1 =1;
      //Int_t c2 = Nb/Ntau;
      // Add the tau->others + all bkg histograms
      //hBSMMmin2pdf->Add(hS_DMTauMmin2pdf,hBMmin2pdf,1,1);
      //hBSMMmin2->Add(hS_DMTauMmin2,hBMmin2,1,1);
	  //Nbsm = hBSMMmin2pdf->GetEntries();
      
      // For Mmin^2
      x = new RooRealVar("x", "x", M_min2Low, M_min2Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hS_DMmin2pdf));
      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hS_DstarGammaMmin2pdf));
      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hS_DstarPi0Mmin2pdf));

      dBkgStagD = new RooDataHist("dBkgStagD", "dBkgStagD", *x, Import(*hBkgStagDMmin2pdf));
      dBkgStagDstarGamma = new RooDataHist("dBkgStagDstarGamma", "dBkgStagDstarGamma", *x, Import(*hBkgStagDstarGammaMmin2pdf));
      dBkgStagDstarPi0 = new RooDataHist("dBkgStagDstarPi0", "dBkgStagDstarPi0", *x, Import(*hBkgStagDstarPi0Mmin2pdf));

    }
  
  if(type==2)
    { 
      /*GetHistogram(tSpdf, hxmaxPdf, 1, 0);  //1 = Mmax^2
      GetHistogram(tSMpdf, hxmaxBPdf, 1, 1);  
      GetHistogram(tSMpdf, hxmaxBOPdf, 1, 2); 
      GetHistogram(tQQpdf, hxmaxQPdf, 1, 3);  
      
      //GetHistogram(tS, hxmax, 1, 0);  //1 = Mmax^2
      GetHistogram(tSM, hxmaxB, 1, 1);  
      GetHistogram(tSM, hxmaxBO, 1, 2); 
      GetHistogram(tQQ, hxmaxQ, 1, 3);  
      
      hxmaxQBOPdf->Add(hxmaxQPdf,hxmaxBOPdf,1,1);
      hxmaxQBBOPdf->Add(hxmaxQBOPdf,hxmaxBPdf,1,1);
      
      hxmaxQBO->Add(hxmaxQ,hxmaxBO,1,1);	    
      hxmaxQBBO->Add(hxmaxQBO,hxmaxB,1,1);
      */
      
      
      
      // Signal D_tag, Dstar_Gamma_tag and Dstar_pi0_tag
      hS_DMmax2pdf = (TH1D*)fSDMmax2pdf->Get("Mmax2");
      hS_DstarGammaMmax2pdf = (TH1D*)fSDstarGammaMmax2pdf->Get("Mmax2");
      hS_DstarPi0Mmax2pdf = (TH1D*)fSDstarPi0Mmax2pdf->Get("Mmax2");
      // Backgrounds D_tag, Dstar_Gamma_tag and Dstar_pi0_tag
      hBkgStagDMmax2pdf = (TH1D*)fBkgStagDMmax2pdf->Get("Mmax2");
      hBkgStagDstarGammaMmax2pdf = (TH1D*)fBkgStagDstarGammaMmax2pdf->Get("Mmax2");
      hBkgStagDstarPi0Mmax2pdf = (TH1D*)fBkgStagDstarPi0Mmax2pdf->Get("Mmax2");
  
     
      
      //hS_DMmin2 = (TH1D*)fSDMmin2->Get("Mmin2");
      //hS_DMPinuMmin2 = (TH1D*)fSDMPinuMmin2->Get("Mmin2");
      //hS_DMTauMmin2 = (TH1D*)fSDMTauMmin2->Get("Mmin2");
      //hBMmin2 = (TH1D*)fBkgStagDstarPi0min2->Get("Mmin2");
      
      
      
      // Number of events

      NsD = hS_DMmax2pdf->GetEntries();
      NsDstarGamma = hS_DstarGammaMmax2pdf->GetEntries();
      NsDstarPi0 = hS_DstarPi0Mmax2pdf->GetEntries();

      NBkgD = hBkgStagDMmax2pdf->GetEntries();
      NBkgDstarGamma = hBkgStagDstarGammaMmax2pdf->GetEntries();
      NBkgDstarPi0 = hBkgStagDstarPi0Mmax2pdf->GetEntries();

      
      //Int_t c1 =1;
      //Int_t c2 = Nb/Ntau;
      // Add the tau->others + all bkg histograms
      //hBSMMmin2pdf->Add(hS_DMTauMmin2pdf,hBMmin2pdf,1,1);
      //hBSMMmin2->Add(hS_DMTauMmin2,hBMmin2,1,1);
	  //Nbsm = hBSMMmin2pdf->GetEntries();
      
      // For Mmin^2
      x = new RooRealVar("x", "x", M_min2Low, M_min2Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hS_DMmax2pdf));
      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hS_DstarGammaMmax2pdf));
      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hS_DstarPi0Mmax2pdf));

      dBkgStagD = new RooDataHist("dBkgStagD", "dBkgStagD", *x, Import(*hBkgStagDMmax2pdf));
      dBkgStagDstarGamma = new RooDataHist("dBkgStagDstarGamma", "dBkgStagDstarGamma", *x, Import(*hBkgStagDstarGammaMmax2pdf));
      dBkgStagDstarPi0 = new RooDataHist("dBkgStagDstarPi0", "dBkgStagDstarPi0", *x, Import(*hBkgStagDstarPi0Mmax2pdf));
    }

 if(type==3)
    {
     
      // Signal D_tag, Dstar_Gamma_tag and Dstar_pi0_tag
      hS_DXpspdf = (TH1D*)fSDXpspdf->Get("Xps");
      hS_DstarGammaXpspdf = (TH1D*)fSDstarGammaXpspdf->Get("Xps");
      hS_DstarPi0Xpspdf = (TH1D*)fSDstarPi0Xpspdf->Get("Xps");
      // Backgrounds D_tag, Dstar_Gamma_tag and Dstar_pi0_tag
      hBkgStagDXpspdf = (TH1D*)fBkgStagDXpspdf->Get("Xps");
      hBkgStagDstarGammaXpspdf = (TH1D*)fBkgStagDstarGammaXpspdf->Get("Xps");
      hBkgStagDstarPi0Xpspdf = (TH1D*)fBkgStagDstarPi0Xpspdf->Get("Xps");
  
     
      
      //hS_DMmin2 = (TH1D*)fSDMmin2->Get("Mmin2");
      //hS_DMPinuMmin2 = (TH1D*)fSDMPinuMmin2->Get("Mmin2");
      //hS_DMTauMmin2 = (TH1D*)fSDMTauMmin2->Get("Mmin2");
      //hBMmin2 = (TH1D*)fBkgStagDstarPi0min2->Get("Mmin2");
      
      
      
      // Number of events

      NsD = hS_DXpspdf->GetEntries();
      NsDstarGamma = hS_DstarGammaXpspdf->GetEntries();
      NsDstarPi0 = hS_DstarPi0Xpspdf->GetEntries();

      NBkgD = hBkgStagDXpspdf->GetEntries();
      NBkgDstarGamma = hBkgStagDstarGammaXpspdf->GetEntries();
      NBkgDstarPi0 = hBkgStagDstarPi0Xpspdf->GetEntries();

      
      //Int_t c1 =1;
      //Int_t c2 = Nb/Ntau;
      // Add the tau->others + all bkg histograms
      //hBSMMmin2pdf->Add(hS_DMTauMmin2pdf,hBMmin2pdf,1,1);
      //hBSMMmin2->Add(hS_DMTauMmin2,hBMmin2,1,1);
	  //Nbsm = hBSMMmin2pdf->GetEntries();
      
      // For Mmin^2
      x = new RooRealVar("x", "x", Xps_Low, Xps_Up);
      x->setBins(nbins);
      
      dSD = new RooDataHist("dSD", "dSD", *x, Import(*hS_DXpspdf));
      dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hS_DstarGammaXpspdf));
      dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hS_DstarPi0Xpspdf));

      dBkgStagD = new RooDataHist("dBkgStagD", "dBkgStagD", *x, Import(*hBkgStagDXpspdf));
      dBkgStagDstarGamma = new RooDataHist("dBkgStagDstarGamma", "dBkgStagDstarGamma", *x, Import(*hBkgStagDstarGammaXpspdf));
      dBkgStagDstarPi0 = new RooDataHist("dBkgStagDstarPi0", "dBkgStagDstarPi0", *x, Import(*hBkgStagDstarPi0Xpspdf));
    }

  //Int_t NevtsGen = lum*(Nb+((0.849421)*NhBkgStagD))/Luminosity;
  //Normalizing for a luminosity of 3.56249e4/ab (1M events of B->piunu)
  //NBkgD = lum*((0.849421)*NBkgD)/Luminosity;
  //NhBkgStagD = lum*NhBkgStagD/Luminosity;

  //NsD = lum*NsD/Luminosity;
  //NsDstarGamma = lum*NsDstarGamma/Luminosity;
  //NsDstarPi0 = lum*NsDstarPi0/Luminosity;

  NBkgD = lum*NBkgD/Luminosity;
  NBkgDstarGamma = lum*NBkgDstarGamma/Luminosity;
  NBkgDstarPi0 = lum*NBkgDstarPi0/Luminosity;
  
  sD_pdf = new RooHistPdf("sD_pdf","signal pdf",*x,*dSD,2);
  sDstarGamma_pdf = new RooHistPdf("sDstarGamma_pdf","signal Gamma pdf",*x,*dSDstarGamma,2);
  sDstarPi0_pdf = new RooHistPdf("sDstarPi0_pdf","signal Pi0 pdf",*x,*dSDstarPi0,2);

  BkgD_pdf = new RooHistPdf("BkgD_pdf","Bkg D pdf",*x,*dBkgStagD,2);
  BkgDstarGamma_pdf = new RooHistPdf("BkgDstarGamma_pdf","Bkg Gamma pdf",*x,*dBkgStagDstarGamma,2);
  BkgDstarPi0_pdf = new RooHistPdf("BkgDstarPi0_pdf","Bkg Pi0 pdf",*x,*dBkgStagDstarPi0,2);



  
  if(type==1){
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);

    c1->Divide(3);
    c1->cd(1);
    hS_DMmin2pdf->Draw();
    hBkgStagDMmin2pdf->Draw("SAME");

    c1->cd(2);
    hS_DstarGammaMmin2pdf->Draw();
    hBkgStagDstarGammaMmin2pdf->Draw("SAME");

    c1->cd(3);
    hS_DstarPi0Mmin2pdf->Draw();
    hBkgStagDstarPi0Mmin2pdf->Draw("SAME");
    c1->Draw();

      /*TCanvas *c1 = new TCanvas("c1","c1",800,600);
      hS_DMmin2pdf->Draw("HIST");
      hBkgStagDMmin2pdf->Draw("HISTsame");
      
      TCanvas *c2 = new TCanvas("c2","c2",800,600);
      hS_DstarGammaMmin2pdf->Draw("HIST");
      hBkgStagDstarGammaMmin2pdf->Draw("HISTsame");

      TCanvas *c3 = new TCanvas("c3","c3",800,600);
      hS_DstarPi0Mmin2pdf->Draw("HIST");
      hBkgStagDstarPi0Mmin2pdf->Draw("HISTsame");*/
  }
    
  if(type==2){

    TCanvas *c1 = new TCanvas("c1","c1",1200,600);

    c1->Divide(3);
    c1->cd(1);
    hS_DMmax2pdf->Draw();
    hBkgStagDMmax2pdf->Draw("SAME");

    c1->cd(2);
    hS_DstarGammaMmax2pdf->Draw();
    hBkgStagDstarGammaMmax2pdf->Draw("SAME");

    c1->cd(3);
    hS_DstarPi0Mmax2pdf->Draw();
    hBkgStagDstarPi0Mmax2pdf->Draw("SAME");
    c1->Draw();

      /*TCanvas *c4 = new TCanvas("c4","c4",800,600);
      hS_DMmax2pdf->Draw("HIST");
      hBkgStagDMmax2pdf->Draw("HISTsame");
      
      TCanvas *c5 = new TCanvas("c5","c5",800,600);
      hS_DstarGammaMmax2pdf->Draw("HIST");
      hBkgStagDstarGammaMmax2pdf->Draw("HIST");

      TCanvas *c6 = new TCanvas("c6","c6",800,600);
      hS_DstarPi0Mmax2pdf->Draw("HIST");
      hBkgStagDstarPi0Mmax2pdf->Draw("HISTsame");*/
  }
  
  if(type==3){
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);

    c1->Divide(3);
    c1->cd(1);
    hS_DXpspdf->Draw();
    hBkgStagDXpspdf->Draw("SAME");

    c1->cd(2);
    hS_DstarGammaXpspdf->Draw();
    hBkgStagDstarGammaXpspdf->Draw("SAME");

    c1->cd(3);
    hS_DstarPi0Xpspdf->Draw();
    hBkgStagDstarPi0Xpspdf->Draw("SAME");
    c1->Draw();

  }
  
  //if(pdf_type == 1){

    //hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x),NhBkgStagD);
    //hist_dataB = b_pdf->generateBinned(RooArgList(*x),1);  
    //hist_dataS = s_pdf->generateBinned(RooArgList(*x),Ns_generated_BR);
    //RooDataHist* generatedDataHist = BkgD_pdf->generateBinned(RooArgList(*x), NBkgD);

    hist_dataBD = BkgD_pdf->generateBinned(RooArgList(*x),NBkgD);  
    hist_dataBDstarGamma = BkgDstarGamma_pdf->generateBinned(RooArgList(*x),NBkgDstarGamma);  
    hist_dataBDstarPi0 = BkgDstarPi0_pdf->generateBinned(RooArgList(*x),NBkgDstarPi0);  

   
    // Create a TH1D histogram from the RooDataHist
    //TH1D* th1dhist_dataBDstarGamma = (TH1D*)hist_dataBDstarGamma->createHistogram("th1dhist_dataBDstarGamma", *x);
    // Clone the TH1D histogram to avoid any potential issues with ownership
    //TH1D* hist_dataBDstarGamma_TH1D = (TH1D*)th1dhist_dataBDstarGamma->Clone("hist_dataBDstarGamma_TH1D");



    

  ///}
  //else{
    //hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x),NhBkgStagD);
    //hist_dataBD = BkgD_pdf->generateBinned(RooArgList(*x),NBkgD);  
//hist_dataBDstarGamma = BkgDstarGamma_pdf->generateBinned(RooArgList(*x),NBkgDstarGamma);  
   // hist_dataBDstarPi0 = BkgDstarPi0_pdf->generateBinned(RooArgList(*x),NBkgDstarPi0);  

    //hist_dataSD = sD_pdf->generateBinned(RooArgList(*x),NsD);
    //hist_dataSDstarGamma = sDstarGamma_pdf->generateBinned(RooArgList(*x),NsDstarGamma);
    //hist_dataSDstarPi0 = sDstarPi0_pdf->generateBinned(RooArgList(*x),NsDstarPi0);

  //}

  
  // Data histograms for Background
  hist_dataBD->Print();
  hist_dataBDstarGamma->Print();
  hist_dataBDstarPi0->Print();
  // Data histograms for signal
  //hist_dataSD->Print();
  //hist_dataSDstarGamma->Print();
  //hist_dataSDstarPi0->Print();

    
  RooDataHist* hist_dataD = hist_dataBD;
  //RooDataHist* hist_dataDstarGamma = hist_dataSDstarGamma;
  //RooDataHist* hist_dataDstarPi0 = hist_dataSDstarPi0;

  hist_dataD->add(*hist_dataBDstarGamma);
  //hist_dataDstarGamma->add(*hist_dataBDstarGamma);
  //hist_dataDstarPi0->add(*hist_dataBDstarOPi0);

  hist_dataD->Print();  
  //hist_dataDstarGamma->Print();  
  //hist_dataDstarPi0->Print();  
  ////////////////////////////
  RooDataHist* hist_data = hist_dataD;
  hist_data->add(*hist_dataBDstarPi0);
    
    
    

  //if(type==1) hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x), Import(*hxminQBBOPdf));
  //if(type==2) hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x), Import(*hxmaxQBBOPdf));
  //if(type==3) hist_data = new RooDataHist("hist_data","hist_data",RooArgList(*x), Import(*hXpsQBBOPdf));
  // hist_data->Print();
  
  //==============================================

  //==============================================
  //  Some tests
  
  // // This may help to check the data by fitting including some signal
  // RooFitResult *fitmass = sb_pdf.fitTo(*hist_data, Extended(),Save(ktRUE), Strategy(0), NumCPU(4)); //For small values of signal, the Likelihood shos some bias
  // //sb_pdf.chi2FitTo(*hist_data);  //This is to test a chi2 fit which show smallest bias than Likelihood 

  // TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  // RooPlot* framer = x->frame() ;
  // hist_data->plotOn(framer) ;
  // sb_pdf.plotOn(framer) ;
  // //sb_pdf.plotOn(framer,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  // sb_pdf.plotOn(framer,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed)) ;
  // //pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed)) ;
  // framer->Draw() ;
  // cR1->Draw();


  //return 0;
  // //=============================================
  // //  Limit
  // //=============================================
    
    

    
    
 /* RooWorkspace w("w");

  if(pdf_type == 1){
    w.import(*sD_pdf);
    w.import(*BkgD_pdf);
    //w.import(*sm_tau_pdf);
    //w.import(*bsm_pdf);
    //w.import(*b_pdf);



    w.factory("effrel[1]");

  }
  if(pdf_type == 2){

    w.import(*sD_pdf);
    w.import(*BkgD_pdf);
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
  }*/

  //w.factory("effenunu[1]");

  Double_t eff_signal_D = 1;
  Double_t eff_signal_DstarGamma = 1;
  Double_t eff_signal_DstarPi0 = 1;
  Double_t eff_bkg_D = 1;
  Double_t eff_bkg_DstarGamma = 1;
  Double_t eff_bkg_DstarPi0 = 1;
  
  //Lets extract the efficiencies
  

    
  if(type==1)
    {
      eff_signal_D = (hS_DMmin2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hS_DstarGammaMmin2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hS_DstarPi0Mmin2pdf->GetEntries())/Ns_generated;
    }
  if(type==2)
    {
      eff_signal_D = (hS_DMmax2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hS_DstarGammaMmax2pdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hS_DstarPi0Mmax2pdf->GetEntries())/Ns_generated;
    }
    
  if(type==3)
    {
      eff_signal_D = (hS_DXpspdf->GetEntries())/Ns_generated;
      eff_signal_DstarGamma = (hS_DstarGammaXpspdf->GetEntries())/Ns_generated;
      eff_signal_DstarPi0 = (hS_DstarPi0Xpspdf->GetEntries())/Ns_generated;
    }



  //cout<<" --> Ns : "<<hxminPdf->GetEntries()<<"  Nnu : "<<hxminBPdf->GetEntries()<<"  Nb : "<<hxminQBOPdf->GetEntries()<<endl;
  cout<<" Signal efficiency for D : "<<eff_signal_D<<" , DstarGamma:  "<<eff_signal_DstarGamma<<"  and DstarPi0 "<<eff_signal_DstarPi0<<endl;
  //cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_sm<<endl;


  cout<<" # bins of Data D: "<<hist_dataBD->numEntries()<<endl;
  cout<<" # bins of Data DstarGamma: "<<hist_dataBDstarGamma->numEntries()<<endl;
  cout<<" # bins of Data DstarPi0: "<<hist_dataBDstarPi0->numEntries()<<endl;


  //cout<<" # events of Data D: "<<hist_dataBDTH1->Integral()<<endl;
  //cout<<" # bins of Data DstarGamma: "<<hist_dataBDstarGammaTH1->Integral()<<endl;
  //cout<<" # bins of Data DstarPi0: "<<hist_dataBDstarPi0->numEntries()<<endl;


  cout<<"NsD: "<<NsD<<endl;
  cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
  cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
  cout<<"NBkgD: "<<NBkgD<<endl;
  cout<<"NBkgDstarGamma: "<<NBkgDstarGamma<<endl;
  cout<<"NBkgDstarPi0: "<<NBkgDstarPi0<<endl;
    // Creating a new variable in order to estimate the branching ratio as a POI
    
  Double_t factorD =  2*lum*CS_ee_BB*Br_tagD*eff_signal_D;
  Double_t factorDstarGamma =  2*lum*CS_ee_BB*Br_tagDstarGamma*eff_signal_DstarGamma;
  Double_t factorDstarPi0 =  2*lum*CS_ee_BB*Br_tagDstarPi0*eff_signal_DstarPi0;

  cout<<" factor: (2*lum*CS_ee_BB*Br_D*eff_D) = "<<factorD<<endl;
  cout<<" factor: (2*lum*CS_ee_BB*Br_DstarGamma*eff_DstarGamma) = "<<factorDstarGamma<<endl;
  cout<<" factor: (2*lum*CS_ee_BB*Br_DstarPi0*eff_DstarPi0) = "<<factorDstarPi0<<endl;

  /*if(pdf_type == 1){
    w.var("effrel")->setVal(effrelVal);
    //w.var("fac")->setVal(factor);
    //w.var("effenunu")->setVal(eff_enunu);


    // For Mmin 
    if (type == 1){
      
      if (lum == 10000){
        w.factory("mu[0,1]");
        w.factory("Snu[2e2,1.5e2,2.6e2]");
        //w.factory("B[2.7e2,1.6e2,2.8e2]");
    


        //w.var("B")->setVal(Nb);
        //w.var("B")->setMin(150);
  
        w.var("Snu")->setVal(NhBkgStagD);
        w.var("Snu")->setMin(0);
      

    
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
          w.factory("Snu[1e3,0,2.5e+03]");
          //w.factory("B[1e2,0,2.0e+02]");

  
          //w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
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

          
          
          w.var("Snu")->setVal(NhBkgStagD);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-4);
          //w.var("mu")->setMin(0);
  
  
          
        
          w.var("mu")->setMax(5e-1);
          // to m_x = 5.0
          //w.var("mu")->setMax(5e-2);
      }
    }
   
    // For Mmax   
    if (type == 2){
        
      if (lum == 10000){
        w.factory("mu[0,1]");
        w.factory("Snu[100,0,300]");
        //w.factory("B[100,0,300]");
    


        //w.var("B")->setVal(Nb);
        //w.var("B")->setMin(0);
  
        w.var("Snu")->setVal(NhBkgStagD);
        w.var("Snu")->setMin(0);
      

    
        //w.var("mu")->setVal(1.0e-5);
        w.var("mu")->setMin(0);
        //w.var("mu")->setMax(1.5e-1);
      
          // to m_x = 0  
        //w.var("mu")->setMax(1.5e-1);
        // to m_x = 1.0
        //w.var("mu")->setMax(2.0e-1);
        // to m_x = 2.0, 3.0
          w.var("mu")->setMax(2.2e-1);
          // to m_x = 4.0
          //w.var("mu")->setMax(1.4e-1);
          // to m_x = 5.0
          //w.var("mu")->setMax(6.0e-2);
          
      }
        
      if (lum == 30000){
      
          w.factory("mu[0,1]");
          w.factory("Snu[1e3,0,2.5e+03]");
          //w.factory("B[1e2,0,2.0e+02]");

  
          //w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
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

  
         // w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-7);
          //w.var("mu")->setMin(0);
  
  
          
          
          w.var("mu")->setMax(2.5e-1);
          // to m_x = 5.0
          //w.var("mu")->setMax(5e-2);
      }

    }

 
  

    w.factory("expr::S('Snu*mu*effrel',Snu, mu,effrel)") ;
    //w.factory("expr::S('fac*mu',fac, mu)") ;

    w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(mu*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;

  }

  if(pdf_type == 2){
    //w.var("effrel")->setVal(effrelVal);
    //w.var("fac")->setVal(factor);
    //w.var("effenunu")->setVal(eff_enunu);


    // For Mmin 
    if (type == 1){
      
      if (lum == 10000){
        w.factory("mu[0,1]");
        w.factory("Snu[2e2,1.5e2,2.6e2]");
        //w.factory("B[2.7e2,1.6e2,2.8e2]");
    


        //w.var("B")->setVal(Nb);
        //w.var("B")->setMin(150);
  
        w.var("Snu")->setVal(NhBkgStagD);
        w.var("Snu")->setMin(0);
      

    
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
          w.factory("Snu[1e3,0,2.5e+03]");
          //w.factory("B[1e2,0,2.0e+02]");

  
          //w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
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

          
          
          w.var("Snu")->setVal(NhBkgStagD);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-7);
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
          w.var("mu")->setMax(2.0e-6);
      }
    }
   
    // For Mmax   
    if (type == 2){
        
      if (lum == 10000){
        w.factory("mu[0,1]");
        w.factory("Snu[100,0,300]");
        //w.factory("B[100,0,300]");
    


        //w.var("B")->setVal(Nb);
        //w.var("B")->setMin(0);
  
        w.var("Snu")->setVal(NhBkgStagD);
        w.var("Snu")->setMin(0);
      

    
        //w.var("mu")->setVal(1.0e-5);
        w.var("mu")->setMin(0);
        //w.var("mu")->setMax(1.5e-1);
      
          // to m_x = 0  
        //w.var("mu")->setMax(1.5e-1);
        // to m_x = 1.0
        //w.var("mu")->setMax(2.0e-1);
        // to m_x = 2.0, 3.0
          w.var("mu")->setMax(2.2e-1);
          // to m_x = 4.0
          //w.var("mu")->setMax(1.4e-1);
          // to m_x = 5.0
          //w.var("mu")->setMax(6.0e-2);
          
      }
        
      if (lum == 30000){
      
          w.factory("mu[0,1]");
          w.factory("Snu[1e3,0,2.5e+03]");
          //w.factory("B[1e2,0,2.0e+02]");

  
          //w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
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

          
          
          w.var("Snu")->setVal(NhBkgStagD);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-7);
          //w.var("mu")->setMin(0);
  
  
          
          
          //w.var("mu")->setMax(2.0e-6);
          // to m_x = 5.0
          w.var("mu")->setMax(8.0e-7);
      }

    }

 
  

    //w.factory("expr::S('Snu*mu*effrel',Snu, mu,effrel)") ;
    //w.factory("expr::S('fac*mu',fac, mu)") ;

    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    w.factory("SUM::model(mu*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;

  }
  
  if(pdf_type == 3){
    
    //w.var("effrel")->setVal(effrelVal);
    w.var("fac")->setVal(factor);
    //w.var("effenunu")->setVal(eff_enunu);


    // For Mmin 
    if (type == 1){
      
      if (lum == 10000){
        w.factory("mu[0,1]");
        w.factory("Snu[2e2,1.5e2,2.6e2]");
        //w.factory("B[2.7e2,1.6e2,2.8e2]");
    


        //w.var("B")->setVal(Nb);
        //w.var("B")->setMin(150);
  
        w.var("Snu")->setVal(NhBkgStagD);
        w.var("Snu")->setMin(0);
      

    
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
          w.factory("Snu[1e3,0,2.5e+03]");
          //w.factory("B[1e2,0,2.0e+02]");

  
          //w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
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

          
          
          w.var("Snu")->setVal(NhBkgStagD);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-7);
          //w.var("mu")->setMin(0);
  
  
          
          // to m_x = 0  
          //w.var("mu")->setMax(1.5e-2);
        // to m_x = 1.0
        //w.var("mu")->setMax(2.0e-1);
        // to m_x = 2.0, 3.0
          //w.var("mu")->setMax(2.2e-1);
          // to m_x = 4.0
          //w.var("mu")->setMax(1.4e-1);
          // to m_x = 5.0
          w.var("mu")->setMax(3.0e-6);
      }
    }
   
    // For Mmax   
    if (type == 2){
        
      if (lum == 10000){
        w.factory("mu[0,1]");
        w.factory("Snu[100,0,300]");
        //w.factory("B[100,0,300]");
    


        //w.var("B")->setVal(Nb);
        //w.var("B")->setMin(0);
  
        w.var("Snu")->setVal(NhBkgStagD);
        w.var("Snu")->setMin(0);
      

    
        //w.var("mu")->setVal(1.0e-5);
        w.var("mu")->setMin(0);
        //w.var("mu")->setMax(1.5e-1);
      
          // to m_x = 0  
        //w.var("mu")->setMax(1.5e-1);
        // to m_x = 1.0
        //w.var("mu")->setMax(2.0e-1);
        // to m_x = 2.0, 3.0
          w.var("mu")->setMax(2.2e-1);
          // to m_x = 4.0
          //w.var("mu")->setMax(1.4e-1);
          // to m_x = 5.0
          //w.var("mu")->setMax(6.0e-2);
          
      }
        
      if (lum == 30000){
      
          w.factory("mu[0,1]");
          w.factory("Snu[1e3,0,2.5e+03]");
          //w.factory("B[1e2,0,2.0e+02]");

  
          //w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(NhBkgStagD);
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

          
          
          w.var("Snu")->setVal(NhBkgStagD);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-7);
          //w.var("mu")->setMin(0);
  
  
          
          
          //w.var("mu")->setMax(2.0e-6);
          // to m_x = 5.0
          w.var("mu")->setMax(8.0e-7);
      }

    }

 
  

    //w.factory("expr::S('Snu*mu*effrel',Snu, mu,effrel)") ;
    w.factory("expr::S('fac*mu',fac, mu)") ;

    //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    //w.factory("SUM::model(mu*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;
    w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf,B*b_pdf)") ;

  }

  
  //w.factory("SUM::model(S*s_pdf,Snu*pi_taunu_pdf)") ;

  //w.factory("expr::S('B*mu*effrel',B, mu[0.00001,0,0.02],effrel)") ;
  //w.factory("SUM::model(S*s_pdf,B*b_pdf)") ;
   
  pdf->fitTo(*hist_data, Extended()) ;
  Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  // w.var("mu")->setMax(100*mu_val);

  if(pdf_type == 2 && pdf_type == 3){

    TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
    RooPlot* frame = w.var("x")->frame() ;
    hist_data->plotOn(frame) ;
    pdf->plotOn(frame) ;
    pdf->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
    pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
    pdf->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
    frame->Draw() ;
    cR11->Draw();

  }
  else
  {
    TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
    RooPlot* frame = w.var("x")->frame() ;
    hist_data->plotOn(frame) ;
    pdf->plotOn(frame) ;
    pdf->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
    //pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
    pdf->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
    frame->Draw() ;
    cR11->Draw();

  } */
  
  // Lets to construct the model (measurement) with the  poi = mu
  HistFactory::Measurement meas("meas", "meas");
  meas.SetExportOnly(true); //true do not perfor fit to data
  meas.SetPOI("mu");
  meas.AddConstantParam("Lumi");
  //meas.SetLumi(9.0506e06);
  meas.SetLumi(0.0000271518);

  //meas.SetLumi(0.0000165039);
  //meas.SetLumi(0.0000513893);
  //meas.SetLumi(0.0000173662);
  
  //meas.SetLumi(0.00001);
  //meas.SetLumiRelErr(0.00000001);
  //meas.SetLumi(1.0);
  //meas.SetLumiRelErr(0.01);

  
  // convert RooDataHist to a TH1 histogram
  TH1 *hist_dataBDTH1 = hist_dataBD->createHistogram("hist_dataBD", *x);
  TH1 *hist_dataBDstarGammaTH1 = hist_dataBDstarGamma->createHistogram("hist_dataDstarGamma", *x);
  TH1 *hist_dataBDstarPi0TH1 = hist_dataBDstarPi0->createHistogram("hist_dataBDstarPi0", *x);
  //TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  //DstarGammachannel.SetData(hist_dataBDstarGamma);
  //TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  //hist_dataBDstarGammaTH1->Draw("HIST");
  //cR1->Draw();
  //TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  //hist_dataBDTH1->Draw("HIST");
  //cR1->Draw();

  //Let's define the channel for D_tag
  // *****************************  D_tag  ****************
  HistFactory::Channel Dchannel("D_channel");
  Dchannel.SetData(hist_dataBDTH1);

  HistFactory::Sample bkgD("bkgD");
  if (type==1 && lum==10000)
  {
    bkgD.SetHisto(hBkgStagDMmin2pdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,600.0,800.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  if (type==1 && lum==30000)
  {
    bkgD.SetHisto(hBkgStagDMmin2pdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,2000.0,2400.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  if (type==1 && lum==50000)
  {
    bkgD.SetHisto(hBkgStagDMmin2pdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,-3200.0,3900.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  
  if (type==2 && lum==10000)
  {
    bkgD.SetHisto(hBkgStagDMmax2pdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,600.0,800.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  if (type==2 && lum==30000)
  {
    bkgD.SetHisto(hBkgStagDMmax2pdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,2000.0,2400.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  if (type==2 && lum==50000)
  {
    bkgD.SetHisto(hBkgStagDMmax2pdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,0.0,4500.0 );
    bkgD.SetNormalizeByTheory(false);
  }

  if (type==3 && lum==10000)
  {
    bkgD.SetHisto(hBkgStagDXpspdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,600.0,800.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  if (type==3 && lum==30000)
  {
    bkgD.SetHisto(hBkgStagDXpspdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,2000.0,2400.0 );
    bkgD.SetNormalizeByTheory(true);
  }
  if (type==3 && lum==50000)
  {
    bkgD.SetHisto(hBkgStagDXpspdf);
    bkgD.AddNormFactor("Nbkg_D",NBkgD,0.0,4000.0 );
    bkgD.SetNormalizeByTheory(false);
  }

  Dchannel.AddSample(bkgD);  

  // Adding the Signal Dtag Pdf to the model
  Char_t factorD_str[5];
  sprintf(factorD_str, "%f",(double)factorD); //convert double to string
  //meas.AddPreprocessFunction("Ns_D", "NsD", "NsD[0.0,0.0,10.0]");
  meas.AddPreprocessFunction("Ns_D", std::string(factorD_str) + " * mu", "mu[0.0,2.0e-7]");
  // Signal sample D_tag
  // Mmin2 and Respectively lumi

  HistFactory::Sample sigD("sigD");
  if (type==1 && lum==10000)
  {
    sigD.SetHisto(hS_DMmin2pdf);
    sigD.AddNormFactor("Ns_D",0.0,-500.0,500.0  );
    sigD.SetNormalizeByTheory(false);
  }
  if (type==1 && lum==30000)
  {
    sigD.SetHisto(hS_DMmin2pdf);
    sigD.AddNormFactor("Ns_D",0.0,-1000.0,1000.0  );
    sigD.SetNormalizeByTheory(false);
  }
  if (type==1 && lum==50000)
  {
    sigD.SetHisto(hS_DMmin2pdf);
    sigD.AddNormFactor("Ns_D",0.0,-2000.0,2000.0 );
    sigD.SetNormalizeByTheory(false);
  }
  
  if (type==2 && lum==10000)
  {
    sigD.SetHisto(hS_DMmax2pdf);
    sigD.AddNormFactor("Ns_D",0.0,-500.0,500.0 );
    sigD.SetNormalizeByTheory(false);
  }
  if (type==2 && lum==30000)
  {
    sigD.SetHisto(hS_DMmax2pdf);
    sigD.AddNormFactor("Ns_D",0.0,-1000.0,1000.0 );
    sigD.SetNormalizeByTheory(false);
  }
  if (type==2 && lum==50000)
  {
    sigD.SetHisto(hS_DMmax2pdf);
    sigD.AddNormFactor("Ns_D",0.0,0.0,10.0 );
    sigD.SetNormalizeByTheory(false);
  }

  if (type==3 && lum==10000)
  {
    sigD.SetHisto(hS_DXpspdf);
    sigD.AddNormFactor("Ns_D",0.0,-500.0,500.0 );
    sigD.SetNormalizeByTheory(false);
  }
  if (type==3 && lum==30000)
  {
    sigD.SetHisto(hS_DXpspdf);
    sigD.AddNormFactor("Ns_D",0.0,-1000.0,1000.00 );
    sigD.SetNormalizeByTheory(false);
  }
  if (type==3 && lum==50000)
  {
    sigD.SetHisto(hS_DXpspdf);
    sigD.AddNormFactor("Ns_D",0.0,0.0,10.0 );
    sigD.SetNormalizeByTheory(false);
  }
  Dchannel.AddSample(sigD);

  // *****************************  DstarGamma_tag  ****************
  RooStats::HistFactory::Channel DstarGammachannel("DstarGamma_channel");
  DstarGammachannel.SetData(hist_dataBDstarGammaTH1);

  //TH1 *hBkgStagDstarGammaMmin2pdfTH1 = hBkgStagDstarGammaMmin2pdf->createHistogram("hBkgStagDstarGammaMmin2pdf", *x);
  //TH1 *hS_DstarGammaMmin2pdfTH1 = hS_DstarGammaMmin2pdf->createHistogram("hS_DstarGammaMmin2pdf", *x);

  //normalize(hist_dataBDstarGammaTH1);
  //normalize(hS_DstarGammaMmin2pdf);
  //normalize(hBkgStagDstarGammaMmin2pdf);

  HistFactory::Sample bkgDstarGamma("bkgDstarGamma");
  if (type==1 && lum==10000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaMmin2pdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,100.0,300.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  if (type==1 && lum==30000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaMmin2pdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,500.0,700.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  if (type==1 && lum==50000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaMmin2pdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,900.0,1100.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  
  if (type==2 && lum==10000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaMmax2pdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,100.0,300.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  if (type==2 && lum==30000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaMmax2pdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,500.0,700.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  if (type==2 && lum==50000)
  {
    bkgDstarGamma.SetNormalizeByTheory(false);
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaMmax2pdf);
    //bkgDstarGamma.AddNormFactor(Nbkg_DstarGamma);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,0.0,1200.0 );
  }

  if (type==3 && lum==10000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaXpspdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,100.0,300.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  if (type==3 && lum==30000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaXpspdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,500.0,700.0 );
    bkgDstarGamma.SetNormalizeByTheory(true);
  }
  if (type==3 && lum==50000)
  {
    bkgDstarGamma.SetHisto(hBkgStagDstarGammaXpspdf);
    bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma",NBkgDstarGamma,0.0,2000.0 );
    bkgDstarGamma.SetNormalizeByTheory(false);
  }
  DstarGammachannel.AddSample(bkgDstarGamma);  

  // Adding the Signal Dtag Pdf to the model
  Char_t factorDstarGamma_str[20];
  sprintf(factorDstarGamma_str, "%f",(double)factorDstarGamma);
  //sprintf(factorDstarGamma_str, "factorDstarGamma = %f",(double)factorDstarGamma);
  //meas.AddPreprocessFunction("Ns_DstarGamma", "Ns_DstarGamma", "Ns_DstarGamma[0.0,10.0]");
  //meas.AddPreprocessFunction("Ns_DstarGamma", std::string(factorDstarGamma_str) +"* mu", "mu[0.0,1.0]");
  meas.AddPreprocessFunction("Ns_DstarGamma", std::string(factorDstarGamma_str) +"* mu", "mu");
  // Signal sample D_tag
  // Mmin2 and Respectively lumi

  HistFactory::Sample sigDstarGamma("sigDstarGamma");
  if (type==1 && lum==10000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaMmin2pdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,00.0,1.0  );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  if (type==1 && lum==30000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaMmin2pdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.0,2.0  );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  if (type==1 && lum==50000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaMmin2pdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0,10.0 );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  
  if (type==2 && lum==10000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaMmax2pdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.0,1.0 );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  if (type==2 && lum==30000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaMmax2pdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.0,3.0 );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  if (type==2 && lum==50000)
  {
    //sigDstarGamma.SetNormalizeByTheory(true);
    sigDstarGamma.SetHisto(hS_DstarGammaMmax2pdf);
    //sigDstarGamma.AddNormFactor("Ns_DstarGamma" );
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.0,10.0 );
    //sigDstarGamma.SetNormalizeByTheory(false);
  }

  if (type==3 && lum==10000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaXpspdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.0,1.0 );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  if (type==3 && lum==30000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaXpspdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.0,2.00 );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  if (type==3 && lum==50000)
  {
    sigDstarGamma.SetHisto(hS_DstarGammaXpspdf);
    sigDstarGamma.AddNormFactor("Ns_DstarGamma",0.0,0.,10.0 );
    sigDstarGamma.SetNormalizeByTheory(false);
  }
  DstarGammachannel.AddSample(sigDstarGamma);


  // *****************************  DstarPi0_tag  ****************
  HistFactory::Channel DstarPi0channel("DstarPi0_channel");
  DstarPi0channel.SetData(hist_dataBDstarPi0TH1);

  HistFactory::Sample bkgDstarPi0("bkgDstarPi0");
  if (type==1 && lum==10000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Mmin2pdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,600.0,800.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  if (type==1 && lum==30000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Mmin2pdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,2000.0,2200.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  if (type==1 && lum==50000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Mmin2pdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,3400.0,3700.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  
  if (type==2 && lum==10000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Mmax2pdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,600.0,800.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  if (type==2 && lum==30000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Mmax2pdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,2000.0,2200.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  if (type==2 && lum==50000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Mmax2pdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,3400.0,3700.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }

  if (type==3 && lum==10000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Xpspdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,600.0,800.0  );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  if (type==3 && lum==30000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Xpspdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,2000.0,2200.0 );
    bkgDstarPi0.SetNormalizeByTheory(true);
  }
  if (type==3 && lum==50000)
  {
    bkgDstarPi0.SetHisto(hBkgStagDstarPi0Xpspdf);
    bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0",NBkgDstarPi0,0.0,3700.0 );
    bkgDstarPi0.SetNormalizeByTheory(false);
  }

  DstarPi0channel.AddSample(bkgDstarPi0);  

  // Adding the Signal Dtag Pdf to the model
  Char_t factorDstarPi0_str[20];
  //sprintf(factorDstarPi0_str, "factorDstarPi0 = %f",(double)factorDstarPi0);
  sprintf(factorDstarPi0_str, "%f",(double)factorDstarPi0);
  //meas.AddPreprocessFunction("Ns_DstarPi0", std::string(factorDstarPi0_str) +"* mu", "mu[0.0,1.0]");
  //meas.AddConstantParam("Lumi");
  //meas.SetLumi(0.0000173662);
  meas.AddPreprocessFunction("Ns_DstarPi0", std::string(factorDstarPi0_str) +"* mu", "mu");
  // Signal sample D_tag
  // Mmin2 and Respectively lumi

  HistFactory::Sample sigDstarPi0("sigDstarPi0");
  if (type==1 && lum==10000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Mmin2pdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-500.0,500.0  );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  if (type==1 && lum==30000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Mmin2pdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-1000.0,1000.0  );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  if (type==1 && lum==50000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Mmin2pdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-2000.0,2000.0 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  
  if (type==2 && lum==10000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Mmax2pdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-500.0,500.0 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  if (type==2 && lum==30000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Mmax2pdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-1000.0,1000.0 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  if (type==2 && lum==50000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Mmax2pdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-2000.0,2000.0 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }

  if (type==3 && lum==10000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Xpspdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-500.0,500.0 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  if (type==3 && lum==30000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Xpspdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,-1000.0,1000.00 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  if (type==3 && lum==50000)
  {
    sigDstarPi0.SetHisto(hS_DstarPi0Xpspdf);
    sigDstarPi0.AddNormFactor("Ns_DstarPi0",0.0,0.0,10.0 );
    sigDstarPi0.SetNormalizeByTheory(false);
  }
  DstarPi0channel.AddSample(sigDstarPi0);

  // perform measurement
  meas.AddChannel(Dchannel);
  meas.AddChannel(DstarGammachannel);
  meas.AddChannel(DstarPi0channel);
  


    /*HistFactory::Channel Dchannel("D_channel");
    Dchannel.SetData(hist_dataBDTH1);

    // Define a background sample for D
    HistFactory::Sample bkgD("bkgD");
    if (type == 3 && lum == 50000) {
        bkgD.SetHisto(hBkgStagDXpspdf);
        bkgD.AddNormFactor("Nbkg_D", NBkgD, 0.0, 3900.0);
        bkgD.SetNormalizeByTheory(false);
    }
    Dchannel.AddSample(bkgD);

    // Define a signal sample for D
    HistFactory::Sample sigD("sigD");
    if (type == 3 && lum == 50000) {
        sigD.SetHisto(hS_DXpspdf);
        sigD.AddNormFactor("Ns_D", 0.0, 0.0, 10.0);
        sigD.SetNormalizeByTheory(false);
    }
    Dchannel.AddSample(sigD);


    HistFactory::Channel DstarGammachannel("DstarGamma_channel");
    DstarGammachannel.SetData(hist_dataBDstarGammaTH1);

    // Define a background sample for DstarGamma
    HistFactory::Sample bkgDstarGamma("bkgDstarGamma");
    if (type == 3 && lum == 50000) {
        bkgDstarGamma.SetHisto(hBkgStagDstarGammaXpspdf);
        bkgDstarGamma.AddNormFactor("Nbkg_DstarGamma", NBkgDstarGamma, 0.0, 2000.0);
        bkgDstarGamma.SetNormalizeByTheory(false);
    }
    DstarGammachannel.AddSample(bkgDstarGamma);

    // Define a signal sample for DstarGamma
    HistFactory::Sample sigDstarGamma("sigDstarGamma");
    if (type == 3 && lum == 50000) {
        sigDstarGamma.SetHisto(hS_DstarGammaXpspdf);
        sigDstarGamma.AddNormFactor("Ns_DstarGamma", 0.0, 0.0, 10.0);
        sigDstarGamma.SetNormalizeByTheory(false);
    }
    DstarGammachannel.AddSample(sigDstarGamma);

    // Define a channel for the DstarPi0 process
    HistFactory::Channel DstarPi0channel("DstarPi0_channel");
    DstarPi0channel.SetData(hist_dataBDstarPi0TH1);

    // Define a background sample for DstarPi0
    HistFactory::Sample bkgDstarPi0("bkgDstarPi0");
    if (type == 3 && lum == 50000) {
        bkgDstarPi0.SetHisto(hBkgStagDstarPi0Xpspdf);
        bkgDstarPi0.AddNormFactor("Nbkg_DstarPi0", NBkgDstarPi0, 0.0, 3700.0);
        bkgDstarPi0.SetNormalizeByTheory(false);
    }
    DstarPi0channel.AddSample(bkgDstarPi0);

    // Define a signal sample for DstarPi0
    HistFactory::Sample sigDstarPi0("sigDstarPi0");
    if (type == 3 && lum == 50000) {
        sigDstarPi0.SetHisto(hS_DstarPi0Xpspdf);
        sigDstarPi0.AddNormFactor("Ns_DstarPi0", 0.0, 0.0, 10.0);
        sigDstarPi0.SetNormalizeByTheory(false);
    }
    DstarPi0channel.AddSample(sigDstarPi0);

    // Add channels to the measurement
    meas.AddChannel(Dchannel);
    meas.AddChannel(DstarGammachannel);
    meas.AddChannel(DstarPi0channel);*/




  cout << "→←ł→←ł→←ł→←ł→←ł→←ł→←ł→←ł→←ł→>|Begining of the RooStats Model|" << endl;
  meas.PrintTree();
  cout << "→←ł→←ł→←ł→←ł→←ł→←ł→←ł→←ł→←ł→>|End of the RooStats Model|" <<endl;



 
  RooWorkspace* w =  MakeModelAndMeasurementFast(meas);
  cout << "→←ł→←ł→←ł→←ł→←ł→←ł→←ł→←ł→←ł→>|End of the Histfactory|" <<endl;



  w->SetName("w");
  w->data("obsData")->SetName("observed_data");
  RooAbsData* data = w->data("observed_data") ;
  w->Print("t");

  
  //data->Print();

  w->pdf("simPdf")->SetName("model");
  RooAbsPdf* model = w->pdf("model");

  model->fitTo(*data, Extended()); 


  TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
  //RooPlot* frame = x->frame();
  RooPlot* frame = w->var("obs_x_D_channel")->frame() ;
  data->plotOn(frame);
  //model->plotOn(frame);
  frame->Draw() ;
  cR11->Draw();

  //RooStats::ModelConfig* modelConf = (RooStats::ModelConfig*) w->obj("ModelConfig");
  //RooSimultaneous* model = (RooSimultaneous*)modelConf->GetPdf();
  //model->fitTo(*data, Extended()) ;

  //w->pdf("model_Dchannel")->SetName("modelD");
  //RooAbsPdf* modelD = w->pdf("modelD") ;
  //w->pdf("model_DstarGamma_channel")->SetName("modelDstar");
  //RooAbsPdf* modelDstar = w->pdf("modelDstar") ;
  /*w->pdf("model_DstarPi0_channel")->SetName("modelDstar");
  RooAbsPdf* modelDstar = w->pdf("modelDstar") ;
  //modelD->fitTo(*data, Extended()); 
  modelDstar->fitTo(*data, Extended()); 
  
  w->pdf("simPdf[")->SetName("model");
  RooAbsPdf* model = w->pdf("model");

  TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
  //RooPlot* frame = w->var("obs_x_Dchannel")->frame() ;
  //data->plotOn(frame, RooFit::Cut("channelCat==channelCat::Dchannel")) ;
  //RooPlot* frame = w->var("obs_x_DstarGamma_channel")->frame() ;
  //data->plotOn(frame, RooFit::Cut("channelCat==channelCat::DstarGamma_channel")) ;
  RooPlot* frame = w->var("obs_x_DstarPi0_channel")->frame() ;
  data->plotOn(frame, RooFit::Cut("channelCat==channelCat::DstarPi0_channel")) ;
  //modelD->plotOn(frame);
  modelDstar->plotOn(frame);
  //w->pdf("model_DstarGamma_channel")->plotOn(frame);
  //model->plotOn(frame, RooFit::LineStyle(kDashed)) ;
  //model->plotOn(frame,RooFit::Components(""),RooFit::LineStyle(kDashed)) ;
  //pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
  //pdf->plotOn(frame,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  frame->Draw() ;
  cR11->Draw();*/

  

  

  
  // s+b model
  RooStats::ModelConfig* SBModel = (RooStats::ModelConfig*) w->obj("ModelConfig");
  SBModel->SetPdf(*model);
  SBModel->SetName("S+B model rhf");
  RooRealVar* poi = (RooRealVar*) SBModel->GetParametersOfInterest()->first();
  poi->setVal(1.0);  
  SBModel->SetSnapshot( *poi );

  // b only model
  RooStats::ModelConfig* BModel = (RooStats::ModelConfig*) SBModel->Clone();
  BModel->SetPdf(*model);
  BModel->SetName("B model rhf");
  poi->setVal(0.0);  
  poi->setError(0.0);
  BModel->SetSnapshot( *poi );

  //BModel->Print();
  //SBModel->Print();

  // Test statistic \lambda(s) = -log L(s,\hat\hat{b})/L(\hat{s},\hat{b})
  RooStats::ProfileLikelihoodTestStat profll(*SBModel->GetPdf());

  // Construct an hypothesis p-value calculator
  RooStats::AsymptoticCalculator asympCalc(*data, *BModel, *SBModel);

  // Configure calculator for a limit (=one-side interval)
  asympCalc.SetOneSided(true);

  // Construct an hypothesis test inverter
  // i.e. a tool that can calculate the POI value for which (in this case) CLS==𝑝(sbModel)/(1 − 𝑝(Model)) takes a certain
  //value. This inversion requires a scan over possible values of 𝜇.
  RooStats::HypoTestInverter inverter(asympCalc);

  // Statistical configuration of hypothesis test inverter
  //inverter.SetConfidenceLevel(0.90);
  inverter.SetConfidenceLevel(0.95);
  //inverter.UseCLs(true);
  bool useCLs = true;
  inverter.UseCLs(useCLs);
  if (useCLs) { profll.SetOneSided(true); }

  // Technical configuration of hypothesis test inverter
  inverter.SetVerbose(false);
  int npoints = 50;  // number of points to scan
  double poimin = poi->getMin();
  double poimax = poi->getMax();
  inverter.SetFixedScan(npoints, 0.0, poimax);  // set number of points, xmin, xmax

  // Perform calculation of limit
  RooStats::HypoTestInverterResult* result = inverter.GetInterval();

  // Print observed limit
  cout << 100*inverter.ConfidenceLevel() << "% upper limit : " << result->UpperLimit() << endl;
  cout << "Expected upper limits, using the B (alternate) model : " << endl;
  cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << endl;
  cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << endl;
  cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << endl;
  cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << endl;
  cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << endl;

    // plot result of the scan 
  RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", result);
  //TCanvas* c2 = new TCanvas("HypoTestInverter Scan"); 
  TCanvas *cR3 = new TCanvas("cR3","cR3",800,600);
  cR3->SetLogy(false);
  plot2->Draw("CLb 2CL");
  //cR3->SaveAs("SimpleCLsLimit.pdf");*/



/*
  // Sets up b-only model 
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*pdf);
  b_modelNM.SetParametersOfInterest(*w.var("mu"));
  //b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("Stau"),*w.var("B")));
  //b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("B")));
  if(pdf_type == 1){
    b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu")));
  }
  else{
    b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("B")));
  }
  b_modelNM.SetObservables(*w.var("x"));  //This is for 1D limit
  w.var("mu")->setVal(0.0);
  //w.var("mu")->setVal(1.5e-5);
  //w.var("mu")->setMin(0);
  //w.var("mu")->setMax(6.0e-7);
  b_modelNM.SetSnapshot(*w.var("mu"));     // sets up b hypothesis as s = 0

  
   
  // create the alternate (s+b) ModelConfig with given value of s
  double s = 1;
  RooStats::ModelConfig sb_modelNM("S+B_modelNM", &w);
  sb_modelNM.SetPdf(*pdf);
  sb_modelNM.SetObservables(*w.var("x"));  //This is for 1D limit
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

  int npoints = 50;  // number of points to scan
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

  // plot result of the scan 
  RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
  //TCanvas* c2 = new TCanvas("HypoTestInverter Scan"); 
  TCanvas *cR3 = new TCanvas("cR3","cR3",800,600);
  cR3->SetLogy(false);
  plot2->Draw("CLb 2CL");
  //cR3->SaveAs("~/tau_lalpha/eChannel/GetLimit/PlotsCLS/SimpleCLsLimitMmax_m16.pdf");

  Double_t limitval = r->GetExpectedUpperLimit(0);*/
  return 0;
  //return limitval;
  

}


/*Double_t getLimit2D(Double_t lum, Int_t pdf_type)
{

  Int_t Ns = 0;
  Int_t NhBkgStagD = 0;
  Int_t Nb = 0;
  //Int_t Ntau =0;
  //Int_t Nbsm =0;

  RooRealVar *x;
  RooRealVar *y;

  RooDataHist *dS;
  //RooDataHist *dBSM;
  RooDataHist *dSMhBkgStagD;
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
    
  
    //TFile *fSDMpinu_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_DM_pinu_2D.root");
    //TFile *fSDMtau_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hS_DM_tau_2D.root");
    //TFile *fB_2D = new TFile("/home/belle2/johancol/tau_lalpha/eChannel/RootFiles/hB_2D.root");

    
    //----------------------------------------------------
    //                         PDF
    //----------------------------------------------------

   
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_2D_m00_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_2D_m10_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_2D_m20_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_2D_m30_pdf.root");
    TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_2D_m40_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_D_2D_m50_pdf.root");
    

  
    TFile *fhBkgStagD_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_D_2D_pdf.root");
    //TFile *fBkgStagDstarGamma_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hS_DM_tau_2D_pdf.root");
    TFile *fB_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hBkgStag_Dstar_gamma_2D_pdf.root");
    
    

    //hS_D2Dpdf = (TH2D*)fS_2DPdf->Get("hS_D_2D_m00");
    //hS_D2Dpdf = (TH2D*)fS_2DPdf->Get("hS_D_2D_m10");
    //hS_D2Dpdf = (TH2D*)fS_2DPdf->Get("hS_D_2D_m20");
    //hS_D2Dpdf = (TH2D*)fS_2DPdf->Get("hS_D_2D_m30");
    hS_D2Dpdf = (TH2D*)fS_2DPdf->Get("hS_D_2D_m40");
    //hS_D2Dpdf = (TH2D*)fS_2DPdf->Get("hS_D_2D_m50");
   
    
    hhBkgStagD2Dpdf = (TH2D*)fhBkgStagD_2DPdf->Get("hBkgStag_D_2D");
    //hS_DMTau2Dpdf = (TH2D*)fSDMtau_2DPdf->Get("hS_DM_tau_2D_1");
    hB2Dpdf = (TH2D*)fB_2DPdf->Get("hBkgStag_Dstar_gamma_2D");
    
      
      
    //hS_D2D = (TH2D*)fSDm00_2D->Get("hS_D_2D_m00_2");
    //hS_DMPinu2D = (TH2D*)fSDMpinu_2D->Get("hS_DM_pinu_2D_2");
    //hS_DMTau2D = (TH2D*)fSDMtau_2D->Get("hS_DM_tau_2D_2");
    //hB2D = (TH2D*)fB_2D->Get("hB_2D_2");
    
    //hBSM2Dpdf->Add(hS_DMTau2Dpdf,hB2Dpdf,1,1);
    //hBSM2D->Add(hS_DMTau2D,hB2D,1,1);
    //hS_DM2D->Add(hS_DMPinu2D,hS_DMTau2D,1,1);
    
    Ns = hS_D2Dpdf->GetEntries();
    NhBkgStagD = hhBkgStagD2Dpdf->GetEntries();
    //Ntau = hS_DMTau2Dpdf->GetEntries();
    Nb = hB2Dpdf->GetEntries();
    //Nbsm = hBSM2Dpdf->GetEntries();
    
  //cout<<" Ratio of Nbsm/Nenunu :" <<Nb/Nenunu<<endl;
  
  //Int_t NevtsGen = lum*(Nb+((7.78571)*NhBkgStagD))/Luminosity;
  //Normalizing for a luminosity of 3.56249e4/ab (1M events of B->piunu)
  NhBkgStagD = lum*((0.849421)*NhBkgStagD)/Luminosity;
  //NhBkgStagD = lum*NhBkgStagD/Luminosity;
 
  Nb = lum*Nb/Luminosity;
  Ns = lum*Ns/Luminosity;
  
    
  //cout<<" Ratio of Nb/Nenunu :" <<Nb/Nenunu<<endl;
  //cout<<"  NevtsGen : "<<NevtsGen<<endl;
  cout<<"  Ns : "<<Ns<<endl;
  cout<<"  NhBkgStagD : "<<NhBkgStagD<<endl;
  cout<<"  Npinunu : "<<Nb<<endl;

 
    
  x = new RooRealVar("x", "x", M_minLow, M_minUp);
  y = new RooRealVar("y", "y", M_maxLow, M_maxUp);
  x->setBins(nbinsL);
  y->setBins(nbinsH);
  
  dS = new RooDataHist("dS", "dS", RooArgList(*x,*y), Import(*hS_D2Dpdf));
  dSMhBkgStagD = new RooDataHist("dSMhBkgStagD", "dSMhBkgStagD", RooArgList(*x,*y), Import(*hhBkgStagD2Dpdf));
  //dSMtau = new RooDataHist("dSMtau", "dSMtau", RooArgList(*x,*y), Import(*hS_DMTau2Dpdf));
  //dBSM = new RooDataHist("dBSM", "dBSM", RooArgList(*x,*y), Import(*hBSM2Dpdf));
  dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hB2Dpdf));
  //dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hMQBBOPdf));
  
  s_pdf = new RooHistPdf("s_pdf","signal pdf",RooArgList(*x,*y),*dS,2);
  pi_taunu_pdf = new RooHistPdf("pi_taunu_pdf","taupinu 2D pdf",RooArgList(*x,*y),*dSMhBkgStagD,2);
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

    hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x,*y),NhBkgStagD);
    hist_dataB = b_pdf->generateBinned(RooArgList(*x,*y),1);  


  }
  else{
    hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x,*y),NhBkgStagD);
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
  // RooFitResult *fitmass = sb_pdf.fitTo(*hist_data, Extended(),Save(ktRUE), Strategy(0), NumCPU(4)); //For small values of signal, the Likelihood shos some bias
  // //sb_pdf.chi2FitTo(*hist_data);  //This is to test a chi2 fit which show smallest bias than Likelihood 

  // TCanvas *cR1 = new TCanvas("cR1","cR1",800,600);
  // RooPlot* framer = x->frame() ;
  // hist_data->plotOn(framer) ;
  // sb_pdf.plotOn(framer) ;
  // //sb_pdf.plotOn(framer,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  // sb_pdf.plotOn(framer,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed)) ;
  // //pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineStyle(kDashed)) ;
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
  
  Double_t eff_signal = (hS_D2Dpdf->GetEntries())/Ns_generated;
  Double_t eff_taupinu = (hhBkgStagD2Dpdf->GetEntries())/NB_generated;
  //Double_t eff_taupinu = (hhBkgStagD2Dpdf->GetEntries()*(7.78571))/NB_generated;
  //Double_t eff_tau = (hS_DMTau2Dpdf->GetEntries())/Ntau_generated;

  //Double_t eff_bkg = (hxminQBPdf->GetEntries())/Ntau_generated;
  Double_t effrelVal = eff_signal/(eff_taupinu);

  //cout<<" --> Ns : "<<hS_D2Dpdf->GetEntries()<<"  Ntaupinu : "<<hhBkgStagD2Dpdf->GetEntries()<<"  Nb : "<<hB2Dpdf->GetEntries()<<endl;
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
 
      w.var("Snu")->setVal(NhBkgStagD);
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

        
        
      w.var("Snu")->setVal(NhBkgStagD);
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

        
        
      w.var("Snu")->setVal(NhBkgStagD);
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
 
      w.var("Snu")->setVal(NhBkgStagD);
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

        
        
      w.var("Snu")->setVal(NhBkgStagD);
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

        
        
      w.var("Snu")->setVal(NhBkgStagD);
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
 
      w.var("Snu")->setVal(NhBkgStagD);
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

        
        
      w.var("Snu")->setVal(NhBkgStagD);
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

        
        
      w.var("Snu")->setVal(NhBkgStagD);
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


   
  pdf->fitTo(*hist_data, Extended()) ;

  //Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  //w.var("mu")->setMax(100*mu_val);
 
  if(pdf_type == 2 && pdf_type == 3){

    TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
    RooPlot* frame = w.var("x")->frame() ;
    hist_data->plotOn(frame) ;
    pdf->plotOn(frame) ;
    pdf->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
    pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
    pdf->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
    frame->Draw() ;
    cR11->Draw();

  

  }
  else{
  
    TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
    RooPlot* frame = w.var("x")->frame() ;
    hist_data->plotOn(frame) ;
    pdf->plotOn(frame) ;
    pdf->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen)) ;
    //pdf->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
    pdf->plotOn(frame,RooFit::Components("pi_taunu_pdf"),RooFit::LineStyle(kDashed)) ;
    frame->Draw() ;
    cR11->Draw();

  }

  // b-only model construction
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*pdf);
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
  sb_modelNM.SetPdf(*pdf);
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

TH1* normalize(TH1* h) {
    double norm = h->Integral();
    if (norm != 0.0) {
        h->Scale(1.0 / norm);
    } else {
        std::cout << "Norm factor not defined" << std::endl;
    }
    return h;
}

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


void GetHistogram( DataGen3x1 *t, TH1D *hS_D, Int_t type, Int_t sample)
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
	  if(type == 1) hS_D->Fill(Ymax);
	  if(type == 2) hS_D->Fill(Ymin);
	  if(type == 3) hS_D->Fill(xPseudo);
	}	  
      
      if(sample==1)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(!Is_enunu) continue;

	  if(type == 1) hS_D->Fill(Ymax);
	  if(type == 2) hS_D->Fill(Ymin);
	  if(type == 3) hS_D->Fill(xPseudo);
	  
	}

      if(sample==2)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(Is_enunu) continue;

	  if(type == 1) hS_D->Fill(Ymax);
	  if(type == 2) hS_D->Fill(Ymin);
	  if(type == 3) hS_D->Fill(xPseudo);
	  
	}

      if(sample==3)
	{
	  if(type == 1) hS_D->Fill(Ymax);
	  if(type == 2) hS_D->Fill(Ymin);
	  if(type == 3) hS_D->Fill(xPseudo);
	}	  
      
    }

  return;
}


void GetHistogram2D( DataGen3x1 *t, TH2D *hS_D, Int_t sample)
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
	  hS_D->Fill(Ymin,Ymax);
	}	  
      
      if(sample==1)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(!Is_enunu) continue;

	  hS_D->Fill(Ymin,Ymax);
	  
	}

      if(sample==2)
	{
	  Bool_t Is_enunu = false;
	  //Let's check the decays
	  if(t->ndaughtersTauN==3 && (t->daughterTauN[0])*(t->daughterTauN[1])*(t->daughterTauN[2])==-16*11*12 ) Is_enunu = true;
	  if(t->ndaughtersTauP==3 && (t->daughterTauP[0])*(t->daughterTauP[1])*(t->daughterTauP[2])==16*11*12 ) Is_enunu = true;

	  if(Is_enunu) continue;

	  hS_D->Fill(Ymin,Ymax);
	  
	}

      if(sample==3)
	{
	  hS_D->Fill(Ymin,Ymax);
	}	  
      

    }

  return;
}
