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
Double_t getLimit(Double_t lum, Int_t type, Int_t pdf_type);
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
Double_t M_min2Up = 27.0;

Double_t M_max2Low = -1.0;
Double_t M_max2Up = 27.0;

  
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
TH2D *hS2D = new TH2D("hS2D","hS2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hSMPinu2D = new TH2D("hSMPinu2D","hSMPinu2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
//TH2D *hSMTau2D = new TH2D("hSMTau2D","hSMTau2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hB2D = new TH2D("hB2D","hB2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
//TH2D *hBSM2D = new TH2D("hBSM2D","hBSM2D",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp); 


//  For Pdf  ------
TH2D *hS2Dpdf = new TH2D("hS2Dpdf","hS2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hPitaunu2Dpdf = new TH2D("hPitaunu2Dpdf","hPitaunu2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
//TH2D *hSMTau2Dpdf = new TH2D("hSMTau2Dpdf","hSMTau2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hB2Dpdf = new TH2D("hB2Dpdf","hB2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp); 
//TH2D *hBSM2Dpdf = new TH2D("hBSM2Dpdf","hBSM2Dpdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);  




// =============================

// For Mmin^2 limit

// For data ------

TH1D *hSMmin2 = new TH1D("hSMmin2","hSMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hSMPinuMmin2 = new TH1D("hSMPinuMmin2","hSMPinuMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hSMTauMmin2 = new TH1D("hSMTauMmin2","hSMTauMmin2",nbins,M_min2Low,M_min2Up);
TH1D *hBMmin2 = new TH1D("hBMmin2","hBMmin2",nbins,M_min2Low,M_min2Up);

TH1D *hBSMMmin2 = new TH1D("hBSMMmin2","hBSMMmin2",nbins,M_min2Low,M_min2Up);

// For pdf   ------
TH1D *hSMmin2pdf = new TH1D("hSMmin2pdf","hSMmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hPitaunuMmin2pdf = new TH1D("hPitaunuMmin2pdf","hPitaunuMmin2pdf",nbins,M_min2Low,M_min2Up);
//TH1D *hSMTauMmin2pdf = new TH1D("hSMTauMmin2pdf","hSMTauMmin2pdf",nbins,M_min2Low,M_min2Up);
TH1D *hBMmin2pdf = new TH1D("hBMmin2pdf","hBMmin2pdf",nbins,M_min2Low,M_min2Up);

//TH1D *hBSMMmin2pdf = new TH1D("hBSMMmin2pdf","hBSMMmin2pdf",nbins,M_min2Low,M_min2Up);
//TH1D *hxminQBOPdf = new TH1D("hxminQBOPdf","hxminQBOPdf",nbinsL,M_minLow,M_minUp);
//TH1D *hxminQBBOPdf = new TH1D("hxminQBBOPdf","hxminQBBOPdf",nbinsL,M_minLow,M_minUp);

// =============================

// For Mmax^2 limit

// For data ------


TH1D *hSMmax2 = new TH1D("hSMmax2","hSMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hSMPinuMmax2 = new TH1D("hSMPinuMmax2","hSMPinuMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hSMTauMmax2 = new TH1D("hSMTauMmax2","hSMTauMmax2",nbins,M_max2Low,M_max2Up);
TH1D *hBMmax2 = new TH1D("hBMmax2","hBMmax2",nbins,M_max2Low,M_max2Up);

TH1D *hBSMMmax2 = new TH1D("hBSMMmax2","hBSMMmax2",nbins,M_max2Low,M_max2Up);

// For pdf ------

TH1D *hSMmax2pdf = new TH1D("hSMmax2pdf","hSMmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hPitaunuMmax2pdf = new TH1D("hPitaunuMmax2pdf","hPitaunuMmax2pdf",nbins,M_max2Low,M_max2Up);
//TH1D *hSMTauMmax2pdf = new TH1D("hSMTauMmax2pdf","hSMTauMmax2pdf",nbins,M_max2Low,M_max2Up);
TH1D *hBMmax2pdf = new TH1D("hBMmax2pdf","hBMmax2pdf",nbins,M_max2Low,M_max2Up);
//TH1D *hxmaxQBOPdf = new TH1D("hxmaxQBOPdf","hxmaxQBOPdf",nbinsH,M_maxLow,M_maxUp);
//TH1D *hBSMMmax2pdf = new TH1D("hBSMMmax2pdf","hBSMMmax2pdf",nbins,M_max2Low,M_max2Up);



//Double_t Br_tau_to_3pi_X = 0.152;
Double_t Br_B_to_l_nu_D0 = 2*0.0235*0.0395;
//Double_t Br_tau_to_h_X = 0.1203;
Double_t Br_B_to_tau_pinu = 0.000109*0.1091;
Double_t Br_B_to_pinunu = 0.000014;
Double_t Luminosity = 4.19402e7;//3.56249e7; //in fb^-1 based on B->taunu(pinu) -> 1M events = 4.19402e7/fb
Double_t CS_ee_BB = 0.54e6; //in fb
Double_t Ns_generated = 1000000;//*(Br_B_to_l_nu_D0);
Double_t NB_generated = 2*Luminosity*CS_ee_BB*Br_B_to_l_nu_D0*Br_B_to_tau_pinu;



void FindLimitCls_3x1_final_v3()
{
  SetData();
  Double_t limit[3];
  //Double_t lum = 10000; 
  //Double_t lum = 30000; 
  Double_t lum = 50000; 
  

  // Depending on the kind of pdf we can estimnate the UL in normalized POI, Number of
  // Events, or directly the branching ratio, the pdfs used are:

  //   PDF1 = Snu*mu*effrel*S_pdf + Snu*pi_taunu_pdf   <- without background
  //   PDF2 = S*S_pdf + Snu*pi_taunu_pdf + B*pi_nunu_pdf  <- Fitting number of events
  //   PDF3 = fac*mu*S_pdf + Snu*pi_taunu_pdf + B*pi_nunu_pdf   <- Branching ratio as POI

  //Int_t pdf_type =1;
  //Int_t pdf_type =2;
  Int_t pdf_type =3;


  // 1:Mmin2, 2:Mmax2  
  //limit[0] = getLimit(lum,1,pdf_type);
  //limit[1] = getLimit(lum,2,pdf_type);
  limit[2] = getLimit2D(lum,pdf_type);

  cout<<" **************************** "<<endl<<endl;
  //cout<<" Limit : "<<limit[0]<<endl;
  //cout<<" Limit : "<<limit[1]<<endl;
  cout<<" Limit : "<<limit[2]<<endl;

  
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

Double_t getLimit(Double_t lum, Int_t type, Int_t pdf_type)
{

  Int_t Ns = 0;
  Int_t Npitaunu = 0;
  Int_t Nb = 0;
  //Int_t Ntau = 0;
  //Int_t Nbsm = 0;

  RooRealVar *x;

  RooDataHist *dS;
  //RooDataHist *dSM;
  RooDataHist *dSMpitaunu;
  //RooDataHist *dSMtau;
  //RooDataHist *dBSM;
  RooDataHist *dB;
    
  
  RooDataHist *hist_dataSM;
  RooDataHist *hist_dataB;
  //RooDataHist *hist_data;
    
    

  RooHistPdf *s_pdf;
  //RooHistPdf *sm_pdf;
  RooHistPdf *pi_taunu_pdf;
  //RooHistPdf *sm_tau_pdf;
  RooHistPdf *b_pdf;
  
      //-----------------------------------------------
    //                              PDF
    //---------------------------------------------
    TFile *fSMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m00_Mmin2_pdf.root");
    TFile *fSMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m00_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m10_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m10_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m20_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m20_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m30_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m30_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m40_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m40_Mmax2_pdf.root");
    //TFile *fSMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m50_Mmin2_pdf.root");
    //TFile *fSMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_m50_Mmax2_pdf.root");
    
    
    

    TFile *fPitaunuMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hpi_taunu_Mmin2_pdf.root");
    TFile *fPitaunuMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hpi_taunu_Mmax2_pdf.root");
    TFile *fPinunuMmin2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hpi_nunu_Mmin2_pdf.root");
    TFile *fPinunuMmax2pdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hpi_nunu_Mmax2_pdf.root");
    //TFile *fBMmin2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_Mmin2_pdf.root");
    //TFile *fBMmax2pdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_Mmax2_pdf.root");

    
    //----------------------------------------------------
    //                         DATA
    //----------------------------------------------------

    
    //TFile *fSMPinuMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_pinu_Mmin2.root");
    //TFile *fSMPinuMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_pinu_Mmax2.root");
    //TFile *fSMTauMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_Mmin2.root");
    //TFile *fSMTauMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_Mmax2.root");
    //TFile *fBMmin2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_Mmin2.root");
    //TFile *fBMmax2 = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hB_Mmax2.root");
    
    
    
    
  if(type==1)
    {
      /*GetHistogram(tSpdf, hxminPdf, 2, 0);  //2 = Mmin^2
      GetHistogram(tSMpdf, hxminBPdf, 2, 1);  
      GetHistogram(tSMpdf, hxminBOPdf, 2, 2); 
      GetHistogram(tQQpdf, hxminQPdf, 2, 3);  
      
      //GetHistogram(tS, hxmin, 2, 0);  //2 = Mmin^2
      GetHistogram(tSM, hxminB, 2, 1);  
      GetHistogram(tSM, hxminBO, 2, 2); 
      GetHistogram(tQQ, hxminQ, 2, 3);  
      
      hxminQBOPdf->Add(hxminQPdf,hxminBOPdf,1,1);
      hxminQBBOPdf->Add(hxminQBOPdf,hxminBPdf,1,1);
      
      hxminQBO->Add(hxminQ,hxminBO,1,1);
      hxminQBBO->Add(hxminQBO,hxminB,1,1);
      */
      
      
      
      hSMmin2pdf = (TH1D*)fSMmin2pdf->Get("Mmin2");
      hPitaunuMmin2pdf = (TH1D*)fPitaunuMmin2pdf->Get("Mmin2");
      //hSMTauMmin2pdf = (TH1D*)fSMTauMmin2pdf->Get("Mmin2");
      hBMmin2pdf = (TH1D*)fPinunuMmin2pdf->Get("Mmin2");
     
      
      //hSMmin2 = (TH1D*)fSMmin2->Get("Mmin2");
      //hSMPinuMmin2 = (TH1D*)fSMPinuMmin2->Get("Mmin2");
      //hSMTauMmin2 = (TH1D*)fSMTauMmin2->Get("Mmin2");
      //hBMmin2 = (TH1D*)fBMmin2->Get("Mmin2");
      
      
      
      
      //Nsm = Nenunu + Ntau
      Ns = hSMmin2pdf->GetEntries();
      Npitaunu = hPitaunuMmin2pdf->GetEntries();
      //Ntau = hSMTauMmin2pdf->GetEntries();
      Nb = hBMmin2pdf->GetEntries();
      //Nsm = hSMMmin2pdf->GetEntries();
      
      //Int_t c1 =1;
      //Int_t c2 = Nb/Ntau;
      // Add the tau->others + all bkg histograms
      //hBSMMmin2pdf->Add(hSMTauMmin2pdf,hBMmin2pdf,1,1);
      //hBSMMmin2->Add(hSMTauMmin2,hBMmin2,1,1);
	  //Nbsm = hBSMMmin2pdf->GetEntries();
      
      // For Mmin^2
      x = new RooRealVar("x", "x", M_min2Low, M_min2Up);
      x->setBins(nbins);
      
      dS = new RooDataHist("dS", "dS", *x, Import(*hSMmin2pdf));
      dSMpitaunu = new RooDataHist("dSMpitaunu", "dSMpitaunu", *x, Import(*hPitaunuMmin2pdf));
      //dSMtau = new RooDataHist("dSMtau", "dSMtau", *x, Import(*hSMTauMmin2pdf));
      //BSM = new RooDataHist("dBSM", "dBSM", *x, Import(*hBSMMmin2pdf));
      dB = new RooDataHist("dB", "dB", *x, Import(*hBMmin2pdf));
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
      
      
      
      hSMmax2pdf = (TH1D*)fSMmax2pdf->Get("Mmax2");
      hPitaunuMmax2pdf = (TH1D*)fPitaunuMmax2pdf->Get("Mmax2");
      //hSMTauMmax2pdf = (TH1D*)fSMTauMmax2pdf->Get("Mmax2");
      hBMmax2pdf = (TH1D*)fPinunuMmax2pdf->Get("Mmax2");
      
      
      //hSMmax2 = (TH1D*)fSMmax2->Get("Mmax2");
      //hSMPinuMmax2 = (TH1D*)fSMPinuMmax2->Get("Mmax2");
      //hSMTauMmax2 = (TH1D*)fSMTauMmax2->Get("Mmax2");
      //hBMmax2 = (TH1D*)fBMmax2->Get("Mmax2");
      
      
      
      Ns = hSMmax2pdf->GetEntries();
      Npitaunu = hPitaunuMmax2pdf->GetEntries();
      //Ntau = hSMTauMmax2pdf->GetEntries();
      Nb = hBMmax2pdf->GetEntries();
      
      
      //Int_t c1 =1;
      //Int_t c2 = Nb/Ntau;
      // Add the tau->others + all bkg histograms
      //hBSMMmax2pdf->Add(hSMTauMmax2pdf,hBMmax2pdf,1,1);
      //hBSMMmax2->Add(hSMTauMmax2,hBMmax2,1,1);
      
      //Nbsm = hBSMMmax2pdf->GetEntries();
      //hSMMmin2->Add(hSMPinuMmin2,hSMTauMmin2,c1,c2);
      
      // For Mmax^2
      x = new RooRealVar("x", "x", M_max2Low, M_max2Up);
      x->setBins(nbins);
      
      dS = new RooDataHist("dS", "dS", *x, Import(*hSMmax2pdf));
      dSMpitaunu = new RooDataHist("dSMpitaunu", "dSMpitaunu", *x, Import(*hPitaunuMmax2pdf));
      //dSMtau = new RooDataHist("dSMtau", "dSMtau", *x, Import(*hSMTauMmax2pdf));
      //dBSM = new RooDataHist("dBSM", "dBSM", *x, Import(*hBSMMmax2pdf));
      dB = new RooDataHist("dB", "dB", *x, Import(*hBMmax2pdf));
    }

 /*if(type==3)
    {
     
       
      hSPdf = (TH2D*)fSMmin2pdf->Get("Mmin2");
      
      Ns = hXpsPdf->GetEntries();
      Nenunu = hXpsBPdf->GetEntries();
      Nb = hXpsQBOPdf->GetEntries();
      
      
      // For x=2E/mtau
      x = new RooRealVar("x", "x", 0, 2);
      x->setBins(nbins);
      
      dS = new RooDataHist("dS", "dS", *x, Import(*hXpsPdf)); //eaplha
      dSM = new RooDataHist("dSM", "dSM", *x, Import(*hXpsBPdf)); //enunu
      dB = new RooDataHist("dB", "dB", *x, Import(*hXpsQBOPdf)); //sm_bkg - enunu
    }*/

  //Int_t NevtsGen = lum*(Nb+((0.849421)*Npitaunu))/Luminosity;
  //Normalizing for a luminosity of 3.56249e4/ab (1M events of B->piunu)
  Npitaunu = lum*((0.849421)*Npitaunu)/Luminosity;
  //Npitaunu = lum*Npitaunu/Luminosity;
    
  Nb = lum*Nb/Luminosity;
  Ns = lum*Ns/Luminosity;
  //Nbsm = lum*Nbsm/Luminosity;
  
  s_pdf = new RooHistPdf("s_pdf","signal pdf",*x,*dS,2);
  pi_taunu_pdf = new RooHistPdf("pi_taunu_pdf","taupinu pdf",*x,*dSMpitaunu,2);
  //sm_tau_pdf = new RooHistPdf("sm_tau_pdf","tau background pdf",*x,*dSMtau,2);
  //bsm_pdf = new RooHistPdf("bsm_pdf","pinunu",*x,*dBSM,2);
  b_pdf = new RooHistPdf("b_pdf","pinunu pdf",*x,*dB,2);


  
  if(type==1){
      TCanvas *c1 = new TCanvas("c1","c1",800,600);

      //hxminPdf->Draw("HIST");
      hSMmin2pdf->Draw("HIST");
      hPitaunuMmin2pdf->Draw("HISTsame");
      hBMmin2pdf->Draw("HISTsame");
      //hPitaunuMmin2pdf->Draw("HISTsame");
  }
    
  if(type==2){
      TCanvas *c2 = new TCanvas("c2","c2",800,600);

      //hxmaxPdf->Draw("HIST");

      hSMmax2pdf->Draw("HIST");
      hPitaunuMmax2pdf->Draw("HISTsame");
      hBMmax2pdf->Draw("HISTsame");
      //hBMmax2pdf->Draw("HISTsame");
  }
  
 /* if(type==3){
      TCanvas *c3 = new TCanvas("c2","c2",800,600);

      
      hXpsPdf->SetLineColor(4);
      hXpsBOPdf->SetLineColor(1);
      hXpsQPdf->SetLineColor(2);
      
      hXpsPdf->Draw("HIST");
      hXpsBOPdf->Draw("HISTsame");
      hXpsQPdf->Draw("HISTsame");
  }*/  
  
  if(pdf_type == 1){

    hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x),Npitaunu);
    hist_dataB = b_pdf->generateBinned(RooArgList(*x),1);  

  }
  else{
    hist_dataSM = pi_taunu_pdf->generateBinned(RooArgList(*x),Npitaunu);
    hist_dataB = b_pdf->generateBinned(RooArgList(*x),Nb);  

  }
  hist_dataSM->Print();
  hist_dataB->Print();
  //hist_dataSM-Print();
  //hist_dataB.Print();
    
  RooDataHist* hist_data = hist_dataSM;

  //hist_data->add(*hist_dataSM);
  hist_data->add(*hist_dataB);
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

  //w.factory("effenunu[1]");

  Double_t eff_signal = 1;
  Double_t eff_taupinu = 1;
  Double_t eff_tau = 1;
  Double_t eff_sm = 1;
  
  //Lets extract the efficiencies
  
  /*if(type==1)
    {
      eff_signal = (hSMmin2pdf->GetEntries())/Ns_generated;
      eff_taupinu = (hPitaunuMmin2pdf->GetEntries()*(0.849421))/NB_generated;
      //eff_sm = (hSMPinuMmin2pdf->GetEntries())/Ntau_generated;
      //eff_tau = (hSMTauMmin2pdf->GetEntries())/Ntau_generated;
    }
  if(type==2)
    {
      eff_signal = (hSMmax2pdf->GetEntries())/Ns_generated;
      eff_taupinu = (hPitaunuMmax2pdf->GetEntries()*(0.849421))/NB_generated;
      //eff_sm = (hSMMmax2pdf->GetEntries())/Ntau_generated;
      //eff_tau = (hSMTauMmax2pdf->GetEntries())/Ntau_generated;
    }*/
    
  if(type==1)
    {
      eff_signal = (hSMmin2pdf->GetEntries())/Ns_generated;
      eff_taupinu = (hPitaunuMmin2pdf->GetEntries())/NB_generated;
      //eff_sm = (hSMPinuMmin2pdf->GetEntries())/Ntau_generated;
      //eff_tau = (hSMTauMmin2pdf->GetEntries())/Ntau_generated;
    }
  if(type==2)
    {
      eff_signal = (hSMmax2pdf->GetEntries())/Ns_generated;
      eff_taupinu = (hPitaunuMmax2pdf->GetEntries())/NB_generated;
      //eff_sm = (hSMMmax2pdf->GetEntries())/Ntau_generated;
      //eff_tau = (hSMTauMmax2pdf->GetEntries())/Ntau_generated;
    }
    
  /*if(type==3)
    {
      eff_signal = (hXpsPdf->GetEntries())/Ns_generated;
      eff_enunu = (hXpsBPdf->GetEntries())/Ntau_generated;
    }
*/
  
  //Double_t eff_bkg = (hxminQBPdf->GetEntries())/Ntau_generated;


  Double_t effrelVal = eff_signal/(eff_taupinu);
  //Double_t effrelVal = eff_signal/(eff_sm);

  //cout<<" --> Ns : "<<hxminPdf->GetEntries()<<"  Nnu : "<<hxminBPdf->GetEntries()<<"  Nb : "<<hxminQBOPdf->GetEntries()<<endl;
  cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_taupinu<<endl;
  //cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_sm<<endl;
  cout<<" Data : "<<hist_data->numEntries()<<endl;

  cout<<"Ns: "<<Ns<<endl;
  cout<<"Ntaupinu: "<<Npitaunu<<endl;
  cout<<"Npinunu: "<<Nb<<endl;

    // Creating a new variable in order to estimate the branching ratio as a POI
    
  Double_t factor =  2*lum*CS_ee_BB*Br_B_to_l_nu_D0;
  
  cout<<" factor: (2*lum*CS_ee_BB*Br_B_to_l_nu_D0) = "<<factor<<endl;

  if(pdf_type == 1){
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
  
        w.var("Snu")->setVal(Npitaunu);
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


          w.var("mu")->setVal(2.5e-4);
          //w.var("mu")->setMin(0);
  
  
          
        
          //w.var("mu")->setMax(2.5e-1);
          // to m_x = 5.0
          w.var("mu")->setMax(5e-2);
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

  
         // w.var("B")->setVal(Nb);
          //w.var("B")->setMin(0);

          
          
          w.var("Snu")->setVal(Npitaunu);
          //w.var("Snu")->setMin(0);


          w.var("mu")->setVal(2.5e-7);
          //w.var("mu")->setMin(0);
  
  
          
          
          //w.var("mu")->setMax(2.5e-1);
          // to m_x = 5.0
          w.var("mu")->setMax(5e-2);
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
  
        w.var("Snu")->setVal(Npitaunu);
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
  
        w.var("Snu")->setVal(Npitaunu);
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
   
  w.pdf("model")->fitTo(*hist_data, Extended()) ;
  Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  // w.var("mu")->setMax(100*mu_val);

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
  else
  {
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
  
  // Sets up b-only model 
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*w.pdf("model"));
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

  Double_t limitval = r->GetExpectedUpperLimit(0);
  return limitval;

}


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

   
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_2D_m00_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_2D_m10_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_2D_m20_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_2D_m30_pdf.root");
    TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_2D_m40_pdf.root");
    //TFile *fS_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hS_2D_m50_pdf.root");
    

  
    TFile *fPitaunu_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hpi_taunu_2D_pdf.root");
    //TFile *fpinunu_2DPdf = new TFile("/home/belle2/johancol/tau_lalpha/muChannel/RootFiles/hSM_tau_2D_pdf.root");
    TFile *fB_2DPdf = new TFile("~/BToInv/PiChannel/Analysis/Files/RootFiles/hpi_nunu_2D_pdf.root");
    
    

    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m00");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m10");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m20");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m30");
    hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m40");
    //hS2Dpdf = (TH2D*)fS_2DPdf->Get("hS_2D_m50");
   
    
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
