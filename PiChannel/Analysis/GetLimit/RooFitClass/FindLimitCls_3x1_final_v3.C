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

using namespace RooFit;

#include "RooTauLeptonInvisible.h"
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
Double_t getLimit(Double_t lum, Int_t type);
Double_t getLimit2D(Double_t lum);
TVector3 calculateThrust(const std::vector<TVector3> momenta);
Double_t get_X_PseudoRestFrame(TLorentzVector piAll, TLorentzVector el);
void GetHistogram( DataGen3x1 *t, TH1D *hS, Int_t type, Int_t sample);
void GetHistogram2D( DataGen3x1 *t, TH2D *hS, Int_t sample);

//========================================
//  Histograms definition
//  Q : qqbar bkg
//  B : tau+tau- bkg
//  QB: Q+B
  
Double_t M_minLow = -6.0;
Double_t M_minUp = 4.0;

Double_t M_maxLow = -1.0;
Double_t M_maxUp = 5.5;

Int_t nbins = 100;
Int_t nbinsH = 100;
Int_t nbinsL = 100;


// =============================
// For 2D limit

//  For data ------
TH2D *hM = new TH2D("hM","hM",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMB = new TH2D("hMB","hMB",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMBO = new TH2D("hMBO","hMBO",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMQ = new TH2D("hMQ","hMQ",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMQBO = new TH2D("hMQBO","hMQBO",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMQBBO = new TH2D("hMQBBO","hMQBBO",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);


//  For Pdf  ------
TH2D *hMPdf = new TH2D("hMPdf","hMPdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMBPdf = new TH2D("hMBPdf","hMBPdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMBOPdf = new TH2D("hMBOPdf","hMBOPdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMQPdf = new TH2D("hMQPdf","hMQPdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);  
TH2D *hMQBOPdf = new TH2D("hMQBOPdf","hMQBOPdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);
TH2D *hMQBBOPdf = new TH2D("hMQBBOPdf","hMQBBOPdf",nbinsL,M_minLow,M_minUp,nbinsH,M_maxLow,M_maxUp);

// =============================

// For Mmin^2 limit

// For data ------

TH1D *hxmin = new TH1D("hxmin","hxmin",nbinsL,M_minLow,M_minUp);
TH1D *hxminB = new TH1D("hxminB","hxminB",nbinsL,M_minLow,M_minUp);
TH1D *hxminBO = new TH1D("hxminBO","hxminBO",nbinsL,M_minLow,M_minUp);
TH1D *hxminQ = new TH1D("hxminQ","hxminQ",nbinsL,M_minLow,M_minUp);
TH1D *hxminQBO = new TH1D("hxminQBO","hxminQBO",nbinsL,M_minLow,M_minUp);
TH1D *hxminQBBO = new TH1D("hxminT","hxminT",nbinsL,M_minLow,M_minUp);

// For pdf   ------
TH1D *hxminPdf = new TH1D("hxminPdf","hxminPdf",nbinsL,M_minLow,M_minUp);
TH1D *hxminBPdf = new TH1D("hxminBPdf","hxminBPdf",nbinsL,M_minLow,M_minUp);
TH1D *hxminBOPdf = new TH1D("hxminBOPdf","hxminBOPdf",nbinsL,M_minLow,M_minUp);
TH1D *hxminQPdf = new TH1D("hxminQPdf","hxminQPdf",nbinsL,M_minLow,M_minUp);
TH1D *hxminQBOPdf = new TH1D("hxminQBOPdf","hxminQBOPdf",nbinsL,M_minLow,M_minUp);
TH1D *hxminQBBOPdf = new TH1D("hxminQBBOPdf","hxminQBBOPdf",nbinsL,M_minLow,M_minUp);

// =============================

// For Mmax^2 limit

// For data ------


TH1D *hxmax = new TH1D("hxmax","hxmax",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxB = new TH1D("hxmaxB","hxmaxB",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxBO = new TH1D("hxmaxBO","hxmaxBO",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxQ = new TH1D("hxmaxQ","hxmaxQ",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxQBO = new TH1D("hxmaxQBO","hxmaxQBO",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxQBBO = new TH1D("hxmaxQBBO","hxmaxBBO",nbinsH,M_maxLow,M_maxUp);

// For pdf ------

TH1D *hxmaxPdf = new TH1D("hxmaxPdf","hxmaxPdf",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxBPdf = new TH1D("hxmaxBPdf","hxmaxBPdf",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxBOPdf = new TH1D("hxmaxBOPdf","hxmaxBOPdf",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxQPdf = new TH1D("hxmaxQPdf","hxmaxQPdf",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxQBOPdf = new TH1D("hxmaxQBOPdf","hxmaxQBOPdf",nbinsH,M_maxLow,M_maxUp);
TH1D *hxmaxQBBOPdf = new TH1D("hxmaxBQBOPdf","hxmaxQBBOPdf",nbinsH,M_maxLow,M_maxUp);

  
// =============================

// For x=2E/mtau in pseudo-rest-frame limit

//  For data ------

TH1D *hXps = new TH1D("hXps","hXps",nbins,0,2);
TH1D *hXpsB = new TH1D("hXpsB","hXpsB",nbins,0,2);
TH1D *hXpsBO = new TH1D("hXpsBO","hXpsBO",nbins,0,2);
TH1D *hXpsQ = new TH1D("hXpsQ","hXps",nbins,0,2);
TH1D *hXpsQBO = new TH1D("hXpsQBO","hXpsQBO",nbins,0,2);
TH1D *hXpsQBBO = new TH1D("hXpsQBBO","hXpsQBBO",nbins,0,2);

//  For pdf ------

TH1D *hXpsPdf = new TH1D("hXpsPdf","hXpsPdf",nbins,0,2);
TH1D *hXpsBPdf = new TH1D("hXpsBPdf","hXpsBPdf",nbins,0,2);
TH1D *hXpsBOPdf = new TH1D("hXpsBOPdf","hXpsBOPdf",nbins,0,2);
TH1D *hXpsQPdf = new TH1D("hXpsQPdf","hXpsPdf",nbins,0,2);
TH1D *hXpsQBOPdf = new TH1D("hXpsQBOPdf","hXpsQBOPdf",nbins,0,2);
TH1D *hXpsQBBOPdf = new TH1D("hXpsQBBOPdf","hXpsQBBOPdf",nbins,0,2);


Double_t Br_tau_to_3pi_X = 0.152;
Double_t Br_tau_to_e_X = 0.1782;
Double_t Luminosity = 100; //in fb^-1
Double_t CS_ee_tautau = 0.919e6; //in fb
Double_t Ns_generated = 1000000*Br_tau_to_3pi_X;
Double_t Ntau_generated = 2*Luminosity*CS_ee_tautau*Br_tau_to_3pi_X*Br_tau_to_e_X;


void FindLimitCls_3x1_final_v3()
{
  SetData();
  Double_t limit[4];
  //Double_t lum = 10; //10/fb
  //Double_t lum = 50; // 50/fb
  //Double_t lum = 100; //100/fb
  //Double_t lum = 500; //500/fb
  //Double_t lum = 1000; //1/ab
  Double_t lum = 5000; //5/ab
  //Double_t lum = 10000; //10/ab
  //Double_t lum = 50000; //50/ab
  //limit[0] = getLimit(lum,1);
  //limit[1] = getLimit(lum,2);
  limit[2] = getLimit(lum,3);
  //limit[3] = getLimit2D(lum);

  cout<<" **************************** "<<endl<<endl;
  //cout<<" Limit : "<<limit[0]<<endl;
  //cout<<" Limit : "<<limit[1]<<endl;
  cout<<" Limit : "<<limit[2]<<endl;
  //cout<<" Limit : "<<limit[3]<<endl;
  
  return;
}


void SetData()
{

    //---------------------------------------------------------
  //                           Data
  //---------------------------------------------------------

  chSpdf = new TChain("tree");
  chSMpdf = new TChain("tree");
  chQQpdf = new TChain("tree");

  chS = new TChain("tree");
  chSM = new TChain("tree");
  chQQ = new TChain("tree");

  // For pdf construction
  


  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_0.000000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_0.200000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_0.400000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_0.600000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_0.800000_r1_1000000_pID_15.root/tree");
  chSpdf->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_inv_m_1.000000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.200000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.400000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.600000_r1_1000000_pID_15.root/tree");

  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.000000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.200000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.400000_r1_1000000_pID_15.root/tree");
  //chSpdf->Add("../../Ntuples/3x1/invisible/pairPYTHIA8_3x1_pi_e_inv_m_1.600000_r1_1000000_pID_15.root/tree");

  chSMpdf->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_taupair.root/tree");
  chQQpdf->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_QQbar.root/tree");
  
  /*
  chSMpdf->Add("../../Ntuples/3x1/taupair/pairPYTHIA8_3x1_pi_e_SM_only_r1_4595000_pID_15.root/tree");
  
  chQQpdf->Add("../../Ntuples/3x1/qqbar/pairPYTHIA8_3x1_pi_e_SM_only_r1_2005000_pID_1.root/tree");
  chQQpdf->Add("../../Ntuples/3x1/qqbar/pairPYTHIA8_3x1_pi_e_SM_only_r1_8025000_pID_2.root/tree");
  chQQpdf->Add("../../Ntuples/3x1/qqbar/pairPYTHIA8_3x1_pi_e_SM_only_r1_1915000_pID_3.root/tree");
  chQQpdf->Add("../../Ntuples/3x1/qqbar/pairPYTHIA8_3x1_pi_e_SM_only_r1_6645000_pID_4.root/tree");
  chQQpdf->Add("../../Ntuples/3x1/qqbar/pairPYTHIA8_3x1_pi_e_SM_only_r1_5500000_pID_5.root/tree");
  */
  
  // Data set for limit calculation

  //chS->Add("../../Ntuples/3x1/pairPYTHIA8_3x1_pi_e_inv_m_1.000000_r1_1000000_pID_15.root/tree");
  
  chSM->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_r2_919000_pID_15.root/tree");      
  
  chQQ->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_r2_401000_pID_1.root/tree");
  chQQ->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_r2_1605000_pID_2.root/tree");
  chQQ->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_r2_383000_pID_3.root/tree");
  chQQ->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_r2_1329000_pID_4.root/tree");
  chQQ->Add("/home/belle2/johancol/tau3x1/LimitPaper/pairPYTHIA8_3x1_pi_e_SM_only_r2_1100000_pID_5.root/tree");
  
  
  treeSpdf = (TTree*) chSpdf;
  treeSMpdf = (TTree*) chSMpdf;
  treeQQpdf = (TTree*) chQQpdf;

  treeS = (TTree*) chS;
  treeSM = (TTree*) chSM;
  treeQQ = (TTree*) chQQ;
  
  tSpdf = new DataGen3x1(treeSpdf);
  tSMpdf = new DataGen3x1(treeSMpdf);
  tQQpdf = new DataGen3x1(treeQQpdf);

  tS = new DataGen3x1(treeS);
  tSM = new DataGen3x1(treeSM);
  tQQ = new DataGen3x1(treeQQ);

  return; 

}

Double_t getLimit(Double_t lum, Int_t type)
{

  Int_t Ns = 0;
  Int_t Nenunu = 0;
  Int_t Nb = 0;

  RooRealVar *x;

  RooDataHist *dS;
  RooDataHist *dSM;
  RooDataHist *dB;
    
  RooDataHist *dataS;
  RooDataHist *dataSM;

  RooHistPdf *s_pdf;
  RooHistPdf *sm_pdf;
  RooHistPdf *b_pdf;
  
  if(type==1)
    {
      GetHistogram(tSpdf, hxminPdf, 2, 0);  //2 = Mmin^2
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
      
      Ns = hxminPdf->GetEntries();
      Nenunu = hxminBPdf->GetEntries();
      Nb = hxminQBOPdf->GetEntries();
      
	  
      // For Mmin^2
      x = new RooRealVar("x", "x", M_minLow, M_minUp);
      x->setBins(nbinsL);
      
      dS = new RooDataHist("dS", "dS", *x, Import(*hxminPdf));
      dSM = new RooDataHist("dSM", "dSM", *x, Import(*hxminBPdf));
      dB = new RooDataHist("dB", "dB", *x, Import(*hxminQBOPdf));
    }
  
  if(type==2)
    {
      GetHistogram(tSpdf, hxmaxPdf, 1, 0);  //1 = Mmax^2
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
      
      Ns = hxmaxPdf->GetEntries();
      Nenunu = hxmaxBPdf->GetEntries();
      Nb = hxmaxQBOPdf->GetEntries();
      
      // For Mmax^2
      x = new RooRealVar("x", "x", M_maxLow, M_maxUp);
      x->setBins(nbinsH);
      
      dS = new RooDataHist("dS", "dS", *x, Import(*hxmaxPdf));
      dSM = new RooDataHist("dSM", "dSM", *x, Import(*hxmaxBPdf));
      dB = new RooDataHist("dB", "dB", *x, Import(*hxmaxQBOPdf));
    }

  if(type==3)
    {
      GetHistogram(tSpdf, hXpsPdf, 3, 0);  //2 = Mmin^2
      GetHistogram(tSMpdf, hXpsBPdf, 3, 1);  
      GetHistogram(tSMpdf, hXpsBOPdf, 3, 2); 
      GetHistogram(tQQpdf, hXpsQPdf, 3, 3);  
      
      //GetHistogram(tS, hXps, 3, 0);  //3 = Xps
      GetHistogram(tSM, hXpsB, 3, 1);  
      GetHistogram(tSM, hXpsBO, 3, 2); 
      GetHistogram(tQQ, hXpsQ, 3, 3);  
      
      hXpsQBOPdf->Add(hXpsQPdf,hXpsBOPdf,1,1);
      hXpsQBBOPdf->Add(hXpsQBOPdf,hXpsBPdf,1,1);
      
      hXpsQBO->Add(hXpsQ,hXpsBO,1,1);
      hXpsQBBO->Add(hXpsQBO,hXpsB,1,1);
      
      Ns = hXpsPdf->GetEntries();
      Nenunu = hXpsBPdf->GetEntries();
      Nb = hXpsQBOPdf->GetEntries();
      
      
      // For x=2E/mtau
      x = new RooRealVar("x", "x", 0, 2);
      x->setBins(nbins);
      
      dS = new RooDataHist("dS", "dS", *x, Import(*hXpsPdf)); //eaplha
      dSM = new RooDataHist("dSM", "dSM", *x, Import(*hXpsBPdf)); //enunu
      dB = new RooDataHist("dB", "dB", *x, Import(*hXpsQBOPdf)); //sm_bkg - enunu
    }

  Int_t NevtsGen = lum*(Nb+Nenunu)/Luminosity;
  Nenunu = lum*Nenunu/Luminosity;
  Nb = lum*Nb/Luminosity;
  Ns = lum*Ns/Luminosity;
  
  //s_pdf = new RooHistPdf("s_pdf","signal pdf",*x,*dS,2);
  //sm_pdf = new RooHistPdf("sm_pdf","signal pdf",*x,*dSM,2);
  b_pdf = new RooHistPdf("b_pdf","background pdf",*x,*dB,2);

  RooRealVar NsR("NsR", "fraction signal", Ns, 0, 3*Ns);
  RooRealVar NenunuR("NenunuR", "fraction enunu", Nenunu, 0, 3*NevtsGen);
  //RooRealVar NbR("NbR", "fraction", NevtsGen , 0 , 3*NevtsGen);
  RooRealVar NbR("NbR", "fraction", Nb , 0 , 3*NevtsGen);
      
  //RooAddPdf sb_pdf("sb_pdf", "Signal+Bkg", RooArgList(*s_pdf, *b_pdf), RooArgList(NsR, NbR));
  //RooAddPdf sb_pdf("sb_pdf", "Signal+Bkg", RooArgList(s_pdf, sm_pdf, b_pdf), RooArgList(NsR, NenunuR, NbR));

  //RooAddPdf smb_pdf("smb_pdf", "Signal+Bkg", RooArgList(*sm_pdf, *b_pdf), RooArgList(NenunuR, NbR));

  
  if(type==1){
      TCanvas *c1 = new TCanvas("c1","c1",800,600);

      //hxminPdf->Draw("HIST");
      hxminBPdf->Draw("HIST");
      hxminBOPdf->Draw("HISTsame");
      hxminQPdf->Draw("HISTsame");
  }
    
  if(type==2){
      TCanvas *c2 = new TCanvas("c2","c2",800,600);

      //hxmaxPdf->Draw("HIST");
      hxmaxBPdf->Draw("HIST");
      hxmaxBOPdf->Draw("HISTsame");
      hxmaxQPdf->Draw("HISTsame");
  }
  
  if(type==3){
      TCanvas *c3 = new TCanvas("c3","c3",800,600);

      
      hXpsBPdf->SetLineColor(4);
      hXpsBOPdf->SetLineColor(1);
      hXpsQPdf->SetLineColor(2);
      
      hXpsBPdf->Draw("HIST");
      hXpsBOPdf->Draw("HISTsame");
      hXpsQPdf->Draw("HISTsame");
  }  
  

  
 
  

  

  //return 0;

  //RooRandom::randomGenerator()->SetSeed(1000);
    
  /*RooDataHist* hist_dataSM;
  RooDataHist* hist_dataB;

  hist_dataSM = sm_pdf->generateBinned(RooArgList(*x),Nenunu);
  hist_dataB = b_pdf->generateBinned(RooArgList(*x),Nb);

  hist_dataSM->Print();
  hist_dataB->Print();
  
  RooDataHist* hist_data = hist_dataSM;
  
  //hist_data->add(*hist_dataSM);
  hist_data->add(*hist_dataB);

  hist_data->Print();*/

  // RooDataHist* hist_data;
  
  // hist_data = smb_pdf.generateBinned(RooArgList(*x),NevtsGen);


  //RooDataHist* hist_data;
  //hist_data = b_pdf->generateBinned(RooArgList(*x),NevtsGen);

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
  
  // We are going to import the RooTauLeptonInvisible class with the purpose of create our 2 and 3 body pdfs
  // s_pdf -> s_analytic_pdf : it is the 2 body pdf (tau->e+alpha)
  // sm_pdf -> sm_analytic_pdf : it is the 3 body pdf (tau->e+nu+nu)  
  
  RooRealVar m("m", "m", 0.0, -1.0, 1.7); //Mass of alpha boson set to 1.0  
  RooTauLeptonInvisible s_anali_pdf("Pdf_e_alpha", "2 body Pdf", *x, m);
  RooTauLeptonInvisible sm_anli_pdf("Pdf_e_nunu", "3 body Pdf", *x);
    
    //Let's generate our data
  m.setVal(0.0);
  RooDataHist* hist_DataS_analytic = s_anali_pdf.generateBinned(RooArgList(*x),Ns);
  RooDataHist* hist_DataSM_analytic = sm_anli_pdf.generateBinned(RooArgList(*x),Nenunu);
    
  TH1 *h_DataS_analytic1 = hist_DataS_analytic->createHistogram("h_DataS_analytic", *x);
  TH1 *h_DataSM_analytic1 = hist_DataSM_analytic->createHistogram("h_DataSM_analytic", *x);
    
  h_DataS_analytic1->Smooth();
  h_DataSM_analytic1->Smooth();
    
    
  RooDataHist h_DataS_analytic("h_DataS_analytic", "h_DataS_analytic", RooArgList(*x), h_DataS_analytic1);
  RooDataHist h_DataSM_analytic("h_DataSM_analytic", "h_DataSM_analytic", RooArgList(*x), h_DataSM_analytic1);
    
  h_DataS_analytic.Print();
  h_DataSM_analytic.Print();
  
    
  dataS = new RooDataHist("dS", "dS", *x, Import(*h_DataS_analytic1)); //eaplha
  dataSM = new RooDataHist("dSM", "dSM", *x, Import(*h_DataSM_analytic1)); //eenunu
  
  s_pdf = new RooHistPdf("s_pdf","signal pdf",*x,*dataS,2);
  sm_pdf = new RooHistPdf("sm_pdf","signal pdf",*x,*dataSM,2);
    
    
  //RooDataHist* hist_data_analytic = hist_dataSM;
  
  //hist_data->add(*hist_dataSM);
  //hist_data->add(*hist_dataB);

  //hist_data->Print();
    
  if(type==3){
      TCanvas *c4 = new TCanvas("c4","c4",800,600);    
      
      hXpsPdf->SetLineColor(4);
      hXpsBPdf->SetLineColor(2);
      h_DataS_analytic1->SetLineColor(1);
      h_DataSM_analytic1->SetLineColor(6);
      
      hXpsBPdf->Draw("HIST");
      h_DataSM_analytic1->Draw("HISTsame");
      hXpsPdf->Draw("HISTsame");
      h_DataS_analytic1->Draw("HISTsame");
      
      //hXpsPdf->Draw("HIST");
      //hXpsBOPdf->Draw("HISTsame");
      //h_DataS_analytic1->Draw("HISTsame");
      //h_DataSM_analytic1->Draw("HISTsame");
  }  
    
  RooDataHist* hist_dataSM;
  RooDataHist* hist_dataB;

  hist_dataSM = sm_pdf->generateBinned(RooArgList(*x),Nenunu);
  hist_dataB = b_pdf->generateBinned(RooArgList(*x),Nb);

  hist_dataSM->Print();
  hist_dataB->Print();
  
  RooDataHist* hist_data = hist_dataSM;
  
  //hist_data->add(*hist_dataSM);
  hist_data->add(*hist_dataB);

  hist_data->Print();
    
  
    
    

    
    
  RooWorkspace w("w");

  w.import(*s_pdf);
  w.import(*sm_pdf);
  w.import(*b_pdf);


   
  w.factory("effrel[1]");
  //w.factory("effenunu[1]");

  Double_t eff_signal = 1;
  Double_t eff_enunu = 1;
  
  //Lets extract the efficiencies
  
  if(type==1)
    {
      eff_signal = (hxminPdf->GetEntries())/Ns_generated;
      eff_enunu = (hxminBPdf->GetEntries())/Ntau_generated;
    }
  if(type==2)
    {
      eff_signal = (hxmaxPdf->GetEntries())/Ns_generated;
      eff_enunu = (hxmaxBPdf->GetEntries())/Ntau_generated;
    }
  if(type==3)
    {
      eff_signal = (hXpsPdf->GetEntries())/Ns_generated;
      eff_enunu = (hXpsBPdf->GetEntries())/Ntau_generated;
    }

  
  //Double_t eff_bkg = (hxminQBPdf->GetEntries())/Ntau_generated;


  Double_t effrelVal = eff_signal/eff_enunu;

  //cout<<" --> Ns : "<<hxminPdf->GetEntries()<<"  Nnu : "<<hxminBPdf->GetEntries()<<"  Nb : "<<hxminQBOPdf->GetEntries()<<endl;
  cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_enunu<<endl;
  cout<<" Data : "<<hist_data->numEntries()<<"  vs  "<<5*(Nb+Nenunu)<<endl;

  cout<<"Nb: "<<Nb<<endl;
  cout<<"Nenun: "<<Nenunu<<endl;
  cout<<"Ns: "<<Ns<<endl;

  w.var("effrel")->setVal(effrelVal);
  //w.var("effenunu")->setVal(eff_enunu);

  if(lum==10){

    w.factory("mu[0,1]");
    w.factory("Snu[1e6,0,5.0e+06]");
    w.factory("B[1e5,0,3.0e+06]");
  

  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(8.0e-3);
    
  }

  if(lum==50){

    w.factory("mu[0,1]");
    w.factory("Snu[1e6,0,5.0e+06]");
    w.factory("B[1e5,0,3.0e+06]");

  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(4.0e-3);
    
  }

  if(lum==100){

    w.factory("mu[0,1]");
    w.factory("Snu[1e6,0,5.0e+06]");
    w.factory("B[1e5,0,3.0e+06]");


  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(4.0e-3);
    
  }

  if(lum==500){

    w.factory("mu[0,1]");
    w.factory("Snu[1e7,0,5.0e+07]");
    w.factory("B[1e5,0,3.0e+06]");
 

  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(2.0e-3);
    
  }

  if(lum==1000){

    w.factory("mu[0,1]");
    w.factory("Snu[1e7,0,4.0e+07]");
    w.factory("B[1e6,0,3.0e+06]");
 

  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(1.0e-3);
    
  }

  if(lum==5000){

    w.factory("mu[0,1]");
    w.factory("Snu[1e8,0,3.0e+08]");
    w.factory("B[1e7,0,2.0e+07]");


  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(4.0e-4);
    
  }

  if(lum==10000){

    w.factory("mu[0,1]");
    w.factory("Snu[1e8,0,5.0e+08]");
    w.factory("B[1e7,0,3.0e+08]");
 

  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb-1e6);
    w.var("B")->setMin(0);
  //w.var("B")->setMax(lum*Nb/Luminosity);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu-1e7);
    w.var("Snu")->setMin(0);
  //w.var("Snu")->setMax(2.0*lum*Nenunu/Luminosity);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(4.0e-4);
    
  }


  if(lum==50000){

    w.factory("mu[0,1]");
    w.factory("Snu[0,4.0e+9]");
    w.factory("B[0,2.0e+8]");
 

  //w.var("B")->setVal(NevtsGen);
  //w.var("B")->setVal(lum*Nb/(2.0*Luminosity));
    w.var("B")->setVal(Nb-1e7);
    w.var("B")->setMin(0);
    //w.var("B")->setMax(2.0e+8);
  //w.var("B")->setMax(100000*Nb);
    
  //w.var("B")->setVal(5*hMQBBO->GetEntries());
  //w.var("Snu")->setVal(lum*Nenunu/(2.0*Luminosity));
    w.var("Snu")->setVal(Nenunu-1e8);
    w.var("Snu")->setMin(0);
    //w.var("Snu")->setMax(4.0e+9);
  //w.var("Snu")->setMax(1000000*Nenunu);

    w.var("mu")->setVal(1.0e-5);
  //w.var("mu")->setVal(lum*Nenunu/Luminosity);
    w.var("mu")->setMin(0);
  //w.var("mu")->setMax((4/100)*1/TMath::Sqrt(lum));
  //w.var("mu")->setMax(1.4e-4);
    w.var("mu")->setMax(1.5e-4);
    
  }
  
  //w.factory("expr::Snu('Bnu*Bre*effenunu',Bnu, Bre[0.25,0,1.0],effenunu)") ;

  //w.factory("expr::S('Snu*mu*effrel',Snu, mu[0.0000001,0,0.00003],effrel)") ;
  w.factory("expr::S('Snu*mu*effrel',Snu, mu,effrel)") ;
  w.factory("SUM::model(S*s_pdf,Snu*sm_pdf,B*b_pdf)") ;

  //w.factory("expr::S('B*mu*effrel',B, mu[0.00001,0,0.02],effrel)") ;
  //w.factory("SUM::model(S*s_pdf,B*b_pdf)") ;
   
  w.pdf("model")->fitTo(*hist_data, Extended()) ;
  Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  // w.var("mu")->setMax(100*mu_val);

  TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
  RooPlot* frame = w.var("x")->frame() ;
  hist_data->plotOn(frame) ;
  w.pdf("model")->plotOn(frame) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed)) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  frame->Draw() ;
  cR11->Draw();

  //return 0;
  
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*w.pdf("model"));
  b_modelNM.SetParametersOfInterest(*w.var("mu"));
  b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("B")));
  b_modelNM.SetObservables(*w.var("x"));  //This is for 1D limit

  w.var("mu")->setVal(0.0);
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
  cout << "med. limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
  cout << "med. limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
  cout << endl;

  // plot result of the scan 
  RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
  //TCanvas* c2 = new TCanvas("HypoTestInverter Scan"); 
  TCanvas *cR3 = new TCanvas("cR3","cR3",800,600);
  cR3->SetLogy(false);
  plot2->Draw("CLb 2CL");
  cR3->SaveAs("SimpleCLsLimit.pdf"); 

  Double_t limitval = r->GetExpectedUpperLimit(0);
  return limitval;

}


Double_t getLimit2D(Double_t lum)
{

  Int_t Ns = 0;
  Int_t Nenunu = 0;
  Int_t Nb = 0;

  RooRealVar *x;
  RooRealVar *y;

  RooDataHist *dS;
  RooDataHist *dSM;
  RooDataHist *dB;

  RooHistPdf *s_pdf;
  RooHistPdf *sm_pdf;
  RooHistPdf *b_pdf;

  //============================================
  //   This is for 2D pdf (Mmin^2, Mmax^2)
  //============================================
  
  GetHistogram2D(tSpdf, hMPdf, 0); 
  GetHistogram2D(tSMpdf, hMBPdf, 1);  
  GetHistogram2D(tSMpdf, hMBOPdf, 2); 
  GetHistogram2D(tQQpdf, hMQPdf, 3);  
  
  //GetHistogram2D(tS, hM, 0); 
  GetHistogram2D(tSM, hMB, 1);  
  GetHistogram2D(tSM, hMBO, 2); 
  GetHistogram2D(tQQ, hMQ, 3);  
  
      
  hMQBOPdf->Add(hMQPdf,hMBOPdf,1,1);
  hMQBBOPdf->Add(hMQBOPdf,hMBPdf,1,1);
  
  hMQBO->Add(hMQ,hMBO,1,1);
  hMQBBO->Add(hMQBO,hMB,1,1);
  
  Ns = hMPdf->GetEntries();
  Nenunu = hMBPdf->GetEntries();
  Nb = hMQBOPdf->GetEntries();

  cout<<" Ratio of Nb/Nenunu :" <<Nb/Nenunu<<endl;
  
  Int_t NevtsGen = lum*(Nb+Nenunu)/Luminosity;
  Nenunu = lum*Nenunu/Luminosity;
  Nb = lum*Nb/Luminosity;
  
  cout<<"  NevtsGen : "<<NevtsGen<<endl;
  
  x = new RooRealVar("x", "x", M_minLow, M_minUp);
  y = new RooRealVar("y", "y", M_maxLow, M_maxUp);
  x->setBins(nbinsL);
  y->setBins(nbinsH);
  
  dS = new RooDataHist("dS", "dS", RooArgList(*x,*y), Import(*hMPdf));
  dSM = new RooDataHist("dSM", "dSM", RooArgList(*x,*y), Import(*hMBPdf));
  dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hMQBOPdf));
  //dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hMQBBOPdf));
  
  s_pdf = new RooHistPdf("s_pdf","signal pdf",RooArgList(*x,*y),*dS,2);
  sm_pdf = new RooHistPdf("sm_pdf","signal sm pdf",RooArgList(*x,*y),*dSM,2);
  b_pdf = new RooHistPdf("b_pdf","background pdf",RooArgList(*x,*y),*dB,2);

  RooRealVar NsR("NsR", "fraction signal", Ns, 0, 3*Ns);
  RooRealVar NenunuR("NenunuR", "fraction enunu", Nenunu, 0, 3*NevtsGen);
  RooRealVar NbR("NbR", "fraction", Nb, 0 , 3*NevtsGen);
  //RooRealVar NbR("NbR", "fraction", NevtsGen, 0 , 5*NevtsGen);
      
  //RooAddPdf sb_pdf("sb_pdf", "Signal+Bkg", RooArgList(*s_pdf, *b_pdf), RooArgList(NsR, NbR));
  //RooAddPdf sb_pdf("sb_pdf", "Signal+Bkg", RooArgList(*s_pdf, *sm_pdf, *b_pdf), RooArgList(NsR, NenunuR, NbR));
  RooAddPdf smb_pdf("smb_pdf", "Signal+Bkg", RooArgList(*sm_pdf, *b_pdf), RooArgList(NenunuR, NbR));

  RooRandom::randomGenerator()->SetSeed(1000);
  
  RooDataHist* hist_dataSM;
  RooDataHist* hist_dataB;

  hist_dataSM = sm_pdf->generateBinned(RooArgList(*x,*y),Nenunu);
  hist_dataB = b_pdf->generateBinned(RooArgList(*x,*y),Nb);

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

  w.import(*s_pdf);
  w.import(*sm_pdf);
  w.import(*b_pdf);
   
  w.factory("effrel[1]");
  //w.factory("effenunu[1]");

  //Lets extract the efficiencies
  
  Double_t eff_signal = (hMPdf->GetEntries())/Ns_generated;
  Double_t eff_enunu = (hMBPdf->GetEntries())/Ntau_generated;

  //Double_t eff_bkg = (hxminQBPdf->GetEntries())/Ntau_generated;
  Double_t effrelVal = eff_signal/eff_enunu;

  cout<<" --> Ns : "<<hMPdf->GetEntries()<<"  Nnu : "<<hMBPdf->GetEntries()<<"  Nb : "<<hMQBOPdf->GetEntries()<<endl;
  cout<<" Relative efficiency = "<<effrelVal<<"  from "<<eff_signal<<"   "<<eff_enunu<<endl;
  cout<<" Data : "<<hist_data->numEntries()<<"  vs  "<<5*(Nb+Nenunu)<<endl;

  w.var("effrel")->setVal(effrelVal);
  //w.var("effenunu")->setVal(eff_enunu);

  w.factory("mu[0,1]");
  w.factory("Snu[10000,0,2500000000]");
  w.factory("B[13000,0,250000000]");
  //w.factory("Snu[100000,0,2000000]");
  //w.factory("B[130,0,200000000]");

  w.var("B")->setVal(Nb);
  w.var("B")->setMin(0);
  w.var("B")->setMax(2*Nb);
    
  w.var("Snu")->setVal(Nenunu);
  w.var("Snu")->setMin(0);
  w.var("Snu")->setMax(2.0*Nenunu);

  w.var("mu")->setVal(0.00001);
  w.var("mu")->setMin(0);
  //w.var("mu")->setMax(15*0.001/TMath::Sqrt(lum));
  //w.var("mu")->setMax(5e-5);
  w.var("mu")->setMax(5e-5);
  
  
  //w.factory("expr::Snu('Bnu*Bre*effenunu',Bnu, Bre[0.25,0,1.0],effenunu)") ;

  w.factory("expr::S('Snu*mu*effrel',Snu, mu,effrel)") ;
  //w.factory("expr::S('Snu*mu*effrel',Snu, mu[0.001,0,0.00001],effrel)") ;
  w.factory("SUM::model(S*s_pdf,Snu*sm_pdf,B*b_pdf)") ;
  //w.factory("expr::S('B*mu*effrel',B, mu[0.00001,0,0.006],effrel)") ;
  //w.factory("SUM::model(S*s_pdf,B*b_pdf)") ;
   
  w.pdf("model")->fitTo(*hist_data, Extended()) ;

  //Double_t mu_val = w.var("mu")->getVal();
  //w.var("mu")->setMin(0);
  //w.var("mu")->setMax(100*mu_val);
 
  
  TCanvas *cR11 = new TCanvas("cR11","cR11",800,600);
  RooPlot* frame = w.var("x")->frame() ;
  hist_data->plotOn(frame) ;
  w.pdf("model")->plotOn(frame) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed)) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed)) ;
  w.pdf("model")->plotOn(frame,RooFit::Components("sm_pdf"),RooFit::LineStyle(kDashed)) ;
  frame->Draw() ;
  cR11->Draw();

  //return 0;
  
  RooStats::ModelConfig b_modelNM("b_modelNM", &w);
  b_modelNM.SetPdf(*w.pdf("model"));
  b_modelNM.SetParametersOfInterest(*w.var("mu"));
  b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("Snu"),*w.var("B")));
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
  cout << "med. limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
  cout << "med. limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
  cout << endl;

  RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
  TCanvas *cR3 = new TCanvas("cR3","cR3",800,600);
  cR3->SetLogy(false);
  plot2->Draw("2CL");
  cR3->SaveAs("SimpleCLsLimit.pdf"); 

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
