//////////////////////////////////////////////////////////////////////
//
//
//  Example of how to use the RooTauLeptonInvisible class
//
//  Author: Eduard De La Cruz Burelo, CINVESTAV IPN, Mexico, July 2019
//
//
///////////////////////////////////////////////////////////////////////
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TVector3.h>

#define DataGen_cxx
#ifndef DataGen_h
#include "DataGen.h"
#endif

#define DataGenSM_cxx
#ifndef DataGenSM_h
#include "DataGenSM.h"
#endif

using namespace RooFit;

#include "RooTauLeptonInvisible.h"
//Adding library: make changes accordingly.
#pragma cling load("/home/eduard/work/BelleII/development/work-repo/roofitclass/RooTauLeptonInvisible.so"):

// This is to use belle2 style plots, make changes as needed
#include "/home/eduard/work/BelleII/belle2style/Belle2Labels.h"
#include "/home/eduard/work/BelleII/belle2style/Belle2Style.h"
#include "/home/eduard/work/BelleII/belle2style/Belle2Style.C"
#include "/home/eduard/work/BelleII/belle2style/Belle2Utils.C"
#include "/home/eduard/work/BelleII/belle2style/Belle2Labels.C"

#ifdef __CINT__
  gROOT->LoadMacro("/home/eduard/work/BelleII/belle2style/Belle2Style.C");
  gROOT->LoadMacro("/home/eduard/work/BelleII/belle2style/Belle2Labels.C");
#endif


void MakeFitNSM()
{
  Double_t Ecm = 10.5794;
  Double_t E_tau = Ecm/2.0;
  Double_t m_tau = 1.776;

  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;

  Double_t xmin = 0.0;
  Double_t xmax = 2.0;

  TChain *chData = new TChain("Z0");
  //chData->Add("myoutput1CMS.root/tree");
  chData->Add("myoutputCMS.root/tree");
  TTree *tree = (TTree*) chData;
  DataGen t(tree);
    
  Long64_t nentries = t.fChain->GetEntries();
  cout<<" Entries : "<<nentries<<endl;

  RooRealVar x("x", "x", xmin, xmax, "");
  //RooDataSet data("data","data",RooArgSet(x));
  TLorentzVector taum,el,invi,elS;
  TLorentzVector piAll,pi1,pi2,pi3;
  TVector3 pboost, pboostS, pbunit;

  TH1D *hXps = new TH1D("hXps","dN/dXps",50,xmin,xmax);
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      for(Int_t i=0;i<t.__ncandidates__;i++)
	{
	  elS.SetPxPyPzE(t.tau_e_px_CMS[i],t.tau_e_py_CMS[i],t.tau_e_pz_CMS[i],t.tau_e_E_CMS[i]);
	  el.SetPxPyPzE(t.tau_e_px_CMS[i],t.tau_e_py_CMS[i],t.tau_e_pz_CMS[i],t.tau_e_E_CMS[i]);
	  invi.SetPxPyPzE(t.tau_invi_px_CMS[i],t.tau_invi_py_CMS[i],t.tau_invi_pz_CMS[i],t.tau_invi_E_CMS[i]);
	  taum=el;taum+=invi;
	  //Due to some energy lost even at gen level (ISR?) the energy of the tau is not exactly half of the energy of
	  // the beam, so this correct this
	  Double_t betav = taum.Beta(); 

	  Double_t XvalCM = 2.0*el.E()/m_tau;
	  
	  pi1.SetPxPyPzE(t.tau_pi_0_px_CMS[i],t.tau_pi_0_py_CMS[i],t.tau_pi_0_pz_CMS[i],t.tau_pi_0_E_CMS[i]);
	  pi2.SetPxPyPzE(t.tau_pi_1_px_CMS[i],t.tau_pi_1_py_CMS[i],t.tau_pi_1_pz_CMS[i],t.tau_pi_1_E_CMS[i]);
	  pi3.SetPxPyPzE(t.tau_pi_2_px_CMS[i],t.tau_pi_2_py_CMS[i],t.tau_pi_2_pz_CMS[i],t.tau_pi_2_E_CMS[i]);
	  piAll = pi1; piAll+=pi2; piAll+=pi3;
	  
	  pboostS=betav*piAll.BoostVector().Unit();

	  elS.Boost(pboostS);
	  Double_t XvalS = 2.0*elS.E()/m_tau;
	  if(XvalS<xmin || XvalS>xmax) continue;
	  hXps->Fill(XvalS);
	  //x = XvalS;
	  //data.add(x);

	}

    }
  
  RooDataHist data("data","data",RooArgSet(x),hXps);
  data.Print();
  // Int_t cSizeX = 600;
  // Int_t cSizeY = 400;

  
  // TCanvas *c5 = new TCanvas("c5","c5",cSizeX,cSizeY);
  // c5->cd();

  //return;
  
  //RooRealVar x("x", "x", 0.0,2.5);
  RooRealVar mb("mb","mb",1.4,-1.0,1.7); //Mass of the alpha boson set to 1.0 
  RooTauLeptonInvisible Pdf_lb("Pdf_lb","Pdf 2 body",x,mb);
  //RooTauLeptonInvisible Pdf_Invisible("Pdf_Invisible","Pdf 3 body",x);
  
  //RooRealVar fsm("fsm", "fsm",0.9,0,1);	// fraction of SM set to 0.5	     
  //RooAddPdf Global("Global","Global",RooArgList(Pdf_Invisible,Pdf_lb),RooArgList(fsm));

  //mb.setVal(0.0);
  //fsm.setVal(0.9);
  
  //RooDataHist* data = Global.generateBinned(x,100000);  //Toy data

  RooFitResult *fit = Pdf_lb.fitTo(data, Save(kTRUE), Strategy(0), NumCPU(4));

  fit->Print("v");

  SetBelle2Style();
   
  Int_t cSizeX = 600;
  Int_t cSizeY = 400;
  TCanvas *c6 = new TCanvas("c6","c6",cSizeX,cSizeY);
  c6->cd();

  // RooPlot* frame3 = x.frame(0.,2.0,100);
  // data.plotOn(frame3, Name("x = E/2m"));
  // Pdf_lb.paramOn(frame3,Layout(0.65,0.95,0.95));
  // Pdf_lb.plotOn(frame3, Name("model"),LineStyle(1), LineColor(1),RooFit::Precision(1.E-12));
  // //data.plotOn(frame3, Name("x = E/2m"));
  // frame3->Draw();

  RooPlot* frame = x.frame(xmin,xmax,100);
  data.plotOn(frame, Name("data"));
  Pdf_lb.paramOn(frame,Layout(0.65,0.95,0.95));
  Pdf_lb.plotOn(frame, Name("model"),LineStyle(1), LineColor(4),RooFit::Precision(1.E-9));
  //Pdf_lb.plotOn(frame, Name("model"),LineStyle(1), LineColor(1));
  //data.plotOn(frame, Name("x = E/2m"));
  //frame->GetXaxis()->SetTitleSize(0);
  frame->GetYaxis()->SetRangeUser(0.,25000.);
  frame->GetYaxis()->SetMaxDigits(2);
  frame->GetYaxis()->SetTitle("#frac{dN}{dx} [N/0.4]");
  frame->GetXaxis()->SetTitle("x [GeV]");
  frame->Draw();

  //BELLE2Label(0.2,0.8);
  myText(0.7,0.8,1,"MC @ gen level");
  // //myText(0.7,0.7,1,"CM frame");
  myText(0.7,0.7,1,"#frac{dN_{x_{ps}}}{dx_{ps}}");



}

void MakeFitSM()
{
  Double_t Ecm = 10.5794;
  Double_t E_tau = Ecm/2.0;
  Double_t m_tau = 1.776;

  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;

  Double_t xmin = 0.0;
  Double_t xmax = 2.0;
  TChain *chData = new TChain("Z0");
  chData->Add("myoutputCMS_SM_1.root/tree");
  TTree *tree = (TTree*) chData;
  DataGenSM t(tree);
    
  Long64_t nentries = t.fChain->GetEntries();
  cout<<" Entries : "<<nentries<<endl;

  RooRealVar x("x", "x", xmin, xmax, "");
  //RooDataSet data("data","data",RooArgSet(x));
  //RooDataHist data("data","data",RooArgSet(x));
  TLorentzVector taum,el,nue,nut,elS;
  TLorentzVector piAll,pi1,pi2,pi3;
  TVector3 pboost, pboostS, pbunit;

  TH1D *hXps = new TH1D("hXps","dN/dXps",50,0.0,2.0);
  
  nentries = 100000;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      for(Int_t i=0;i<t.__ncandidates__;i++)
	{

	  elS.SetPxPyPzE(t.tau_e_px_CMS[i],t.tau_e_py_CMS[i],t.tau_e_pz_CMS[i],t.tau_e_E_CMS[i]);
	  el.SetPxPyPzE(t.tau_e_px_CMS[i],t.tau_e_py_CMS[i],t.tau_e_pz_CMS[i],t.tau_e_E_CMS[i]);
	  nut.SetPxPyPzE(t.tau_nu_tau_px_CMS[i],t.tau_nu_tau_py_CMS[i],t.tau_nu_tau_pz_CMS[i],t.tau_nu_tau_E_CMS[i]);
	  nue.SetPxPyPzE(t.tau_nu_e_px_CMS[i],t.tau_nu_e_py_CMS[i],t.tau_nu_e_pz_CMS[i],t.tau_nu_e_E_CMS[i]);
	  taum=el;taum+=nut;taum+=nue;
	  //Due to some energy lost even at gen level (ISR?) the energy of the tau is not exactly half of the energy of
	  // the beam, so this correct this
	  Double_t betav = taum.Beta(); 

	  Double_t XvalCM = 2.0*el.E()/m_tau;
	  
	  pi1.SetPxPyPzE(t.tau_pi_0_px_CMS[i],t.tau_pi_0_py_CMS[i],t.tau_pi_0_pz_CMS[i],t.tau_pi_0_E_CMS[i]);
	  pi2.SetPxPyPzE(t.tau_pi_1_px_CMS[i],t.tau_pi_1_py_CMS[i],t.tau_pi_1_pz_CMS[i],t.tau_pi_1_E_CMS[i]);
	  pi3.SetPxPyPzE(t.tau_pi_2_px_CMS[i],t.tau_pi_2_py_CMS[i],t.tau_pi_2_pz_CMS[i],t.tau_pi_2_E_CMS[i]);
	  piAll = pi1; piAll+=pi2; piAll+=pi3;
	  
	  pboostS=betav*piAll.BoostVector().Unit();

	  elS.Boost(pboostS);
	  Double_t XvalS = 2.0*elS.E()/m_tau;
	  if(XvalS<xmin || XvalS>xmax) continue;
	  //x = XvalS;	  
	  //data.add(x);
	  hXps->Fill(XvalS);
	}

    }
  RooDataHist data("data","data",RooArgSet(x),hXps);
  data.Print();
  // Int_t cSizeX = 600;
  // Int_t cSizeY = 400;

  
  // TCanvas *c5 = new TCanvas("c5","c5",cSizeX,cSizeY);
  // c5->cd();

  //return;
  
  //RooRealVar x("x", "x", 0.0,2.5);
  //RooRealVar mb("mb","mb",1.4,-1.0,1.7); //Mass of the alpha boson set to 1.0 
  //RooTauLeptonInvisible Pdf_lb("Pdf_lb","Pdf 2 body",x,mb);
  RooTauLeptonInvisible Pdf_lb("Pdf_lb","Pdf 3 body",x);
  
  //RooRealVar fsm("fsm", "fsm",0.9,0,1);	// fraction of SM set to 0.5	     
  //RooAddPdf Global("Global","Global",RooArgList(Pdf_Invisible,Pdf_lb),RooArgList(fsm));

  //mb.setVal(0.0);
  //fsm.setVal(0.9);
  
  //RooDataHist* data = Global.generateBinned(x,100000);  //Toy data

  //RooFitResult *fit = Pdf_lb.fitTo(data, Save(kTRUE), Strategy(0), NumCPU(4));

  //fit->Print("v");

  SetBelle2Style();
  
  Int_t cSizeX = 600;
  Int_t cSizeY = 400;
  TCanvas *c6 = new TCanvas("c6","c6",cSizeX,cSizeY);
  c6->cd();

  // TCanvas *c0 = new TCanvas("c0","c0",800,800);
  // c0->SetTitle("Massust");
  // TPad *p1 = new TPad("p1","p1",0.0,0.22,1.,1.);
  // p1->Draw();
  
  // TPad *p2 = new TPad("p2","p2",0.0,0.,1.,0.3);
  // p2->Draw();

  // p1->cd();  
  RooPlot* frame = x.frame(xmin,xmax,100);
  data.plotOn(frame, Name("data"));
  Pdf_lb.paramOn(frame,Layout(0.65,0.95,0.95));
  Pdf_lb.plotOn(frame, Name("model"),LineStyle(1), LineColor(4),RooFit::Precision(1.E-9));
  //Pdf_lb.plotOn(frame, Name("model"),LineStyle(1), LineColor(1));
  //data.plotOn(frame, Name("x = E/2m"));
  //frame->GetXaxis()->SetTitleSize(0);
  frame->GetYaxis()->SetRangeUser(0.,25000.);
  frame->GetYaxis()->SetMaxDigits(2);
  frame->GetYaxis()->SetTitle("#frac{dN}{dx} [N/0.4]");
  frame->GetXaxis()->SetTitle("x [GeV]");
  frame->Draw();

  //BELLE2Label(0.2,0.8);
  myText(0.7,0.8,1,"MC @ gen level");
  // //myText(0.7,0.7,1,"CM frame");
  myText(0.7,0.7,1,"#frac{dN_{x_{ps}}}{dx_{ps}}");

  //p2->cd();

  // Residual (Histogram, curve, pull=true)
  // RooHist* hres = frame->residHist("data", "model", true);
  // hres->SetFillColor(kGray+2);
  // RooPlot* ResFrame = x.frame(Title("Residual"));
  // ResFrame->addPlotable(hres,"B");
  
  // TAxis* X = ResFrame->GetXaxis();
  // TAxis* Y = ResFrame->GetYaxis();
  // X->SetLabelSize(0.13);
  // X->SetTitleSize(0.13);
  // X->SetTitleOffset(1.2);
  // Y->SetNdivisions(5,5,0);
  // Y->SetTickSize(0.02);
  // Y->SetLabelSize(0.13);
  // Y->SetTitleSize(0.13);
  // Y->SetTitleOffset(0.4);
  // Y->SetTitle("Pull");
  
  // p2->SetBottomMargin(0.4);
  
  // ResFrame->Draw();


}

void MakeFit()
{
  MakeFitNSM();
  //MakeFitSM();
}
