#include "TROOT.h"
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
#include "TH2D.h"
#include "TRandom.h"
#include <TF1.h>
#include <TFile.h>
#include <iostream>
#include <utility>
#include <RooAbsPdf.h>
#include <TMath.h>
#include <RooFormulaVar.h>
#include <cmath> // Include this for cos and sin functions

using namespace RooFit;
using namespace std; // Added for cout

// Function declarations
void CLsLimit_BToInvCopy();
void loadSignalFiles(const std::string& baseDir, const std::string& mass, const std::string& variable);
//std::pair<TH1D*, TH1D*> ConstructHist1D(TFile*& fS, TFile*& fB, const std::string& variable);
//std::pair<TH2D*, TH2D*> ConstructHist2D(TFile* fS, TFile* fB, const std::string& variable, bool is2D);
void loadAllHistograms();
void SetData(const std::string& mass, const std::string& variable);
void SetInsideEllipseLabel(TH2D hist, RooDataHist dataHist);
Double_t getLimit(Double_t lum, const std::string& variable);
void setupModel(RooDataHist* hist_data, RooHistPdf& s1_pdf, RooHistPdf& s2_pdf, RooHistPdf& s3_pdf, RooHistPdf& b_pdf, Double_t lum, const std::string& variable);

// Histogram limits and bins
Double_t BDTLow = -0.05, BDTUp = 1.05;

Double_t MminLow = -5.0, MminUp = 23.0;
Double_t MmaxLow = -2.0, MmaxUp = 23.0;
Double_t XpsLow = 0.2, XpsUp = 1.2;
Double_t q2Low = -4.0, q2Up = 22;

/*Double_t center_x = -6.649126147447433e-16;
Double_t center_y = -2.5308192913854242e-17;
Double_t width = 3.456188977481231;
Double_t height = 0.1946340391172905;
Double_t angle_degrees = 0.0;*/

/*Double_t center_x = 13.043160254320464;
Double_t center_y =  1.5735658560518488;
Double_t width = 17.849120804919746;
Double_t height = 0.999622685359045;
Double_t angle_degrees = -0.012507487593706729;*/

//Double_t center_x = 8.542308879590326;
//Double_t center_y = 9.98143766843682;
//Double_t width = 27.104025767572782;
//Double_t height = 1.5179346544820214;
//Double_t angle_degrees = 42.55084814068665;



Double_t nbinsX = 150, nbinsY = 150;
Double_t X_min = -5.5, X_max = 33.5;
Double_t Y_min = -0.1, Y_max = 2.6;


//Double_t X_min = -5, X_max = 24;
//Double_t Y_min = -1.2, Y_max = -0.8;

//Double_t X_min = -5, X_max = 23;
//Double_t Y_min = -2, Y_max = 23;

//Double_t XpsLow = -3.0, XpsUp = 22;
//Double_t q2Low = -1.0, q2Up = 22;

Int_t nbins = 250;

// Histograms
// BDT histograms
TH1D *hSD_BDT = new TH1D("hSD_BDT","hSD_BDT",nbins,BDTLow,BDTUp);
TH1D *hSDstarGamma_BDT = new TH1D("hSDstarGamma_BDT","hSDstarGamma_BDT",nbins,BDTLow,BDTUp);
TH1D *hSDstarPi0_BDT = new TH1D("hSDstarPi0_BDT","hSDstarPi0_BDT",nbins,BDTLow,BDTUp);

TH1D *hBkgD_BDT = new TH1D("hBkgD_BDT","hBkgD_BDT",nbins,BDTLow,BDTUp);
TH1D *hBkgDstarGamma_BDT = new TH1D("hBkgDstarGamma_BDT","hBkgDstarGamma_BDT",nbins,BDTLow,BDTUp);
TH1D *hBkgDstarPi0_BDT = new TH1D("hBkgDstarPi0_BDT","hBkgDstarPi0_BDT",nbins,BDTLow,BDTUp);
TH1D *hB1_BDT = new TH1D("hB1_BDT","hB1_BDT",nbins,BDTLow,BDTUp);
TH1D *hB_BDT = new TH1D("hB_BDT","hB_BDT",nbins,BDTLow,BDTUp);

// Xps histograms
TH1D *hSD_Xps = new TH1D("hSD_Xps","hSD_Xps",nbins,XpsLow,XpsUp);
TH1D *hSDstarGamma_Xps= new TH1D("hSDstarGamma_Xps","hSDstarGamma_Xps",nbins,XpsLow,XpsUp);
TH1D *hSDstarPi0_Xps= new TH1D("hSDstarPi0_Xps","hSDstarPi0_Xps",nbins,XpsLow,XpsUp);

TH1D *hBkgD_Xps = new TH1D("hBkgD_Xps","hBkgD_Xps",nbins,XpsLow,XpsUp);
TH1D *hBkgDstarGamma_Xps= new TH1D("hBkgDstarGamma_Xps","hBkgDstarGamma_Xps",nbins,XpsLow,XpsUp);
TH1D *hBkgDstarPi0_Xps = new TH1D("hBkgDstarPi0_Xps","hBkgDstarPi0_Xps",nbins,XpsLow,XpsUp);
TH1D *hB1_Xps = new TH1D("hB1_Xps","hB1_Xps",nbins,XpsLow,XpsUp);
TH1D *hB_Xps = new TH1D("hB_Xps","hB_Xps",nbins,XpsLow,XpsUp);

// q2 histograms
TH1D *hSD_q2 = new TH1D("hSD_q2 ","hSD_q2 ",nbins,q2Low,q2Up);
TH1D *hSDstarGamma_q2= new TH1D("hSDstarGamma_q2","hSDstarGamma_q2",nbins,q2Low,q2Up);
TH1D *hSDstarPi0_q2= new TH1D("hSDstarPi0_q2","hSDstarPi0_q2",nbins,q2Low,q2Up);

TH1D *hBkgD_q2 = new TH1D("hBkgD_q2","hBkgD_q2",nbins,q2Low,q2Up);
TH1D *hBkgDstarGamma_q2= new TH1D("hBkgDstarGamma_q2","hBkgDstarGamma_q2",nbins,q2Low,q2Up);
TH1D *hBkgDstarPi0_q2 = new TH1D("hBkgDstarPi0_q2","hBkgDstarPi0_q2",nbins,q2Low,q2Up);
TH1D *hB1_q2 = new TH1D("hB1_q2","hB1_q2",nbins,q2Low,q2Up);
TH1D *hB_q2 = new TH1D("hB_q2","hB_q2",nbins,q2Low,q2Up);

// Mmin2
TH1D *hSD_Mmin2 = new TH1D("hSD_Mmin2 ","hSD_Mmin2 ",nbins,MminLow,MminUp);
TH1D *hSDstarGamma_Mmin2= new TH1D("hSDstarGamma_Mmin2","hSDstarGamma_Mmin2",nbins,MminLow,MminUp);
TH1D *hSDstarPi0_Mmin2= new TH1D("hSDstarPi0_Mmin2","hSDstarPi0_Mmin2",nbins,MminLow,MminUp);

TH1D *hBkgD_Mmin2 = new TH1D("hBkgD_Mmin2","hBkgD_Mmin2",nbins,MminLow,MminUp);
TH1D *hBkgDstarGamma_Mmin2= new TH1D("hBkgDstarGamma_Mmin2","hBkgDstarGamma_Mmin2",nbins,MminLow,MminUp);
TH1D *hBkgDstarPi0_Mmin2 = new TH1D("hBkgDstarPi0_Mmin2","hBkgDstarPi0_Mmin2",nbins,MminLow,MminUp);
TH1D *hB1_Mmin2 = new TH1D("hB1_Mmin2","hB1_Mmin2",nbins,MminLow,MminUp);
TH1D *hB_Mmin2 = new TH1D("hB_Mmin2","hB_Mmin2",nbins,MminLow,MminUp);

// Mmax2
TH1D *hSD_Mmax2 = new TH1D("hSD_Mmax2 ","hSD_Mmax2 ",nbins,MmaxLow,MmaxUp);
TH1D *hSDstarGamma_Mmax2= new TH1D("hSDstarGamma_Mmax2","hSDstarGamma_Mmax2",nbins,MmaxLow,MmaxUp);
TH1D *hSDstarPi0_Mmax2= new TH1D("hSDstarPi0_Mmax2","hSDstarPi0_Mmax2",nbins,MmaxLow,MmaxUp);

TH1D *hBkgD_Mmax2 = new TH1D("hBkgD_Mmax2","hBkgD_Mmax2",nbins,MmaxLow,MmaxUp);
TH1D *hBkgDstarGamma_Mmax2= new TH1D("hBkgDstarGamma_Mmax2","hBkgDstarGamma_Mmax2",nbins,MmaxLow,MmaxUp);
TH1D *hBkgDstarPi0_Mmax2 = new TH1D("hBkgDstarPi0_Mmax2","hBkgDstarPi0_Mmax2",nbins,MmaxLow,MmaxUp);
TH1D *hB1_Mmax2 = new TH1D("hB1_Mmax2","hB1_Mmax2",nbins,MmaxLow,MmaxUp);
TH1D *hB_Mmax2 = new TH1D("hB_Mmax2","hB_Mmax2",nbins,MmaxLow,MmaxUp);

// 2D histograms
TH2D *hSD_2D = new TH2D("hSD_2D","hSD_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);
TH2D *hSDstarGamma_2D = new TH2D("hSDstar_gamma_2D","hSDstar_gamma_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);
TH2D *hSDstarPi0_2D = new TH2D("hSDstar_pi0_2D","hSDstar_pi0_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);

TH2D *hBkgD_2D = new TH2D("hBkgD_2D","hBkgD_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);
TH2D *hBkgDstarGamma_2D = new TH2D("hBkgDstarGamma_2D","hBkgDstarGamma_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);
TH2D *hBkgDstarPi0_2D = new TH2D("hBkgDstarPi0_2D","hBkgDstarPi0_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);
TH2D *hB1_2D = new TH2D("hB1_2D","hB1_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);
TH2D *hB_2D = new TH2D("hB_2D","hB_2D",nbinsX,X_min,X_max,nbinsY,Y_min,Y_max);

// Declare the TFiles in the global scope
TFile* fSD = nullptr;
TFile* fSDstarGamma = nullptr;
TFile* fSDstarPi0 = nullptr;
TFile* fBkgD = nullptr;
TFile* fBkgDstarGamma = nullptr;
TFile* fBkgDstarPi0 = nullptr;

RooRealVar *x;
RooRealVar *y;

// Creation of the workspace
RooWorkspace w("w");
RooStats::LikelihoodInterval* interval;
RooStats::HypoTestInverterResult* r;
Double_t limitval; 

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
Double_t Ns_generated = 100000;//

Int_t NsD = 0;
Int_t NsDstarGamma = 0;
Int_t NsDstarPi0 = 0;
Int_t Nb = 0;



int main()
{
    CLsLimit_BToInvCopy();
    return 0;
}

void CLsLimit_BToInvCopy()
{
    std::string mass;
    std::cout << "Enter the mass for the signal file: ";
    std::cin >> mass;

     // Get user input for variable
    std::string variable;
    std::cout << "Enter the name of the variable (BDT, Xps, q2, Mmin2, Mmax2, 2D = {Xps, q2}): ";
    std::cin >> variable;
    SetData(mass, variable);
    Double_t limit[4];
    //Double_t lum = 362; 
    //Double_t lum = 403.9; //(361.6 + 42.3 (off 4S)) 
    //Double_t lum = 1000;
    //Double_t lum = 5000;  
    //Double_t lum = 10000; 
    //Double_t lum = 30000; 
    Double_t lum = 50000; 
    
    

    limit[0] = getLimit(lum, variable);

    // 1:BDT, 2:2D,  3:q2,  
    //limit[0] = getLimit(lum);
    //limit[1] = getLimit(lum,2);
    //limit[2] = getLimit(lum,3);
    //limit[3] = getLimit2D(lum);

    cout<<" **************************** "<<endl<<endl;
    cout<<" Limit : "<<limit[0]<<endl;
    //cout<<" Limit : "<<limit[1]<<endl;
    //cout<<" Limit : "<<limit[2]<<endl;
    //cout<<" Limit : "<<limit[3]<<endl;
    //cout<<" Limit : "<<limit[4]<<endl;
}

void loadSignalFiles(const std::string& baseDir, const std::string& mass, const std::string& variable)
{   
    std::string fileNameSD;
    std::string fileNameSDstarGamma;
    std::string fileNameSDstarPi0;
    std::string fileNameBkgD;
    std::string fileNameBkgDstarGamma;
    std::string fileNameBkgDstarPi0;

    if (variable == "BDT")
    {
        fileNameSD = baseDir + "hS_D_m";
        fileNameSDstarGamma = baseDir + "hS_Dstar_gamma_m";
        fileNameSDstarPi0 = baseDir + "hS_Dstar_pi0_m";
        fileNameBkgD = baseDir + "hBkgStag_D_m";
        fileNameBkgDstarGamma = baseDir + "hBkgStag_Dstar_gamma_m";
        fileNameBkgDstarPi0 = baseDir + "hBkgStag_Dstar_pi0_m";
    }
    else if (variable == "2D" || variable == "Xps" || variable == "q2" || variable == "Mmin2" || variable == "Mmax2")
    {
        fileNameSD = baseDir + "hS_D_m";
        fileNameSDstarGamma = baseDir + "hS_Dstar_gamma_m";
        fileNameSDstarPi0 = baseDir + "hS_Dstar_pi0_m";
        fileNameBkgD = baseDir + "hBkgStag_D_";
        fileNameBkgDstarGamma = baseDir + "hBkgStag_Dstar_gamma_";
        fileNameBkgDstarPi0 = baseDir + "hBkgStag_Dstar_pi0_";
    }

    if (variable == "BDT")
    {
        fSD = new TFile((fileNameSD + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
        fSDstarGamma = new TFile((fileNameSDstarGamma + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
        fSDstarPi0 = new TFile((fileNameSDstarPi0 + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());

        fBkgD = new TFile((fileNameBkgD + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
        fBkgDstarGamma = new TFile((fileNameBkgDstarGamma + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
        fBkgDstarPi0 = new TFile((fileNameBkgDstarPi0 + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
    }
    else if ( variable == "Xps" || variable == "q2" || variable == "Mmin2" || variable == "Mmax2")
    {
        fSD = new TFile((fileNameSD + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
        fSDstarGamma = new TFile((fileNameSDstarGamma + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());
        fSDstarPi0 = new TFile((fileNameSDstarPi0 + mass.c_str() + "_" +  variable.c_str() + "_pdf.root").c_str());

        fBkgD = new TFile((fileNameBkgD + variable.c_str() + "_pdf.root").c_str());
        fBkgDstarGamma = new TFile((fileNameBkgDstarGamma + variable.c_str() + "_pdf.root").c_str());
        fBkgDstarPi0 = new TFile((fileNameBkgDstarPi0 + variable.c_str() + "_pdf.root").c_str());

        cout<<fileNameSD + mass.c_str() + "_" +  variable.c_str() + "_pdf.root" <<endl;
    }
    else if (variable == "2D")
    {
        fSD = new TFile((fileNameSD + mass.c_str() + "_" +  variable.c_str() + "_ROT_pdf.root").c_str());
        fSDstarGamma = new TFile((fileNameSDstarGamma + mass.c_str() + "_" +  variable.c_str() + "_ROT_pdf.root").c_str());
        fSDstarPi0 = new TFile((fileNameSDstarPi0 + mass.c_str() + "_" +  variable.c_str() + "_ROT_pdf.root").c_str());

        fBkgD = new TFile((fileNameBkgD + variable.c_str() + "_ROT_pdf.root").c_str());
        fBkgDstarGamma = new TFile((fileNameBkgDstarGamma + variable.c_str() + "_ROT_pdf.root").c_str());
        fBkgDstarPi0 = new TFile((fileNameBkgDstarPi0 + variable.c_str() + "_ROT_pdf.root").c_str());

        cout<<fileNameSD + mass.c_str() + "_" +  variable.c_str() + "_ROT_pdf.root" <<endl;
    }

    
    if (!fSD || !fSDstarGamma|| !fSDstarPi0 || !fBkgD || !fBkgDstarGamma || !fBkgDstarPi0) {
        std::cerr << "Error opening one or more files." << std::endl;
        
        // Handle the error or exit the program
        return;
    }  
}

template<class HISTO>
std::pair<HISTO*, HISTO*> ConstructHist(TFile* fS, TFile* fB, const std::string& variable) {
    HISTO* hS = (HISTO*)fS->Get(variable.c_str());
    if (!hS) {
        // handle error 
        std::cerr << "Error: Unable to retrieve histogram " << variable << " from file S." << std::endl;
        exit(1); // or handle the error in an appropriate way
    }

    HISTO* hB = (HISTO*)fB->Get(variable.c_str());
    if (!hB) {
        // handle error
        std::cerr << "Error: Unable to retrieve histogram " << variable << " from file B." << std::endl;
        exit(1); // or handle the error in an appropriate way
    }

    return std::make_pair(hS, hB);
}


void loadAllHistograms(const std::string& variable)
{
    if (variable == "BDT") {
        auto histPairD_BDT = ConstructHist<TH1D>(fSD, fBkgD, variable);
        hSD_BDT = histPairD_BDT.first;
        hBkgD_BDT = histPairD_BDT.second;

        auto histPairDstarGamma_BDT = ConstructHist<TH1D>(fSDstarGamma, fBkgDstarGamma, variable);
        hSDstarGamma_BDT = histPairDstarGamma_BDT.first;
        hBkgDstarGamma_BDT = histPairDstarGamma_BDT.second;

        auto histPairDstarPi0_BDT = ConstructHist<TH1D>(fSDstarPi0, fBkgDstarPi0, variable);
        hSDstarPi0_BDT = histPairDstarPi0_BDT.first;
        hBkgDstarPi0_BDT = histPairDstarPi0_BDT.second;

        hB1_BDT->Add(hBkgD_BDT, hBkgDstarGamma_BDT, 1, 1);
        hB_BDT->Add(hB1_BDT, hBkgDstarPi0_BDT, 1, 1);
    }
    else if (variable == "Xps") {
        auto histPairD_Xps = ConstructHist<TH1D>(fSD, fBkgD, variable);
        hSD_Xps = histPairD_Xps.first;
        hBkgD_Xps = histPairD_Xps.second;

        auto histPairDstarGamma_Xps = ConstructHist<TH1D>(fSDstarGamma, fBkgDstarGamma, variable);
        hSDstarGamma_Xps = histPairDstarGamma_Xps.first;
        hBkgDstarGamma_Xps = histPairDstarGamma_Xps.second;

        auto histPairDstarPi0_Xps = ConstructHist<TH1D>(fSDstarPi0, fBkgDstarPi0, variable);
        hSDstarPi0_Xps = histPairDstarPi0_Xps.first;
        hBkgDstarPi0_Xps = histPairDstarPi0_Xps.second;

        hB1_Xps->Add(hBkgD_Xps, hBkgDstarGamma_Xps, 1, 1);
        hB_Xps->Add(hB1_Xps, hBkgDstarPi0_Xps, 1, 1);
    }
    else if (variable == "q2") {
        auto histPairD_q2 = ConstructHist<TH1D>(fSD, fBkgD, variable);
        hSD_q2 = histPairD_q2.first;
        hBkgD_q2 = histPairD_q2.second;

        auto histPairDstarGamma_q2 = ConstructHist<TH1D>(fSDstarGamma, fBkgDstarGamma, variable);
        hSDstarGamma_q2 = histPairDstarGamma_q2.first;
        hBkgDstarGamma_q2 = histPairDstarGamma_q2.second;

        auto histPairDstarPi0_q2 = ConstructHist<TH1D>(fSDstarPi0, fBkgDstarPi0, variable);
        hSDstarPi0_q2 = histPairDstarPi0_q2.first;
        hBkgDstarPi0_q2 = histPairDstarPi0_q2.second;

        hB1_q2->Add(hBkgD_q2, hBkgDstarGamma_q2, 1, 1);
        hB_q2->Add(hB1_q2, hBkgDstarPi0_q2, 1, 1);
    }
    else if (variable == "Mmin2") {
        auto histPairD_Mmin2 = ConstructHist<TH1D>(fSD, fBkgD, variable);
        hSD_Mmin2 = histPairD_Mmin2.first;
        hBkgD_Mmin2 = histPairD_Mmin2.second;

        auto histPairDstarGamma_Mmin2 = ConstructHist<TH1D>(fSDstarGamma, fBkgDstarGamma, variable);
        hSDstarGamma_Mmin2 = histPairDstarGamma_Mmin2.first;
        hBkgDstarGamma_Mmin2 = histPairDstarGamma_Mmin2.second;

        auto histPairDstarPi0_Mmin2 = ConstructHist<TH1D>(fSDstarPi0, fBkgDstarPi0, variable);
        hSDstarPi0_Mmin2 = histPairDstarPi0_Mmin2.first;
        hBkgDstarPi0_Mmin2 = histPairDstarPi0_Mmin2.second;

        hB1_Mmin2->Add(hBkgD_Mmin2, hBkgDstarGamma_Mmin2, 1, 1);
        hB_Mmin2->Add(hB1_Mmin2, hBkgDstarPi0_Mmin2, 1, 1);
    }
    else if (variable == "Mmax2") {
        auto histPairD_Mmax2 = ConstructHist<TH1D>(fSD, fBkgD, variable);
        hSD_Mmax2 = histPairD_Mmax2.first;
        hBkgD_Mmax2 = histPairD_Mmax2.second;

        auto histPairDstarGamma_Mmax2 = ConstructHist<TH1D>(fSDstarGamma, fBkgDstarGamma, variable);
        hSDstarGamma_Mmax2 = histPairDstarGamma_Mmax2.first;
        hBkgDstarGamma_Mmax2 = histPairDstarGamma_Mmax2.second;

        auto histPairDstarPi0_Mmax2 = ConstructHist<TH1D>(fSDstarPi0, fBkgDstarPi0, variable);
        hSDstarPi0_Mmax2 = histPairDstarPi0_Mmax2.first;
        hBkgDstarPi0_Mmax2 = histPairDstarPi0_Mmax2.second;

        hB1_Mmax2->Add(hBkgD_Mmax2, hBkgDstarGamma_Mmax2, 1, 1);
        hB_Mmax2->Add(hB1_Mmax2, hBkgDstarPi0_Mmax2, 1, 1);
    }
    else if (variable == "2D") {
        auto histPairD_2D = ConstructHist<TH2D>(fSD, fBkgD, variable);
        hSD_2D = histPairD_2D.first;
        hBkgD_2D = histPairD_2D.second;

        auto histPairDstarGamma_2D = ConstructHist<TH2D>(fSDstarGamma, fBkgDstarGamma, variable);
        hSDstarGamma_2D = histPairDstarGamma_2D.first;
        hBkgDstarGamma_2D = histPairDstarGamma_2D.second;

        auto histPairDstarPi0_2D = ConstructHist<TH2D>(fSDstarPi0, fBkgDstarPi0, variable);
        hSDstarPi0_2D = histPairDstarPi0_2D.first;
        hBkgDstarPi0_2D = histPairDstarPi0_2D.second;

        hB1_2D->Add(hBkgD_2D, hBkgDstarGamma_2D, 1, 1);
        hB_2D->Add(hB1_2D, hBkgDstarPi0_2D, 1, 1);
    }
}

void SetData(const std::string& mass, const std::string& variable)
{
    // Define the base directory for the signal files
    std::string BaseDir = "~/BToInv/KChannel/Analysis/Files/RootFiles/";

    loadSignalFiles(BaseDir, mass, variable);

    // Load all histograms
    loadAllHistograms(variable);
}

void setupModel(RooDataHist* hist_data, RooHistPdf* s1_pdf, RooHistPdf* s2_pdf, RooHistPdf* s3_pdf, RooHistPdf* b_pdf, Double_t lum, const std::string& variable)
{
    Double_t eff_signal_D = 1;
    Double_t eff_signal_DstarGamma = 1;
    Double_t eff_signal_DstarPi0 = 1;

    Double_t factor = 1;
    Double_t f1 = 1;
    Double_t f2 = 1;
    Double_t f3 = 1;
    Double_t scal_fac = 2.0;

    RooDataHist *data_ellipse;

    //Lets extract the efficiencies
    if (variable == "BDT")
    {
        eff_signal_D = ((scal_fac) * hSD_BDT->GetEntries())/Ns_generated;
        eff_signal_DstarGamma = ((scal_fac) * hSDstarGamma_BDT->GetEntries())/Ns_generated;
        eff_signal_DstarPi0 = ((scal_fac) * hSDstarPi0_BDT->GetEntries())/Ns_generated;
     
        
        cout<<" Data : "<<hist_data->numEntries()<<endl;

        cout<<"NsD: "<<NsD<<endl;
        cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
        cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
        cout<<"Nb: "<<scal_fac * Nb<<endl;
    }
    else if (variable == "Xps")
    {
        eff_signal_D = (hSD_Xps->GetEntries())/Ns_generated;
        eff_signal_DstarGamma = (hSDstarGamma_Xps->GetEntries())/Ns_generated;
        eff_signal_DstarPi0 = (hSDstarPi0_Xps->GetEntries())/Ns_generated;
        
        cout<<" Data : "<<hist_data->numEntries()<<endl;

        cout<<"NsD: "<<NsD<<endl;
        cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
        cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
        cout<<"Nb: "<<Nb<<endl;
    }
    else if (variable == "q2")
    {
        eff_signal_D = (hSD_q2->GetEntries())/Ns_generated;
        eff_signal_DstarGamma = (hSDstarGamma_q2->GetEntries())/Ns_generated;
        eff_signal_DstarPi0 = (hSDstarPi0_q2->GetEntries())/Ns_generated;
        
        cout<<" Data : "<<hist_data->numEntries()<<endl;

        cout<<"NsD: "<<NsD<<endl;
        cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
        cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
        cout<<"Nb: "<<Nb<<endl;
    }
    else if (variable == "Mmin2")
    {
        eff_signal_D = (hSD_Mmin2->GetEntries())/Ns_generated;
        eff_signal_DstarGamma = (hSDstarGamma_Mmin2->GetEntries())/Ns_generated;
        eff_signal_DstarPi0 = (hSDstarPi0_Mmin2->GetEntries())/Ns_generated;
        
        cout<<" Data : "<<hist_data->numEntries()<<endl;

        cout<<"NsD: "<<NsD<<endl;
        cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
        cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
        cout<<"Nb: "<<Nb<<endl;
    }
    else if (variable == "Mmax2")
    {
        eff_signal_D = (hSD_Mmax2->GetEntries())/Ns_generated;
        eff_signal_DstarGamma = (hSDstarGamma_Mmax2->GetEntries())/Ns_generated;
        eff_signal_DstarPi0 = (hSDstarPi0_Mmax2->GetEntries())/Ns_generated;
        
        cout<<" Data : "<<hist_data->numEntries()<<endl;

        cout<<"NsD: "<<NsD<<endl;
        cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
        cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
        cout<<"Nb: "<<Nb<<endl;
    }
    else if (variable == "2D")
    {
        eff_signal_D = (hSD_2D->GetEntries())/Ns_generated;
        eff_signal_DstarGamma = (hSDstarGamma_2D->GetEntries())/Ns_generated;
        eff_signal_DstarPi0 = (hSDstarPi0_2D->GetEntries())/Ns_generated;
        
        cout<<" Data : "<<hist_data->numEntries()<<endl;

        cout<<"NsD: "<<NsD<<endl;
        cout<<"NsDstarGamma: "<<NsDstarGamma<<endl;
        cout<<"NsDstarPi0: "<<NsDstarPi0<<endl;
        cout<<"Nb: "<<Nb<<endl;
    }

    // Creation of the workspace
    //RooWorkspace w("w");

    w.import(*s1_pdf);
    w.import(*s2_pdf);
    w.import(*s3_pdf);
    w.import(*b_pdf);

    w.factory("fac[1]");
    w.factory("f1[1]");
    w.factory("f2[1]");
    w.factory("f3[1]");

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
    cout<<" factor: (2*lum*CS_ee_BB*[Br_tagD*eff_signal_D+Br_tagDstarGamma*eff_signal_DstarGamma+Br_tagDstarPi0*eff_signal_DstarPi0) = "<<factor<<endl;

    // The respective factors for each tag
    f1 = 1.0 *NsD/(NsD+NsDstarGamma+NsDstarPi0);
    f2 = 1.0 *NsDstarGamma/(NsD+NsDstarGamma+NsDstarPi0);
    f3 = 1.0 *NsDstarPi0/(NsD+NsDstarGamma+NsDstarPi0);

    cout<<" f1 = "<<f1<<endl;
    cout<<" f2 = "<<f2<<endl;
    cout<<" f3 = "<<f3<<endl;
    
    // Setting the values
    w.var("fac")->setVal(factor);
    w.var("f1")->setVal(f1);
    w.var("f2")->setVal(f2);
    w.var("f3")->setVal(f3);

    if (variable == "BDT"){
        if (lum == 362){
            w.factory("mu[4.0e-6,-1.0,1.0]");
            w.factory("B[43,2.0,35.0]");

            w.var("B")->setVal(Nb);

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
        if (lum == 50000){
            
            w.factory("mu[1e-8,-1.0.0,1.0]");
            //w.factory("S[0.0,0.0,10.0]");
            w.factory("B[3048.0,3000.0,4000.0]");
        
            w.var("B")->setVal(scal_fac*Nb);

            w.var("mu")->setMin(-0.5e-6);
            w.var("mu")->setMax(0.5e-6);
        }
    }
    if (variable == "Xps"){

        if (lum == 362){
            w.factory("mu[4.0e-6,-1.0,1.0]");
            w.factory("B[43,2.0,35.0]");

            w.var("B")->setVal(Nb);

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
        if (lum == 50000){
            
            w.factory("mu[1e-8,-1.0.0,1.0]");
            //w.factory("S[0.0,0.0,10.0]");
            w.factory("B[3048.0,2000.0,5000.0]");

            w.var("B")->setVal(Nb);
            //w.var("B")->setMin(0);

            // to m_x = 0.0    
            w.var("mu")->setMin(-3e-7);

            // to m_x = 0.0    
            w.var("mu")->setMax(7e-7);
            // to m_x = 4.0
            //w.var("mu")->setMax(4e-7);
            // to m_x = 5.0
            //w.var("mu")->setMax(5e-7);
        }
    }
    if (variable == "q2" || variable == "Mmin2" || variable == "Mmax2"){

        if (lum == 362){
            w.factory("mu[4.0e-6,-1.0,1.0]");
            w.factory("B[43,2.0,35.0]");

            w.var("B")->setVal(Nb);

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
        if (lum == 50000){
            
            w.factory("mu[1e-8,-1.0.0,1.0]");
            //w.factory("S[0.0,0.0,10.0]");
            w.factory("B[3048.0,2000.0,5000.0]");

            w.var("B")->setVal(Nb);
            //w.var("B")->setMin(0);

            // to m_x = 0.0    
            w.var("mu")->setMin(-2e-7);

            // to m_x = 0.0    
            w.var("mu")->setMax(7e-7);
            // to m_x = 4.0
            //w.var("mu")->setMax(4e-7);
            // to m_x = 5.0
            //w.var("mu")->setMax(5e-7);
        }
    }
    else if (variable == "2D"){

        if (lum == 362){
            w.factory("mu[4.0e-6,-1.0,1.0]");
            w.factory("B[43,2.0,35.0]");

            w.var("B")->setVal(Nb);

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
        if (lum == 50000){
            
            w.factory("mu[1e-8,-1.0.0,1.0]");
            //w.factory("S[0.0,0.0,10.0]");
            w.factory("B[3048.0,300.0,4000.0]");
        
            w.var("B")->setVal(Nb);

            w.var("mu")->setMin(-0.09e-6);
            w.var("mu")->setMax(0.5e-6);

            //w.var("mu")->setMin(-0.01e-7);
            //w.var("mu")->setMax(2e-7);


        }
    }

    //RooRealVar meanX("meanX", "Mean X", center_x);
    //RooRealVar meanY("meanY", "Mean Y", center_y);
    //RooRealVar angle("angle", "Angle", angle_degrees);
    //RooRealVar sigmaX("sigmaX", "Sigma X", width);
    //RooRealVar sigmaY("sigmaY", "Sigma Y", height);

    // Now, 'pdf' is composite PDF:
    // pdf(x) = Ns * (f1 * S1(x) + f2 * S2(x) + (1 - f1 - f2) * S3(x)) + Nb * B(x)
    
    w.factory("expr::S('fac*mu',fac, mu)") ;
    w.factory("SUM::s_pdf(f1*s1_pdf,f2*s2_pdf,f3*s3_pdf)");
    w.factory("SUM::model(S*s_pdf, B*b_pdf)");

    /*// Create RooDataSet to hold the original data points
    RooDataSet originalData("originalData", "Original Data", RooArgSet(*x, *y));
    // Loop through each data point in hist_data and add it to originalData
    for (int i = 0; i < hist_data->numEntries(); ++i) {
        const RooArgSet* dataPoint = hist_data->get(i);
        originalData.add(*dataPoint);
    }

    cout << "Original Data Entries: " << originalData.numEntries() << endl;


    // Create another RooDataSet to hold the data points within the ellipse region
    RooDataSet data("data", "Data within Ellipse Region", RooArgSet(*x, *y));

    // Fit the model with the data
    // Select data within the ellipse region

    // Define the equation of the rotated ellipse as a RooFormulaVar
    RooFormulaVar ellipse_constraint("ellipse_constraint", 
        "((((x - meanX) * cos(angle) - (y - meanY) * sin(angle)) / sigmaX)^2 + (((x - meanX) * sin(angle) + (y - meanY) * cos(angle)) / sigmaY)^2) <= 1.0", 
       RooArgList(*x, *y, meanX, meanY, sigmaX, sigmaY, angle));

    // Loop through each data point in hist_data and add it to originalData
    for (int i = 0; i < hist_data->numEntries(); ++i) {
        const RooArgSet* dataPoint = hist_data->get(i);
        double xVal = static_cast<RooRealVar*>(dataPoint->find("x"))->getVal();
        double yVal = static_cast<RooRealVar*>(dataPoint->find("y"))->getVal();
        x->setVal(xVal);
        y->setVal(yVal);
        double constraintVal = ellipse_constraint.getVal();
        //cout << "Data Point " << i << ": x = " << xVal << ", y = " << yVal << ", Constraint = " << constraintVal << endl; // Debug print
        if (constraintVal) {
            //cout << "Adding data point to 'data' dataset." << endl; // Debug print
            data.add(*dataPoint);
        }
    }
    cout << "Data Entries within Ellipse Region: " << data.numEntries() << endl;


    // Create a canvas and draw the 2D histogram
    
    // Create a TGraph representing the ellipse contour
    //TGraph* graph_ellipse = new TGraph();
    //for (int i = 0; i < 360; ++i) {
      //  double angle_rad = TMath::DegToRad() * angle.getVal();
      //  double x_ellipse = meanX.getVal() + sigmaX.getVal() * cos(angle_rad) * cos(i * TMath::DegToRad()) - sigmaY.getVal() * sin(angle_rad) * sin(i * TMath::DegToRad());
      //  double y_ellipse = meanY.getVal() + sigmaY.getVal() * sin(angle_rad) * cos(i * TMath::DegToRad()) + sigmaY.getVal() * cos(angle_rad) * sin(i * TMath::DegToRad());
      //  graph_ellipse->SetPoint(i, x_ellipse, y_ellipse);
    //}
    //exit(0);

    // Create a canvas
    TCanvas* c2D = new TCanvas("c2D", "2D Data and Ellipse", 800, 600);

    // Create a 2D histogram of the data
    //TH2F* hist_2D_data = new TH2F("hist_2D_data", "2D Data Histogram", 50, x->getMin(), x->getMax(), 50, y->getMin(), y->getMax());
    //data.fillHistogram(hist_2D_data, RooArgList(*x, *y));

    // Create a 2D histogram of the originalData
    TH2F* hist_2D_originalData = new TH2F("hist_2D_originalData", "Original Data Histogram", 50, x->getMin(), x->getMax(), 50, y->getMin(), y->getMax());
    originalData.fillHistogram(hist_2D_originalData, RooArgList(*x, *y));
    hist_2D_originalData->SetMarkerStyle(20);
    hist_2D_originalData->SetMarkerSize(0.5);

    // Draw the 2D histogram of the data
    //hist_2D_data->Draw("colz");

    // Draw the 2D histogram of the originalData
    hist_2D_originalData->Draw("colz");
    //hist_2D_originalData->Draw("same text");

    // Create and draw the ellipse
    TEllipse* ellipse = new TEllipse(meanX.getVal(), meanY.getVal(), sigmaX.getVal(), sigmaY.getVal(), 0., 360., angle.getVal());
    ellipse->SetLineColor(kRed);
    ellipse->SetLineWidth(2);
    ellipse->Draw("same");

    // Draw the canvas
    c2D->Draw();
    

    // Create custom PDF for ellipse
    //EllipsePdf ellipsePdf("ellipsePdf", "ellipsePdf", *x, *y, meanX, meanY, sigmaX, sigmaY, angle);

    // Combine ellipse and background PDF
    RooAbsPdf *model_pdf = w.pdf("model");
    RooProdPdf model_scaled("model_scaled", "model_scaled", RooArgSet(*model_pdf, ellipse_constraint));
    w.import(model_scaled);
    //exit(0);
    // Calculate the normalization factor
    //RooAbsReal *norm = model_scaled.createIntegral(RooArgSet(x, y), NormSet(RooArgSet(x, y)), Range("ellipse_region"));

    // Get the value of the normalization factor
    //double norm_value = norm->getVal();

    // Create a scaling factor
    //RooRealVar scale("scale", "scale", 1.0 / norm_value);

    // Scale the model PDF
    //RooRealVar *scaled_model_pdf = new RooRealVar("scaled_model_pdf", "scaled_model_pdf", 0.0, 1.0);
    //RooGenericPdf model_pdf_generic("model_pdf_generic", "model_pdf_generic", "@0*@1", RooArgList(*model_pdf, scale));

    // Construct the normalized model
    //RooAddPdf model_scaled_normalized("model_scaled_normalized", "model_scaled_normalized", RooArgList(model_pdf_generic), RooArgList(*scaled_model_pdf), true);

    // Fit the model to data
    RooFitResult *result = model_scaled.fitTo(*hist_data, Constrain(ellipse_constraint), Save());

    // Plot the result
    RooPlot* frame = x->frame(Title("2D Fit"));
    hist_data->plotOn(frame);
    w.pdf("model").plotOn(frame);

    // Draw the plot
    TCanvas c("2D Fit", "2D Fit", 800, 600);
    frame->Draw();

    // Plot the fitted model
    RooPlot *frameX = x->frame(Title("X Projection"));
    data.plotOn(frameX, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(kFullCircle)); // Plot data points within the circular region
    originalData.plotOn(frameX, MarkerColor(kRed), LineColor(kRed), MarkerStyle(kOpenCircle)); // Plot original data points
    pdf.plotOn(frameX);

    RooPlot *frameY = y.frame(Title("Y Projection"));
    data.plotOn(frameY, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(kFullCircle)); // Plot data points within the circular region
    originalData.plotOn(frameY, MarkerColor(kRed), LineColor(kRed), MarkerStyle(kOpenCircle)); // Plot original data points
    pdf.plotOn(frameY);

    // Create a canvas and draw the plots
    TCanvas *canvas = new TCanvas("canvas", "Fit Result", 800, 600);
    canvas->Divide(2, 1);
    canvas->cd(1);
    frameX->Draw();
    canvas->cd(2);
    frameY->Draw();

    // Save the fit result
    RooWorkspace workspace("workspace");
    workspace.import(pdf);
    workspace.import(*result);
    workspace.writeToFile("workspace.root");

    // Keep the program running to view the plot
    canvas->Draw();

    //exit(0);

    //RooFormulaVar cut_condition("cut_condition", "constraint > 0", RooArgList(constraint));

    // Reduce the dataset using the cut condition
    RooDataSet* reducedData = dynamic_cast<RooDataSet*>(hist_data->reduce(RooArgSet(*x, *y), cut_condition));
    RooDataHist* hist_data_reduced = new RooDataHist("hist_data_reduced", "Reduced histogram data", RooArgList(*x, *y), *reducedData);
    delete reducedData; // Clean up the dynamically allocated reducedData



    // Create a 2D histogram of the data
    TH2F* hist_2D_data = dynamic_cast<TH2F*>(hist_data_reduced->createHistogram("hist_2D_data", *x, RooFit::Binning(50), RooFit::YVar(*y, RooFit::Binning(50))));

    // Create a TGraph representing the ellipse contour
    // Create a TGraph representing the rotated ellipse contour
    TGraph* graph_ellipse = new TGraph();
    for (int i = 0; i < 360; ++i) {
        double angle_rad = TMath::DegToRad() * angle.getVal();
        double x_ellipse = meanX.getVal() + sigmaX.getVal() * cos(angle_rad) * cos(i * TMath::DegToRad()) - sigmaY.getVal() * sin(angle_rad) * sin(i * TMath::DegToRad());
        double y_ellipse = meanY.getVal() + sigmaX.getVal() * sin(angle_rad) * cos(i * TMath::DegToRad()) + sigmaY.getVal() * cos(angle_rad) * sin(i * TMath::DegToRad());
        graph_ellipse->SetPoint(i, x_ellipse, y_ellipse);
    }

    // Create a canvas and draw the 2D histogram
    TCanvas* c2D = new TCanvas("c2D", "2D Data and Ellipse", 800, 600);
    hist_2D_data->Draw("colz");

    // Draw the ellipse contour
    graph_ellipse->SetLineColor(kRed);
    graph_ellipse->SetLineWidth(2);
    graph_ellipse->Draw("same");

    // Draw the canvas
    c2D->Draw();*/
    // Define the category for the inside/outside ellipse
   

    // Perform the fit based on the condition
    if (variable == "2D") {   

        // Fit the PDF model to the original data within the ellipse region
        w.pdf("model")->fitTo(*hist_data, Extended(), Minimizer("Minuit2", "Migrad"), Strategy(2), Save(true), NumCPU(2), Optimize(true));
        //model_scaled.fitTo(data, Constrain(ellipse_constraint), Extended(), Minimizer("Minuit2", "Migrad"), Strategy(2), Minos(true), Save(true), NumCPU(2), Optimize(true));
    } else {
        // Fit the PDF model to the original data
        w.pdf("model")->fitTo(*hist_data, Extended(), Minimizer("Minuit2", "Migrad"), Strategy(2), Minos(true), Save(true), NumCPU(2), Optimize(true), Offset(true));
    }
    //exit(0);


    // Plot the fitted model
    /*RooPlot *frameX = x->frame(Title("X Projection"));
    data.plotOn(frameX, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(kFullCircle)); // Plot data points within the circular region
    originalData.plotOn(frameX, MarkerColor(kRed), LineColor(kRed), MarkerStyle(kOpenCircle)); // Plot original data points
    model_scaled.plotOn(frameX);

    RooPlot *frameY = y->frame(Title("Y Projection"));
    data.plotOn(frameY, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(kFullCircle)); // Plot data points within the circular region
    originalData.plotOn(frameY, MarkerColor(kRed), LineColor(kRed), MarkerStyle(kOpenCircle)); // Plot original data points
    model_scaled.plotOn(frameY);

    // Create a canvas and draw the plots
    TCanvas *canvas = new TCanvas("canvas", "Fit Result", 800, 600);
    canvas->Divide(2, 1);
    canvas->cd(1);
    frameX->Draw();
    canvas->cd(2);
    frameY->Draw();

    // Keep the program running to view the plot
    canvas->Draw();
    // Create a 2D histogram of the data
    TH2F* hist_2D_data = dynamic_cast<TH2F*>(hist_data->createHistogram("hist_2D_data", *x, RooFit::Binning(50), RooFit::YVar(*y, RooFit::Binning(50))));

    // Create a TGraph representing the ellipse contour
    TGraph* graph_ellipse = new TGraph();
    for (int i = 0; i < 360; ++i) {
        double angle_rad = TMath::DegToRad() * angle.getVal();
        double x_ellipse = meanX.getVal() + sigmaX.getVal() * cos(angle_rad) * cos(i * TMath::DegToRad()) - sigmaY.getVal() * sin(angle_rad) * sin(i * TMath::DegToRad());
        double y_ellipse = meanY.getVal() + sigmaX.getVal() * sin(angle_rad) * cos(i * TMath::DegToRad()) + sigmaY.getVal() * cos(angle_rad) * sin(i * TMath::DegToRad());
        graph_ellipse->SetPoint(i, x_ellipse, y_ellipse);
    }

    // Create a canvas and draw the 2D histogram
    TCanvas* c2D = new TCanvas("c2D", "2D Data and Ellipse", 800, 600);
    hist_2D_data->Draw("colz");

    // Draw the ellipse contour
    graph_ellipse->SetLineColor(kRed);
    graph_ellipse->SetLineWidth(2);
    graph_ellipse->Draw("same");

    // Draw the canvas
    c2D->Draw();*/

  

    Double_t mu_val = w.var("mu")->getVal();
    
    // Sets up b-only model 
    RooStats::ModelConfig b_modelNM("b_modelNM", &w);
    if (variable == "BDT" || variable == "Xps" || variable == "q2"  || variable == "Mmin2" || variable == "Mmax2")
    {
        b_modelNM.SetPdf(*w.pdf("model"));
    }
    else if (variable == "2D")
    {
        b_modelNM.SetPdf(*w.pdf("model"));
    }
    b_modelNM.SetParametersOfInterest(*w.var("mu"));
    b_modelNM.SetNuisanceParameters(RooArgSet(*w.var("B")));

    if (variable == "BDT" || variable == "Xps" || variable == "q2"  || variable == "Mmin2" || variable == "Mmax2")
    {
        b_modelNM.SetObservables(*w.var("x")); //This is for 1D limit
    }
    else if (variable == "2D")
    {
        b_modelNM.SetObservables(RooArgSet(*x,*y));   //This is for 2D limit
    }
   
    w.var("mu")->setVal(0.0);
    b_modelNM.SetSnapshot(*w.var("mu"));     // sets up b hypothesis as s = 0
    
    // create the alternate (s+b) ModelConfig with given value of s
    double s = 1;
    RooStats::ModelConfig sb_modelNM("S+B_modelNM", &w);
    if (variable == "BDT" || variable == "Xps" || variable == "q2" || variable == "Mmin2" || variable == "Mmax2")
    {
        sb_modelNM.SetPdf(*w.pdf("model"));
    }
    else if (variable == "2D")
    {
        sb_modelNM.SetPdf(*w.pdf("model"));
    }

    // Replace the original model with the scaled model in sb_modelNM
    if (variable == "2D")
    {
        sb_modelNM.SetObservables(RooArgSet(*x, *y));   // This is for 2D limit
    }
    else
    {
        sb_modelNM.SetObservables(*w.var("x"));  // This is for 1D limit
    }

    sb_modelNM.SetParametersOfInterest(*w.var("mu"));
    RooRealVar* poi = (RooRealVar*)sb_modelNM.GetParametersOfInterest()->first();
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
    //RooStats::LikelihoodInterval* interval = plCalc.GetInterval();
    interval = plCalc.GetInterval();

    if (interval) {
        double lowerLimitL = interval->LowerLimit(*poi);
        double upperLimitL = interval->UpperLimit(*poi);
        
        std::cout << "Lower limit: " << lowerLimitL << std::endl;
        std::cout << "Upper limit: " << upperLimitL << std::endl;

        cout << "RESULT: " << 100*plCalc.ConfidenceLevel() << "% interval is : ["<< lowerLimitL << ", "<< upperLimitL << "]" <<endl;
    } else {
        std::cerr << "Failed to retrieve interval." << std::endl;
    }
    
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
    //RooStats::HypoTestInverterResult* r = calc.GetInterval();
    r = calc.GetInterval();

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


    limitval = r->GetExpectedUpperLimit(0);
    return;
}




Double_t getLimit(Double_t lum, const std::string& variable)
{
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

    if (variable == "BDT")
    {
        //Number of events
        NsD = hSD_BDT->GetEntries();
        NsDstarGamma = hSDstarGamma_BDT->GetEntries();
        NsDstarPi0 = hSDstarPi0_BDT->GetEntries();
        Nb = hB_BDT->GetEntries();
        
        // Set the variable and create the dataset
        x = new RooRealVar("x", "x", BDTLow, BDTUp);
        x->setBins(nbins);
        
        dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_BDT));  // B->Pi+alpha Dtag
        dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_BDT));  // B->Pi+alpha DstarGamma tag
        dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_BDT));  // B->Pi+alpha  Dstar Pi0 tag
        dB = new RooDataHist("dB", "dB", *x, Import(*hB_BDT));  // BKG

        Nb = lum*Nb/Luminosity;
    
        s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",*x,*dSD,2);
        s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",*x,*dSDstarGamma,2);
        s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",*x,*dSDstarPi0,2);
        b_pdf = new RooHistPdf("b_pdf","bkg pdf",*x,*dB,2);
    }
    else if (variable == "Xps")
    {
        //Number of events
        NsD = hSD_Xps->GetEntries();
        NsDstarGamma = hSDstarGamma_Xps->GetEntries();
        NsDstarPi0 = hSDstarPi0_Xps->GetEntries();
        Nb = hB_Xps->GetEntries();

        // Set the variable and create the dataset
        x = new RooRealVar("x", "x", XpsLow,XpsUp);
        x->setBins(nbins);
        
        dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_Xps));  // B->Pi+alpha Dtag
        dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_Xps));  // B->Pi+alpha DstarGamma tag
        dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_Xps));  // B->Pi+alpha  Dstar Pi0 tag
        dB = new RooDataHist("dB", "dB", *x, Import(*hB_Xps));  // BKG

        Nb = lum*Nb/Luminosity;
    
        s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",*x,*dSD,2);
        s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",*x,*dSDstarGamma,2);
        s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",*x,*dSDstarPi0,2);
        b_pdf = new RooHistPdf("b_pdf","bkg pdf",*x,*dB,2);
    }
    else if (variable == "q2")
    {
        //Number of events
        NsD = hSD_q2->GetEntries();
        NsDstarGamma = hSDstarGamma_q2->GetEntries();
        NsDstarPi0 = hSDstarPi0_q2->GetEntries();
        Nb = hB_q2->GetEntries();

        Nb = lum*Nb/Luminosity;
        
        // Set the variable and create the dataset
        x = new RooRealVar("x", "x", q2Low,q2Up);
        x->setBins(nbins);
        
        dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_q2));  // B->Pi+alpha Dtag
        dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_q2));  // B->Pi+alpha DstarGamma tag
        dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_q2));  // B->Pi+alpha  Dstar Pi0 tag
        dB = new RooDataHist("dB", "dB", *x, Import(*hB_q2));  // BKG

    
        s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",*x,*dSD,2);
        s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",*x,*dSDstarGamma,2);
        s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",*x,*dSDstarPi0,2);
        b_pdf = new RooHistPdf("b_pdf","bkg pdf",*x,*dB,2);
    }
    else if (variable == "Mmin2")
    {
        //Number of events
        NsD = hSD_Mmin2->GetEntries();
        NsDstarGamma = hSDstarGamma_Mmin2->GetEntries();
        NsDstarPi0 = hSDstarPi0_Mmin2->GetEntries();
        Nb = hB_Mmin2->GetEntries();

        Nb = lum*Nb/Luminosity;
        
        // Set the variable and create the dataset
        x = new RooRealVar("x", "x", MminLow,MminUp);
        x->setBins(nbins);
        
        dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_Mmin2));  // B->Pi+alpha Dtag
        dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_Mmin2));  // B->Pi+alpha DstarGamma tag
        dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_Mmin2));  // B->Pi+alpha  Dstar Pi0 tag
        dB = new RooDataHist("dB", "dB", *x, Import(*hB_Mmin2));  // BKG

    
        s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",*x,*dSD,2);
        s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",*x,*dSDstarGamma,2);
        s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",*x,*dSDstarPi0,2);
        b_pdf = new RooHistPdf("b_pdf","bkg pdf",*x,*dB,2);
    }
    else if (variable == "Mmax2")
    {
        //Number of events
        NsD = hSD_Mmax2->GetEntries();
        NsDstarGamma = hSDstarGamma_Mmax2->GetEntries();
        NsDstarPi0 = hSDstarPi0_Mmax2->GetEntries();
        Nb = hB_Mmax2->GetEntries();

        Nb = lum*Nb/Luminosity;
        
        // Set the variable and create the dataset
        x = new RooRealVar("x", "x", MmaxLow,MmaxUp);
        x->setBins(nbins);
        
        dSD = new RooDataHist("dSD", "dSD", *x, Import(*hSD_Mmax2));  // B->Pi+alpha Dtag
        dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", *x, Import(*hSDstarGamma_Mmax2));  // B->Pi+alpha DstarGamma tag
        dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", *x, Import(*hSDstarPi0_Mmax2));  // B->Pi+alpha  Dstar Pi0 tag
        dB = new RooDataHist("dB", "dB", *x, Import(*hB_Mmax2));  // BKG

    
        s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",*x,*dSD,2);
        s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",*x,*dSDstarGamma,2);
        s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",*x,*dSDstarPi0,2);
        b_pdf = new RooHistPdf("b_pdf","bkg pdf",*x,*dB,2);
    }
    else if (variable == "2D")
    {   
        //Number of events
        NsD = hSD_2D->GetEntries();
        NsDstarGamma = hSDstarGamma_2D->GetEntries();
        NsDstarPi0 = hSDstarPi0_2D->GetEntries();
        Nb = hB_2D->GetEntries();

        Nb = lum*Nb/Luminosity;

        /*// Identify density regions (for illustration purposes)
        double highDensityThreshold = 20; // Adjust based on your data
        double lowDensityThreshold = 5;   // Adjust based on your preference
        int nbinsX = hB_2D->GetNbinsX();
        int nbinsY = hB_2D->GetNbinsY();

       // Loop over all bins to identify low-density regions
        for (int ix = 1; ix <= nbinsX; ++ix) {
            for (int iy = 1; iy <= nbinsY; ++iy) {
                double binContent = hB_2D->GetBinContent(ix, iy);
                if (binContent == 0) {
                    // Set content to one event
                    hB_2D->SetBinContent(ix, iy, 1); 
                } else if (binContent <= lowDensityThreshold) {
                    // Merge adjacent bins in low-density region
                    int ix_end = ix;
                    while (ix_end <= nbinsX && hB_2D->GetBinContent(ix_end, iy) <= lowDensityThreshold) {
                        ++ix_end;
                    }
                    // Calculate the new range
                    double newRangeMin = hB_2D->GetXaxis()->GetBinLowEdge(ix);
                    double newRangeMax = hB_2D->GetXaxis()->GetBinUpEdge(ix_end);
                    // Set the new range for the X-axis
                    hB_2D->GetXaxis()->SetRangeUser(newRangeMin, newRangeMax);
                    // Move to the next bin after the merged region
                    ix = ix_end;
                } else if (binContent > highDensityThreshold) {
                    // Increase number of bins in high-density region
                    // You can choose the number of additional bins based on your preference
                    int additionalBins = 5; // Adjust based on your preference
                    // Calculate the bin width
                    double binWidth = hB_2D->GetXaxis()->GetBinWidth(ix);
                    // Expand the binning
                    for (int i = 0; i < additionalBins; ++i) {
                        hB_2D->GetXaxis()->SetRangeUser(hB_2D->GetXaxis()->GetBinLowEdge(ix), hB_2D->GetXaxis()->GetBinUpEdge(ix) + binWidth);
                        hB_2D->RebinX(2); // Double the number of bins
                    }
                    // Move to the next bin after the expanded region
                    ix += additionalBins;
                }
            }
        }*/
        // Loop over all bins to identify empty bins and set them to 1 event
        int nbinsX = hB_2D->GetNbinsX();
        int nbinsY = hB_2D->GetNbinsY();

        for (int ix = 1; ix <= nbinsX; ++ix) {
            for (int iy = 1; iy <= nbinsY; ++iy) {
                double binContent = hB_2D->GetBinContent(ix, iy);
                if (binContent == 0) {
                    // Set content to one event for empty bins
                    hB_2D->SetBinContent(ix, iy, 1);
                }
            }
        }


        x = new RooRealVar("x", "x", X_min, X_max);
        y = new RooRealVar("y", "y", Y_min, Y_max);
        //x->setBins(nbinsX);
        //y->setBins(nbinsY);
        x->setBins(hB_2D->GetNbinsX());
        y->setBins(hB_2D->GetNbinsY());


        
        // Import histograms along with the "Inside_Ellipse" label
        dSD = new RooDataHist("dSD", "dSD", RooArgList(*x,*y), Import(*hSD_2D));
        dSDstarGamma = new RooDataHist("dSDstarGamma", "dSDstarGamma", RooArgList(*x,*y), Import(*hSDstarGamma_2D));
        dSDstarPi0 = new RooDataHist("dSDstarPi0", "dSDstarPi0", RooArgList(*x,*y), Import(*hSDstarPi0_2D));
        dB = new RooDataHist("dB", "dB", RooArgList(*x,*y), Import(*hB_2D));

        // Create RooHistPdf objects
        s1_pdf = new RooHistPdf("s1_pdf","signal D tag pdf",RooArgList(*x,*y),*dSD,2);
        s2_pdf = new RooHistPdf("s2_pdf","signal DstarGamma tag pdf",RooArgList(*x,*y),*dSDstarGamma,2);
        s3_pdf = new RooHistPdf("s3_pdf","signal DstarPi0 tag pdf",RooArgList(*x,*y),*dSDstarPi0,2);
        b_pdf = new RooHistPdf("b_pdf","bkg pdf",RooArgList(*x,*y),*dB,2); 
    }
    
    // Plot the signal and backgraound histograms
    TCanvas *c = new TCanvas("c","c",800,600);
    c->Divide(2, 2);

  
    if(variable=="BDT"){
        //TCanvas *c1 = new TCanvas("c1","c1",800,600);
        c->cd(1);
        // Customize histogram attributes
        hSD_BDT->SetLineColor(kBlue);
        hSD_BDT->SetFillColor(kBlue);
        hSD_BDT->SetFillStyle(3001);  // Diagonal lines
      
        hSDstarGamma_BDT->SetLineColor(kRed);
        hSDstarGamma_BDT->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
        hSDstarGamma_BDT->SetMarkerStyle(21);
        hSDstarGamma_BDT->SetMarkerSize(0.7);
      
        hSDstarPi0_BDT->SetLineColor(kGreen);

        // Draw histograms
        hSD_BDT->Draw("HIST");
        hSDstarGamma_BDT->Draw("HISTsame");
        hSDstarPi0_BDT->Draw("HISTsame");
        hB_BDT->Draw("HISTsame");
        //c1->Draw();
        c->Draw();
    }
    else if(variable=="Xps"){
        //TCanvas *c1 = new TCanvas("c1","c1",800,600);
        c->cd(1);
        // Customize histogram attributes
        hSD_Xps->SetLineColor(kBlue);
        hSD_Xps->SetFillColor(kBlue);
        hSD_Xps->SetFillStyle(3001);  // Diagonal lines
      
        hSDstarGamma_Xps->SetLineColor(kRed);
        hSDstarGamma_Xps->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
        hSDstarGamma_Xps->SetMarkerStyle(21);
        hSDstarGamma_Xps->SetMarkerSize(0.7);
      
        hSDstarPi0_Xps->SetLineColor(kGreen);

        // Draw histograms
        hSD_Xps->Draw("HIST");
        hSDstarGamma_Xps->Draw("HISTsame");
        hSDstarPi0_Xps->Draw("HISTsame");
        hB_Xps->Draw("HISTsame");
        //c1->Draw();
        c->Draw();
    }
    else if(variable=="q2"){
        //TCanvas *c1 = new TCanvas("c1","c1",800,600);
        c->cd(1);
        // Customize histogram attributes
        hSD_q2->SetLineColor(kBlue);
        hSD_q2->SetFillColor(kBlue);
        hSD_q2->SetFillStyle(3001);  // Diagonal lines
      
        hSDstarGamma_q2->SetLineColor(kRed);
        hSDstarGamma_q2->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
        hSDstarGamma_q2->SetMarkerStyle(21);
        hSDstarGamma_q2->SetMarkerSize(0.7);
      
        hSDstarPi0_q2->SetLineColor(kGreen);

        // Draw histograms
        hSD_q2->Draw("HIST");
        hSDstarGamma_q2->Draw("HISTsame");
        hSDstarPi0_q2->Draw("HISTsame");
        hB_q2->Draw("HISTsame");
        //c1->Draw();
        c->Draw();
    }
    else if(variable=="Mmin2"){
        //TCanvas *c1 = new TCanvas("c1","c1",800,600);
        c->cd(1);
        // Customize histogram attributes
        hSD_Mmin2->SetLineColor(kBlue);
        hSD_Mmin2->SetFillColor(kBlue);
        hSD_Mmin2->SetFillStyle(3001);  // Diagonal lines
      
        hSDstarGamma_Mmin2->SetLineColor(kRed);
        hSDstarGamma_Mmin2->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
        hSDstarGamma_Mmin2->SetMarkerStyle(21);
        hSDstarGamma_Mmin2->SetMarkerSize(0.7);
      
        hSDstarPi0_Mmin2->SetLineColor(kGreen);

        // Draw histograms
        hSD_Mmin2->Draw("HIST");
        hSDstarGamma_Mmin2->Draw("HISTsame");
        hSDstarPi0_Mmin2->Draw("HISTsame");
        hB_Mmin2->Draw("HISTsame");
        //c1->Draw();
        c->Draw();
    }
    else if(variable=="Mmax2"){
        //TCanvas *c1 = new TCanvas("c1","c1",800,600);
        c->cd(1);
        // Customize histogram attributes
        hSD_Mmax2->SetLineColor(kBlue);
        hSD_Mmax2->SetFillColor(kBlue);
        hSD_Mmax2->SetFillStyle(3001);  // Diagonal lines
      
        hSDstarGamma_Mmax2->SetLineColor(kRed);
        hSDstarGamma_Mmax2->SetFillColorAlpha(kRed, 0.5);  // Transparent fill
        hSDstarGamma_Mmax2->SetMarkerStyle(21);
        hSDstarGamma_Mmax2->SetMarkerSize(0.7);
      
        hSDstarPi0_Mmax2->SetLineColor(kGreen);

        // Draw histograms
        hSD_Mmax2->Draw("HIST");
        hSDstarGamma_Mmax2->Draw("HISTsame");
        hSDstarPi0_Mmax2->Draw("HISTsame");
        hB_Mmax2->Draw("HISTsame");
        //c1->Draw();
        c->Draw();
    }

    // Let's generate the pseudo-data
    if (variable == "BDT")
    {
        hist_dataB = b_pdf->generateBinned(RooArgList(*x),(2.0*Nb)); // the factor of 2.0 it is because in the creation of the BDT variable, we used the 50%
    }
    else if (variable == "Xps" || variable == "q2" || variable == "Mmin2" || variable == "Mmax2")
    {
        hist_dataB = b_pdf->generateBinned(RooArgList(*x),(Nb)); // the factor of 2.0 it is because in the creation of the BDT variable, we used the 50%
    }
    else if (variable == "2D")
    {
        hist_dataB = b_pdf->generateBinned(RooArgList(*x,*y),Nb);
    }
    
    hist_dataB->Print();
    RooDataHist* hist_data = hist_dataB;
    hist_data->Print();  
    //exit(0);

    // Construct the model and estimate the UL
    setupModel(hist_data, s1_pdf, s2_pdf, s3_pdf, b_pdf, lum, variable);

    if (variable == "2D")
    {
        /*c->cd(1);
        RooPlot* frame1 = w.var("x")->frame();
        hist_data->plotOn(frame1);
        w.pdf("model")->plotOn(frame1);
        w.pdf("model")->plotOn(frame1,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        w.pdf("model")->plotOn(frame1,RooFit::Components("b_pdf"),RooFit::LineColor(kRed));
        frame1->Draw();

        c->cd(2);
        RooPlot* frame2 = w.var("y")->frame();
        hist_data->plotOn(frame2);
        w.pdf("model")->plotOn(frame2);
        w.pdf("model")->plotOn(frame2,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        w.pdf("model")->plotOn(frame2,RooFit::Components("b_pdf"),RooFit::LineColor(kRed));
        frame2->Draw();*/

        c->cd(1);
        RooPlot* frame1 = w.var("x")->frame();
        hist_data->plotOn(frame1);
        w.pdf("model")->plotOn(frame1);
        w.pdf("model")->plotOn(frame1,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        w.pdf("model")->plotOn(frame1,RooFit::Components("b_pdf"),RooFit::LineColor(kRed));
        frame1->Draw();

        c->cd(2);
        RooPlot* frame2 = w.var("y")->frame();
        hist_data->plotOn(frame2);
        w.pdf("model")->plotOn(frame2);
        w.pdf("model")->plotOn(frame2,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        w.pdf("model")->plotOn(frame2,RooFit::Components("b_pdf"),RooFit::LineColor(kRed));
        frame2->Draw();
    }
    else 
    {
        c->cd(2);
        RooPlot* frame = w.var("x")->frame();
        hist_data->plotOn(frame);
        w.pdf("model")->plotOn(frame);
        w.pdf("model")->plotOn(frame,RooFit::Components("s_pdf"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        w.pdf("model")->plotOn(frame,RooFit::Components("b_pdf"),RooFit::LineColor(kRed));
        frame->Draw();
    }


    // Profile Likelihood plot
    RooStats::LikelihoodIntervalPlot *plot = new RooStats::LikelihoodIntervalPlot(interval);
    plot->SetNPoints(50); // Use this to reduce sampling granularity (trades speed for precis
    //TCanvas *cR2 = new TCanvas("cR2","cR2",800,600);
    c->cd(3);
    //plot->SetRange(-0.2e-6, 0.3e-6);
    //plot->SetRange(-5e-6,2.0e-6);
    plot->Draw("TF1"); gPad->Draw();

    // plot result of the scan 
    RooStats::HypoTestInverterPlot* plot2 = new RooStats::HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
    c->cd(4);
    plot2->Draw("CLb 2CL");

    return limitval;
    //return 0;
}




