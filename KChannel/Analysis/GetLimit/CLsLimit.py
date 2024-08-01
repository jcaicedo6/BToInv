import ROOT
from ROOT import RooFit, RooStats, RooRealVar, RooArgSet, RooWorkspace, RooDataHist, RooHistPdf, RooArgList, RooDataSet, RooCmdArg
#from ROOT.RooFit import Import
from ROOT import TFile, TH1D, TH2D, TCanvas, TPad, TChain, TLegend
from ROOT import kGreen, kBlue, kRed, kBlack, kOrange, kGreen, kCyan, kMagenta, kAzure, kViolet
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
# Global variables
# Branching ratios extracted from PDG

# For the semileptonic tag
Br_B_to_l_nu_D0 = 2 * 0.0235
Br_B_to_l_nu_D0star = 2 * 0.0558
Br_D0_to_pi_K = 0.03947

Br_D0star_to_D0_gamma = 0.353
Br_D0star_to_D0_pi0 = 0.647
Br_pi0_to_2gamma = 1

# For the three different tags
Br_tagD = Br_B_to_l_nu_D0 * Br_D0_to_pi_K
Br_tagDstarGamma = Br_B_to_l_nu_D0star * Br_D0star_to_D0_gamma * Br_D0_to_pi_K
Br_tagDstarPi0 = Br_B_to_l_nu_D0star * Br_D0star_to_D0_pi0 * Br_D0_to_pi_K * Br_pi0_to_2gamma

Luminosity = 1.08e6 # Normalized to the SemiLepTag_Dstar_pi0 
CS_ee_BB = 0.565e6 # in fb
Ns_generated = 100000



#hSD_BDT, hSDstarGamma_BDT, hSDstarPi0_BDT = None, None, None  # Histograms for BDT
#hSD_Xps, hSDstarGamma_Xps, hSDstarPi0_Xps = None, None, None  # Histograms for Xps
#hSD_q2, hSDstarGamma_q2, hSDstarPi0_q2 = None, None, None  # Histograms for q2
#hSD_Mmin2, hSDstarGamma_Mmin2, hSDstarPi0_Mmin2 = None, None, None  # Histograms for Mmin2
#hSD_Mmax2, hSDstarGamma_Mmax2, hSDstarPi0_Mmax2 = None, None, None  # Histograms for Mmax2
#hSD_2D, hSDstarGamma_2D, hSDstarPi0_2D = None, None, None  # Histograms for 2D

# Histogram limits and bins
BDTLow, BDTUp =-0.05, 1.05;

Mmin2Low, Mmin2Up = -5.0, 23.0
Mmax2Low, Mmax2Up = -2.0, 23.0
#XpsLow = 0.2, XpsUp = 1.2
#q2Low = -4.0, q2Up = 22

XpsLow, XpsUp = 0.2, 1.2
q2Low, q2Up = -4.0, 22

X_min, X_max = -5.5, 33.5;
Y_min, Y_max = -0.1, 2.6;

nbins = 200
nbinsX, nbinsY = 150, 150


# Scale factor for the BDT histograms
scal_fac = 2.0



# BDT histograms
'''hSD_BDT = TH1D("hSD_BDT", "hSD_BDT", nbins, BDTLow, BDTUp)
hSDstarGamma_BDT = TH1D("hSDstarGamma_BDT", "hSDstarGamma_BDT", nbins, BDTLow, BDTUp)
hSDstarPi0_BDT = TH1D("hSDstarPi0_BDT", "hSDstarPi0_BDT", nbins, BDTLow, BDTUp)

hBkgD_BDT = TH1D("hBkgD_BDT", "hBkgD_BDT", nbins, BDTLow, BDTUp)
hBkgDstarGamma_BDT = TH1D("hBkgDstarGamma_BDT", "hBkgDstarGamma_BDT", nbins, BDTLow, BDTUp)
hBkgDstarPi0_BDT = TH1D("hBkgDstarPi0_BDT", "hBkgDstarPi0_BDT", nbins, BDTLow, BDTUp)
hB1_BDT = TH1D("hB1_BDT", "hB1_BDT", nbins, BDTLow, BDTUp)
hB_BDT = TH1D("hB_BDT", "hB_BDT", nbins, BDTLow, BDTUp)

# Xps histograms
hSD_Xps = TH1D("hSD_Xps", "hSD_Xps", nbins, XpsLow, XpsUp)
hSDstarGamma_Xps = TH1D("hSDstarGamma_Xps", "hSDstarGamma_Xps", nbins, XpsLow, XpsUp)
hSDstarPi0_Xps = TH1D("hSDstarPi0_Xps", "hSDstarPi0_Xps", nbins, XpsLow, XpsUp)

hBkgD_Xps = TH1D("hBkgD_Xps", "hBkgD_Xps", nbins, XpsLow, XpsUp)
hBkgDstarGamma_Xps = TH1D("hBkgDstarGamma_Xps", "hBkgDstarGamma_Xps", nbins, XpsLow, XpsUp)
hBkgDstarPi0_Xps = TH1D("hBkgDstarPi0_Xps", "hBkgDstarPi0_Xps", nbins, XpsLow, XpsUp)
hB1_Xps = TH1D("hB1_Xps", "hB1_Xps", nbins, XpsLow, XpsUp)
hB_Xps = TH1D("hB_Xps", "hB_Xps", nbins, XpsLow, XpsUp)

# q2 histograms
hSD_q2 = TH1D("hSD_q2", "hSD_q2", nbins, q2Low, q2Up)
hSDstarGamma_q2 = TH1D("hSDstarGamma_q2", "hSDstarGamma_q2", nbins, q2Low, q2Up)
hSDstarPi0_q2 = TH1D("hSDstarPi0_q2", "hSDstarPi0_q2", nbins, q2Low, q2Up)

hBkgD_q2 = TH1D("hBkgD_q2", "hBkgD_q2", nbins, q2Low, q2Up)
hBkgDstarGamma_q2 = TH1D("hBkgDstarGamma_q2", "hBkgDstarGamma_q2", nbins, q2Low, q2Up)
hBkgDstarPi0_q2 = TH1D("hBkgDstarPi0_q2", "hBkgDstarPi0_q2", nbins, q2Low, q2Up)
hB1_q2 = TH1D("hB1_q2", "hB1_q2", nbins, q2Low, q2Up)
hB_q2 = TH1D("hB_q2", "hB_q2", nbins, q2Low, q2Up)

# Mmin2 histograms
hSD_Mmin2 = TH1D("hSD_Mmin2", "hSD_Mmin2", nbins, MminLow, MminUp)
hSDstarGamma_Mmin2 = TH1D("hSDstarGamma_Mmin2", "hSDstarGamma_Mmin2", nbins, MminLow, MminUp)
hSDstarPi0_Mmin2 = TH1D("hSDstarPi0_Mmin2", "hSDstarPi0_Mmin2", nbins, MminLow, MminUp)

hBkgD_Mmin2 = TH1D("hBkgD_Mmin2", "hBkgD_Mmin2", nbins, MminLow, MminUp)
hBkgDstarGamma_Mmin2 = TH1D("hBkgDstarGamma_Mmin2", "hBkgDstarGamma_Mmin2", nbins, MminLow, MminUp)
hBkgDstarPi0_Mmin2 = TH1D("hBkgDstarPi0_Mmin2", "hBkgDstarPi0_Mmin2", nbins, MminLow, MminUp)
hB1_Mmin2 = TH1D("hB1_Mmin2", "hB1_Mmin2", nbins, MminLow, MminUp)
hB_Mmin2 = TH1D("hB_Mmin2", "hB_Mmin2", nbins, MminLow, MminUp)

# Mmax2 histograms
hSD_Mmax2 = TH1D("hSD_Mmax2", "hSD_Mmax2", nbins, MmaxLow, MmaxUp)
hSDstarGamma_Mmax2 = TH1D("hSDstarGamma_Mmax2", "hSDstarGamma_Mmax2", nbins, MmaxLow, MmaxUp)
hSDstarPi0_Mmax2 = TH1D("hSDstarPi0_Mmax2", "hSDstarPi0_Mmax2", nbins, MmaxLow, MmaxUp)

hBkgD_Mmax2 = TH1D("hBkgD_Mmax2", "hBkgD_Mmax2", nbins, MmaxLow, MmaxUp)
hBkgDstarGamma_Mmax2 = TH1D("hBkgDstarGamma_Mmax2", "hBkgDstarGamma_Mmax2", nbins, MmaxLow, MmaxUp)
hBkgDstarPi0_Mmax2 = TH1D("hBkgDstarPi0_Mmax2", "hBkgDstarPi0_Mmax2", nbins, MmaxLow, MmaxUp)
hB1_Mmax2 = TH1D("hB1_Mmax2", "hB1_Mmax2", nbins, MmaxLow, MmaxUp)
hB_Mmax2 = TH1D("hB_Mmax2", "hB_Mmax2", nbins, MmaxLow, MmaxUp)

# 2D histograms
hSD_2D = TH2D("hSD_2D", "hSD_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)
hSDstarGamma_2D = TH2D("hSDstarGamma_2D", "hSDstarGamma_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)
hSDstarPi0_2D = TH2D("hSDstarPi0_2D", "hSDstarPi0_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)

hBkgD_2D = TH2D("hBkgD_2D", "hBkgD_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)
hBkgDstarGamma_2D = TH2D("hBkgDstarGamma_2D", "hBkgDstarGamma_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)
hBkgDstarPi0_2D = TH2D("hBkgDstarPi0_2D", "hBkgDstarPi0_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)

hB1_2D = TH2D("hB1_2D", "hB1_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)
hB_2D = TH2D("hB_2D", "hB_2D", nbinsX, X_min, X_max, nbinsY, Y_min, Y_max)'''

# Read the histograms, first load the files
def load_signal_files(base_dir, mass, variable):
    fileNameSD = ""
    fileNameSDstarGamma = ""
    fileNameSDstarPi0 = ""
    fileNameBkgD = ""
    fileNameBkgDstarGamma = ""
    fileNameBkgDstarPi0 = ""

    if variable == "BDT":
        fileNameSD = base_dir + "hS_D_m"
        fileNameSDstarGamma = base_dir + "hS_Dstar_gamma_m"
        fileNameSDstarPi0 = base_dir + "hS_Dstar_pi0_m"
        fileNameBkgD = base_dir + "hBkgStag_D_m"
        fileNameBkgDstarGamma = base_dir + "hBkgStag_Dstar_gamma_m"
        fileNameBkgDstarPi0 = base_dir + "hBkgStag_Dstar_pi0_m"
    elif variable in ["2D", "Xps", "q2", "Mmin2", "Mmax2"]:
        fileNameSD = base_dir + "hS_D_m"
        fileNameSDstarGamma = base_dir + "hS_Dstar_gamma_m"
        fileNameSDstarPi0 = base_dir + "hS_Dstar_pi0_m"
        fileNameBkgD = base_dir + "hBkgStag_D_"
        fileNameBkgDstarGamma = base_dir + "hBkgStag_Dstar_gamma_"
        fileNameBkgDstarPi0 = base_dir + "hBkgStag_Dstar_pi0_"

    if variable == "BDT":
        fSD = ROOT.TFile(fileNameSD + mass + "_" + variable + "_pdf.root")
        fSDstarGamma = ROOT.TFile(fileNameSDstarGamma + mass + "_" + variable + "_pdf.root")
        fSDstarPi0 = ROOT.TFile(fileNameSDstarPi0 + mass + "_" + variable + "_pdf.root")

        fBkgD = ROOT.TFile(fileNameBkgD + mass + "_" + variable + "_pdf.root")
        fBkgDstarGamma = ROOT.TFile(fileNameBkgDstarGamma + mass + "_" + variable + "_pdf.root")
        fBkgDstarPi0 = ROOT.TFile(fileNameBkgDstarPi0 + mass + "_" + variable + "_pdf.root")
    elif variable in ["Xps", "q2", "Mmin2", "Mmax2"]:
        fSD = ROOT.TFile(fileNameSD + mass + "_" + variable + "_pdf.root")
        fSDstarGamma = ROOT.TFile(fileNameSDstarGamma + mass + "_" + variable + "_pdf.root")
        fSDstarPi0 = ROOT.TFile(fileNameSDstarPi0 + mass + "_" + variable + "_pdf.root")

        fBkgD = ROOT.TFile(fileNameBkgD + variable + "_pdf.root")
        fBkgDstarGamma = ROOT.TFile(fileNameBkgDstarGamma + variable + "_pdf.root")
        fBkgDstarPi0 = ROOT.TFile(fileNameBkgDstarPi0 + variable + "_pdf.root")
        print(fileNameSD + mass + "_" + variable + "_pdf.root")
    elif variable == "2D":
        fSD = ROOT.TFile(fileNameSD + mass + "_" + variable + "_ROT_pdf.root")
        fSDstarGamma = ROOT.TFile(fileNameSDstarGamma + mass + "_" + variable + "_ROT_pdf.root")
        fSDstarPi0 = ROOT.TFile(fileNameSDstarPi0 + mass + "_" + variable + "_ROT_pdf.root")

        fBkgD = ROOT.TFile(fileNameBkgD + variable + "_ROT_pdf.root")
        fBkgDstarGamma = ROOT.TFile(fileNameBkgDstarGamma + variable + "_ROT_pdf.root")
        fBkgDstarPi0 = ROOT.TFile(fileNameBkgDstarPi0 + variable + "_ROT_pdf.root")
        print(fileNameSD + mass + "_" + variable + "_ROT_pdf.root")

    if not (fSD and fSDstarGamma and fSDstarPi0 and fBkgD and fBkgDstarGamma and fBkgDstarPi0):
        raise RuntimeError("Error opening one or more files.")
    
    return fSD, fSDstarGamma, fSDstarPi0, fBkgD, fBkgDstarGamma, fBkgDstarPi0

def construct_hist(fS, fB, variable, histo_type):
    hS = fS.Get(variable)
    if not hS:
        raise RuntimeError(f"Error: Unable to retrieve histogram {variable} from file S.")
    
    hB = fB.Get(variable)
    if not hB:
        raise RuntimeError(f"Error: Unable to retrieve histogram {variable} from file B.")
    
    return hS, hB

def load_all_histograms(variable):
    global hSD_BDT, hBkgD_BDT, hSDstarGamma_BDT, hBkgDstarGamma_BDT, hSDstarPi0_BDT, hBkgDstarPi0_BDT
    global hSD_Xps, hBkgD_Xps, hSDstarGamma_Xps, hBkgDstarGamma_Xps, hSDstarPi0_Xps, hBkgDstarPi0_Xps
    global hSD_q2, hBkgD_q2, hSDstarGamma_q2, hBkgDstarGamma_q2, hSDstarPi0_q2, hBkgDstarPi0_q2
    global hSD_Mmin2, hBkgD_Mmin2, hSDstarGamma_Mmin2, hBkgDstarGamma_Mmin2, hSDstarPi0_Mmin2, hBkgDstarPi0_Mmin2
    global hSD_Mmax2, hBkgD_Mmax2, hSDstarGamma_Mmax2, hBkgDstarGamma_Mmax2, hSDstarPi0_Mmax2, hBkgDstarPi0_Mmax2
    global hSD_2D, hBkgD_2D, hSDstarGamma_2D, hBkgDstarGamma_2D, hSDstarPi0_2D, hBkgDstarPi0_2D
    global hB_BDT, hB1_BDT, hB_Xps, hB1_Xps, hB_q2, hB1_q2, hB_Mmin2, hB1_Mmin2, hB_Mmax2, hB1_Mmax2, hB_2D, hB1_2D
    
    if variable == "BDT":
        hSD_BDT, hBkgD_BDT = construct_hist(fSD, fBkgD, variable, ROOT.TH1D)
        hSDstarGamma_BDT, hBkgDstarGamma_BDT = construct_hist(fSDstarGamma, fBkgDstarGamma, variable, ROOT.TH1D)
        hSDstarPi0_BDT, hBkgDstarPi0_BDT = construct_hist(fSDstarPi0, fBkgDstarPi0, variable, ROOT.TH1D)
        
        hB1_BDT = hBkgD_BDT.Clone()
        hB1_BDT.Add(hBkgDstarGamma_BDT)
        hB_BDT = hB1_BDT.Clone()
        hB_BDT.Add(hBkgDstarPi0_BDT)
        
    elif variable == "Xps":
        hSD_Xps, hBkgD_Xps = construct_hist(fSD, fBkgD, variable, ROOT.TH1D)
        hSDstarGamma_Xps, hBkgDstarGamma_Xps = construct_hist(fSDstarGamma, fBkgDstarGamma, variable, ROOT.TH1D)
        hSDstarPi0_Xps, hBkgDstarPi0_Xps = construct_hist(fSDstarPi0, fBkgDstarPi0, variable, ROOT.TH1D)
        
        hB1_Xps = hBkgD_Xps.Clone()
        hB1_Xps.Add(hBkgDstarGamma_Xps)
        hB_Xps = hB1_Xps.Clone()
        hB_Xps.Add(hBkgDstarPi0_Xps)
        
    elif variable == "q2":
        hSD_q2, hBkgD_q2 = construct_hist(fSD, fBkgD, variable, ROOT.TH1D)
        hSDstarGamma_q2, hBkgDstarGamma_q2 = construct_hist(fSDstarGamma, fBkgDstarGamma, variable, ROOT.TH1D)
        hSDstarPi0_q2, hBkgDstarPi0_q2 = construct_hist(fSDstarPi0, fBkgDstarPi0, variable, ROOT.TH1D)
        
        hB1_q2 = hBkgD_q2.Clone()
        hB1_q2.Add(hBkgDstarGamma_q2)
        hB_q2 = hB1_q2.Clone()
        hB_q2.Add(hBkgDstarPi0_q2)
        
    elif variable == "Mmin2":
        hSD_Mmin2, hBkgD_Mmin2 = construct_hist(fSD, fBkgD, variable, ROOT.TH1D)
        hSDstarGamma_Mmin2, hBkgDstarGamma_Mmin2 = construct_hist(fSDstarGamma, fBkgDstarGamma, variable, ROOT.TH1D)
        hSDstarPi0_Mmin2, hBkgDstarPi0_Mmin2 = construct_hist(fSDstarPi0, fBkgDstarPi0, variable, ROOT.TH1D)
        
        hB1_Mmin2 = hBkgD_Mmin2.Clone()
        hB1_Mmin2.Add(hBkgDstarGamma_Mmin2)
        hB_Mmin2 = hB1_Mmin2.Clone()
        hB_Mmin2.Add(hBkgDstarPi0_Mmin2)
        
    elif variable == "Mmax2":
        hSD_Mmax2, hBkgD_Mmax2 = construct_hist(fSD, fBkgD, variable, ROOT.TH1D)
        hSDstarGamma_Mmax2, hBkgDstarGamma_Mmax2 = construct_hist(fSDstarGamma, fBkgDstarGamma, variable, ROOT.TH1D)
        hSDstarPi0_Mmax2, hBkgDstarPi0_Mmax2 = construct_hist(fSDstarPi0, fBkgDstarPi0, variable, ROOT.TH1D)
        
        hB1_Mmax2 = hBkgD_Mmax2.Clone()
        hB1_Mmax2.Add(hBkgDstarGamma_Mmax2)
        hB_Mmax2 = hB1_Mmax2.Clone()
        hB_Mmax2.Add(hBkgDstarPi0_Mmax2)
        
    elif variable == "2D":
        hSD_2D, hBkgD_2D = construct_hist(fSD, fBkgD, variable, ROOT.TH2D)
        hSDstarGamma_2D, hBkgDstarGamma_2D = construct_hist(fSDstarGamma, fBkgDstarGamma, variable, ROOT.TH2D)
        hSDstarPi0_2D, hBkgDstarPi0_2D = construct_hist(fSDstarPi0, fBkgDstarPi0, variable, ROOT.TH2D)
        
        hB1_2D = hBkgD_2D.Clone()
        hB1_2D.Add(hBkgDstarGamma_2D)
        hB_2D = hB1_2D.Clone()
        hB_2D.Add(hBkgDstarPi0_2D)

def set_data(mass, variable):
    base_dir = "~/BToInv/KChannel/Analysis/Files/RootFiles/"
    
    global fSD, fSDstarGamma, fSDstarPi0, fBkgD, fBkgDstarGamma, fBkgDstarPi0
    fSD, fSDstarGamma, fSDstarPi0, fBkgD, fBkgDstarGamma, fBkgDstarPi0 = load_signal_files(base_dir, mass, variable)
    
    load_all_histograms(variable)

def extract_efficiencies(variable):
    # Extract efficiencies
    if variable == "BDT":
        eff_signal_D = (scal_fac * hSD_BDT.GetEntries()) / Ns_generated
        eff_signal_DstarGamma = (scal_fac * hSDstarGamma_BDT.GetEntries()) / Ns_generated
        eff_signal_DstarPi0 = (scal_fac * hSDstarPi0_BDT.GetEntries()) / Ns_generated
    elif variable == "Xps":
        eff_signal_D = (hSD_Xps.GetEntries()) / Ns_generated
        eff_signal_DstarGamma = (hSDstarGamma_Xps.GetEntries()) / Ns_generated
        eff_signal_DstarPi0 = (hSDstarPi0_Xps.GetEntries()) / Ns_generated
    elif variable == "q2":
        eff_signal_D = (hSD_q2.GetEntries()) / Ns_generated
        eff_signal_DstarGamma = (hSDstarGamma_q2.GetEntries()) / Ns_generated
        eff_signal_DstarPi0 = (hSDstarPi0_q2.GetEntries()) / Ns_generated
    elif variable == "Mmin2":
        eff_signal_D = (hSD_Mmin2.GetEntries()) / Ns_generated
        eff_signal_DstarGamma = (hSDstarGamma_Mmin2.GetEntries()) / Ns_generated
        eff_signal_DstarPi0 = (hSDstarPi0_Mmin2.GetEntries()) / Ns_generated
    elif variable == "Mmax2":
        eff_signal_D = (hSD_Mmax2.GetEntries()) / Ns_generated
        eff_signal_DstarGamma = (hSDstarGamma_Mmax2.GetEntries()) / Ns_generated
        eff_signal_DstarPi0 = (hSDstarPi0_Mmax2.GetEntries()) / Ns_generated
    elif variable == "2D":
        eff_signal_D = (hSD_2D.GetEntries()) / Ns_generated
        eff_signal_DstarGamma = (hSDstarGamma_2D.GetEntries()) / Ns_generated
        eff_signal_DstarPi0 = (hSDstarPi0_2D.GetEntries()) / Ns_generated
    else:
        raise ValueError("Unknown variable: {variable}")

    print("Efficiencies: {eff_signal_D}, {eff_signal_DstarGamma}, {eff_signal_DstarPi0}".format(eff_signal_D=eff_signal_D, eff_signal_DstarGamma=eff_signal_DstarGamma, eff_signal_DstarPi0=eff_signal_DstarPi0))

    return eff_signal_D, eff_signal_DstarGamma, eff_signal_DstarPi0

# Create the workspace
def setup_workspace():
    w = RooWorkspace("w")
    return w


# Fit the model
def fit_model(w, hist_data):
    model = w.pdf("model")
    model.fitTo(hist_data, RooFit.Extended(), RooFit.Minimizer("Minuit2", "Migrad"), RooFit.Strategy(2), RooFit.Minos(True), RooFit.Save(True), RooFit.NumCPU(2), RooFit.Optimize(True), RooFit.Offset(True))

    mu_val = w.var("mu").getVal()


# Set up background model
def setup_b_model(w, variable, x=None, y=None):
    b_modelNM = RooStats.ModelConfig("b_modelNM", w)
    b_modelNM.SetPdf(w.pdf("model"))
    b_modelNM.SetParametersOfInterest(w.var("mu"))
    b_modelNM.SetNuisanceParameters(RooArgSet(w.var("B")))

    if variable in ["BDT", "Xps", "q2", "Mmin2", "Mmax2"]:
        b_modelNM.SetObservables(w.var("x"))  # 1D limit
    elif variable == "2D":
        b_modelNM.SetObservables(RooArgSet(x, y))  # 2D limit

    w.var("mu").setVal(0.0)
    b_modelNM.SetSnapshot(w.var("mu"))  # b hypothesis as s = 0

    b_modelNM.Print();

    return b_modelNM

# Set up background + signal model
def setup_sb_model(w, variable, x=None, y=None):
    sb_modelNM = RooStats.ModelConfig("S+B_modelNM", w)
    sb_modelNM.SetPdf(w.pdf("model"))

    if variable in ["BDT", "Xps", "q2", "Mmin2", "Mmax2"]:
        sb_modelNM.SetObservables(w.var("x"))  # 1D limit
    elif variable == "2D":
        sb_modelNM.SetObservables(RooArgSet(x, y))  # 2D limit

    sb_modelNM.SetParametersOfInterest(w.var("mu"))
    poi = sb_modelNM.GetParametersOfInterest().first()
    # Setting s = 1, in order to have signal in the model
    s = 1
    w.var("mu").setVal(s)
    sb_modelNM.SetSnapshot(w.var("mu"))  # sb hypothesis with given s

    return sb_modelNM, poi

def get_limit(w, hist_data, variable, x=None, y=None):
    if variable in ["BDT", "Xps", "q2", "Mmin2", "Mmax2"]:
        b_modelNM = setup_b_model(w, variable)
        sb_modelNM, poi = setup_sb_model(w, variable)
    elif variable == "2D":
        b_modelNM = setup_b_model(w, variable, x, y)
        sb_modelNM, poi = setup_sb_model(w, variable, x, y)
        
    # Print both models
    b_modelNM.Print()
    sb_modelNM.Print()
    
    # Proceding with the statistical test
    # Statistic test \lambda(s) = -log L(s,\hat\hat{b})/L(\hat{s},\hat{b})
    profll = RooStats.ProfileLikelihoodTestStat(sb_modelNM.GetPdf())

    # The interval calculator
    plCalc = RooStats.ProfileLikelihoodCalculator(hist_data, sb_modelNM)

    # Seting confidence level for the Likelihood Ratio
    plCalc.SetConfidenceLevel(0.68)

    # The interval
    interval = plCalc.GetInterval()
    print("Interval: ", interval)

    if interval:
        lowerLimitL = interval.LowerLimit(poi)
        upperLimitL = interval.UpperLimit(poi)
        print("Lower Limit: ", lowerLimitL)
        print("Upper Limit: ", upperLimitL)
        print("RESULT: ", 100 * plCalc.ConfidenceLevel(), "%The interval is: [", lowerLimitL, ", ", upperLimitL, "]")
    else:
        print("Failed to retreive interval.")

    # Now compute using asymptotic formula; b is alt, sb is null
    ac = RooStats.AsymptoticCalculator(hist_data, b_modelNM, sb_modelNM)
    ac.SetOneSided(True)  # Should want 1-side (true) for the Limit
    RooStats.AsymptoticCalculator.SetPrintLevel(-1) # No printout
    # Get Hypotesis test
    asymResult = ac.GetHypoTest()
    asymResult.SetPValueIsRightTail(False)
    asym_pb = 1 - asymResult.AlternatePValue()
    asymResult.SetPValueIsRightTail(True)
    asym_psb = asymResult.NullPValue()

    print("Results based on asymptotic formulae:")
    print("Asymptotic p-value: ", asym_pb)
    print("Asymptotic p-value sb: ", asym_psb)

    asymp_clb = 1 - asym_pb
    asymp_clsb = asym_psb
    asymp_cls = 9999

    if asymp_clb > 0:
        asymp_cls = asymp_clb / asymp_clsb

    print("Asymptotic CLs: ", asymp_cls)

    # Create hypotest inverter passing the desired calculator (hc or ac)
    calc = RooStats.HypoTestInverter(ac)
    calc.SetVerbose(False)
    calc.SetConfidenceLevel(0.95)
    useCLs = True
    calc.UseCLs(useCLs)
    if useCLs:
        profll.SetOneSided(True)

    npoints = 150 # Number of points to scan
    # min and max for scan (better to choose smaller intervals)
    poimin = poi.getMin()
    poimax = poi.getMax()
    print("Doing a fixed scan in interval: [", poimin, ", ", poimax, "]")

    calc.SetFixedScan(npoints, poimin, poimax)  # Fixed scan
    r = calc.GetInterval()

    upperLimit = r.UpperLimit()
    ulError = r.UpperLimitEstimatedError()
    print("The computed limit is: ", upperLimit, " +/- ", ulError)

    # Compute the expected upper limit
    print("Expected upper limits using b-only model:")
    print("median limit: ", r.GetExpectedUpperLimit(0))
    print("limit (-2 sigma): ", r.GetExpectedUpperLimit(-2))
    print("limit (+2 sigma): ", r.GetExpectedUpperLimit(2))
    print("limit (-1 sigma): ", r.GetExpectedUpperLimit(-1))
    print("limit (+1 sigma): ", r.GetExpectedUpperLimit(1))

def generate_pseudo_data(variable, b_pdf, Nb, x, y=None):
    """
    Generate pseudo-data based on the variable type.
    
    Parameters:
        variable (str): Variable name.
        b_pdf (RooHistPdf): Background PDF.
    
    Returns:
        RooDataHist: Generated pseudo-data histogram.
    """
    hist_dataB = None
    # Generate pseudo-data logic here based on variable type
    if variable == "BDT":
        hist_dataB = b_pdf.generateBinned(RooArgList(x), 2.0 * Nb)
    elif variable in ["Xps", "q2", "Mmin2", "Mmax2"]:
        hist_dataB = b_pdf.generateBinned(RooArgList(x), Nb)
    elif variable == "2D":
        hist_dataB = b_pdf.generateBinned(RooArgList(x, y), Nb)

    hist_dataB.Print()

    if variable in ["Xps", "q2", "Mmin2", "Mmax2", "BDT"]:
        hist_data = ROOT.RooDataHist("hist_data", "hist_data", ROOT.RooArgList(x), hist_dataB)
    else:
        hist_data = ROOT.RooDataHist("hist_data", "hist_data", ROOT.RooArgList(x, y), hist_dataB)
    hist_data.Print()
    

    return hist_data

def create_model(w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb, lum, variable):
    
    factor = 1;
    f1 = 1;
    f2 = 1;
    f3 = 1;
    scal_fac = 2.0;

    eff_signal_D, eff_signal_DstarGamma, eff_signal_DstarPi0 = extract_efficiencies(variable)

    # Print number of events
    print(f"Data: {hist_data.numEntries()} events")
    print(f"Signal D: {NsD} events")
    print(f"Signal DstarGamma: {NsDstarGamma} events")
    print(f"Signal DstarPi0: {NsDstarPi0} events")
    print(f"Background: {Nb} events")


    # Setup the workspace

    w.factory("fac[1]")
    w.factory("f1[1]")
    w.factory("f2[1]")
    w.factory("f3[1]")

    # Weighted branch fractions
    B_eff1 = Br_tagD * eff_signal_D
    B_eff2 = Br_tagDstarGamma * eff_signal_DstarGamma
    B_eff3 = Br_tagDstarPi0 * eff_signal_DstarPi0
    lum_CS = lum * CS_ee_BB

    # Defining the factor = 2 * L * sigma(ee->BB) * (Br_tagD * eff_signal_D + Br_tagDstarGamma * eff_signal_DstarGamma + Br_tagDstarPi0 * eff_signal_DstarPi0)
    factor = 2 * lum_CS * (B_eff1 + B_eff2 + B_eff3)
    # Defining the fractions of D, DstarGamma, and DstarPi0 in the total signal sample
    f1 = NsD / (NsD + NsDstarGamma + NsDstarPi0)
    f2 = NsDstarGamma / (NsD + NsDstarGamma + NsDstarPi0)
    f3 = NsDstarPi0 / (NsD + NsDstarGamma + NsDstarPi0)

    # Set values of the parameters
    w.var("fac").setVal(factor)
    w.var("f1").setVal(f1)
    w.var("f2").setVal(f2)
    w.var("f3").setVal(f3)

    if variable == "BDT":
        if lum == 362:
            w.factory("mu[4.0e-6,-1.0,1.0]")
            w.factory("B[43,2.0,35.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(0)
            w.var("mu").setMax(4.0e-5)
        elif lum == 50000:
            w.factory("mu[1e-8,-1.0,1.0]")
            w.factory("B[3048.0,3000.0,5000.0]")
            w.var("B").setVal(scal_fac * Nb)
            w.var("mu").setMin(-0.5e-6)
            w.var("mu").setMax(0.5e-5)
    elif variable == "Xps":
        if lum == 362:
            w.factory("mu[4.0e-6,-1.0,1.0]")
            w.factory("B[43,2.0,35.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(0)
            w.var("mu").setMax(4.0e-5)
        elif lum == 50000:
            w.factory("mu[1e-8,-1.0,1.0]")
            w.factory("B[3048.0,2000.0,5000.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(-3e-7)
            w.var("mu").setMax(7e-7)
    elif variable in ["q2", "Mmin2", "Mmax2"]:
        if lum == 362:
            w.factory("mu[4.0e-6,-1.0,1.0]")
            w.factory("B[43,2.0,35.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(0)
            w.var("mu").setMax(4.0e-5)
        elif lum == 50000:
            w.factory("mu[1e-8,-1.0,1.0]")
            w.factory("B[3048.0,2000.0,5000.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(-2e-7)
            w.var("mu").setMax(7e-7)
    elif variable == "2D":
        if lum == 362:
            w.factory("mu[4.0e-6,-1.0,1.0]")
            w.factory("B[43,2.0,35.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(0)
            w.var("mu").setMax(4.0e-5)
        elif lum == 50000:
            w.factory("mu[1e-8,-1.0,1.0]")
            w.factory("B[3048.0,2000.0,5000.0]")
            w.var("B").setVal(Nb)
            w.var("mu").setMin(-0.09e-6)
            w.var("mu").setMax(0.5e-6)

    # Composing the model
    w.factory("expr::S('fac*mu',fac, mu)") 
    w.factory("SUM::s_pdf(f1*s1_pdf,f2*s2_pdf,f3*s3_pdf)")
    w.factory("SUM::model(S*s_pdf, B*b_pdf)")
    

    
def plot_histograms(histograms, variable):
    c = TCanvas("c","c",800,600)
    c.Divide(2, 2)

    if variable == "BDT":
        c.cd(1)
        # Customize histogram attributes
        hSD_BDT.SetLineColor(kBlue)
        hSD_BDT.SetFillColor(kBlue)
        hSD_BDT.SetFillStyle(3001)  # Diagonal lines
      
        hSDstarGamma_BDT.SetLineColor(kRed)
        hSDstarGamma_BDT.SetFillColorAlpha(kRed, 0.5);  # Transparent fill
        hSDstarGamma_BDT.SetMarkerStyle(21)
        hSDstarGamma_BDT.SetMarkerSize(0.7)
      
        hSDstarGamma_BDT.SetLineColor(kGreen)

        # Draw histograms
        hSD_BDT.Draw("HIST")
        hSDstarGamma_BDT.Draw("HISTsame")
        hSDstarPi0_BDT.Draw("HISTsame")
        hB_BDT.Draw("HISTsame")
        c.Draw()
    elif variable == "Xps":
        c.cd(1)
        # Customize histogram attributes
        hSD_Xps.SetLineColor(kBlue)
        hSD_Xps.SetFillColor(kBlue)
        hSD_Xps.SetFillStyle(3001)  # Diagonal lines
      
        hSDstarGamma_Xps.SetLineColor(kRed)
        hSDstarGamma_Xps.SetFillColorAlpha(kRed, 0.5);  # Transparent fill
        hSDstarGamma_Xps.SetMarkerStyle(21)
        hSDstarGamma_Xps.SetMarkerSize(0.7)
      
        hSDstarGamma_Xps.SetLineColor(kGreen)

        # Draw histograms
        hSD_Xps.Draw("HIST")
        hSDstarGamma_Xps.Draw("HISTsame")
        hSDstarPi0_Xps.Draw("HISTsame")
        hB_Xps.Draw("HISTsame")
        c.Draw()
    elif variable == "q2":
        c.cd(1)
        # Customize histogram attributes
        hSD_q2.SetLineColor(kBlue)
        hSD_q2.SetFillColor(kBlue)
        hSD_q2.SetFillStyle(3001)  # Diagonal lines
      
        hSDstarGamma_q2.SetLineColor(kRed)
        hSDstarGamma_q2.SetFillColorAlpha(kRed, 0.5);  # Transparent fill
        hSDstarGamma_q2.SetMarkerStyle(21)
        hSDstarGamma_q2.SetMarkerSize(0.7)
      
        hSDstarGamma_q2.SetLineColor(kGreen)

        # Draw histograms
        hSD_q2.Draw("HIST")
        hSDstarGamma_q2.Draw("HISTsame")
        hSDstarPi0_q2.Draw("HISTsame")
        hB_q2.Draw("HISTsame")
        c.Draw()
    elif variable == "Mmin2":
        c.cd(1)
        # Customize histogram attributes
        hSD_Mmin2.SetLineColor(kBlue)
        hSD_Mmin2.SetFillColor(kBlue)
        hSD_Mmin2.SetFillStyle(3001)  # Diagonal lines
      
        hSDstarGamma_Mmin2.SetLineColor(kRed)
        hSDstarGamma_Mmin2.SetFillColorAlpha(kRed, 0.5);  # Transparent fill
        hSDstarGamma_Mmin2.SetMarkerStyle(21)
        hSDstarGamma_Mmin2.SetMarkerSize(0.7)
      
        hSDstarGamma_Mmin2.SetLineColor(kGreen)

        # Draw histograms
        hSD_Mmin2.Draw("HIST")
        hSDstarGamma_Mmin2.Draw("HISTsame")
        hSDstarPi0_Mmin2.Draw("HISTsame")
        hB_Mmin2.Draw("HISTsame")
        c.Draw()




# Creating pdfs and pseudo data
def create_data_pdfs(lum, mass, variable):
    """
    Create data PDFs from histograms.
    
    Parameters:
        histograms (dict): Dictionary of histograms.
        lum (float): Luminosity value.
        variable (str): Variable name.
    
    Returns:
        tuple: Contains the generated data and PDFs.
    """
    # Initialize variables
    dSD = None
    dSDstarGamma = None
    dSDstarPi0 = None
    dB = None
    hist_dataB = None
    s1_pdf = None
    s2_pdf = None
    s3_pdf = None
    b_pdf = None
    hist_dataB = RooDataHist()
 
    # Define histograms and variables based on 'variable'
    if variable == "BDT":
        # Number of events
        NsD = hSD_BDT.GetEntries()
        NsDstarGamma = hSDstarGamma_BDT.GetEntries()
        NsDstarPi0 = hSDstarPi0_BDT.GetEntries()
        Nb = hB_BDT.GetEntries()

        # Set the variable and create the dataset
        x = ROOT.RooRealVar("x", "x", BDTLow, BDTUp)
        x.setBins(nbins)

        # Create RooDataHist
        dSD = ROOT.RooDataHist("dSD", "dSD", ROOT.RooArgList(x), ROOT.RooFit.Import(hSD_BDT))
        dSDstarGamma = ROOT.RooDataHist("dSDstarGamma", "dSDstarGamma", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarGamma_BDT))
        dSDstarPi0 = ROOT.RooDataHist("dSDstarPi0", "dSDstarPi0", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarPi0_BDT))
        dB = ROOT.RooDataHist("dB", "dB", ROOT.RooArgList(x), ROOT.RooFit.Import(hB_BDT))

        Nb = lum * Nb / Luminosity

        s1_pdf = ROOT.RooHistPdf("s1_pdf", "signal D tag pdf", ROOT.RooArgSet(x), dSD, 2)
        s2_pdf = ROOT.RooHistPdf("s2_pdf", "signal DstarGamma tag pdf", ROOT.RooArgSet(x), dSDstarGamma, 2)
        s3_pdf = ROOT.RooHistPdf("s3_pdf", "signal DstarPi0 tag pdf", ROOT.RooArgSet(x), dSDstarPi0, 2)
        b_pdf = ROOT.RooHistPdf("b_pdf", "bkg pdf", ROOT.RooArgSet(x), dB, 2)

    elif variable == "Xps":
        # Define histogram names and binning for Xps
        # Assuming hSD_Xps, hSDstarGamma_Xps, hSDstarPi0_Xps, hB_Xps are defined
        NsD = hSD_Xps.GetEntries()
        NsDstarGamma = hSDstarGamma_Xps.GetEntries()
        NsDstarPi0 = hSDstarPi0_Xps.GetEntries()
        Nb = hB_Xps.GetEntries()

        x = ROOT.RooRealVar("x", "x", XpsLow, XpsUp)
        x.setBins(nbins)

        dSD = ROOT.RooDataHist("dSD", "dSD", ROOT.RooArgList(x), ROOT.RooFit.Import(hSD_Xps))
        dSDstarGamma = ROOT.RooDataHist("dSDstarGamma", "dSDstarGamma", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarGamma_Xps))
        dSDstarPi0 = ROOT.RooDataHist("dSDstarPi0", "dSDstarPi0", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarPi0_Xps))
        dB = ROOT.RooDataHist("dB", "dB", ROOT.RooArgList(x), hB_Xps)

        Nb = lum * Nb / Luminosity

        s1_pdf = ROOT.RooHistPdf("s1_pdf", "signal D tag pdf", ROOT.RooArgSet(x), dSD, 2)
        s2_pdf = ROOT.RooHistPdf("s2_pdf", "signal DstarGamma tag pdf", ROOT.RooArgSet(x), dSDstarGamma, 2)
        s3_pdf = ROOT.RooHistPdf("s3_pdf", "signal DstarPi0 tag pdf", ROOT.RooArgSet(x), dSDstarPi0, 2)
        b_pdf = ROOT.RooHistPdf("b_pdf", "bkg pdf", ROOT.RooArgSet(x), dB, 2)

    elif variable == "q2":
        # Define histogram names and binning for q2
        # Assuming hSD_q2, hSDstarGamma_q2, hSDstarPi0_q2, hB_q2 are defined
        NsD = hSD_q2.GetEntries()
        NsDstarGamma = hSDstarGamma_q2.GetEntries()
        NsDstarPi0 = hSDstarPi0_q2.GetEntries()
        Nb = hB_q2.GetEntries()

        x = ROOT.RooRealVar("x", "x", q2Low, q2Up)
        x.setBins(nbins)

        dSD = ROOT.RooDataHist("dSD", "dSD", ROOT.RooArgList(x), ROOT.RooFit.Import(hSD_q2))
        dSDstarGamma = ROOT.RooDataHist("dSDstarGamma", "dSDstarGamma", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarGamma_q2))
        dSDstarPi0 = ROOT.RooDataHist("dSDstarPi0", "dSDstarPi0", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarPi0_q2))
        dB = ROOT.RooDataHist("dB", "dB", ROOT.RooArgList(x), ROOT.RooFit.Import(hB_q2))

        Nb = lum * Nb / Luminosity

        s1_pdf = ROOT.RooHistPdf("s1_pdf", "signal D tag pdf", ROOT.RooArgSet(x), dSD, 2)
        s2_pdf = ROOT.RooHistPdf("s2_pdf", "signal DstarGamma tag pdf", ROOT.RooArgSet(x), dSDstarGamma, 2)
        s3_pdf = ROOT.RooHistPdf("s3_pdf", "signal DstarPi0 tag pdf", ROOT.RooArgSet(x), dSDstarPi0, 2)
        b_pdf = ROOT.RooHistPdf("b_pdf", "bkg pdf", ROOT.RooArgSet(x), dB, 2)

    elif variable == "Mmin2":
        # Define histogram names and binning for Mmin2
        # Assuming hSD_Mmin2, hSDstarGamma_Mmin2, hSDstarPi0_Mmin2, hB_Mmin2 are defined
        NsD = hSD_Mmin2.GetEntries()
        NsDstarGamma = hSDstarGamma_Mmin2.GetEntries()
        NsDstarPi0 = hSDstarPi0_Mmin2.GetEntries()
        Nb = hB_Mmin2.GetEntries()

        x = ROOT.RooRealVar("x", "x", Mmin2Low, Mmin2Up)
        x.setBins(nbins)

        dSD = ROOT.RooDataHist("dSD", "dSD", ROOT.RooArgList(x), ROOT.RooFit.Import(hSD_Mmin2))
        dSDstarGamma = ROOT.RooDataHist("dSDstarGamma", "dSDstarGamma", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarGamma_Mmin2))
        dSDstarPi0 = ROOT.RooDataHist("dSDstarPi0", "dSDstarPi0", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarPi0_Mmin2))
        dB = ROOT.RooDataHist("dB", "dB", ROOT.RooArgList(x), hB_Mmin2)

        Nb = lum * Nb / Luminosity

        s1_pdf = ROOT.RooHistPdf("s1_pdf", "signal D tag pdf", ROOT.RooArgSet(x), dSD, 2)
        s2_pdf = ROOT.RooHistPdf("s2_pdf", "signal DstarGamma tag pdf", ROOT.RooArgSet(x), dSDstarGamma, 2)
        s3_pdf = ROOT.RooHistPdf("s3_pdf", "signal DstarPi0 tag pdf", ROOT.RooArgSet(x), dSDstarPi0, 2)
        b_pdf = ROOT.RooHistPdf("b_pdf", "bkg pdf", ROOT.RooArgSet(x), dB, 2)

    elif variable == "Mmax2":
        # Define histogram names and binning for Mmax2
        # Assuming hSD_Mmax2, hSDstarGamma_Mmax2, hSDstarPi0_Mmax2, hB_Mmax2 are defined
        NsD = hSD_Mmax2.GetEntries()
        NsDstarGamma = hSDstarGamma_Mmax2.GetEntries()
        NsDstarPi0 = hSDstarPi0_Mmax2.GetEntries()
        Nb = hB_Mmax2.GetEntries()

        x = ROOT.RooRealVar("x", "x", Mmax2Low, Mmax2Up)
        x.setBins(nbins)

        dSD = ROOT.RooDataHist("dSD", "dSD", ROOT.RooArgList(x), ROOT.RooFit.Import(hSD_Mmax2))
        dSDstarGamma = ROOT.RooDataHist("dSDstarGamma", "dSDstarGamma", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarGamma_Mmax2))
        dSDstarPi0 = ROOT.RooDataHist("dSDstarPi0", "dSDstarPi0", ROOT.RooArgList(x), ROOT.RooFit.Import(hSDstarPi0_Mmax2))
        dB = ROOT.RooDataHist("dB", "dB", ROOT.RooArgList(x), hB_Mmax2)

        Nb = lum * Nb / Luminosity

        s1_pdf = ROOT.RooHistPdf("s1_pdf", "signal D tag pdf", ROOT.RooArgSet(x), dSD, 2)
        s2_pdf = ROOT.RooHistPdf("s2_pdf", "signal DstarGamma tag pdf", ROOT.RooArgSet(x), dSDstarGamma, 2)
        s3_pdf = ROOT.RooHistPdf("s3_pdf", "signal DstarPi0 tag pdf", ROOT.RooArgSet(x), dSDstarPi0, 2)
        b_pdf = ROOT.RooHistPdf("b_pdf", "bkg pdf", ROOT.RooArgSet(x), dB, 2)

    elif variable == "2D":
        # Define histogram names and binning for 2D
        # Assuming hSD_2D, hSDstarGamma_2D, hSDstarPi0_2D, hB_2D are defined
        NsD = hSD_2D.GetEntries()
        NsDstarGamma = hSDstarGamma_2D.GetEntries()
        NsDstarPi0 = hSDstarPi0_2D.GetEntries()
        Nb = hB_2D.GetEntries()

        nbinsX = hB_2D.GetNbinsX()
        nbinsY = hB_2D.GetNbinsY()
        # Filling empty bins with 1 at least
        for ix in range(nbinsX):
            for iy in range(nbinsY):
                bincontent = hB_2D.GetBinContent(ix, iy)
                if bincontent == 0:
                    hB_2D.SetBinContent(ix, iy, 1)
        
        x = ROOT.RooRealVar("x", "x", X_min, X_max)
        y = ROOT.RooRealVar("y", "y", Y_min, Y_max)

        x.setBins(hB_2D.GetNbinsX())
        y.setBins(hB_2D.GetNbinsY())

        dSD = ROOT.RooDataHist("dSD", "dSD", ROOT.RooArgList(x, y), ROOT.RooFit.Import(hSD_2D))
        dSDstarGamma = ROOT.RooDataHist("dSDstarGamma", "dSDstarGamma", ROOT.RooArgList(x, y), ROOT.RooFit.Import(hSDstarGamma_2D))
        dSDstarPi0 = ROOT.RooDataHist("dSDstarPi0", "dSDstarPi0", ROOT.RooArgList(x, y), ROOT.RooFit.Import(hSDstarPi0_2D))
        dB = ROOT.RooDataHist("dB", "dB", ROOT.RooArgList(x, y), hB_2D)

        Nb = lum * Nb / Luminosity

        s1_pdf = ROOT.RooHistPdf("s1_pdf", "signal D tag pdf", ROOT.RooArgSet(x, y), dSD, 2)
        s2_pdf = ROOT.RooHistPdf("s2_pdf", "signal DstarGamma tag pdf", ROOT.RooArgSet(x, y), dSDstarGamma, 2)
        s3_pdf = ROOT.RooHistPdf("s3_pdf", "signal DstarPi0 tag pdf", ROOT.RooArgSet(x, y), dSDstarPi0, 2)
        b_pdf = ROOT.RooHistPdf("b_pdf", "bkg pdf", ROOT.RooArgSet(x, y), dB, 2)

    else:
        raise ValueError(f"Variable '{variable}' is not recognized.")

    # Let's generate pseudo-data
    if variable in ["BDT", "Xps", "q2", "Mmin2", "Mmax2"]:
        hist_data = generate_pseudo_data(variable, b_pdf, Nb, x)
    else:
        hist_data = generate_pseudo_data(variable, b_pdf, Nb, x, y)

    w = setup_workspace()

    # Import the pdf
    w.Import(s1_pdf)
    w.Import(s2_pdf)
    w.Import(s3_pdf)
    w.Import(b_pdf)

    #return
    if variable in ["BDT", "Xps", "q2", "Mmin2", "Mmax2"]:
        return w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb
    else:
        return w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb, x, y


def main(lum, mass, variable):
    # Set up workspace
    w = setup_workspace()
    # Set up data
    set_data(mass, variable)
    # Create the pdfs based on templates
    if variable in ["BDT", "Xps", "q2", "Mmin2", "Mmax2"]:
        w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb = create_data_pdfs(lum, mass, variable)
        print("Pseudo data created ...")

        # Create the model in the workspace
        create_model(w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb, lum, variable)

        # Fit the model with data
        fit_model(w, hist_data)
        print("Model fitted")

        get_limit(w, hist_data, variable)
    else:
        w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb, x, y = create_data_pdfs(lum, mass, variable)
        print("Pseudo data created ...")

        # Create the model in the workspace
        create_model(w, hist_data, NsD, NsDstarGamma, NsDstarPi0, Nb, lum, variable)

        # Fit the model with data
        fit_model(w, hist_data)
        print("Model fitted")

        get_limit(w, hist_data, variable, x, y)



import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--luminosity', type=float, help='Luminosity', required=True)
    parser.add_argument('-m', '--mass', type=str, help='Mass', required=True)
    parser.add_argument('-v', '--variable', type=str, help='Variable', required=True)
    
    args = parser.parse_args()
    main(lum = args.luminosity, mass = args.mass, variable =args.variable)