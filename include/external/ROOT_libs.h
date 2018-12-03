/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/
 
 /*-------------------------------------
 * Include classes from ROOT to fill
 * and draw histograms and plots.
 *-------------------------------------
 */
#include <TROOT.h>
#include <TLegend.h>
#include <TCut.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TSystemFile.h>
#include <TChain.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TLine.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TString.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TDirectory.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooCBShape.h"
#include "RooDstD0BG.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"

//Minuit2
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include "Minuit2/VariableMetricMinimizer.h"


