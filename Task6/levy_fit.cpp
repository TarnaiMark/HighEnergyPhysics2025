#include <iostream>
#include <vector>
#include <string>
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TFile.h"
#include "TStyle.h"

// To calculate the 1D Levy-function
#include "levy_calc.h"

using namespace std;

TGraphErrors* gr;

double rmin;
double rmax;
const int NPARS = 3;
int NDF;

// Levy function to fit
double FitFunction(const double *x, const double *par)
{
	double r = x[0];
	double N = par[0];
    double R = par[1];
    double alpha = par[2];
    double Rcc = R*pow(2.,1. / alpha);
    return N*levy_calc(r/Rcc,Rcc,alpha);
}

// Chi-square error function
double MyChi2(const double *par)
{
    double chi2 = 0;
    NDF = 0;
    for(int ix=0;ix<gr->GetN();ix++)
    {
        double r = gr->GetX()[ix];
        if(r<rmin || r>rmax) continue;
        double exp = gr->GetY()[ix];
        double theor = FitFunction(&r,par);
        double err = gr->GetEY()[ix];
        if(err==0) continue;
        double chi = (exp-theor)/err;
        chi2 += chi*chi;
        NDF++;
    }
    NDF -= NPARS; // NDF = (number of fitted data points) - (number of fitting parameters)
    return chi2;
}

// Fitting and plotting one rho distribution
vector<double> levy_fit(const string& root_filename) {
    // Opening the root file and read the histogram
    TFile *file = new TFile(root_filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open the file: " << root_filename << std::endl;
        return {};
    }
    TH1F *hist = (TH1F*)file->Get("rho_hist");
   
    // Creating TGraphErrors object from the histogram data
    int nPoints = hist->GetNbinsX();
    gr = new TGraphErrors(nPoints);
    for (int i = 1; i <= nPoints; i++) {
        double x = hist->GetBinCenter(i);
        double y = hist->GetBinContent(i);
        double err = hist->GetBinError(i);
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, 0, err);
    }

    // Set fitting and plotting boundaries
    rmin = 1.0;
    rmax = 200.0;
    double rminplot = 0.3;
    double rmaxplot = 10000;
    
    // Set up the minimization:
    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kCombined ); // Choose method upon creation between: kMigrad, kSimplex, kCombined, kScan, kFumili
    
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(0.001);
    min.SetPrintLevel(2);
    
    ROOT::Math::Functor f(&MyChi2,NPARS);
    min.SetFunction(f);
    
    min.SetLimitedVariable(0,"N",      1.0 ,0.01, 0.4, 1.0);
    min.SetLimitedVariable(1,"R",      8.0 ,0.01, 3.0, 15.0);
    min.SetLimitedVariable(2,"alpha",  1.1 ,0.01, 0.8, 1.8);
    
    // Actual minimization
    min.Minimize();
    const double *par = min.X();
    const double *err = min.Errors();
    double chi2 = MyChi2(par);
    
    // Plotting the results
    /* Plot the results of the fit. Show the fitted parameters and their errors on the plot with TLatex) */
    TCanvas *c = new TCanvas("c","c",800,600);
    c->SetLogx();
    c->SetLogy();

    gr->Draw("AP");

    //plotting the result of the fit on the same plot
    TF1 *fit = new TF1("fit",FitFunction,rminplot,rmaxplot,NPARS);
    fit->SetParameters(par);
    fit->SetLineColor(kRed);
    
    fit->Draw("same");




    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.2,0.85,Form("N = %.3f #pm %.3f",par[0],err[0]));
    tex->DrawLatex(0.2,0.80,Form("R = %.3f #pm %.3f",par[1],err[1]));
    tex->DrawLatex(0.2,0.75,Form("#alpha = %.3f #pm %.3f",par[2],err[2]));
    tex->DrawLatex(0.2,0.70,Form("#chi^{2}/NDF = %.3f / %d",chi2,NDF));

    //saving all 10 graphs separatly
    std::string number = root_filename.substr(root_filename.find("my59a-")+6,1);
    std::string output_name = "out/levy_fit_" + number + ".pdf";

    c->Print(output_name.c_str());

    c->Close();
    return {par[0], err[0], par[1], err[1], par[2], err[2]};
}

int main() {
    // D_rho histo file names with centralities
    std::vector<std::string> root_files = {
        "histos/z-gg2my59a-1-DrhoHist.root",
        "histos/z-gg2my59a-2-DrhoHist.root",
        "histos/z-gg2my59a-3-DrhoHist.root",
        "histos/z-gg2my59a-4-DrhoHist.root",
        "histos/z-gg2my59a-5-DrhoHist.root",
        "histos/z-gg2my59a-6-DrhoHist.root",
        "histos/z-gg2my59a-7-DrhoHist.root",
        "histos/z-gg2my59a-8-DrhoHist.root",
        "histos/z-gg2my59a-9-DrhoHist.root",
        "histos/z-gg2my59a-10-DrhoHist.root"
    };

    vector<double> N_part = {351.4, 299, 253.9, 215.3, 181.6 , 151.5, 125.7, 102.7, 82.9, 65.9};
    
    vector<double> N_vals, N_errs, R_vals, R_errs, alpha_vals, alpha_errs, chi2_probs;
    
    // Iterating through the histograms and perform the fitting on all of them
    for (size_t i = 0; i < root_files.size(); i++) {
        vector<double> params = levy_fit(root_files[i]);
        if (params.empty()) continue;
        
        N_vals.push_back(params[0]);
        N_errs.push_back(params[1]);
        R_vals.push_back(params[2]);
        R_errs.push_back(params[3]);
        alpha_vals.push_back(params[4]);
        alpha_errs.push_back(params[5]);
    }
    
    // Plotting the parameters
    /* Plot the fitted parameters with their errors in function of <N_part> */

    TCanvas* c1 = new TCanvas("c1", "Error Bars", 800, 600);
    TGraphErrors* N_graph = new TGraphErrors(N_part.size(), N_part.data(), N_vals.data(), 0, N_errs.data());
    N_graph->SetTitle("N vs <N_{part}>; <N_{part}>; N");
    N_graph->SetMarkerStyle(20);
    N_graph->SetMarkerSize(1);
    N_graph->SetMarkerColor(kBlue);
    N_graph->GetYaxis()->SetRangeUser(0.5, 1.0);
    N_graph->Draw("AP");
    c1->Print("out/N_vs_Npart.pdf");

    TGraphErrors* R_graph = new TGraphErrors(N_part.size(), N_part.data(), R_vals.data(), 0, R_errs.data());
    R_graph->SetTitle("R vs <N_{part}>; <N_{part}>; R");
    R_graph->SetMarkerStyle(20);
    R_graph->SetMarkerSize(1);
    R_graph->SetMarkerColor(kGreen);
    R_graph->GetYaxis()->SetRangeUser(5, 15);
    R_graph->Draw("AP");
    c1->Print("out/R_vs_Npart.pdf");

    TGraphErrors* alpha_graph = new TGraphErrors(N_part.size(), N_part.data(), alpha_vals.data(), 0, alpha_errs.data());
    alpha_graph->SetTitle("#alpha vs <N_{part}>; <N_{part}>; #alpha");
    alpha_graph->SetMarkerStyle(20);
    alpha_graph->SetMarkerSize(1);
    alpha_graph->SetMarkerColor(kRed);
    alpha_graph->GetYaxis()->SetRangeUser(0.5, 1.3);
    alpha_graph->Draw("AP");
    c1->Print("out/alpha_vs_Npart.pdf");


    c1->Close();

    
    return 0;
}

