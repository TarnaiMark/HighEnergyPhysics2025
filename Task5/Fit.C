using namespace std;
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TString.h>
#include <vector>
#include <TColor.h> 
#include <TLegend.h>

const int NCENT = 5;  //number of centrality bins
const int NPT = 18;   //number of pt bins
TGraphErrors *v2_graph[NCENT];
TH1F *phiDist[NCENT][NPT];
double fitFunc(double* x, double* par);
//centrality classes:
std::vector<double> cent_bins = {0.0, 10.0, 20.0, 30.0 ,40.0, 50.0};
//pt classes:
std::vector<double> pt_bins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
double R_value;


double fitFunc(double* x, double* par)
{
  double a = par[0];
  double b = par[1]; 
  double phi = x[0];
  return a + b * cos(2 * phi);
}

double R(double cent)
{
  double R = 1e-9 * TMath::Power(cent, 5) - 2.335e-7 * TMath::Power(cent, 4) + 2.4645e-05 * TMath::Power(cent, 3)
   - 0.001539 * TMath::Power(cent, 2) +  0.042389 * cent + 0.36115;
  return R;
}

void Fit()
{
  //READ THE INPUT FILE
  TFile *inputfile = new TFile("out/result.root");
  for(int icent=0; icent<NCENT; icent++)
    for(int ipt=0; ipt<NPT; ipt++)
    {
      phiDist[icent][ipt] = (TH1F*)inputfile->Get(Form("phidist_icent%i_ipt%i",icent,ipt));
    }
  
  //CREATING GRAPHS TO STORE V2 VALUES
  for(int icent=0; icent<NCENT; icent++)
  {
    v2_graph[icent] = new TGraphErrors();
    v2_graph[icent]->SetName(Form("v2_graph_icent%i",icent));
  }

  //FIT THE HISTOGRAMS, CALCULATE V2 AND ITS ERRORS
  TF1* f_fitFunc = new TF1("f_fitFunc", fitFunc, -M_PI/2, M_PI/2, 2);
  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);


  for(int icent=0; icent<NCENT; icent++){
    //Calculating R for each centrality bin
    R_value = R(cent_bins[icent] + 5.0);

    for(int ipt=0; ipt<NPT; ipt++)
    {
      //Fitting the histogram
      phiDist[icent][ipt]->Fit("f_fitFunc","","",-M_PI/2,M_PI/2);

      //Calculating v2 and its error
      double a = f_fitFunc->GetParameter(0);
      double b = f_fitFunc->GetParameter(1);
      double a_err = f_fitFunc->GetParError(0);
      double b_err = f_fitFunc->GetParError(1);


      double v2 = b / (2*a* R_value);
      double v2_err = v2 * sqrt(pow(a_err/a, 2) + pow(b_err/b, 2));

      //Filling the graph
      v2_graph[icent]->SetPoint(ipt, pt_bins[ipt], v2);
      v2_graph[icent]->SetPointError(ipt, 0, v2_err);


    
      //PLOT THE HISTOS AND THE FITS, with displaying the fit parameters and their errors in the legend
      gStyle->SetOptStat(0);
      phiDist[icent][ipt]->GetXaxis()->SetTitle("#phi");
      phiDist[icent][ipt]->GetYaxis()->SetTitle("Counts");
      phiDist[icent][ipt]->SetTitle(Form("Centrality: %i-%i%%, p_{T}: %.1f-%.1f GeV/c", (int)cent_bins[icent], (int)cent_bins[icent+1], pt_bins[ipt], pt_bins[ipt+1]));
      // Create a legend
      TLegend *legend = new TLegend(0.15, 0.75, 0.45, 0.9);
      legend->SetBorderSize(1);
      legend->SetFillColor(0);
      // Add parameters with errors
      legend->AddEntry((TObject*)0, Form("a = %.3f #pm %.3f", a, a_err), "");
      legend->AddEntry((TObject*)0, Form("b = %.3f #pm %.3f", b, b_err), "");
      
      phiDist[icent][ipt]->Draw("e");
      legend->Draw();
      c1->Print(Form("figs/phidist_icent%i_ipt%i_.pdf",icent,ipt));
      c1->Clear();
    }
  }
  //PLOT THE RESULTS, every centrality bin in a different color
  TLegend *legend = new TLegend(0.15, 0.75, 0.45, 0.9);
  for(int icent=0; icent<NCENT; icent++)
  {
    
    if (icent == 0 ){
      
      v2_graph[icent]->SetMarkerStyle(20);
      v2_graph[icent]->SetMarkerSize(1);
      v2_graph[icent]->SetMarkerColor(TColor::GetColorPalette(icent*50));
      v2_graph[icent]->SetLineColor(TColor::GetColorPalette(icent*50));
      v2_graph[icent]->GetHistogram()->GetYaxis()->SetRangeUser(0, 0.25);  // Set Y-axis range
      
      v2_graph[icent]->Draw("AP");
      
    } 
    else{
      
      v2_graph[icent]->SetMarkerStyle(20);
      v2_graph[icent]->SetMarkerSize(1);
      v2_graph[icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      v2_graph[icent]->GetYaxis()->SetTitle("v_{2}");
      v2_graph[icent]->SetTitle("v_{2} vs p_{T}");
      v2_graph[icent]->SetMarkerColor(TColor::GetColorPalette(icent*50));
      v2_graph[icent]->SetLineColor(TColor::GetColorPalette(icent*50));
      v2_graph[icent]->Draw("P SAME");
      
    }
    
  }
  for(int icent=0; icent<NCENT; icent++)
  {
    legend->AddEntry(v2_graph[icent], Form("%i-%i%%", (int)cent_bins[icent], (int)cent_bins[icent+1]), "ep");
  }
  

  legend->Draw();
  c1->Print("figs/v2_graphs.pdf");
  delete c1;
  inputfile->Close();
  delete inputfile;
}
