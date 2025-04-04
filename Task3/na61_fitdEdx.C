#include <TH1D.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <TLegend.h>

#include <stdlib.h>


#include <iostream>
#include <cstring>


void na61_fitdEdx()
{
  //loading input file and setting output file 
  TFile *inputFile;
  inputFile = new TFile("out/output.root","read");
  if(inputFile->IsZombie()) {
    std::cout << "File " << " not opened" << std::endl;
    return;
  }
  TString inputFileName = inputFile->GetName(); 
  inputFileName.Remove(inputFileName.Sizeof()-6,inputFileName.Sizeof());
  inputFileName.Remove(0,4);
  

  
  // load histograms into containers and change their graphics

  TH2D *dhdEdx = (TH2D*)inputFile->Get("dhdEdx");
  TH1D *dhP = (TH1D*)inputFile->Get("dhP");
  TH2D *dhselecteddEdx = (TH2D*)inputFile->Get("dhselecteddEdx");
  

  
  //Projecting of selected particles
  int bin_x_min =  dhdEdx ->GetXaxis()->FindBin(-0.6);
  int bin_x_max =  dhdEdx ->GetXaxis()->FindBin(-0.4);
  TH1D *projecteddEdx = dhdEdx->ProjectionY("projectedSelected", bin_x_min, bin_x_max);
  
  
  // fitting two gaussians on the two peaks
  TF1 *f1 = new TF1("f1","gaus",0,1.3); 
  TF1 *f2 = new TF1("f2","gaus",1.4,2);
  f1->SetParameters(100, 2, 1);
  f2->SetParameters(100, 8, 1);
  f1->SetParNames("Constant","Mean","Sigma");
  f2->SetParNames("Constant","Mean","Sigma");
  f1->SetLineColor(2);
  f2->SetLineColor(3);
  f1->SetLineStyle(2);
  f2->SetLineStyle(2);
  f1->SetLineWidth(1);
  f2->SetLineWidth(1);
  projecteddEdx->Fit("f1","R");
  projecteddEdx->Fit("f2","R+");

  // plot the histogram
  TCanvas *c1 = new TCanvas("c1","c1", 600,450);
  gStyle->SetOptStat(0);
  // dEdx
  c1->SetLogz(1);
  dhdEdx->Draw("colz");
  dhdEdx->GetYaxis()->SetTitle("dEdx (arb. units)");
  dhdEdx->GetXaxis()->SetTitle("log(p/(1 GeV/c))");
  c1->Print("figs/dEdx.png");
  c1->Clear();

  // p  
  dhP->Draw();
  dhP->GetXaxis()->SetTitle("p [GeV/c]");
  c1->Print("figs/p.png");
  c1->Clear();
  
  // dEdx of selected particles
  dhselecteddEdx->Draw("colz");
  dhselecteddEdx->GetYaxis()->SetTitle("dEdx (arb. units)");
  c1->Print("figs/selecteddEdx.png");

  // projected dEdx
  c1 -> SetLogy(1);
  projecteddEdx->Draw();
  projecteddEdx->GetXaxis()->SetTitle("dEdx (arb. units)");
  projecteddEdx->GetYaxis()->SetTitle("Entries");
  projecteddEdx->SetTitle("Projected dEdx for selected particles");
  projecteddEdx->GetXaxis()->SetRangeUser(0.2, 4);

  // legend in the upper right corner
  TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);
  legend->AddEntry(f1,Form("Mean, Sigm:  %.2f , %.2f",f1->GetParameter(1),f1->GetParameter(2)),"l");
  legend->AddEntry(f2,Form("Mean, Sigm: %.2f,  %.2f",f2->GetParameter(1),f2->GetParameter(2)),"l");
  legend->Draw();

  c1->Print("figs/projecteddEdx.png");
  

  return;
}
