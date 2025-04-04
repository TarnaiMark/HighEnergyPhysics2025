#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>




int barWidth = 70;
void progressbar(Int_t event, Int_t nEvents)
{
  /*
     Simple progressbar
  */
  float progress = event*1.0 / nEvents * 1.0;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << "% ( " << event << " / " << nEvents << " )" << "\r";
  std::cout.flush();
  
}

double BB_GausfitMean(const double LogP)
{
  /*
    Expected values of negative pion dEdx
  */
  double BB_GFitMeanpar1 = 3.10453;
  double BB_GFitMeanpar2 = -5.50479;
  double BB_GFitMeanpar3 = 3.88143;
  double BB_GFitMeanpar4 = -0.477965;
  return BB_GFitMeanpar1 + BB_GFitMeanpar2 / sqrt(LogP+2.) + BB_GFitMeanpar3 / (LogP+2.) + BB_GFitMeanpar4 / (LogP+2.)/(LogP+2.);
}


double BB_GausfitSigma(const double LogP)
{
  /*
     1 sigma deviation from the expected values of negative pion dEdx
  */
  double BB_GFitSigmapar1 = 0.348025;
  double BB_GFitSigmapar2 = -1.3454;
  double BB_GFitSigmapar3 = 1.61325;
  double BB_GFitSigmapar4 = -0.52485;

  return BB_GFitSigmapar1 + BB_GFitSigmapar2 / sqrt(LogP+2.) + BB_GFitSigmapar3 / (LogP+2.) + BB_GFitSigmapar4 / (LogP+2.)/(LogP+2.);
}




void na61dEdx()
{
  //loading input file and setting output file 
  TFile *inputFile;
  inputFile = new TFile("/project/femtoscopy/data/NA61_data/ArSc_150_DOKK_DataLike.root","read");
  if(inputFile->IsZombie()) {
    std::cout << "File " << " not opened" << std::endl;
    return;
  }
  
  TFile* outFile = new TFile("out/output.root","RECREATE");
  
  //reading event tree from the input file
  TTree *event = (TTree*) inputFile->Get("ArSc_simpleTree");

  //create histograms containers and set their attributes
  TH1D *dhP = new TH1D("dhP", "Momentum;p [GeV/c];Entries", 500, 0.0, 50.0);
  TH2D *dhdEdx = new TH2D("dhdEdx", "particle identification; log(p/(1 GeV/c)); dEdx", 500, -5.0, 5.0,500, 0, 15.0);
  TH2D *dhselecteddEdx = new TH2D("dhselecteddEdx", "particle identification; log(p/(1 GeV/c)); dEdx", 500, -5.0, 5.0,500, 0, 15.0);
  

  // Setting up variables to be read in from input file
  int evID = 0;
  int nTrk = 0;
  double xVtx = 0; 
  double yVtx = 0; 
  double zVtx = 0; 
  int numPi = 0;
  std::vector<int> *charge = 0; 
  std::vector<double> *PX = 0;
  std::vector<double> *PY = 0;
  std::vector<double> *PZ = 0;
  std::vector<double> *dEdx = 0;

  //Read into the variables from different branches
  event->SetBranchAddress("evID",&evID);
  event->SetBranchAddress("nTrk",&nTrk);
  event->SetBranchAddress("charge",&charge);
  event->SetBranchAddress("xVtx",&xVtx);
  event->SetBranchAddress("yVtx",&yVtx);
  event->SetBranchAddress("zVtx",&zVtx);
  event->SetBranchAddress("trkPx",&PX);
  event->SetBranchAddress("trkPy",&PY);
  event->SetBranchAddress("trkPz",&PZ);
  event->SetBranchAddress("dEdx",&dEdx);

  // Get number of events in the file and loop through
  Int_t nEvent = event->GetEntries();
  // If you want to run for fewer events uncomment
  //nEvent = 100000;
  for (Int_t ievent = 0; ievent < nEvent; ievent++) {
    event->GetEntry(ievent);

    // remove during testing
    progressbar(ievent+1.0, nEvent);

    int eventsize = nTrk;
    // Loop through particles of the given event
    for(int ipart = 0; ipart < eventsize; ipart++)
    {
      if(charge->at(ipart) == -1) {
        // values
        double dEdxValue = dEdx->at(ipart);
        double p = sqrt(PX->at(ipart)*PX->at(ipart) + PY->at(ipart)*PY->at(ipart) + PZ->at(ipart)*PZ->at(ipart));
        double logp = TMath::Log(p/1.0);

        // Fill histograms
        dhP->Fill(p);
        dhdEdx->Fill(logp, dEdxValue);

        // Select particles
        if(dEdxValue > BB_GausfitMean(logp) - 2.0 * BB_GausfitSigma(logp) && dEdxValue < BB_GausfitMean(logp) + 3.0 * BB_GausfitSigma(logp)) {
          dhselecteddEdx->Fill(logp, dEdxValue);
          numPi++; 
        }
      }
    }

  }
  std::cout << std::endl;

  std::stringstream title;
  title << "number of selected Pions:  " << numPi;
  dhselecteddEdx->SetTitle(title.str().c_str());


  // Open outputfile and write out histograms
  outFile->cd();

  dhP->Write();
  dhdEdx->Write();
  dhselecteddEdx->Write();
  

  outFile->Close();

  return;
}

