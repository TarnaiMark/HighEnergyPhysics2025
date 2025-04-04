#define particle_tree_cxx
#include "particle_tree.h"
#include <iostream>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <numeric>  
#include <complex>
#include <iostream> 
#include <fstream>  
#include "TGraphErrors.h"
#include <vector>
#include "TMath.h"
#include "TF1.h"
using namespace std;

const int NCENT = 5;  //number of centrality bins
const int NPT = 18;   //number of pt bins
//centrality classes:
std::vector<double> cent_bins = {0.0, 10.0, 20.0, 30.0 ,40.0, 50.0};
//pt classes:
std::vector<double> pt_bins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};


TH1F *phiDist[NCENT][NPT];

int findBinIndex(const std::vector<double>& binLimits, int numBins, double value) {
  if (numBins <= 0 || binLimits.size() != static_cast<size_t>(numBins + 1)) {
      throw std::invalid_argument("Invalid bin limits or number of bins.");
  }
  for (int i = 0; i < numBins; ++i) {
      if (value >= binLimits[i] && value < binLimits[i + 1]) {
          return i;
      }
  }
  if (value == binLimits[numBins]) {
      return numBins - 1;
  }
  return -1; 
}

double pt(double px, double py) {
  return sqrt(px*px + py*py);
}



int main(int argc, const char** argv)
{
  string infilename = "/project/femtoscopy/data/PHENIX_AuAu_200GeV/measure_createtree_01234.root";
  string outfilename = "out/result.root"; 

  int Nmaxevt = -1;
  if(argc==2) Nmaxevt = atoi(argv[1]);
  if(Nmaxevt<1) Nmaxevt = -1;
  
  //INITIALIZE particle_tree OBJECT
  particle_tree p(infilename.c_str());
  if(p.fChain) cerr << "Tree initialized" << endl;
  else { cerr << "No tree found." << endl; return 1; }

  //DETERMINE HOW MANY EVENTS TO RUN ON
  long int Nevents = p.fChain->GetEntries();
  if(Nmaxevt>0 && Nmaxevt<(int)Nevents) Nevents=Nmaxevt;
  cerr << "Will run on " << Nevents << " events (out of " << p.fChain->GetEntries()  << ")." << endl;

  //CREATE HISTOGRAMS
  for(int icent=0; icent<NCENT; icent++)
    for(int ipt=0; ipt<NPT; ipt++)
    {
      phiDist[icent][ipt] = new TH1F(Form("phidist_icent%i_ipt%i",icent,ipt),"Phi Distribution",200,-M_PI/2, M_PI/2);
    }
  //LOOP THROUGH EVENTS IN THE GIVEN DATASET
  for(int ievent=0;ievent<Nevents;ievent++)
  {
    //MONITOR PROGRESS THROUGH STDERR OUTPUT
    if(ievent>0&&ievent%1000==0) cerr << ".";
    if(ievent>0&&ievent%10000==0) cerr <<", Analyzing event #" << ievent << endl;

    //LOAD EVENT DATA INTO particle_tree OBJECT
    p.GetEntry(ievent);
    int cent = p.Centrality;  //centrality of the current event
    
    //Determining which centrality bin every track from the event belongs to

    int cent_bin = findBinIndex(cent_bins, NCENT, cent);
    if (cent_bin == -1) {
      continue;
    }


    //LOOP THROUGH TRACKS
    for(int ipart=0;ipart<p.Ntracks;ipart++)
    {
      //determining pt bin of track
      double pt_track = pt(p.px[ipart], p.py[ipart]);
      int pt_bin = findBinIndex(pt_bins, NPT, pt_track);
      if (pt_bin == -1) {
        continue;
      }


      //CALCULATE phiR.
      double phi0 = atan2(p.py[ipart],p.px[ipart]);
      double phiR = phi0 - p.ReactionPlane;
      if (phiR < -M_PI) phiR += 2*M_PI;
      if (phiR > M_PI) phiR -= 2*M_PI;
      if (phiR < -M_PI/2) phiR += M_PI;
      if (phiR > M_PI/2) phiR -= M_PI;
      
      //FILL HISTOGRAMS
      
      phiDist[cent_bin][pt_bin]->Fill(phiR);
      
    }
  }
  
  //WRITE OUT HISTOGRAMS TO ROOT FILE
  TFile *f = new TFile(outfilename.c_str(),"RECREATE");
  if(!f->IsWritable()) cerr << "File " << outfilename.c_str() << " was not opened!" << endl;
  else cerr << "Analysis done, writing histos to " << outfilename.c_str() << endl;
  f->cd();
  for(int icent=0; icent<NCENT; icent++)
    for(int ipt=0; ipt<NPT; ipt++)
      phiDist[icent][ipt]->Write();
  f->Close();

  return 0; 
}


