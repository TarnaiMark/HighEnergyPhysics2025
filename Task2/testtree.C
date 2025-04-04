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

#define MAXTRACKS 1000

using namespace std;

void testtree(const char* infilename)
{
    //Open file and create a TChain
    ifstream infile(infilename);
    if (!infile) {
        cerr << "Error: Could not open input file " << infilename << endl;
    }

    //Declare TChain
    TChain fChain("demo/TrackTree");

    // Check if input is a .root file or a list of ROOT files
    string filename(infilename);
    if (filename.size() >= 5 && filename.substr(filename.size()-5) == ".root") {
        TFile* f = new TFile(filename.c_str(), "READ");
        TTree* tree = (TTree*)f->Get("TrackTree");
        fChain.Add(filename.c_str());  // Add single file to TChain
        cerr << "Initialized TTree from " << filename << endl;
    }
    else {
        string rootfilename;    
        while (infile >> rootfilename) {
            fChain.Add(rootfilename.c_str());
        }
        cerr << "Initialized TChain with multiple ROOT files." << endl;
    } 
    cerr << "Initiaized input file from " << infilename << endl;

    // Define event necessary variables
    Int_t nTrk;
    Float_t zVtx;
    Float_t EP2MbyPi;
    
    //...

    //Define track variables
    Float_t trkPt[MAXTRACKS];
    Int_t trkCharge[MAXTRACKS];
    Float_t trkEta[MAXTRACKS];


    // Set branch addresses
    
    
    fChain.SetBranchAddress("nTrk", &nTrk);
    fChain.SetBranchAddress("zVtx", &zVtx);
    fChain.SetBranchAddress("EP3MbyPi", &EP2MbyPi);
    

    fChain.SetBranchAddress("trkPt", &trkPt);
    fChain.SetBranchAddress("trkCharge", &trkCharge);
    fChain.SetBranchAddress("trkEta", &trkEta);
    //...
    
    // Create histograms
    TH1D* zvtxdist = new TH1D("zvtxdist", "zVtx distribution", 140, -35, 35);
    TH1D* ntrkdist = new TH1D("ntrkdist", "Track multiplicity distribution", 1000, 1, 1000);
    TH1D* eventPlane2M = new TH1D("eventPlane2M", "Event Plane", 400, -3.14, 3.14);
    

    TH1D* trkPtHist = new TH1D("trkPtHist", "Track pT distribution", 100, 0, 4);
    TH1D* trkChargeHist = new TH1D("trkChargeHist", "Track charge distribution", 300, -1.5, 1.5);
    TH1D* trkEtaHist = new TH1D("trkEtaHist", "Track eta distribution", 400, -4, 4);
    


    //...

    // Get number of events
    Long64_t Nevents = fChain.GetEntries();
    cerr << "Processing  " << Nevents << " events..." << endl;

    // Loop through events
    for(Long64_t ievent = 0; ievent < Nevents; ievent++){
        if(ievent % 10000 == 0) cerr << ".";
        if(ievent % 100000 == 0) cerr << "Event " << ievent << endl;

        fChain.GetEntry(ievent);
        zvtxdist->Fill(zVtx);
        ntrkdist->Fill(nTrk);
        eventPlane2M->Fill(EP2MbyPi * 3.14);
        
        //loop through tracks
        for (int i = 0; i < nTrk; i++) {
            trkPtHist->Fill(trkPt[i]);
            trkChargeHist->Fill(trkCharge[i]);
            trkEtaHist->Fill(trkEta[i]);
        }
        
        //...   
    }

    std::cerr << "\nFinished processing events." << std::endl;

    // Vertex distribution
    TCanvas* canvas = new TCanvas();
    canvas->SetLogy();  // Enable log scale on Y-axis
    zvtxdist ->SetStats(0);  // Disable statistics box
    zvtxdist ->GetXaxis()->SetTitle("Number of tracks");
    zvtxdist->GetYaxis()->SetTitle("Number of events");
    zvtxdist->Draw();
    canvas->Print("zvtxdist.png");
    
    
    // Track multiplicity distribution
    TCanvas* canvas2 = new TCanvas();
    canvas2->SetLogy();  // Enable log scale on Y-axis
    canvas2->SetLogx();  // Enable log scale on X-axis
    ntrkdist->SetStats(0);  // Disable statistics box
    ntrkdist->GetXaxis()->SetTitle("Number of tracks");
    ntrkdist->GetYaxis()->SetTitle("Number of events");
    ntrkdist->Draw();
    canvas2->Print("ntrkdist.png");
    

    // Event Plane M distribution
    TCanvas* canvas3 = new TCanvas();  
    eventPlane2M->SetStats(0);  // Disable statistics box
    eventPlane2M->GetXaxis()->SetTitle("Event Plane M");
    eventPlane2M->GetYaxis()->SetTitle("Number of events");
    eventPlane2M->Draw();
    canvas3->Print("eventPlane2M.png");


    // Track pT distribution
    TCanvas* canvas5 = new TCanvas();
    trkPtHist->SetStats(0);  // Disable statistics box
    trkPtHist->GetXaxis()->SetTitle("Track pT");
    trkPtHist->GetYaxis()->SetTitle("Number of tracks");
    canvas5->SetLogy();  // Enable log scale on Y-axis
    trkPtHist->Draw();
    canvas5->Print("trkPtHist.png");
    

    // Track charge distribution
    TCanvas* canvas6 = new TCanvas();
    trkChargeHist->SetStats(0);  // Disable statistics box
    trkChargeHist->GetXaxis()->SetTitle("Track charge");
    trkChargeHist->GetYaxis()->SetTitle("Number of tracks");
    trkChargeHist->Draw();
    canvas6->Print("trkChargeHist.png");
    

    // Track eta distribution
    TCanvas* canvas7 = new TCanvas();
    trkEtaHist->SetStats(0);  // Disable statistics box
    trkEtaHist->GetXaxis()->SetTitle("Track eta");
    trkEtaHist->GetYaxis()->SetTitle("Number of tracks");
    trkEtaHist->Draw();
    canvas7->Print("trkEtaHist.png");
    
    




    

    //...
}



