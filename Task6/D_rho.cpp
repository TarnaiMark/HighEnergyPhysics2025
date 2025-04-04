#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>

const int npMax = 84000;

int main(int argc, char* argv[]) {
    // File name to process in the argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    std::string filename = "/mnt/st04pool/project/femtoscopy/data/EPOS3_AuAu_200GeV/" + std::string(argv[1]);
    
    // Open ROOT file and obtain tree
    TFile* inputfile = new TFile(filename.c_str());
    TTree* tree = (TTree*)inputfile->Get("teposevent");
    if(!tree){ std::cerr << "Could not open root file." << std::endl; return 1;}
    
    // Variables to read from the tree
    int np;
    int ist[npMax], id[npMax];
    float t[npMax], x[npMax], y[npMax], z[npMax], mass[npMax], px[npMax], py[npMax], pz[npMax];

    tree->SetBranchAddress("np", &np); // number of particles in the event
    tree->SetBranchAddress("ist", &ist);  // Determines if the particle is final (detected) or not
    tree->SetBranchAddress("id", &id); // Particle id
    tree->SetBranchAddress("t", &t);  // |>Freeze out coordinates
    tree->SetBranchAddress("x", &x);  // |
    tree->SetBranchAddress("y", &y);  // |
    tree->SetBranchAddress("z", &z);  // |
    tree->SetBranchAddress("e", &mass); // !!!e is mass instead of energy!!! Later calculating energy from it
    tree->SetBranchAddress("px", &px);  // |>Impulse components 
    tree->SetBranchAddress("py", &py);  // |
    tree->SetBranchAddress("pz", &pz);  // |

    // Creating histogram for D(rho) with log-binning
    const int nbins = 60;
    double rho_min = 0.8;
    double rho_max = 1000.0;

    double binfactor = TMath::Power(rho_max/rho_min,1.0/nbins);

    double rho_bins[nbins+1];
    rho_bins[0] = rho_min;
    for(int ibin = 1; ibin <= nbins; ibin++)
        rho_bins[ibin] = rho_min * TMath::Power(binfactor,ibin);

    TH1D *rho_hist = new TH1D("rho_hist", "D(#rho); #rho (fm); D(#rho)", nbins, rho_bins);
    rho_hist->Sumw2(); // To calculete errors correctly

    // To store selected particles
    std::vector<double> x_vals, y_vals, z_vals, t_vals, px_vals, py_vals, pz_vals, E_vals;
    
    // Looping through events
    int nEntries = tree->GetEntries();
    for (int ievt = 0; ievt < nEntries; ievt++) {
        tree->GetEntry(ievt);
        
        // Apply one particle cuts
        for (int iPart = 0; iPart < np; iPart++)
        {
            
            if (id[iPart] != 120 or ist[iPart] != 0) continue;
            
            double p_abs = TMath::Sqrt(px[iPart]*px[iPart] + py[iPart]*py[iPart] + pz[iPart]*pz[iPart]);
            double energy = TMath::Sqrt(p_abs*p_abs + mass[iPart]*mass[iPart]);
            double pt = TMath::Sqrt(px[iPart]*px[iPart] + py[iPart]*py[iPart]);
            double eta = 0.5 * TMath::Log((p_abs + pz[iPart]) / (p_abs - pz[iPart]));
            
            if (0.15 > pt or pt > 1.0 or eta > 1 ) continue;
            
            x_vals.push_back(x[iPart]);
            y_vals.push_back(y[iPart]);
            z_vals.push_back(z[iPart]);
            t_vals.push_back(t[iPart]);
            px_vals.push_back(px[iPart]);
            py_vals.push_back(py[iPart]);
            pz_vals.push_back(pz[iPart]);
            E_vals.push_back(energy);
        }
        
        // Looping through particle pairs and calculate rho
        for (size_t i = 0; i < x_vals.size(); i++) {
            for (size_t j = i + 1; j < x_vals.size(); j++) {
                
                double K_x = (px_vals[i] + px_vals[j]) / 2;
                double K_y = (py_vals[i] + py_vals[j]) / 2;
                double K_z = (pz_vals[i] + pz_vals[j]) / 2;
                double K_0 = (E_vals[i] + E_vals[j]) / 2;
                double k_T = TMath::Sqrt(K_x*K_x + K_y*K_y);

                double r_x = x_vals[i] - x_vals[j];
                double r_y = y_vals[i] - y_vals[j];
                double r_z = z_vals[i] - z_vals[j];
                double t = t_vals[i] - t_vals[j];

                /* Calculate rhoOut, rhoSide, rhoLong */
                double rhoOut =  (r_x * K_x / k_T) + (r_y * K_y / k_T) - k_T * (K_0 * t - K_z *r_z) / (K_0 * K_0  - K_z * K_z);
                double rhoSide = -r_x * K_y / k_T + r_y * K_x / k_T;
                double rhoLong = (K_0 * r_z - K_z * t) / TMath::Sqrt(K_0 * K_0 - K_z * K_z);
                
                double rho = TMath::Sqrt(rhoOut*rhoOut + rhoSide*rhoSide + rhoLong*rhoLong);
                
                // To get a normalized histogram dividing with 4*pi*rho^2 and (number of particle pairs) is needed
                double weight = 1.0 / (4 * TMath::Pi() * rho * rho) / (nEntries*x_vals.size()*(x_vals.size()-1)/2);
                rho_hist->Fill(rho, weight);
            }
        }
        
        // Clear the vectors storing the selected particles
        x_vals.clear();
        y_vals.clear();
        z_vals.clear();
        t_vals.clear();
        px_vals.clear();
        py_vals.clear();
        pz_vals.clear();
        E_vals.clear();
    }
    
    // Normalize the histogram by bin widths
    for (int i = 1; i <= rho_hist->GetNbinsX(); i++) {
        double bin_content = rho_hist->GetBinContent(i);
        double bin_width = rho_hist->GetBinWidth(i);
        rho_hist->SetBinContent(i, bin_content / bin_width);
        rho_hist->SetBinError(i, rho_hist->GetBinError(i) / bin_width); // Error normalization
    }
    
    // Saving the histogram to a root file
    std::string input_filename = argv[1];
    /* Need to create histos directory for this to work */
    std::string output_filename = "histos/" + input_filename.substr(0, input_filename.find(".root")) + "-DrhoHist.root";

    TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
    rho_hist->Write();
    outfile->Close();

    std::cout << "Histogram saved: " << output_filename << std::endl;
    
    /* Check the result histogram by plotting it */
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    rho_hist->Draw();
    c1->Print("DrhoHist.png");

    
    return 0;
}