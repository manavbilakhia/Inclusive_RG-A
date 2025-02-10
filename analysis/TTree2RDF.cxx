#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>

const std::string OUTPUT_FOLDER = "../analysis_out/";

ROOT::RDataFrame convert_ttrees_to_rdataframe(const std::string &root_file_path) {
    TFile *file = TFile::Open(root_file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file " << root_file_path << std::endl;
        return ROOT::RDataFrame(0);
    }

    std::vector<std::string> keys;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (std::string(key->GetClassName()) == "TTree") {
            keys.push_back(key->GetName());
        }
    }

    if (keys.empty()) {
        std::cerr << "No TTrees found in the ROOT file." << std::endl;
        return ROOT::RDataFrame(0);
    }

    std::string tree_name = keys[0];
    std::cout << "Processing TTree: " << tree_name << std::endl;

    ROOT::RDataFrame rdf(tree_name, root_file_path);
    file->Close();
    return rdf;
}

// Use ROOT::RDF::RNode instead of RDataFrame& to fix type mismatch
void plot_1d_W(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c1", "W Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("W_distribution", "Invariant Mass W; W (GeV); Events", 100, 0, 5), "W");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "W_distribution.png").c_str());
    std::cout << "Saved 1D histogram as W_distribution.png" << std::endl;
}

void plot_1d_QSquared(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c3", "Q^{2} Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("QSquared", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "Q2_distribution.png").c_str());
    std::cout << "Saved 1D histogram as QSquared_distribution.png" << std::endl;
}

void plot_2d_W_vs_QSquared(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c4", "W vs QSquared", 800, 600);
    auto hist2D = rdf.Histo2D(ROOT::RDF::TH2DModel("W_vs_Q2", "W vs Q^{2}; W (GeV); Q^2 (GeV^2)", 100, 0, 3, 100, 0, 10), "W", "Q2");
    hist2D->Draw("COLZ");
    canvas.SaveAs((OUTPUT_FOLDER + "W_vs_Q2.png").c_str());
    std::cout << "Saved 2D histogram as W_vs_Q2.png" << std::endl;
}

int TTree2RDF() {
    std::string root_file_path = "../data/outH2R_test/allRunsP1NickPart_2023.dat_QADBtest.root";
    auto rdf = convert_ttrees_to_rdataframe(root_file_path);
    if (rdf.GetColumnNames().empty()) {
        std::cerr << "Error: Could not create RDataFrame." << std::endl;
        return 1;
    }

    // Define necessary variables in RDataFrame
    auto rdf2 = rdf.Define("el_final", "return TLorentzVector(p4_ele_px[0], p4_ele_py[0], p4_ele_pz[0], p4_ele_E[0]);")
                   .Define("el_initial", "return TLorentzVector(0, 0, 10.6, 10.6);")
                   .Define("proton_initial", "return TLorentzVector(0, 0, 0, 0.938);")
                   .Define("Q2", "-(el_initial - el_final).M2()")
                   .Define("W", "(el_initial - el_final + proton_initial).M()");

    // Print column names
    std::cout << "Columns in RDataFrame:" << std::endl;
    for (const auto &col : rdf2.GetColumnNames()) {
        std::cout << col << std::endl;
    }

    // Generate plots
    plot_1d_W(rdf2);
    plot_1d_QSquared(rdf2);
    plot_2d_W_vs_QSquared(rdf2);

    return 0;
}
