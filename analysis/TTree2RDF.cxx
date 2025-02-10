// to run, use  g++ TTree2RDF.cxx -o executable `root-config --cflags --glibs`
// ./executable

#include "cuts.cxx"
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
#include <fstream>


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

void plot_1d_abs_mom(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c1", "Absolute Momenta", 800, 600);
    auto hist1 = rdf.Histo1D(ROOT::RDF::TH1DModel("Absolute Momenta", "Absolute Momenta; |P| (GeV); Events", 200, 0, 10), "el_abs_mom");
    auto hist2 = rdf.Histo1D(ROOT::RDF::TH1DModel("Absolute Momenta corr", "Absolute Momenta corr; |P| (GeV); Events", 200, 0, 10), "el_abs_mom_corr");
    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);

    hist1->Draw("SAME");
    hist2->Draw("SAME");
    canvas.BuildLegend();
    canvas.SaveAs((OUTPUT_FOLDER + "el_abs.png").c_str());
    std::cout << "Saved 1D histogram as el_abs.png" << std::endl;
}
void plot_1d_W(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c1", "W Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("W_distribution", "Invariant Mass W; W (GeV); Events", 100, 0, 5), "W");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "W_distribution.png").c_str());

    TCanvas canvas2("c2", "W Distribution corr", 800, 600);
    auto hist2 = rdf.Histo1D(ROOT::RDF::TH1DModel("W_distributioncorr", "Invariant Mass W; W (GeV); Events", 100, 0, 5), "W_corr");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "W_corr_distribution.png").c_str());
    std::cout << "Saved 1D histogram as Wcorr_distribution.png" << std::endl;
}

void plot_1d_QSquared(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c3", "Q^{2} Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("QSquared", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "Q2_distribution.png").c_str());

    TCanvas canvas2("c4", "Q^{2}corr Distribution", 800, 600);
    auto hist2 = rdf.Histo1D(ROOT::RDF::TH1DModel("Q2_corr", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2_corr");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "Q2corr_distribution.png").c_str());
    std::cout << "Saved 1D histogram as Q2corr_distribution.png" << std::endl;
}

void plot_2d_W_vs_QSquared(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c5", "W vs QSquared", 800, 600);
    auto hist2D = rdf.Histo2D(ROOT::RDF::TH2DModel("W_vs_Q2", "W vs Q^{2}; W (GeV); Q^2 (GeV^2)", 100, 0, 3, 100, 0, 10), "W", "Q2");
    hist2D->Draw("COLZ");
    canvas.SaveAs((OUTPUT_FOLDER + "W_vs_Q2.png").c_str());
    std::cout << "Saved 2D histogram as W_vs_Q2.png" << std::endl;
}

int main() {
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
                   .Define("W", "(el_initial - el_final + proton_initial).M()")
                   .Define("el_px", "return p4_ele_px[0];")
                   .Define("el_py", "return p4_ele_py[0];")
                   .Define("el_pz", "return p4_ele_pz[0];")
//                   .Define("el_phi", "return atan2(p4_ele_py[0], p4_ele_px[0])*57.2958;")
//                   .Define("el_theta", "return p4_ele_theta[0];")
                   .Define("el_sector", "return double(sectorE[0]);")
                   .Define("el_final_corr", +Get4mom_corr, {"el_px", "el_py", "el_pz", "el_sector"})
                   .Define("Q2_corr", "-(el_initial - el_final_corr).M2()")
                   .Define("W_corr", "(el_initial - el_final_corr + proton_initial).M()")
                   .Define("el_abs_mom", "return sqrt(p4_ele_px[0]*p4_ele_px[0] + p4_ele_py[0]*p4_ele_py[0] + p4_ele_pz[0]*p4_ele_pz[0]);")
                   .Define("el_abs_mom_corr", "return sqrt(el_final_corr.Px()*el_final_corr.Px() + el_final_corr.Py()*el_final_corr.Py() + el_final_corr.Pz()*el_final_corr.Pz());")
//                   .Define("phiSpikeCut", +phiSpikeCut, {"el_phi", "el_theta", "el_sector", "1"})
                   .Define("el_vz", "return p4_ele_vz[0];")
                   .Define("el_vz_corr", [](double el_vz) { return CutVz(el_vz, 1); }, {"el_vz"});
                   
//    rdf2.Display({"phiSpikeCut"}, 10)->Print();
    rdf2.Display({"el_vz_corr"}, 10)->Print();

    // Print column names
    std::cout << "Columns in RDataFrame:" << std::endl;
    for (const auto &col : rdf2.GetColumnNames()) {
        std::cout << col << std::endl;
    }

    // Generate plots
    plot_1d_abs_mom(rdf2);
    plot_1d_W(rdf2);
    plot_1d_QSquared(rdf2);
    plot_2d_W_vs_QSquared(rdf2);

    return 0;
}
