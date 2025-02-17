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
#include <TLine.h> 
#include <TLegend.h> 
#include <string>

// Define the root file path
std::string root_file_path = "../data/outH2R_test/allRunsP1NickPart_2023_short.dat_QADBtest_trigger_electron.root";

// Define the suffix manually using substr()
std::string target_suffix = "trigger_electron.root";
std::string suffix = 
    (root_file_path.size() >= target_suffix.size() && 
     root_file_path.substr(root_file_path.size() - target_suffix.size()) == target_suffix) 
    ? "trigger_electron" 
    : "first_electron";

// Define the output folder as a constant
const std::string OUTPUT_FOLDER = "../analysis_out_" + suffix + "/";



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
}

void plot_1d_el_final_corr_P_with_filter(ROOT::RDF::RNode rdf) {
    // Apply a filter for values greater than 10.5
    auto filtered_rdf = rdf.Filter("el_final_corr_P > 10.5");

    TCanvas canvas("c3", "el_final_corr_P Distribution", 800, 600);
    
    // Plot histogram after filtering
    auto hist = filtered_rdf.Histo1D(
        ROOT::RDF::TH1DModel("el_final_corr_P", "el_final_corr_P; el_final_corr_P (GeV^); Events", 100, 10.5, 250), 
        "el_final_corr_P"
    );
    
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "el_final_corr_P_distribution.png").c_str());
}

void plot_1d_el_final_corr_P(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c3", "el_final_corr_P Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("el_final_corr_P", "el_final_corr_P; el_final_corr_P (GeV^); Events", 100, 0, 300), "el_final_corr_P");
    hist->Draw();
    canvas.SaveAs((OUTPUT_FOLDER + "el_final_corr_P_distribution.png").c_str());
}

void plot_1d_Q2corr(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c4", "Q^{2}corr Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("Q2_corr", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2_corr");
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

void plot_2d_pcal_sf_vs_ecin_sf(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c6", "pcal_sf vs ecin_sf", 800, 600);
    auto hist2D = rdf.Histo2D(ROOT::RDF::TH2DModel("pcal_sf_vs_ecin_sf", "pcal_sf vs ecin_sf; ecin_sf; pcal_sf", 100, 0, 0.2, 100, 0, 0.25), "ecin_sf", "pcal_sf");
    hist2D->Draw("COLZ");
    canvas.SaveAs((OUTPUT_FOLDER + "pcal_sf_vs_ecin_sf_afterW-Trianglecut.png").c_str());
    std::cout << "Saved 2D histogram as pcal_sf_vs_ecin_sf.png" << std::endl;
}

void plot_2d_pcalsf_vs_ecinsf_afterW_Trianglecut_bin_sector(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c7", "pcal_sf vs ecin_sf after W and Trianglecut", 800, 600);
    auto hist2D = rdf.Filter("el_sector == 2 && pBin == 7")
                .Histo2D(ROOT::RDF::TH2DModel("pcal_sf_vs_ecin_sf_sector1", "pcal_sf vs ecin_sf; ecin_sf; pcal_sf", 100, 0, 0.2, 100, 0, 0.25), "ecin_sf", "pcal_sf");
    hist2D->Draw("COLZ");
    canvas.SaveAs((OUTPUT_FOLDER + "pcal_sf_vs_ecin_sf_sector2_pBin7.png").c_str());
    std::cout << "Saved 2D histogram as pcal_sf_vs_ecin_sf_sector1.png" << std::endl;
}


void plot_elastic_W_sector(ROOT::RDF::RNode rdf) {
    for (int sector = 1; sector <= 6; ++sector) {
        TCanvas c1(Form("c1_sector%d", sector), Form("Elastic W Sector %d", sector), 800, 600);

        auto hist1 = rdf.Filter(Form("el_sector == %d", sector))
                         .Histo1D(ROOT::RDF::TH1DModel(Form("W_sector%d", sector), 
                                                       Form("Elastic W Sector %d; W (GeV); Events", sector), 
                                                       100, 0.8, 1.2), "W");
        auto hist2 = rdf.Filter(Form("el_sector == %d", sector))
                         .Histo1D(ROOT::RDF::TH1DModel(Form("W_sector%d_corr", sector), 
                                                       Form("Elastic W Sector %d; W (GeV); Events", sector), 
                                                       100, 0.8, 1.2), "W_corr");
        
        hist1->SetLineColor(kRed);
        hist2->SetLineColor(kBlue);
        hist1->Draw("SAME");
        hist2->Draw("SAME");

        // Draw vertical line at proton mass (0.938 GeV)
        TLine *line = new TLine(0.938, 0, 0.938, hist1->GetMaximum());
        line->SetLineColor(kBlack);
        line->SetLineStyle(2); // Dashed line
        line->SetLineWidth(2);
        line->Draw("SAME");

        // Save the canvas as an image
        std::string filename = OUTPUT_FOLDER + "W_sector" + std::to_string(sector) + ".png";
        c1.SaveAs(filename.c_str());
        std::cout << "Saved 1D histogram as " << filename << std::endl;
    }
}


std::pair<double, double> calculate_phi_theta(TLorentzVector el_final, int el_sector) {
    float toRD = 57.2958; // 360/2pi
    double thetaEMeasured = el_final.Theta() * toRD;

    // To get Phi:
    double phiDC = el_final.Phi() * toRD;
    if (phiDC < 0 && el_sector > 0) {
        phiDC = phiDC + 360. - el_sector * 60.;
    } else {
        phiDC = phiDC - el_sector * 60.;
    }

    return {phiDC, thetaEMeasured}; // Returning a pair
}

int main() {
    // Ensure triangleCutParams is initialized
    initializeTriangleCut();

    std::cout << "triangleCutParams check:\n";
    for (int s = 0; s < 6; s++) {
        for (int t = 0; t < 10; t++) {
            std::cout << "Sector: " << s << ", Bin: " << t
                      << " -> Zero: " << triangleCutParams[s][t][0]
                      << ", Slope: " << triangleCutParams[s][t][1] << std::endl;
        }
    }

    // Load ROOT file and convert TTrees to RDataFrame
    auto rdf = convert_ttrees_to_rdataframe(root_file_path);
    if (rdf.GetColumnNames().empty()) {
        std::cerr << "Error: Could not create RDataFrame." << std::endl;
        return 1;
    }

    // Define necessary variables in RDataFrame

    auto init_rdf = rdf.Define("el_initial", "return TLorentzVector(0, 0, 10.6, 10.6);")
                        .Define("el_final", "return TLorentzVector(p4_ele_px[0], p4_ele_py[0], p4_ele_pz[0], p4_ele_E[0]);")
                        .Define("proton_initial", "return TLorentzVector(0, 0, 0, 0.938);")
                        .Define("Q2", "-(el_initial - el_final).M2()")
                        .Define("W", "(el_initial - el_final + proton_initial).M()")
                        .Define("el_px", "return p4_ele_px[0];")
                        .Define("el_py", "return p4_ele_py[0];")
                        .Define("el_pz", "return p4_ele_pz[0];")
                        .Define("el_sector", "return int(sectorE[0]);")
                        .Define("el_final_corr", +Get4mom_corr, {"el_px", "el_py", "el_pz", "el_sector"})
                        .Define("Q2_corr", "-(el_initial - el_final_corr).M2()")
                        .Define("W_corr", "(el_initial - el_final_corr + proton_initial).M()");

    auto W_filtered_rdf = init_rdf//.Filter("W > 0.7 && W < 3.0") //W cut
                    .Define("el_phi", [](const TLorentzVector& el_final, int el_sector) {
                        return calculate_phi_theta(el_final, el_sector).first;}, {"el_final", "el_sector"})
                    .Define("el_theta", [](const TLorentzVector& el_final, int el_sector) {
                        return calculate_phi_theta(el_final, el_sector).second;}, {"el_final", "el_sector"})
                    .Define("el_final_corr_P", "return (sqrt(el_final_corr.Px()*el_final_corr.Px() + el_final_corr.Py()*el_final_corr.Py() + el_final_corr.Pz()*el_final_corr.Pz()));")
                    .Define("el_abs_mom", "return sqrt(p4_ele_px[0]*p4_ele_px[0] + p4_ele_py[0]*p4_ele_py[0] + p4_ele_pz[0]*p4_ele_pz[0]);")
                    .Define("el_abs_mom_corr", "return sqrt(el_final_corr.Px()*el_final_corr.Px() + el_final_corr.Py()*el_final_corr.Py() + el_final_corr.Pz()*el_final_corr.Pz());")
                    .Define("phiSpikeCut", [](double el_phi, double el_theta, int el_sector) {
                        return phiSpikeCut(el_phi, el_theta, el_sector, 1); }, {"el_phi", "el_theta", "el_sector"})
                    .Define("el_vz", "return p4_ele_vz[0];")
                    .Define("el_vz_cut", [](double el_vz) { return CutVz(el_vz, 1); }, {"el_vz"})
                    .Define("Hx_pcal", "return pcalHX[0];")
                    .Define("Hy_pcal", "return pcalHY[0];")
                    .Define("Hx_ecin", "return ecinHX[0];")
                    .Define("Hy_ecin", "return ecinHY[0];")
                    .Define("Hx_ecout", "return ecoutHX[0];")
                    .Define("Hy_ecout","return ecoutHY[0];")
                    .Define("BadElementKnockOut", [](double Hx_pcal, double Hy_pcal, double Hx_ecin, double Hy_ecin, double Hx_ecout, double Hy_ecout, int el_sector) {
                        return BadElementKnockOut(Hx_pcal, Hy_pcal, Hx_ecin, Hy_ecin, Hx_ecout, Hy_ecout, el_sector, 1); }, {"Hx_pcal", "Hy_pcal", "Hx_ecin", "Hy_ecin", "Hx_ecout", "Hy_ecout", "el_sector"})
                    .Define("ecin_Energy", "return ecinE[0];")
                    .Define("pcal_Energy", "return pcalE[0];")
                    .Define("ecin_sf", "return ecin_Energy / el_final_corr.P();")
                    .Define("pcal_sf", "return pcal_Energy / el_final_corr.P();")
                    .Define("pBin", "return static_cast<int>(el_final_corr.P()/1);")
        .Define("SFTriangleCut", [](double ecinE, double pcalE, const TLorentzVector& el_final_corr, int sector) {
            float shift = 0.0f;
            int pBin = static_cast<int>(el_final_corr.P()/1);
            if (el_final_corr.P() > 10.5){
            std::cout << "el_final_corr.P()=" << el_final_corr.P() << std::endl;}
            return SFTriangleCut(ecinE, pcalE, triangleCutParams, sector, pBin, shift);
        }, {"ecin_sf", "pcal_sf", "el_final_corr", "el_sector"})
        .Filter("SFTriangleCut == true");

//    W_filtered_rdf.Display({"pBin"}, 10)->Print();
//    W_filtered_rdf.Display({"ecin_sf"}, 10)->Print();
//    W_filtered_rdf.Display({"pcal_sf"}, 10)->Print();
//    W_filtered_rdf.Display({"el_sector"}, 10)->Print();
    //W_filtered_rdf.Filter("el_final_corr_P > 10.5").Display({"el_final_corr_P"}, 10000000)->Print();

//    W_filtered_rdf.Filter("pBin ==4 && el_final_corr_P <4").Display({"el_final_corr_P"}, 100)->Print();
    //std::cout << "Finished filtering" << std::endl;
    //W_filtered_rdf.Filter("SFTriangleCut == true").Display({"SFTriangleCut"}, 100)->Print();

    // Print column names
//    std::cout << "Columns in RDataFrame:" << std::endl;
//    for (const auto &col : W_filtered_rdf.GetColumnNames()) {
//        std::cout << col << std::endl;
//    }

    // Generate plots
    plot_1d_abs_mom(W_filtered_rdf);
    plot_1d_W(W_filtered_rdf);
    plot_1d_QSquared(W_filtered_rdf);
    plot_2d_W_vs_QSquared(W_filtered_rdf);
    plot_elastic_W_sector(W_filtered_rdf);
    plot_1d_el_final_corr_P(W_filtered_rdf);
      plot_2d_pcal_sf_vs_ecin_sf(W_filtered_rdf);
plot_2d_pcalsf_vs_ecinsf_afterW_Trianglecut_bin_sector(W_filtered_rdf);


    return 0;
}
