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
#include <cmath>
#include <chrono>


int isData = 1;  // 1 for real data, 0 for MC
bool isBigStatistics = false;
bool toFarm = false;

std::string farm_out = (toFarm == true) ? "/farm_out/" : "/";
std::string root_file_path = (isBigStatistics == true) 
? "../data/outH2R/allRunsP1NickPart_2023.dat_QADB_Valerii_runs_first_electron.root"
: "../data/outH2R_test/allRunsP1NickPart_2023_short.dat_QADBtest_first_electron.root";


// Define the suffix manually using substr()
std::string target_suffix = "trigger_electron.root";
std::string suffix = 
    (root_file_path.size() >= target_suffix.size() && 
     root_file_path.substr(root_file_path.size() - target_suffix.size()) == target_suffix) 
    ? "trigger_electron" 
    : "first_electron";

// Define the output folder as a constant
const std::string OUTPUT_FOLDER = "../analysis_out_" + suffix + farm_out ;



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

void rotatedY_vs_rotated_x_sectorwise(ROOT::RDF::RNode rdf) {
    for (int sector = 1; sector <= 6; ++sector) {
        // Create histogram
        TH2D hist2D(Form("rotatedY_vs_rotatedX_sector%d", sector),
                    Form("Rotated Y vs Rotated X Sector %d; Rotated X; Rotated Y", sector),
                    100, -200, 200, 100, -200, 200);

        // Use Foreach to extract data and fill histogram
        rdf.Filter(Form("el_sector == %d", sector)).Foreach(
            [&](const TVector3& vec) {
                hist2D.Fill(vec.X(), vec.Y()); // Fill histogram
            },
            {"dcR1Vector3_rot"} // Extract the vector without defining extra columns
        );

        // Draw the histogram
        TCanvas canvas(Form("c8_sector%d", sector), Form("Rotated Y vs Rotated X Sector %d", sector), 800, 600);
        hist2D.Draw("COLZ");
        canvas.SaveAs((OUTPUT_FOLDER + Form("rotatedY_vs_rotatedX_sector%d.png", sector)).c_str());
        std::cout << "Saved 2D histogram as rotatedY_vs_rotatedX_sector" << sector << ".png" << std::endl;
    }
}

void rotatedY_vs_rotated_x_all_sectors(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c8", "Rotated Y vs Rotated X for All Sectors", 1200, 800);
    canvas.Divide(3, 2); // Create 6 subpads (3 columns × 2 rows)

    std::vector<TH2D*> histograms;

    for (int sector = 1; sector <= 6; ++sector) {
        histograms.push_back(new TH2D(Form("rotatedY_vs_rotatedX_sector%d", sector),
                                      Form("Sector %d; Rotated X; Rotated Y", sector),
                                      100, -200, 200, 100, -200, 200));

        rdf.Filter(Form("el_sector == %d", sector)).Foreach(
            [&](const TVector3& vec) {
                histograms[sector - 1]->Fill(vec.X(), vec.Y());
            },
            {"dcR1Vector3_rot"} // Extract vector directly
        );

        canvas.cd(sector); // Move to subpad for the sector
        histograms[sector - 1]->Draw("COLZ");
    }

    canvas.SaveAs((OUTPUT_FOLDER + "rotatedY_vs_rotatedX_all_sectors.png").c_str());
    std::cout << "Saved 2D histograms for all sectors on one canvas." << std::endl;
}



void unrotatedY_vs_unrotated_x_all_sectors(ROOT::RDF::RNode rdf) {
    TCanvas canvas("c8", "unRotated Y vs unRotated X for All Sectors", 1200, 800);
    canvas.Divide(3, 2); // Create 6 subpads (3 columns × 2 rows)

    std::vector<TH2D*> histograms;

    for (int sector = 1; sector <= 6; ++sector) {
        histograms.push_back(new TH2D(Form("unrotatedY_vs_unrotatedX_sector%d", sector),
                                      Form("Sector %d; unRotated X; unRotated Y", sector),
                                      100, -200, 200, 100, -200, 200));

        rdf.Filter(Form("el_sector == %d", sector)).Foreach(
            [&](const TVector3& vec) {
                histograms[sector - 1]->Fill(vec.X(), vec.Y());
            },
            {"dcR1Vector3"} // Extract vector directly
        );

        canvas.cd(sector); // Move to subpad for the sector
        histograms[sector - 1]->Draw("COLZ");
    }

    canvas.SaveAs((OUTPUT_FOLDER + "unrotatedY_vs_unrotatedX_all_sectors.png").c_str());
    std::cout << "Saved 2D histograms for all sectors on one canvas." << std::endl;
}




void unrotatedY_vs_unrotated_x_sectorwise(ROOT::RDF::RNode rdf) {
    for (int sector = 1; sector <= 6; ++sector) {
        // Create histogram
        TH2D hist2D(Form("unrotatedY_vs_unrotatedX_sector%d", sector),
                    Form("unRotated Y vs unRotated X Sector %d; unRotated X; unRotated Y", sector),
                    100, -200, 200, 100, -200, 200);

        // Use Foreach to extract data and fill histogram
        rdf.Filter(Form("el_sector == %d", sector)).Foreach(
            [&](const TVector3& vec) {
                hist2D.Fill(vec.X(), vec.Y()); // Fill histogram
            },
            {"dcR1Vector3"} // Extract the vector without defining extra columns
        );

        // Draw the histogram
        TCanvas canvas(Form("c8_sector%d", sector), Form("unRotated Y vs unRotated X Sector %d", sector), 800, 600);
        hist2D.Draw("COLZ");
        canvas.SaveAs((OUTPUT_FOLDER + Form("unrotatedY_vs_unrotatedX_sector%d.png", sector)).c_str());
        std::cout << "Saved 2D histogram as unrotatedY_vs_unrotatedX_sector" << sector << ".png" << std::endl;
    }
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
        line->SetLineStyle(2); // Dashed line
        line->SetLineWidth(2);
        line->Draw("SAME");

        // Save the canvas as an image
        std::string filename = OUTPUT_FOLDER + "W_sector" + std::to_string(sector) + ".png";
        c1.SaveAs(filename.c_str());
        std::cout << "Saved 1D histogram as " << filename << std::endl;
    }
}



void W_for_each_Q2_bin(ROOT::RDF::RNode rdf, const std::string& root_filename) {
    const float wBinSize = 0.05;
    const unsigned nWBins = 30;
    const float lowBorderW = 1.025;
    const float highBorderW = lowBorderW + wBinSize * nWBins;

    // Open a ROOT file to store histograms
    TFile rootFile((OUTPUT_FOLDER + root_filename).c_str(), "RECREATE");

    for (int sector = 1; sector <= 6; ++sector) {
        TCanvas canvas("c8", Form("W_for_each_Q2_bin_sector %d", sector), 1200, 800);
        canvas.Divide(5, 2); // Create 10 subpads (5 columns × 2 rows)

        std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms; // Store histograms

        for (int q2bin = 6; q2bin <= 14; ++q2bin) {
            // Define filtered dataframe
            auto filtered_df = rdf.Filter(Form("el_sector == %d && Q2_bin == %d", sector, q2bin));

            // Create histogram using Histo1D
            auto hist1D = filtered_df.Histo1D(
                {Form("W_for_each_Q2_bin_%d_sector_%d", q2bin, sector),
                 Form("W for Q2 bin %d, Sector %d; W (GeV); Events", q2bin, sector),
                 nWBins, lowBorderW, highBorderW},
                "W_corr"
            );

            histograms.push_back(hist1D); // Store histogram

            // Force execution
            hist1D->GetEntries();  // Forces computation

            // Draw histogram on correct subpad
            canvas.cd(q2bin - 4);  // Adjust indexing so q2bin=5 maps to pad 1
            hist1D->Draw();

            // Save histogram to ROOT file
            hist1D->Write();  
        }

        // Save the canvas as PNG
        std::string png_filename = OUTPUT_FOLDER + Form("W_for_each_Q2_bin_sector_%d.png", sector);
        canvas.SaveAs(png_filename.c_str());
        std::cout << "Saved histogram for sector " << sector << " as " << png_filename << std::endl;
    }

    // Close the ROOT file
    rootFile.Close();
    std::cout << "Histograms saved in ROOT file: " << OUTPUT_FOLDER + root_filename<<std::endl;

}

void plotWvsQ2andSector_SaveROOT(ROOT::RDF::RNode rdf, const std::string& outputFileName) {
    const float wBinSize = 0.05;
    const unsigned nWBins = 30;
    const float lowBorderW = 1.025;
    const float highBorderW = lowBorderW + wBinSize * nWBins;

    // Create a 3D histogram
    auto h3 = rdf.Histo3D(
        {"h3", "W vs Q2_bin vs Sector;W_{corr};Q2_bin;Sector",
         nWBins, lowBorderW, highBorderW, // W axis
         9, 6, 15,                        // Q2_bin axis (bins for 6–14)
         6, 1, 7},                        // Sector axis (bins for 1–6)
        "W_corr", "Q2_bin", "el_sector");
    TH3* h3_ptr = h3.GetPtr();

    // Open output ROOT file for saving histograms
    TFile outFile((OUTPUT_FOLDER + outputFileName).c_str(), "RECREATE");
    if (outFile.IsZombie()) {
        std::cerr << "Error: Cannot open output file " << outputFileName << std::endl;
        return;
    }

    // Loop over sectors (1–6)
    for (int sector = 1; sector <= 6; ++sector) {
        TCanvas* c = new TCanvas(Form("c_sector%d", sector), Form("Sector %d", sector), 1200, 800);
        c->Divide(3, 3);

        // Loop over Q2_bin values (6–14)
        for (int i = 0; i < 9; ++i) {
            int q2_bin = 6 + i;
            c->cd(i + 1);

            // Project W distribution for given sector & Q2_bin
            TH1D* h_proj = (TH1D*)h3_ptr->ProjectionX(
                Form("h_W_sec%d_Q2bin%d", sector, q2_bin),
                q2_bin - 5, q2_bin - 5,   // Correct Q2_bin indexing
                sector, sector);         // Correct sector indexing

            // Histogram appearance
            h_proj->SetTitle(Form("Sector %d, Q2_bin %d;W_{corr};Counts", sector, q2_bin));
            h_proj->Draw("HIST");

            // === SAVE HISTOGRAM INTO ROOT FILE ===
            outFile.cd();
            h_proj->Write();  // Writes each histogram with unique name
        }

        c->Update();
        // === SAVE CANVAS AS PNG FILE ===
        std::string png_filename = OUTPUT_FOLDER + Form("NEW_W_for_each_Q2_bin_sector_%d.png", sector);
        c->SaveAs(png_filename.c_str());
        delete c;  // Clean up canvas after saving
    }

    // Close the ROOT file after writing all histograms
    outFile.Close();
    std::cout << "All histograms saved in: " << outputFileName << std::endl;
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
    auto start = std::chrono::high_resolution_clock::now(); // STRAT
    // Ensure triangleCutParams is initialized
    initializeTriangleCut(isData);

    // Load ROOT file and convert TTrees to RDataFrame
    ROOT::EnableImplicitMT(); // Enable multi-threading
    auto rdf = convert_ttrees_to_rdataframe(root_file_path);
    if (rdf.GetColumnNames().empty()) {
        std::cerr << "Error: Could not create RDataFrame." << std::endl;
        return 1;
    }

    // Define necessary variables in RDataFrame

    auto init_rdf = rdf.Filter("p4_ele_px.size() >0")
                        .Define("el_initial", "return TLorentzVector(0, 0, 10.6, 10.6);")
                        .Define("el_final", "return TLorentzVector(p4_ele_px[0], p4_ele_py[0], p4_ele_pz[0], p4_ele_E[0]);")
                        .Define("proton_initial", "return TLorentzVector(0, 0, 0, 0.938);")
                        .Define("el_px", "return p4_ele_px[0];")
                        .Define("el_py", "return p4_ele_py[0];")
                        .Define("el_pz", "return p4_ele_pz[0];")
                        .Define("el_sector", "return int(sectorE[0]);")
                        //.Define("el_final_corr", +Get4mom_corr, {"el_px", "el_py", "el_pz", "el_sector"})
                        .Define("el_final_corr", [=](double el_px, double el_py, double el_pz, int el_sector) {
                            return Get4mom_corr(el_px, el_py, el_pz, el_sector, isData);
                        }, {"el_px", "el_py", "el_pz", "el_sector"})
                        .Define("el_phi", [](const TLorentzVector& el_final, int el_sector) {
                            return calculate_phi_theta(el_final, el_sector).first;}, {"el_final", "el_sector"})
                        .Define("el_theta", [](const TLorentzVector& el_final, int el_sector) {
                            return calculate_phi_theta(el_final, el_sector).second;}, {"el_final", "el_sector"})
                        .Define("phiSpikeCut", [](double el_phi, double el_theta, int el_sector) {
                            return phiSpikeCut(el_phi, el_theta, el_sector, 1); }, {"el_phi", "el_theta", "el_sector"})
                        .Filter("phiSpikeCut == true")
                        .Define("el_vz", "return p4_ele_vz[0];")
                        .Define("el_vz_cut", [](double el_vz) { return CutVz(el_vz, 1); }, {"el_vz"})
                        .Filter("el_vz_cut == true")
                        .Define("Hx_pcal", "return pcalHX[0];")
                        .Define("Hy_pcal", "return pcalHY[0];")
                        .Define("Hx_ecin", "return ecinHX[0];")
                        .Define("Hy_ecin", "return ecinHY[0];")
                        .Define("Hx_ecout", "return ecoutHX[0];")
                        .Define("Hy_ecout","return ecoutHY[0];")
                        .Define("BadElementKnockOut", [](double Hx_pcal, double Hy_pcal, double Hx_ecin, double Hy_ecin, double Hx_ecout, double Hy_ecout, int el_sector) {
                            return BadElementKnockOut(Hx_pcal, Hy_pcal, Hx_ecin, Hy_ecin, Hx_ecout, Hy_ecout, el_sector, 1); }, {"Hx_pcal", "Hy_pcal", "Hx_ecin", "Hy_ecin", "Hx_ecout", "Hy_ecout", "el_sector"})
                        .Filter("BadElementKnockOut == true")
                        .Define("ecin_Energy", "return ecinE[0];")
                        .Define("pcal_Energy", "return pcalE[0];")
                        .Define("ecout_Energy", "return ecoutE[0];")
                        .Define("ecin_sf", "return ecin_Energy / el_final_corr.P();")
                        .Define("pcal_sf", "return pcal_Energy / el_final_corr.P();")
                        .Define("SFTriangleCut", [](double ecinE, double pcalE, const TLorentzVector& el_final_corr, int sector) {
                            float shift = 0.0f;
                            int pBin = static_cast<int>(el_final_corr.P()/1);
                            return SFTriangleCut(ecinE, pcalE, triangleCutParams, sector, pBin, shift);}, {"ecin_sf", "pcal_sf", "el_final_corr", "el_sector"})
                        .Filter("SFTriangleCut == true")
                        .Define("dcR1Vector3", "return TVector3(dcXR1[0], dcYR1[0], dcZR1[0]);")
                        .Define("dcR2Vector3", "return TVector3(dcXR2[0], dcYR2[0], dcZR2[0]);")
                        .Define("dcR3Vector3", "return TVector3(dcXR3[0], dcYR3[0], dcZR3[0]);")
                        .Define("dcR1Vector3_rot", [](TVector3 vec, int sec) {
                            sec = sec-1;
                            vec.RotateZ(-60 * sec / 57.2958);
                            vec.RotateY(-25 / 57.2958);
                            return vec;}, {"dcR1Vector3", "el_sector"})
                        .Define("dcR2Vector3_rot", [](TVector3 vec, int sec) {
                            sec = sec-1;
                            vec.RotateZ(-60 * sec / 57.2958);
                            vec.RotateY(-25 / 57.2958);
                            return vec;}, {"dcR2Vector3", "el_sector"})
                        .Define("dcR3Vector3_rot", [](TVector3 vec, int sec) {
                            sec = sec-1;
                            vec.RotateZ(-60 * sec / 57.2958);
                            vec.RotateY(-25 / 57.2958);
                            return vec;}, {"dcR3Vector3", "el_sector"})
                        .Define("dcXY",
                            [](const TVector3& dcR1, const TVector3& dcR2, const TVector3& dcR3) {
                                return DCXY{dcR1.X(), dcR1.Y(), dcR2.X(), dcR2.Y(), dcR3.X(), dcR3.Y()};}, {"dcR1Vector3_rot", "dcR2Vector3_rot", "dcR3Vector3_rot"})
                        .Define("CutDCfid",
                            [](const DCXY& dcXY, int sec) {
                                return CutDCfid(dcXY, sec, 1);}, {"dcXY", "el_sector"})
                        .Filter("CutDCfid == true")
                        .Define("pcalLu_f", "return pcalLu[0];")
                        .Define("pcalLv_f", "return pcalLv[0];")
                        .Define("pcalLw_f", "return pcalLw[0];")
                        .Define ("PCALFid_VW", [](double pcalLu, double pcalLv, double pcalLw) {
                            return PCALFid_VW(pcalLv, pcalLw, pcalLu, 1);  }, {"pcalLu_f", "pcalLv_f", "pcalLw_f"})
                        .Filter("PCALFid_VW == true")
                        .Define("Edep", "return ecout_Energy + ecin_Energy + pcal_Energy;")
                        .Define("sf", "return Edep / el_final_corr.P();")
                        .Define("SfCut", [](double sf, double Edep, int sec) {
                            sec = sec-1;
                            return SfCutValerii_Edepos(sf, Edep, sec, 1, isData); }, {"sf", "Edep", "el_sector"})
                        .Filter("SfCut == true")
                        .Define("Q2", "-(el_initial - el_final).M2()")
                        .Define("W", "(el_initial - el_final + proton_initial).M()")
                        .Define("Q2_corr", "-(el_initial - el_final_corr).M2()")
                        .Define("W_corr", "(el_initial - el_final_corr + proton_initial).M()")
                        .Define("W_bin", [](double W_corr) {
                            double low_bin = 1.025;
                            double high_bin = 2.525;
                            double bin_width = 0.05;
                            if (W_corr < low_bin || W_corr >= high_bin) return -1; // Out of range
                                return static_cast<int>((W_corr - low_bin) / bin_width);}, {"W_corr"})
                        .Define("Q2_bin", [](double Q2_corr) {
                            double low_bin = 1;
                            double high_bin = 2500;
                            double bin_number = 50;
                            double delta_Q2 = std::log(high_bin/low_bin)/bin_number;
                            int q2bin = std::log(Q2_corr/low_bin)/delta_Q2;
                                return q2bin;}, {"Q2_corr"})
                        .Define("Q2W_bin","30*Q2_bin + W_bin");


    
    plotWvsQ2andSector_SaveROOT(init_rdf,"qwerty.root");
    //W_for_each_Q2_bin(init_rdf,"JOPA.root");
    //rotatedY_vs_rotated_x_sectorwise(init_rdf);

    //W_for_each_Q2_bin(init_rdf, "W_for_each_Q2_bin.root");  
    //plot_elastic_W_sector(init_rdf);      
    //plot_1d_abs_mom(init_rdf);  
    //unrotatedY_vs_unrotated_x_all_sectors(init_rdf);   
    //rotatedY_vs_rotated_x_all_sectors(init_rdf);  
    //plot_1d_W(init_rdf);
    //plot_1d_Q2corr(init_rdf);  
    //plot_2d_W_vs_QSquared(init_rdf);
    
    //init_rdf.Display({"Q2_bin","W_bin","Q2W_bin"},100)->Print();
    //rotatedY_vs_rotated_x_all_sectors(init_rdf);
    // Print column names
    //std::cout << "Columns in RDataFrame:" << std::endl;
    //for (const auto &col : init_rdf.GetColumnNames()) {
    //    std::cout << col << std::endl;
    //}
     
    //plot_1d_W(init_rdf);
    //plot_2d_W_vs_QSquared(init_rdf);


    auto end = std::chrono::high_resolution_clock::now(); // END

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time of execution: " << elapsed.count() << " sec" << std::endl;

    return 0;
}
