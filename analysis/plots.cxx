#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <TVector3.h>
#include <ROOT/RDataFrame.hxx>
#include <iostream>

void plot_1d_abs_mom(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    TCanvas canvas("c1", "Absolute Momenta", 800, 600);
    auto hist1 = rdf.Histo1D(ROOT::RDF::TH1DModel("Absolute Momenta", "Absolute Momenta; |P| (GeV); Events", 200, 0, 10), "el_abs_mom");
    auto hist2 = rdf.Histo1D(ROOT::RDF::TH1DModel("Absolute Momenta corr", "Absolute Momenta corr; |P| (GeV); Events", 200, 0, 10), "el_abs_mom_corr");
    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);

    hist1->Draw("SAME");
    hist2->Draw("SAME");
    canvas.BuildLegend();
    canvas.SaveAs((output_folder + "el_abs.png").c_str());
    std::cout << "Saved 1D histogram as el_abs.png" << std::endl;
}
void plot_1d_W(ROOT::RDF::RNode rdf, const std::string& output_folder) {
    TCanvas canvas("c1", "W Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("W_distribution", "Invariant Mass W; W (GeV); Events", 100, 0, 5), "W");
    hist->Draw();
    canvas.SaveAs((output_folder + "W_distribution.png").c_str());

    TCanvas canvas2("c2", "W Distribution corr", 800, 600);
    auto hist2 = rdf.Histo1D(ROOT::RDF::TH1DModel("W_distributioncorr", "Invariant Mass W; W (GeV); Events", 100, 0, 5), "W_corr");
    hist->Draw();
    canvas.SaveAs((output_folder + "W_corr_distribution.png").c_str());
    std::cout << "Saved 1D histogram as Wcorr_distribution.png" << std::endl;
}

void plot_1d_QSquared(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c3", "Q^{2} Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("QSquared", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2");
    hist->Draw();
    canvas.SaveAs((output_folder + "Q2_distribution.png").c_str());
}

void plot_1d_el_final_corr_P_with_filter(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    // Apply a filter for values greater than 10.5
    auto filtered_rdf = rdf.Filter("el_final_corr_P > 10.5");

    TCanvas canvas("c3", "el_final_corr_P Distribution", 800, 600);
    
    // Plot histogram after filtering
    auto hist = filtered_rdf.Histo1D(
        ROOT::RDF::TH1DModel("el_final_corr_P", "el_final_corr_P; el_final_corr_P (GeV^); Events", 100, 10.5, 250), 
        "el_final_corr_P"
    );
    
    hist->Draw();
    canvas.SaveAs((output_folder + "el_final_corr_P_distribution.png").c_str());
}

void plot_1d_el_final_corr_P(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c3", "el_final_corr_P Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("el_final_corr_P", "el_final_corr_P; el_final_corr_P (GeV^); Events", 100, 0, 300), "el_final_corr_P");
    hist->Draw();
    canvas.SaveAs((output_folder + "el_final_corr_P_distribution.png").c_str());
}

void plot_1d_Q2corr(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c4", "Q^{2}corr Distribution", 800, 600);
    auto hist = rdf.Histo1D(ROOT::RDF::TH1DModel("Q2_corr", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2_corr");
    hist->Draw();
    canvas.SaveAs((output_folder + "Q2corr_distribution.png").c_str());
    std::cout << "Saved 1D histogram as Q2corr_distribution.png" << std::endl;
}

void plot_2d_W_vs_QSquared(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c5", "W vs QSquared", 800, 600);
    auto hist2D = rdf.Histo2D(ROOT::RDF::TH2DModel("W_vs_Q2", "W vs Q^{2}; W (GeV); Q^2 (GeV^2)", 100, 0, 3, 100, 0, 10), "W", "Q2");
    hist2D->Draw("COLZ");
    canvas.SaveAs((output_folder + "W_vs_Q2.png").c_str());
    std::cout << "Saved 2D histogram as W_vs_Q2.png" << std::endl;
}

void plot_2d_pcal_sf_vs_ecin_sf(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c6", "pcal_sf vs ecin_sf", 800, 600);
    auto hist2D = rdf.Histo2D(ROOT::RDF::TH2DModel("pcal_sf_vs_ecin_sf", "pcal_sf vs ecin_sf; ecin_sf; pcal_sf", 100, 0, 0.2, 100, 0, 0.25), "ecin_sf", "pcal_sf");
    hist2D->Draw("COLZ");
    canvas.SaveAs((output_folder + "pcal_sf_vs_ecin_sf_afterW-Trianglecut.png").c_str());
    std::cout << "Saved 2D histogram as pcal_sf_vs_ecin_sf.png" << std::endl;
}

void plot_2d_pcalsf_vs_ecinsf_afterW_Trianglecut_bin_sector(ROOT::RDF::RNode rdf,const std::string& output_folder) {
    TCanvas canvas("c7", "pcal_sf vs ecin_sf after W and Trianglecut", 800, 600);
    auto hist2D = rdf.Filter("el_sector == 2 && pBin == 7")
                .Histo2D(ROOT::RDF::TH2DModel("pcal_sf_vs_ecin_sf_sector1", "pcal_sf vs ecin_sf; ecin_sf; pcal_sf", 100, 0, 0.2, 100, 0, 0.25), "ecin_sf", "pcal_sf");
    hist2D->Draw("COLZ");
    canvas.SaveAs((output_folder + "pcal_sf_vs_ecin_sf_sector2_pBin7.png").c_str());
    std::cout << "Saved 2D histogram as pcal_sf_vs_ecin_sf_sector1.png" << std::endl;
}

void rotatedY_vs_rotated_x_sectorwise(ROOT::RDF::RNode rdf,const std::string& output_folder) {
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
        canvas.SaveAs((output_folder + Form("rotatedY_vs_rotatedX_sector%d.png", sector)).c_str());
        std::cout << "Saved 2D histogram as rotatedY_vs_rotatedX_sector" << sector << ".png" << std::endl;
    }
}

void rotatedY_vs_rotated_x_all_sectors(ROOT::RDF::RNode rdf,const std::string& output_folder) {
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

    canvas.SaveAs((output_folder + "rotatedY_vs_rotatedX_all_sectors.png").c_str());
    std::cout << "Saved 2D histograms for all sectors on one canvas." << std::endl;
}



void unrotatedY_vs_unrotated_x_all_sectors(ROOT::RDF::RNode rdf,const std::string& output_folder) {
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

    canvas.SaveAs((output_folder + "unrotatedY_vs_unrotatedX_all_sectors.png").c_str());
    std::cout << "Saved 2D histograms for all sectors on one canvas." << std::endl;
}




void unrotatedY_vs_unrotated_x_sectorwise(ROOT::RDF::RNode rdf,const std::string& output_folder) {
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
        canvas.SaveAs((output_folder + Form("unrotatedY_vs_unrotatedX_sector%d.png", sector)).c_str());
        std::cout << "Saved 2D histogram as unrotatedY_vs_unrotatedX_sector" << sector << ".png" << std::endl;
    }
}


void plot_elastic_W_sector(ROOT::RDF::RNode rdf,const std::string& output_folder) {
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
        std::string filename = output_folder + "W_sector" + std::to_string(sector) + ".png";
        c1.SaveAs(filename.c_str());
        std::cout << "Saved 1D histogram as " << filename << std::endl;
    }
}



void W_for_each_Q2_bin(ROOT::RDF::RNode rdf, const std::string& root_filename,const std::string& output_folder) {
    const float wBinSize = 0.05;
    const unsigned nWBins = 30;
    const float lowBorderW = 1.025;
    const float highBorderW = lowBorderW + wBinSize * nWBins;

    // Open a ROOT file to store histograms
    TFile rootFile((output_folder + root_filename).c_str(), "RECREATE");

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
        std::string png_filename = output_folder + Form("W_for_each_Q2_bin_sector_%d.png", sector);
        canvas.SaveAs(png_filename.c_str());
        std::cout << "Saved histogram for sector " << sector << " as " << png_filename << std::endl;
    }

    // Close the ROOT file
    rootFile.Close();
    std::cout << "Histograms saved in ROOT file: " << output_folder + root_filename<<std::endl;

}

void plotWvsQ2andSector_SaveROOT(ROOT::RDF::RNode rdf, const std::string& outputFileName,const std::string& output_folder) {
    
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
    TFile outFile((output_folder + outputFileName).c_str(), "RECREATE");
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
        std::string png_filename = output_folder + Form("NEW_W_for_each_Q2_bin_sector_%d.png", sector);
        c->SaveAs(png_filename.c_str());
        delete c;  // Clean up canvas after saving
    }

    // Close the ROOT file after writing all histograms
    outFile.Close();
    std::cout << "All histograms saved in: " << outputFileName << std::endl;
}