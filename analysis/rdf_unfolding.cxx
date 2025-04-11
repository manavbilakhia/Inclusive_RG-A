// rdf_unfolding.cxx
// Compile with: g++ rdf_unfolding.cxx -o rdf_unfolding_executable `root-config --cflags --libs` -I/w/hallb-scshelf2102/clas12/manavb/RooUnfold/src -L/w/hallb-scshelf2102/clas12/manavb/RooUnfold/build -lRooUnfold

// setenv LD_LIBRARY_PATH /w/hallb-scshelf2102/clas12/manavb/RooUnfold/build:$LD_LIBRARY_PATH

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "cuts.cxx"
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int main() {
    ROOT::EnableImplicitMT();

    string filePath = "/w/hallb-scshelf2102/clas12/manavb/grad/Inclusive_RG-A/analysis_out_MC_first_electron/farm_out/rooUnfold.root";
    TFile* file = TFile::Open(filePath.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Failed to open file: " << filePath << endl;
        return 1;
    }

    ROOT::RDataFrame df("T", filePath);

    // Define bin edges from your table
    std::vector<double> Q2_bin_edges = {2.557, 2.990, 3.497, 4.089, 4.782, 5.592, 6.539, 7.646, 8.942, 10.456};
    std::vector<double> W_bin_edges = {
        1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425,
        1.475, 1.525, 1.575, 1.625, 1.675, 1.725,
        1.775, 1.825, 1.875, 1.925, 1.975, 2.025,
        2.075, 2.125, 2.175, 2.225, 2.275, 2.325,
        2.375, 2.425, 2.475, 2.525
    };

    const int W_bins = W_bin_edges.size() - 1;
    const int Q2_bins = Q2_bin_edges.size() - 1;
    const int total_bins = W_bins * Q2_bins;

    // Define binning function for RDataFrame
    auto get_bin_index = [](double val, const std::vector<double>& edges) {
        auto it = std::upper_bound(edges.begin(), edges.end(), val);
        int bin = std::distance(edges.begin(), it) - 1;
        return (bin >= 0 && bin < edges.size() - 1) ? bin : -1;
    };

    auto df_binned = df
        .Redefine("W_bin", [W_bin_edges, get_bin_index](double W_corr) {
            return get_bin_index(W_corr, W_bin_edges);
        }, {"W_corr"})
        .Redefine("Q2_bin", [Q2_bin_edges, get_bin_index](double Q2_corr) {
            return get_bin_index(Q2_corr, Q2_bin_edges);
        }, {"Q2_corr"})
        .Define("W_gen_bin", [W_bin_edges, get_bin_index](double W_gen) {
            return get_bin_index(W_gen, W_bin_edges);
        }, {"W_gen"})
        .Define("Q2_gen_bin", [Q2_bin_edges, get_bin_index](double Q2_gen) {
            return get_bin_index(Q2_gen, Q2_bin_edges);
        }, {"Q2_gen"})
        .Redefine("Q2W_bin", "Q2_bin * " + std::to_string(W_bins) + " + W_bin")
        .Define("Q2W_gen_bin", "Q2_gen_bin * " + std::to_string(W_bins) + " + W_gen_bin")
        .Filter("Q2W_bin >= 0 && Q2W_bin < " + std::to_string(total_bins) +
                " && Q2W_gen_bin >= 0 && Q2W_gen_bin < " + std::to_string(total_bins));

    // Build response matrix
    RooUnfoldResponse response(total_bins, 0, total_bins);
    df_binned.Foreach([&response](int recoBin, int truthBin, bool passCut) {
        if (passCut) response.Fill(recoBin, truthBin);
        else response.Miss(truthBin);
    }, {"Q2W_bin", "Q2W_gen_bin", "passesCut"});

    // Reco histogram
    auto hReco = df_binned.Filter("passesCut")
    .Histo1D({"hReco", "Reco Q2-W Bin", total_bins, 0.0, static_cast<double>(total_bins)}, "Q2W_bin");

    // Unfold
    RooUnfoldBayes unfold(&response, hReco.GetPtr(), 4);
    TH1D* hUnfold = (TH1D*) unfold.Hunfold(RooUnfolding::kErrors);

    // Response matrix cropped view
    TH2D* hResponse = (TH2D*)response.Hresponse();
    hResponse->GetXaxis()->SetRangeUser(0, 600);
    hResponse->GetYaxis()->SetRangeUser(0, 600);
    TCanvas* cResponse = new TCanvas("cResponse", "Response Matrix", 600, 600);
    hResponse->Draw("COLZ");
    cResponse->SaveAs("response_matrix_2D.png");

    // 1D unfolded W for Q2_bin == 0
    TH1D* hUnfolded1D_W = new TH1D("hUnfolded1D_W", "Unfolded W for Q^{2} bin 0;W (GeV);Events", W_bins, &W_bin_edges[0]);
    for (int wbin = 0; wbin < W_bins; ++wbin) {
        hUnfolded1D_W->SetBinContent(wbin + 1, hUnfold->GetBinContent(wbin + 1));
        hUnfolded1D_W->SetBinError(wbin + 1, hUnfold->GetBinError(wbin + 1));
    }

    TCanvas* cW = new TCanvas("cW", "Unfolded W in Q2 bin 0", 800, 600);
    hUnfolded1D_W->Draw("E1");
    cW->SaveAs("unfolded_W_Q2bin0.png");

    // W_rec vs W_gen for Q2_bin == 0
    auto df_q2bin0 = df_binned.Filter("Q2_bin == 0 && W_bin >= 0 && W_gen_bin >= 0");
    TH2D* hWrecVsWgen = new TH2D("hWrecVsWgen", "W_{rec} vs W_{gen} for Q^{2} bin 0;W_{gen} (GeV);W_{rec} (GeV)",
                                 52, 1.0, 2.6, 52, 1.0, 2.6);
    df_q2bin0.Foreach([&hWrecVsWgen](double wrec, double wgen) {
        hWrecVsWgen->Fill(wgen, wrec);
    }, {"W_corr", "W_gen"});
    TCanvas* cWrecGen = new TCanvas("cWrecGen", "W_{rec} vs W_{gen}", 700, 600);
    hWrecVsWgen->Draw("COLZ");
    cWrecGen->SaveAs("Wrec_vs_Wgen_Q2bin0.png");

    // Reco vs truth bin index (0–300)
    TH2D* hRecoVsTruthBin = new TH2D("hRecoVsTruthBin", "Reco vs Truth Bins;Truth bin;Reco bin",
                                     300, 0, 300, 300, 0, 300);
    df_binned.Foreach([&hRecoVsTruthBin](int reco, int truth) {
        hRecoVsTruthBin->Fill(truth, reco);
    }, {"Q2W_bin", "Q2W_gen_bin"});
    TCanvas* cRecoTruth = new TCanvas("cRecoTruth", "Reco vs Truth Bins", 700, 600);
    hRecoVsTruthBin->Draw("COLZ");
    cRecoTruth->SaveAs("RecoVsTruth_Q2Wbin_AllQ2.png");

    // Save everything
    TFile out("unfolded_output_2D.root", "RECREATE");
    hUnfold->Write("hUnfolded");
    hReco->Write("hReco");
    hResponse->Write("hResponse");
    hUnfolded1D_W->Write("hUnfolded1D_W_Q2bin0");
    hWrecVsWgen->Write("hWrec_vs_Wgen_Q2bin0");
    hRecoVsTruthBin->Write("hRecoVsTruthBin_AllQ2");
    out.Close();

    cout << "Unfolding and plotting complete!" << endl;
    return 0;
}