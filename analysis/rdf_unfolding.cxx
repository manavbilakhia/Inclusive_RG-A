// rdf_unfolding.cxx
// Compile with: g++ rdf_unfolding.cxx -o rdf_unfolding `root-config --cflags --libs`


/*
g++ rdf_unfolding.cxx -o rdf_unfolding `root-config --cflags --libs` -I/w/hallb-scshelf2102/clas12/manavb/RooUnfold/src -L/w/hallb-scshelf2102/clas12/manavb/RooUnfold/build -lRooUnfold
*/

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "cuts.cxx"
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <iostream>

using namespace std;

int main() {
    ROOT::EnableImplicitMT();

    // Load ROOT file
    string filePath = "/w/hallb-scshelf2102/clas12/manavb/grad/Inclusive_RG-A/analysis_out_MC_first_electron/rooUnfold.root";  // Replace with your input
    TFile *file = TFile::Open(filePath.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Failed to open file: " << filePath << endl;
        return 1;
    }

    ROOT::RDataFrame df("T", filePath);  // replace "T" with your actual TTree name

    // Create response matrix
    const int nbins = 20;
    double min = 0.0, max = 5.0;
    RooUnfoldResponse response(nbins, min, max);

    df.Foreach([&response](double reco, double truth, bool passCut) {
        if (passCut) {
            response.Fill(reco, truth);
        } else {
            response.Miss(truth);
        }
    }, {"W_corr", "W_gen", "passesCut"});

    // Fill reco histogram
    auto hReco = df.Filter("passesCut").Histo1D({"hReco", "Reco W", nbins, min, max}, "W_corr");

    // Unfold
    RooUnfoldBayes unfold(&response, hReco.GetPtr(), 4); // 4 iterations
    TH1D* hUnfold = (TH1D*) unfold.Hunfold(RooUnfolding::kErrors);

    // Save result
    TFile out("unfolded_output.root", "RECREATE");
    hUnfold->Write("hUnfolded");
    hReco->Write("hReco");
    out.Close();

    cout << "Unfolding complete!" << endl;
    return 0;
}