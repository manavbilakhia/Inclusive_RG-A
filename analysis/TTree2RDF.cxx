// to run, use  g++ TTree2RDF.cxx -o executable `root-config --cflags --glibs`
// ./executable

#include "plots.cxx"
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
#include <TEntryList.h>
#include <string>
#include <cmath>
#include <chrono>


// Global variables
int isData = 0;  // 1 for real data, 0 for MC
bool isBigStatistics = true;
bool toFarm = true;

// Determine the output folder
const std::string farm_out = toFarm ? "/farm_out/" : "/";

// Determine the root file path using ternary operators
const std::string root_file_path = (isData == 1) ? 
    (isBigStatistics ? 
        "../data/outH2R/allRunsP1NickPart_2023.dat_QADB_Valerii_runs_first_electron.root" :
        "../data/outH2R_test/allRunsP1NickPart_2023_short.dat_QADBtest_first_electron.root") 
    : (isBigStatistics ? 
        "/w/hallb-scshelf2102/clas12/valerii/data/outH2R/simElasFix_iter3/SimIter3_10.datAna.root" :
        "/w/hallb-scshelf2102/clas12/valerii/data/outH2R/simElasFix_iter3/SimIter3_10.datAna.root");


std::string target_suffix = "trigger_electron.root";
std::string suffix = (root_file_path.size() >= target_suffix.size() && 
     root_file_path.substr(root_file_path.size() - target_suffix.size()) == target_suffix) 
    ? "trigger_electron" 
    : "first_electron";
// Define the output folder as a constant
//const std::string OUTPUT_FOLDER = "../analysis_out_" + suffix + farm_out ;

const std::string OUTPUT_FOLDER = (isData == 1) ? "../analysis_out_" + suffix + farm_out : "../analysis_out_MC_" + suffix + farm_out;



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
    fastMC();
    initializeTriangleCut(isData);

    // Load ROOT file and convert TTrees to RDataFrame
    ROOT::EnableImplicitMT(); // Enable multi-threading
    auto rdf = convert_ttrees_to_rdataframe(root_file_path);
    if (rdf.GetColumnNames().empty()) {
        std::cerr << "Error: Could not create RDataFrame." << std::endl;
        return 1;
    }

auto rdf_with_entry = rdf.Define("entry", "rdfentry_"); // Define entry column
auto rdf_filtered = (isData == 0 && !isBigStatistics) 
                    ? rdf_with_entry.Filter("entry < 1000000") 
                    : ROOT::RDF::RNode(rdf_with_entry);


if (rdf_filtered.HasColumn("dcYZ1")) {
    rdf_filtered = rdf_filtered.Define("dcZR1", "dcYZ1");
}
else {
    std::cerr << "Error: Column dcYZ1 not found." << std::endl;
}
 
    auto init_rdf = rdf_filtered
                        //.Filter("p4_ele_px.size() >0")
                        .Define("el_initial", "return TLorentzVector(0, 0, 10.6, 10.6);")
                        .Define("el_final", "return TLorentzVector(p4_ele_px[0], p4_ele_py[0], p4_ele_pz[0], p4_ele_E[0]);")
                        .Define("proton_initial", "return TLorentzVector(0, 0, 0, 0.938);")
                        .Define("el_final_gen", "return TLorentzVector(gen_px[0], gen_py[0], gen_pz[0], sqrt(gen_px[0]*gen_px[0] + gen_py[0]*gen_py[0] + gen_pz[0]*gen_pz[0]));")
                        .Define("el_px", "return p4_ele_px[0];")
                        .Define("el_py", "return p4_ele_py[0];")
                        .Define("el_pz", "return p4_ele_pz[0];")
                        .Define("el_sector", "return int(sectorE[0]);")
                        .Define("el_final_corr", [=](double el_px, double el_py, double el_pz, int el_sector) {
                            return Get4mom_corr(el_px, el_py, el_pz, el_sector, isData);
                        }, {"el_px", "el_py", "el_pz", "el_sector"})
                        .Define("el_phi", [](const TLorentzVector& el_final, int el_sector) {
                            return calculate_phi_theta(el_final, el_sector).first;}, {"el_final_corr", "el_sector"})
                        .Define("el_theta", [](const TLorentzVector& el_final, int el_sector) {
                            return calculate_phi_theta(el_final, el_sector).second;}, {"el_final_corr", "el_sector"})
                        .Define("phiSpikeCut", [](double el_phi, double el_theta, int el_sector) {
                            return phiSpikeCut(el_phi, el_theta, el_sector, 1); }, {"el_phi", "el_theta", "el_sector"})
                        //.Filter("phiSpikeCut == true")
                        .Define("el_vz", "return p4_ele_vz[0];")
                        .Define("el_vz_cut", [](double el_vz) { return CutVz(el_vz, 1); }, {"el_vz"})
                        //.Filter("el_vz_cut == true")
                        .Define("Hx_pcal", "return pcalHX[0];")
                        .Define("Hy_pcal", "return pcalHY[0];")
                        .Define("Hx_ecin", "return ecinHX[0];")
                        .Define("Hy_ecin", "return ecinHY[0];")
                        .Define("Hx_ecout", "return ecoutHX[0];")
                        .Define("Hy_ecout","return ecoutHY[0];")
                        .Define("BadElementKnockOut", [](double Hx_pcal, double Hy_pcal, double Hx_ecin, double Hy_ecin, double Hx_ecout, double Hy_ecout, int el_sector) {
                            return BadElementKnockOut(Hx_pcal, Hy_pcal, Hx_ecin, Hy_ecin, Hx_ecout, Hy_ecout, el_sector, 1); }, {"Hx_pcal", "Hy_pcal", "Hx_ecin", "Hy_ecin", "Hx_ecout", "Hy_ecout", "el_sector"})
                        //.Filter("BadElementKnockOut == true")
                        .Define("ecin_Energy", "return ecinE[0];")
                        .Define("pcal_Energy", "return pcalE[0];")
                        .Define("ecout_Energy", "return ecoutE[0];")
                        .Define("ecin_sf", "return ecin_Energy / el_final_corr.P();")
                        .Define("pcal_sf", "return pcal_Energy / el_final_corr.P();")
                        .Define("SFTriangleCut", [](double ecinE, double pcalE, const TLorentzVector& el_final_corr, int sector) {
                            float shift = 0.0f;
                            int pBin = static_cast<int>(el_final_corr.P()/1);
                            return SFTriangleCut(ecinE, pcalE, triangleCutParams, sector, pBin, shift);}, {"ecin_sf", "pcal_sf", "el_final_corr", "el_sector"})
                        //.Filter("SFTriangleCut == true")
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
                        //.Filter("CutDCfid == true")
                        .Define("pcalLu_f", "return pcalLu[0];")
                        .Define("pcalLv_f", "return pcalLv[0];")
                        .Define("pcalLw_f", "return pcalLw[0];")
                        .Define ("PCALFid_VW", [](double pcalLu, double pcalLv, double pcalLw) {
                            return PCALFid_VW(pcalLv, pcalLw, pcalLu, 1);  }, {"pcalLu_f", "pcalLv_f", "pcalLw_f"})
                        //.Filter("PCALFid_VW == true")
                        .Define("Edep", "return ecout_Energy + ecin_Energy + pcal_Energy;")
                        .Define("sf", "return Edep / el_final_corr.P();")
                        .Define("SfCut", [](double sf, double Edep, int sec) {
                            sec = sec-1;
                            return SfCutValerii_Edepos(sf, Edep, sec, 1, isData); }, {"sf", "Edep", "el_sector"})
                        //.Filter("SfCut == true")
                        .Define("Q2", "-(el_initial - el_final).M2()")
                        .Define("W", "(el_initial - el_final + proton_initial).M()")
                        .Define("W_gen", "(el_initial - el_final_gen + proton_initial).M()")
                        .Define("Q2_gen", "-(el_initial - el_final_gen).M2()")
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
                        .Define("Q2W_bin","30*Q2_bin + W_bin")
                        .Define("el_final_P", "el_final.P()")  // Define new column
                        .Define("el_final_corr_P", "el_final_corr.P()")  // Define new column
                         // this is true when every single defined cut is true
                         .Define("passesCut", [](bool phiSpikeCut, bool el_vz_cut, bool BadElementKnockOut, bool SFTriangleCut, bool CutDCfid, bool PCALFid_VW, bool SfCut) {
                            return phiSpikeCut && el_vz_cut && BadElementKnockOut && SFTriangleCut && CutDCfid && PCALFid_VW && SfCut;}, 
                            {"phiSpikeCut", "el_vz_cut", "BadElementKnockOut", "SFTriangleCut", "CutDCfid", "PCALFid_VW", "SfCut"});

    //save the rdf as a root file
    // Only save standard types
    std::vector<std::string> columns_to_save = {
    "el_final_P", "el_final_corr_P",
    "Q2", "Q2_corr", "Q2_gen", "W", "W_corr", "W_gen",
    "Q2_bin", "W_bin", "Q2W_bin","passesCut"
    // add more as needed, just avoid DCXY or any other struct
    };

    init_rdf.Snapshot("T", OUTPUT_FOLDER + "rooUnfold.root", columns_to_save);


    //std::cout<<OUTPUT_FOLDER<<std::endl;
    //plotWvsQ2andSector_SaveROOT(init_rdf_MC,"qwerty_MC.root" ,OUTPUT_FOLDER);
    //W_for_each_Q2_bin(init_rdf,"JOPA.root",OUTPUT_FOLDER);
    //rotatedY_vs_rotated_x_sectorwise(init_rdf, OUTPUT_FOLDER);

    //W_for_each_Q2_bin(init_rdf, "W_for_each_Q2_bin.root", OUTPUT_FOLDER);  
    //plot_elastic_W_sector(init_rdf, OUTPUT_FOLDER);      
    //plot_1d_abs_mom(init_rdf, OUTPUT_FOLDER);  
    //unrotatedY_vs_unrotated_x_all_sectors(init_rdf, OUTPUT_FOLDER);   
    //rotatedY_vs_rotated_x_all_sectors(init_rdf,   OUTPUT_FOLDER);  
    //plot_1d_W(init_rdf, OUTPUT_FOLDER);
    //plot_1d_Q2corr(init_rdf, OUTPUT_FOLDER);  
    //plot_2d_W_vs_QSquared(init_rdf, OUTPUT_FOLDER);
    
    //init_rdf.Display({"Q2_bin","W_bin","Q2W_bin"},100)->Print();
    //rotatedY_vs_rotated_x_all_sectors(init_rdf, OUTPUT_FOLDER);
    // Print column names
    //rdf_filtered.Display( "dcZR1", 100)->Print();
    //std::cout << "Columns in RDataFrame:" << std::endl;
    //for (const auto &col : rdf.GetColumnNames()) {
    //    std::cout << col << std::endl;
    //}
     
     // show me the first 10 events of genpx, genpy, genpz and W_gen
    //init_rdf.Display({"gen_px", "gen_py", "gen_pz", "W_gen"}, 10)->Print();

    //plot_1d_W(init_rdf, OUTPUT_FOLDER);
    //plot_2d_W_vs_QSquared(init_rdf, OUTPUT_FOLDER);

    //init_rdf.Display({"el_final_P", "el_final_corr_P"}, 10)->Print();


    auto end = std::chrono::high_resolution_clock::now(); // END

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time of execution: " << elapsed.count() << " sec" << std::endl;

    return 0;
}
