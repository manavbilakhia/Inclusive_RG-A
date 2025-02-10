import ROOT
import math

M_P = 0.938  # Proton mass in GeV/c^2
E_BEAM = 10.6  # Beam energy in GeV
OUTPUT_FOLDER = "/w/hallb-scshelf2102/clas12/manavb/grad/Inclusive_RG-A/analysis_out/"

def main():
    root_file_path = "/w/hallb-scshelf2102/clas12/manavb/grad/Inclusive_RG-A/data/outH2R_test/allRunsP1NickPart_2023.dat_QADBtest.root"
    
    rdf = convert_ttrees_to_rdataframe(root_file_path)
    if rdf is None:
        print("Error: Could not create RDataFrame.")
        return
    
    # Declare vectorized function in ROOT
    ROOT.gInterpreter.Declare("""
    #include <vector>
    #include <cmath>
    #include <utility>
    using namespace ROOT::VecOps;

    std::pair<RVec<double>, RVec<double>> calc_W_Q2(const RVec<double>& p4_ele_px, 
                                                    const RVec<double>& p4_ele_py, 
                                                    const RVec<double>& p4_ele_pz, 
                                                    const RVec<double>& p4_ele_E) 
    {
        RVec<double> W;
        RVec<double> QSquared;
        double M_P = 0.938;
        double E_BEAM = 10.6;

        for (size_t i = 0; i < p4_ele_px.size(); ++i) {
            double qx = -p4_ele_px[i];
            double qy = -p4_ele_py[i];
            double qz = E_BEAM - p4_ele_pz[i];
            double q0 = E_BEAM - p4_ele_E[i];

            double Q2 = - q0 * q0 + (qx * qx + qy * qy + qz * qz);
            double W2 = (M_P * M_P + 2 * M_P * q0 - Q2);
            
            W.push_back(W2 > 0 ? sqrt(W2) : 0);
            QSquared.push_back(Q2);
        }
        return std::make_pair(W, QSquared);
    }
    """)

    # Use Define to store both W and QSquared in the RDataFrame
    rdf = rdf.Define("W_Q2", "calc_W_Q2(p4_ele_px, p4_ele_py, p4_ele_pz, p4_ele_E)")
    rdf = rdf.Define("W", "W_Q2.first")       # Extract W from the pair
    rdf = rdf.Define("QSquared", "W_Q2.second")  # Extract QSquared from the pair

    print("Columns in RDataFrame:", rdf.GetColumnNames())

    # Plot W, QSquared, and W vs QSquared
    plot_1d_W(rdf)
    plot_1d_QSquared(rdf)
    plot_2d_W_vs_QSquared(rdf)

def convert_ttrees_to_rdataframe(root_file_path):
    file = ROOT.TFile.Open(root_file_path, "READ")
    if not file or file.IsZombie():
        print(f"Error: Cannot open ROOT file {root_file_path}")
        return None

    keys = [key.GetName() for key in file.GetListOfKeys() if key.GetClassName() == "TTree"]
    if not keys:
        print("No TTrees found in the ROOT file.")
        return None

    tree_name = keys[0]
    print(f"Processing TTree: {tree_name}")
    
    rdf = ROOT.RDataFrame(tree_name, root_file_path)
    file.Close()
    
    return rdf

def plot_1d_W(rdf):
    """
    Creates and saves a 1D histogram of W.
    """
    canvas = ROOT.TCanvas("c1", "W Distribution", 800, 600)
    
    hist = rdf.Histo1D(("W_distribution", "Invariant Mass W; W (GeV); Events", 100, 0, 5), "W")
    hist.Draw()
    
    canvas.SaveAs(OUTPUT_FOLDER+"W_distribution.png")
    print("Saved 1D histogram as W_distribution.png")

def plot_1d_QSquared(rdf):
    """
    Creates and saves a 1D histogram of QSquared.
    """
    canvas = ROOT.TCanvas("c3", "QSquared Distribution", 800, 600)
    
    hist = rdf.Histo1D(("QSquared_distribution", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "QSquared")
    hist.Draw()
    
    canvas.SaveAs(OUTPUT_FOLDER+"QSquared_distribution.png")
    print("Saved 1D histogram as QSquared_distribution.png")

def plot_2d_W_vs_QSquared(rdf):
    """
    Creates and saves a 2D histogram of W vs QSquared.
    """
    canvas = ROOT.TCanvas("c4", "W vs QSquared", 800, 600)
    
    hist2D = rdf.Histo2D(("W_vs_QSquared", "W vs QSquared; W (GeV); Q^2 (GeV^2)", 100, 0, 3, 100, 0, 10), "W", "QSquared")
    hist2D.Draw("COLZ")  # COLZ option shows the color gradient
    
    canvas.SaveAs("W_vs_QSquared.png")
    print(OUTPUT_FOLDER+"Saved 2D histogram as W_vs_QSquared.png")

if __name__ == "__main__":
    main()
