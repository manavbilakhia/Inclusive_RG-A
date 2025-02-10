import ROOT
import math

M_P = 0.938  # Proton mass in GeV/c^2
E_BEAM = 10.6  # Beam energy in GeV
OUTPUT_FOLDER = "/w/hallb-scshelf2102/clas12/bulgakov/projects/Inclusive_RG-A/analysis_out/"

def main():
    root_file_path = "/w/hallb-scshelf2102/clas12/bulgakov/projects/Inclusive_RG-A/data/outH2R_test/allRunsP1NickPart_2023.dat_QADBtest.root"
    rdf = convert_ttrees_to_rdataframe(root_file_path)
    if rdf is None:
        print("Error: Could not create RDataFrame.")
        return
    
   
    # Use Define to store both W and QSquared in the RDataFrame
    
    rdf = rdf.Define("el_final","return (TLorentzVector) {p4_ele_px[0],p4_ele_py[0],p4_ele_pz[0],p4_ele_E[0]};")
    
    rdf = rdf.Define("el_initial","return (TLorentzVector) {0,0,10.6,10.6};")
    
    rdf = rdf.Define("proton_initial","return (TLorentzVector) {0,0,0,0.938};")
    
    rdf = rdf.Define("Q2", "-(el_initial - el_final).M2()")
    
    rdf = rdf.Define("W", "(el_initial - el_final + proton_initial).M()")

    print("Columns in RDataFrame:", rdf.GetColumnNames())
    
    # # Retrieve the vector of TLorentzVectors
    # electron_vectors = rdf.Take[ROOT.TLorentzVector]("el_final")

    # # Convert it to a standard Python list
    # electron_list = electron_vectors.GetValue()

    # # Print the first 10 elements
    # for i, vec in enumerate(electron_list[:10]):
    #     print(f"Element {i}: (E={vec.E()}, px={vec.Px()}, py={vec.Py()}, pz={vec.Pz()})")
    

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
    canvas = ROOT.TCanvas("c3", "Q^{2} Distribution", 800, 600)
    
    hist = rdf.Histo1D(("Q^{2}", "Momentum Transfer Q^2; Q^2 (GeV^2); Events", 100, 0, 10), "Q2")
    hist.Draw()
    
    canvas.SaveAs(OUTPUT_FOLDER+"Q2_distribution.png")
    print("Saved 1D histogram as QSquared_distribution.png")

def plot_2d_W_vs_QSquared(rdf):
    """
    Creates and saves a 2D histogram of W vs QSquared.
    """
    canvas = ROOT.TCanvas("c4", "W vs QSquared", 800, 600)
    
    hist2D = rdf.Histo2D(("W_vs_Q2", "W vs Q^{2}; W (GeV); Q^2 (GeV^2)", 100, 0, 3, 100, 0, 10), "W", "Q2")
    hist2D.Draw("COLZ")  # COLZ option shows the color gradient
    
    canvas.SaveAs(OUTPUT_FOLDER+"W_vs_Q2.png")
    print("Saved 2D histogram as W_vs_Q2.png")

if __name__ == "__main__":
    main()
