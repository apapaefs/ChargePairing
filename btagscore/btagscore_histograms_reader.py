import numpy as np
import ROOT
import matplotlib.pyplot as plt

###################################
# Read the btagscore histograms.  # 
###################################

# the bin edges have been constructed with:
# const vector<int> pt_bins = {20, 30, 40, 50, 60, 70, 80, 90, 100,
#                                         120, 140, 160, 180, 200,
#                                         250, 300, 400, 600, 1000};
#    const vector<double> abs_eta_bins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
# see: https://github.com/jiashu-huang/belphes/blob/master/belphes-examples/btagscore_analysis/PlotHistogram.cpp
# The naming convention of each histogram in the .root file is
#hist_eta[X]_pt[Y]
# where X and Y indicate the n-th number of bins.

# function to read the btagscore histograms from a ROOT file
def read_btagscore_histograms(root_file_path):
    # Define the pt and abs(eta) bin edges
    pt_bins = [20, 30, 40, 50, 60, 70, 80, 90, 100,
               120, 140, 160, 180, 200,
               250, 300, 400, 600, 1000]
    abs_eta_bins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]

    # Open the ROOT file
    root_file = ROOT.TFile.Open(root_file_path)

    # Initialize a dictionary to hold the histograms
    btagscore_histograms = {}

    # Loop over eta and pt bins to read histograms
    for i in range(len(abs_eta_bins) - 1):
        for j in range(len(pt_bins) - 1):
            hist_name = f"hist_eta{i}_pt{j}"
            histogram = root_file.Get(hist_name)
            if histogram:
                # Convert ROOT histogram to numpy array
                bin_edges = [histogram.GetBinLowEdge(k) for k in range(1, histogram.GetNbinsX() + 2)]
                bin_contents = [histogram.GetBinContent(k) for k in range(1, histogram.GetNbinsX() + 1)]
                btagscore_histograms[(i, j)] = (np.array(bin_edges), np.array(bin_contents))
            else:
                print(f"Warning: Histogram {hist_name} not found in the ROOT file.")

    # Close the ROOT file
    root_file.Close()

    return btagscore_histograms

# function that returns the btagscore for given pt and abs(eta) of a jet:
def jet_btagscore(jet_pt, jet_abs_eta, btagscore_histograms):
    # Define the pt and abs(eta) bin edges
    pt_bins = [20, 30, 40, 50, 60, 70, 80, 90, 100,
               120, 140, 160, 180, 200,
               250, 300, 400, 600, 1000]
    abs_eta_bins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]

    # Find the appropriate bin indices for pt and abs(eta)
    pt_bin_index = next((i for i in range(len(pt_bins) - 1) if pt_bins[i] <= jet_pt < pt_bins[i + 1]), None)
    eta_bin_index = next((i for i in range(len(abs_eta_bins) - 1) if abs_eta_bins[i] <= jet_abs_eta < abs_eta_bins[i + 1]), None)

    if pt_bin_index is None or eta_bin_index is None:
        raise ValueError("Jet pt or abs(eta) is out of the defined bin ranges.")

    # Retrieve the corresponding histogram
    bin_edges, bin_contents = btagscore_histograms.get((eta_bin_index, pt_bin_index), (None, None))
    if bin_edges is None or bin_contents is None:
        raise ValueError("No histogram found for the given pt and abs(eta) bin.")

    # Sample a btagscore from the histogram
    btagscore = np.random.choice(bin_edges[:-1], p=bin_contents / np.sum(bin_contents))
    return btagscore

# plot all histograms (for testing purposes)
def plot_all_histograms(btagscore_histograms, b_or_nonb='B'):
    for (eta_idx, pt_idx), (bin_edges, bin_contents) in btagscore_histograms.items():
        plt.figure()
        plt.hist(bin_edges[:-1], bins=bin_edges, weights=bin_contents, histtype='stepfilled', alpha=0.7)
        plt.title(f'Histogram for eta bin {eta_idx}, pt bin {pt_idx}')
        plt.xlabel('B-tag Score')
        plt.ylabel('Entries')
        plt.grid()
        # print pdfs in subdirectory "histograms_plots"
        plt.savefig(f'histograms_plots/' + b_or_nonb + f'_hist_eta{eta_idx}_pt{pt_idx}.pdf')
        plt.close()

############################################## 
# Read the ROOT files and use the functions: #
##############################################

# Read the histograms for B and Non-B jets:
# B jets
file_path_B = "JetBtagDeepFlavB_B_Distributions.root"
JetBtagDeepFlavB_B_Distributions = read_btagscore_histograms(file_path_B)
# non-B jets
file_path_NonB = "JetBtagDeepFlavB_B_Distributions.root"
JetBtagDeepFlavB_NonB_Distributions = read_btagscore_histograms(file_path_NonB)

# plot all histograms:
plot_all_histograms(JetBtagDeepFlavB_B_Distributions,b_or_nonb='B')
plot_all_histograms(JetBtagDeepFlavB_NonB_Distributions,b_or_nonb='NonB')


# ###################################################### 
# Example usage on a jet                               #
########################################################

jet_pt = 85
jet_abs_eta = 1.2
# if B-jet:
btagscore_B = jet_btagscore(jet_pt, jet_abs_eta, JetBtagDeepFlavB_B_Distributions)
# if non-B-jet:
btagscore_nonB = jet_btagscore(jet_pt, jet_abs_eta, JetBtagDeepFlavB_NonB_Distributions)

print("Example: For a jet with:")
print("Selected random btagscore from the histograms:")
print(f"Jet pT: {jet_pt} GeV, abs(eta): {jet_abs_eta}")
print(f"B-jet btagscore: {btagscore_B}")
print(f"Non-B-jet btagscore: {btagscore_nonB}")
