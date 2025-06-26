import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import ROOT

# Create pic folder
output_dir = "pic"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# File groups
base_path = "/eos/home-p/pelai/HZa/ALP/gridpacks/check_LHE/run3/rootfile"
ma_0p1_0p9_files = [
    "ALP_M0p1.root", "ALP_M0p2.root", "ALP_M0p3.root", "ALP_M0p4.root", "ALP_M0p5.root",
    "ALP_M0p6.root", "ALP_M0p7.root", "ALP_M0p8.root", "ALP_M0p9.root"
]
ma_1_30_files = [
    "ALP_M1.root", "ALP_M2.root", "ALP_M3.root", "ALP_M4.root", "ALP_M5.root",
    "ALP_M6.root", "ALP_M7.root", "ALP_M8.root", "ALP_M9.root", "ALP_M10.root",
    "ALP_M15.root", "ALP_M20.root", "ALP_M25.root", "ALP_M30.root"
]

# Set global tick font size
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

# Colors
colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896'
]

def plot_gamma_dr_events_and_bins(file_list, ma_range, events_output_filename, events_extended_output_filename, bin_output_filename):
    gamma_dr_list = []
    total_valid_events_list = []
    labels = [f.replace("ALP_", "").replace(".root", "") for f in file_list]

    # Process files
    for filename in file_list:
        file_path = os.path.join(base_path, filename)
        gamma_dr = []
        try:
            with uproot.open(file_path) as file:
                tree = file["events"]
                px = tree["px"].array(library="np")
                py = tree["py"].array(library="np")
                pz = tree["pz"].array(library="np")
                energy = tree["energy"].array(library="np")
                total_events = len(px)
                print(f"Processing {filename}: {total_events} events")

                for i, (event_px, event_py, event_pz, event_energy) in enumerate(zip(px, py, pz, energy)):
                    if len(event_px) in [8, 9]:
                        v1 = ROOT.TLorentzVector()
                        v2 = ROOT.TLorentzVector()
                        v1.SetPxPyPzE(event_px[7], event_py[7], event_pz[7], event_energy[7])
                        if len(event_px) == 9:
                            v2.SetPxPyPzE(event_px[8], event_py[8], event_pz[8], event_energy[8])
                            gamma_dr.append(v1.DeltaR(v2))
                        else:
                            gamma_dr.append(np.nan)  # Events with 8 particles have no second photon
                    else:
                        print(f"Warning: Event {i} in {filename} has {len(event_px)} particles, skipping.")
                        continue
            gamma_dr_array = np.array(gamma_dr)
            total_valid_events = np.sum(~np.isnan(gamma_dr_array))  # Count events with valid gamma_dr
            print(f"Valid events with gamma_dr in {filename}: {total_valid_events}")
            gamma_dr_list.append(gamma_dr_array)
            total_valid_events_list.append(total_valid_events)
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            gamma_dr_list.append(np.array([]))
            total_valid_events_list.append(0)

    # 1. Line plot for efficiency with points (ΔR range 0 to 1.1)
    dr_bins = np.arange(0, 1.1, 0.1)  # 0 to 1 with step 0.1
    plt.figure(figsize=(10, 6))
    handles = []
    for i, (gamma_dr, label, total_valid) in enumerate(zip(gamma_dr_list, labels, total_valid_events_list)):
        if len(gamma_dr) > 0 and total_valid > 0:
            hist, bin_edges = np.histogram(gamma_dr, bins=dr_bins, density=False)
            events_remaining = total_valid - np.cumsum(hist)  # Events remaining after cut
            efficiency = events_remaining / total_valid  # Efficiency = remaining events / valid events
            # Prepend point at ΔR = 0 with efficiency = 1
            plot_x = np.concatenate(([0], bin_edges[1:]))
            plot_y = np.concatenate(([1], efficiency))
            plt.plot(plot_x, plot_y, marker='o', linestyle='-', color=colors[i], linewidth=2.0, markersize=6)
            handles.append(Line2D([0], [0], color=colors[i], linewidth=2.0, marker='o', markersize=6, label=label))
    plt.xlabel(r"$\Delta R(\gamma_1, \gamma_2)$ Cut", fontsize=24)
    plt.ylabel("Efficiency", fontsize=24)
    plt.legend(handles=handles, fontsize=14, ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, events_output_filename))
    plt.close()

    # 2. Line plot for efficiency with points (ΔR range 0 to 3.0)
    dr_bins_extended = np.arange(0, 3.1, 0.1)  # 0 to 3 with step 0.1
    plt.figure(figsize=(10, 6))
    handles = []
    for i, (gamma_dr, label, total_valid) in enumerate(zip(gamma_dr_list, labels, total_valid_events_list)):
        if len(gamma_dr) > 0 and total_valid > 0:
            hist, bin_edges = np.histogram(gamma_dr, bins=dr_bins_extended, density=False)
            events_remaining = total_valid - np.cumsum(hist)  # Events remaining after cut
            efficiency = events_remaining / total_valid  # Efficiency = remaining events / valid events
            # Prepend point at ΔR = 0 with efficiency = 1
            plot_x = np.concatenate(([0], bin_edges[1:]))
            plot_y = np.concatenate(([1], efficiency))
            plt.plot(plot_x, plot_y, marker='o', linestyle='-', color=colors[i], linewidth=2.0, markersize=6)
            handles.append(Line2D([0], [0], color=colors[i], linewidth=2.0, marker='o', markersize=6, label=label))
    plt.xlabel(r"$\Delta R(\gamma_1, \gamma_2)$ Cut", fontsize=24)
    plt.ylabel("Efficiency", fontsize=24)
    plt.legend(handles=handles, fontsize=14, ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, events_extended_output_filename))
    plt.close()

    # 3. Bar plot for proportions
    ranges = [(0, 0.1), (0.1, 0.3), (0.3, np.inf)]
    range_labels = ["ΔR ≤ 0.1", "0.1 < ΔR ≤ 0.3", "ΔR > 0.3"]
    proportions = []
    for gamma_dr in gamma_dr_list:
        if len(gamma_dr) > 0:
            valid_dr = gamma_dr[~np.isnan(gamma_dr)]
            props = []
            for r_min, r_max in ranges:
                count = np.sum((valid_dr > r_min) & (valid_dr <= r_max)) if r_max != np.inf else np.sum(valid_dr > r_min)
                props.append(count / len(valid_dr))
            proportions.append(props)
        else:
            proportions.append([0, 0, 0])

    proportions = np.array(proportions).T  # [range, file]

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(labels))
    width = 0.25
    for i, (props, label) in enumerate(zip(proportions, range_labels)):
        ax.bar(x + i * width, props, width, label=label, color=colors[i % len(colors)])
    ax.set_xlabel("ALP Mass (GeV)", fontsize=24)
    ax.set_ylabel("Proportion", fontsize=24)
    ax.set_xticks(x + width)
    ax.set_xticklabels(labels, rotation=45)
    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, bin_output_filename))
    plt.close()

# Generate plots
plot_gamma_dr_events_and_bins(
    ma_0p1_0p9_files,
    ma_range="0.1-0.9 GeV",
    events_output_filename="gamma_dr_efficiency_ma_0p1_0p9.pdf",
    events_extended_output_filename="gamma_dr_efficiency_extended_ma_0p1_0p9.pdf",
    bin_output_filename="gamma_dr_bins_ma_0p1_0p9.pdf"
)
plot_gamma_dr_events_and_bins(
    ma_1_30_files,
    ma_range="1-30 GeV",
    events_output_filename="gamma_dr_efficiency_ma_1_30.pdf",
    events_extended_output_filename="gamma_dr_efficiency_extended_ma_1_30.pdf",
    bin_output_filename="gamma_dr_bins_ma_1_30.pdf"
)