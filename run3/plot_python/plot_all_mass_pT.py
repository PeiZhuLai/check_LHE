import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import os

# 檢查並創建 pic 文件夾
output_dir = "pic"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 定義文件分組
base_path = "/eos/home-p/pelai/HZa/ALP/gridpacks/check_LHE/run3/rootfile"
ma_0p1_0p9_files = [
    "ALP_M0p1.root", "ALP_M0p2.root", "ALP_M0p3.root", "ALP_M0p5.root",
    "ALP_M0p6.root", "ALP_M0p7.root", "ALP_M0p8.root", "ALP_M0p9.root"
]
ma_1_30_files = [
    "ALP_M1.root", "ALP_M2.root", "ALP_M3.root", "ALP_M4.root", "ALP_M5.root",
    "ALP_M6.root", "ALP_M7.root", "ALP_M8.root", "ALP_M9.root", "ALP_M10.root",
    "ALP_M15.root", "ALP_M20.root", "ALP_M25.root", "ALP_M30.root"
]

# 設置全局刻度字體大小
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

# 定義顏色循環：使用 tab20
# colors = plt.cm.tab20(np.linspace(0, 1, 20))  # 20 種顏色

colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896'
]

# 定義繪圖函數
def plot_mass_and_pt_distributions(file_list, ma_range, mass_output_filename, pt_output_filename, alp_mass_range):
    # 初始化每個文件的質量和 pT 數據
    higgs_mass_list = []
    alp_mass_list = []
    z_mass_list = []
    electron_mass_list = []
    muon_mass_list = []
    tau_mass_list = []
    gamma_mass_list = []
    higgs_pt_list = []
    alp_pt_list = []
    z_pt_list = []
    electron_pt_list = []
    muon_pt_list = []
    tau_pt_list = []
    gamma_pt_list = []

    # 輕子質量範圍（電子和 μ子使用 MeV/c²，τ子使用 GeV/c²）
    electron_mass_range = (0.4, 0.6)       # ~0.511 MeV/c²
    muon_mass_range = (90, 120)            # ~105.66 MeV/c²
    tau_mass_range = (1.7, 1.9)            # ~1776.86 MeV/c² (GeV/c²)

    # 遍歷文件
    total_events = 0
    skipped_events = 0
    unknown_lepton_count = 0
    for filename in file_list:
        file_path = os.path.join(base_path, filename)
        # 為該文件初始化質量和 pT 列表
        higgs_mass = []
        alp_mass = []
        z_mass = []
        electron_mass = []
        muon_mass = []
        tau_mass = []
        gamma_mass = []
        higgs_pt = []
        alp_pt = []
        z_pt = []
        electron_pt = []
        muon_pt = []
        tau_pt = []
        gamma_pt = []
        try:
            with uproot.open(file_path) as file:
                tree = file["events"]
                masses = tree["mass"].array(library="np")
                px = tree["px"].array(library="np")
                py = tree["py"].array(library="np")
                total_events += len(masses)
                print(f"Processing {filename}: {len(masses)} events")

                # 遍歷每個事件
                for i, (event_masses, event_px, event_py) in enumerate(zip(masses, px, py)):
                    if len(event_masses) != len(event_px) or len(event_masses) != len(event_py):
                        print(f"Warning: Event {i} in {filename} has mismatched lengths: masses={len(event_masses)}, px={len(event_px)}, py={len(event_py)}")
                        skipped_events += 1
                        continue
                    if len(event_masses) == 9:
                        higgs_mass.append(event_masses[2])  # instance 2: Higgs
                        alp_mass.append(event_masses[3])    # instance 3: ALP
                        z_mass.append(event_masses[4])      # instance 4: Z
                        higgs_pt.append(np.sqrt(event_px[2]**2 + event_py[2]**2))  # Higgs pT
                        alp_pt.append(np.sqrt(event_px[3]**2 + event_py[3]**2))    # ALP pT
                        z_pt.append(np.sqrt(event_px[4]**2 + event_py[4]**2))      # Z pT
                        # 處理 instance 5 和 6（輕子質量和 pT）
                        for j in [5, 6]:
                            lepton_mass = event_masses[j]
                            lepton_pt = np.sqrt(event_px[j]**2 + event_py[j]**2)
                            if tau_mass_range[0] <= lepton_mass <= tau_mass_range[1]:
                                tau_mass.append(lepton_mass)  # τ子保持 GeV/c²
                                tau_pt.append(lepton_pt)
                            else:
                                lepton_mass_mev = lepton_mass * 1000  # 轉換為 MeV/c²
                                if electron_mass_range[0] <= lepton_mass_mev <= electron_mass_range[1]:
                                    electron_mass.append(lepton_mass_mev)
                                    electron_pt.append(lepton_pt)
                                elif muon_mass_range[0] <= lepton_mass_mev <= muon_mass_range[1]:
                                    muon_mass.append(lepton_mass_mev)
                                    muon_pt.append(lepton_pt)
                                else:
                                    print(f"Warning: Event {i} in {filename} has unknown lepton mass {lepton_mass} GeV/c² ({lepton_mass_mev} MeV/c²)")
                                    unknown_lepton_count += 1
                        gamma_mass.extend([event_masses[7], event_masses[8]])  # instance 7, 8: Gamma
                        gamma_pt.extend([np.sqrt(event_px[7]**2 + event_py[7]**2), np.sqrt(event_px[8]**2 + event_py[8]**2)])
                    elif len(event_masses) == 8:
                        print(f"Warning: Event {i} in {filename} has 8 masses: {event_masses}")
                        higgs_mass.append(event_masses[2])
                        alp_mass.append(event_masses[3])
                        z_mass.append(event_masses[4])
                        higgs_pt.append(np.sqrt(event_px[2]**2 + event_py[2]**2))
                        alp_pt.append(np.sqrt(event_px[3]**2 + event_py[3]**2))
                        z_pt.append(np.sqrt(event_px[4]**2 + event_py[4]**2))
                        for j in [5, 6]:
                            lepton_mass = event_masses[j]
                            lepton_pt = np.sqrt(event_px[j]**2 + event_py[j]**2)
                            if tau_mass_range[0] <= lepton_mass <= tau_mass_range[1]:
                                tau_mass.append(lepton_mass)
                                tau_pt.append(lepton_pt)
                            else:
                                lepton_mass_mev = lepton_mass * 1000
                                if electron_mass_range[0] <= lepton_mass_mev <= electron_mass_range[1]:
                                    electron_mass.append(lepton_mass_mev)
                                    electron_pt.append(lepton_pt)
                                elif muon_mass_range[0] <= lepton_mass_mev <= muon_mass_range[1]:
                                    muon_mass.append(lepton_mass_mev)
                                    muon_pt.append(lepton_pt)
                                else:
                                    print(f"Warning: Event {i} in {filename} has unknown lepton mass {lepton_mass} GeV/c² ({lepton_mass_mev} MeV/c²)")
                                    unknown_lepton_count += 1
                        gamma_mass.extend([event_masses[7], 0.0])
                        gamma_pt.extend([np.sqrt(event_px[7]**2 + event_py[7]**2), 0.0])
                    else:
                        print(f"Warning: Event {i} in {filename} has {len(event_masses)} masses, skipping.")
                        skipped_events += 1
            # 將該文件的數據添加到列表
            higgs_mass_list.append(np.array(higgs_mass))
            alp_mass_list.append(np.array(alp_mass))
            z_mass_list.append(np.array(z_mass))
            electron_mass_list.append(np.array(electron_mass))
            muon_mass_list.append(np.array(muon_mass))
            tau_mass_list.append(np.array(tau_mass))
            gamma_mass_list.append(np.array(gamma_mass))
            higgs_pt_list.append(np.array(higgs_pt))
            alp_pt_list.append(np.array(alp_pt))
            z_pt_list.append(np.array(z_pt))
            electron_pt_list.append(np.array(electron_pt))
            muon_pt_list.append(np.array(muon_pt))
            tau_pt_list.append(np.array(tau_pt))
            gamma_pt_list.append(np.array(gamma_pt))
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    # 打印總結
    print(f"Total events processed: {total_events}")
    print(f"Total skipped events: {skipped_events}")
    print(f"Total unknown lepton masses: {unknown_lepton_count}")

    # 繪製質量分布圖
    fig_mass, axes_mass = plt.subplots(3, 3, figsize=(24, 18))
    axes_mass = axes_mass.flatten()

    # Higgs mass
    handles_mass0 = []
    for i, (higgs_mass, filename) in enumerate(zip(higgs_mass_list, file_list)):
        if len(higgs_mass) > 0:
            axes_mass[0].hist(higgs_mass, bins=100, range=(120, 130), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass0.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[0].set_title("Higgs Mass Distribution", fontsize=18)
    axes_mass[0].set_xlabel("Mass (GeV)", fontsize=18)
    axes_mass[0].set_ylabel("A.U.", fontsize=18)
    axes_mass[0].legend(handles=handles_mass0, fontsize=14, handlelength=2)

    # ALP mass
    handles_mass1 = []
    for i, (alp_mass, filename) in enumerate(zip(alp_mass_list, file_list)):
        if len(alp_mass) > 0:
            axes_mass[1].hist(alp_mass, bins=100, range=alp_mass_range, histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass1.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[1].set_title(f"ALP Mass Distribution (ma {ma_range})", fontsize=18)
    axes_mass[1].set_xlabel("Mass (GeV)", fontsize=18)
    axes_mass[1].set_ylabel("A.U.", fontsize=18)
    axes_mass[1].legend(handles=handles_mass1, fontsize=14, handlelength=2)

    # Z mass
    handles_mass2 = []
    for i, (z_mass, filename) in enumerate(zip(z_mass_list, file_list)):
        if len(z_mass) > 0:
            axes_mass[2].hist(z_mass, bins=25, range=(70, 110), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass2.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[2].set_title("Z Mass Distribution", fontsize=18)
    axes_mass[2].set_xlabel("Mass (GeV)", fontsize=18)
    axes_mass[2].set_ylabel("A.U.", fontsize=18)
    axes_mass[2].legend(handles=handles_mass2, fontsize=14, handlelength=2)

    # Electron mass
    handles_mass3 = []
    for i, (electron_mass, filename) in enumerate(zip(electron_mass_list, file_list)):
        if len(electron_mass) > 0:
            axes_mass[3].hist(electron_mass, bins=100, range=(0.4, 0.6), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass3.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[3].set_title("Electron Mass Distribution", fontsize=18)
    axes_mass[3].set_xlabel("Mass (MeV)", fontsize=18)
    axes_mass[3].set_ylabel("A.U.", fontsize=18)
    axes_mass[3].legend(handles=handles_mass3, fontsize=14, handlelength=2)

    # Muon mass
    handles_mass4 = []
    for i, (muon_mass, filename) in enumerate(zip(muon_mass_list, file_list)):
        if len(muon_mass) > 0:
            axes_mass[4].hist(muon_mass, bins=100, range=(90, 120), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass4.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[4].set_title("Muon Mass Distribution", fontsize=18)
    axes_mass[4].set_xlabel("Mass (MeV)", fontsize=18)
    axes_mass[4].set_ylabel("A.U.", fontsize=18)
    axes_mass[4].legend(handles=handles_mass4, fontsize=14, handlelength=2)

    # Tau mass
    handles_mass5 = []
    for i, (tau_mass, filename) in enumerate(zip(tau_mass_list, file_list)):
        if len(tau_mass) > 0:
            axes_mass[5].hist(tau_mass, bins=100, range=(1.7, 1.9), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass5.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[5].set_title("Tau Mass Distribution", fontsize=18)
    axes_mass[5].set_xlabel("Mass (GeV)", fontsize=18)
    axes_mass[5].set_ylabel("A.U.", fontsize=18)
    axes_mass[5].legend(handles=handles_mass5, fontsize=14, handlelength=2)

    # Gamma mass
    handles_mass6 = []
    for i, (gamma_mass, filename) in enumerate(zip(gamma_mass_list, file_list)):
        if len(gamma_mass) > 0:
            axes_mass[6].hist(gamma_mass, bins=100, range=(0, 1), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_mass6.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_mass[6].set_title("Gamma Mass Distribution", fontsize=18)
    axes_mass[6].set_xlabel("Mass (GeV)", fontsize=18)
    axes_mass[6].set_ylabel("A.U.", fontsize=18)
    axes_mass[6].legend(handles=handles_mass6, fontsize=14, handlelength=2)

    # 移除空白子圖
    axes_mass[7].axis('off')
    axes_mass[8].axis('off')

    # 保存質量圖
    plt.figure(fig_mass.number)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, mass_output_filename))
    plt.close(fig_mass)

    # 繪製 pT 分布圖
    fig_pt, axes_pt = plt.subplots(3, 3, figsize=(24, 18))
    axes_pt = axes_pt.flatten()

    # Higgs pT
    handles_pt0 = []
    for i, (higgs_pt, filename) in enumerate(zip(higgs_pt_list, file_list)):
        if len(higgs_pt) > 0:
            axes_pt[0].hist(higgs_pt, bins=25, range=(0, 200), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt0.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[0].set_title("Higgs pT Distribution", fontsize=18)
    axes_pt[0].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[0].set_ylabel("A.U.", fontsize=18)
    axes_pt[0].legend(handles=handles_pt0, fontsize=14, handlelength=2)

    # ALP pT
    handles_pt1 = []
    for i, (alp_pt, filename) in enumerate(zip(alp_pt_list, file_list)):
        if len(alp_pt) > 0:
            axes_pt[1].hist(alp_pt, bins=25, range=(0, 50), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt1.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[1].set_title(f"ALP pT Distribution (ma {ma_range})", fontsize=18)
    axes_pt[1].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[1].set_ylabel("A.U.", fontsize=18)
    axes_pt[1].legend(handles=handles_pt1, fontsize=14, handlelength=2)

    # Z pT
    handles_pt2 = []
    for i, (z_pt, filename) in enumerate(zip(z_pt_list, file_list)):
        if len(z_pt) > 0:
            axes_pt[2].hist(z_pt, bins=25, range=(0, 80), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt2.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[2].set_title("Z pT Distribution", fontsize=18)
    axes_pt[2].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[2].set_ylabel("A.U.", fontsize=18)
    axes_pt[2].legend(handles=handles_pt2, fontsize=14, handlelength=2)

    # Electron pT
    handles_pt3 = []
    for i, (electron_pt, filename) in enumerate(zip(electron_pt_list, file_list)):
        if len(electron_pt) > 0:
            axes_pt[3].hist(electron_pt, bins=25, range=(0, 90), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt3.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[3].set_title("Electron pT Distribution", fontsize=18)
    axes_pt[3].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[3].set_ylabel("A.U.", fontsize=18)
    axes_pt[3].legend(handles=handles_pt3, fontsize=14, handlelength=2)

    # Muon pT
    handles_pt4 = []
    for i, (muon_pt, filename) in enumerate(zip(muon_pt_list, file_list)):
        if len(muon_pt) > 0:
            axes_pt[4].hist(muon_pt, bins=25, range=(0, 90), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt4.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[4].set_title("Muon pT Distribution", fontsize=18)
    axes_pt[4].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[4].set_ylabel("A.U.", fontsize=18)
    axes_pt[4].legend(handles=handles_pt4, fontsize=14, handlelength=2)

    # Tau pT
    handles_pt5 = []
    for i, (tau_pt, filename) in enumerate(zip(tau_pt_list, file_list)):
        if len(tau_pt) > 0:
            axes_pt[5].hist(tau_pt, bins=25, range=(0, 90), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt5.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[5].set_title("Tau pT Distribution", fontsize=18)
    axes_pt[5].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[5].set_ylabel("A.U.", fontsize=18)
    axes_pt[5].legend(handles=handles_pt5, fontsize=14, handlelength=2)

    # Gamma pT
    handles_pt6 = []
    for i, (gamma_pt, filename) in enumerate(zip(gamma_pt_list, file_list)):
        if len(gamma_pt) > 0:
            axes_pt[6].hist(gamma_pt, bins=25, range=(0, 30), histtype='step', density=True, color=colors[i], linewidth=1.5)
            handles_pt6.append(Line2D([0], [0], color=colors[i], linewidth=1.5, label=filename.replace(".root", "")))
    axes_pt[6].set_title("Gamma pT Distribution", fontsize=18)
    axes_pt[6].set_xlabel("pT (GeV)", fontsize=18)
    axes_pt[6].set_ylabel("A.U.", fontsize=18)
    axes_pt[6].legend(handles=handles_pt6, fontsize=14, handlelength=2)

    # 移除空白子圖
    axes_pt[7].axis('off')
    axes_pt[8].axis('off')

    # 保存 pT 圖
    plt.figure(fig_pt.number)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pt_output_filename))
    plt.close(fig_pt)

# 繪製質量和 pT 圖
plot_mass_and_pt_distributions(
    ma_0p1_0p9_files,
    ma_range="0.1-0.9 GeV",
    mass_output_filename="mass_distributions_ma_0p1_0p9.pdf",
    pt_output_filename="pt_distributions_ma_0p1_0p9.pdf",
    alp_mass_range=(0, 1)
)
plot_mass_and_pt_distributions(
    ma_1_30_files,
    ma_range="1-30 GeV",
    mass_output_filename="mass_distributions_ma_1_30.pdf",
    pt_output_filename="pt_distributions_ma_1_30.pdf",
    alp_mass_range=(0, 35)
)