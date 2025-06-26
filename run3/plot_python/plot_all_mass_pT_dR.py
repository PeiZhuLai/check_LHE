import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import os
import ROOT  # 導入 PyROOT 用於 TLorentzVector

# 檢查並創建 pic 文件夾
output_dir = "pic"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 定義文件分組
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

# 設置全局刻度字體大小
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20

# 定義顏色
colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896'
]

# 定義繪圖函數
def plot_mass_and_pt_distributions(file_list, ma_range, mass_output_filename, pt_output_filename, dr_output_filename, alp_mass_range):
    # 初始化每個文件的質量、pT 和 ΔR 數據
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
    gamma_dr_list = []
    electron_dr_list = []
    muon_dr_list = []

    # 輕子質量範圍（電子和 μ子使用 MeV/c²，τ子使用 GeV/c²）
    electron_mass_range = (0.4, 0.6)  # ~0.511 MeV/c²
    muon_mass_range = (90, 120)       # ~105.66 MeV/c²
    tau_mass_range = (1.7, 1.9)       # ~1776.86 MeV/c² (GeV/c²)

    # 遍歷文件
    total_events = 0
    skipped_events = 0
    unknown_lepton_count = 0
    for filename in file_list:
        file_path = os.path.join(base_path, filename)
        # 為該文件初始化數據列表
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
        gamma_dr = []
        electron_dr = []
        muon_dr = []
        try:
            with uproot.open(file_path) as file:
                tree = file["events"]
                masses = tree["mass"].array(library="np")
                px = tree["px"].array(library="np")
                py = tree["py"].array(library="np")
                pz = tree["pz"].array(library="np")
                energy = tree["energy"].array(library="np")  # 讀取 energy 分支
                total_events += len(masses)
                print(f"Processing {filename}: {len(masses)} events")

                # 遍歷每個事件
                for i, (event_masses, event_px, event_py, event_pz, event_energy) in enumerate(zip(masses, px, py, pz, energy)):
                    if len(event_masses) != len(event_px) or len(event_masses) != len(event_py) or len(event_masses) != len(event_pz) or len(event_masses) != len(event_energy):
                        print(f"Warning: Event {i} in {filename} has mismatched lengths: masses={len(event_masses)}, px={len(event_px)}, py={len(event_py)}, pz={len(event_pz)}, energy={len(event_energy)}")
                        skipped_events += 1
                        continue
                    if len(event_masses) == 9:
                        # 質量和 pT 數據
                        higgs_mass.append(event_masses[2])  # instance 2: Higgs
                        alp_mass.append(event_masses[3])    # instance 3: ALP
                        z_mass.append(event_masses[4])      # instance 4: Z
                        higgs_pt.append(np.sqrt(event_px[2]**2 + event_py[2]**2))
                        alp_pt.append(np.sqrt(event_px[3]**2 + event_py[3]**2))
                        z_pt.append(np.sqrt(event_px[4]**2 + event_py[4]**2))
                        # 處理 instance 5 和 6（輕子）
                        lepton_types = []
                        for j in [5, 6]:
                            lepton_mass = event_masses[j]
                            lepton_pt = np.sqrt(event_px[j]**2 + event_py[j]**2)
                            if tau_mass_range[0] <= lepton_mass <= tau_mass_range[1]:
                                tau_mass.append(lepton_mass)
                                tau_pt.append(lepton_pt)
                                lepton_types.append("tau")
                            else:
                                lepton_mass_mev = lepton_mass * 1000
                                if electron_mass_range[0] <= lepton_mass_mev <= electron_mass_range[1]:
                                    electron_mass.append(lepton_mass_mev)
                                    electron_pt.append(lepton_pt)
                                    lepton_types.append("electron")
                                elif muon_mass_range[0] <= lepton_mass_mev <= muon_mass_range[1]:
                                    muon_mass.append(lepton_mass_mev)
                                    muon_pt.append(lepton_pt)
                                    lepton_types.append("muon")
                                else:
                                    print(f"Warning: Event {i} in {filename} has unknown lepton mass {lepton_mass} GeV/c² ({lepton_mass_mev} MeV/c²)")
                                    unknown_lepton_count += 1
                                    lepton_types.append("unknown")
                        # Gamma 質量和 pT
                        gamma_mass.extend([event_masses[7], event_masses[8]])
                        gamma_pt.extend([np.sqrt(event_px[7]**2 + event_py[7]**2), np.sqrt(event_px[8]**2 + event_py[8]**2)])
                        # 計算 Gamma ΔR (instance 7 和 8)
                        v1 = ROOT.TLorentzVector()
                        v2 = ROOT.TLorentzVector()
                        v1.SetPxPyPzE(event_px[7], event_py[7], event_pz[7], event_energy[7])
                        v2.SetPxPyPzE(event_px[8], event_py[8], event_pz[8], event_energy[8])
                        gamma_dr.append(v1.DeltaR(v2))
                        # 計算 Electron 或 Muon ΔR
                        if lepton_types == ["electron", "electron"]:
                            v1.SetPxPyPzE(event_px[5], event_py[5], event_pz[5], event_energy[5])
                            v2.SetPxPyPzE(event_px[6], event_py[6], event_pz[6], event_energy[6])
                            electron_dr.append(v1.DeltaR(v2))
                        elif lepton_types == ["muon", "muon"]:
                            v1.SetPxPyPzE(event_px[5], event_py[5], event_pz[5], event_energy[5])
                            v2.SetPxPyPzE(event_px[6], event_py[6], event_pz[6], event_energy[6])
                            muon_dr.append(v1.DeltaR(v2))
                    elif len(event_masses) == 8:
                        print(f"Warning: Event {i} in {filename} has 8 masses: {event_masses}")
                        higgs_mass.append(event_masses[2])
                        alp_mass.append(event_masses[3])
                        z_mass.append(event_masses[4])
                        higgs_pt.append(np.sqrt(event_px[2]**2 + event_py[2]**2))
                        alp_pt.append(np.sqrt(event_px[3]**2 + event_py[3]**2))
                        z_pt.append(np.sqrt(event_px[4]**2 + event_py[4]**2))
                        lepton_types = []
                        for j in [5, 6]:
                            lepton_mass = event_masses[j]
                            lepton_pt = np.sqrt(event_px[j]**2 + event_py[j]**2)
                            if tau_mass_range[0] <= lepton_mass <= tau_mass_range[1]:
                                tau_mass.append(lepton_mass)
                                tau_pt.append(lepton_pt)
                                lepton_types.append("tau")
                            else:
                                lepton_mass_mev = lepton_mass * 1000
                                if electron_mass_range[0] <= lepton_mass_mev <= electron_mass_range[1]:
                                    electron_mass.append(lepton_mass_mev)
                                    electron_pt.append(lepton_pt)
                                    lepton_types.append("electron")
                                elif muon_mass_range[0] <= lepton_mass_mev <= muon_mass_range[1]:
                                    muon_mass.append(lepton_mass_mev)
                                    muon_pt.append(lepton_pt)
                                    lepton_types.append("muon")
                                else:
                                    print(f"Warning: Event {i} in {filename} has unknown lepton mass {lepton_mass} GeV/c² ({lepton_mass_mev} MeV/c²)")
                                    unknown_lepton_count += 1
                                    lepton_types.append("unknown")
                        gamma_mass.extend([event_masses[7], 0.0])
                        gamma_pt.extend([np.sqrt(event_px[7]**2 + event_py[7]**2), 0.0])
                        gamma_dr.append(np.nan)  # 8 粒子事件無第二個光子
                        if lepton_types == ["electron", "electron"]:
                            v1 = ROOT.TLorentzVector()
                            v2 = ROOT.TLorentzVector()
                            v1.SetPxPyPzE(event_px[5], event_py[5], event_pz[5], event_energy[5])
                            v2.SetPxPyPzE(event_px[6], event_py[6], event_pz[6], event_energy[6])
                            electron_dr.append(v1.DeltaR(v2))
                        elif lepton_types == ["muon", "muon"]:
                            v1 = ROOT.TLorentzVector()
                            v2 = ROOT.TLorentzVector()
                            v1.SetPxPyPzE(event_px[5], event_py[5], event_pz[5], event_energy[5])
                            v2.SetPxPyPzE(event_px[6], event_py[6], event_pz[6], event_energy[6])
                            muon_dr.append(v1.DeltaR(v2))
                    else:
                        print(f"Warning: Event {i} in {filename} has {len(event_masses)} masses, skipping.")
                        skipped_events += 1
                        continue
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
            gamma_dr_list.append(np.array(gamma_dr))
            electron_dr_list.append(np.array(electron_dr))
            muon_dr_list.append(np.array(muon_dr))
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
            axes_mass[0].hist(higgs_mass, bins=100, range=(120, 130), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass0.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[0].set_title("Higgs Mass Distribution", fontsize=18)
    axes_mass[0].set_xlabel("Higgs Mass (GeV)", fontsize=24)
    axes_mass[0].set_ylabel("A.U.", fontsize=24)
    axes_mass[0].legend(handles=handles_mass0, fontsize=14, handlelength=2, ncol=2)

    # ALP mass
    handles_mass1 = []
    for i, (alp_mass, filename) in enumerate(zip(alp_mass_list, file_list)):
        if len(alp_mass) > 0:
            axes_mass[1].hist(alp_mass, bins=100, range=alp_mass_range, histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass1.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[1].set_title(f"ALP Mass Distribution (ma {ma_range})", fontsize=18)
    axes_mass[1].set_xlabel("ALP Mass (GeV)", fontsize=24)
    axes_mass[1].set_ylabel("A.U.", fontsize=24)
    axes_mass[1].legend(handles=handles_mass1, fontsize=14, handlelength=2, ncol=2)

    # Z mass
    handles_mass2 = []
    for i, (z_mass, filename) in enumerate(zip(z_mass_list, file_list)):
        if len(z_mass) > 0:
            axes_mass[2].hist(z_mass, bins=25, range=(70, 110), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass2.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[2].set_title("Z Mass Distribution", fontsize=18)
    axes_mass[2].set_xlabel("Z boson Mass (GeV)", fontsize=24)
    axes_mass[2].set_ylabel("A.U.", fontsize=24)
    axes_mass[2].legend(handles=handles_mass2, fontsize=14, handlelength=2, ncol=2)

    # Electron mass
    handles_mass3 = []
    for i, (electron_mass, filename) in enumerate(zip(electron_mass_list, file_list)):
        if len(electron_mass) > 0:
            axes_mass[3].hist(electron_mass, bins=100, range=(0.4, 0.6), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass3.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[3].set_title("Electron Mass Distribution", fontsize=18)
    axes_mass[3].set_xlabel("e Mass (MeV)", fontsize=24)
    axes_mass[3].set_ylabel("A.U.", fontsize=24)
    axes_mass[3].legend(handles=handles_mass3, fontsize=14, handlelength=2, ncol=2)

    # Muon mass
    handles_mass4 = []
    for i, (muon_mass, filename) in enumerate(zip(muon_mass_list, file_list)):
        if len(muon_mass) > 0:
            axes_mass[4].hist(muon_mass, bins=100, range=(90, 120), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass4.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[4].set_title("Muon Mass Distribution", fontsize=18)
    axes_mass[4].set_xlabel(r"$\mu \ Mass (MeV)$", fontsize=24)
    axes_mass[4].set_ylabel("A.U.", fontsize=24)
    axes_mass[4].legend(handles=handles_mass4, fontsize=14, handlelength=2, ncol=2)

    # Tau mass
    handles_mass5 = []
    for i, (tau_mass, filename) in enumerate(zip(tau_mass_list, file_list)):
        if len(tau_mass) > 0:
            axes_mass[5].hist(tau_mass, bins=100, range=(1.7, 1.9), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass5.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[5].set_title("Tau Mass Distribution", fontsize=24)
    axes_mass[5].set_xlabel(r"$\tau \ Mass (GeV)$", fontsize=24)
    axes_mass[5].set_ylabel("A.U.", fontsize=18)
    axes_mass[5].legend(handles=handles_mass5, fontsize=14, handlelength=2, ncol=2)

    # Gamma mass
    handles_mass6 = []
    for i, (gamma_mass, filename) in enumerate(zip(gamma_mass_list, file_list)):
        if len(gamma_mass) > 0:
            axes_mass[6].hist(gamma_mass, bins=100, range=(0, 1), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_mass6.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_mass[6].set_title("Gamma Mass Distribution", fontsize=18)
    axes_mass[6].set_xlabel(r"$\gamma \ Mass (GeV)$", fontsize=24)
    axes_mass[6].set_ylabel("A.U.", fontsize=24)
    axes_mass[6].legend(handles=handles_mass6, fontsize=14, handlelength=2, ncol=2)

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
            axes_pt[0].hist(higgs_pt, bins=25, range=(0, 200), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt0.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[0].set_title("Higgs pT Distribution", fontsize=18)
    axes_pt[0].set_xlabel("Higgs pT (GeV)", fontsize=24)
    axes_pt[0].set_ylabel("A.U.", fontsize=24)
    axes_pt[0].legend(handles=handles_pt0, fontsize=14, handlelength=2, ncol=2)

    # ALP pT
    handles_pt1 = []
    for i, (alp_pt, filename) in enumerate(zip(alp_pt_list, file_list)):
        if len(alp_pt) > 0:
            axes_pt[1].hist(alp_pt, bins=25, range=(0, 50), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt1.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[1].set_title(f"ALP pT Distribution (ma {ma_range})", fontsize=18)
    axes_pt[1].set_xlabel("ALP pT (GeV)", fontsize=24)
    axes_pt[1].set_ylabel("A.U.", fontsize=24)
    axes_pt[1].legend(handles=handles_pt1, fontsize=14, handlelength=2, ncol=2)

    # Z pT
    handles_pt2 = []
    for i, (z_pt, filename) in enumerate(zip(z_pt_list, file_list)):
        if len(z_pt) > 0:
            axes_pt[2].hist(z_pt, bins=25, range=(0, 80), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt2.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[2].set_title("Z pT Distribution", fontsize=18)
    axes_pt[2].set_xlabel("Z boson pT (GeV)", fontsize=24)
    axes_pt[2].set_ylabel("A.U.", fontsize=24)
    axes_pt[2].legend(handles=handles_pt2, fontsize=14, handlelength=2, ncol=2)

    # Electron pT
    handles_pt3 = []
    for i, (electron_pt, filename) in enumerate(zip(electron_pt_list, file_list)):
        if len(electron_pt) > 0:
            axes_pt[3].hist(electron_pt, bins=25, range=(0, 90), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt3.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[3].set_title("Electron pT Distribution", fontsize=18)
    axes_pt[3].set_xlabel("e pT (GeV)", fontsize=24)
    axes_pt[3].set_ylabel("A.U.", fontsize=24)
    axes_pt[3].legend(handles=handles_pt3, fontsize=14, handlelength=2, ncol=2)

    # Muon pT
    handles_pt4 = []
    for i, (muon_pt, filename) in enumerate(zip(muon_pt_list, file_list)):
        if len(muon_pt) > 0:
            axes_pt[4].hist(muon_pt, bins=25, range=(0, 90), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt4.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[4].set_title("Muon pT Distribution", fontsize=18)
    axes_pt[4].set_xlabel(r"$\mu \ pT (GeV)$", fontsize=24)
    axes_pt[4].set_ylabel("A.U.", fontsize=24)
    axes_pt[4].legend(handles=handles_pt4, fontsize=14, handlelength=2, ncol=2)

    # Tau pT
    handles_pt5 = []
    for i, (tau_pt, filename) in enumerate(zip(tau_pt_list, file_list)):
        if len(tau_pt) > 0:
            axes_pt[5].hist(tau_pt, bins=25, range=(0, 90), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt5.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[5].set_title("Tau pT Distribution", fontsize=18)
    axes_pt[5].set_xlabel(r"$\tau \ pT (GeV)$", fontsize=24)
    axes_pt[5].set_ylabel("A.U.", fontsize=24)
    axes_pt[5].legend(handles=handles_pt5, fontsize=14, handlelength=2, ncol=2)

    # Gamma pT
    handles_pt6 = []
    for i, (gamma_pt, filename) in enumerate(zip(gamma_pt_list, file_list)):
        if len(gamma_pt) > 0:
            axes_pt[6].hist(gamma_pt, bins=25, range=(0, 30), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_pt6.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_pt[6].set_title("Gamma pT Distribution", fontsize=18)
    axes_pt[6].set_xlabel(r"$\gamma \ pT (GeV)$", fontsize=24)
    axes_pt[6].set_ylabel("A.U.", fontsize=24)
    axes_pt[6].legend(handles=handles_pt6, fontsize=14, handlelength=2, ncol=2)

    # 移除空白子圖
    axes_pt[7].axis('off')
    axes_pt[8].axis('off')

    # 保存 pT 圖
    plt.figure(fig_pt.number)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, pt_output_filename))
    plt.close(fig_pt)

    # 繪製 ΔR 分布圖
    fig_dr, axes_dr = plt.subplots(1, 3, figsize=(24, 6))
    axes_dr = axes_dr.flatten()

    # Gamma ΔR
    handles_dr0 = []
    for i, (gamma_dr, filename) in enumerate(zip(gamma_dr_list, file_list)):
        if len(gamma_dr) > 0:
            axes_dr[0].hist(gamma_dr, bins=25, range=(0, 5), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_dr0.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_dr[0].set_title(f"Gamma ΔR Distribution (ma {ma_range})", fontsize=18)
    axes_dr[0].set_xlabel(r"$\Delta R(\gamma1, \gamma2)$", fontsize=24)
    axes_dr[0].set_ylabel("A.U.", fontsize=24)
    axes_dr[0].legend(handles=handles_dr0, fontsize=14, handlelength=2, ncol=2)

    # Electron ΔR
    handles_dr1 = []
    for i, (electron_dr, filename) in enumerate(zip(electron_dr_list, file_list)):
        if len(electron_dr) > 0:
            axes_dr[1].hist(electron_dr, bins=25, range=(0, 5), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_dr1.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_dr[1].set_title(f"Electron ΔR Distribution (ma {ma_range})", fontsize=18)
    axes_dr[1].set_xlabel(r"$\Delta R(e1, e2)$", fontsize=24)
    axes_dr[1].set_ylabel("A.U.", fontsize=24)
    axes_dr[1].legend(handles=handles_dr1, fontsize=14, handlelength=2, ncol=2)

    # Muon ΔR
    handles_dr2 = []
    for i, (muon_dr, filename) in enumerate(zip(muon_dr_list, file_list)):
        if len(muon_dr) > 0:
            axes_dr[2].hist(muon_dr, bins=25, range=(0, 5), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_dr2.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # axes_dr[2].set_title(f"Muon ΔR Distribution (ma {ma_range})", fontsize=18)
    axes_dr[2].set_xlabel(r"$\Delta R(\mu1, \mu2)$", fontsize=24)
    axes_dr[2].set_ylabel("A.U.", fontsize=24)
    axes_dr[2].legend(handles=handles_dr2, fontsize=14, handlelength=2, ncol=2)

    # 保存 ΔR 圖
    plt.figure(fig_dr.number)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, dr_output_filename))
    plt.close(fig_dr)

    # 繪製 Gamma ΔR 放大圖 (範圍 0-1)
    fig_gamma_dr_zoom, ax_gamma_dr_zoom = plt.subplots(figsize=(8, 6))
    handles_gamma_dr_zoom = []
    for i, (gamma_dr, filename) in enumerate(zip(gamma_dr_list, file_list)):
        if len(gamma_dr) > 0:
            ax_gamma_dr_zoom.hist(gamma_dr, bins=25, range=(0, 1), histtype='step', density=True, color=colors[i], linewidth=2.0)
            handles_gamma_dr_zoom.append(Line2D([0], [0], color=colors[i], linewidth=2.0, label=filename.replace("ALP_", "").replace(".root", "")))
    # ax_gamma_dr_zoom.set_title(f"Gamma ΔR Distribution (ma {ma_range}, Zoomed)", fontsize=18)
    ax_gamma_dr_zoom.set_xlabel(r"$\Delta R(\gamma1, \gamma2)$", fontsize=24)
    ax_gamma_dr_zoom.set_ylabel("A.U.", fontsize=24)
    ax_gamma_dr_zoom.legend(handles=handles_gamma_dr_zoom, fontsize=14, handlelength=2, ncol=2)
    plt.tight_layout()
    gamma_zoom_output_filename = dr_output_filename.replace(".pdf", "_gamma_zoom.pdf")
    plt.savefig(os.path.join(output_dir, gamma_zoom_output_filename))
    plt.close(fig_gamma_dr_zoom)

# 繪製質量、pT 和 ΔR 圖
plot_mass_and_pt_distributions(
    ma_0p1_0p9_files,
    ma_range="0.1-0.9 GeV",
    mass_output_filename="mass_distributions_ma_0p1_0p9.pdf",
    pt_output_filename="pt_distributions_ma_0p1_0p9.pdf",
    dr_output_filename="dr_distributions_ma_0p1_0p9.pdf",
    alp_mass_range=(0, 1)
)
plot_mass_and_pt_distributions(
    ma_1_30_files,
    ma_range="1-30 GeV",
    mass_output_filename="mass_distributions_ma_1_30.pdf",
    pt_output_filename="pt_distributions_ma_1_30.pdf",
    dr_output_filename="dr_distributions_ma_1_30.pdf",
    alp_mass_range=(0, 35)
)