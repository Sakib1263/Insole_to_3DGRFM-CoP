import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.signal import resample, medfilt, find_peaks
from scipy.signal import butter, lfilter, freqz

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 18,
    "text.color": "#212121",
    "axes.edgecolor": "#212121",
    "xtick.color": "#212121",
    "ytick.color": "#212121",
    "axes.labelcolor": "#212121",
    'legend.frameon': False,
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ROOT_DIR = "../../"
DATA_DIR = os.path.join(ROOT_DIR, "dataset", "dataset3")
PROCESSED_DATA_DIR = os.path.join(ROOT_DIR, "dataset", "processed")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")
FORCE_DIR = "Force"
IMU_DIR = "IMU"

SUBJECTS = ["Subj04", "Subj05", "Subj06", "Subj07", "Subj08", "Subj09",
            "Subj10", "Subj11", "Subj12"]

WEIGHTS = [75.2, 65.1, 56.0, 80.8, 61.5, 81.0, 61.6, 62.4, 69.0]

TRIALS = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54"]
        #   "run_63", "run_81", "run_99", "jump", "squat", "land", "lunge"]

# TRIALS = ["walk_27", "walk_36", "walk_45", "walk_54"]

FORCE_COLUMNS = ["RFx", "RFy", "RFz", "RMx", "RMy", "RMz",
                 "LFx", "LFy", "LFz", "LMx", "LMy", "LMz"]

IMU_COLUMNS = ["HAx", "HAy", "HAz", "RAx", "RAy", "RAz",
               "LAx", "LAy", "LAz"]

SEGMENT_LEN = 1024
FORCE_THRESHOLD = 20

FS_FORCE = 2048
FS_IMU = 240
OFFSET = 200

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def normalize(x: np.ndarray) -> np.ndarray:
    return x
    # return (x - np.min(x)) / (np.max(x) - np.min(x))


def apply_low_pass_filter(
    x: np.ndarray,
    order: int = 4,
    cutoff: int = 15,
    fs: int = FS_FORCE
) -> np.ndarray:
    b, a = butter(order, cutoff, fs=fs, btype='low', analog=False)
    return lfilter(b, a, x)


def preprocess_data(x: np.ndarray, fs: int) -> np.ndarray:
    return normalize(apply_low_pass_filter(resample(x, SEGMENT_LEN), fs=fs))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plt.figure(figsize=(8, 6))
ax = plt.gca()


for trial in TRIALS:
    for subject, weight in zip(SUBJECTS, WEIGHTS):
        fx_segments = []
        fy_segments = []
        fz_segments = []

        ax_segments = []
        ay_segments = []
        az_segments = []

        file_name = subject + "_" + trial + ".csv"
        force_path = os.path.join(DATA_DIR, FORCE_DIR, file_name)
        imu_path = os.path.join(DATA_DIR, IMU_DIR, file_name)

        force = pd.read_csv(force_path, names=FORCE_COLUMNS)
        imu = pd.read_csv(imu_path, names=IMU_COLUMNS)

        delta = imu.shape[0] / force.shape[0]

        rfx = force["RFx"].to_numpy()
        rfy = force["RFy"].to_numpy()
        rfz = force["RFz"].to_numpy()

        rmx = force["RMx"].to_numpy()
        rmy = force["RMy"].to_numpy()
        rmz = force["RMz"].to_numpy()

        rax = imu["LAx"].to_numpy()
        ray = imu["LAy"].to_numpy()
        raz = imu["LAz"].to_numpy()

        ram = np.sqrt(rax ** 2 + ray ** 2 + raz ** 2)

        mask = (rfy > FORCE_THRESHOLD).astype("int")
        _, config = find_peaks(mask, width=500)

        for (start, end) in zip(config["left_bases"], config["right_bases"]):
            start = start + OFFSET
            end = end + OFFSET
            fx_segments.append(preprocess_data(
                rfx[start: end] / weight, FS_FORCE))
            fy_segments.append(preprocess_data(
                rfy[start: end] / weight, FS_FORCE))
            fz_segments.append(preprocess_data(
                rfz[start: end] / weight, FS_FORCE))

            start = int(np.round(start * delta))
            end = int(np.round(end * delta))

            ax_segments.append(preprocess_data(rax[start: end], FS_IMU))
            ay_segments.append(preprocess_data(ray[start: end], FS_IMU))
            az_segments.append(preprocess_data(raz[start: end], FS_IMU))
        # break

        segments = np.array(fy_segments)

    try:
        mean = np.mean(segments, axis=0)
        std = np.std(segments, axis=0)
        t = np.arange(0, len(mean)) / SEGMENT_LEN * 100
        ax.plot(t, mean, label=trial)
        ax.fill_between(t, mean + std, mean - std, alpha=0.3)
        ax.set_xlabel("Gait Cycle (\%)")
        ax.set_ylabel("Force ($N/kg$)")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlim([0, 150])
        plt.title(f"$F_y$")
        plt.tight_layout()

    except Exception as e:
        print(f"Not enoogh samples for {trial}")

    # break

plt.legend()
plt.savefig(os.path.join(FIGURES_DIR, "Fy_all_trials.svg"))
plt.savefig(os.path.join(FIGURES_DIR, "Fy_all_trials.png"))
plt.show()
