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
    "legend.frameon": False,
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ROOT_DIR = "../../"
DATA_DIR = os.path.join(ROOT_DIR, "dataset", "dataset3")
PROCESSED_DATA_DIR = os.path.join(ROOT_DIR, "dataset", "processed")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")
INSOLE_DIR = "Insole"

SUBJECTS = ["Subj04", "Subj05", "Subj06", "Subj07", "Subj08", "Subj09",
            "Subj10", "Subj11", "Subj12"]

WEIGHTS = [75.2, 65.1, 56.0, 80.8, 61.5, 81.0, 61.6, 62.4, 69.0]

TRIALS = ["walk_09", "walk_18", "walk_27", "walk_36", "walk_45", "walk_54"]
        #   "run_63", "run_81", "run_99", "jump", "squat", "land", "lunge"]

# TRIALS = ["walk_27", "walk_36", "walk_45", "walk_54"]

INSOLE_COLUMNS = [
    "LP0", "LP1", "LP2", "LP3", "LP4", "LP5", "LP6", "LP7", "LP8", "LP9",
    "LP10", "LP11","LP12", "LP13", "LP14", "LP15", "LAX", "LAY", "LAZ",
    "LGX", "LGY", "LGZ", "LTF", "LCX", "LCY",
    "RP0", "RP1", "RP2", "RP3", "RP4", "RP5", "RP6", "RP7", "RP8", "RP9",
    "RP10", "RP11","RP12", "RP13", "RP14", "RP15", "RAX", "RAY", "RAZ",
    "RGX", "RGY", "RGZ", "RTF", "RCX", "RCY"
]

SEGMENT_LEN = 1024
FORCE_THRESHOLD = 20
FS_INSOLE = 100

OFFSET = 0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def normalize(x: np.ndarray) -> np.ndarray:
    return x
    # return (x - np.min(x)) / (np.max(x) - np.min(x))


def apply_low_pass_filter(
    x: np.ndarray,
    order: int = 4,
    cutoff: int = 15,
    fs: int = FS_INSOLE
) -> np.ndarray:
    b, a = butter(order, cutoff, fs=fs, btype="low", analog=False)
    return lfilter(b, a, x)


def preprocess_data(x: np.ndarray, fs: int) -> np.ndarray:
    return normalize(apply_low_pass_filter(resample(x, SEGMENT_LEN), fs=fs))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plt.figure(figsize=(8, 6))
ax = plt.gca()


for subject, weight in zip(SUBJECTS, WEIGHTS):
    for trial in TRIALS:
        f_segments = []

        ax_segments = []
        ay_segments = []
        az_segments = []

        file_name = subject + "_" + trial + ".csv"
        insole_path = os.path.join(DATA_DIR, INSOLE_DIR, file_name)

        insole = pd.read_csv(insole_path, names=INSOLE_COLUMNS)

        rf = insole["RTF"].to_numpy()
        lf = insole["LTF"].to_numpy()

        # plt.plot(rf)
        # plt.show()

        # break

        rax = insole["RAX"].to_numpy()
        ray = insole["RAY"].to_numpy()
        raz = insole["RAZ"].to_numpy()

        ram = np.sqrt(rax ** 2 + ray ** 2 + raz ** 2)

        mask = (lf > FORCE_THRESHOLD).astype("int")
        _, config = find_peaks(mask, width=[50, 300])

        for (start, end) in zip(config["left_bases"], config["right_bases"]):
            start = start + OFFSET
            end = end + OFFSET
            f_segments.append(preprocess_data(
                lf[start: end] / weight, FS_INSOLE))

            ax_segments.append(preprocess_data(rax[start: end], FS_INSOLE))
            ay_segments.append(preprocess_data(ray[start: end], FS_INSOLE))
            az_segments.append(preprocess_data(raz[start: end], FS_INSOLE))
        # break

        segments = np.array(az_segments)

    try:
        # for i in range(segments.shape[0]):
        #     ax.plot(segments[i, :].ravel())
        mean = np.mean(segments, axis=0)
        std = np.std(segments, axis=0)
        t = np.arange(0, len(mean)) / SEGMENT_LEN * 100
        ax.plot(t, mean, label=subject)
        ax.fill_between(t, mean + std, mean - std, alpha=0.3)
        ax.set_xlabel("Gait Cycle (\%)")
        ax.set_ylabel(r"Acceleration ($ms^-2$)")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xlim([0, 150])
        plt.title(r"$A_z$")
        plt.tight_layout()

    except Exception as e:
        print(f"Not enoogh samples for {trial}")

        # break
    # break

plt.legend()
plt.savefig(os.path.join(FIGURES_DIR, "Az_insole_all_subjects.svg"))
plt.savefig(os.path.join(FIGURES_DIR, "Az_insole_all_subjects.png"))
plt.show()
