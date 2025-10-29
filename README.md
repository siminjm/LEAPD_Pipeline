# LEAPD Pipeline (MATLAB)

Unified training, testing, and preprocessing framework for EEG-based classification and correlation using the **LEAPD** method.  
Developed by **Simin Jamshidi**, *University of Iowa*, Departments of **Electrical & Computer Engineering (ECE)** and **Neurology**.

---

## Overview
This repository provides a full **MATLAB workflow** for EEG analysis using the **LEAPD (Linear Predictive Coding-based EEG Analysis for Prognosis and Diagnosis)** approach.

It now includes:
- **EEG Preprocessing module** (EEGLAB + ICLabel) for artifact detection and cleaning  
- **Training and testing modules** for LEAPD-based classification and correlation analyses  

The pipeline supports:
- **Binary classification** (e.g., Deceased vs Living subjects)  
- **Correlation analysis** (e.g., LEAPD indices vs clinical scores such as MoCA or UPDRS)

---

## Core Features

### Preprocessing
- Optional preprocessing using EEGLAB and ICLabel  
- Automatic detection and removal of noisy channels  
- ICA-based artifact removal (blink, muscle, heart)  
- Optional 60 Hz notch filter and harmonics suppression  
- Visualization of raw vs. cleaned EEG signals  
- Automatic saving of cleaned EEG to `/cleaned_data/`

To run the standalone preprocessing demo:
```matlab
cd Preprocessing
Demo
```

---

### LEAPD Pipeline
- Automatic parameter search over filter bands and LPC orders  
- Single- and multi-channel evaluation (1–10 channels)  
- Cross-validation and out-of-sample testing  
- Correlation polarity alignment  
- Comprehensive metrics: ACC, AUC, SEN, SPC, PPV, NPV, OR, LR⁺, ρ, p-value  

---

## Repository Structure
```
LEAPD_Pipeline_MATLAB/
├── Preprocessing/                   # EEG artifact removal module
│   ├── Demo.m
│   ├── pipeline_preprocessing.m
│   ├── detect_noisy_channels.m
│   ├── remove_line_noise.m
│   ├── README.md
│   └── cleaned_data/
│
├── main_train.m                     # LEAPD training script
├── main_test.m                      # LEAPD testing script
├── plot_results.m                   # Visualization
│
├── +utils/                          # Helper functions
│   ├── load_data.m
│   ├── filter_data.m
│   ├── create_filter.m
│   ├── compute_yw.m
│   ├── build_hyperplanes.m
│   ├── compute_leapd_scores.m
│   ├── evaluate_classification.m
│   ├── evaluate_correlation.m
│   ├── combine_scores.m
│   ├── generate_combinations.m
│   ├── pick_polarity_and_rho.m
│   ├── read_labels_table.m
│   ├── fetch_targets.m
│   ├── count_subjects.m
│   └── save_results.m
│
├── results/
│   ├── train_results/
│   └── test_results/
│
├── README.md
├── LICENSE
└── .gitignore
```

---

## Quick Start

### 1. (Optional) Preprocessing
Run the preprocessing demo to remove noise and artifacts:
```matlab
cd Preprocessing
Demo
```

You can also call it manually:
```matlab
[X_clean, labels_clean, report, savePath] = pipeline_preprocessing(X_raw, Fs, labels);
```

---

### 2. Training (LEAPD)
```matlab
cfg = struct;
cfg.mode        = "correlation";              % or "classification"
cfg.data_train  = "data/EEG_train.mat";
cfg.labels_file = "data/ClinicalLabels.xlsx"; % columns: ID, Target
cfg.save_dir    = "results/train_results";

results_train = main_train(cfg);
```

---

### 3. Testing (Out-of-Sample)
```matlab
cfg2 = struct;
cfg2.mode            = "correlation";                 
cfg2.data_test       = "data/EEG_test.mat";
cfg2.trained_model   = "results/train_results/BestParamsAll.mat";
cfg2.labels_file     = "data/ClinicalLabels_Test.xlsx";
cfg2.combo_sizes     = 1:10;
cfg2.max_full_combos = 5;
cfg2.save_dir        = "results/test_results";

results_test = main_test(cfg2);
```

---

## Data Privacy Notice
Due to clinical confidentiality agreements, the original EEG datasets used in this research cannot be shared publicly.  
However, the **entire analysis pipeline**, algorithms, and reproducible code structure are fully provided to ensure transparency.  
All stages—from preprocessing to LEAPD evaluation—can be executed using **synthetic or anonymized EEG data**.

---

## Skills Demonstrated
- **MATLAB** — advanced EEG signal processing and algorithm design  
- **EEG Preprocessing** — ICA, ICLabel, and noise/artifact removal  
- **Feature Extraction** — LPC coefficients and hyperplane distances  
- **Statistical Analysis** — correlation, classification metrics, AUC, OR, LR⁺  
- **Reproducible Research** — modular, documented, version-controlled workflow  

---

## Authors
**Simin Jamshidi** — Ph.D. Candidate  
Departments of **Electrical & Computer Engineering (ECE)** and **Neurology**  
**University of Iowa**

Supervisors:  
**Prof. Soura Dasgupta** and **Dr. Nandakumar Narayanan**

---

## Citation
If you use this pipeline in academic or research work, please cite:

> Jamshidi, S., et al. (Year).  
> *EEG-Based Mortality and Cognitive Decline Prediction in Parkinson’s Disease using the LEAPD Method.*  
> Departments of Electrical & Computer Engineering (ECE) and Neurology, University of Iowa.  
> (Preprint or Journal Reference — to be updated)

---

## License
Released under the **MIT License (with Citation Request)**.  
See the [LICENSE](./LICENSE) file for details.
