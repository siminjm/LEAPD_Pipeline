## LEAPD Pipeline (MATLAB)

Unified training and testing framework for EEG-based classification and correlation using the LEAPD method.  
Developed by Simin Jamshidi (University of Iowa).

---

## Overview
This repository implements a full MATLAB workflow for LEAPD-based analysis of EEG data.  
It supports both:

- Binary classification (e.g., Deceased vs Living subjects)  
- Correlation analysis (e.g., LEAPD indices vs clinical scores such as MoCA or UPDRS)

Core features:
- Automatic parameter search over filter bands and LPC orders  
- Single- and multi-channel evaluation (1–10 channels)  
- Cross-validation and out-of-sample testing  
- Correlation polarity alignment  
- Comprehensive performance metrics (ACC, AUC, SEN, SPC, PPV, NPV, OR, LR⁺, ρ, p-value)

---

## Repository Structure
LEAPD_Pipeline/
├── main_train.m
├── main_test.m
├── +utils/
│ ├── load_data.m
│ ├── filter_data.m
│ ├── create_filter.m
│ ├── compute_yw.m
│ ├── build_hyperplanes.m
│ ├── compute_leapd_scores.m
│ ├── evaluate_classification.m
│ ├── evaluate_correlation.m
│ ├── combine_scores.m
│ ├── generate_combinations.m
│ ├── pick_polarity_and_rho.m
│ ├── read_labels_table.m
│ ├── fetch_targets.m
│ ├── count_subjects.m
│ └── save_results.m
├── results/
│ ├── train_results/
│ └── test_results/
├── README.md
├── LICENSE
└── .gitignore


---

## Quick Start
```matlab
%% === Training ===
cfg = struct;
cfg.mode        = "correlation";              % or "classification"
cfg.data_train  = "data/EEG_train.mat";
cfg.labels_file = "data/ClinicalLabels.xlsx"; % columns: ID, Target
cfg.save_dir    = "results/train_results";
results_train   = main_train(cfg);

%% === Testing ===
cfg2 = struct;
cfg2.mode          = "correlation";           % or "classification"
cfg2.data_test     = "data/EEG_test.mat";
cfg2.trained_model = "results/train_results/BestParamsAll.mat";
cfg2.combo_sizes   = 1:10; 
cfg2.max_full_combos = 5;
results_test = main_test(cfg2);


---

## Data Privacy Notice
Due to clinical confidentiality agreements, the original EEG datasets used for this research cannot be shared publicly.  
However, the **entire analysis pipeline**, algorithms, and reproducible code structure are provided here to ensure transparency and demonstrate the complete workflow.  
All processing, feature extraction, and evaluation steps can be executed using synthetic or anonymized EEG data.

---

## Skills Demonstrated
- **MATLAB** — advanced signal processing, data structures, and parallelization  
- **EEG Preprocessing** — filtering, normalization, and artifact exclusion  
- **Feature Extraction** — Linear Predictive Coding (LPC), LEAPD indices, and hyperplane distance computation  
- **Statistical Analysis** — correlation (Spearman), classification metrics (ACC, AUC, SEN, SPC, PPV, NPV, OR, LR⁺)  
- **Algorithm Design** — parameter search grids, LOOCV, and progressive multi-channel combination logic  
- **Reproducible Research** — modular pipeline design with documentation, version control (Git), and open-source licensing  

---

