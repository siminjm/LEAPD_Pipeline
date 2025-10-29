# EEG Preprocessing Pipeline (MATLAB + EEGLAB)

This repository contains a **fully automated EEG preprocessing pipeline** built in MATLAB, designed for both synthetic and real EEG datasets.  
It demonstrates how to clean EEG data by removing line noise, detecting noisy channels, and performing **ICA-based artifact rejection** using the [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) toolbox and **ICLabel** plugin.

---

## Features

- **Automatic line-noise removal** (60 Hz notch filter, configurable)
- **Noisy-channel detection** using standard deviation
- **ICA + ICLabel** for automatic artifact component rejection  
  (blink, muscle, heart, and line-noise ICs)
- **Automatic EEG reconstruction** after IC removal
- **Before/after visualization** with aligned time axes and power spectra
- **Automatic saving** of cleaned EEG (`.mat` file)

---

## Folder Structure

Preprocessing_MATLAB/
│
├── eeglab/ # EEGLAB toolbox folder (do not rename subfolders)
│
├── Demo.m # Main demonstration script
├── pipeline_preprocessing.m # Core preprocessing pipeline
├── detect_noisy_channels.m # Identifies bad EEG channels
├── remove_line_noise.m # 60 Hz notch filter (configurable)
│
└── cleaned_data/ # Automatically generated cleaned EEG files


---

## Configuration Options

The pipeline can be customized via a simple configuration structure inside `pipeline_preprocessing.m`.

```matlab
% Example configuration
cfg.notch60Hz = true;    % Enable (true) or disable (false) 60 Hz notch filter
cfg.plotResults = true;  % Show before/after plots
cfg.saveCleaned = true;  % Automatically save cleaned EEG to /cleaned_data