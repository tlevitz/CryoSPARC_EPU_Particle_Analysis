# CryoSPARC_EPU_Particle_Analysis

This script (or Jupyter notebook) will analyze particles from a given CryoSPARC job and its associated EPU-collected micrographs. The default EPU and CryoSPARC directory styles with metadata are expected. 

# Getting Ready
1. The script requires a python3.9 environment with numpy and matplotlib
  ```bash
  conda create --name particle_analysis python=3.9 numpy matplotlib
  ```
2. You also must have cs2star available and callable as cs2star: https://github.com/brisvag/cs2star

# Running
```bash
python particle_distribution.py {job} {rawdatapath (optional)}")
```

Examples (running from CryoSPARC project directory):
```bash
python particle_distribution.py J88
```
```bash
python particle_distribution.py J88 /path/to/raw/data/symlinks
```

_The script will attempt to find the raw data automatically based on the CryoSPARC (Live or traditional) input, but the directory can be given directly as input if needed_

# Outputs

The script (or notebook) will put all outputs in a directory named particle_stats_JXX.

_"Non-empty micrographs" refer to micrographs that contributed zero particles to the input CryoSPARC job_

## Average transmission per gridsquare, for all micrographs and for all non-empty micrographs
<img width="800" height="600" alt="transmission_vs_gridsquare_allmics" src="https://github.com/user-attachments/assets/1d1a0283-861d-4d1c-a022-2db7b4f5cdd3" />

## Total particles aggregated per grid square
<img width="1000" height="600" alt="total_particles" src="https://github.com/user-attachments/assets/2d54cfb8-b182-4e9d-9ed5-871044f6d4f3" />

## Percentage of micrographs per grid square with zero particles
<img width="640" height="480" alt="percent_empty" src="https://github.com/user-attachments/assets/18d4ed88-8c5c-41f1-b79f-e3d447b03830" />

## Scatterplot of number of particles vs percent transmission of parent micrograph, color-coded by grid square
<img width="1000" height="600" alt="particles_vs_transmission" src="https://github.com/user-attachments/assets/e9795fa8-79cc-4289-9ecc-96a3960dd265" />

## Average number of particles per micrograph at each defocus value, for all micrographs and for all non-empty micrographs
<img width="800" height="600" alt="particles_vs_defocus_allmics" src="https://github.com/user-attachments/assets/a91b1ffa-5ae2-4c2e-99a0-6b9861fe93eb" />

## Number of particles per micrograph, for all micrographs and for all non-empty micrographs, color-coded by grid square
<img width="1200" height="600" alt="particles_per_mic_allmics" src="https://github.com/user-attachments/assets/5ce19f6f-0228-4a7a-88cd-f4c72c209e06" />

## Distribution of particles per micrograph across all micrographs
<img width="1000" height="600" alt="distribution_particles_per_mic" src="https://github.com/user-attachments/assets/c465cbac-c5a3-4fd3-96c8-89824c0374fb" />

## Average particles per micrograph for each grid square, for all micrographs and for all non-empty micrographs
<img width="640" height="480" alt="avg_particles_allmics" src="https://github.com/user-attachments/assets/9e13f5c2-ccd6-451d-a7e6-c00f5ed42d85" />

There will also be one output subdirectory (in addition to script housekeeping) named particles_vs_transmission_gridsquares. This directory contains 
## Individual scatterplots of number of particles per micrograph vs transmission for each gridsquare. 
<img width="1000" height="600" alt="particles_vs_transmission_sq1" src="https://github.com/user-attachments/assets/850c8deb-a689-48f3-8297-4fa5cd380b48" />

This script was generated with the assistance of GPT4DFCI, a private, HIPAA-secure endpoint to GPT-4o provided by Dana-Farber Cancer Institute
