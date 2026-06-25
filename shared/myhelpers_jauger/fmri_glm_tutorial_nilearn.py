# Title: fmri_glm_tutorial_nilearn.py

# Description:
# Tutorial for performing a general linear model (GLM) analysis on fMRI data using the Nilearn library in Python. 
#

# Created on: June 2026
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital

# %%
from nilearn.datasets import fetch_spm_auditory
subject_data = fetch_spm_auditory()
print("Subject data keys:", subject_data.keys())

# %%
from nilearn.image import mean_img
from nilearn.plotting import plot_anat, plot_img, plot_stat_map, show

fmri_img = subject_data.func
mean_img = mean_img(subject_data.func[0])
# Display mean functional image
plot_img(mean_img, cbar_tick_format="%i")
# Display anatomical image
plot_anat(subject_data.anat, cbar_tick_format="%i")
show()

# %%
# Specify experiment paradigm with timings of the auditory stimulation and rest periods
import pandas as pd
events = pd.read_table(subject_data.events)
events

# %%
# Create a first-level GLM model to fit the experiment paradigm
from nilearn.glm.first_level import FirstLevelModel

fmri_glm = FirstLevelModel(
    t_r=subject_data.t_r,   # repetition time (TR) of acquisition
    noise_model="ar1",      # noise covariance model
    standardize=False,      # do not rescale time series to mean 0, variance 1
    hrf_model="spm",        # use canonical HRF model from SPM
    drift_model="cosine",   # model signal drifts as slow oscillating cosine function
    high_pass=0.01,         # high-pass filter cutoff frequency (Hz)
    verbose=1,
    minimize_memory=False,
)

# %%
# Run the GLM model on the fMRI data
fmri_glm = fmri_glm.fit(fmri_img, events)
fmri_glm


# %%
# Inspect the design matrix (rows = time, columns = predictors)
design_matrix = fmri_glm.design_matrices_[0]

from nilearn.plotting import plot_design_matrix
plot_design_matrix(design_matrix)
show()


# %%
# Save design matrix
from pathlib import Path
output_dir = Path.cwd() / "results" / "plot_single_subject_single_run"
output_dir.mkdir(exist_ok=True, parents=True)
print(f"Output will be saved to: {output_dir}")
plot_design_matrix(design_matrix, output_file=output_dir / "design_matrix.png")


# %%
# Plot first column of design matrix, the expected response profile of auditory stimulation
import matplotlib.pyplot as plt
plt.plot(design_matrix["listening"])
plt.xlabel("scan")
plt.title("Expected Auditory Response")
show()


# %%
# Define contrast to test the effect of stimulation vs rest
import numpy as np
n_regressors = design_matrix.shape[1]
activation = np.zeros(n_regressors)
activation[0] = 1  # Set the first regressor (listening) to be active


# %%
# Plot coefficients of the contrast
from nilearn.plotting import plot_contrast_matrix
plot_contrast_matrix(contrast_def=activation, design_matrix=design_matrix)

# %%
# Compute the estimated effect size in BOLD signal
eff_map = fmri_glm.compute_contrast(activation, output_type="effect_size")
# Compute statistical significance of the contrast, using z-scale
z_map = fmri_glm.compute_contrast(activation, output_type="z_score")

# %%
# Plot thresholded Z-scores
plotting_config = {
    "bg_img": mean_img,
    "display_mode": "z",
    "cut_coords": 3,
    "black_bg": True,
}
plot_stat_map(
    z_map,
    threshold=3,
    title="listening > rest (|Z|>3)",
    figure=plt.figure(figsize=(10, 4)),
    **plotting_config,
)
show()

# %%
# Threshold to control the false positive rate (fpr) at 0.001
from nilearn.glm import threshold_stats_img

clean_map, threshold = threshold_stats_img(
    z_map,
    alpha=0.001,
    height_control="fpr",
    two_sided=False,    # using a one-sided test
)
plotting_config["cmap"] = "inferno"
plot_stat_map(
    clean_map,
    threshold=threshold,
    title=f"listening > rest (Uncorrected p<0.001; threshold={threshold:.3f})",
    figure=plt.figure(figsize=(10, 4)),
    **plotting_config,
)
show()

# %%
# Bonferroni correction to control family wise error rate at 0.05 (very conservative)
clean_map, threshold = threshold_stats_img(
    z_map, 
    alpha=0.05, 
    height_control="bonferroni", 
    two_sided=False
)
plot_stat_map(
    clean_map,
    threshold=threshold,
    title=(
        "listening > rest (p<0.05 Bonferroni-corrected, "
        f"threshold: {threshold:.3f})"
    ),
    figure=plt.figure(figsize=(10, 4)),
    **plotting_config,
)
show()

# %%
# False discovery rate (FDR) correction to control the expected proportion of false positives
# Popular alternative to Bonferroni correction
clean_map, threshold = threshold_stats_img(
    z_map, 
    alpha=0.05, 
    height_control="fdr", 
    two_sided=False
)
plot_stat_map(
    clean_map,
    threshold=threshold,
    title=(
        f"listening > rest (p<0.05 FDR-corrected; threshold: {threshold:.3f})"
    ),
    figure=plt.figure(figsize=(10, 4)),
    **plotting_config,
)
show()

# %%
# Discard isolated voxels (small clusters < 10 voxels) that are likely false positives
clean_map, threshold = threshold_stats_img(
    z_map,
    alpha=0.05,
    height_control="fdr",
    cluster_threshold=10,
    two_sided=False,
)
plot_stat_map(
    clean_map,
    threshold=threshold,
    title=(
        "listening > rest "
        f"(p<0.05 FDR-corrected; threshold: {threshold:.3f}; "
        "clusters > 10 voxels)"
    ),
    figure=plt.figure(figsize=(10, 4)),
    **plotting_config,
)
show()

# %%
# Save effect and zscore maps
z_map.to_filename(output_dir / "listening_gt_rest_z_map.nii.gz")
eff_map.to_filename(output_dir / "listening_gt_rest_eff_map.nii.gz")

# %%
# Extract and report the found positions in a table
from nilearn.reporting import get_clusters_table

table = get_clusters_table(
    z_map, 
    stat_threshold=threshold, 
    cluster_threshold=20
)
table

# %%
table.to_csv(output_dir / "table.csv")

# %%
# F-test of which voxels are well explained by combo of more active or less active than rest
z_map = fmri_glm.compute_contrast(
    activation,
    output_type="z_score",
    stat_type="F",  # set stat_type to 'F' to perform an F test
)

# %%
# Apply same corrections to F-test results
clean_map, threshold = threshold_stats_img(
    z_map,
    alpha=0.05,
    height_control="fdr",
    cluster_threshold=10,
    two_sided=False,
)
plot_stat_map(
    clean_map,
    threshold=threshold,
    title="Effects of interest (fdr=0.05), clusters > 10 voxels",
    figure=plt.figure(figsize=(10, 4)),
    **plotting_config,
)
show()