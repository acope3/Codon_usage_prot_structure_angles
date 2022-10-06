#!/bin/bash

# Edit these to suit your needs
PROCESSES=48
DATASET_DIR="dataverse_files_sim_w_ena/"
TAG="simulated"

# Values used in the paper results
MIN_GROUP=1
KDE_NBINS=128
KDE_WIDTH=200
DDIST_BS_NITER=25
DDIST_K=200
DDIST_K_MIN=100
DDIST_K_TH=50
DDIST_NMAX=200
DDIST_STATISTIC="kde_g"
DDIST_KERNEL_SIZE=2.0
FDR=0.05

set -eux
pp5 -p="$PROCESSES" \
 analyze-pointwise \
 --dataset-dir="$DATASET_DIR" \
 --min-group-size="$MIN_GROUP" \
 --kde-width="$KDE_WIDTH" \
 --kde-nbins="$KDE_NBINS" \
 --ddist-statistic="$DDIST_STATISTIC" \
 --ddist-k="$DDIST_K" \
 --ddist-k-min="$DDIST_K_MIN" \
 --ddist-k-th="$DDIST_K_TH" \
 --ddist-bs-niter="$DDIST_BS_NITER" \
 --ddist-n-max="$DDIST_NMAX" \
 --ddist-kernel-size="$DDIST_KERNEL_SIZE" \
 --fdr="$FDR" \
 --comparison-types aa cc \
 --ignore-omega \
 --out-tag="$TAG"
