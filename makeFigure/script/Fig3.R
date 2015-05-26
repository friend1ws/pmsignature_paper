# Script for generating figures in SText 1 (simulation study)
# First, we need the results ("../../analysis/UCUT/result"),

if (!file.exists("../result/Fig3")) {
  dir.create("../result/Fig3")
}


# script for generating (full and independent) mutation signatures obtained in UCUT analysis.
source("subscript_Fig3/UCUT_signature.R")

# script for generating (full and independent) the result of down-sampling analysis.
source("subscript_Fig3/UCUT_downsampling.R")