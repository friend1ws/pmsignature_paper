# Script for generating figures in SText2.
# We need several results ("../../analysis/UCUT/result"),
# that are generated the script (""../../analysis/UCUT/script/runall.sh").

if (!file.exists("../result/SText2")) {
  dir.create("../result/SText2")
}

# You may need to modify the size of device...
source("subscript_SText2/UCUT_variousK.R")

# generating statistics for the estimation results such as
# likelihood, bootstrap errors and correlations.
source("subscript_SText2/UCUT_stat.R")