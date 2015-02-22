
# Hoang et al. Science Trandlational Medicine, 2013
# 3006200s.txt: downloaded from the science web site.

echo "perl subscript/proc.hoang.pl ../../data/UTUC/raw/3006200s.txt > ../../data/UTUC/raw/mutInfo.hoang.txt"
perl subscript/proc.hoang.pl ../../data/UTUC/raw/3006200s.txt > ../../data/UTUC/raw/mutInfo.hoang.txt

echo "sh ../lib/procMutForPMS/batch_proc_mutInfo.sh ../../data/UTUC/raw/mutInfo.hoang.txt ../../data/UTUC/input hoang"
sh ../lib/procMutForPMS/batch_proc_mutInfo.sh ../../data/UTUC/raw/mutInfo.hoang.txt ../../data/UTUC/input hoang


