#!/bin/bash

for m in {0..40}
do
for z in {5..30}
do
	fn_out=$z"chi""$m.out"
        jb_out=$z"chi""$m.o"
        fn_err=$z"chi""$m.e"
		subCommand="Rscript euler_scripts/cluster_euler_k"$z"_chi""$m.R"
        bsubCommand="bsub -W 1400 -n 1 -J $z"tcga"$m -o $jb_out -e $fn_err -M 64000 -R \"rusage[mem=64000]\" \"$subCommand >> $fn_out\""
		command="echo $subCommand > $fn_out; $bsubCommand"
		echo "$command"
		eval "$command"
done
done

