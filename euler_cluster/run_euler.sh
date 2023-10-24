#!/bin/bash

for m in {1..100}
do
#for z in {1..4}
#do
	fn_out="seed""$m.out"
        jb_out="seed""$m.o"
        fn_err="seed""$m.e"
		subCommand="Rscript euler_scripts/cluster_euler_seed""$m.R"
        bsubCommand="bsub -W 1400 -n 1 -J "tcga"$m -o $jb_out -e $fn_err -M 64000 -R \"rusage[mem=64000]\" \"$subCommand >> $fn_out\""
		command="echo $subCommand > $fn_out; $bsubCommand"
		echo "$command"
		eval "$command"
#done
done

