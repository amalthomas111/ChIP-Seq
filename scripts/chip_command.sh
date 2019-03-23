snakemake -s chip_preprocess_fromfastq.Snakemake --cluster-config chip_cluster.yaml --jobs 50  --jobname 'chip.{rulename}.{jobid}.sh' --restart-times 22 --cluster "sbatch -p {cluster.p} --time={cluster.time} --mem-per-cpu={cluster.mem} --ntasks={cluster.ntasks}" 