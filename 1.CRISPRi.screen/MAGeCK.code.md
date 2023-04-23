library="lncRNA_CRISPRi_sgRNA_library_MAGeCK.txt"
ctrl="lncRNA_CRISPRi_sgRNA_CTRL_MAGeCK.txt"
Output="U87" # (U87/U251)
OutPackage="MAGeCK_result_analysis_secondbest"

mageck count -l /rsrch3/scratch/bcb/ywei4/sgRNA/MAGeCK-library/"$library" --control-sgrna /rsrch3/scratch/bcb/ywei4/sgRNA/MAGeCK-library/"$ctrl" --norm-method control -n "$Output" --sample-label ${labels#*,} --fastq ${files#* }

mageck test -k "$Output".count.txt -t ${groupD21#*,} -c ${groupD0#*,} --control-sgrna /rsrch3/scratch/bcb/ywei4/sgRNA/MAGeCK-library/"$ctrl" --norm-method control --keep-tmp -n "$Output"_D21_D0 --gene-lfc-method secondbest

