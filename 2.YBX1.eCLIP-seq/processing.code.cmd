### Initial Job Request #############################################################################################################################################################
dir="eCLIPseq/111519_YC_4CLIP-144151011/"
fastqDIR="fastq"
fastqSuffix=".fastq.gz"
outTO="yeo.lab.out"

starIdxRepeat="RepBase27.07/hum.rep_sub/"
starIdxGenome="ENCODE/GRCh38_no_alt_analysis_set_GCA_000001405.15/"

clipperSpecies="GRCh38"
#####################################################################################################################################################################################

function="YeoLab.pipeline/script/"
cd "$dir"

# 01.extract demux (not needed in this project)

# 02.cutadapt.round1
sh "$function"/02_cutadapt_round1.sh "$outTO"/02.cutadapt.round1/"$sampleID".Tr.fastq.gz "$fastqDIR"/"$filename" "$outTO"/02.cutadapt.round1/"$sampleID".Tr.log

# 03.cutadapt.round2
sh "$function"/03_cutadapt_round2.sh "$outTO"/03.cutadapt.round2/"$sampleID".TrTr.fastq.gz "$outTO"/02.cutadapt.round1/"$sampleID".Tr.fastq.gz "$outTO"/03.cutadapt.round2/"$sampleID".TrTr.log

# 04.fastq sort
sh "$function"/04_fastq_sort.sh "$outTO"/03.cutadapt.round2/"$sampleID".TrTr.fastq.gz "$outTO"/04.fastq.sort/"$sampleID".TrTr.sort.fastq

# 05.mapping to repeat
sh "$function"/05_star_repeat.sh "$starIdxRepeat" "$outTO"/05.star.repeat/"$sampleID".TrTr.sort.repeat-map "$outTO"/04.fastq.sort/"$sampleID".TrTr.sort.fastq.gz

# 06.mapping to genome
sh "$function"/06_star_genome.sh "$starIdxGenome" "$outTO"/06.star.genome/"$sampleID".TrTr.sort.genome-mapped "$outTO"/05.star.repeat/"$sampleID".TrTr.sort.repeat-mapUnmapped.out.mate1

# 07.sort align bam
sh "$function"/07_sort.sh "$outTO"/06.star.genome/"$sampleID".TrTr.sort.genome-mappedAligned.out.bam "$outTO"/07.bam.sort/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.bam "$outTO"/07.bam.sort/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.bam

# 08.remove duplication
sh "$function"/08_umi_tools_dedup.sh "$outTO"/07.bam.sort/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.bam "$outTO"/08.remove.duplication/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost "$outTO"/08.remove.duplication/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.bam

# 09.clipper.peaks
sh "$function"/09_clipper.sh "$clipperSpecies" "$outTO"/08.remove.duplication/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.bam "$outTO"/09.clipper.peaks/"$sampleID".TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.peakClusters.bed

# 10.IP's peak normalization based on input
sh "$function"/10_normalize.sh yeo.lab.out/8/YINPUT1_S7_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.bam yeo.lab.out/10/YINPUT1_S7_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.sort.bam yeo.lab.out/8/YIP1_S5_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.bam yeo.lab.out/10/YIP1_S5_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.sort.bam yeo.lab.out/10/YINPUT1_S7_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.sort.mapped_readnum.txt yeo.lab.out/10/YIP1_S5_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.sort.mapped_readnum.txt yeo.lab.out/9/YIP1_S5_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.peakClusters.bed yeo.lab.out/10/YIP1_S5_L002_R1_001.TrTr.sort.genome-mappedAligned.out.sortByName.thenByLeftmost.rmDup.renameHeader.peakClusters.normed.bed
