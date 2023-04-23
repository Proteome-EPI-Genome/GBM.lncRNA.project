### Initial Job Request ######################################################################################################
dir="RNAseq/U251/190228_YC_PM968" # (DARS-AS1: RNAseq/U251/190228_YC_PM968 or YBX1: RNAseq/U251/190613_YC531-1_PM1027)
FQsubDIR="fastq"

index="hg38/GRCh38.v25.primary_assembly.genome"
annotation="hg38/gencode.v25.primary_assembly.annotation.gff3"

package="star_result_hg38v25"
CountOutFile="U251.190228_YC_PM968.RNA-seq.gc25.htseqcount"
#############################################################################################################################
header="Gid     Symbol  Type"

#Step 1, align fastq reads to refgenome using STAR
~/Programs/STAR-2.6.1b/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ~/star-index/"$index"/ --readFilesCommand zcat --readFilesIn "$dir"fastq/"${filename%.fastq*}".fastq.gz --outFileNamePrefix "$dir$package"/bam/"${filename%.fastq*}" --outSAMunmapped Within --outFilterType BySJout --twopassMode Basic --outSAMtype BAM SortedByCoordinate

#Step 2, quantifying Genes
echo "$header" > "$dir$package"/htseq-count/"$CountOutFile".reverse.out
python -m HTSeq.scripts.count --stranded reverse --additional-attr gene_name gene_type -f bam $bamfiles ~/refgenomes/"$annotation" >> "$dir$package"/htseq-count/"$CountOutFile".reverse.out

