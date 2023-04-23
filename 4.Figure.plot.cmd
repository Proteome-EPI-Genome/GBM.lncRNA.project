setwd("GBM.lncRNA.project")

## Figure 1C & 1D. Volcanoplot for screen (gene-level) ======================================
GenelevelData = lncRNA_CRISPRi_MAGeCK_sgrna_summary2genelevel_secondbest.1182.20220828
GenelevelData$Symbol = substr(GenelevelData$Gene, unlist(gregexpr(pattern ='\\|', GenelevelData$Gene))+1, nchar(GenelevelData$Gene))

GenelevelData = as.data.frame(GenelevelData)

SGnumCutoff = 2
LFC = log2(1.5)
pvalue = 0.01

GenelevelData$U251Group = "non-significant"
GenelevelData$U251Group[which(
    GenelevelData$U251_down_fc2 <= -LFC & 
    GenelevelData$U251_down_pvalue2 <= pvalue
)] = "significant"
GenelevelData$U251Group[which(
  GenelevelData$U251_down <= -SGnumCutoff & 
    GenelevelData$U251_down_fc2 <= -LFC & 
    GenelevelData$U251_down_pvalue2 <= pvalue
)] = ">= 2 sgRNAs"
GenelevelData$U251Group[which(
  (
    GenelevelData$U251_down <= -SGnumCutoff & 
      GenelevelData$U251_down_fc2 <= -LFC & 
      GenelevelData$U251_down_pvalue2 <= pvalue
  ) & 
    (
      GenelevelData$U87_down <= -SGnumCutoff & 
        GenelevelData$U87_down_fc2 <= -LFC & 
        GenelevelData$U87_down_pvalue2 <= pvalue
    )
)] = "common"

GenelevelData$U87Group = "non-significant"
GenelevelData$U87Group[which(
  GenelevelData$U87_down_fc2 <= -LFC & 
    GenelevelData$U87_down_pvalue2 <= pvalue
)] = "significant"
GenelevelData$U87Group[which(
  GenelevelData$U87_down <= -SGnumCutoff & 
    GenelevelData$U87_down_fc2 <= -LFC & 
    GenelevelData$U87_down_pvalue2 <= pvalue
)] = ">= 2 sgRNAs"
GenelevelData$U87Group[which(
  (
    GenelevelData$U251_down <= -SGnumCutoff & 
      GenelevelData$U251_down_fc2 <= -LFC & 
      GenelevelData$U251_down_pvalue2 <= pvalue
  ) & 
    (
      GenelevelData$U87_down <= -SGnumCutoff & 
        GenelevelData$U87_down_fc2 <= -LFC & 
        GenelevelData$U87_down_pvalue2 <= pvalue
    )
)] = "common"

table(GenelevelData$U251Group)
table(GenelevelData$U87Group)

# numeric
GenelevelData$U251_down_fc2 = as.numeric(unlist(GenelevelData$U251_down_fc2))
GenelevelData$U251_down_pvalue2 = as.numeric(unlist(GenelevelData$U251_down_pvalue2))
GenelevelData$U87_down_fc2 = as.numeric(unlist(GenelevelData$U87_down_fc2))
GenelevelData$U87_down_pvalue2 = as.numeric(unlist(GenelevelData$U87_down_pvalue2))
# log10 transformation of pvalue
GenelevelData$U251_down_pvalue2 = -log10(GenelevelData$U251_down_pvalue2)
GenelevelData$U87_down_pvalue2 = -log10(GenelevelData$U87_down_pvalue2)

Selected_lncRNAs <- read_excel("Selected-lncRNAs.xlsx", sheet = "Sheet1", col_names = FALSE)

pdf("1.genelevel.VolcanoPlot.1182.lncRNA_CRISPRi.pdf", width=8.5, height=8.5)
ggscatter(GenelevelData, 
          x = "U251_down_fc2", 
          y = "U251_down_pvalue2", 
          ylab = "-log10(U251_down_pvalue2)",
          xlab = "log2FC",
          main = "U251",
          #xlim = c(-2, 0),
          #ylim = c(0, 15),
          color = "U251Group",
          size = ifelse(GenelevelData$Gene %in% Selected_lncRNAs$...2, 3, 1.5),
          #shape = ifelse(GenelevelData$Gene %in% Selected_lncRNAs$...2, 8, 19),
          label = "Symbol", 
          repel = T,
          palette = c("#FFB000", "#DC267F", "#999999", "#648FFF"),
          label.select = GenelevelData$Symbol[which(GenelevelData$Gene %in% Selected_lncRNAs$...2)],
          font.label = c(14, "bold")
) +
  geom_hline(yintercept = -log10(pvalue), colour="gray", linetype="dashed") +
  geom_vline(xintercept = -LFC, colour="gray", linetype="dashed")

ggscatter(GenelevelData, 
          x = "U87_down_fc2", 
          y = "U87_down_pvalue2", 
          ylab = "-log10(U87_down_pvalue2)",
          xlab = "log2FC",
          main = "U87",
          #xlim = c(-2, 0),
          #ylim = c(0, 15),
          color = "U87Group",
          size = ifelse(GenelevelData$Gene %in% Selected_lncRNAs$...2, 3, 1.5),
          #shape = ifelse(GenelevelData$Gene %in% Selected_lncRNAs$...2, 8, 19),
          label = "Symbol", 
          repel = T,
          palette = c("#FFB000", "#DC267F", "#999999", "#648FFF"),
          label.select = GenelevelData$Symbol[which(GenelevelData$Gene %in% Selected_lncRNAs$...2)],
          font.label = c(14, "bold")
) +
  geom_hline(yintercept = -log10(pvalue), colour="gray", linetype="dashed") +
  geom_vline(xintercept = -LFC, colour="gray", linetype="dashed")
dev.off()

## Figure 2E & 4E. Volcanoplot for screen (gene-level) ======================================
# Fig4Publication/TCGA-GBM.GTEx-Brain.exp.deseq2-normalized.log2.tsv

# DARS-AS1 (ENSG00000231890)
# LINC01063 (ENSG00000232065)
# RPPH1 (ENSG00000277209)
# ENSG00000256128.4|LINC00944 (Version 2: 2022.09.01 Added)
# ENSG00000065978.18|YBX1 (Version 2: 2022.09.01 Added)

DARSAS14box <- data.frame(
  DARSAS1 = as.numeric(TCGA.GBM.GTEx.Brain.exp[grep("ENSG00000231890", row.names(TCGA.GBM.GTEx.Brain.exp)), c(GTEx, GBM.Normal, GBM.Tumor)]),
  LINC01063 = as.numeric(TCGA.GBM.GTEx.Brain.exp[grep("ENSG00000232065", row.names(TCGA.GBM.GTEx.Brain.exp)), c(GTEx, GBM.Normal, GBM.Tumor)]),
  RPPH1 = as.numeric(TCGA.GBM.GTEx.Brain.exp[grep("ENSG00000277209", row.names(TCGA.GBM.GTEx.Brain.exp)), c(GTEx, GBM.Normal, GBM.Tumor)]),
  LINC00944 = as.numeric(TCGA.GBM.GTEx.Brain.exp[grep("ENSG00000256128", row.names(TCGA.GBM.GTEx.Brain.exp)), c(GTEx, GBM.Normal, GBM.Tumor)]),
  YBX1 = as.numeric(TCGA.GBM.GTEx.Brain.exp[grep("ENSG00000065978", row.names(TCGA.GBM.GTEx.Brain.exp)), c(GTEx, GBM.Normal, GBM.Tumor)]),
  SampleType = c(rep('Normal',(length(GTEx) + length(GBM.Normal))),rep('Tumor',length(GBM.Tumor))),
  SampleID = colnames(TCGA.GBM.GTEx.Brain.exp)[c(GTEx, GBM.Normal, GBM.Tumor)]
)

pdf("2.GTEx-Brain.TCGA-GBM.DARS-AS1.YBX1.seperate.pdf", width=4.5, height=5.5)
ggboxplot(DARSAS14box, x = "SampleType",
          y = c("DARSAS1"),
          title = "DARS-AS1",
          color = "SampleType", palette = c("#0066ff","#ff3399"),
          ylab = "Expression\nlog2NormalizedCount",
          add = "dotplot",                              # Add dotplot
          add.params = list(binwidth = 0.1, dotsize = 0.2)
) + stat_compare_means(comparisons = my_comparisons)
ggboxplot(DARSAS14box, x = "SampleType",
          y = c("YBX1"),
          title = "YBX1",
          color = "SampleType", palette = c("#0066ff","#ff3399"),
          ylab = "Expression\nlog2NormalizedCount",
          add = "dotplot",                              # Add dotplot
          add.params = list(binwidth = 0.1, dotsize = 0.2)
) + stat_compare_means(comparisons = my_comparisons)
dev.off()

## Figure 2D. GSEA plot ======================================
gsea.out = "RNA-seq/U251/YBX1_genelist/GSEA_results/DARSAS1_All.MSigDB.v6.2.gene.sets.NOFILTER.Gsea.1570023107777"

replotGSEA(path = gsea.out, gene.set = "BASAKI_YBX1_TARGETS_UP", class.name = "DARS-AS1_SI / NC")
replotGSEA(path = gsea.out, gene.set = "FUJII_YBX1_TARGETS_DN", class.name = "DARS-AS1_SI / NC")
replotGSEA(path = gsea.out, gene.set = "HALLMARK_E2F_TARGETS", class.name = "DARS-AS1_SI / NC")
replotGSEA(path = gsea.out, gene.set = "HALLMARK_DNA_REPAIR", class.name = "DARS-AS1_SI / NC")

## Figure 5B. Enrichment plot ======================================
SigPval = 0.01

# Initial dataset
enrichChart = read.delim(paste(c(DavidFileName, "txt"), collapse="."))

SigEnrichItem = arrange(enrichChart[which(enrichChart$PValue <= SigPval),],Category, PValue)
SigEnrichItem$logP <- -log10(SigEnrichItem$PValue)

pdf(paste(c(DavidFileName, "pdf"), collapse="."), width=18, height=8.5)
ggbarplot(SigEnrichItem, x="Term", y="logP", ylab = "-log10(PValue)", fill = "Category", color = "white", 
          main = title, 
          palette = "npg",
          sort.val = "desc",
          sort.by.groups = T,
          x.text.angle = 0,
          rotate = TRUE
          )
dev.off()
