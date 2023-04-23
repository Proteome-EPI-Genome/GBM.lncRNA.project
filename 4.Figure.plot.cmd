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

## Figure 7B. Motif Logo plot ======================================
MF.CATC <- read.table("eCLIP/111519_YC_4CLIP/yeo.lab.out/111519_YC_4CLIP.motifAnalysis.CATC.txt", header=T)
colnames(MF.CATC) = c("A", "C", "G", "U")
CATC_motif <- new("pcm", mat = as.matrix(t(MF.CATC)), name="CAUC")
pdf("111519_YC_4CLIP.motifAnalysis.CATC.pdf", width=6, height=4)
plot(CATC_motif)
dev.off()

MF.TTACCATC <- read.table("eCLIP/111519_YC_4CLIP/yeo.lab.out/111519_YC_4CLIP.motifAnalysis.TTACCATC.txt", header=T)
colnames(MF.TTACCATC) = c("A", "C", "G", "U")
TTACCATC_motif <- new("pcm", mat = as.matrix(t(MF.TTACCATC)), name="UUACCAUC")
pdf("111519_YC_4CLIP.motifAnalysis.TTACCATC.pdf", width=6, height=4)
plot(TTACCATC_motif)
dev.off()

## Figure 7D. signal track plot ======================================
UCSC

## Figure S1C & S1D. Gini index ======================================
gi = rep(0,dim(InputData_counts)[2])
for (i in 1:dim(InputData_counts)[2]) {
  gi[i] = gini(log2(InputData_counts[,i] + 1))
}
pdf(paste(OutputPrefix, "_Gini.pdf", sep = ""), width=11, height=8.5)
par(mar=c(7,6,2,2)+0.1,mgp=c(5,1,0))
barplot(gi, names.arg = colnames(InputData_counts), ylab = "Gini Index", col='steelblue', las=2)
dev.off()

## Figure S1E & S1F. Screen Compare (QC LFC boxplot for control & lncRNAs) ======================================
#U251_180529_D21_D0.sgrna_summary
OutputPrefix = "U251"
sgRNA.level.data = U251_180529_D21_D0.sgrna_summary
sgRNA.level.data$positive = ifelse(startsWith(sgRNA.level.data$sgrna, "POS"), "POS_CTRL", ifelse(startsWith(sgRNA.level.data$sgrna, "NEG"), "NEG_CTRL", "lncRNAs"))
table(sgRNA.level.data$positive)

LFC.box1 <- ggboxplot(sgRNA.level.data,
                      x = "positive",
                      y = c("LFC"),
                      ylim = c(-7.5, 5),
                      xlab = OutputPrefix,
                      ylab = "log2(Fold Change)",
                      order = c("lncRNAs", "NEG_CTRL", "POS_CTRL"),
                      color = "positive", palette = rainbow(length(unique(sgRNA.level.data$positive))),
                      add = "median_iqr") +
  theme(legend.position = "none")
#U87_180529_D21_D0.sgrna_summary
LFC.box2

pdf("3.U251.U87.lncRNA.CRISPRi.Screen.LFC.compare.pdf", width=7.5, height=5)
grid.arrange(LFC.box1, LFC.box2, ncol=2)
dev.off()


## Additional Analysis
# subtype
CLINICAL = "TCGA.GBM.3Class_subtype[1-s2.0-S1535610817302532-mmc5].txt"
ClinicalInfo <- read.delim(CLINICAL, check.names = FALSE)
ID2SubType = unique(ClinicalInfo[, c("sampleId", "Group")])
row.names(ID2SubType) = ID2SubType$sampleId

GeneExpr4box <- data.frame(
  GeneExpr = as.numeric(TCGA.GBM.GTEx.Brain.exp[grep(GID, row.names(TCGA.GBM.GTEx.Brain.exp)), ]),
  SampleID = colnames(TCGA.GBM.GTEx.Brain.exp)
)
GeneExpr4box$SampleType = ID2SubType[GeneExpr4box$SampleID, "Group"]
GeneExpr4box$SampleType[grep("^GTEX", GeneExpr4box$SampleID)] = "Normal"
GeneExpr4box = GeneExpr4box[-which(is.na(GeneExpr4box$SampleType)), ]
head(GeneExpr4box)
table(GeneExpr4box$SampleType)

compare.out = compare_means(GeneExpr~SampleType, data = GeneExpr4box, paired = FALSE)
compare.out
sig.pairs = compare.out[which(compare.out$p.signif != "ns"), c("group1", "group2")]
my_comparisons <- vector("list", nrow(sig.pairs))
for (i in 1:nrow(sig.pairs)) {
  my_comparisons[[i]] <- c(as.character(sig.pairs[i, "group1"]), as.character(sig.pairs[i, "group2"]))
}

pdf(paste(c("Boxplot", GeneName, "GBM.GTEx.3subtypes.v2.pdf"), collapse="."), width=11.5, height=7.5)
ggboxplot(GeneExpr4box, x = "SampleType",
          y = "GeneExpr",
          color = "SampleType", palette = "jco",
          ylab = "Expression\nlog2NormalizedCount",
          add = "dotplot",                              # Add "dotplot"
          add.params = list(binwidth = 0.1, dotsize = 0.3)
) + stat_compare_means(comparisons = my_comparisons)
dev.off()

## Purity Analysis (GBM & LGG)
#InputData_counts
#PurityPloidy
#GeneName = "DARS-AS1"
#GID = "ENSG00000231890"

#Purity Filter (> 0.7, 2023.0303): PurityPloidy = PurityPloidy[-which(PurityPloidy$purity < 0.7), ]

CorData = as.data.frame(t(InputData_counts[GID, which(substr(colnames(InputData_counts), 1, 12) %in% row.names(PurityPloidy))]))
CorData$purity = PurityPloidy[substr(row.names(CorData), 1, 12), "purity"]
CorData$ploidy = PurityPloidy[substr(row.names(CorData), 1, 12), "ploidy"]

head(CorData)

cor.test(CorData$ENSG00000231890, CorData$purity, alternative = "two.sided", method = "pearson")

data4boxLGG <- data.frame(
  ENSG00000231890 = as.numeric(CorData$ENSG00000231890),
  group = rep("LGG", length(CorData$ENSG00000231890))
)
data4boxGBM <- data.frame(
  ENSG00000231890 = as.numeric(CorData$ENSG00000231890),
  group = rep("GBM", length(CorData$ENSG00000231890))
)
data4boxNorm <- data.frame(
  ENSG00000231890 = as.numeric(InputData_counts["ENSG00000231890", grep("-11A.rep", colnames(InputData_counts))]),
  group = rep("Normal", length(grep("-11A.rep", colnames(InputData_counts))))
)
data4box = rbind(data4boxNorm, data4boxLGG, data4boxGBM)
data4box$ENSG00000231890 = log2(data4box$ENSG00000231890)

pdf("purity.gt.0.7.v3.pdf", width=11.5, height=8.5)
ggboxplot(data4box2,
          x = "group",
          y = "ENSG00000231890",
          xlab = "Sample Types",
          ylab = "Normalized count",
          color = "group",
          order = c("Normal", "LGG", "GBM"),
          add = "dotplot",                              # Add dotplot
          add.params = list(binwidth = 0.1, dotsize = 0.2)
) +
  theme(legend.position = "top")
dev.off()

compare.out = compare_means(ENSG00000231890~group, data = data4box2, paired = FALSE)
compare.out

fold.change(
  median(data4box2$ENSG00000231890[which(data4box2$group == "Normal")]),
  median(data4box2$ENSG00000231890[which(data4box2$group == "GBM")])
  )

## 4 COX multivariate analysis
ClinicalInfoTCGA
ClinicalInfo

ID2SurvivalData = unique(ClinicalInfoTCGA[, c("sampleID", "_EVENT", "_OS", "gender", "age_at_initial_pathologic_diagnosis")])
row.names(ID2SurvivalData) = ID2SurvivalData$sampleID
ID2SurvivalData[, c("IDH1_status", "Group", "OS", "Age", "Gender")] = ClinicalInfo[row.names(ID2SurvivalData), c("IDH1_status", "Group", "OS", "Age", "Gender")]

coxData = as.data.frame(t(InputData_counts[GID, ]))
coxData$age = ID2SurvivalData[substr(row.names(coxData), 1, 15), "age_at_initial_pathologic_diagnosis"]
coxData$gender = ID2SurvivalData[substr(row.names(coxData), 1, 15), "gender"]
coxData$IDH1status = ID2SurvivalData[substr(row.names(coxData), 1, 15), "IDH1_status"]
coxData$Group = ID2SurvivalData[substr(row.names(coxData), 1, 15), "Group"]
coxData$days = ID2SurvivalData[substr(row.names(coxData), 1, 15), "_OS"]
coxData$status = ID2SurvivalData[substr(row.names(coxData), 1, 15), "_EVENT"]

# REMOVE survival time "NA"
if(length(which(is.na(coxData$days))) > 0){
  coxData = coxData[-which(is.na(coxData$days)), ]
}
# REMOVE survival LiveStatus "NULL"
if(length(which(is.na(coxData$status))) > 0){
  coxData = coxData[-which(is.na(coxData$status)), ]
}
# REMOVE survival subtype.Group "NULL"
if(length(which(is.na(coxData$Group))) > 0){
  coxData = coxData[-which(is.na(coxData$Group)), ]
}

coxData$gender = ifelse(toupper(coxData$gender) == "MALE", 1, 0)
coxData$IDH1status = ifelse(toupper(coxData$IDH1status) == "WT", 0, 1)
coxData$Group = ifelse(toupper(coxData$Group) == "CL", 1, ifelse(toupper(coxData$Group) == "MS", 2, ifelse(toupper(coxData$Group) == "PN", 3, 0)))

for(i in 1:(dim(coxData)[2]-2)){
  result=coxph(Surv(days,status)~coxData[,i],coxData)
  p<-summary(result)[[7]][5]
  coef<-summary(result)[[7]][1]
  write.table(matrix(c(colnames(coxData)[i],p,coef),1,3), "coxresult.v2.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=T)
}

coxData$regression = 0.00975363602355131*coxData$ENSG00000231890 + 0.0252234718110747*coxData$age
coxData$regressionGroup = ifelse(coxData$regression >= median(coxData$regression), "High", "Low")
sfit<-survfit(Surv(days,status)~regressionGroup,data=coxData)
sfit
ggsurv <- ggsurvplot(sfit,
                     conf.int = T, conf.int.style = "step",
                     censor.size = 6, size = 1,
                     pval = T,
                     pval.method = T,
                     log.rank.weights = "1",
                     surv.median.line ="hv",
                     legend = c(0.8, 0.95), legend.title = "")
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 0, y = 0.4, # x and y coordinates of the text
                    label = paste(paste(c("High", "Low"), table(coxData$regressionGroup), sep = ":"), collapse = ", "),
                    hjust = "left",
                    size = 5)
res <- ggpar(ggsurv, 
             font.main = c(12, "bold"),
             font.x = c(12, "bold"),
             font.y = c(12, "bold"),
             font.caption = c(12, "bold"), 
             font.legend = c(12), 
             font.tickslab = c(12))
pdf("Cox.regression.Median.High-Low.survival.pdf", width=6.5, height=6.5)
print(res, newpage = FALSE)
dev.off()

## infiltration
infiltration_estimation_for_tcga
GBM_sample_sheet.2020.02.25
LGG_sample_sheet.2020.02.25

# CIBERSORT, CIBERSORT-ABS, EPIC, MCPCOUNTER, QUANTISEQ, TIMER, XCELL
methodType = "QUANTISEQ"
data4box = infiltration_estimation_for_tcga[
  c(
    substr(GBM_sample_sheet.2020.02.25$`Sample ID`[which(endsWith(GBM_sample_sheet.2020.02.25$`Sample Type`, "Tumor"))], 1, 15),
    substr(LGG_sample_sheet.2020.02.25$`Sample ID`[which(endsWith(LGG_sample_sheet.2020.02.25$`Sample Type`, "Tumor"))], 1, 15)
  ),
  which(endsWith(colnames(infiltration_estimation_for_tcga), methodType))
]

data4box[
  substr(GBM_sample_sheet.2020.02.25$`Sample ID`[which(endsWith(GBM_sample_sheet.2020.02.25$`Sample Type`, "Tumor"))], 1, 15),
  "TumorType"
] = "GBM"
data4box[
  substr(LGG_sample_sheet.2020.02.25$`Sample ID`[which(endsWith(LGG_sample_sheet.2020.02.25$`Sample Type`, "Tumor"))], 1, 15),
  "TumorType"
] = "LGG"


df <- data.frame(
  infiltration = as.vector(as.matrix(data4box[, -dim(data4box)[2]])),
  CellType = rep(
    substr(
      colnames(data4box[, -dim(data4box)[2]]), 1 , 
      unlist(gregexpr(pattern ='_', colnames(data4box[, -dim(data4box)[2]])))-1), 
    rep(dim(data4box)[1],dim(data4box)[2]-1)),
  TumorType = rep(data4box$TumorType, dim(data4box)[2]-1)
)
df = df[-which(is.na(df$TumorType)) ,]

library(ggpubr)

pdf(paste("GBM.LGG.infiltration.", methodType, ".BoxPlot.pdf", sep = ""), width=11, height=8.5)
ggboxplot(df, "CellType", "infiltration", color = "TumorType",
          palette = c("#00AFBB", "#E7B800")) + 
  rotate_x_text(45)
dev.off()
