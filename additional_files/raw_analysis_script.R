library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)
library(rtracklayer)
library(edgeR)
library(dplyr)
library(GenomicFeatures)
library(diffloop)
library(Sushi)
library(GenomicInteractions)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plotgardener)
library(RColorBrewer)
library(FactoMineR)
library(factoextra)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(rstatix)
library(ChIPpeakAnno)
library(chipseq)
library(ShortRead)
library(readr)
library(tidyr)
library(tidyverse)
library(BSgenome)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(EnsDb.Hsapiens.v86)
annoData <- toGRanges(EnsDb.Hsapiens.v86, feature = "gene")
library(stringr)
library(enrichplot)

#loading diffloop matrix

#reading HiChIP data

data <- read.delim("/Volumes/ucmm/GrpSR/Itzel&Chaitali/HiChIP_H3K4me3/hichipper_loops/diffloops_lost_gain_common/Allstats_interactions_rosum2loops-all.significant.loops_withoutU3047.txt")

data_chr15 <- data %>% dplyr::filter(chr.x == "chr15") #8972 loopa
head(data_chr15)
"left right     name_inter loopWidth    mango.FDR      mango.P hAstro U3008 U3009 U3013 U3028 U3031 U3039 U3042 U3054 U3073 U3078 U3085 U3086 U3118 U3121 GBM_counts       FC
184974 51609 51610 51609_to_51610     20302 9.926786e-04 6.096431e-04      3     0     0    14     0     0     0     0     0     0     0     0     0     0     0         14 4.666667
184975 51609 51612 51609_to_51612    377977 9.263936e-05 4.832287e-05      0     0     0     2     0     0     0     0     0     0     0     0     0     0     0          2      Inf
184976 51611 51613 51611_to_51613    362341 1.001896e-04 5.268039e-05      0     0     0     2     0     0     0     0     0     0     0     0     0     0     0          2      Inf
184977 51612 51615 51612_to_51615    792619 5.180832e-08 1.517130e-08      0     0     0     1     2     0     0     0     0     0     0     0     0     0     0          3      Inf
184978 51613 51615 51613_to_51615    771853 2.378743e-05 1.058543e-05      0     0     0     0     2     0     0     0     0     0     0     0     0     0     0          2      Inf
184979 51614 51615 51614_to_51615    750126 2.531266e-05 1.135395e-05      0     0     0     2     0     0     0     0     0     0     0     0     0     0     0          2      Inf
         log2FC seqnames.x  start.x    end.x width.x chr.x              anchorID.L seqnames.y  start.y    end.y width.y chr.y              anchorID.R distance
184974 2.222392         15 17050582 17070200   19619 chr15 chr15:17050582-17070200         15 17075230 17086156   10927 chr15 chr15:17075230-17086156    15956
184975      Inf         15 17050582 17070200   19619 chr15 chr15:17050582-17070200         15 17435126 17441610    6485 chr15 chr15:17435126-17441610   371410
184976      Inf         15 17096047 17097538    1492 chr15 chr15:17096047-17097538         15 17457775 17460493    2719 chr15 chr15:17457775-17460493   362955
184977      Inf         15 17435126 17441610    6485 chr15 chr15:17435126-17441610         15 18230127 18231848    1722 chr15 chr15:18230127-18231848   790238
184978      Inf         15 17457775 17460493    2719 chr15 chr15:17457775-17460493         15 18230127 18231848    1722 chr15 chr15:18230127-18231848   771355
184979      Inf         15 17480098 17481625    1528 chr15 chr15:17480098-17481625         15 18230127 18231848    1722 chr15 chr15:18230127-18231848   750223"
colnames(data_chr15)
#[1] "left"       "right"      "name_inter" "loopWidth"  "mango.FDR"  "mango.P"    "hAstro"     "U3008"      "U3009"      "U3013"      "U3028"      "U3031"      "U3039"      "U3042"     
#[15] "U3054"      "U3073"      "U3078"      "U3085"      "U3086"      "U3118"      "U3121"      "GBM_counts" "FC"         "log2FC"     "seqnames.x" "start.x"    "end.x"      "width.x"   
#[29] "chr.x"      "anchorID.L" "seqnames.y" "start.y"    "end.y"      "width.y"    "chr.y"      "anchorID.R" "distance" 


#adding loop IDs
data_chr15$loopID <- paste(data_chr15$chr.x, data_chr15$start.x, sep = ":")
data_chr15$loopID <- paste(data_chr15$loopID, data_chr15$end.x, sep = "-")
data_chr15$loopID <- paste(data_chr15$loopID, data_chr15$chr.y, sep = "_")
data_chr15$loopID <- paste(data_chr15$loopID, data_chr15$start.y, sep = ":")
data_chr15$loopID <- paste(data_chr15$loopID, data_chr15$end.y, sep = "-")

selected_cols <- c("chr.x", "start.x", "end.x", "chr.y","start.y", "end.y", "" )

data_chr15.1 <- data_chr15[, c(29,26,27,35,32,33,7,11,4,5,6,38) ] #make a simple data frame
data_chr15.1$total_counts <- rowSums(data_chr15.1[, c(7,8)]) # sum counts
data_chr15.2 <- data_chr15.1[!(data_chr15.1$total_counts == 0), ] #3511

write.table(data_chr15.1, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_full.txt", sep = "\t", quote = F)

write.table(data_chr15.2[, -c(12,13)], "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe.txt", sep = "\t", quote = F)



#overlap_loops


hAstro_loops <- data_chr15.2 %>% dplyr::filter(hAstro != 0) #1088

hAstro <- hAstro_loops$loopID


U3028_loops <- data_chr15.2 %>% dplyr::filter(U3028 != 0) #3085

U3028 <- U3028_loops$loopID

#Overlap matrix with loop IDs

lt1 = list("hAstro"=hAstro, "U3028"=U3028)
mat1 <- list_to_matrix(lt1)
ltdf1 <- as.data.frame(mat1) #a matrix with binaries for TRUE and FALSe on presense and absence of loops
write.csv(ltdf1, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe_upsetmatrix.csv") 

mat2 <- make_comb_mat(lt1)


UpSet(mat2, top_annotation = upset_top_annotation(mat2, add_numbers= TRUE))

UpSet(mat2, set_order = c("hAstro", "U3028"), top_annotation = upset_top_annotation(mat2, add_numbers= TRUE), right_annotation = upset_right_annotation(mat2, add_numbers=TRUE))



png("/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe_upsetmatrix_UPSETPLOT.png", width = 14, height = 4, res = 600, units = 'in')
UpSet(mat2, set_order = c("hAstro", "U3028"), top_annotation = upset_top_annotation(mat2, add_numbers= TRUE), right_annotation = upset_right_annotation(mat2, add_numbers=TRUE))
par(mar=c(6,5,5,6) + 0.1)
dev.off()


p <- ggVennDiagram(lt1) + scale_x_continuous(expand = expansion(mult = .2)) + scale_fill_distiller(palette = "Set3")
p 


png("/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe_upsetmatrix_VENNDIAGRAM.png", width = 14, height = 4, res = 600, units = 'in')
p
par(mar=c(6,5,5,6) + 0.1)
dev.off()


#assigning loop types

head(ltdf1)
"                                                     hAstro U3028
chr15:100047761-100064820_chr15:100166768-100177483      1     0
chr15:100047761-100064820_chr15:100327465-100347878      1     0
chr15:100047761-100064820_chr15:100599116-100606695      1     0
chr15:100111670-100120379_chr15:100599116-100606695      0     1
chr15:100221330-100246394_chr15:100379020-100387966      1     0
chr15:100327465-100347878_chr15:100348754-100351622      0     1"

ltdf1$presence <- rowSums(ltdf1[, c(1,2)])

#let's annotate
ltdf1 <- ltdf1 %>% dplyr::mutate(annotation = case_when((presence == 2) ~ "common",
                                                        (presence == 1 & hAstro == 1) ~"hAstro.specific",
                                                        (presence == 1 & U3028 == 1)~ "GB.specific"))

table(ltdf1$annotation)

#common     GB.specific hAstro.specific 
#662            2423             426 

write.csv(ltdf1, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe_upsetmatrix_annotated.csv") 

#now annotate your loop matrix

data_chr15.2$annotation <- ltdf1$annotation[match(data_chr15.2$loopID, rownames(ltdf1))]
data_chr15.2 <- data_chr15.2 %>% dplyr::mutate(colour = case_when(data_chr15.2$annotation == "common" ~ "#a0a0a0",
                                                                  data_chr15.2$annotation ==  "hAstro.specific" ~ "#1f78b4",
                                                                  data_chr15.2$annotation == "GB.specific"~"#e31a1c"))


write.table(data_chr15.2, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe_annotated.txt", sep = "\t", quote = F)


#prepare a interact file for visualization in the UCSC browser

#first three columns are 
"interaction between two regions"
#( 
#  string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records" 
#  uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region" 
#  uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
#  string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
#  uint score;          "Score (0-1000)"
#  double value;        "Strength of interaction or other data value. Typically basis for score"
#  string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
#  string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
#  string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
#  uint sourceStart;    "Start position in chromosome of source/lower/this region"
#  uint sourceEnd;      "End position in chromosome of source/lower/this region"
#  string sourceName;   "Identifier of source/lower/this region"
#  string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
#  string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
#  uint targetStart;    "Start position in chromosome of target/upper/this region"
#  uint targetEnd;      "End position in chromosome of target/upper/this region"
#  string targetName;   "Identifier of target/upper/this region"
#  string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
  
#)



df <- data_chr15.2[, c(1,2,3)]
df$name <- paste("Chaitali", 1:nrow(df), sep = "_")
df$score <- data_chr15.2$total_counts
df$value <- -log10(data_chr15.2$mango.FDR)
df$value <- ifelse(df$value == "Inf", 20, df$value)
class(df$value)
#[1] "numeric"

df$exp <- data_chr15.2$annotation
df$color <- data_chr15.2$colour
df$sourceChrom <- df$chr.x
df$sourceStart <- df$start.x
df$sourceEnd <- df$end.x
df$sourceName <- "anchor1"
df$sourceStrand <- "+"
df$targetChrom <- data_chr15.2$chr.y
df$targetStart <- data_chr15.2$start.y
df$targetEnd <- data_chr15.2$end.y
df$targetName <- "anchor2"
df$targetStrand <- "+"

head(df)

write.table(df, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/data_chr15.2_studentsdataframe_annotated_UCSCinteract.txt", sep = "\t", quote = F)

#on terminal
#tail -n +2 data_chr15.2_studentsdataframe_annotated_UCSCinteract.txt > tmpfile
#less tmpfile
#cut -f 2-19 tmpfile > tmpfile1
#less tmpfile1
#echo -e " track type=interact name="H3K4Me3.HiChIP" description="Chr15 subset" interactDirectional=true maxHeightPixels=200:100:50 visibility=full itemRgb="On" " | cat - tmpfile1  > data_chr15.2_studentsdataframe_annotated_UCSCbrowser.interact


#look at annotations for the specific identified groups

#prepare bedpe files  #8 coloumns - chr.x. start.x, end.x, chr.y, start.y, end.y, name,counts/score/fdr
#astro specific 
astro_specific <- data_chr15.2[data_chr15.2$annotation == "hAstro.specific", ]

head(astro_specific)

bedpe <- astro_specific[, c(1,2,3,4,5,6, 12,11)]
write.table(bedpe, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/hAstro_sp_chr15_studentsdataframet.txt", sep = "\t", quote = F)


GB_specific <- data_chr15.2[data_chr15.2$annotation == "GB.specific", ]

head(GB_specific)

bedpe <- GB_specific[, c(1,2,3,4,5,6, 12,11)]
write.table(bedpe, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/GB_sp_chr15_studentsdataframet.txt", sep = "\t", quote = F)

#less GB_sp_chr15_studentsdataframet.txt
#cat GB_sp_chr15_studentsdataframet.txt | wc -l #2424
#tail -n +2 GB_sp_chr15_studentsdataframet.txt > tmpfile
#cut -f 2-9 tmpfile > GB_sp_chr15_studentsdataframet.bedpe
#less GB_sp_chr15_studentsdataframet.bedpe
#cat GB_sp_chr15_studentsdataframet.bedpe | wc -l #2423


common <- data_chr15.2[data_chr15.2$annotation == "common", ]

head(common)

bedpe <- common[, c(1,2,3,4,5,6, 12,11)]
write.table(bedpe, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/common_chr15_studentsdataframet.txt", sep = "\t", quote = F)


#Run Genomic annotations
dfloops <- file.path("/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/GB_sp_chr15_studentsdataframet.bedpe")
dfinter <- makeGenomicInteractionsFromFile(dfloops, type = 'bedpe', experiment_name = "GB specific")
head(interactionCounts(dfinter))
#[1] 1.517130e-08 1.058543e-05 1.445453e-05 6.808408e-05 0.000000e+00 6.242519e-09
median_distance_interactions <- median(calculateDistances(dfinter, method = "midpoint"))
median_distance_interactions
#[1] 335226


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
refseq.genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
refseq.transcripts = transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")
non_pseudogene = names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) 
refseq.transcripts = refseq.transcripts[non_pseudogene] 

refseq.promoters = promoters(refseq.transcripts, upstream=2000, downstream=2000)
# unlist object so "strand" is one vector
refseq.transcripts.ul = unlist(refseq.transcripts) 
# terminators can be called as promoters with the strand reversed
strand(refseq.transcripts.ul) = ifelse(strand(refseq.transcripts.ul) == "+", "-", "+") 
refseq.terminators.ul = promoters(refseq.transcripts.ul, upstream=1000, downstream=1000) 
# change back to original strand
strand(refseq.terminators.ul) = ifelse(strand(refseq.terminators.ul) == "+", "-", "+") 
# `relist' maintains the original names and structure of the list
refseq.terminators = relist(refseq.terminators.ul, refseq.transcripts)


annotation.features = list(promoter=refseq.promoters, 
                           terminator=refseq.terminators, 
                           gene.body=refseq.transcripts)
annotateInteractions(dfinter, annotation.features)

annotationFeatures(dfinter)
#[1] "distal"     "promoter"   "gene.body"  "terminator"

categoriseInteractions(dfinter)
"              category count
1          distal-distal    26
2        distal-promoter   316
3       distal-gene.body    21
4      distal-terminator    12
5      promoter-promoter  1530
6     promoter-gene.body   346
7    promoter-terminator   132
8    gene.body-gene.body    29
9   gene.body-terminator     9
10 terminator-terminator     2"

dfinter_df <- data.frame(dfinter)
head(dfinter_df[, c(1,2,3,6,10,11,12,15,19,20)])
"seqnames1   start1     end1 node.class1 seqnames2   start2     end2 node.class2       counts                                            name
1     chr15 17435127 17441610      distal     chr15 18230128 18231848      distal 1.517130e-08 chr15:17435126-17441610_chr15:18230127-18231848
2     chr15 17457776 17460493      distal     chr15 18230128 18231848      distal 1.058543e-05 chr15:17457775-17460493_chr15:18230127-18231848
3     chr15 18465421 18523463      distal     chr15 19167039 19170276      distal 1.445453e-05 chr15:18465420-18523463_chr15:19167038-19170276
4     chr15 18754403 18770002      distal     chr15 19070205 19108832      distal 6.808408e-05 chr15:18754402-18770002_chr15:19070204-19108832
5     chr15 22491926 22497810      distal     chr15 22836948 22842739    promoter 0.000000e+00 chr15:22491925-22497810_chr15:22836947-22842739
6     chr15 22491926 22497810      distal     chr15 23560791 23572771    promoter 6.242519e-09 chr15:22491925-22497810_chr15:23560790-23572771"

dfinter_df$start1 <- dfinter_df$start1 -1
dfinter_df$start2 <- dfinter_df$start2 -1

head(dfinter_df[, c(1,2,3,6,10,11,12,15,19,20)])

"seqnames1   start1     end1 node.class1 seqnames2   start2     end2 node.class2       counts                                            name
1     chr15 17435126 17441610      distal     chr15 18230127 18231848      distal 1.517130e-08 chr15:17435126-17441610_chr15:18230127-18231848
2     chr15 17457775 17460493      distal     chr15 18230127 18231848      distal 1.058543e-05 chr15:17457775-17460493_chr15:18230127-18231848
3     chr15 18465420 18523463      distal     chr15 19167038 19170276      distal 1.445453e-05 chr15:18465420-18523463_chr15:19167038-19170276
4     chr15 18754402 18770002      distal     chr15 19070204 19108832      distal 6.808408e-05 chr15:18754402-18770002_chr15:19070204-19108832
5     chr15 22491925 22497810      distal     chr15 22836947 22842739    promoter 0.000000e+00 chr15:22491925-22497810_chr15:22836947-22842739
6     chr15 22491925 22497810      distal     chr15 23560790 23572771    promoter 6.242519e-09 chr15:22491925-22497810_chr15:23560790-23572771"

dfinter_df <- dfinter_df %>% dplyr::mutate(interaction.type =
                                             case_when((dfinter_df$node.class1 == "promoter" &dfinter_df$node.class2 == "promoter"~"pp"),
                                                       (dfinter_df$node.class1 == "promoter" &dfinter_df$node.class2 == "distal"~"pe"),
                                                       (dfinter_df$node.class1 == "distal" &dfinter_df$node.class2 == "promoter"~"pe"),
                                                       (dfinter_df$node.class1 == "promoter" &dfinter_df$node.class2 == "gene.body"~"pe"),
                                                       (dfinter_df$node.class1 == "gene.body" &dfinter_df$node.class2 == "promoter"~"pe"),
                                                       (dfinter_df$node.class1 == "promoter" &dfinter_df$node.class2 == "terminator"~"pe"),
                                                       (dfinter_df$node.class1 == "terminator" &dfinter_df$node.class2 == "promoter"~"pe"),
                                                       (dfinter_df$node.class1 == "distal" &dfinter_df$node.class2 == "distal"~"ee"),
                                                       (dfinter_df$node.class1 == "distal" &dfinter_df$node.class2 == "gene.body"~"others"),
                                                       (dfinter_df$node.class1 == "gene.body" &dfinter_df$node.class2 == "distal"~"others"),
                                                       (dfinter_df$node.class1 == "gene.body" &dfinter_df$node.class2 == "gene.body"~"others"),
                                                       (dfinter_df$node.class1 == "distal" &dfinter_df$node.class2 == "terminator"~"others"),
                                                       (dfinter_df$node.class1 == "terminator" &dfinter_df$node.class2 == "distal"~"others"),
                                                       (dfinter_df$node.class1 == "terminator" &dfinter_df$node.class2 == "terminator"~"others"),
                                                       (dfinter_df$node.class1 == "gene.body" &dfinter_df$node.class2 == "terminator"~"others"),
                                                       (dfinter_df$node.class1 == "terminator" &dfinter_df$node.class2 == "gene.body"~"others")))




dfinter_df_sh <- dfinter_df[, c(1,2,3,6,10,11,12,15,19,20,21)]
head(dfinter_df_sh)
"seqnames1   start1     end1 node.class1 seqnames2   start2     end2 node.class2       counts                                            name interaction.type
1     chr15 17435126 17441610      distal     chr15 18230127 18231848      distal 1.517130e-08 chr15:17435126-17441610_chr15:18230127-18231848               ee
2     chr15 17457775 17460493      distal     chr15 18230127 18231848      distal 1.058543e-05 chr15:17457775-17460493_chr15:18230127-18231848               ee
3     chr15 18465420 18523463      distal     chr15 19167038 19170276      distal 1.445453e-05 chr15:18465420-18523463_chr15:19167038-19170276               ee
4     chr15 18754402 18770002      distal     chr15 19070204 19108832      distal 6.808408e-05 chr15:18754402-18770002_chr15:19070204-19108832               ee
5     chr15 22491925 22497810      distal     chr15 22836947 22842739    promoter 0.000000e+00 chr15:22491925-22497810_chr15:22836947-22842739               pe
6     chr15 22491925 22497810      distal     chr15 23560790 23572771    promoter 6.242519e-09 chr15:22491925-22497810_chr15:23560790-23572771               pe"
write.csv(dfinter_df_sh, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/GB_sp_chr15_GA.csv")

dfinter_pp <- dfinter[is.pp(dfinter)]
dfinter_pp_genes <- dfinter_pp@regions@elementMetadata@listData$promoter.id
dfinter_pp_genes_un <- unlist(dfinter_pp_genes)
dfinter_pp_genes_un <- na.omit(dfinter_pp_genes_un)
dfinter_pp_genes_dedup <- dfinter_pp_genes_un[!duplicated(dfinter_pp_genes_un)]
converted <- bitr(dfinter_pp_genes_dedup, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db, drop = FALSE) #14054
converted_pp <- converted
converted_pp$annotation <- "pp"

converted_pp <- converted_pp[!duplicated(converted_pp$ENTREZID),] 

dfinter_pd <- dfinter[is.pd(dfinter)]
dfinter_pd_genes <- dfinter_pd@regions@elementMetadata@listData$promoter.id
dfinter_pd_genes_un <- unlist(dfinter_pd_genes)
dfinter_pd_genes_un <- na.omit(dfinter_pd_genes_un)
dfinter_pd_genes_dedup <- dfinter_pd_genes_un[!duplicated(dfinter_pd_genes_un)] 
converted_pd <- bitr(dfinter_pd_genes_dedup, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db, drop = FALSE)

dfinter_pt <- dfinter[is.pt(dfinter)]
dfinter_pt_genes <- dfinter_pt@regions@elementMetadata@listData$promoter.id
dfinter_pt_genes_un <- unlist(dfinter_pt_genes)
dfinter_pt_genes_un <- na.omit(dfinter_pt_genes_un)
dfinter_pt_genes_dedup <- dfinter_pt_genes_un[!duplicated(dfinter_pt_genes_un)] 
converted_pt <- bitr(dfinter_pt_genes_dedup, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db, drop = FALSE)

dfinter_pg <- dfinter[isInteractionType(dfinter, "promoter", "gene.body")]
dfinter_pg_genes <- dfinter_pg@regions@elementMetadata@listData$promoter.id
dfinter_pg_genes_un <- unlist(dfinter_pg_genes)
dfinter_pg_genes_un <- na.omit(dfinter_pg_genes_un)
dfinter_pg_genes_dedup <- dfinter_pg_genes_un[!duplicated(dfinter_pg_genes_un)] #2715
converted_pg <- bitr(dfinter_pg_genes_dedup, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db, drop = FALSE)

converted_pe <- rbind(converted_pd, converted_pg, converted_pt)
converted_pe <- converted_pe[!duplicated(converted_pe$ENTREZID),]
converted_pe$annotation <- "pe"

converted <- rbind(converted_pp, converted_pe)

#unique genes annotated to anchors 
length(unique(converted$ENTREZID)) #[1] 622


write.csv(converted, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/GB_sp_chr15_GA_genes.csv")


#chromatin compaction

data_summary <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}


data_chr15.2$log2.loopWidth <- log2(data_chr15.2$loopWidth)



hist(data_chr15.2$loopWidth)
png("/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/df_chr15_group_loopdistribution.png", width = 3, height = 3, res = 600, units = 'in')
hist(data_chr15.2$loopWidth, main = "lop length dist.", xlab = "loop length (bp)")
par(mar=c(8, 6, 6, 5) + 0.1)
dev.off()

hist(data_chr15.2$log2.loopWidth)
png("/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/df_chr15_group_log2loopdistribution.png", width = 3, height = 3, res = 600, units = 'in')
hist(data_chr15.2$log2.loopWidth, main = "log2 loop length dist.", xlab = "log2 of loop length (bp)")
par(mar=c(8, 6, 6, 5) + 0.1)
dev.off()


data_chr15.2$annotation <- factor(data_chr15.2$annotation, levels = c("hAstro.specific", "GB.specific", "common") )

stat.test1 <- data_chr15.2 %>%
  wilcox_test(log2.loopWidth~annotation, p.adjust.method = "BH") %>%
  add_significance()
stat.test1


stat.test1 <- stat.test1 %>% add_xy_position(x = "annotation")
stattest1df <- as.data.frame(stat.test1)
stattest1df$groups <- as.character(stattest1df$groups)
write.csv(stattest1df, "/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/Wilcoxon.test.csv")

p1 <- ggplot(df, aes(x=data_chr15.2$annotation, y=data_chr15.2$loopWidth)) + geom_violin(fill = "gray", size = 0.5, linewidth=0.3)
p1

p2 <- p1 + stat_summary(fun.data=data_summary, size = 0.1)
p2

p3 <- p2 + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
p3
p4 <- p3  + ggtitle("Loop length") +
  theme(plot.title = element_text(hjust = 0.5, size= 10)) + xlab("Annotations") + ylab("loop length (bp)") + theme(axis.text.x = element_text(size = 6, angle = 0),
                                                                                                                          axis.text.y = element_text(size = 8),axis.title = element_text(size = 10),
                                                                                                                          axis.title.x = element_text(size = 8),  axis.title.y = element_text(size = 8))

p4


p5 <- p4 + stat_pvalue_manual(
  stat.test1, label = "Wilcoxon-test, p.adj = {p.adj} {p.adj.signif}",
  vjust = -1.3, bracket.nudge.y = 0.5, size = 1.5,y.position = seq(from=2050000,by=400000,length.out=3)
) + ylim(5000, max(seq(from=2050000,by=400000,length.out=4)))

p5


png("/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/looplength.png", width = 4, height = 6, res = 600, units = 'in')
p5
par(mar=c(8, 6, 6, 5) + 0.1)
dev.off()


#fetch anchors for analysis
#for GB specific loops
#anchor.x

anchor.x <- GRanges(seqnames = GB_specific$chr.x, IRanges(start = GB_specific$start.x, end =GB_specific$end.x))

anchor.x
"GRanges object with 2423 ranges and 0 metadata columns:
         seqnames              ranges strand
            <Rle>           <IRanges>  <Rle>
     [1]    chr15   17435126-17441610      *
     [2]    chr15   17457775-17460493      *
     [3]    chr15   18465420-18523463      *
     [4]    chr15   18754402-18770002      *
     [5]    chr15   22491925-22497810      *
     ...      ...                 ...    ...
  [2419]    chr15 101292364-101299094      *
  [2420]    chr15 101292364-101299094      *
  [2421]    chr15 101307065-101336013      *
  [2422]    chr15 101389337-101395268      *
  [2423]    chr15 101648802-101656417      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths"


anchor.y <- GRanges(seqnames = GB_specific$chr.y, IRanges(start = GB_specific$start.y, end =GB_specific$end.y))

anchor.y

"GRanges object with 2423 ranges and 0 metadata columns:
         seqnames              ranges strand
            <Rle>           <IRanges>  <Rle>
     [1]    chr15   18230127-18231848      *
     [2]    chr15   18230127-18231848      *
     [3]    chr15   19167038-19170276      *
     [4]    chr15   19070204-19108832      *
     [5]    chr15   22836947-22842739      *
     ...      ...                 ...    ...
  [2419]    chr15 101722364-101725660      *
  [2420]    chr15 101953134-101970476      *
  [2421]    chr15 101648802-101656417      *
  [2422]    chr15 101608901-101621083      *
  [2423]    chr15 101722364-101725660      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths"

bed_file <- GRangesList(anchor.x, anchor.y)

bed_file <- unlist(bed_file)

bed_file

"GRanges object with 4846 ranges and 0 metadata columns:
         seqnames              ranges strand
            <Rle>           <IRanges>  <Rle>
     [1]    chr15   17435126-17441610      *
     [2]    chr15   17457775-17460493      *
     [3]    chr15   18465420-18523463      *
     [4]    chr15   18754402-18770002      *
     [5]    chr15   22491925-22497810      *
     ...      ...                 ...    ...
  [4842]    chr15 101722364-101725660      *
  [4843]    chr15 101953134-101970476      *
  [4844]    chr15 101648802-101656417      *
  [4845]    chr15 101608901-101621083      *
  [4846]    chr15 101722364-101725660      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths"


bed_unique <- unique(bed_file)

bed_unique

"GRanges object with 789 ranges and 0 metadata columns:
        seqnames              ranges strand
           <Rle>           <IRanges>  <Rle>
    [1]    chr15   17435126-17441610      *
    [2]    chr15   17457775-17460493      *
    [3]    chr15   18465420-18523463      *
    [4]    chr15   18754402-18770002      *
    [5]    chr15   22491925-22497810      *
    ...      ...                 ...    ...
  [785]    chr15 101608901-101621083      *
  [786]    chr15 101722364-101725660      *
  [787]    chr15 101751634-101768881      *
  [788]    chr15 101595158-101600078      *
  [789]    chr15 101953134-101970476      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths"


export.bed(bed_unique,"/Volumes/ucmm/GrpSR/ChaitaliChakraborty/conferences/PALS2024/Gb.specific.anchors.bed")




