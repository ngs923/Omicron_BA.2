library(ggtree)
library(ggplot2)
library(ape)
library(castor)
library(ggnewscale)
library(gridExtra)
library(phytools)
require(cowplot)
library(RColorBrewer)
library(dplyr)


## Those input files will be available upon request.
folder <- "Omicron624_run_treefile"
info_file <- "Omicron624.metadata.tsv"
gene_loc_file <- "Ben637.snp.txt.snp_pos.txt"

##info
info <- read.delim(info_file, header = T, sep = "\t") ;dim(info)
names(info)
year_month <- length(unique(info$year_month))
rownames(info) <- info$gisaid_epi_isl
dim(info)
names(info)
info[is.na(info$year),"year"] <- info[is.na(info$year),"date"]
info[is.na(info$year_month),"year_month"] <- info[is.na(info$year_month),"date"]

##tree
all.trees <- list.files(folder)

##colour
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
cols <- col_vector

##mutation location
read.delim(gene_loc_file) -> snp
head(snp)

##VOC
voc.info <- read.delim("lineages.voc.txt")
unique(voc.info$VOC)
row.names(voc.info) <- voc.info$Pango
info$voc <- voc.info[info$pango_lineage,2]

names(info)
info.lineage <- data.frame(info[,c(15)])
rownames(info.lineage) <- info$gisaid_epi_isl
colnames(info.lineage) <- "Lineage"
head(info.lineage)

lineage_colour <- read.delim("lineage.hex_color.txt", sep = "\t")
lineage_colour
dim(lineage_colour)
cols <- c(lineage_colour$col, cols)


t <- all.trees[5]

for (t in all.trees){
  tree_name <- paste(folder,"/",t,sep = "")
  tree <- read.tree(tree_name)

  tips <- tree$tip.label
  
  if ("EPI_ISL_402125" %in% tips){
    tree <- root(tree,"EPI_ISL_402125")
  }
  
  adjust0 <- -0.6
  size0 <- 0.8
  cut0 <- 0.5
  

  info.sub <- info[tree$tip.label,]
  write.table(info.sub, file = "Omicron624.metadata.tsv", quote = F, sep = "\t", row.names = F)
  value.sub <- unique(info.sub$year_month)
  Omicron <- c("BA.2","BA.3","BA.1.1","BA.1")
  list0 <- rownames(info.sub[info.sub$voc %in% Omicron,])
  Omicron_anc <- findMRCA(tree, list0)
  Omicron_anc
  off_spring <- tree$edge[tree$edge[,1]==Omicron_anc,]
  
  tree$edge[Omicron_anc,]

  name0 <- rownames(lineage_colour)

  
  ##genome
  head(snp)
  anno0 <- setdiff(unique(snp$gene),"non_coding")
  name1 <- c(name0,anno0)
  name1
  
  shape0 <- "."
  
  data1 <- as.matrix(info.sub)[,"voc"]
  
  td0 <- data.frame(node = nodeid(tree, names(data1)),
                   Lineage = data1)
  d0 <- td0
  
  tree <- full_join(tree, d0, by = 'node')

  size0 <- 0.5 
  dot.size <- 0.8
  p <- ggtree(tree,aes(color=Lineage), size = size0) + theme_tree2() +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 79), size = dot.size, color = "red", shape = 20) +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50 & as.numeric(label) < 80), size = dot.size, color = "blue", shape = 20)
  
  p1 <- p %<+% info.sub +
    scale_color_manual(breaks=name0, values=cols, name="Lineages&Genes")
  p2 <- p1 + geom_facet(panel = "genome", data = snp, 
                        geom = geom_point, 
                        mapping=aes(x = pos, color = gene), 
                        shape = shape0)+ 
    scale_color_manual(breaks=name1, values=cols, name="Lineages&Genes")+
    ggtitle(t)+
    theme(legend.position = 'right',
          legend.background = element_rect(),
          legend.key = element_blank(), # removes the border
          legend.key.size = unit(0.5, 'cm'), # sets overall area/size of the legend
          legend.text = element_text(size = 8), # text size
          title = element_text(size = 10))

  output_file <- paste(t,".pdf",sep="")
  pdf(output_file,width=11.69,height=8.27)
  print(p1)
  print(p2)
  
  #S
  gene <- "S"
  gene
  snp.sub <- snp[snp$gene == gene,];dim(snp.sub);head(snp.sub)
  
  p <- ggtree(tree,aes(color=Lineage), size = size0) + theme_tree2() +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 79), size = dot.size, color = "red", shape = 20) +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50 & as.numeric(label) < 80), size = dot.size, color = "blue", shape = 20)
  p1 <- p %<+% info.sub +
    scale_color_manual(breaks=name0, values=cols, name="Lineages&Genes")
  p3 <- p1 + geom_facet(panel = gene, data = snp.sub, 
                        geom = geom_point, 
                        mapping=aes(x = pos, color = voc), 
                        shape = shape0) +
    scale_fill_manual(breaks=name0, 
                      values=cols, name="Lineages&Genes")
  
  print(p3)
  
  ##
  gene <- "ORF1ab"
  gene
  snp.sub <- snp[snp$gene == gene,];dim(snp.sub);head(snp.sub)
  
  p <- ggtree(tree,aes(color=Lineage), size = size0) + theme_tree2() +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 79), size = dot.size, color = "red", shape = 20) +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50 & as.numeric(label) < 80), size = dot.size, color = "blue", shape = 20)
  
  p1 <- p %<+% info.sub +
    scale_color_manual(breaks=name0, values=cols, name="Lineages&Genes")
  p4 <- p1 + geom_facet(panel = gene, data = snp.sub, 
                        geom = geom_point, 
                        mapping=aes(x = pos, color = voc), 
                        shape = shape0) +
    scale_fill_manual(breaks=name0, 
                      values=cols, name="Lineages&Genes")
  
  print(p4)
  
  ##
  gene <- "S_RBD"
  gene
  snp.sub <- snp[snp$gene == gene,];dim(snp.sub);head(snp.sub)
  
  p <- ggtree(tree,aes(color=Lineage), size = size0) + theme_tree2() +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 79), size = dot.size, color = "red", shape = 20) +
    geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50 & as.numeric(label) < 80), size = dot.size, color = "blue", shape = 20)
  
  p1 <- p %<+% info.sub +
    scale_color_manual(breaks=name0, values=cols, name="Lineages&Genes")
  p5 <- p1 + geom_facet(panel = gene, data = snp.sub, 
                        geom = geom_point, 
                        mapping=aes(x = pos, color = voc), 
                        shape = shape0) +
    scale_fill_manual(breaks=name0, 
                      values=cols, name="Lineages&Genes")
  
  print(p5)
  
  dev.off()
}
