## SpaceRanger

```sh
#!/bin/bash

#=====0. General Information=====
# GEO Accession: GSE208654
# Purpose: Download data, process samples, and run spaceranger pipeline
# Author: PengXia

#=====1. Configuration=====
# Define GEO sample IDs
cat > config <<EOF
SRR20330028
SRR20330029
SRR20330030
SRR20330031
EOF

#=====2. Software Download=====
# Download spaceranger software
curl -o spaceranger-3.0.1.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.0.1.tar.gz?Expires=1721576115&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=io~5pxarBXu0cnJnFte7zAN00SmQqTzIEfz7KZ7QxKaGJvs27b5PyIsWiV1sPFLxMgysv3tRhKiD0wHFfrHgEbtRYaHFSSx3nZd~aD8oCd4~9rMCjIFzIMMuW--R1ZrdU5WLwx2~FJdrYq9RnMzEFqNqzIerFcT4X3z21Q9Y0uYAWX7YBueMg8Effrl5Gn6gCS99VvB4~kZA5F24baLYgJ7JuAWgoUVuh8CTlTkqbSoC2dtmzZbuQrj7FUAy8~Hf7AM95dQXYerQ3-6RdHFUk7tZD6zJzHYaC9TjrwKDzNnx6uVPz7YL1GgCmUyXQAsqLosPKpYVNKZzrhGs67SQlw__"

#=====3. Reference Genome Download=====
# Download reference genome
curl -O "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"
wget "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"

#=====4. Data Download=====
# Create directories for data
mkdir -p 1.sra 2.output figure

# Download SRA files using prefetch
nohup cat config | while read id; do
    prefetch $id \
    --max-size 100G \
    -O 1.sra;
done > download.log 2>&1 &

#=====5. Data Download Using Kingfisher=====
# Create and activate a conda environment for Kingfisher
conda create -n kingfisher -y
conda activate kingfisher
conda install mamba -y
mamba install -c conda-forge -c bioconda kingfisher -y

# Download SRA files using Kingfisher
nohup cat config | while read id; do
    kingfisher get -r $id \
    -m aws-http prefetch aws-cp gcp-cp ena-ascp ena-ftp \
    -t 35 > download.log 2>&1 &
done

# Extract SRA files to FASTQ format
nohup kingfisher extract --sra SRR20330028.sra -t 16 -f fastq > trans.log 2>&1 &
nohup kingfisher extract --sra SRR20330029.sra -t 16 -f fastq > trans.log 2>&1 &
nohup kingfisher extract --sra SRR20330030.sra -t 16 -f fastq > trans.log 2>&1 &
nohup kingfisher extract --sra SRR20330031.sra -t 16 -f fastq > trans.log 2>&1 &

#=====6. Spaceranger Pipeline=====
# Define sample names and corresponding image files
sample_names=("SRR20330028" "SRR20330029" "SRR20330030" "SRR20330031")
image_files=("GSM6360689_N.jpg" "GSM6360690_HSIL.jpg" "GSM6360691_SCC.jpg" "GSM6360692_ADC.jpg")

# Loop through each sample and process
for ((i=0; i<${#sample_names[@]}; i++)); do
    # Rename FASTQ files
    mv ${sample_names[i]}_1.fastq GSM6360689_ST1_S1_L001_R1_001.fastq
    mv ${sample_names[i]}_2.fastq GSM6360689_ST1_S1_L001_R2_001.fastq

    # Run spaceranger count
    nohup spaceranger count --id=${sample_names[i]} \
                            --description=${sample_names[i]} \
                            --transcriptome=/home/somnus/database/refdata-gex-GRCh38-2020-A \
                            --fastqs=/home/somnus/project \
                            --image=/home/somnus/project/${image_files[i]} \
                            --unknown-slide visium-1 \
                            --create-bam false \
                            --localcores=20 \
                            --localmem=200 \
                            > spaceranger_${sample_names[i]}.log 2>&1 &
done

```



## Config

```R
pal <- c( "#F1788D",  "#54990F","#E6550D","#843C39", "#3182BD","#8C6D31",
                  "#BD9E39", "#E7BA52", "#31A354", "#E41A1C", "#6BAED6",
                  "#9ECAE1", "#AD494A", "#A1D99B", "#C7E9C0", "#99600F",
                  "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B", 
                  "#66A61E","#E6550D", "#E7969C",'#53A85F')
                  
pal1 <- c('#EFD4B4','#9BCFB9','#F0C0BD','#FB6264','#4E9C99','#70C055','#E98A27','#FEBC28',
                   "#F1788D",  "#54990F","#E6550D","#843C39","#66A61E","#E6550D", "#E7969C",'#53A85F')
                   


pal2 <-  c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")
library(data.table)
library(Seurat)
library(scRNAtoolVis)
library(qs)
library(export)
library(harmony)
library(future)
library(future.apply)
library(viridis)
library(paletteer)
options(future.globals.maxSize = 8000 * 1024^2)


run_normalize <- function(seurat_obj,dims = 1:30){
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  scale_genes <- VariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = scale_genes)
  seurat_obj <- RunPCA(seurat_obj, features = scale_genes)
  scRNA_harmony <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")
  seurat_obj <- FindNeighbors(scRNA_harmony, dims = dims, reduction = "harmony")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction = "harmony") %>% 
    RunTSNE(., dims = dims, reduction = "harmony")
}



feature_plot <- function(data, feature, reduction = "umap", pt.size = 0.0001, max.cutoff = 1.5, cols=pal2, title) {
  plot <- FeaturePlot(
    object = data,
    features = feature,
    reduction = reduction,
    pt.size = pt.size,
    max.cutoff = max.cutoff,
    cols = cols
  ) +
    scale_x_continuous("") +
    scale_y_continuous("") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    ggtitle(title)
  
  return(plot)
}



run_normalize <- function(seurat_obj,dims = 1:30){
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  scale_genes <- VariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = scale_genes)
  seurat_obj <- RunPCA(seurat_obj, features = scale_genes)
  scRNA_harmony <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")
  seurat_obj <- FindNeighbors(scRNA_harmony, dims = dims, reduction = "harmony")
  seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction = "harmony")
}




plot_umap <- function(scdata, reduction = "umap", group.by = 'celltype', label = TRUE, repel = TRUE, pt.size = 0.01, 
                      cols = c('#6C67AC', '#FDBF6F', '#B49D99', '#9DCAE1', '#E31A1C', '#FF7F00', '#DEEDF9', '#3F93C7', '#CAB2D6')) {
  # 绘制 UMAP 图
  Seurat::DimPlot(scdata, reduction = reduction, group.by = group.by, label = label, repel = repel, pt.size = pt.size, cols = cols)
}





create_spatial_plot <- function(RDS, feature = NULL, ad = 1, limits = 1, min.cutoff = 0, 
                                plot_type = c("feature", "spatial"), title = "") {
  # 生成颜色调色板
  viridis_plasma_light_high <- as.vector(x = paletteer_c(palette = "viridis::inferno", n = 250, direction = 1))
  viridis_plasma_light_high <- c(rep("black", ad), viridis_plasma_light_high)
  
  # 设置隐藏轴的主题
  hide_axis <- theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
  
  plot_type <- match.arg(plot_type)
  
  if (plot_type == "feature" & !is.null(feature)) {
    # 生成特征图
    p <- FeaturePlot(RDS, features = feature, reduction = 'spatial', min.cutoff = min.cutoff)
    p <- p + theme_void() + theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "white", fill = NA, size = 2)
    ) + scale_y_reverse() +
      ggtitle(paste0(feature,title)) +
      DarkTheme() +
      xlab(NULL) + 
      ylab(NULL)
    
    if (length(limits) == 1) {
      p <- p + scale_colour_gradientn(colours = viridis_plasma_light_high, na.value = "black")
    } else {
      p <- p + scale_colour_gradientn(colours = viridis_plasma_light_high, na.value = "black", limits = limits)
    }
    
    return(p)
    
  } else if (plot_type == "spatial") {
    # 生成空间图
    p <- SpatialPlot(RDS, repel = FALSE, label = FALSE, image.alpha = 1, alpha = c(0, 0), pt.size.factor = 0.000001) + 
      geom_point(alpha = 0) +
      NoLegend() +
      DarkTheme() +
      hide_axis + 
      ggtitle(title) +  
      theme(text = element_text(size = 14, face = "bold"))
    
    return(p)
  } else {
    stop("For 'feature' plot type, 'feature' parameter must be provided.")
  }
}



### define function
Rcpp::sourceCpp(code='
                
                #include <Rcpp.h>
                
                using namespace Rcpp;
                
                // [[Rcpp::export]]
                
                IntegerMatrix asMatrix(NumericVector rp,
                                       
                                       NumericVector cp,
                                       
                                       NumericVector z,
                                       
                                       int nrows,
                                       
                                       int ncols){
                  
                  int k = z.size() ;
                  
                  IntegerMatrix  mat(nrows, ncols);
                  
                  for (int i = 0; i < k; i++){
                    
                    mat(rp[i],cp[i]) = z[i];
                    
                  }
                  
                  return mat;
                  
                }
                
                ' )

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}



do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


# 加载必要的包
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

# 定义主分析函数
analyze_tissue_dist <- function(meta_data, output_prefix, pdf_width = 8, pdf_height = 4, verbose = 1) {
  
  # 调用 do.tissueDist 函数进行主要分析
  OR_immune_list <- do.tissueDist(cellInfo.tb = meta_data,
                                  out.prefix = sprintf("%s.Immune_cell", output_prefix),
                                  pdf.width = pdf_width, pdf.height = pdf_height, verbose = verbose)
  
  # 返回分析结果
  return(OR_immune_list)
}

# 定义绘图函数
plot_heatmap <- function(OR_list) {
  # 提取 OR 值结果
  a <- OR_list[["OR.dist.tb"]] %>%
    as.data.frame() %>%
    column_to_rownames(var = "rid") %>%
    na.omit()
  
  # 提取 P 值结果
  b <- OR_list$count.dist.melt.ext.tb[, c(1, 2, 6)] %>%
    spread(key = "cid", value = "adj.p.value") %>%
    column_to_rownames(var = "rid")
  
  # 只选择在a中的行
  b <- b[rownames(a),]
  
  # 调整 P 值符号表示
  col <- viridis(11, option = "D")
  b <- ifelse(b >= 0.05 & (a > 1.5 | a < 0.5), "",
              ifelse(b < 0.0001 & (a > 1.5 | a < 0.5), "****",
                     ifelse(b < 0.001 & (a > 1.5 | a < 0.5), "***",
                            ifelse(b < 0.01 & (a > 1.5 | a < 0.5), "**",
                                   ifelse(b < 0.05 & (a > 1.5 | a < 0.5), "*", "")))))
  
  bk <- c(seq(0, 0.99, by = 0.01), seq(1, 2, by = 0.01))
  
  # 绘制热图
  pheatmap(a, border_color = NA, fontsize = 9, cellheight = 12, cellwidth = 20,
           clustering_distance_rows = "correlation", display_numbers = b,
           number_color = "black", fontsize_number = 10, cluster_col = FALSE,
           cluster_rows = TRUE, breaks = bk, treeheight_row = 20, treeheight_col = 20,
           color = c(colorRampPalette(colors = col[1:6])(length(bk) / 2),
                     colorRampPalette(colors = col[6:11])(length(bk) / 2)))
}



# ======== 绘制monocle2 热图
generate_branch_enrichment_plot <- function(sce, cell_type_col, branch_point, num_clusters = 3, top_n_markers = 50, 
                                            cores = 1, pvalue_cutoff = 0.05, topn_enrich = 8, seed = 5201314, 
                                            num_mark_genes = 25, pdf_path = './output_figure/Figure3/branch-enrich.pdf', 
                                            pdf_height = 9, pdf_width = 16, plot_type = "both", 
                                            column_names_rot = 45, show_row_dend = FALSE, 
                                            markGenes_side = "left", go_colors = jjAnno::useMyCol("calm", n = 3)) {
  
  library(org.Hs.eg.db)
  library(ClusterGVis)
  library(dplyr)
  library(scutilsR)
  library(monocle)
  library(ggplot2)
  
  # 设置细胞类型
  Idents(sce) <- cell_type_col
  
  # 找到所有标记基因
  cell_marker <- scutilsR::mcFindAllMarkers(sce)
  
  # 选择每个cluster的top n个标记基因
  top <- cell_marker %>% 
    group_by(cluster) %>% 
    top_n(n = top_n_markers, wt = avg_log2FC)
  
  # 生成分支热图数据
  df <- plot_genes_branched_heatmap2(mycds[unique(top$Gene.name.uniq),],
                                     branch_point = branch_point,
                                     num_clusters = num_clusters,
                                     cores = cores,
                                     use_gene_short_name = TRUE,
                                     show_rownames = TRUE)
  
  # 富集分析
  enrich <- enrichCluster(object = df, OrgDb = org.Hs.eg.db, 
                          type = "BP", organism = "hsa", 
                          pvalueCutoff = pvalue_cutoff, topn = topn_enrich, 
                          seed = seed)
  
  # 随机选择标记基因
  markGenes <- sample(unique(df$wide.res$gene), num_mark_genes, replace = FALSE)
  
  # 绘制并保存PDF
  pdf(pdf_path, height = pdf_height, width = pdf_width, onefile = FALSE)
  visCluster(object = df, plot.type = plot_type, column_names_rot = column_names_rot, 
             show_row_dend = show_row_dend, markGenes = markGenes, 
             markGenes.side = markGenes_side, annoTerm.data = enrich, 
             go.col = c(rep(go_colors, each = topn_enrich)), 
             add.bar = TRUE, line.side = markGenes_side)
  dev.off()
  
  # 调用函数示例
  # generate_branch_enrichment_plot(sce = sce, cell_type_col = 'cell_type', branch_point = 1)
  
}


# ======== gsea
run_gsea_analysis <- function(sce, ident1, ident2, logfc_threshold = 0.1, 
                              gmt_paths = list("./input_data/h.all.v7.4.symbols.gmt", 
                                               "./input_data/c5.go.bp.v7.4.symbols.gmt", 
                                               "./input_data/c2.cp.kegg.v7.4.symbols.gmt"), 
                              pvalue_cutoff = 0.05, p_adjust_method = 'none') {
  
  # 设置细胞类型
  Idents(sce) <- ident1
  
  # 找到差异表达基因
  sce_edg <- FindMarkers(sce, ident.1 = ident1, ident.2 = ident2, logfc.threshold = logfc_threshold)
  
  # 添加基因符号
  sce_edg$SYMBOL <- rownames(sce_edg)
  names(sce_edg)[2] <- "logFC"
  sce_edg <- sce_edg %>% arrange(desc(logFC))
  
  # 创建基因列表
  geneList <- sce_edg$logFC 
  names(geneList) <- sce_edg$SYMBOL 
  
  # 读取GMT文件
  gmt_list <- lapply(gmt_paths, read.gmt)
  
  # 运行GSEA分析
  gsea_results <- lapply(gmt_list, function(gmt) {
    GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_adjust_method)
  })
  
  # 获取GSEA结果
  gsea_results_list <- lapply(gsea_results, function(gsea) gsea@result)
  
  return(gsea_results_list)
}
```



## Basic Analysis

```R
# Code by PengXia
# ================= 1. Create Directories =================
# Create output directories for each figure
for (i in 1:12) {
  dir.create(paste0('./output_data/Figure', i))
  dir.create(paste0('./output_figure/Figure', i))
}

# ================= 2. Data Loading =================
# Load configuration file
source('./script/config.R')

# Get the list of files in the input data directory
files <- list.files('./input_data/GSE208653', pattern = '^GSM.*', full.names = TRUE)
files

# Define sample groups
group <- data.frame(
  sample_id = c("GSM6360680", "GSM6360681", "GSM6360682", 
                "GSM6360683", "GSM6360684", 
                "GSM6360685", "GSM6360686", "GSM6360687", "GSM6360688"),
  condition = c("NO_HPV", "NO_HPV", "N_HPV", "N_HPV", "HSIL_HPV",
                "HSIL_HPV", "CA_HPV", "CA_HPV", "CA_HPV")
)

# Read and create Seurat objects for each dataset
sce <- mclapply(1:length(files), function(x) {
  data <- Read10X(data.dir = files[x])
  CreateSeuratObject(counts = data, project = group$sample_id[x],
                     min.cells = 3, min.features = 250)
}, mc.cores = length(files))

# Merge all Seurat objects
sce <- Reduce(merge, sce)

# Add group information to Seurat object
sce$group <- group$condition[match(sce$orig.ident, group$sample_id)]

# Save raw Seurat object
qs::qsave(sce, file = './output_data/Figure1/scRNA_raw.qs')

# ================= 3. Dimensionality Reduction and Clustering =================
rm(list = ls())
source('./script/config.R')

# Load raw Seurat object
sce <- qs::qread('./output_data/scRNA/scRNA_raw.qs')

# Calculate mitochondrial gene percentage
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

# Plot quality control metrics
p1 <- VlnPlot(sce, assay = 'RNA', features = c("nCount_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
p2 <- VlnPlot(sce, assay = 'RNA', features = c("nFeature_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
p3 <- VlnPlot(sce, assay = 'RNA', features = c("percent.mt"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
p1 + p2 + p3

# Filter cells based on QC metrics
minGene <- 200
maxGene <- 5000
pctMT <- 15
maxCounts <- 20000

scRNA_filtered <- subset(sce, subset = nFeature_RNA > minGene & 
                           nFeature_RNA < maxGene & 
                           percent.mt < pctMT & nCount_RNA < maxCounts)

# Plot filtered QC metrics
p1 <- VlnPlot(scRNA_filtered, assay = 'RNA', features = c("nCount_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
p2 <- VlnPlot(scRNA_filtered, assay = 'RNA', features = c("nFeature_RNA"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
p3 <- VlnPlot(scRNA_filtered, assay = 'RNA', features = c("percent.mt"), cols = pal, pt.size = 0, ncol = 1) + NoLegend()
p1 + p2 + p3

# Save filtered QC plots
output_fig_dir <- './output_figure/scRNA/'
graph2pdf(p1, file = paste0(output_fig_dir, 'AfterQC-nCount_RNA.pdf'), width = 6, height = 4)
graph2pdf(p2, file = paste0(output_fig_dir, 'AfterQC-nFeature.pdf'), width = 6, height = 4)
graph2pdf(p3, file = paste0(output_fig_dir, 'AfterQC-percent-mt.pdf'), width = 6, height = 4)

# Normalize and perform dimensionality reduction
options(future.globals.maxSize = 8000 * 1024^2)
plan(multisession, workers = 10)

run_normalize <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  scale_genes <- VariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = scale_genes)
  seurat_obj <- RunPCA(seurat_obj, features = scale_genes)
  scRNA_harmony <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")
  seurat_obj <- FindNeighbors(scRNA_harmony, dims = 1:30, reduction = "harmony")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony") %>% 
    RunTSNE(., dims = 1:30, reduction = "harmony")
}

scRNA_normalized <- run_normalize(scRNA_filtered)
plan(sequential)

# Visualize PCA and Harmony
DimPlot(scRNA_normalized, reduction = 'pca')
DimPlot(scRNA_normalized, reduction = 'harmony', group.by = 'orig.ident')

# Determine clustering resolution and run clustering
resolutions <- c(0.4, 0.6, 0.8, 1, 1.2)
scRNA_test <- FindClusters(scRNA_normalized, resolution = resolutions)
index <- colnames(scRNA_test@meta.data)[str_detect(colnames(scRNA_test@meta.data), pattern = '^RNA_snn.*')]
DimPlot(scRNA_test, reduction = 'tsne', group.by = paste0("RNA_snn_res.", resolutions))
DimPlot(scRNA_test, reduction = 'umap', group.by = 'orig.ident', label = TRUE, repel = TRUE)

# Save final clustered object
scRNA_normalized <- FindClusters(scRNA_normalized, resolution = 1)
rm(scRNA_test)
qsave(scRNA_normalized, file = './output_data/Figure1/scRNA_merge_reduction.qs')

# ================= 4. Cell Annotation =================
rm(list = ls())
source('./script/config.R')

# Load clustered object
sce <- qread('./output_data/Figure1/scRNA_merge_reduction.qs')

# Identify marker genes for all clusters
plan(multisession, workers = 10)
degs <- FindAllMarkers(sce)
save(degs, file = './output_data/Figure1/deg_clusters.Rdata')
top20 <- degs %>% group_by(cluster) %>% top_n(20)
plan(sequential)

# Plot UMAP with labels
DimPlot(sce, label = TRUE, repel = TRUE)

# Plot feature markers for specific cell types
FeaturePlot(sce, c("EPCAM", "KLF5", "MKI67")) # Epithelial cells
FeaturePlot(sce, c("DCN", "COL1A1", "COL3A1")) # Fibroblasts
FeaturePlot(sce, c("PECAM1", "CDH5", "VWF")) # Endothelial cells
FeaturePlot(sce, c("ACTA2", "RGS5")) # Smooth muscle cells
FeaturePlot(sce, c("CD68", "CSF1R", "CD163", "LYZ")) # Myeloid cells
FeaturePlot(sce, c("NKG7", "CCL5", "GZMA", "CD3G", "CD3E", "CD3D")) # NK/T cells
FeaturePlot(sce, c("NCF1", "SORL1")) # Neutrophils
FeaturePlot(sce, c("CD19", "BANK1", "MS4A1")) # B cells
FeaturePlot(sce, c("TPSAB1", "CPA3")) # Mast cells
FeaturePlot(sce, c("MZB1", "IGHG1", "IGKC", "IGHG3", "XBP1", "JCHAIN")) # Plasma cells

# Annotate clusters with cell types
sce_Anno <- RenameIdents(
  sce,
  '0' = 'NK/T',
  '1' = 'Neutrophil',
  '2' = 'NK/T',
  '3' = 'NK/T',
  '4' = 'Epithelial',
  '5' = 'Epithelial',
  '6' = 'NK/T',
  '7' = 'Fibroblast',
  '8' = 'Myeloid',
  '9' = 'NK/T',
  '10' = 'Plasma cell',
  '11' = 'Epithelial',
  '12' = 'Mast cell',
  '13' = 'Epithelial',
  '14' = 'Endothelial',
  '15' = 'Myeloid',
  '16' = 'NK/T',
  '17' = 'B cell',
  '18' = 'Epithelial',
  '19' = 'Neutrophil',
  '20' = 'Epithelial',
  '21' = 'Smooth muscle cells',
  '22' = 'Epithelial',
  '23' = 'Epithelial',
  '24' = 'Fibroblast',
  '25' = 'Epithelial',
  '26' = 'Epithelial',
  '27' = 'Epithelial',
  '28' = 'Epithelial',
  '29' = 'Epithelial'
)

# Assign cell type to metadata
sce_Anno$cell_type <- Idents(sce_Anno)

# Plot UMAP with cell type annotations
DimPlot(sce_Anno, group.by = 'seurat_clusters', reduction = "umap", label = TRUE)
DimPlot(sce_Anno, group.by = 'cell_type', reduction = "umap", label = TRUE)

# Save annotated Seurat object
qsave(sce_Anno, file = './output_data/Figure1/scRNA_anno.qs')
```

## InferCNV

```R
# Code by PengXia
# Set working directory and load configuration
setwd(here::here())
source('./script/config.R')
library(infercnv)
library(parallel)

# ================== 1. Input ==================
# Load annotated Seurat object
sce <- qread('./output_data/Figure1/scRNA_anno.qs')

# Create output directory for inferCNV analysis
dir.create('./output_data/inferCNV', showWarnings = FALSE)
path <- './output_data/inferCNV/'

# ================== 2. Split by Sample ==================
# Check cell type and group distribution
table(sce$cell_type, sce$group)

# Subset epithelial cells
sce.epi <- subset(sce, cell_type == 'epithelial')

# Separate reference (normal) and tumor groups
sce.refer <- subset(sce.epi, group == 'NO_HPV')
sce.tumor <- subset(sce.epi, group == 'CA_HPV')

# Split tumor samples by sample ID
sce.split <- SplitObject(sce.tumor, split.by = 'orig.ident')

# Merge reference cells into each tumor sample
name <- names(sce.split)
sce.split <- lapply(1:length(sce.split), function(i) {
  sce.new <- merge(sce.split[[i]], sce.refer)
  sce.new$cell_type <- paste0(sce.new$cell_type, '_', sce.new$group)
  return(sce.new)
})
names(sce.split) <- name

# Check cell type and sample distribution
table(sce.split[[1]]$cell_type)
table(sce.split[[1]]$orig.ident)

# ================== 3. Generate Matrix and Cell Labels ==================
# Save count matrices and cell type labels
mclapply(1:length(sce.split), function(i) {
  qsave(as.matrix(sce.split[[i]][["RNA"]]@counts), file = paste0(path, names(sce.split)[i], '.qs'))
  write.table(sce.split[[i]]$cell_type, 
              file = paste0(path, names(sce.split)[i], '.celltype.label.txt'), 
              sep = "\t", quote = FALSE, col.names = FALSE)
}, mc.cores = length(sce.split))

# ================== 4. Run InferCNV ==================
# Define input parameters
gene_order_file <- './output_data/inferCNV/hg38_gencode_v27.txt'
ref_group_names <- c("epithelial_NO_HPV")  # Reference group (normal cells)
name <- 'infer_run'

# Run inferCNV for each sample
lapply(1:length(sce.split), function(i) {
  matrix_counts <- qread(paste0(path, names(sce.split)[i], '.qs'))  # Counts matrix
  annotations_file <- paste0(path, names(sce.split)[i], '.celltype.label.txt')  # Cell type labels
  out_path <- paste0(path, name, "_", names(sce.split)[i])  # Output directory
  
  # Create inferCNV object
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = matrix_counts,
    annotations_file = annotations_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = ref_group_names,
    chr_exclude = c("chrY", "chrM")  # Exclude unwanted chromosomes
  )
  
  # Run inferCNV
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,  # Suitable for 10x Genomics
    out_dir = out_path,
    no_prelim_plot = TRUE,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    min_cells_per_gene = 10,
    num_threads = 1,
    write_expr_matrix = TRUE
  )
})

# ================== 5. Visualization ==================
gc()
library(ComplexHeatmap)

# Define input files
infercnv_res <- paste0('./output_data/inferCNV/', names(sce.split), '/run.final.infercnv_obj')
gene_order_file <- "./output_data/inferCNV/hg38_gencode_v27.txt"
ref_cell_name <- c("epithelial_NO_HPV")
obs_cell_name <- c("epithelial_CA_HPV")
ref_group <- 'NO_HPV'
obs_group <- 'CA_HPV'

# Generate heatmaps for each sample
mclapply(1:length(infercnv_res), function(i) {
  # Load inferCNV results
  infercnv_obj <- read_rds(infercnv_res[i])
  expr <- infercnv_obj@expr.data
  
  # Load chromosome annotation
  gene_pos <- read.delim(gene_order_file, header = FALSE)
  gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr), ]
  new_cluster <- unique(gene_pos$V2)
  
  # Top annotation for chromosomes
  top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                      gp = gpar(col = "white"),
                                                      labels_gp = gpar(cex = 1, col = "black"),
                                                      height = unit(5, "mm")))
  
  # Extract reference and observation cell indices
  ref_cell <- c(infercnv_obj@reference_grouped_cell_indices[ref_cell_name][[1]])
  obs_cell <- c(infercnv_obj@observation_grouped_cell_indices[obs_cell_name][[1]])
  cell_anno <- data.frame(
    cell_id = c(colnames(expr)[ref_cell], colnames(expr)[obs_cell]),
    group = c(rep(ref_group, length(ref_cell)), rep(obs_group, length(obs_cell)))
  )
  
  # Perform k-means clustering
  set.seed(123)
  kmeans.result <- kmeans(t(expr), 5)
  kmeans_df <- data.frame(kmeans.result$cluster)
  colnames(kmeans_df) <- "k_cluster"
  kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
  kmeans_df <- kmeans_df %>% inner_join(cell_anno, by = "cell_id") %>% arrange(k_cluster)
  kmeans_df$k_cluster <- as.factor(kmeans_df$k_cluster)
  
  # Row annotation
  annotation_row <- data.frame(
    k_cluster = kmeans_df$k_cluster,
    group = kmeans_df$group
  )
  row.names(annotation_row) <- kmeans_df$cell_id
  saveRDS(annotation_row, file = paste0('./output_data/inferCNV/', names(sce.split)[i], '/Kmeans.rds'))
  
  # Define colors for clusters
  color_cluster <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
  names(color_cluster) <- as.character(1:5)
  
  # Left annotation
  left_anno <- rowAnnotation(
    df = annotation_row,
    col = list(group = c("NO_HPV" = "#00A0877F", "CA_HPV" = "#E64B357F"),
               k_cluster = color_cluster),
    show_annotation_name = FALSE
  )
  
  # Draw heatmap
  pdf(paste0(path, names(sce.split)[i], '/CNV_heatmap.pdf'), width = 9.5, height = 4)
  ht <- Heatmap(
    t(log2(expr))[rownames(annotation_row), ],
    col = colorRamp2::colorRamp2(c(-0.5, 0, 0.5), c("#2166ac", "white", "#b2182b")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_split = factor(gene_pos$V2, new_cluster),
    heatmap_legend_param = list(title = "inferCNV",
                                direction = "vertical",
                                title_position = "leftcenter-rot",
                                legend_height = unit(3, "cm")),
    left_annotation = left_anno,
    row_title = NULL,
    column_title = NULL,
    top_annotation = top_color,
    border = TRUE
  )
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}, mc.cores = 10)

# ================== 6. Identify Tumor Cells ==================
# Define tumor clusters for each sample
sample_tumor_k1 <- c(2, 4)
sample_tumor_k2 <- c(3, 4)
sample_tumor_k3 <- c(3, 4)
k_check <- list(sample_tumor_k1, sample_tumor_k2, sample_tumor_k3)

# Load k-means results and identify tumor cells
files <- paste0('./output_data/inferCNV/', names(sce.split), '/Kmeans.rds')
tumor.table <- lapply(1:length(files), function(i) {
  kmeans_res <- read_rds(files[i])
  tumor.table <- kmeans_res %>% 
    filter(group == 'CA_HPV') %>% 
    mutate(celltype = ifelse(.$k_cluster %in% k_check[[i]], 'tumor_cell', 'normal_cell'))
  message(table(tumor.table$group, tumor.table$celltype))
  return(tumor.table)
})

# Combine results and save
tumor.table <- Reduce(rbind, tumor.table)
table(sce$cell_type, sce$group)
dim(tumor.table)
table(tumor.table$celltype, tumor.table$group)
openxlsx::write.xlsx(rownames_to_column(tumor.table, var = 'id'), file = './output_data/Figure2/tumor.table.xlsx')
```

## scMetabolism

```R
# ========================== Code by PengXia ==========================
# Date: 2024-04-26
# Author: PengXia

# ========================== scMetabolism Analysis ==========================

# Load necessary libraries
library(scMetabolism)
library(Seurat)
library(ggplot2)
library(rsvd)
library(pheatmap)
library(qs)

# ========================== Step 1: Load Data ==========================
# Load single-cell data (same input as scFEA)
human_data <- qread('./output_data/Figure2/scRNA_epi.qs')

# ========================== Step 2: Calculate Metabolic Pathway Scores ==========================
# Compute metabolic scores using KEGG pathways via AUCell (without imputation)
human_countexp_Seurat <- sc.metabolism.Seurat(
  obj = human_data,
  method = "AUCell",
  imputation = FALSE,
  ncores = 10,
  metabolism.type = "KEGG"
)

# Save the Seurat object with metabolism scores
save(human_countexp_Seurat, file = './output_data/Figure4/scMetabolism.Rdata')

# ========================== Step 3: Visualization ==========================
# Generate DotPlot for top 20 KEGG metabolic pathways by score
DotPlot.metabolism(
  obj = human_countexp_Seurat,
  pathway = rownames(human_countexp_Seurat@assays[["METABOLISM"]][["score"]])[1:20],
  phenotype = "cell_type",
  norm = "y"
) + xlab('')

# Save DotPlot as PDF
graph2pdf(
  file = './output_figure/Figure2/scMetabolism.pdf',
  width = 6,
  height = 10
)
```

## scFEA

```R
# Code by PengXia
# Date: 2024-04-26
# Author: PengXia

# ========================== 1. scFEA ==========================
source('./script/config.R')
getwd()

# ========================== Step 1: Prepare scFEA Input Files ==========================
# Create output directory for scFEA
dir.create('./scFEA', showWarnings = FALSE)

# Load single-cell data
sce <- qread('./output_data/Figure2/scRNA_epi.qs')

# Export gene expression data for scFEA
ExportGeneExpr <- sce@assays$RNA@counts
write.csv(ExportGeneExpr, file = './scFEA/input/Seurat_geneExpr.csv', row.names = TRUE)

# Visualize cell types using DimPlot
DimPlot(sce, group.by = 'cell_type')

# ========================== Step 2: Run scFEA (Python) ==========================
# Note: This step assumes scFEA is run externally using Python
# Example command:
# python src/scFEA.py --data_dir data --input_dir input \
# --test_file Seurat_geneExpr.csv \
# --moduleGene_file module_gene_m168.csv \
# --stoichiometry_matrix cmMat_c70_m168.csv \
# --output_flux_file output/Seurat_geneExpr.csv \
# --output_balance_file output/Cervix_flux_balance.csv \
# --sc_imputation True

# ========================== Step 3: Load scFEA Results ==========================
# Load predicted metabolic fluxes
predFlux <- read.csv('./scFEA/output/Cervix_flux.csv', header = TRUE, row.names = 1)
predFlux <- as.data.frame(t(predFlux))

# Load pathway annotations
anno <- read.csv('./scFEA/data/Human_M168_information.symbols.csv')
anno$pathway <- paste0(anno$Compound_IN_name, ' -> ', anno$Compound_OUT_name)

# Ensure consistent row names
if (!identical(anno$X, rownames(predFlux))) {
  stop("Row names of annotations and predicted fluxes do not match.")
}
rownames(predFlux) <- anno$pathway

# Clean up column names
colnames(predFlux) <- gsub('\\.', '-', colnames(predFlux))

# ========================== Step 4: Merge scFEA Results into Seurat ==========================
# Add metabolic fluxes to Seurat object
sce[["FLUX"]] <- CreateAssayObject(counts = predFlux)

# Set default assay to FLUX
DefaultAssay(sce) <- 'FLUX'

# Normalize, identify variable features, and scale data
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)

# Define color palette
colors <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", 
            "#56B4E9", "#E69F00", "#00ADA9", "#D0E429", "#ED008C", "#68217A")

# Find differentially active metabolic pathways between tumor and normal cells
Idents(sce) <- 'cell_type'
df <- FindAllMarkers(sce, only.pos = TRUE, logfc.threshold = 0.5)

# Generate heatmap of average pathway activity
scRNAtoolVis::AverageHeatmap(sce, rownames(df), assays = 'FLUX')

# Save heatmap as PDF
graph2pdf(file = './output_figure/Figure4/Epi_metabolism_scFEA_heatmap.pdf',
          width = 6, height = 7)

# Generate heatmap of top 50 pathways
DoHeatmap(sce, features = rownames(df), size = 5, group.colors = 
            c("#C77CFF", "#7CAE00", "#00BFC4", "#F8766D", "#AB82FF", "#90EE90", 
              "#00CD00", "#008B8B", "#FFA500")) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))

graph2pdf(file = './output_figure/scRNA/Epithelial_metabolism_scFEA.pdf',
          width = 14, height = 7)

# ========================== Step 5: Oxygen Phosphorylation Analysis ==========================
# Define pathways for oxidative phosphorylation
oxp <- c("Pyruvate -> Acetyl-Coa", "Succinyl-CoA -> Succinate",
         "Acetyl-CoA-in -> Acetyl-CoA", "Pyruvate -> Oxaloacetate")

# Generate violin plots for oxidative phosphorylation pathways
p1 <- lapply(1:length(oxp), function(i) {
  tmp <- VlnPlot(sce, features = oxp[i], pt.size = 0, ncol = 1, 
                 cols = colors, stack = FALSE, flip = TRUE) +
    theme(legend.position = "none") +
    ggpubr::stat_compare_means(aes(label = ..p.signif..), label.x = 1.5)
  
  graph2pdf(tmp, file = paste0('./output_figure/Figure4/OXPHOS_Vln_', i, '.pdf'),
            width = 6, height = 5)
  return(tmp)
})

aplot::gglist(p1)

# ========================== Step 6: Glycolysis Analysis ==========================
# Define pathways for glycolysis
gly <- c("Glucose-in -> Glucose", "Glucose -> G6P")

# Generate violin plots for glycolysis pathways
p2 <- lapply(1:length(gly), function(i) {
  tmp <- VlnPlot(sce, features = gly[i], pt.size = 0, ncol = 1, 
                 cols = colors, stack = FALSE, flip = TRUE) +
    theme(legend.position = "none") +
    ggpubr::stat_compare_means(aes(label = ..p.signif..), label.x = 1.5)
  
  graph2pdf(tmp, file = paste0('./output_figure/Figure4/Gly_Vln_', i, '.pdf'),
            width = 6, height = 5)
  return(tmp)
})

aplot::gglist(p2)

# ========================== Step 7: Lipid Metabolism Analysis ==========================
# Define pathways for lipid metabolism
fat <- c("Fatty Acid-in -> Fatty Acid", "Fatty Acid -> Acetyl-CoA",
         "(E,E)-Farnesyl-PP -> Cholesterol")

# Generate violin plots for lipid metabolism pathways
p3 <- lapply(1:length(fat), function(i) {
  tmp <- VlnPlot(sce, features = fat[i], pt.size = 0, ncol = 1, 
                 cols = colors, stack = FALSE, flip = TRUE) +
    theme(legend.position = "none") +
    ggpubr::stat_compare_means(aes(label = ..p.signif..), label.x = 1.5)
  
  graph2pdf(tmp, file = paste0('./output_figure/Figure4/Fatty_Vln_', i, '.pdf'),
            width = 6, height = 5)
  return(tmp)
})

aplot::gglist(p3)

# ========================== Step 8: Save Seurat Object ==========================
# Save updated Seurat object with metabolic fluxes
qs::qsave(sce, file = './output_data/Figure4/scRNA_epi_scFEA.qs')

```

## CellChat

```R
# Code by PengXia
# Date: 2024-04-26
# Author: PengXia
# ========================== Required Libraries ==========================
source('./script/config.R') # Load configuration
library(Seurat)
library(dplyr)
library(CellChat)
library(parallel)
library(circlize)
library(reshape2)
library(ggrepel)
library(patchwork)
library(aplot)

# ========================== 1. Run CellChat Analysis ==========================

# Load scRNA data
sce.all <- qread('./output_data/Figure1/scRNA_anno.qs')

# Visualize UMAP by cell type
plot_umap(sce.all, group.by = 'cell_type', cols = pal)

# Split the dataset by group
sce <- SplitObject(sce.all, split.by = 'group')
name <- names(sce)

# Normalize and cluster each group in parallel
sce <- mclapply(1:length(sce), function(i) {
  sce[[i]] %>%
    run_normalize() %>%
    FindClusters()
}, mc.cores = 4)

# Run CellChat analysis for each group in parallel
cellchat_res <- mclapply(1:length(sce), function(i) {
  data.input <- sce[[i]]@assays$RNA@data
  identity <- sce[[i]]@meta.data
  Idents(sce[[i]]) <- "cell_type"
  
  # Initialize CellChat object
  cellchat <- createCellChat(data.input, meta = identity, group.by = 'cell_type')
  cellchat <- addMeta(cellchat, meta = identity)
  cellchat <- setIdent(cellchat, ident.use = "cell_type")
  
  # Use human CellChat database
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # Subset data and identify overexpressed genes and interactions
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 5)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  # Compute communication probabilities and pathways
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}, mc.cores = 4)

# Save results
names(cellchat_res) <- name
save(sce, cellchat_res, file = './output_data/Figure11/cellchat_raw.Rdata')

# ========================== 2. Visualization ==========================

# Load CellChat results
load('./output_data/Figure11/cellchat_raw.Rdata')

# Extract CellChat objects for specific groups
cellchat_No_HPV <- cellchat_res[[1]]
cellchat_N_HPV <- cellchat_res[[2]]
cellchat_HSIL_HPV <- cellchat_res[[3]]
cellchat_CA_HPV <- cellchat_res[[4]]

# Calculate weight differences between CA_HPV and HSIL_HPV
weight_M <- cellchat_CA_HPV@net$weight %>% data.frame()
weight_NM <- cellchat_HSIL_HPV@net$weight %>% data.frame()

weight <- cellchat_CA_HPV@net$weight
weight2 <- cellchat_HSIL_HPV@net$weight

# Prepare data for chord diagram
weight_diff <- data.frame(weight - weight2)
weight_diff$name <- rownames(weight_diff)

# Reshape and clean data
weight_diff <- reshape2::melt(weight_diff)
weight_diff$name <- gsub(' cell', '', weight_diff$name)
weight_diff$variable <- gsub('.cell', '', as.character(weight_diff$variable))

# Plot upregulated interactions
chordDiagram(filter(weight_diff, value >= 0))
graph2pdf(file = './output_figure/Figure11/cellchat_CA_upregulated.pdf', width = 6, height = 6)

# Plot downregulated interactions
chordDiagram(filter(weight_diff, value <= 0))
graph2pdf(file = './output_figure/Figure11/cellchat_CA_downregulated.pdf', width = 6, height = 6)

# ========================== 3. Ligand Analysis ==========================

# Extract ligands for CA_HPV and HSIL_HPV
ligand_CA_HPV <- cellchat_CA_HPV@LR$LRsig[, 2]
ligand_HSIL_HPV <- cellchat_HSIL_HPV@LR$LRsig[, 2]

# Subset communication data for ligands
ligand_CA_HPV <- subsetCommunication(cellchat_CA_HPV)[which(subsetCommunication(cellchat_CA_HPV)$ligand %in% ligand_CA_HPV), ]
ligand_HSIL_HPV <- subsetCommunication(cellchat_HSIL_HPV)[which(subsetCommunication(cellchat_HSIL_HPV)$ligand %in% ligand_HSIL_HPV), ]

# Summarize ligand probabilities
ligand_CA_HPV <- ligand_CA_HPV %>% group_by(ligand) %>% summarize(prob = median(prob)) %>% data.frame()
ligand_HSIL_HPV <- ligand_HSIL_HPV %>% group_by(ligand) %>% summarize(prob = median(prob)) %>% data.frame()

# Add class labels
ligand_CA_HPV$class <- 'CA_HPV'
ligand_HSIL_HPV$class <- 'HSIL_HPV'

# Combine ligand data
ligand <- rbind(ligand_CA_HPV, ligand_HSIL_HPV)
ligand$y <- factor(rownames(ligand))

# Plot ligand comparison
ggplot(data = ligand, aes(x = y, y = -log(prob), color = as.factor(class), size = -log(prob))) +
  geom_point() +  
  geom_text_repel(label = ligand$ligand, color = 'black') +
  scale_color_manual(values = pal) +
  theme_classic() +
  theme(text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(title = "Ligand Comparison: CA_HPV vs HSIL_HPV", 
       x = "", 
       y = "-log10(communication probability)")

graph2pdf(file = './output_figure/Figure11/cellchat_ligand_comparison.pdf', width = 8, height = 5.5)

# ========================== 4. Merge and Compare ==========================

# Merge CellChat results for comparison
cellchat <- mergeCellChat(cellchat_res, add.names = names(cellchat_res), cell.prefix = TRUE)

# Compute centrality for each group
object.list <- lapply(cellchat_res, function(x) {
  netAnalysis_computeCentrality(x)
})

# Visualize signaling roles for each group
num.link <- sapply(object.list, function(x) { rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count) })
weight.MinMax <- c(min(num.link), max(num.link))
gg <- lapply(1:length(object.list), function(i) {
  netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
})
patchwork::wrap_plots(plots = gg)

graph2pdf(file = './output_figure/Figure11/cellchat_signaling_comparison.pdf', width = 8, height = 5.5)

# ========================== 5. Heatmap Comparison ==========================

# Generate heatmaps for communication comparison
gg1 <- netVisual_heatmap(cellchat, comparison = c(3, 4))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(3, 4))
gg1 + gg2

graph2pdf(file = './output_figure/Figure11/cellchat_heatmap_comparison.pdf', width = 8, height = 4)

# ========================== 6. SPP1 Expression Analysis ==========================

# Plot SPP1 expression by cell type and group
VlnPlot(sce.all, group.by = 'cell_type', features = 'SPP1', pt.size = 0, split.by = 'group')

# Split dataset by group and plot SPP1 expression
sce1 <- SplitObject(sce.all, split.by = 'group')
sce1 <- sce1[c(3, 4)] # Select specific groups
p <- lapply(1:length(sce1), function(x) {
  feature_plot(data = sce1[[x]], feature = 'SPP1', title = paste0('SPP1 ', names(sce1)[x]))
})
aplot::gglist(p)

graph2pdf(file = './output_figure/Figure11/SPP1_expression.pdf', width = 5, height = 2.5)

```

## mistyR

```R
# 加载必要的库和配置文件
library(argparse)
library(tidyverse)
library(Seurat)
library(mistyR)
library(qs)
library(future)
source("./script/stRNA/misty_utilities.R")

# 定义函数以运行共定位分析
run_colocalization <- function(slide, assay, useful_features, out_label, misty_out_alias) {
  # 定义每个视图的assay
  view_assays <- list(
    "main" = assay,
    "juxta" = assay,
    "para" = assay
  )
  
  # 定义每个视图的特征
  view_features <- list(
    "main" = useful_features, 
    "juxta" = useful_features,
    "para" = useful_features
  )
  
  # 定义每个视图的空间上下文
  view_types <- list(
    "main" = "intra", 
    "juxta" = "juxta",
    "para" = "para"
  )
  
  # 定义额外的参数
  view_params <- list(
    "main" = NULL, 
    "juxta" = 2,
    "para" = 5
  )
  
  # 输出路径
  misty_out <- paste0(misty_out_alias, out_label, "_", assay)
  
  # 运行Misty Seurat分析
  run_misty_seurat(
    visium.slide = slide,
    view.assays = view_assays,
    view.features = view_features,
    view.types = view_types,
    view.params = view_params,
    spot.ids = NULL,
    out.alias = misty_out
  )
  
  return(misty_out)
}

# ============1. 加载空间数据和RCTD去卷积的结果============

# 加载单细胞数据
sc <- qs::qread('./output_data/Figure11/sce.tumor.所有细胞注释.qs')

# 加载空间数据
load('./output_data/stRNA/ST_all.Rdata')

spatial_samples = filtered_results[c(3,4)]
# 为空间数据添加坐标信息
for (i in seq_along(filtered_results)) {
  spatial_coords <- filtered_results[[i]]@images$slice1@coordinates[, c('col', 'row')]
  colnames(spatial_coords) <- c('s_1', 's_2')
  spatial_coords <- as.matrix(spatial_coords)
  filtered_results[[i]][["spatial"]] <- CreateDimReducObject(embeddings = spatial_coords, key = "s_", assay = "Spatial")
}

cell2loc1 <- read.csv('./output_data/Figure12/st_cell2location_res_SRR20330030.csv',row.names = 1)
row.names(cell2loc1) <- sapply(strsplit(row.names(cell2loc1),'_'),'[',2)
colnames(cell2loc1) <- sapply(strsplit(colnames(cell2loc1),'w_sf_'),'[',2)

cell2loc2 <- read.csv('./output_data/Figure12/st_cell2location_res_SRR20330031.csv',row.names = 1)
row.names(cell2loc2) <- sapply(strsplit(row.names(cell2loc2),'_'),'[',2)
colnames(cell2loc2) <- sapply(strsplit(colnames(cell2loc2),'w_sf_'),'[',2)

cell2loc <- list(cell2loc1,cell2loc2)
identical(colnames(spatial_samples[[1]]),row.names(cell2loc1))

# 加载RCTD结果
# load('./output_data/stRNA/RCTD_res_all.Rdata')

# ============2. 运行MistyR分析============

# 遍历所有样本（例如HCC1T, HCC2T, HCC3T）
for (i in 1:2) {
  future::plan(future::multisession, workers = 5)
  
  # 获取当前样本的名称
  slide_id <- names(spatial_samples)[i]
  outdir <- './output_data/stRNA/'
  
  # 获取当前样本的RCTD结果和空间样本
  predictions <- cell2loc[[i]]
  
  # 获取空间样本
  slide <- spatial_samples[[i]]
  
  # 检查并处理不匹配的条形码
  unmatched_barcodes <- setdiff(colnames(slide), rownames(predictions))
  if (length(unmatched_barcodes) > 0) {
    slide$barcode <- colnames(slide)
    slide <- subset(slide, barcode %in% setdiff(colnames(slide), unmatched_barcodes))
  }
  
  # 创建新的Assay对象并设置默认Assay
  slide[["predictions"]] <- CreateAssayObject(t(predictions))
  DefaultAssay(slide) <- 'predictions'
  
  # 定义有用的特征
  useful_features <- rownames(slide)
  useful_features <- useful_features[!useful_features %in% "prolif"]
  
  # 运行共定位分析
  misty_output <- run_colocalization(
    slide = slide,
    assay = 'predictions',
    useful_features = useful_features,
    out_label = slide_id,
    misty_out_alias = outdir
  )
  
  # 收集和可视化Misty结果
  misty_results <- collect_results(misty_output)
  plot_folder <- paste0(misty_output, "/plots")
  
  # 创建目录（若不存在）
  if (!dir.exists(plot_folder)) {
    dir.create(plot_folder, recursive = TRUE)
  }
  
  # 创建PDF文件并生成各种图表
  pdf(file = paste0(plot_folder, "/", slide_id, "_summary_plots.pdf"))
  
  mistyR::plot_improvement_stats(misty_results)
  mistyR::plot_view_contributions(misty_results)
  
  mistyR::plot_interaction_heatmap(misty_results, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_results, "intra", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_results, "juxta_2", cutoff = 0)
  mistyR::plot_interaction_communities(misty_results, "juxta_2", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_results, "para_5", cutoff = 0)
  mistyR::plot_interaction_communities(misty_results, "para_5", cutoff = 0)
  
  dev.off()
  future::plan(future::sequential)
}






```

## Cell2loc

```python
import os
os.getcwd()
# 单细胞部分
import scanpy as sc
import scipy.sparse as sp
import pandas as pd
import cell2location
import warnings
warnings.filterwarnings('ignore')
THEANO_FLAGS='force_device=True'
sc_adata = sc.read_h5ad('./scrna_mix.h5ad')
# 输出每个基因在每种细胞类型中的估计表达
if 'means_per_cluster_mu_fg' in sc_adata.varm.keys():
    inf_aver = sc_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in sc_adata.uns['mod']['factor_names']]].copy()
else:
    inf_aver = sc_adata.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in sc_adata.uns['mod']['factor_names']]].copy()
inf_aver.columns = sc_adata.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:13]
from pathlib import Path
import scanpy as sc

adata_dct = {}
for i in Path("./rawdata_st").glob("SRR20330030"):
    _s = i.name  # 直接使用目录名作为样本ID
    _a = sc.read_visium(i, library_id=_s)
    _a.obs.index = [_s + "_" + bc for bc in _a.obs.index.tolist()]
    _a.var_names_make_unique()
    adata_dct[_s] = _a

st_adata = sc.concat(adata_dct, label="sample", uns_merge="unique")
import numpy as np
# 提取共享基因并准备anndata
intersect = np.intersect1d(st_adata.var_names, inf_aver.index)
st_adata = st_adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=st_adata, batch_key="sample")
# 构建训练模型，注意 N_cells_per_location（每个spot的细胞数）和detection_alpha（试验受技术影响从成都）参数
mod = cell2location.models.Cell2location(
    st_adata, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=20,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
# 训练模型
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )
mod.plot_history(200)
# 保存训练数据于空转对象中
st_adata = mod.export_posterior(
    st_adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)
# 质控，强对角线为优
mod.plot_QC()
# 输出结果，每个spot的细胞组成情况
pd.DataFrame(st_adata.obsm['q05_cell_abundance_w_sf']).to_csv("./rawdata_st/st_cell2location_res_SRR20330030.csv")
```







## Commot

```python
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


st_adata  = sc.read_visium('./rawdata_st/SRR20330030/')
st_adata.var_names_make_unique()
sc.pp.normalize_total(st_adata, inplace=True)
sc.pp.log1p(st_adata)
LR=np.array([['SPP1', 'SPP1_CD44', 'SPP1_pathway'],
    ['FN1', 'FN1_CD44', 'FN1_pathway']],dtype=str)
df_ligrec = pd.DataFrame(data=LR)
df_ligrec
ct.tl.spatial_communication(st_adata,
                            database_name='user_database', 
                            df_ligrec=df_ligrec, dis_thr=500,
                            heteromeric=True, pathway_sum=True)

st_adata.obsm['commot-user_database-sum-sender']

#Plot the amount of sent and received signal, for example, of the LR pair Fgf1-Fgfr1.
pts = st_adata.obsm['spatial']
s =  st_adata.obsm['commot-user_database-sum-sender']['s-FN1-FN1_CD44']
r =  st_adata.obsm['commot-user_database-sum-receiver']['r-FN1-FN1_CD44']
fig, ax = plt.subplots(1,2, figsize=(10,4))

# 发送信号的散点图
ax[0].scatter(pts[:,0], pts[:,1], c=s, s=5, cmap='Blues')
ax[0].set_title('Sender')
ax[0].invert_yaxis()  # 上下翻转发送信号图像

# 接收信号的散点图
ax[1].scatter(pts[:,0], pts[:,1], c=r, s=5, cmap='Reds')
ax[1].set_title('Receiver')
ax[1].invert_yaxis()  # 上下翻转接收信号图像

plt.show()
ct.tl.communication_direction(st_adata, database_name='user_database', lr_pair=('FN1','FN1_CD44'), k=5)
# 调用 plot_cell_communication 函数生成图像
ct.pl.plot_cell_communication(st_adata, database_name='user_database', lr_pair=('FN1','FN1_CD44'), plot_method='grid', background_legend=True,
    scale=0.00003, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',
    normalize_v = True, normalize_v_quantile=0.995)

# 获取当前轴对象并翻转 y 轴
plt.gca().invert_yaxis()

# 显示图像
plt.show()
ct.pl.plot_cell_communication(st_adata, database_name='user_database', lr_pair=('FN1','FN1_CD44'), plot_method='grid', background_legend=True,
    scale=0.00003, ndsize=8, grid_density=0.4, summary='receiver', background='summary', clustering='leiden', cmap='Reds',
    normalize_v = True, normalize_v_quantile=0.995)

# 获取当前轴对象并翻转 y 轴
plt.gca().invert_yaxis()

# 显示图像
plt.show()

```

