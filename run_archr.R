library(ArchR)
library(Seurat)
library(patchwork)
library(gridExtra)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Mm.eg.db)
library(repr)
library(purrr)
library(presto)
library(ggseqlogo)
library(chromVARmotifs)


data_species <- 'mm10'
num_threads <- 50
tile_size <- 5000
genomeSize = 3.0e+09
min_TSS <- 3.5 
min_Frags <- 1e+04 
set.seed(42)
#Path to fragments.tsv.gz located in <sample>_cellranger_outs/
inputFile <- "../data/melanoma/D3G5/UCLA_D3G5_ATAC2_fragments.tsv.gz" # path of input fragments file
project_name <- "melanoma_D3G5"
# Path to spatial folder located in <sample>
spatialFolder <- "../data/melanoma/D3G5/UCLA_D3G5_2_spatial"


addArchRGenome(data_species)
geneAnnotation <- getGeneAnnotation()
genomeAnnotation <- getGenomeAnnotation()
addArchRThreads(threads = num_threads)

ArrowFiles <- createArrowFiles(
   inputFiles = inputFile,
   sampleNames = project_name,
   geneAnnotation = geneAnnotation,
   genomeAnnotation = genomeAnnotation,
   minTSS = min_TSS,
   minFrags = min_Frags,
   maxFrags = 1e+07,
   addTileMat = TRUE,
   addGeneScoreMat = TRUE,
   offsetPlus = 0,
   offsetMinus = 0,
   force = TRUE,
   TileMatParams = list(tileSize = tile_size)
)


proj <- ArchRProject(
   ArrowFiles = ArrowFiles,
   outputDirectory = project_name,
   geneAnnotation = geneAnnotation,
   genomeAnnotation = genomeAnnotation,
   copyArrows = TRUE
)


############### Prepare meta.data
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

############### Filtering off-tissue tixels using image data
image = Read10X_Image(image.dir = spatialFolder, filter.matrix = TRUE)
sequenced_tixels <- row.names(meta.data)
image <- image[sequenced_tixels, ]
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]

############### Dimension reduction, clustering, and add UMAP embedding
proj_in_tissue <- addIterativeLSI(
   ArchRProj = proj_in_tissue,
   useMatrix = "TileMatrix",
   name = "IterativeLSI",
   iterations = 2,
   clusterParams = list(
   resolution = c(0.5),
   sampleCells = 10000,
   n.start = 10
   ),
   varFeatures = 50000,
   dimsToUse = 1:30,
   force = TRUE
)

proj_in_tissue <- addClusters(
   input = proj_in_tissue,
   reducedDims = "IterativeLSI",
   method = "Seurat",
   name = "Clusters",
   resolution = 0.2,
   force = TRUE
)

proj_in_tissue <- addUMAP(
   ArchRProj = proj_in_tissue,
   reducedDims = "IterativeLSI",
   name = "UMAP",
   nNeighbors = 30,
   minDist = 0.5,
   metric = "cosine",
   force = TRUE
)

############## Creating Seurat object
gene_score <- getMatrixFromProject(proj_in_tissue)
rownames(gene_score) <- rowData(gene_score)$name
#proj_in_tissue <- addImputeWeights(proj_in_tissue)
#gene_score <- imputeMatrix(assay(gene_score), getImputeWeights(proj_in_tissue))
#gene_score <- log(gene_score+1, base = 2)
gene_score <- log(assay(gene_score)+1, base = 2)
colnames(gene_score) <- gsub(pattern = paste0(project_name, "#|-1"), replacement = "", x= colnames(gene_score))

object <- CreateSeuratObject(counts = gene_score, assay = "Spatial", meta.data = meta.data)

image <- image[Cells(x = object)]
DefaultAssay(object = image) <- "Spatial"
object[["slice1"]] <- image
spatial_in_tissue.obj <- object

spatial_in_tissue.obj$orig.ident = as.factor(project_name)
Idents(spatial_in_tissue.obj) = "orig.ident"
spatial_in_tissue.obj = AddMetaData(spatial_in_tissue.obj, spatial_in_tissue.obj@images$slice1@coordinates)


############## Define aesthetic parameters
n_clusters <- length(unique(proj_in_tissue$Clusters))
palette  = c("navyblue", "turquoise2", "tomato", "tan2", "pink", "mediumpurple1", "steelblue", "springgreen2","violetred", "orange", "violetred", "slateblue1",  "violet", "purple", "purple3","blue2",  "pink", "coral2", "palevioletred", "red2", "yellowgreen", "palegreen4", "wheat2", "tan", "tan3", "brown", "grey70", "grey50", "grey30")
cols <- palette[seq_len(n_clusters)]
names(cols) <- names(proj_in_tissue@sampleMetadata)
names(cols) <- paste0('C', seq_len(n_clusters))
cols_hex <- lapply(X = cols, FUN = function(x){
    do.call(rgb, as.list(col2rgb(x)/255))
})
cols <- unlist(cols_hex)
pt_size_factor <- 1

############## Plotting UMAP/cluster identities to spatial histology
spatial_in_tissue.obj@meta.data$Clusters = proj_in_tissue$Clusters
plot_spatial = Seurat::SpatialDimPlot(
    spatial_in_tissue.obj,
    label = FALSE, label.size = 3,
    group.by = "Clusters",
    pt.size.factor = pt_size_factor, cols = cols, stroke = 0) +
    theme(
       plot.title = element_blank(),
       legend.position = "right",
       text=element_text(size=21)) +
       ggtitle(project_name) + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=21))

plot_spatial$layers[[1]]$aes_params <- c(plot_spatial$layers[[1]]$aes_params, shape=22)

plot_umap = plotEmbedding(
  ArchRProj = proj_in_tissue,
  pal = cols,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP",
  size = 2) +
  theme(
    plot.title = element_blank(),
    legend.position = "none",
    text=element_text(size=21))

cluster_plots <- plot_spatial + plot_umap
cluster_plots

# we do basic clustering analysis based on LSI and UMAP. Also draw a plot in a spatial domain

## QC 
############## Plotting quality control metrics to spatial histology
spatial_in_tissue.obj@meta.data$log10_nFrags <- log10(spatial_in_tissue.obj@meta.data$nFrags)
plot_metadata = SpatialFeaturePlot(
  object = spatial_in_tissue.obj,
  features = c("log10_nFrags", "NucleosomeRatio", "TSSEnrichment"),
  alpha = c(0.2, 1), pt.size.factor = pt_size_factor) +
  theme(plot.title = element_text(hjust = 0.5), text=element_text(size=10))
plot_metadata$layers[[1]]$aes_params <-c(plot_metadata$layers[[1]]$aes_params, shape=22)

plot_metadata



## Finding marker genes 
markersGS <- getMarkerFeatures(ArchRProj = proj_in_tissue, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
plot_marker_heatmap = plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25", log2Norm=FALSE, nPrint=100)
plotMarkerHeatmap = ComplexHeatmap::draw(plot_marker_heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")



## Peak calling
pathToMacs2 <- findMacs2()

proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Clusters")
proj_in_tissue = addReproduciblePeakSet(ArchRProj = proj_in_tissue, groupBy = "Clusters", pathToMacs2 = pathToMacs2)
proj_in_tissue <- addPeakMatrix(proj_in_tissue)

getPeakSet(proj_in_tissue)

markersPeaks <- getMarkerFeatures(ArchRProj = proj_in_tissue, useMatrix = "PeakMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.05", transpose = TRUE)
plot_heatmap_peaks <- draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")



## Motif enrichment 
proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif")
proj_in_tissue <- addBgdPeaks(proj_in_tissue)
proj_in_tissue <- addDeviationsMatrix(ArchRProj = proj_in_tissue, peakAnnotation = "Motif")

markersMotifs <- getMarkerFeatures(ArchRProj = proj_in_tissue, useMatrix = "MotifMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon", useSeqnames = 'z')
markerMotifsList <- getMarkers(markersMotifs, cutOff = "FDR <= 0.05")
plot_marker_motif_heatmap = plotMarkerHeatmap(seMarker = markersMotifs, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")



## Footprinting
seFoot <- getFootprints(ArchRProj = proj_in_tissue, positions = motifPositions[markerMotifs], groupBy = "Clusters")
plotFootprints(seFoot = seFoot, ArchRProj = proj_in_tissue, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias-melanoma-d3g5", addDOC = FALSE, smoothWindow = 5)
plotFootprints(seFoot = seFoot, ArchRProj = proj_in_tissue, normMethod = "Divide", plotName = "Footprints-Divide-Bias-melanoma-d3g5", addDOC = FALSE, smoothWindow = 5)


