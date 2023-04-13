library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggplotify)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(grid)
library(shinycssloaders)
library(patchwork)

# load font for plot
if(!system.file(package="showtext")=="") sysfonts::font_add_google(name = "Lato", family = "Lato")
if(!system.file(package="showtext")=="") showtext::showtext_auto()

if(!exists("sc1conf")) sc1conf = readRDS("sc1conf.rds")
if(!exists("sc1def")) sc1def  = readRDS("sc1def.rds")
if(!exists("sc1gene")) sc1gene = readRDS("sc1gene.rds")
if(!exists("sc1meta")) sc1meta = readRDS("sc1meta.rds")

### Initialise variables and functions ----

# Colour palette
cList <- list(
  c(
    "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
  ),
  c(
    "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027"
  )[c(1, 1:9, 9)],
  c(
    "#FDE725", "#AADC32", "#5DC863", "#27AD81", "#21908C",
    "#2C728E", "#3B528B", "#472D7B", "#440154"
  )
)
names(cList) <- c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

# Panel sizes
pList <- c("400px", "600px", "800px")
names(pList) <- c("Small", "Medium", "Large")
pList2 <- c("500px", "700px", "900px")
names(pList2) <- c("Small", "Medium", "Large")
pList3 <- c("600px", "800px", "1000px")
names(pList3) <- c("Small", "Medium", "Large")
# baseplot font size
sList <- c(12, 14, 18, 22)
names(sList) <- c("Smaller", "Small", "Medium", "Large")
# ggrepel font size
lList <- c(5, 6, 7)
names(lList) <- c("Small", "Medium", "Large")

# Function to extract legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# progress indicator
show_progress <- function(...){
  return(shinycssloaders::withSpinner(..., type = 7, color = "#95a5a6"))
}

# Plot theme
# @description Custom ggplot theme
# @param base_size (Numeric) Base font size
# @param XYval (Logical) Show XY axes text?
# @param Xang (Numeric) X axis text angle
# @param XjusH (Numeric) X axis horizontal justification
# @param lpos (Character) Position of Legend
# @param font (Character) Google font
# @param col_text (Character) Text colour
# @param col_line (Character) Line colour
#
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5, lpos = "bottom", font = "Lato", col_text = "grey30", col_line = "grey60") {
  oupTheme <- theme(
    text = element_text(size = base_size, family = font),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(colour = col_line),
    axis.ticks = element_line(colour = col_line),
    axis.title = element_text(colour = col_text),
    axis.text = element_text(size = base_size, colour = col_text),
    axis.text.x = element_text(angle = Xang, hjust = XjusH),
    strip.background = element_rect(colour = "white"),
    legend.position = lpos,
    legend.key = element_rect(colour = NA, fill = NA)
  )
  if (!XYval) {
    oupTheme <- oupTheme +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
      )
  }
  
  return(oupTheme)
}

### Plotting functions ----

# @description DR scatterplot for gene expression
# @param dtab (Data.table) with columns X, Y, geneName and val
# @param bgCells (Logical) Background points
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scScatter <- function(dtab, bgCells = FALSE, inpdrX, inpdrY, inpsiz, inpcol, inpfsz, inpasp, inptxt){
  
  if(any(!c("X", "Y", "geneName", "val") %in% colnames(dtab))) "Input missing one or more columns: X, Y, geneName, val."
  
  ggOut <- ggplot(dtab, aes(X, Y, color = val))
  rat <- (max(dtab$X) - min(dtab$X)) / (max(dtab$Y) - min(dtab$Y))
  ltitle <- dtab$geneName[1]
  
  if (bgCells) {
    ggOut <- ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }
  
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) + xlab(inpdrX) + ylab(inpdrY) +
    scale_color_gradientn(ltitle, colours = cList[[inpcol]]) +
    #guides(color = guide_colorbar(barwidth = 20)) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt, lpos = "right")
  
  if (inpasp == "Square") {
    ggOut <- ggOut + coord_fixed(ratio = rat)
  } else if (inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }
  
  return(ggOut)
}

# @description Plot gene expression on dimred for any number of input genes
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to expression h5
# @param inpGene (Numeric) Named gene expression vector
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
# @param inpncol (Integer) Number of rows of plots
# @details 
# Config table contains columns ID (Character, columns name in metadata), UI (Character, UI id), fID (Character, Levels for categorical data, | separated), fCL (Character, Colours for categorical data, | separated), fRow (Integer, number of rows), grp (Logical), dimred (Logical)
#
scFeature <- function(inpConf, inpMeta, inpdrX, inpdrY, inp, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inpncol = 0){
  
  if (is.null(inpsub1)) inpsub1 <- inpConf$UI[1]
  
  # Identify genes that are in our dataset
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 36, "More than 36 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 0, "Please input at least 1 gene to plot!"))
  
  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData <- data.table()
  for (iGene in geneList$gene) {
    tmp <- inpMeta[, c("sampleID",inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, inpConf[UI == inpsub1]$ID),with = FALSE]
    colnames(tmp) <- c("sampleID", "X", "Y", "sub")
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr = )))
    ggData <- rbindlist(list(ggData, tmp))
  }
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) ggData <- ggData[sub %in% inpsub2]
  
  bgCells <- FALSE
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }
  
  if (inpord == "Max") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min") {
    ggData <- ggData[order(-val)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }
  
  ggDataSplit <- split(ggData, by=c("geneName"), flatten=FALSE)
  plist <- vector("list", length = length(ggDataSplit))
  for(i in seq_along(ggDataSplit)){
    plist[[i]] <- scScatter(dtab = ggDataSplit[[i]], bgCells, inpdrX, inpdrY, inpsiz, inpcol, inpfsz, inpasp, inptxt)
  }

  if(inpncol < 1) inpncol <- floor(sqrt(length(plist)))
  ggOut <- wrap_plots(plist) + plot_layout(ncol = inpncol)
  return(ggOut)
}

# @description Plot cell information on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
# @param inplab 
#
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "val", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
  bgCells <- FALSE
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }

  if (inpord == "Max") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min") {
    ggData <- ggData[order(-val)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }

  # Do factoring if required
  if (!is.na(inpConf[UI == inp1]$fCL)) {
    ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
    names(ggCol) <- levels(ggData$val)
    ggLvl <- levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)]
    ggData$val <- factor(ggData$val, levels = ggLvl)
    ggCol <- ggCol[ggLvl]
  }

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y, color = val))
  if (bgCells) {
  ggOut <- ggOut +
    geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) +
    xlab(inpdrX) +
    ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt)

  if (is.na(inpConf[UI == inp1]$fCL)) {
    ggOut <- ggOut +
    scale_color_gradientn("", colours = cList[[inpcol]]) +
    guides(color = guide_colorbar(barwidth = 20)) 
  } else {
    sListX <- min(nchar(paste0(levels(ggData$val), collapse = "")), 200)
    sListX <- 0.75 * (sList - (1.5 * floor(sListX / 50)))
    ggOut <- ggOut + scale_color_manual("", values = ggCol) +
      guides(color = guide_legend(
        override.aes = list(size = 5),
        nrow = inpConf[UI == inp1]$fRow
      )) +
      theme(legend.text = element_text(size = sListX[inpfsz]))

    if (inplab) { ggData3 <- ggData[, .(X = mean(X), Y = mean(Y)), by = "val"]
      lListX <- min(nchar(paste0(ggData3$val, collapse = "")), 200)
      lListX <- lList - (0.25 * floor(lListX / 50))
      ggOut <- ggOut +
        geom_text_repel(
          data = ggData3, aes(X, Y, label = val),
          color = "grey10", bg.color = "grey95", bg.r = 0.15,
          size = lListX[inpfsz], seed = 42
        ) }
    }

  if (inpasp == "Square") { ggOut <- ggOut + coord_fixed(ratio = rat) } else if (inpasp == "Fixed") { ggOut <- ggOut + coord_fixed() }

  if(is.numeric(inpMeta[[inp1]])) ggOut <- ggOut + guides(color = guide_colorbar(barwidth = 25))
  return(ggOut)
}

# @description Cell number, stats data
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) Name of metadata column for coloring points
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Name/path to h5 object
# @param inpGene (Character) Gene to use
# @param inpsplt (Character) Whether to split table by Quartile or Decile
#
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inpsplt) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
    with = FALSE
  ]

  colnames(ggData) <- c("group", "sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = T)))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { ggData <- ggData[sub %in% inpsub2] }

  # Split inp1 if necessary
  if (is.na(inpConf[UI == inp1]$fCL)) { if (inpsplt == "Quartile") { nBk <- 4 }
    if (inpsplt == "Decile") { nBk <- 10 }
    ggData$group <- cut(ggData$group, breaks = nBk) }

  # Actual data.table
  ggData$express <- FALSE
  ggData[val2 > 0]$express <- TRUE
  ggData1 <- ggData[express == TRUE, .(nExpress = .N), by = "group"]
  ggData <- ggData[, .(nCells = .N), by = "group"]
  ggData <- ggData1[ggData, on = "group"]
  ggData <- ggData[, c("group", "nCells", "nExpress"), with = FALSE]
  ggData[is.na(nExpress)]$nExpress <- 0
  ggData$pctExpress <- 100 * ggData$nExpress / ggData$nCells
  ggData <- ggData[order(group)]
  colnames(ggData)[3] <- paste0(colnames(ggData)[3], "_", inp2)

  return(ggData)
}

# @description Plot gene expression on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val < 0]$val <- 0
  h5file$close_all()
  bgCells <- FALSE

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2] }

  if (inpord == "Max") {
    ggData <- ggData[order(val)]
  } else if (inpord == "Min") { 
    ggData <- ggData[order(-val)] 
  } else if (inpord == "Random") { 
    ggData <- ggData[sample(nrow(ggData))] 
  }

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y, color = val))
  if (bgCells) { 
    ggOut <- ggOut +
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }

  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20) + xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) +
    guides(color = guide_colorbar(barwidth = 20))

  if (inpasp == "Square") { 
    ggOut <- ggOut + coord_fixed(ratio = rat) 
  } else if (inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }

  return(ggOut)
}

# Plot gene coexpression on dimred
bilinear <- function(x, y, xy, Q11, Q21, Q12, Q22) {
  oup <- (xy - x) * (xy - y) * Q11 + x * (xy - y) * Q21 + (xy - x) * y * Q12 + x * y * Q22
  oup <- oup / (xy * xy)
  return(oup)
}

# @description Gene Co-expression on dimred with legend
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRcoexFull <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                     inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                     inpasp, inptxt) {
  g1 <- scDRcoex(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                       inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                       inpasp, inptxt)
  g2 <- scDRcoexLeg(inp1, inp2, inpcol, inpfsz)
  g <- wrap_plots(g1, g2) + plot_layout(ncol = 2, nrow = 1, widths = c(10, 2))
  return(g)
}

# @description Gene Co-expression on dimred
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inpdrX (Character) X axis variable for DR
# @param inpdrY (Character) Y axis variable for DR
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpsiz (Numeric) Point size
# @param inpcol (Character) Custom colour label
# @param inpord (Character) Custom plotting order
# @param inpfsz (Character) Custom font size
# @param inppasp (Character) Custom aspect ratio
# @param inptxt (Logical) Show XY labels
#
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, inpsub1,
                     inpsub2, inpH5, inpGene, inpsiz, inpcol, inpord, inpfsz,
                     inpasp, inptxt) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "Y", "sub")
  rat <- (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))

  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val1 <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val1 < 0]$val1 <- 0
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  bgCells <- FALSE

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells <- TRUE
    ggData2 <- ggData[!sub %in% inpsub2]
    ggData <- ggData[sub %in% inpsub2]
  }

  # Generate coex color palette
  cInp <- strsplit(inpcol, "; ")[[1]]
  if (cInp[1] == "Red (Gene1)") {
    c10 <- c(255, 0, 0) 
  } else if (cInp[1] == "Orange (Gene1)") {
    c10 <- c(255, 140, 0)
  } else {
    c10 <- c(0, 255, 0)
  }

  if (cInp[2] == "Green (Gene2)") { c01 <- c(0, 255, 0) } else { c01 <- c(0, 0, 255) }

  c00 <- c(217, 217, 217)
  c11 <- c10 + c01
  nGrid <- 16
  nPad <- 2
  nTot <- nGrid + nPad * 2
  gg <- data.table(v1 = rep(0:nTot, nTot + 1), v2 = sort(rep(0:nTot, nTot + 1)))
  gg$vv1 <- gg$v1 - nPad
  gg[vv1 < 0]$vv1 <- 0
  gg[vv1 > nGrid]$vv1 <- nGrid
  gg$vv2 <- gg$v2 - nPad
  gg[vv2 < 0]$vv2 <- 0
  gg[vv2 > nGrid]$vv2 <- nGrid
  gg$cR <- bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG <- bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB <- bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix <- rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg <- gg[, c("v1", "v2", "cMix")]

  # Map colours
  ggData$v1 <- round(nTot * ggData$val1 / max(ggData$val1))
  ggData$v2 <- round(nTot * ggData$val2 / max(ggData$val2))
  ggData$v0 <- ggData$v1 + ggData$v2
  ggData <- gg[ggData, on = c("v1", "v2")]

  if (inpord == "Max") {
    ggData <- ggData[order(v0)]
  } else if (inpord == "Min") {
    ggData <- ggData[order(-v0)]
  } else if (inpord == "Random") {
    ggData <- ggData[sample(nrow(ggData))]
  }

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, Y))
  if (bgCells) {
    ggOut <- ggOut +
    geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 20)
  }
  
  ggOut <- ggOut +
    geom_point(size = inpsiz, shape = 20, color = ggData$cMix) +
    xlab(inpdrX) + ylab(inpdrY) +
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    scale_color_gradientn(inp1, colours = cList[[1]]) +
    guides(color = guide_colorbar(barwidth = 20))

  if (inpasp == "Square") {
    ggOut <- ggOut + coord_fixed(ratio = rat)
  } else if (inpasp == "Fixed") {
    ggOut <- ggOut + coord_fixed()
  }
  
  return(ggOut)
}
  
# Co-exp plot legend
# @param inp1
# @param inp2
# @param inpcol Colour
# @param inpfsz Font size
#
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz) {
  
  # Generate coex color palette
  cInp <- strsplit(inpcol, "; ")[[1]]
  if (cInp[1] == "Red (Gene1)") {
    c10 <- c(255, 0, 0)
  } else if (cInp[1] == "Orange (Gene1)") { 
    c10 <- c(255, 140, 0) 
  } else { 
    c10 <- c(0, 255, 0) 
  }

  if (cInp[2] == "Green (Gene2)") { 
    c01 <- c(0, 255, 0) 
  } else { 
    c01 <- c(0, 0, 255) 
  }

  c00 <- c(217, 217, 217)
  c11 <- c10 + c01
  nGrid <- 16
  nPad <- 2
  nTot <- nGrid + nPad * 2
  gg <- data.table(v1 = rep(0:nTot, nTot + 1), v2 = sort(rep(0:nTot, nTot + 1)))
  gg$vv1 <- gg$v1 - nPad
  gg[vv1 < 0]$vv1 <- 0
  gg[vv1 > nGrid]$vv1 <- nGrid
  gg$vv2 <- gg$v2 - nPad
  gg[vv2 < 0]$vv2 <- 0
  gg[vv2 > nGrid]$vv2 <- nGrid
  gg$cR <- bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1])
  gg$cG <- bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2])
  gg$cB <- bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3])
  gg$cMix <- rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255)
  gg <- gg[, c("v1", "v2", "cMix")]

  # Actual ggplot
  ggOut <- ggplot(gg, aes(v1, v2)) +
    geom_tile(fill = gg$cMix) +
    xlab(inp1) +
    ylab(inp2) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) +
    sctheme(base_size = sList[inpfsz], XYval = TRUE)
  return(ggOut)
}


scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene) {
  
  if (is.null(inpsub1)) {
  inpsub1 <- inpConf$UI[1]
  }

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE]
  colnames(ggData) <- c("sub")
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData$val1 <- h5data$read(args = list(inpGene[inp1], quote(expr = )))
  ggData[val1 < 0]$val1 <- 0
  ggData$val2 <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
  ggData[val2 < 0]$val2 <- 0
  h5file$close_all()
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { ggData <- ggData[sub %in% inpsub2] }

  # Actual data.table
  ggData$express <- "none"
  ggData[val1 > 0]$express <- inp1
  ggData[val2 > 0]$express <- inp2
  ggData[val1 > 0 & val2 > 0]$express <- "both"
  ggData$express <- factor(ggData$express, levels = unique(c("both", inp1, inp2, "none")))
  ggData <- ggData[, .(nCells = .N), by = "express"]
  ggData$percent <- 100 * ggData$nCells / sum(ggData$nCells)
  ggData <- ggData[order(express)]
  colnames(ggData)[1] <- "expression > 0"
  
  return(ggData)
}

# Plot violin / boxplot / lineplot
# @description Violin plot gene expression
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) X axis cell info
# @param inp2 (Character) Y axis cell info / gene
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inptyp (Character) Plot type. "violin", "boxplot" or "lineplot".
# @param inppts (Logical) Should points be displayed?
# @param inpsiz (Numeric) Point size
# @param inpfsz (Character) Custom font size
# @param inpbarsz (Numeric) Bar size for lineplot
#
scVioBox <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inpH5, inpGene, inptyp, inppts, inpsiz, inpfsz, inpbarsz) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }
  if (is.null(inpbarsz)) inpbarsz <- 0.3

  # Prepare ggData
  ggData <- inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),
    with = FALSE
  ]
  colnames(ggData) <- c("X", "sub")

  # Load in either cell meta or gene expr
  if (inp2 %in% inpConf$UI) {
    ggData$val <- inpMeta[[inpConf[UI == inp2]$ID]]
  } else {
    h5file <- H5File$new(inpH5, mode = "r")
    h5data <- h5file[["grp"]][["data"]]
    ggData$val <- h5data$read(args = list(inpGene[inp2], quote(expr = )))
    ggData[val < 0]$val <- 0
    set.seed(42)
    tmpNoise <- rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000
    ggData$val <- ggData$val + tmpNoise
    h5file$close_all()
  }

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2]
  }

  # Do factoring
  ggCol <- strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]]
  names(ggCol) <- levels(ggData$X)
  ggLvl <- levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)]
  ggData$X <- factor(ggData$X, levels = ggLvl)
  ggCol <- ggCol[ggLvl]

  # Actual ggplot
  ggOut <- ggplot(ggData, aes(X, val, fill = X))
  
  if (inptyp == "violin") {
    ggOut <- ggOut + geom_violin(scale = "width")
  } else if (inptyp == "boxplot")  {
    ggOut <- ggOut + geom_boxplot()
  } else if (inptyp == "lineplot") {
    ggOut <- ggOut + geom_col(aes(col = X), position = position_dodge2(preserve = "single"), size = inpbarsz)
  }
  
  if (inppts) ggOut <- ggOut + geom_jitter(size = inpsiz, shape = 20, alpha = 0.4)

  ggOut <- ggOut + xlab(inp1) + ylab(inp2) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_colour_manual("", values = ggCol) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
    
  return(ggOut)
  }

# Plot proportion plot
# @description Proportion barplot
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp1 (Character) Gene name to use
# @param inp2 (Character) Gene name to use
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inptyp (Character) Plot type. "violin" else boxplot.
# @param inpflp (Logical) Flip coordinates?
# @param inpfsz (Character) Custom font size
#
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, inptyp, inpflp, inpfsz) {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Prepare ggData
  ggData <- inpMeta[, c(
    inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID,
    inpConf[UI == inpsub1]$ID
  ),
  with = FALSE
  ]
  colnames(ggData) <- c("X", "grp", "sub")

  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) { ggData <- ggData[sub %in% inpsub2] }
  ggData <- ggData[, .(nCells = .N), by = c("X", "grp")]
  ggData <- ggData[, { tot <- sum(nCells)
    .SD[, .(
      pctCells = 100 * sum(nCells) / tot,
      nCells = nCells
    ), by = "grp"] }, by = "X"]

  # Do factoring
  ggCol <- strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]]
  names(ggCol) <- levels(ggData$grp)
  ggLvl <- levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)]
  ggData$grp <- factor(ggData$grp, levels = ggLvl)
  ggCol <- ggCol[ggLvl]

  # Actual ggplot
  if (inptyp == "Proportion") {
    ggOut <- ggplot(ggData, aes(X, pctCells, fill = grp)) +
      geom_col() +
      ylab("Cell Proportion (%)") 
  } else { 
    ggOut <- ggplot(ggData, aes(X, nCells, fill = grp)) +
      geom_col() +
      ylab("Number of Cells")
  }

  if (inpflp) ggOut <- ggOut + coord_flip()

  ggOut <- ggOut + xlab(inp1) +
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "right")

  return(ggOut)
}

# Get gene list
scGeneList <- function(inp, inpGene) {
  geneList <- data.table(
  gene = inp,
  present = TRUE
  )
  geneList[!gene %in% names(inpGene)]$present <- FALSE
  return(geneList)
}

# Plot gene expression bubbleplot / heatmap
# @description dotplot / heatmap gene expression
# @param inpConf (data.frame) Configuration table
# @param inpMeta (data.frame) Metadata table
# @param inp (Character) Gene names
# @param inpGrp
# @param inpPlt
# @param inpsub1 (Character) Name of metadata column for subsetting
# @param inpsub2 (Character/Vector) Levels under metadata column for subsetting
# @param inpH5 (Character) Path to gene expression h5 file (sc1gexpr.h5)
# @param inpGene (integer) Named integer vector of gene expression values (sc1gene.rds)
# @param inpScl
# @param inpRow
# @param inpCol
# @param inpcols
# @param inpfsz (Character) Custom font size
# @param col_line (Character) Line colour
#
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, inpcols, inpfsz, col_line = "grey60") {
  
  if (is.null(inpsub1)) { inpsub1 <- inpConf$UI[1] }

  # Identify genes that are in our dataset
  geneList <- scGeneList(inp, inpGene)
  geneList <- geneList[present == TRUE]
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))

  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData <- data.table()
  for (iGene in geneList$gene) {
    tmp <- inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) <- c("sampleID", "sub")
    tmp$grpBy <- inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName <- iGene
    tmp$val <- h5data$read(args = list(inpGene[iGene], quote(expr = )))
    ggData <- rbindlist(list(ggData, tmp)) 
  }
  h5file$close_all()
  
  if (length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    ggData <- ggData[sub %in% inpsub2] 
  }
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))

  # Aggregate
  ggData$val <- expm1(ggData$val)
  ggData <- ggData[, .(val = mean(val), prop = sum(val > 0) / length(sampleID)),
    by = c("geneName", "grpBy")
  ]
  ggData$val <- log1p(ggData$val)

  # Scale if required
  colRange <- range(ggData$val)
  if (inpScl) {
    ggData[, val := scale(val), keyby = "geneName"]
    colRange <- c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  }

  # hclust row/col if necessary
  ggMat <- dcast.data.table(ggData, geneName ~ grpBy, value.var = "val")
  tmp <- ggMat$geneName
  ggMat <- as.matrix(ggMat[, -1])
  rownames(ggMat) <- tmp
  if (inpRow) {
    hcRow <- dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow <- ggplot() +
      coord_flip() +
      geom_segment(data = hcRow$segments, aes(x = x, y = y, xend = xend, yend = yend), col = col_line) +
      scale_y_continuous(
        breaks = rep(0, uniqueN(ggData$grpBy)),
        labels = unique(ggData$grpBy), expand = c(0, 0)
      ) +
      scale_x_continuous(
        breaks = seq_along(hcRow$labels$label),
        labels = hcRow$labels$label, expand = c(0, 0.5)
      ) +
      sctheme(base_size = sList[inpfsz]) +
      theme(
        axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "white", angle = 45, hjust = 1)
      )
    ggData$geneName <- factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName <- factor(ggData$geneName, levels = rev(geneList$gene)) 
  }

  if (inpCol) {
    hcCol <- dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol <- ggplot() +
      geom_segment(data = hcCol$segments, aes(x = x, y = y, xend = xend, yend = yend), col = col_line) +
      scale_x_continuous(
        breaks = seq_along(hcCol$labels$label),
        labels = hcCol$labels$label, expand = c(0.05, 0)
      ) +
      scale_y_continuous(
        breaks = rep(0, uniqueN(ggData$geneName)),
        labels = unique(ggData$geneName), expand = c(0, 0)
      ) +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      theme(
        axis.title = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(color = "white")
      )
    ggData$grpBy <- factor(ggData$grpBy, levels = hcCol$labels$label)
    }

  # Actual plot according to plottype
  if (inpPlt == "Bubbleplot") {
    # Bubbleplot
    ggOut <- ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) +
      geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_size_continuous("proportion",
        range = c(0, 8),
        limits = c(0, 1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)
      ) +
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 20)) +
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else {
    # Heatmap
    ggOut <- ggplot(ggData, aes(grpBy, geneName, fill = val)) +
      geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(fill = guide_colorbar(barwidth = 20)) +
      theme(axis.title = element_blank()) 
  }

  # Final tidy
  ggLeg <- g_legend(ggOut)
  ggOut <- ggOut + theme(legend.position = "none")
  
  if (inpRow & inpCol) {
    ggOut <- wrap_plots(ggCol, plot_spacer(), ggOut, ggRow, ggplotify::as.ggplot(ggLeg))+
              plot_layout(ncol = 2, widths = c(7,1), heights = c(1,7,2))
  } else if (inpRow) {
    ggOut <- wrap_plots(ggOut, ggRow, ggplotify::as.ggplot(ggLeg))+
              plot_layout(ncol = 2, widths = c(7,1), heights = c(7,2))
  } else if (inpCol) {
    ggOut <- wrap_plots(ggCol, ggOut, ggplotify::as.ggplot(ggLeg))+
               plot_layout(ncol = 1, heights = c(1,7,2))
  } else {
    ggOut <- wrap_plots(ggOut, ggplotify::as.ggplot(ggLeg))+
      plot_layout(ncol = 1, heights = c(7,2))
  }
  
  return(ggOut)
}


### Server code ----
shinyServer(function(input, output, session) {

### For all tags and Server-side selectize
observe_helpers()
optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
updateSelectizeInput(session, "sc1_civge_inp2", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_gevge_inp1", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_gevge_inp2", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_gec_inp1", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$gene1, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_gec_inp2", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$gene2, options = list(
                       maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_gem_inp", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$genes[1:9], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_hea_inp", choices = names(sc1gene), server = TRUE,
                     selected = sc1def$genes[1:12], options = list(
                     create = TRUE, persist = TRUE, render = I(optCrt)))
updateSelectizeInput(session, "sc1_vio_inp2", server = TRUE,
                     choices = c(sc1conf[is.na(fID)]$UI,names(sc1gene)),
                     selected = sc1conf[is.na(fID)]$UI[1], options = list(
                       maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
                       create = TRUE, persist = TRUE, render = I(optCrt)))  
### Tab civge cell info vs gene exp ----

  output$sc1_civge_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_civge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_civge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_civge_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_civge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_civge_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_civge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_civge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_civge_oup1 <- renderPlot({
  req(input$sc1_civge_inp1)
  scDRcell(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp1, input$sc1_civge_sub1, input$sc1_civge_sub2, input$sc1_civge_siz, input$sc1_civge_col1, input$sc1_civge_ord1, input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt, input$sc1_civge_lab1)
})

output$sc1_civge_oup1.ui <- renderUI({
  show_progress(imageOutput("sc1_civge_oup1", height = pList[input$sc1_civge_psz]))
})

output$sc1_civge_oup1.png <- downloadHandler(
 filename = function() { tolower(paste0("sc1", "_", input$sc1_civge_drX, "_", input$sc1_civge_drY, "_", input$sc1_civge_inp1, ".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$sc1_civge_oup1.res, bg = "white",
   plot = scDRcell(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp1,   input$sc1_civge_sub1, input$sc1_civge_sub2, input$sc1_civge_siz, input$sc1_civge_col1, input$sc1_civge_ord1,  input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt, input$sc1_civge_lab1)
   )
})

output$sc1_civge_oup1.pdf <- downloadHandler(
 filename = function() { tolower(paste0("sc1", "_", input$sc1_civge_drX,"_", input$sc1_civge_drY,"_", input$sc1_civge_inp1,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRcell(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp1,   input$sc1_civge_sub1, input$sc1_civge_sub2, input$sc1_civge_siz, input$sc1_civge_col1, input$sc1_civge_ord1,  input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt, input$sc1_civge_lab1)
   )
})

output$sc1_civge_oup1.svg <- downloadHandler(
 filename = function() { tolower(paste0("sc1", "_", input$sc1_civge_drX,"_", input$sc1_civge_drY,"_", input$sc1_civge_inp1,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRcell(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp1,   input$sc1_civge_sub1, input$sc1_civge_sub2, input$sc1_civge_siz, input$sc1_civge_col1, input$sc1_civge_ord1,  input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt, input$sc1_civge_lab1)
   )
})

output$sc1_civge_.dt <- renderDataTable({
 req(input$sc1_civge_inp2)
 ggData = scDRnum(sc1conf, sc1meta, input$sc1_civge_inp1, input$sc1_civge_inp2, input$sc1_civge_sub1, input$sc1_civge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_civge_splt)
 datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
   formatRound(columns = c("pctExpress"), digits = 2)
})

output$sc1_civge_oup2 <- renderPlot({
 req(input$sc1_civge_inp2)
 scDRgene(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp2, input$sc1_civge_sub1, input$sc1_civge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_civge_siz, input$sc1_civge_col2, input$sc1_civge_ord2, input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt)
})

output$sc1_civge_oup2.ui <- renderUI({
 show_progress(imageOutput("sc1_civge_oup2", height = pList[input$sc1_civge_psz]))
})

output$sc1_civge_oup2.png <- downloadHandler(
 filename = function() { tolower(paste0("sc1", "_", input$sc1_civge_drX,"_",input$sc1_civge_drY,"_", input$sc1_civge_inp2,".png")) },
 content = function(file) {
   ggsave(
   file, device = "png", dpi = input$sc1_civge_oup2.res, bg = "white",
   plot = scDRgene(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp2, input$sc1_civge_sub1, input$sc1_civge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_civge_siz, input$sc1_civge_col2, input$sc1_civge_ord2, input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt)
   )
})

output$sc1_civge_oup2.pdf <- downloadHandler(
 filename = function() { tolower(paste0("sc1", "_", input$sc1_civge_drX,"_",input$sc1_civge_drY,"_", input$sc1_civge_inp2,".pdf")) },
 content = function(file) {
   ggsave(
   file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
   plot = scDRgene(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp2,  input$sc1_civge_sub1, input$sc1_civge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_civge_siz, input$sc1_civge_col2, input$sc1_civge_ord2, input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt)
   )
}) 

output$sc1_civge_oup2.svg <- downloadHandler(
 filename = function() { tolower(paste0("sc1", "_", input$sc1_civge_drX,"_",input$sc1_civge_drY,"_", input$sc1_civge_inp2,".svg")) },
 content = function(file) {
   ggsave(
   file, device = "svg", bg = "white",
   plot = scDRgene(sc1conf, sc1meta, input$sc1_civge_drX, input$sc1_civge_drY, input$sc1_civge_inp2,  input$sc1_civge_sub1, input$sc1_civge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_civge_siz, input$sc1_civge_col2, input$sc1_civge_ord2, input$sc1_civge_fsz, input$sc1_civge_asp, input$sc1_civge_txt)
   )
}) # End of tab civge



### Tab civci cell info vs cell info ----
  
  output$sc1_civci_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_civci_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_civci_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_civci_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_civci_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_civci_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_civci_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_civci_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_civci_oup1 <- renderPlot({
  req(input$sc1_civci_inp1)
  scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp1, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col1, input$sc1_civci_ord1, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab1)
})

output$sc1_civci_oup1.ui <- renderUI({
  show_progress(imageOutput("sc1_civci_oup1", height = pList[input$sc1_civci_psz]))
})

output$sc1_civci_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_civci_drX, "_", input$sc1_civci_drY, "_", input$sc1_civci_inp1, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", dpi = input$sc1_civci_oup1.res, bg = "white",
    plot = scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp1, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col1, input$sc1_civci_ord1, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab1)
    )
})

output$sc1_civci_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_civci_drX, "_", input$sc1_civci_drY, "_", input$sc1_civci_inp1, ".pdf")) },
  content = function(file) { ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp1, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col1, input$sc1_civci_ord1, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab1) )
})

output$sc1_civci_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_civci_drX, "_", input$sc1_civci_drY, "_", input$sc1_civci_inp1, ".svg")) },
  content = function(file) { ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp1, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col1, input$sc1_civci_ord1, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab1) )
})

output$sc1_civci_oup2 <- renderPlot({
  req(input$sc1_civci_inp2)
  scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp2, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col2, input$sc1_civci_ord2, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab2)
})

output$sc1_civci_oup2.ui <- renderUI({
  show_progress(imageOutput("sc1_civci_oup2", height = pList[input$sc1_civci_psz]))
})

output$sc1_civci_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_civci_drX,"_",input$sc1_civci_drY,"_", input$sc1_civci_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$sc1_civci_oup2.res,
    plot = scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp2, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col2, input$sc1_civci_ord2, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab2)
    )
}) 

output$sc1_civci_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_civci_drX,"_",input$sc1_civci_drY,"_", input$sc1_civci_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp2, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col2, input$sc1_civci_ord2, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab2) 
    )
})

output$sc1_civci_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_civci_drX,"_",input$sc1_civci_drY,"_", input$sc1_civci_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcell(sc1conf, sc1meta, input$sc1_civci_drX, input$sc1_civci_drY, input$sc1_civci_inp2, input$sc1_civci_sub1, input$sc1_civci_sub2, input$sc1_civci_siz, input$sc1_civci_col2, input$sc1_civci_ord2, input$sc1_civci_fsz, input$sc1_civci_asp, input$sc1_civci_txt, input$sc1_civci_lab2) 
    )
}) # End of tab civci



### Tab gevge gene exp vs gene exp ----

  output$sc1_gevge_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_gevge_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_gevge_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_gevge_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_gevge_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_gevge_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_gevge_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_gevge_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_gevge_oup1 <- renderPlot({
  req(input$sc1_gevge_inp1)
  scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp1, input$sc1_gevge_sub1, input$sc1_gevge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gevge_siz, input$sc1_gevge_col1, input$sc1_gevge_ord1, input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
})

output$sc1_gevge_oup1.ui <- renderUI({
  show_progress(imageOutput("sc1_gevge_oup1", height = pList[input$sc1_gevge_psz]))
})

output$sc1_gevge_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gevge_drX, "_", input$sc1_gevge_drY, "_", input$sc1_gevge_inp1,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$sc1_gevge_oup1.res,
    plot = scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp1, 
                    input$sc1_gevge_sub1, input$sc1_gevge_sub2,
                    "sc1gexpr.h5", sc1gene,
                    input$sc1_gevge_siz, input$sc1_gevge_col1, input$sc1_gevge_ord1,
                    input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
    )
})

output$sc1_gevge_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gevge_drX, "_", input$sc1_gevge_drY, "_", input$sc1_gevge_inp1, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp1, input$sc1_gevge_sub1, input$sc1_gevge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gevge_siz, input$sc1_gevge_col1, input$sc1_gevge_ord1, input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
    )
})

output$sc1_gevge_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gevge_drX, "_", input$sc1_gevge_drY, "_", input$sc1_gevge_inp1, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp1, input$sc1_gevge_sub1, input$sc1_gevge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gevge_siz, input$sc1_gevge_col1, input$sc1_gevge_ord1, input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
    )
})

output$sc1_gevge_oup2 <- renderPlot({
  req(input$sc1_gevge_inp2)
  scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp2, input$sc1_gevge_sub1, input$sc1_gevge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gevge_siz, input$sc1_gevge_col2, input$sc1_gevge_ord2, input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
})

output$sc1_gevge_oup2.ui <- renderUI({
  show_progress(imageOutput("sc1_gevge_oup2", height = pList[input$sc1_gevge_psz]))
})

output$sc1_gevge_oup2.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gevge_drX, "_", input$sc1_gevge_drY, "_", input$sc1_gevge_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$sc1_gevge_oup2.res,
    plot = scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp2, input$sc1_gevge_sub1, input$sc1_gevge_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gevge_siz, input$sc1_gevge_col2, input$sc1_gevge_ord2, input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
    )
}) 

output$sc1_gevge_oup2.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gevge_drX, "_", input$sc1_gevge_drY, "_", input$sc1_gevge_inp2,".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp2, 
                    input$sc1_gevge_sub1, input$sc1_gevge_sub2,
                    "sc1gexpr.h5", sc1gene,
                    input$sc1_gevge_siz, input$sc1_gevge_col2, input$sc1_gevge_ord2,
                    input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
    )
})

output$sc1_gevge_oup2.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gevge_drX, "_", input$sc1_gevge_drY, "_", input$sc1_gevge_inp2,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRgene(sc1conf, sc1meta, input$sc1_gevge_drX, input$sc1_gevge_drY, input$sc1_gevge_inp2, 
                    input$sc1_gevge_sub1, input$sc1_gevge_sub2,
                    "sc1gexpr.h5", sc1gene,
                    input$sc1_gevge_siz, input$sc1_gevge_col2, input$sc1_gevge_ord2,
                    input$sc1_gevge_fsz, input$sc1_gevge_asp, input$sc1_gevge_txt)
    )
}) # End of tab gevge




### Tab gem gene expression multi ----

output$sc1_gem_sub1.ui <- renderUI({
  sub = strsplit(sc1conf[UI == input$sc1_gem_sub1]$fID, "\\|")[[1]]
  checkboxGroupInput("sc1_gem_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
})
observeEvent(input$sc1_gem_sub1non, {
  sub = strsplit(sc1conf[UI == input$sc1_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "sc1_gem_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
})
observeEvent(input$sc1_gem_sub1all, {
  sub = strsplit(sc1conf[UI == input$sc1_gem_sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "sc1_gem_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
})

output$sc1_gem_oup1 <- renderPlot({
  req(input$sc1_gem_inp)
  
  scFeature(sc1conf, sc1meta, input$sc1_gem_drX, input$sc1_gem_drY, input$sc1_gem_inp, input$sc1_gem_sub1, input$sc1_gem_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gem_siz, input$sc1_gem_col, input$sc1_gem_ord, input$sc1_gem_fsz, input$sc1_gem_asp, input$sc1_gem_txt, input$sc1_gem_ncol)
})

output$sc1_gem_oup1.ui <- renderUI({
  show_progress(imageOutput("sc1_gem_oup1", height = pList[input$sc1_gem_psz]))
})

output$sc1_gem_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gem_drX, "_", input$sc1_gem_drY, "_expression.png")) },
  content = function(file) {
    ggsave(
      file, device = "png", height = input$sc1_gem_oup1.height, width = input$sc1_gem_oup1.width, dpi = input$sc1_gem_oup1.res, units = "cm", bg = "white",
      plot = scFeature(sc1conf, sc1meta, input$sc1_gem_drX, input$sc1_gem_drY, input$sc1_gem_inp, input$sc1_gem_sub1, input$sc1_gem_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gem_siz, input$sc1_gem_col, input$sc1_gem_ord, input$sc1_gem_fsz, input$sc1_gem_asp, input$sc1_gem_txt, input$sc1_gem_ncol)
    )
}) 

output$sc1_gem_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gem_drX, "_", input$sc1_gem_drY, "_expression.pdf")) },
  content = function(file) {
  pdf(file, useDingbats = FALSE, height = input$sc1_gem_oup1.height/2.54, width = input$sc1_gem_oup1.width/2.54, bg = "white", onefile = TRUE)
  if(!system.file(package="showtext")=="") showtext::showtext_begin()
  print(scFeature(sc1conf, sc1meta, input$sc1_gem_drX, input$sc1_gem_drY, input$sc1_gem_inp, input$sc1_gem_sub1, input$sc1_gem_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gem_siz, input$sc1_gem_col, input$sc1_gem_ord, input$sc1_gem_fsz, input$sc1_gem_asp, input$sc1_gem_txt, input$sc1_gem_ncol))
  if(!system.file(package="showtext")=="") showtext::showtext_end()
  dev.off()
})

output$sc1_gem_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gem_drX, "_", input$sc1_gem_drY, "_expression.svg")) },
  content = function(file) { ggsave(
    file, device = "svg", height = input$sc1_gem_oup1.height, width = input$sc1_gem_oup1.width, units = "cm", bg = "white",
    plot = scFeature(sc1conf, sc1meta, input$sc1_gem_drX, input$sc1_gem_drY, input$sc1_gem_inp, input$sc1_gem_sub1, input$sc1_gem_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gem_siz, input$sc1_gem_col, input$sc1_gem_ord, input$sc1_gem_fsz, input$sc1_gem_asp, input$sc1_gem_txt, input$sc1_gem_ncol))
}) # End of tab gem



### Tab gec gene co-expression ----

  output$sc1_gec_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_gec_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_gec_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_gec_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_gec_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_gec_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_gec_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_gec_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_gec_oup1 <- renderPlot({
  scDRcoexFull(sc1conf, sc1meta, input$sc1_gec_drX, input$sc1_gec_drY, input$sc1_gec_inp1, input$sc1_gec_inp2, input$sc1_gec_sub1, input$sc1_gec_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gec_siz, input$sc1_gec_col1, input$sc1_gec_ord1, input$sc1_gec_fsz, input$sc1_gec_asp, input$sc1_gec_txt)
})

output$sc1_gec_oup1.ui <- renderUI({
  show_progress(imageOutput("sc1_gec_oup1", height = pList2[input$sc1_gec_psz]))
})

output$sc1_gec_oup1.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gec_drX, "_", input$sc1_gec_drY, "_", input$sc1_gec_inp1, "_", input$sc1_gec_inp2, ".png")) },
  content = function(file) { ggsave(
    file, device = "png", bg = "white", dpi = input$sc1_gec_oup1.res,
    plot = scDRcoexFull(sc1conf, sc1meta, input$sc1_gec_drX, input$sc1_gec_drY, input$sc1_gec_inp1, input$sc1_gec_inp2, input$sc1_gec_sub1, input$sc1_gec_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gec_siz, input$sc1_gec_col1, input$sc1_gec_ord1, input$sc1_gec_fsz, input$sc1_gec_asp, input$sc1_gec_txt) )
})

output$sc1_gec_oup1.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gec_drX, "_", input$sc1_gec_drY, "_", input$sc1_gec_inp1, "_", input$sc1_gec_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scDRcoexFull(sc1conf, sc1meta, input$sc1_gec_drX, input$sc1_gec_drY, input$sc1_gec_inp1, input$sc1_gec_inp2, input$sc1_gec_sub1, input$sc1_gec_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gec_siz, input$sc1_gec_col1, input$sc1_gec_ord1, input$sc1_gec_fsz, input$sc1_gec_asp, input$sc1_gec_txt)
    )
})

output$sc1_gec_oup1.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_gec_drX, "_", input$sc1_gec_drY, "_", input$sc1_gec_inp1, "_", input$sc1_gec_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scDRcoexFull(sc1conf, sc1meta, input$sc1_gec_drX, input$sc1_gec_drY, input$sc1_gec_inp1, input$sc1_gec_inp2, input$sc1_gec_sub1, input$sc1_gec_sub2, "sc1gexpr.h5", sc1gene, input$sc1_gec_siz, input$sc1_gec_col1, input$sc1_gec_ord1, input$sc1_gec_fsz, input$sc1_gec_asp, input$sc1_gec_txt)
    )
})

output$sc1_gec_.dt <- renderDataTable({
  ggData = scDRcoexNum(sc1conf, sc1meta, input$sc1_gec_inp1, input$sc1_gec_inp2, input$sc1_gec_sub1, input$sc1_gec_sub2, "sc1gexpr.h5", sc1gene)
  datatable(ggData, rownames = FALSE, extensions = "Buttons", options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
            formatRound(columns = c("percent"), digits = 2)
}) # End of tab gec



### Tab vio violinplot / boxplot / lineplot ----

  output$sc1_vio_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_vio_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_vio_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_vio_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_vio_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_vio_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_vio_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_vio_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_vio_oup <- renderPlot({
  scVioBox(sc1conf, sc1meta, input$sc1_vio_inp1, input$sc1_vio_inp2, input$sc1_vio_sub1, input$sc1_vio_sub2, "sc1gexpr.h5", sc1gene, input$sc1_vio_typ, input$sc1_vio_pts, input$sc1_vio_siz, input$sc1_vio_fsz, input$sc1_vio_barsz)
})

output$sc1_vio_oup.ui <- renderUI({
  show_progress(imageOutput("sc1_vio_oup", height = pList2[input$sc1_vio_psz]))
})

output$sc1_vio_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_vio_typ, "_", input$sc1_vio_inp1, "_", input$sc1_vio_inp2,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$sc1_vio_oup.res,
    plot = scVioBox(sc1conf, sc1meta, input$sc1_vio_inp1, input$sc1_vio_inp2, input$sc1_vio_sub1, input$sc1_vio_sub2, "sc1gexpr.h5", sc1gene, input$sc1_vio_typ, input$sc1_vio_pts, input$sc1_vio_siz, input$sc1_vio_fsz, input$sc1_vio_barsz)
    )
})

output$sc1_vio_oup.pdf <- downloadHandler(
  filename = function() { paste0("sc1", input$sc1_vio_typ, "_", input$sc1_vio_inp1, "_", input$sc1_vio_inp2, ".pdf") },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white",
    plot = scVioBox(sc1conf, sc1meta, input$sc1_vio_inp1, input$sc1_vio_inp2, input$sc1_vio_sub1, input$sc1_vio_sub2, "sc1gexpr.h5", sc1gene, input$sc1_vio_typ, input$sc1_vio_pts, input$sc1_vio_siz, input$sc1_vio_fsz, input$sc1_vio_barsz)
    )
})

output$sc1_vio_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_vio_typ, "_", input$sc1_vio_inp1, "_", input$sc1_vio_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scVioBox(sc1conf, sc1meta, input$sc1_vio_inp1, input$sc1_vio_inp2, input$sc1_vio_sub1, input$sc1_vio_sub2, "sc1gexpr.h5", sc1gene, input$sc1_vio_typ, input$sc1_vio_pts, input$sc1_vio_siz, input$sc1_vio_fsz, input$sc1_vio_barsz)
    )
}) # End of tab vio




### Tab pro proportion plot ----

  output$sc1_pro_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_pro_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_pro_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_pro_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_pro_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_pro_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_pro_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_pro_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_pro_oup <- renderPlot({
  scProp(sc1conf, sc1meta, input$sc1_pro_inp1, input$sc1_pro_inp2, input$sc1_pro_sub1, input$sc1_pro_sub2, input$sc1_pro_typ, input$sc1_pro_flp, input$sc1_pro_fsz)
})

output$sc1_pro_oup.ui <- renderUI({
  show_progress(imageOutput("sc1_pro_oup", height = pList2[input$sc1_pro_psz]))
})

output$sc1_pro_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_pro_typ, "_", input$sc1_pro_inp1, "_", input$sc1_pro_inp2, ".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", dpi = input$sc1_pro_oup.res,
    plot = scProp(sc1conf, sc1meta, input$sc1_pro_inp1, input$sc1_pro_inp2, input$sc1_pro_sub1, input$sc1_pro_sub2, input$sc1_pro_typ, input$sc1_pro_flp, input$sc1_pro_fsz)
    )
  }) 
  
output$sc1_pro_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_pro_typ, "_", input$sc1_pro_inp1, "_", input$sc1_pro_inp2, ".pdf")) },
  content = function(file) {
    ggsave(
    file, device = "pdf", useDingbats = FALSE, bg = "white", onefile = TRUE,
    plot = scProp(sc1conf, sc1meta, input$sc1_pro_inp1, input$sc1_pro_inp2, input$sc1_pro_sub1, input$sc1_pro_sub2, input$sc1_pro_typ, input$sc1_pro_flp, input$sc1_pro_fsz)
    )
  })
  
output$sc1_pro_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_pro_typ, "_", input$sc1_pro_inp1, "_", input$sc1_pro_inp2, ".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white",
    plot = scProp(sc1conf, sc1meta, input$sc1_pro_inp1, input$sc1_pro_inp2, input$sc1_pro_sub1, input$sc1_pro_sub2, input$sc1_pro_typ, input$sc1_pro_flp, input$sc1_pro_fsz)
    )
  }) # End of tab pro



### Tab hea heatmap / dotplot ----

  output$sc1_hea_sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1_hea_sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1_hea_sub2", "Groups to display:", inline = TRUE, choices = sub, selected = sub)
  })
  observeEvent(input$sc1_hea_sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_hea_sub2", label = "Groups to display:", choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc1_hea_sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1_hea_sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc1_hea_sub2", label = "Groups to display:", choices = sub, selected = sub, inline = TRUE)
  })

output$sc1_hea_oupTxt <- renderUI({
  geneList = scGeneList(input$sc1_hea_inp, sc1gene)
  if(nrow(geneList) > 50){
    HTML("More than 50 input genes! Please reduce the gene list!")
  }
})

output$sc1_hea_oup <- renderPlot({
  scBubbHeat(sc1conf, sc1meta, input$sc1_hea_inp, input$sc1_hea_grp, input$sc1_hea_plt, input$sc1_hea_sub1, input$sc1_hea_sub2, "sc1gexpr.h5", sc1gene, input$sc1_hea_scl, input$sc1_hea_row, input$sc1_hea_col, input$sc1_hea_cols, input$sc1_hea_fsz)
})

output$sc1_hea_oup.ui <- renderUI({
  show_progress(imageOutput("sc1_hea_oup", height = pList3[input$sc1_hea_psz]))
})

output$sc1_hea_oup.png <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_hea_plt,"_",input$sc1_hea_grp,".png")) },
  content = function(file) {
    ggsave(
    file, device = "png", bg = "white", height = input$sc1_hea_oup.height, width = input$sc1_hea_oup.width, units = "cm", dpi = input$sc1_hea_oup.res, plot = scBubbHeat(sc1conf, sc1meta, input$sc1_hea_inp, input$sc1_hea_grp, input$sc1_hea_plt, input$sc1_hea_sub1, input$sc1_hea_sub2, "sc1gexpr.h5", sc1gene, input$sc1_hea_scl, input$sc1_hea_row, input$sc1_hea_col, input$sc1_hea_cols, input$sc1_hea_fsz)
    )
})

output$sc1_hea_oup.pdf <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_hea_plt,"_",input$sc1_hea_grp,".pdf")) },
  content = function(file) {
    pdf(file, useDingbats = FALSE, bg = "white", height = input$sc1_hea_oup.height/2.54, width = input$sc1_hea_oup.width/2.54, onefile = TRUE)
    if(!system.file(package="showtext")=="") showtext::showtext_begin()
    print(scBubbHeat(sc1conf, sc1meta, input$sc1_hea_inp, input$sc1_hea_grp, input$sc1_hea_plt, input$sc1_hea_sub1, input$sc1_hea_sub2, "sc1gexpr.h5", sc1gene, input$sc1_hea_scl, input$sc1_hea_row, input$sc1_hea_col, input$sc1_hea_cols, input$sc1_hea_fsz))
    if(!system.file(package="showtext")=="") showtext::showtext_end()
    dev.off()
})

output$sc1_hea_oup.svg <- downloadHandler(
  filename = function() { tolower(paste0("sc1", "_", input$sc1_hea_plt,"_",input$sc1_hea_grp,".svg")) },
  content = function(file) {
    ggsave(
    file, device = "svg", bg = "white", height = input$sc1_hea_oup.height, width = input$sc1_hea_oup.width, units = "cm", 
    plot = scBubbHeat(sc1conf, sc1meta, input$sc1_hea_inp, input$sc1_hea_grp, input$sc1_hea_plt, input$sc1_hea_sub1, input$sc1_hea_sub2, "sc1gexpr.h5", sc1gene, input$sc1_hea_scl, input$sc1_hea_row, input$sc1_hea_col, input$sc1_hea_cols, input$sc1_hea_fsz)
    )
}) # End of tab hea      
       
})


