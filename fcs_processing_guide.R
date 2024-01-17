# install packages
if (!require("BiocManager")) install.packages('BiocManager') 
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("PeacoQC")
BiocManager::install("flowDensity")
devtools::install_url("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/flowType_2.20.1.tar.gz")
BiocManager::install("aya49/flowGraph")

# load packages
library("flowCore")
library("PeacoQC")
library("flowWorkspace")
library("flowDensity")
# library("ggpointdensity") # for easier density scatterplots

# directory to save results in
res_dir <- "/home/maruko/projects/gating"
dir.create(res_dir)

gateplot_dir <- paste0(res_dir, "/gs_plots")
dir.create(gateplot_dir)

# path to raw fcs file
# file from: http://flowrepository.org/id/FR-FCM-ZYXN
# download file here: https://drive.google.com/file/d/1PpSM93GTj9zejVDZzD89_k3sx7Lc-TQl/view?usp=sharing
fcs_path <- paste0(res_dir, "/sangerP2.fcs")

# load fcs file
f <- flowCore::read.FCS(fcs_path)

# explore fcs file
# flowCore::exprs(f) # raw cell x marker matrix
head(flowCore::exprs(f))
dim(flowCore::exprs(f))
flowCore::markernames(f) # marker names (excluding morphology columns)
f@parameters@data


## 1.1 compensate ####
spillover <- flowCore::keyword(f)$SPILL
f <- flowCore::compensate(f, spillover=spillover)

## 1.2 logicle transform ####
transformList <- flowCore::estimateLogicle(f, channels=colnames(spillover))
f <- flowWorkspace::transform(f, transformList)

## 1.3 cleaning; see res_path for plot ####
fmr <- PeacoQC::RemoveMargins(f, channels=c(1,4), output="full")
pQC <- PeacoQC::PeacoQC(fmr[["flowframe"]], channels=colnames(spillover),
                        plot=TRUE, save_fcs=FALSE, report=FALSE,
                        output_directory=res_dir)
f <- pQC[["FinalFF"]]

# clean memory by removing variables
rm("fmr")
rm("pQC")
gc()


## 2 GATING: cell population identification based on a gating strategy ####
# This section is to extract the lymphocytes because we are only interested in those.

## function to rotate 2D frame
## input: 2D matrix and angle
## output: rotated 2D matrix
rotate_data <- function(data, theta=pi/2 - atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])) {
    data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
}

# initialize gating set (containing 1 file)
fs <- flowCore::flowSet(list(f))
# fs <- flowCore::flowSet(list(f1, f2, f3)) # can add more files if needed
gs <- flowWorkspace::GatingSet(fs)


## 2.1 gating all > singlets ####

# get threshold gates
gate_singlets <- flowDensity::deGate(
    f, channel="SSC-W", 
    use.upper=TRUE, upper=TRUE, tinypeak.removal=0.95)
gate_singlets_low <- flowDensity::deGate(
    f, channel="SSC-W", 
    use.percentile=TRUE, percentile=0.0001)

# gate
temp <- flowDensity::flowDensity(
    f, channels=c("FSC-A", "SSC-W"), 
    position=c(NA,FALSE), gates=c(NA, gate_singlets))
fd_singlets <- flowDensity::flowDensity(
    temp, channels=c("FSC-A", "SSC-W"), 
    position=c(NA,TRUE), gates=c(NA, gate_singlets_low))

# plot
png(paste0(gateplot_dir, "/01_all_singlets.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(f)[,"FSC-A"])
plot(d, ylab="", axes=FALSE, 
     main="all events (cells, particles, etc.) > singlets")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(f)[,"SSC-W"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=c(gate_singlets, gate_singlets_low), lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(f, channels=c("FSC-A", "SSC-W"), main="")
lines(fd_singlets@filter)
abline(h=c(gate_singlets, gate_singlets_low), lty="dashed", col="red")
graphics.off()


## 2.2 gating singlets > live ####

# get threshold gates
temp <- flowDensity::flowDensity(
    fd_singlets, channels=c("APC-Cy7-A", "SSC-A"), 
    position=c(NA,F), gates=c(NA,50000))
gate_live <- flowDensity::deGate(
    flowDensity::getflowFrame(temp), channel="APC-Cy7-A")

# gate
fd_live <- flowDensity::flowDensity(
    fd_singlets, channels=c("APC-Cy7-A", "SSC-A"), 
    position = c(FALSE,NA), gates=c(gate_live,NA))

# register gate into gs
gate <- flowCore::rectangleGate(
    filterId="live", 
    "APC-Cy7-A"=c(-1, gate_live), # include all
    "SSC-W"=c(gate_singlets_low, gate_singlets))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="root")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "live")

# plot and save as png
png(paste0(gateplot_dir, "/02_singlets_live.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(
    flowDensity::getflowFrame(fd_singlets))[,"APC-Cy7-A"])
plot(d, ylab="", axes=FALSE, main="singlets > live")
abline(v=gate_live, lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(
    flowDensity::getflowFrame(fd_singlets))[,"SSC-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")

par(mar=c(5,5,1,1))
flowDensity::plotDens(fd_singlets, channels=c("APC-Cy7-A", "SSC-A"), main="")
abline(v=gate_live, lty="dashed", col="red")
graphics.off()

# clean memory by removing variables
rm(fd_singlets)
gc()


## 2.3 gating live > lymphocytes ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    as(flowWorkspace::gs_pop_get_data(
        gs, "live"), Class="list")[[1]] )

# get upper limit gates
gate_ssca_high <- flowDensity::deGate(
    ff, channel="SSC-A", 
    use.percentile=TRUE, percentile=0.9999999)
gate_fsca_high <- flowDensity::deGate(
    ff, channel="FSC-A", 
    use.percentile=TRUE, percentile=0.99999999)

# get threshold gate
gate_fsca <- flowDensity::deGate(
    ff, channel="FSC-A")

# register gate into gs
gate <- flowCore::rectangleGate(
    filterId="lymphocytes", 
    "FSC-A"=c(gate_fsca, gate_fsca_high), # include all
    "SSC-A"=c(0, gate_ssca_high))
node <- flowWorkspace::gs_pop_add(gs, gate, parent="live")
flowWorkspace::recompute(gs)
# flowWorkspace::gs_pop_remove(gs, "lymphocytes")

# plot
png(paste0(gateplot_dir, "/03_live_lymphocytes.png"))
layout(matrix(c(1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 1,3,3,3,3, 0,2,2,2,2), nrow=5))

par(mar=c(0,5,3,1))
d <- density(flowCore::exprs(ff)[,"FSC-A"])
plot(d, ylab="", axes=FALSE, main="live > lymphocytes")
abline(v=c(gate_fsca, gate_fsca_high), lty="dashed", col="red")

par(mar=c(5,0,1,3))
d <- density(flowCore::exprs(ff)[,"SSC-A"])
plot(d$y, d$x, type="l", xlab="", axes=FALSE, main="")
abline(h=c(0, gate_ssca_high), lty="dashed", col="red")

par(mar=c(5,5,1,1))
flowDensity::plotDens(ff, channels=c("FSC-A", "SSC-A"), main="")
abline(v=c(gate_fsca, gate_fsca_high), 
       h=c(0, gate_ssca_high), lty="dashed", col="red")
graphics.off()


## save .fcs file with only lymphocytes ####
ff <- flowWorkspace::cytoframe_to_flowFrame(
    flowWorkspace::gs_pop_get_data(gs, "lymphocytes")[[1]])
flowCore::write.FCS(ff, "sangerP2_new.fcs")

# let's see how many events we removed
dim(flowCore::exprs(f))
dim(flowCore::exprs(ff))

## get thresholds for each marker
marker_indices <- c(7:16)
marker_thres <- list()
# prepare plot
png(paste0(gateplot_dir, "/thresholds.png"), width=900, height=700)
plot_dim <- ceiling(sqrt(length(marker_indices)))
plot_dim_total <- plot_dim^2
layout(matrix(seq_len(plot_dim_total), nrow=plot_dim, ncol=plot_dim))
# for each marker...
for (mi in seq_len(length(marker_indices))) {
    mii <- marker_indices[mi]
    marker_thres[[mi]] <- flowDensity::deGate(
        ff, channel=mii)
    
    # plot
    plot(density(flowCore::exprs(ff)[,mii]), 
         main=ff@parameters@data$name[mii])
    abline(v=marker_thres[[mi]])
}
graphics.off()

## get all possible cell populations and their counts ####
# create one row of the sample x cell population matrix for your current file 

# function to get cell population counts vector
get_counts <- function(
        f, # flowframe
        marker_indices, # vector of numeric indices
        markers, # marker names to use
        thresholds, # list of threshold vectors
        max_markers_per_pop=min(length(thresholds, 4)) # 
) {
    if (!(length(marker_indices) == length(markers) &
          length(markers) == length(thresholds))) {
        print("marker_indices, markers, and thresholds need to be the same length.")
        return()
    }
    paritions_per_marker <- lapply(thresholds, function(t) {
        seq_len(length(t) + 2) - 2
    })
    # enumerate combos
    combos <- expand.grid(paritions_per_marker)
    combos <- combos[apply(combos, 1, function(x) {
        sum(x != -1) <= max_markers_per_pop
    }), , drop=FALSE]
    # create marker labels
    markers <- gsub("[ -+]", "_", markers)
    marker_labels <- apply(combos, 1, function(x) {
        marker_label <- ""
        for (xi in which(x != -1)) {
            marker_label <- paste0(
                marker_label, 
                markers[xi], ifelse(x[xi] == 0, "-", rep("+", x[xi]))
            )
        }
        return(marker_label)
    })
    # get counts
    counts <- apply(combos, 1, function(x) {
        tf <- rep(TRUE, nrow(flowCore::exprs(f)))
        for (xi in which(x != -1)) {
            if (x[xi] == 0) {
                tfx <- flowCore::exprs(f)[, marker_indices[xi]] < thresholds[[xi]][1]
            } else {
                tfx <- flowCore::exprs(f)[, marker_indices[xi]] >= thresholds[[xi]][x[xi]]
            }
            tf <- tf & tfx
        }
        return(sum(tf))
    })
    names(counts) <- marker_labels
    return(counts)
}

ft <- get_counts(
    ff, 
    marker_indices=marker_indices, 
    markers=c("a","b","c","d","e","f","g","h","i","j"),
    thresholds=marker_thres,
    max_markers_per_pop=3
)
ft <- matrix(ft, nrow=1)
rownames(ft) <- "sample_id"

# if you have other files, rbind their `ft`s together
# fts <- Reduce(list(ft1, ft2, ft3), rbind)

# create meta-data
# ft_meta <- data.frame(id=c("sample_id1", "sample_id2", "sample_id3"),
#                       class=c("control", "experiment", "experiment"))
