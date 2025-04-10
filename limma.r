#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
theme_set(theme_classic())
options(stringsAsFactors=FALSE)

# helper variables
cols.group <- c(Baseline="white", Early="#fdd49e", Middle="#e41a1c",
                Late="#7f0000", Transitional="grey40", Recovered="grey80")
shapes.treat <- c(Control=21, "4' FIU"=22)
cols.treat <- c(Control="#b2df8a", "4' FIU"="#33a02c")
cols.reg <- c(Down="#377eb8", None="grey60", Up="#e41a1c")
samp.days <- c(0, 3, 4, 6, 7, 10, 12, 15, 21, 28, 35)
supl.fig <- list()
main.fig <- list()

# helper functions
normalize.counts <- function(counts.matrix, min.samples=3, threshold=20) {
  # how many samples per gene have >threshold counts?
  keepgenes <- rowSums(counts.matrix > threshold)
  # only keep the genes with >= min.samples
  keepgenes <- (keepgenes >= min.samples)
  # subset the counts matrix 
  counts.matrix <- counts.matrix[keepgenes, ]
  # format with edgeR and return
  counts.matrix %>%
    DGEList() %>%
    calcNormFactors()
}

run.pca <- function(norm.matrix, metadata) {
  # take normalized counts and transform to log2CPM before running PCA
  pca <- norm.matrix %>%
         cpm() %>%
         log2() %>%
         t() %>%
         prcomp()
  # pull out the principal component % variance explained
  pcs <- summary(pca)$importance["Proportion of Variance", c("PC1", "PC2")]
  pcs <- round(100*pcs)
  pcs <- paste0(names(pcs), " (", pcs, "%)")
  # format PC matrix with metadata
  pca <- pca$x %>%
         as.data.frame() %>% 
         rownames_to_column("ID") %>%
         select(ID, PC1, PC2) %>%
         left_join(metadata, by="ID")
  # return a list with both
  list(pcs=pca,
       percentvar=pcs)
}

diffexpr <- function(norm.matrix, model.matrix, contrasts, plim=0.05, llim=1) {
  # run differential expression and format output
  x <- norm.matrix %>%
       voom(model.matrix) %>%
       lmFit(design=model.matrix) %>%
       contrasts.fit(contrasts) %>%
       eBayes() %>%
       topTable(number=Inf) %>%
       rownames_to_column("gene") %>%
       dplyr::rename(padj=adj.P.Val,
                     pvalue=P.Value,
                     lfc=logFC,
                     avexpr=AveExpr) %>%
       select(gene, avexpr, lfc, padj, pvalue)
  # extract name of comparison
  x$Comparison <- names(which(contrasts[, 1]==1))
  # add significance and annotate regulation
  x$Significant <- (x$padj < plim & abs(x$lfc) > llim)
  x$Regulation <- "None"
  x[x$Significant & x$lfc > llim, "Regulation"] <- "Up"
  x[x$Significant & x$lfc < -llim, "Regulation"] <- "Down"
  # return full matrix
  return(x)
}

## read inputs -----------------------------------------------------------------
# metadata
meta <- read.csv("samplesheet.csv") %>%
        mutate(Group=factor(Group, levels=names(cols.group)),
               Treatment=factor(Treatment, levels=names(shapes.treat)))

# counts
cmat <- read.csv("counts-thresholded.csv",
                 row.names=1)
# remove control probes
controls <- str_detect(rownames(cmat), "^POS_|^NEG_")
cmat <- cmat[!controls, ]
rm(controls)

# plot samples (supplemental A)
supl.fig$A <- meta %>%
              ggplot(aes(DPI, NHP, fill=Group, shape=Treatment)) +
              geom_line(aes(group=NHP)) +
              geom_point(size=3) +
              scale_shape_manual(values=shapes.treat, ) +
              scale_fill_manual(values=cols.group) +
              facet_wrap(~Treatment, scales="free_y") +
              scale_x_continuous(breaks=samp.days) +
              guides(fill=guide_legend(override.aes=list(pch=21)),
                     shape="none") +
              labs(x="Days postinfection",
                   y=element_blank(),
                   fill=element_blank()) +
              theme(legend.position=c(0.4, 0.5))
supl.fig$A
ggsave("analysis/samples.png",
       units="in", width=5, height=3)

## identify K for k-means of all samples ---------------------------------------
# get log2CPM
logcpm <- cmat %>%
          normalize.counts() %>%
          cpm() %>%
          log2()

# calculate the total within sum of squares (WSS) for K=1 to K=10
# iterate 500 times since this is a stochastic process and plot the st. dev.
km <- lapply(1:100, function(i) {
        elbow <- 1:20 %>%
                 sapply(function(k) {
                   logcpm %>%
                     t() %>%
                     kmeans(centers=k) %>%
                     magrittr::extract2("tot.withinss")
                 })
        data.frame(WSS=elbow, k=1:20)
      }) %>%
        do.call(rbind, .) %>% 
        group_by(k) %>%
        summarise(mean=mean(WSS),
                  stdev=sd(WSS),
                  .groups="drop")
km %>%
  ggplot(aes(k, mean)) +
  geom_ribbon(aes(ymin=mean-stdev, ymax=mean+stdev), fill="grey60") +
  geom_line() +
  geom_point() +
  geom_point(data=filter(km, k==7), pch=21, fill="red", size=3) +
  labs(x="k-means",
       y="Total within sum of squares")
ggsave("analysis/kmeans-wss.png",
       units="in", width=3, height=3)

# calculate K-means when K=7
km <- logcpm %>%
      t() %>%
      kmeans(centers=7) %>%
      magrittr::extract2("cluster")
km <- data.frame(Cluster=km,
                 ID=names(km))

## PCA and clustering ----------------------------------------------------------
# all samples
pca <- cmat %>%
       normalize.counts() %>%
       run.pca(meta)
# PCA of k-means clusters
pca$pcs %>%
  left_join(km, by="ID") %>%
  ggplot(aes(PC1, PC2, shape=Treatment, fill=as.factor(Cluster)))+
  geom_hline(yintercept=0, linetype=2, col="lightgrey") +
  geom_vline(xintercept=0, linetype=2, col="lightgrey") +
  geom_point(size=3) +
  scale_shape_manual(values=c(21, 22))+
  scale_fill_brewer(palette="Paired") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x=pca$percentvar[1],
       y=pca$percentvar[2],
       fill="Cluster",
       title="K-means clusters (K=7)")
ggsave("analysis/pca-kmeans.png",
       units="in", width=4, height=3)
# PCA of groups (main figure A)
main.fig$A <- pca$pcs %>%
              ggplot(aes(PC1, PC2, shape=Treatment, fill=Group)) +
              geom_hline(yintercept=0, linetype=2, col="lightgrey") +
              geom_vline(xintercept=0, linetype=2, col="lightgrey") +
              geom_point(size=3) +
              scale_shape_manual(values=c(21, 22)) +
              scale_fill_manual(values=cols.group) +
              guides(fill=guide_legend(override.aes=list(shape=21))) +
              labs(x=pca$percentvar[1],
                   y=pca$percentvar[2],
                   fill="Group")
main.fig$A
ggsave("analysis/pca-groups.png",
       units="in", width=3.75, height=3)
# clean up
rm(km, pca, logcpm)

# controls only: symptom groups
m <- filter(meta, Treatment=="Control")
pca <- cmat[, m$ID] %>%
       normalize.counts() %>%
       run.pca(m)
pca$pcs %>%
  ggplot(aes(PC1, PC2, fill=Group)) +
  geom_hline(yintercept=0, linetype=2, col="lightgrey") +
  geom_vline(xintercept=0, linetype=2, col="lightgrey") +
  geom_point(size=3, shape=shapes.treat["Control"]) +
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values=cols.group) +
  labs(x=pca$percentvar[1],
       y=pca$percentvar[2],
       fill="Group",
       title="Disease state: controls only")
ggsave("analysis/pca-controls.png",
       units="in", width=4, height=3)

# controls only: symptom groups
m <- filter(meta, Treatment=="4' FIU")
pca <- cmat[, m$ID] %>%
       normalize.counts() %>%
       run.pca(m)
pca$pcs %>%
  ggplot(aes(PC1, PC2, fill=Group)) +
  geom_hline(yintercept=0, linetype=2, col="lightgrey") +
  geom_vline(xintercept=0, linetype=2, col="lightgrey") +
  geom_point(size=3, shape=shapes.treat["4' FIU"]) +
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values=cols.group) +
  labs(x=pca$percentvar[1],
       y=pca$percentvar[2],
       fill="Group",
       title="Disease state: treated only")
ggsave("analysis/pca-treated.png",
       units="in", width=4, height=3)

# clean up
rm(m, pca)

## controls DE analysis: does this look like LASV? -----------------------------
# filter metadata and counts
m <- filter(meta, Treatment=="Control")
c <- cmat[, m$ID] %>%
     normalize.counts()

# format model matrix
x <- droplevels(m$Group)
m <- model.matrix(~0+x)
colnames(m) <- levels(x)

# run contrasts
dexp <- list(diffexpr(c, m, makeContrasts(Early - Baseline, 
                                          levels=colnames(m))),
             diffexpr(c, m, 
                      makeContrasts(Middle - Baseline, 
                                    levels=colnames(m))),
             diffexpr(c, m, 
                      makeContrasts(Late - Baseline, 
                                    levels=colnames(m)))) %>%
        do.call(rbind, .) %>%
        mutate(Comparison=factor(Comparison, 
                                 levels=c("Early", "Middle", "Late")),
               Regulation=factor(Regulation, 
                                 levels=c("Down", "None", "Up")))
write.csv(dexp, "analysis/results-controls.csv",
          row.names=FALSE, na="")
rm(m, c)

# write IPA input
dexp %>%
  select(gene, lfc, padj, Comparison) %>%
  reshape2::melt(id.vars=c("gene", "Comparison")) %>%
  mutate(Header=paste0(Comparison, ".", variable)) %>%
  reshape2::dcast(gene ~ Header, value.var="value") %>%
  write.csv("analysis/ipa-input-controls.csv", row.names=FALSE)

# total DE genes
tot.degs <- dexp %>%
            filter(Significant) %>%
            group_by(Comparison, Regulation) %>%
            summarise(Genes=n(),
                      .groups="drop") %>%
            # add in zeros, if necessary
            right_join(expand.grid(Comparison=c("Baseline", "Early", 
                                                "Middle", "Late"),
                                   Regulation=c("Down", "Up")),
                       by=c("Comparison", "Regulation")) %>%
            replace_na(list(Genes=0)) %>%
              mutate(Treatment="Control")

# volcano plot (supplemental figure C)
topgenes <- dexp %>%
            filter(Significant) %>%
            group_by(Comparison) %>%
            top_n(n=10, wt=abs(lfc)) %>%
            ungroup()
supl.fig$C <- dexp %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
              ggrepel::geom_text_repel(data=topgenes, aes(label=gene),
                                       size=2, force=20, max.overlaps=15) +
              scale_size_manual(values=c(1, 0.5, 1)) +
              scale_color_manual(values=cols.reg) +
              facet_wrap(~Comparison) +
              labs(x="Fold change (log2)",
                   y="FDR-adjusted p-value (-log10)")
supl.fig$C
ggsave("analysis/volcano-controls.png",
       units="in", width=7.5, height=3)
rm(topgenes, dexp)

## 4' FIU treated DE analysis: resolution by endpoint? -------------------------
# filter metadata and counts
m <- meta %>%
     filter(Treatment=="4' FIU")
c <- cmat[, m$ID] %>%
     normalize.counts()

# format model matrix
x <- droplevels(m$Group)
m <- model.matrix(~0+x)
colnames(m) <- levels(x)

# run contrasts by DPI
dexp <- list(diffexpr(c, m, makeContrasts(Middle - Baseline, levels=colnames(m))),
             diffexpr(c, m, makeContrasts(Transitional - Baseline, levels=colnames(m))),
             diffexpr(c, m, makeContrasts(Recovered - Baseline, levels=colnames(m)))) %>%
        do.call(rbind, .) %>%
        mutate(Comparison=factor(Comparison, levels=names(cols.group)),
               Regulation=factor(Regulation, 
                                 levels=c("Down", "None", "Up")))
write.csv(dexp, "analysis/results-treated.csv",
          row.names=FALSE, na="")
rm(x, m, c)

# write IPA input
dexp %>%
  select(gene, lfc, padj, Comparison) %>%
  reshape2::melt(id.vars=c("gene", "Comparison")) %>%
  mutate(Header=paste0(Comparison, ".", variable)) %>%
  reshape2::dcast(gene ~ Header, value.var="value") %>%
  write.csv("analysis/ipa-input-treated.csv", row.names=FALSE)

# total DE genes (main figure B)
main.fig$B <- dexp %>%
              filter(Significant) %>%
              group_by(Comparison, Regulation) %>%
              summarise(Genes=n(),
                        .groups="drop") %>%
              # add in zeros, if necessary
              right_join(expand.grid(Comparison=c("Baseline", "Middle",
                                                  "Transitional", "Recovered"),
                                     Regulation=c("Down", "Up")),
                         by=c("Comparison", "Regulation")) %>%
              replace_na(list(Genes=0)) %>%
              # add in DE genes from control
              mutate(Treatment="4' FIU") %>%
              rbind(tot.degs) %>% 
              mutate(Treatment=factor(Treatment, 
                                      levels=levels(meta$Treatment))) %>%
              # adjust levels to make the transitional-vs-late comparison
              mutate(Comparison=factor(Comparison, levels=names(cols.group),
                                       labels=c("Baseline", "Early", "Middle", 
                                                "Transitional/Late", 
                                                "Transitional/Late", 
                                                "Recovered"))) %>%
              # plot it
              ggplot(aes(Comparison, Genes)) +
              geom_line(aes(linetype=Regulation, 
                            group=paste0(Regulation, Treatment))) +
              geom_point(aes(fill=Regulation, shape=Treatment), size=3) +
              scale_linetype_manual(values=c(Down=2, Up=1)) +
              scale_fill_manual(values=cols.reg) +
              scale_shape_manual(values=shapes.treat) +
              guides(fill=guide_legend(override.aes=list(pch=21))) +
              ylim(0, 150) +
              labs(x="Disease state",
                   y="Total DE genes") +
              theme(axis.text.x=element_text(angle=45, hjust=1))
main.fig$B
ggsave("analysis/degenes-totals.png",
       units="in", width=3.75, height=3)

# volcano plot (supplemental figure D)
topgenes <- dexp %>%
            filter(Significant) %>%
            group_by(Comparison) %>%
            top_n(n=10, wt=abs(lfc)) %>%
            ungroup()
supl.fig$D <- dexp %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
              ggrepel::geom_text_repel(data=topgenes, aes(label=gene),
                                       size=2, force=20, max.overlaps=15) +
              scale_size_manual(values=c(1, 0.5, 1)) +
              scale_color_manual(values=cols.reg) +
              facet_wrap(~Comparison, nrow=1) +
              labs(x="Fold change (log2)",
                   y="FDR-adjusted p-value (-log10)")
supl.fig$D
ggsave("analysis/volcano-treated.png",
       units="in", width=7.5, height=3)

# clean up
rm(dexp, topgenes, tot.degs)

## compare treated vs. controls ------------------------------------------------
# baseline
m <- filter(meta, Group=="Baseline") %>%
     mutate(Treatment=factor(Treatment, levels=names(shapes.treat),
                             labels=c("Control", "Treated")))
c <- cmat[, m$ID] %>%
     normalize.counts()
x <- m$Treatment
m <- model.matrix(~0+x)
colnames(m) <- levels(x)
dexp <- list(Baseline=diffexpr(c, m, 
                               makeContrasts(Treated - Control, 
                                            levels=colnames(m))))
dexp$Baseline <- mutate(dexp$Baseline, Comparison="Baseline")

# middle
m <- filter(meta, Group=="Middle") %>%
     mutate(Treatment=factor(Treatment, levels=names(shapes.treat),
                             labels=c("Control", "Treated")))
c <- cmat[, m$ID] %>%
     normalize.counts()
x <- m$Treatment
m <- model.matrix(~0+x)
colnames(m) <- levels(x)
dexp$Middle <- diffexpr(c, m, 
                        makeContrasts(Treated - Control, 
                                      levels=colnames(m))) %>%
               mutate(Comparison="Middle")

# late/intermediate
m <- filter(meta, Group %in% c("Late", "Transitional")) %>%
     mutate(Treatment=factor(Treatment, levels=names(shapes.treat),
                             labels=c("Control", "Treated")))
c <- cmat[, m$ID] %>%
     normalize.counts()
x <- m$Treatment
m <- model.matrix(~0+x)
colnames(m) <- levels(x)
dexp$LateIntermediate <- diffexpr(c, m, 
                                  makeContrasts(Treated - Control, 
                                                levels=colnames(m))) %>%
                         mutate(Comparison="Transitional/Late")

# combine results
dexp <- dexp %>% 
        do.call(rbind, .) %>%
        mutate(Comparison=factor(Comparison, 
                                 levels=c("Baseline", "Middle", 
                                          "Transitional/Late")),
               Regulation=factor(Regulation, 
                                 levels=c("Down", "None", "Up"),
                                 labels=c("Control", "Neither", "4' FIU")))
write.csv(dexp, "analysis/results-comparison.csv",
          row.names=FALSE, na="")
rm(x, m, c)

# write IPA input
dexp %>%
  select(gene, lfc, padj, Comparison) %>%
  reshape2::melt(id.vars=c("gene", "Comparison")) %>%
  mutate(Header=paste0(Comparison, ".", variable)) %>%
  reshape2::dcast(gene ~ Header, value.var="value") %>% 
  write.csv("analysis/ipa-input-comparison.csv", row.names=FALSE)

# total DE genes (main figure C)
main.fig$C <- dexp %>%
              filter(Significant) %>%
              group_by(Comparison, Regulation) %>%
              summarise(Genes=n(),
                        .groups="drop") %>%
              # add in zeros, if necessary
              right_join(expand.grid(Comparison=levels(dexp$Comparison),
                                     Regulation=c("Control", "4' FIU")),
                         by=c("Comparison", "Regulation")) %>%
              replace_na(list(Genes=0)) %>%
              # plot it
              ggplot(aes(Comparison, Genes, group=Regulation)) +
              geom_col(aes(fill=Regulation), col="black", position="dodge") +
              geom_text(aes(label=Genes, y=Genes+3), 
                        position=position_dodge(width=0.9)) +
              scale_fill_manual(values=cols.treat) +
              labs(x="Disease state",
                   y="Total significantly DE genes",
                   fill="Higher in") +
              theme(legend.position=c(0.25, 0.8))
main.fig$C
ggsave("analysis/degenes-comparison.png",
       units="in", width=3.75, height=3) 

# volcano plot
topgenes <- dexp %>%
            filter(Significant) %>%
            group_by(Comparison, Regulation) %>%
            top_n(n=10, wt=abs(lfc)) %>%
            ungroup()
# full volcano
dexp %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=gene), max.iter=1e8,
                           size=2, force=30, max.overlaps=15) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  scale_color_manual(values=c(cols.treat, Neither="grey60")) +
  facet_wrap(~Comparison, nrow=1) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value (-log10)",
       col="Higher in",
       size="Higher in")
ggsave("analysis/volcano-comparison-all.png",
       units="in", width=7.5, height=3)
# baseline only (supplemental figure B)
supl.fig$B <- dexp %>%
              filter(Comparison=="Baseline") %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
              ggrepel::geom_text_repel(data=filter(topgenes, 
                                                   Comparison=="Baseline"), 
                                       aes(label=gene), max.iter=1e8,
                                       size=2, force=30, max.overlaps=15) +
              scale_size_manual(values=c(1, 0.5, 1)) +
              scale_color_manual(values=c(cols.treat, Neither="grey60")) +
              ylim(0, 6) +
              labs(x="Fold change (log2)",
                   y="FDR-adjusted p-value (-log10)",
                   col="Higher in",
                   size="Higher in") +
              theme(legend.position=c(0.25, 0.8))
supl.fig$B
ggsave("analysis/volcano-comparison-baseline.png",
       units="in", width=2.5, height=3)
# transitional vs. late only (main figure D)
main.fig$D <- dexp %>%
              filter(Comparison=="Transitional/Late") %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(size=Regulation, col=Regulation), alpha=0.8) +
              ggrepel::geom_text_repel(data=filter(topgenes, 
                                                   Comparison=="Transitional/Late"), 
                                       aes(label=gene), max.iter=1e8,
                                       size=2, force=30, max.overlaps=15) +
              scale_size_manual(values=c(1, 0.5, 1)) +
              scale_color_manual(values=c(cols.treat, Neither="grey60")) +
              ylim(0, 6) +
              labs(x="Fold change (log2)",
                   y="FDR-adjusted p-value (-log10)",
                   col="Higher in",
                   size="Higher in")
main.fig$D
ggsave("analysis/volcano-comparison-transitionallate.png",
       units="in", width=3.75, height=3)
rm(dexp, topgenes)

## cell type scoring -----------------------------------------------------------
# load annotations
anno <- read.csv("annotations.csv", na.strings="") %>%
        filter(!is.na(Cell.Type)) %>%
        select(Gene, Cell.Type)

# get log2CPM 
logcpm <- cmat %>%
          normalize.counts() %>%
          cpm() %>%
          log2()

# subset genes -- some may be filtered out in the normalization process
genelist <- intersect(rownames(logcpm), anno$Gene)
anno <- filter(anno, Gene %in% genelist)
logcpm <- logcpm[genelist, ]
rm(genelist)

# get mean log2CPM per cell type
anno <- logcpm %>%
        as.data.frame() %>%
        rownames_to_column("Gene") %>%
        # melt to long
        reshape2::melt(id.vars="Gene",
                       variable.name="ID",
                       value.name="log2CPM") %>%
        # add in cell types
        left_join(anno, by="Gene") %>%
        # compile to mean logCPM by cell type
        group_by(ID, Cell.Type) %>%
        summarise(log2CPM=mean(log2CPM),
                  .groups="drop") %>%
        # add in metadata
        left_join(meta, by="ID")

# plot cell types: controls only
anno %>%
  filter(Treatment=="Control") %>%
  ggplot(aes(Group, log2CPM)) +
  geom_boxplot(aes(fill=Group), alpha=0.5, outliers=FALSE, col="black") +
  geom_jitter(aes(fill=Group), pch=21, width=0.2, height=0, size=0.5) +
  scale_fill_manual(values=cols.group) +
  facet_wrap(~Cell.Type, ncol=5, scales="free_y") +
  labs(x="Disease state",
       y="Mean log2CPM") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))
ggsave("analysis/celltypes-controls.png", 
       units="in", width=7.5, height=5)

# plot cell types: treated only
anno %>%
  filter(Treatment=="4' FIU") %>%
  ggplot(aes(Group, log2CPM)) +
  geom_boxplot(aes(fill=Group), alpha=0.5, outliers=FALSE, col="black") +
  geom_jitter(aes(fill=Group), pch=21, width=0.2, height=0, size=0.5) +
  scale_fill_manual(values=cols.group) +
  facet_wrap(~Cell.Type, ncol=5, scales="free_y") +
  labs(x="Disease state",
       y="Mean log2CPM") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))
ggsave("analysis/celltypes-treated.png", 
       units="in", width=7.5, height=5)

# plot cell types: treated vs. controls
pvalue <- anno %>%
          filter(Group!="Recovered",
                 Group!="Early") %>%
          mutate(Group=factor(Group, 
                              levels=c("Baseline", "Middle", 
                                       "Transitional", "Late"),
                              labels=c("Baseline", "Middle", 
                                       "Transitional/Late",
                                       "Transitional/Late"))) %>%
          group_by(Cell.Type, Group) %>%
          summarise(pvalue=wilcox.test(log2CPM ~ Treatment)$p.value,
                    y.position=max(log2CPM)+0.2,
                    .groups="drop") %>%
          filter(pvalue <= 0.05) %>%
          mutate(xmin=as.integer(Group)-0.25,
                 xmax=as.integer(Group)+0.25) %>%
          rstatix::add_significance(p.col="pvalue", output.col="pvalue")
supl.fig$E <- anno %>%
              filter(Group!="Recovered",
                     Group!="Early") %>%
              mutate(Group=factor(Group, 
                                  levels=c("Baseline", "Middle", 
                                           "Transitional", "Late"),
                                  labels=c("Baseline", "Middle", 
                                           "Transitional/Late",
                                           "Transitional/Late"))) %>%
              ggplot(aes(Group, log2CPM)) +
              geom_boxplot(aes(fill=Treatment), alpha=0.5, 
                           outliers=FALSE, col="black") +
              geom_point(aes(fill=Treatment, group=Treatment), 
                          pch=21, size=0.5, 
                         position=position_jitterdodge(jitter.height=0, 
                                                       jitter.width=0.2)) +
              scale_fill_manual(values=cols.treat) +
              ggpubr::geom_bracket(data=pvalue, aes(label=pvalue)) +
              facet_wrap(~Cell.Type, ncol=5, scales="free_y") +
              # modulator to increase the y-axis limit of plots 
              # with significance bars
              geom_blank(data=pvalue, aes(y=y.position+0.6)) +
              labs(x="Disease state",
                   y="Mean log2CPM") +
              theme(legend.position=c(0.9, -0.1),
                    axis.text.x=element_text(angle=45, hjust=1))
supl.fig$E
ggsave("analysis/celltypes-comparison.png", 
       units="in", width=7.5, height=5)

# subset figure
typelist <- c("B-cells", "NK cells", "Neutrophils", "Macrophages")
pvalue <- filter(pvalue, Cell.Type %in% typelist)
main.fig$E <- anno %>%
              filter(Group!="Recovered",
                     Group!="Early",
                     Cell.Type %in% typelist) %>%
              mutate(Group=factor(Group, 
                                  levels=c("Baseline", "Middle", 
                                           "Transitional", "Late"),
                                  labels=c("Baseline", "Middle", 
                                           "Transitional/Late",
                                           "Transitional/Late"))) %>%
              ggplot(aes(Group, log2CPM)) +
              geom_boxplot(aes(fill=Treatment), alpha=0.5, 
                           outliers=FALSE, col="black") +
              geom_point(aes(fill=Treatment, group=Treatment), 
                         pch=21, size=0.5, 
                         position=position_jitterdodge(jitter.height=0, 
                                                       jitter.width=0.2)) +
              scale_fill_manual(values=cols.treat) +
              ggpubr::geom_bracket(data=pvalue, aes(label=pvalue)) +
              facet_wrap(~Cell.Type, nrow=1, scales="free_y") +
              # modulator to increase the y-axis limit of plots 
              # with significance bars
              geom_blank(data=pvalue, aes(y=y.position+0.6)) +
              labs(x="Disease state",
                   y="Mean log2CPM") +
              theme(legend.position="right",
                    axis.text.x=element_text(angle=45, hjust=1))
main.fig$E
ggsave("analysis/celltypes-highlight.png", 
       units="in", width=7.5, height=3)

## pathways --------------------------------------------------------------------
# what are the main differences between transitional vs late?
pathways <- read.csv("analysis/ipa-output-controls.csv") %>%
            select(Type, Name, Middle, Late) %>%
            mutate(Treatment="Control")
pathways <- read.csv("analysis/ipa-output-treated.csv") %>%
            select(Type, Name, Middle, Transitional) %>%
            rename(Late=Transitional) %>%
            mutate(Treatment="FIU") %>%
            rbind(pathways) %>% 
            reshape2::melt(id.vars=c("Type", "Name", "Treatment"),
                           variable.name="Group",
                           value.name="zscore") %>%
            reshape2::dcast(Type + Name + Group ~ Treatment, 
                            value.var="zscore", fill=0) %>%
            # filter to late/transitional time points only
            filter(Group=="Late")

# what has the largest differences?
pathways %>%
  mutate(Difference=FIU-Control) %>%
  filter(Difference != 0) 

# subset to shortlist with nicknames
main.fig$F <- read.csv("analysis/ipa-output-shortlist.csv") %>%
              left_join(pathways, by="Name") %>%
              # make a "difference" column for sorting
              mutate(Difference=FIU-Control) %>%
              arrange(FIU, abs(Difference)) %>%
              mutate(Nickname=factor(Nickname, levels=Nickname)) %>%
              reshape2::melt(id.vars=c("Nickname", "Type"),
                             measure.vars=c("Control", "FIU"),
                             variable.name="Treatment",
                             value.name="zscore") %>%
              mutate(Treatment=factor(Treatment, levels=c("Control", "FIU"),
                                      labels=c("Control, late", 
                                               "4' FIU, transitional")),
                     Type=factor(Type, 
                                 levels=c("Canonical", "Functional", 
                                          "Upstream"),
                                 labels=c("Pathways", "Pathways", 
                                          "Regulators"))) %>%
              arrange(Type) %>%
              mutate(Nickname=factor(Nickname, levels=unique(Nickname))) %>%
              ggplot(aes(zscore, Nickname)) +
              geom_vline(xintercept=0, linetype=2, col="grey") +
              geom_line(aes(group=Nickname, linetype=Type)) +
              geom_point(aes(fill=Treatment, shape=Treatment), size=2) +
              scale_fill_manual(values=as.character(cols.treat)) +
              scale_shape_manual(values=as.numeric(shapes.treat)) +
              labs(x="z-score",
                   y=element_blank(),
                   fill=element_blank(),
                   shape=element_blank(),
                   linetype=element_blank()) +
              theme(legend.position="right")
main.fig$F
ggsave("analysis/ipa-late.png",
       units="in", width=7.5, height=3)

## assemble figures ------------------------------------------------------------
# main figure
x <- cowplot::plot_grid(plotlist=main.fig[c("A", "B", "C", "D")], labels="AUTO")
cowplot::plot_grid(x, main.fig$E, main.fig$F, labels=c(NA, "E", "F"), 
                   ncol=1, rel_heights=c(2, 1, 1))
ggsave("analysis/figure-main.png", units="in", width=7.5, height=12)

# supplemental figure requires more paneling
x <- cowplot::plot_grid(supl.fig$A, supl.fig$B, nrow=1, rel_widths=c(2, 1), 
                        labels="AUTO")
cowplot::plot_grid(x, supl.fig$C, supl.fig$D, supl.fig$E, 
                   ncol=1, labels=c(NA, "C", "D", "E"),
                   rel_heights=c(3, 3, 3, 4))
ggsave("analysis/figure-supplemental.png", units="in", width=7.5, height=13)

# clean up
rm(x)

## fin -------------------------------------------------------------------------
sessionInfo()
