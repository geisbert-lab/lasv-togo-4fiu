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
rm(dexp, topgenes)

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

## assemble figures ------------------------------------------------------------
# main figure is a 2x2 grid
cowplot::plot_grid(plotlist=main.fig, labels="AUTO")
ggsave("analysis/figure-main.png", units="in", width=7.5, height=6)



# plot samples (supplemental A)
supl.fig$A <- meta %>%
  ggplot(aes(DPI, NHP, fill=Group, shape=Treatment)) +
  geom_line(aes(group=NHP)) +
  geom_point(size=3) +
  scale_shape_manual(values=shapes.treat, ) +
  scale_fill_manual(values=cols.group) +
  facet_wrap(~Treatment, scales="free_y", ncol=1) +
  scale_x_continuous(breaks=samp.days) +
  guides(fill=guide_legend(ncol=2, override.aes=list(pch=21)),
         shape="none") +
  labs(x="Days postinfection",
       y=element_blank(),
       fill=element_blank()) +
  theme(legend.position=c(0.8, 0.77))
supl.fig$A

# supplemental figure requires more paneling
x <- cowplot::plot_grid(supl.fig$A, supl.fig$B, nrow=1, rel_widths=c(2, 1), 
                        labels="AUTO")
cowplot::plot_grid(x, supl.fig$C, supl.fig$D, ncol=1, labels=c(NA, "C", "D"))
ggsave("analysis/figure-supplemental.png", units="in", width=7.5, height=9)

## fin -------------------------------------------------------------------------
sessionInfo()
