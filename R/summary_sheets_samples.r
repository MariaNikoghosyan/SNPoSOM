pipeline.summarySheetsSamples <- function(env)
{
  if(ncol(env$indata) >= 1000) return()
    
  dir.create(env$output.paths["Summary Sheets Samples"], showWarnings=FALSE)

  #### Summary Sheets ####

  ylim.max <- 0

  for (m in 1:ncol(env$indata))
  {
    h <- hist(env$p.g.m[,m], bre=20, plot=FALSE)
    y.max <- max(h$density)

    if (y.max > ylim.max)
    {
      ylim.max <- y.max
    }
  }

  util.info("Writing:", file.path(env$output.paths["Summary Sheets Samples"], "*.pdf"))

  for (m in 1:ncol(env$indata))
  {
    basename <- paste(make.names(make.unique(colnames(env$indata))[m]), ".pdf", sep="")
    pdf(file.path(env$output.paths["Summary Sheets Samples"], basename), 29.7/2.54, 21/2.54, useDingbats=FALSE)

    ## Global Sheet
    #layout(matrix(c(1,2,0,1,3,4,5,5,6,7,7,8), 3, 4), widths=c(1,1,2,2), heights=c(2,1,1))
    layout(matrix(c(1,2,0,1,3,4,5,5,6,7,7,6), 3, 4), widths=c(1,1,2,2), heights=c(2,1,1))
    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

    text(0.1, 0.94, colnames(env$indata)[m] , cex=3, adj=0)
    text(0.1, 0.8, "Global Summary" , cex=1.8, adj=0)
    #text(0.1, 0.7,  paste("%DE =", round(env$perc.DE.m[colnames(env$indata)[m]], 2)), adj=0)

    # Infor related to genes
    # all.fdr.genes <- which(env$fdr.g.m[,m] < 0.2)
    # plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 0)
    # minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 0)
    # 
    # text(0.1, 0.65, paste("# genes with fdr < 0.2  =",
    #                       length(all.fdr.genes), " (",
    #                       length(plus.fdr.genes), "+ /",
    #                       length(minus.fdr.genes), " -)"), adj=0)
    # 
    # all.fdr.genes <- which(env$fdr.g.m[,m] < 0.1)
    # plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 0)
    # minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 0)
    # 
    # text(0.1, 0.6, paste("# genes with fdr < 0.1  =",
    #                      length(all.fdr.genes)," (",
    #                      length(plus.fdr.genes), "+ /",
    #                      length(minus.fdr.genes), " -)"), adj=0)
    # 
    # all.fdr.genes <- which(env$fdr.g.m[,m] < 0.05)
    # plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 0)
    # minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 0)
    # 
    # text(0.1, 0.55, paste("# genes with fdr < 0.05  =",
    #                       length(all.fdr.genes)," (",
    #                       length(plus.fdr.genes), "+ /",
    #                       length(minus.fdr.genes), " -)"), adj=0)
    # 
    # all.fdr.genes <- which(env$fdr.g.m[,m] < 0.01)
    # plus.fdr.genes <- which(env$indata[all.fdr.genes, m] > 0)
    # minus.fdr.genes <- which(env$indata[all.fdr.genes, m] <= 0)
    # 
    # text(0.1, 0.5, paste("# genes with fdr < 0.01 =",
    #                      length(all.fdr.genes)," (",
    #                      length(plus.fdr.genes), "+ /",
    #                      length(minus.fdr.genes), " -)"), adj=0)

    
    if(env$preferences$group.maf)
    {
      
      text(0.1, 0.65, paste("# SNPs with homozygous minor allele =",
                            env$indata.temp[[m]]$minor), adj=0)
      text(0.1, 0.6, paste("# SNPs with homozygous major allele =",
                           env$indata.temp[[m]]$major), adj=0)
      text(0.1, 0.55, paste("# Heterozygous SNPs =",
                            env$indata.temp[[m]]$hetero), adj=0)

      text(0.1, 0.35,  paste("<MAF> =", 
                             round(env$indata.temp[[m]]$minor/sum(env$indata.temp[[m]]$minor, env$indata.temp[[m]]$major, env$indata.temp[[m]]$hetero), 2)), adj=0)
    }else
    {
      minor.allele <- which(env$indata.temp[,m] == 2)
      ref.allele <- which(env$indata.temp[,m] == 0)
      hetero.allele <- which(env$indata.temp[,m] == 1)
      
      text(0.1, 0.65, paste("# SNPs with homozygous minor allele =",
                            length(minor.allele)), adj=0)
      text(0.1, 0.6, paste("# SNPs with homozygous major allele =",
                           length(ref.allele)), adj=0)
      text(0.1, 0.55, paste("# Heterozygous SNPs =",
                            length(hetero.allele)), adj=0)
      
      text(0.1, 0.35,  paste("<MAF> =", round(sum(env$indata.temp[,m])/(nrow(env$indata.temp)*2), 2)), adj=0)
    }

    #text(0.1, 0.35,  paste("<FC> =", round(mean(env$indata[,m]), 2)), adj=0)
    
    text(0.1, 0.3,  paste("<p-value> =", round(10 ^ mean(log10(env$p.g.m[,m])), 2)), adj=0)
    text(0.1, 0.25, paste("<fdr> =", round(mean(env$fdr.g.m[,m]), 2)), adj=0)

    
    # portrait
    
    par(mar=c(2,3,3,1))

    image(matrix(env$metadata[,m], env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom),
          axes=FALSE, col = env$color.palette.portraits(1000), main="Portrait", cex.main=1.5)

    axis(1,
         seq(0, 1, length.out=env$preferences$dim.1stLvlSom/10+1),
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out=env$preferences$dim.1stLvlSom/10)),
         cex.axis=1.0)

    axis(2,
         seq(0, 1, length.out=env$preferences$dim.1stLvlSom/10+1),
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out=env$preferences$dim.1stLvlSom/10)),
         cex.axis=1.0, las=1)

    box()
    
    
    # top 100 differentially expressed genes
    
    n.genes <- min( 100, nrow(env$indata) )
    n.map <- matrix(0,env$preferences$dim.1stLvlSom,env$preferences$dim.1stLvlSom)
    set.genes <- names(sort(env$p.g.m[,m]))[1:n.genes]
    gs.nodes <- env$som.result$feature.BMU[set.genes]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
    
    lim <- c(1,env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
    colr <- env$color.palette.heatmaps(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                                               max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                                               999 + 1]
    
    par(mar=c(2,3,3,1))
    plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
         xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="Top 100 SNPs", cex.main=1.5,
         cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
    
    axis(1,
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out=env$preferences$dim.1stLvlSom/10)),
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out=env$preferences$dim.1stLvlSom/10)),
         cex.axis=1.0)
    
    axis(2,
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out=env$preferences$dim.1stLvlSom/10)),
         c(1, seq(10, env$preferences$dim.1stLvlSom, length.out=env$preferences$dim.1stLvlSom/10)),
         cex.axis=1.0, las=1)
    
    box()
    
    
    # volcano plot
    
    par(mar=c(3,3,3,1))
    
    plot( env$indata[set.genes,m], -log10(env$p.g.m[set.genes,m]), xlim=c(-1,1)*max(abs(env$indata[set.genes,m])), xlab="", ylab="", pch=16,col="gray30", las=1)
      mtext("log FC", 1, line=2, cex=.6)
      mtext("-log10(p)", 2, line=2, cex=.6)
      abline(v=0,lty=3,col="gray80",lwd=1.5)
      abline(h=-log10(0.05),lty=3,col="gray80",lwd=1.5)
    
      
    # differentially expressed genes list
    
    n.genes <- 20
    
    de.genes <- names(sort(env$p.g.m[,m]))
    de.genes <- de.genes[ which(env$indata[de.genes,m]>0) ][1:n.genes]
    
    
    # de.genes.labels <- env$gene.info$names[de.genes]
    # de.genes.labels[which(de.genes.labels=="")] <- de.genes[which(de.genes.labels=="")]
    
    # change to only snps 
    de.genes.labels <- de.genes

    
    par(mar=c(0,0,0,0))

    #x.coords <- c(0, 0.06, 0.2, 0.28, 0.36, 0.44, 0.52)
    x.coords <- c(0, 0.06, 0.2, 0.28, 0.36, 0.44, 0.52)
    y.coords <- seq(0.75, 0.4, length.out=n.genes)

    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

    # remove log(FC) values add gene name
    text(0, 0.88, "SNPs with overrepresented minor allele", cex=1.8, adj=0)
    text(x.coords, rep(c(0.82, 0.80), 4)[1:7],
         c("Rank", "ID", "Gene", "p-value", "fdr", "Metagene", "Description"),
         cex=1, adj=0)
    #text(x.coords[1], 0.77, "Overexpressed", cex=0.8, adj=0, font=3)
    
    text(x.coords[1], y.coords, c(1:n.genes), adj=0)
    text(x.coords[2], y.coords, de.genes.labels, cex=0.6, adj=0)
    rect(x.coords[3]-0.02, y.coords[1]+0.01, 1, 0, border="white", col="white")
    text(x.coords[3], y.coords,env$gene.info$names[de.genes.labels], cex=0.6, adj=0)
    text(x.coords[4], y.coords, format(env$p.g.m[de.genes, m], digits=1), cex=0.6, adj=0)
    text(x.coords[5], y.coords, format(env$fdr.g.m[de.genes, m], digits=1), cex=0.6, adj=0)
    text(x.coords[6], y.coords, env$gene.info$coordinates[de.genes], cex=0.6, adj=0)
    text(x.coords[7], y.coords, env$gene.info$descriptions[de.genes], cex=0.6, adj=0)

    
    # ###
    # text(x.coords[1], y.coords, c(1:n.genes), adj=0)
    # text(x.coords[2], y.coords, de.genes.labels, cex=0.6, adj=0)
    # rect(x.coords[3]-0.02, y.coords[1]+0.01, 1, 0, border="white", col="white")
    # text(x.coords[3], y.coords, round(env$indata[de.genes, m], 2), cex=0.6, adj=0)
    # text(x.coords[4], y.coords, format(env$p.g.m[de.genes, m], digits=1), cex=0.6, adj=0)
    # text(x.coords[5], y.coords, format(env$fdr.g.m[de.genes, m], digits=1), cex=0.6, adj=0)
    # text(x.coords[6], y.coords, env$gene.info$coordinates[de.genes], cex=0.6, adj=0)
    # text(x.coords[7], y.coords, env$gene.info$descriptions[de.genes], cex=0.6, adj=0)
    # ###
    
    # 
    # 
    # de.genes <- names(sort(env$p.g.m[,m]))
    # de.genes <- de.genes[ which(env$indata[de.genes,m]<0) ][1:n.genes]
    # 
    # # de.genes.labels <- env$gene.info$names[de.genes]
    # # de.genes.labels[which(de.genes.labels=="")] <- de.genes[which(de.genes.labels=="")]
    # 
    # de.genes.labels <- de.genes
    # 
    # y.coords <- seq(0.35, 0.02, length.out=n.genes)
    # 
    # text(x.coords[1], 0.37, "Underexpressed", cex=0.8, adj=0, font=3)
    # 
    # text(x.coords[1], y.coords, c(1:n.genes), adj=0)
    # text(x.coords[2], y.coords, de.genes.labels, cex=0.6, adj=0)
    # rect(x.coords[3]-0.02, y.coords[1]+0.01, 1, 0, border="white", col="white")
    # text(x.coords[3], y.coords,env$gene.info$names[de.genes.labels], cex=0.6, adj=0)
    # text(x.coords[4], y.coords, format(env$p.g.m[de.genes, m], digits=1), cex=0.6, adj=0)
    # text(x.coords[5], y.coords, format(env$fdr.g.m[de.genes, m], digits=1), cex=0.6, adj=0)
    # text(x.coords[6], y.coords, env$gene.info$coordinates[de.genes], cex=0.6, adj=0)
    # text(x.coords[7], y.coords, env$gene.info$descriptions[de.genes], cex=0.6, adj=0)
    
    
    # p-value histogram
    
    # par(mar=c(3,6,2,6))
    # 
    # hist(env$p.g.m[,m], bre=20, freq=FALSE, xlab="p-value", ylab="", main="p-values",
    #      ylim=c(0,ylim.max), las=1, cex.main=1.5, cex.lab=1, cex.axis=1)
    # 
    # box()
    # mtext("Density", side=2, line=3, cex=1)
    # mtext("FDR", side=4, line=3, cex=1)
    # mtext(paste("%DE =", round(env$perc.DE.m[colnames(env$indata)[m]] ,2)), line=-1.2, cex=0.5)
    # 
    # abline(h=env$n.0.m[colnames(env$indata)[m]], col="gray", lwd=2)
    # 
    # par(new=TRUE)
    # plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
    # axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)
    # o <- order(env$p.g.m[,m])
    # lines(env$p.g.m[o,m], env$fdr.g.m[o,m], lty=3, lwd=3)
    # 
    # legend("topright", c("p", expression(env$eta[0]), "fdr"),
    #        col=c("black","gray","black"), lty=c(1,1,3),
    #        lwd=c(1,1,2), cex=0.7)

    #  manhattan snps plot
    par(mar=c(5,5,0,3))
    chr.mat <- data.frame(ids = names(env$gene.info$chr.name), chr = env$gene.info$chr.name, position = env$gene.info$chr.start)
    chr.mat[(is.na(chr.mat$chr) | is.na(chr.mat$position)),] <- ""
    chr.mat <- chr.mat[-which(chr.mat$chr == "" | chr.mat$position == ""),]
    #filter SNPs on autosomal chromosoms
    chr.mat <- chr.mat[which(chr.mat$chr %in% as.character(1:23)),]
    ps <- env$p.g.m[chr.mat$ids,m]
    ps[which(ps == 0)] <- min(ps[which(ps != 0)])/2
    d = data.frame(CHR = as.numeric(chr.mat$chr), BP = as.numeric(chr.mat$position),
                   pos = NA, index = NA, 
                   SNP =  as.character(chr.mat$ids), stringsAsFactors = FALSE,
                   logp =-log10(ps))

    d <- d[order(d$CHR, d$BP), ]
    col = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
            "#FCCDE5", '#D9D9D9', "#BC80BD","#CCEBC5", "#FFED6F", "#8DD3C7",
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
            "#FCCDE5", '#D9D9D9', "#BC80BD")
    d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                               d$CHR, length))
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
      d$pos = d$BP
      xlabel = paste("Chromosome", unique(d$CHR), "position")
    }  else {
      lastbase = 0
      ticks = NULL
      for (i in unique(d$index)) {
        if (i == 1) {
          d[d$index == i, ]$pos = d[d$index == i, ]$BP
        }
        else {
          lastbase = lastbase + max(d[d$index == (i - 
                                                    1), "BP"])
          d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
            min(d[d$index == i, "BP"]) + 1
          d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
            lastbase
        }
      }
      ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
      xlabel = "Chromosome"
      labs <- unique(d$CHR)
    }
    
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    # plot
    plot(NA,xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
         las = 1, pch = 20, xlim = c(xmin, xmax), 
         ylim = c(0,ceiling(max(d$logp))), 
         xlab = xlabel, ylab = expression(-log[10](italic(p))))
    #mtext(expression(-log[10](italic(p))), 2, line=2, cex=.6)
    # add chromosomes
    if (nchr == 1) {
      axis(1)
    }else 
    {
      for (chr in 1:length(labs))
        {
        axis(1, at = ticks[chr], labels = labs[chr], cex.axis = 0.7, lwd = 0,
             line = ifelse(chr%%2 > 0, 0, -0.5))
      }
      axis(1, at = ticks, labels = NA)
      
    }
    # add colored porints
    col = rep_len(col, max(d$index))
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"],col = col[icol], pch = 20)
      icol = icol + 1
    }
    # line for significant level
    abline(h = -log10(0.05), col = "blue")
    # annotate significant SNPs
    # topHits = subset(d,  logp >= -log10(0.05))
    # if(nrow(topHits) > 0)
    # {
    #   text(logp ~ pos, labels = topHits$SNP,data = topHits, cex=0.5, font=1, pos=1)
    # }
   
    box()
    
    if (env$preferences$activated.modules$geneset.analysis)
    {
      # differentially expressed gene sets list
      
      n.sets <- 20

      top.gs.score <- sort(env$spot.list.samples[[m]]$GSZ.score, decreasing=TRUE)[1:n.sets]
      top.gs.p <- env$spot.list.samples[[m]]$GSZ.p.value[names(top.gs.score)]

      par(mar=c(0,0,0,0))

      #x.coords <- c(0, 0.1, 0.18, 0.30, 0.39, 0.47)
      #x.coords <- c(0, 0.09, 0.15, 0.25, 0.34, 0.40, 0.47)
      x.coords <- c(0, 0.05, 0.12, 0.19, 0.25, 0.34,0.40, 0.47)
      y.coords <- seq(0.75, 0.4, length.out=n.sets)

      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

      text(0, 0.88, "Overrepresented SNP sets", cex=1.8, adj=0)
      text(x.coords, rep(c(0.82, 0.80), 4)[1:8],
           c("Rank", "SSZ", "p-value","S.score","p-value", "#all", "Geneset", ""), cex=1, adj=0)
      #text(x.coords[1], 0.77, "Overexpressed", cex=0.8, adj=0, font=3)

      text(x.coords[1], y.coords, c(1:n.genes), adj=0)
      text(x.coords[2], y.coords, round(top.gs.score, 2), cex=0.6, adj=0)
      text(x.coords[3], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)
      
      text(x.coords[4], y.coords, format(env$spot.list.samples[[m]]$SNP.score[names(top.gs.score)], digits = 1), cex=0.6, adj=0)
      text(x.coords[5], y.coords, format(env$spot.list.samples[[m]]$SNP.p.value[names(top.gs.score)], digits = 3), cex=0.6, adj=0)

      text(x.coords[6], y.coords, sapply(env$gs.def.list[names(top.gs.score)],
                                         function(x) { length(x$Snps) }), cex=0.6, adj=0)

      text(x.coords[7], y.coords, sapply(env$gs.def.list, function(x) { x$Type })[names(top.gs.score)], cex=0.6, adj=0)
      text(x.coords[8], y.coords, names(top.gs.score), cex=0.6, adj=0)

      # top.gs.score <- sort(env$spot.list.samples[[m]]$GSZ.score, decreasing=FALSE)[1:n.sets]
      # top.gs.p <- env$spot.list.samples[[m]]$GSZ.p.value[names(top.gs.score)]
      # 
      # y.coords <- seq(0.35, 0.02, length.out=n.sets)
      # 
      # text(x.coords[1], 0.37, "Underexpressed", cex=0.8, adj=0, font=3)
      # text(x.coords[1], y.coords, c(1:n.genes), adj=0)
      # text(x.coords[2], y.coords, round(top.gs.score, 2), cex=0.6, adj=0)
      # text(x.coords[3], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)
      # 
      # text(x.coords[4], y.coords, sapply(env$gs.def.list[names(top.gs.score)],
      #                                    function(x) { length(x$Genes) }), cex=0.6, adj=0)
      # 
      # text(x.coords[5], y.coords, sapply(env$gs.def.list, function(x) { x$Type })[names(top.gs.score)], cex=0.6, adj=0)
      # text(x.coords[6], y.coords, names(top.gs.score), cex=0.6, adj=0)

      
      # p-value histogram
      
      # p <- env$spot.list.samples[[m]]$GSZ.p.value
      # 
      # fdrtool.result <- suppressWarnings(fdrtool(p, statistic="pvalue", verbose=FALSE, plot=FALSE))
      # fdr.spot.list.samples <- fdrtool.result$lfdr
      # Fdr.spot.list.samples <- fdrtool.result$qval
      # 
      # n.0.spot.list.samples <- fdrtool.result$param[1,"eta0"]
      # perc.DE.spot.list.samples <- 1 - n.0.spot.list.samples
      # 
      # par(mar=c(3,6,2,6))
      # 
      # hist(p, bre=20, freq=FALSE, xlab="p-value", ylab="", main="p-values",
      #      las=1, cex.main=1.5, cex.lab=1, cex.axis=1)
      # 
      # box()
      # mtext("Density", side=2, line=3, cex=1)
      # mtext("FDR", side=4, line=3, cex=1)
      # mtext(paste("%DE =", round(perc.DE.spot.list.samples ,2)), line=-1.2, cex=0.5)
      # 
      # abline(h=n.0.spot.list.samples , col="gray", lwd=2)
      # 
      # par(new=TRUE)
      # plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
      # axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)
      # o = order(p)
      # lines(p[o], Fdr.spot.list.samples[o], lty=2, lwd=2)
      # lines(p[o], fdr.spot.list.samples[o], lty=3, lwd=3)
      # 
      # legend("topright",
      #        c("p", expression(env$eta[0]), "Fdr", "fdr"),
      #        col=c("black","gray","black","black"),
      #        lty=c(1,1,2,3), lwd=c(1,1,1,2), cex=0.7)
      

    }
   
    dev.off()
    }

  }

