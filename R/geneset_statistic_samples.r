pipeline.genesetStatisticSamples <- function(env)
{

  ### perform GS analysis ###
  
  # divide env$gs.def.list into 2 parts for genes and snps 
  # option 2 make enrichment only for SNPs
  # ind.snps <- c()
  #   
  # for (i in seq_along(env$gs.def.list)) 
  #   {
  #   if(env$gs.def.list[[i]]$Type == "Clinical trait" | env$gs.def.list[[i]]$Type == "Trait" | env$gs.def.list[[i]]$Type == "Disease")
  #   {
  #     ind.snps <- c(ind.snps, i)
  #   }
  # }
  # 
  # gs.def.list.snp <- env$gs.def.list[ind.snps] 
  # gs.def.list.gene <- env$gs.def.list[-ind.snps] 
  
  ## GS for genes ##
  # env$indata.ensID.m <- env$indata[env$gene.info$ensembl.mapping[,1],]
  # env$indata.ensID.m <- do.call(rbind, by(env$indata.ensID.m, env$gene.info$ensembl.mapping[,2], colMeans))
  # 
  # mean.ex.all <- colMeans( env$indata.ensID.m )
  # sd.ex.all <- apply( env$indata.ensID.m, 2, sd )
  # 
  # gs.null.list <- list()
  # for (i in seq_along(gs.def.list.gene))
  # {
  #   gs.null.list[[i]] <-
  #     list(Genes=sample(unique(env$gene.info$ensembl.mapping$ensembl_gene_id), length(gs.def.list.gene[[i]]$Genes)))
  # }
  # 
  # null.scores <- sapply( gs.null.list, Sample.GSZ, env$indata.ensID.m, mean.ex.all, sd.ex.all )
  # null.culdensity <- ecdf(abs(unlist(null.scores)))
  # 
  # env$samples.GSZ.scores <- t( sapply( gs.def.list.gene, Sample.GSZ, env$indata.ensID.m, mean.ex.all, sd.ex.all ) )
  # 
  # env$spot.list.samples <- lapply(seq_along(env$spot.list.samples) , function(m)
  # {
  #   x <- env$spot.list.samples[[m]]
  # 
  #   x$GSZ.score <- env$samples.GSZ.scores[,m]
  #   x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
  #   names(x$GSZ.p.value) <- names(x$GSZ.score)
  # 
  #   return(x)
  # })
  # names(env$spot.list.samples) <- colnames(env$indata)
  
  ## GS for snps ##

  mean.ex.all.snp <- colMeans( env$indata )
  sd.ex.all.snp <- apply( env$indata, 2, sd )
  
  gs.null.list.snp <- list()
  for (i in seq_along(env$gs.def.list))
  {
    gs.null.list.snp[[i]] <-
      list(Snps=sample(rownames(env$indata), length(env$gs.def.list[[i]]$Snps)))
  }
  
  null.scores <- sapply( gs.null.list.snp, Sample.SSZ, env$indata, mean.ex.all.snp, sd.ex.all.snp )
  null.culdensity <- ecdf(abs(unlist(null.scores)))
  
  env$samples.GSZ.scores <- t( sapply( env$gs.def.list, Sample.SSZ, env$indata, mean.ex.all.snp, sd.ex.all.snp ) )
  
  
  env$spot.list.samples <- lapply(seq_along(env$spot.list.samples) , function(m)
  {
    x <- env$spot.list.samples[[m]]
    
    x$GSZ.score <- env$samples.GSZ.scores[,m]
    x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
    names(x$GSZ.p.value) <- names(x$GSZ.score)
    
    return(x)
  })
  names(env$spot.list.samples) <- colnames(env$indata)

  
  ### perform SNPs analysis ###

  for (i in seq_along(env$spot.list.samples)) 
    {
    env$spot.list.samples[[i]]$SNP.score <- sapply(env$gs.def.list, function(x)
      {
      x <- sum(env$indata.temp[x$Snps,i])/(length(x$Snps)*2)
      return(x)
    })
  }

  # create null distribution
  for (i in seq_along(env$spot.list.samples)) 
    {
    null.snp.score <- sapply(gs.null.list.snp, function(x)
      {
      x <- sum(env$indata.temp[x$Snps, i])/(length(x$Snps)*2)
    })
    null.snp.culdensity <- ecdf(abs(unlist(null.snp.score)))
    env$spot.list.samples[[i]]$SNP.p.value <-  null.snp.culdensity(abs(env$spot.list.samples[[i]]$SNP.score))
    names(env$spot.list.samples[[i]]$SNP.p.value) <- names(env$spot.list.samples[[i]]$SNP.score)
  }
  
  return(env)
}
