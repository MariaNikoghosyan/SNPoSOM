pipeline.genesetStatisticModules <- function(env)
{
  # spot.fisher.p <- function(spot)
  # {
  #   spot$Fisher.p <- GeneSet.Fisher(unique(env$gene.info$ensembl.mapping$ensembl_gene_id[ which(env$gene.info$ensembl.mapping[,1]%in%spot$genes) ]),
  #                                   unique(env$gene.info$ensembl.mapping$ensembl_gene_id), env$gs.def.list, sort=TRUE)
  # 
  #   return(spot)
  # }
  
  ## modify for snps
  # divide env$gs.def.list into 2 parts for genes and snps 
  # ind.snps <- c()
  # 
  # for (i in seq_along(env$gs.def.list))
  # {
  #   if(env$gs.def.list[[i]]$Type == "Clinical trait" | env$gs.def.list[[i]]$Type == "Trait" | env$gs.def.list[[i]]$Type == "Disease")
  #   {
  #     ind.snps <- c(ind.snps, i)
  #   }
  # }
  # 
  # gs.def.list.snp <- env$gs.def.list[ind.snps]
  # gs.def.list.gene <- env$gs.def.list[-ind.snps]
  
  # option 2 make enrichment
  spot.fisher.p <- function(spot)
  {
    #genes
    # spot$Fisher.p <- GeneSet.Fisher(unique(env$gene.info$ensembl.mapping$ensembl_gene_id[ which(env$gene.info$ensembl.mapping[,1]%in%spot$genes) ]),
    #                                 unique(env$gene.info$ensembl.mapping$ensembl_gene_id), gs.def.list.gene, sort=TRUE)
    # #snps
    # Fisher.p.snp <- GeneSet.Fisher.Snp(unique(names(env$gene.info$ids)[ which(names(env$gene.info$ids)%in%spot$genes) ]),
    #                                     names(env$gene.info$ids), gs.def.list.snp, sort=TRUE)
    # spot$Fisher.p <- c(spot$Fisher.p, Fisher.p.snp)

    spot$Fisher.p <- GeneSet.Fisher.Snp(unique(names(env$gene.info$ids)[ which(names(env$gene.info$ids)%in%spot$genes) ]),
                                        names(env$gene.info$ids),  env$gs.def.list, sort=TRUE)
    
    return(spot)
  }

  list.ids <- unique(names(env$gene.info$ids)[ which(names(env$gene.info$ids)%in%env$spot.list.overexpression$spots$C$genes) ])
  all.ids <-  names(env$gene.info$ids)
  gs.def.list <- env$gs.def.list
  x <- gs.def.list$`actin binding`
  

  
  
  
  
  # overexpression spot
  env$spot.list.overexpression$spots <- lapply( env$spot.list.overexpression$spots, spot.fisher.p)
  
  # underexpression spot
  env$spot.list.underexpression$spots <- lapply( env$spot.list.underexpression$spots, spot.fisher.p)
  # correlation spot
  if (length(env$spot.list.correlation$spots) > 0)
  {
    env$spot.list.correlation$spots <- lapply( env$spot.list.correlation$spots, spot.fisher.p)
  }
  # kmeans spot
  env$spot.list.kmeans$spots <- lapply( env$spot.list.kmeans$spots, spot.fisher.p)
  
  if (length(unique(env$group.labels)) > 1)
  {
    env$spot.list.group.overexpression$spots <- lapply( env$spot.list.group.overexpression$spots, spot.fisher.p)
  }  
  # dmap spot
  env$spot.list.dmap$spots <- lapply( env$spot.list.dmap$spots, spot.fisher.p)

  return(env)
}
