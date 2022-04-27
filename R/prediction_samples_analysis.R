pipeline.prediction.Samples.Analysis <- function(env)
{
  if(env$prediction$preferences$annotation.analysis &  env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    env$prediction$gene.info <-lapply(env$gene.info[c("ids","names","descriptions","chr.name",
                                                      "chr.band","chr.start","coordinates")] , function(x)
                                                      {
                                                        x <- x[names(x) %in% rownames(env$prediction$prediction.indata)]
                                                        return(x)
                                                      })
    
    env$prediction$gs.def.list <- lapply(env$gs.def.list, function(x)
    {
      x$Genes <- x$Genes[x$Genes %in% env$prediction$gene.info$ids]
      x$Snps <- x$Snps[x$Snps %in% names(env$prediction$gene.info$ids)]
      return(x)
    })
  }
  
  
  #define sample spot
  local.env <- new.env()
  local.env$indata <- env$prediction$prediction.indata
  local.env$gene.info <- env$prediction$gene.info
  local.env$gs.def.list <- env$prediction$gs.def.list
  local.env$metadata <- env$prediction$predicted.metadata
  local.env$preferences$dim.1stLvlSom <-  env$preferences$dim.1stLvlSom
  local.env$preferences$dim.2ndLvlSom <-  env$preferences$dim.2ndLvlSom
  local.env$som.result$feature.BMU <- env$prediction$feature.BMU
  local.env$indata.temp <- env$prediction$prediction.indata.temp
  local.env$preferences$group.maf <- ifelse(env$preferences$snp.analysis, yes = env$preferences$group.maf, no = F)
  local.env$color.palette.portraits <- env$color.palette.portraits
  local.env$color.palette.heatmaps <- env$color.palette.heatmaps
  local.env$preferences$activated.modules$geneset.analysis <- env$preferences$activated.modules$geneset.analysis
  local.env$output.paths <- paste(env$files.name, "- Results/Prediction/Summary Sheets - Samples")
  names(local.env$output.paths) <- "Summary Sheets Samples"
  
  
  local.env <- pipeline.detectSpotsSamples(local.env)
  local.env <- pipeline.diffExpressionStatistics(local.env)
  if(env$prediction$preferences$annotation.analysis & env$preferences$activated.modules$geneset.analysis)
  {
    local.env <- pipeline.genesetStatisticSamples(local.env)
  }

  
  pipeline.summarySheetsSamples(local.env)
  
  #save the results
  env$prediction$fdr.g.m <- local.env$fdr.g.m
  env$prediction$n.0.m <- local.env$n.0.m
  env$prediction$p.g.m <- local.env$p.g.m
  env$prediction$p.m <- local.env$p.m
  env$prediction$perc.DE.m <- local.env$perc.DE.m
  env$prediction$samples.GSZ.scores <- local.env$samples.GSZ.scores
  env$prediction$spot.list.samples <- local.env$spot.list.samples
  
  return(env)

}