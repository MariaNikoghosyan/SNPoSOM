pipeline.prediction.data.preprocessing <- function(env)
{
  #### check input prediction.indata ####
  if (is.null(env$prediction$prediction.indata))
  {
    util.fatal("No prediction.indata supplied!")
    env$prediction$passedInputChecking <- FALSE
    return(env)
  }
  
  if (!is(env$prediction$prediction.indata,"matrix") && (is.null(dim(env$prediction$prediction.indata)) || dim(env$prediction$prediction.indata) < 1))
  {
    util.fatal("Invalid prediction.indata! Provide a two-dimensional numerical matrix.")
    env$passedInputChecking <- FALSE
    return(env)
  }
  
  if (!is(env$prediction$prediction.indata,"matrix") ||
      mode(env$prediction$prediction.indata) != "numeric" )
    # storage.mode(indata) != "numeric")
  {
    rn <- rownames(env$prediction$prediction.indata)
    num.mode <- sapply(seq(ncol(env$prediction$prediction.indata)), function(i){ all(grepl("^-?[0-9\\.]+$",env$prediction$prediction.indata[,i])) })
    
    if( num.mode[1]==FALSE && all(num.mode[-1]==TRUE) ) # check if IDs are contained as first row
    {
      rn <- env$prediction$prediction.indata[,1]
      env$prediction$prediction.indata <- env$prediction$prediction.indata[,-1]
      num.mode <- num.mode[-1]
      util.warn("Gene IDs adopted from first data column.")
    }
    
    if( any(num.mode!=TRUE) ) # check if all columns contain numbers or convertible characters
    {
      util.fatal("Invalid train.indata! Provide a two-dimensional numerical matrix.")
      env$passedInputChecking <- FALSE
      return(env)
      
    } else
    {
      env$prediction$prediction.indata <- apply(env$prediction$prediction.indata, 2, function(x){ as.numeric(as.vector(x)) })
      rownames(env$prediction$prediction.indata) <- rn
      storage.mode( env$prediction$prediction.indata) <- "numeric"
      util.warn("Indata converted to two-dimensional numerical matrix.")    
    }
  }
  
  
  # check for the constant columns and rows
  const.cols <- which(apply(env$prediction$prediction.indata, 2, function(col) { diff(range(col)) == 0 }))
  
  if (length(const.cols) > 0)
  {
    env$prediction$prediction.indata <- env$prediction$prediction.indata[,-const.cols]
    env$group.labels <- env$group.labels[-const.cols]
    env$group.colors <- env$group.colors[-const.cols]
    util.warn("Removed",length(const.cols),"constant columns from data set.")
  }
  
  const.rows <- which(apply(env$prediction$prediction.indata, 1, function(row) { diff(range(row)) == 0 }))
  if (length(const.rows) > 0)
  {
    env$prediction$prediction.indata <- env$prediction$prediction.indata[-const.rows,]
    util.warn("Removed",length(const.rows),"constant rows from data set.")
  }
  
  # check for NAs and infinite values
  na.rows <- which( apply(env$prediction$prediction.indata, 1, function(x) sum( is.na(x) | is.infinite(x) ) ) > 0 )
  
  if (length(na.rows) > 0)
  {
    env$prediction$prediction.indata <- env$prediction$prediction.indata[-na.rows,]
    util.warn("Removed NAs or infinite values from data set")
  }
  
  # compare rownames of indata and prediction indata
  common_ids <- intersect(rownames(env$prediction$train.indata), rownames(env$prediction$prediction.indata))

  if(length(common_ids) < nrow(env$prediction$prediction.indata) |
     length(common_ids) < nrow(env$prediction$train.indata))
  {
    env$prediction$train.indata <-  env$prediction$train.indata[common_ids,]
    env$prediction$prediction.indata <- env$prediction$prediction.indata[common_ids,]
    util.warn(length(common_ids),"SNPs are matched in training and prediction data")
    
  }
  
  # check group.labels and group.colors
  if ((!is.null(env$prediction$prediction.group.labels) && length(env$prediction$prediction.group.labels) != ncol(env$prediction$prediction.indata)) ||
      (!is.null(env$prediction$prediction.group.labels) && length(env$prediction$prediction.group.colors) != ncol(env$prediction$prediction.indata)))
  {
    env$prediction$prediction.group.labels <- NULL
    env$prediction$prediction.group.labels <- NULL
    util.warn("Group assignment doesnt fit number of samples")
  }
  
  if (!is.null(env$prediction$prediction.group.labels) && max(table(env$prediction$prediction.group.colors)) == 1)
  {
    env$prediction$prediction.group.labels <- NULL
    env$prediction$prediction.group.colors <- NULL
    util.warn("Each sample has an own group")
  }
  
  if (!is.null(env$prediction$prediction.group.labels))
  {
    for (sample in unique(colnames(env$prediction$prediction.indata)))
    {
      if (length(unique(env$prediction$prediction.group.labels[which(colnames(env$prediction$prediction.indata) == sample)])) > 1)
      {
        util.warn("Sample is in multiple groups:", sample)
        env$prediction$prediction.group.labels <- NULL
        env$prediction$prediction.group.colors <- NULL
        break
      }
    }
  }
  
  if (!is.null(env$prediction$prediction.group.labels))
  {
    env$prediction$prediction.group.labels <- as.character(env$prediction$prediction.group.labels)
    names(env$prediction$prediction.group.labels) <- colnames(env$prediction$prediction.indata)
    
    if (is.null(env$prediction$prediction.group.colors))
    {
      env$prediction$prediction.group.colors <- rep("", ncol(env$prediction$prediction.indata))
      
      for (i in seq_along(unique(env$prediction$prediction.group.labels)))
      {
        env$prediction$prediction.group.colors[which(env$prediction$prediction.group.labels == unique(env$prediction$prediction.group.labels)[i])] <- color.palette.discrete(length(unique(env$prediction$prediction.group.labels)))[i]
      }
    }
    
    # catch userdefined group.colors --> convert to #hex
    if (length(unique(substr(env$prediction$prediction.group.colors, 1, 1)) > 1) || unique(substr(env$prediction$prediction.group.colors, 1, 1))[1] != "#")
    {
      env$prediction$prediction.group.colors <- apply(col2rgb(env$prediction$prediction.group.colors), 2, function(x) { rgb(x[1]/255, x[2]/255, x[3]/255) })
    }
    names(env$prediction$prediction.group.colors) <- colnames(env$prediction$prediction.indata)
  } else
  {
    env$group.labels <- rep("auto",ncol(env$prediction$prediction.indata))
    names(env$group.labels) <- colnames(env$prediction$prediction.indata)
    
    env$prediction$prediction.group.colors <- rep("#000000", ncol(env$prediction$prediction.indata))
    names(env$prediction$prediction.group.colors) <- colnames(env$prediction$prediction.indata)
  }
  
  return(env)
    
  
}

pipeline.prepare.Prediction.Indata <- function(env)
{
  # save original indata in temarory file
  env$prediction$prediction.indata.temp <- env$prediction$prediction.indata
  #env$indata.sample.mean <- colMeans(env$indata)
  
  if (env$prediction$preferences$sample.quantile.normalization)
  {
    env$prediction$prediction.indata <- Quantile.Normalization(env$prediction$prediction.indata)
  }
  
  colnames(env$prediction$prediction.indata) <- make.unique(colnames(env$prediction$prediction.indata))
  names(env$prediction$prediction.group.labels) <- make.unique(names(env$prediction$prediction.group.labels))
  names(env$prediction$prediction.group.colors) <- make.unique(names(env$prediction$prediction.group.colors))
  
  
  env$prediction$prediction.indata.gene.mean <- rowMeans(env$prediction$prediction.indata)
  
  if (env$prediction$preferences$feature.centralization)
  {
    env$prediction$prediction.indata<- env$prediction$prediction.indata - env$prediction$prediction.indata.gene.mean
  }
  
  #batch effect correction
  if(env$prediction$preferences$batch.effect.correction)
  {
    batches <- as.factor(c(rep(0, ncol(env$prediction$train.indata)),
                           rep(1, ncol(env$prediction$prediction.indata))))
    batch_data <- sva::ComBat(cbind(env$prediction$train.indata, env$prediction$prediction.indata), 
                              batches)
    
    env$prediction$train.indata <- batch_data[,batches == 0]
    env$prediction$prediction.indata <- batch_data[,batches == 1]
  }
  
  
  return(env)
}
