# get the SNPs from neighvor metagenes
get.neighbors.radius <- function(x,dim=env$preferences$dim.1stLvlSom, r=env$prediction$metagene.radius)
{
  y <- (x - 1)%/%dim + 1
  x <- (x - 1)%%dim + 1
  
  ret <- c()
  
  for(x.i in 1:dim)
    for(y.i in 1:dim)
    {
      if( sqrt( (x.i-x)^2 + (y.i-y)^2  ) < r ) ret[length(ret)+1] = (y.i - 1) * dim + x.i
    }
  
  return(ret)
}

pipeline.svmr.training <- function(env)
{
  util.info("Start SVMR training with radius = ", env$prediction$metagene.radius)
  
  env$prediction$predicted.metadata <- matrix(0,nrow = nrow(env$metadata), ncol = ncol(env$prediction$prediction.indata),
                                              dimnames = list(rownames(env$metadata),
                                                              colnames(env$prediction$prediction.indata)))
  # run the svmr training
  pb <- txtProgressBar(0,nrow(env$prediction$predicted.metadata),style=3)
  for( i in 1:nrow(env$prediction$predicted.metadata) )
  {
    env$prediction$feature.BMU <- env$som.result$feature.BMU[rownames(env$prediction$prediction.indata)]
    
    metagene.determinants = get.neighbors.radius(i)
    metagene.determinants = do.call(c,lapply( metagene.determinants, 
                                              function(x) names(which( env$prediction$feature.BMU == x )) ) )
    if(length(metagene.determinants)>0)
    {
      f = as.formula( paste( "M ~", paste(metagene.determinants,collapse="+") ) )
      df = data.frame( t(env$prediction$train.indata[metagene.determinants,,drop=F]), 
                       t(env$metadata [i,,drop=F]) )
      colnames(df)[ncol(df)] = "M"
      
      env$prediction$svmr.model.list[[i]] = e1071:::svm(formula = f, data = df, type = "eps-regression")
      env$prediction$metagene.determinant.list[[i]] = metagene.determinants
      
    }
    setTxtProgressBar(pb, i)    
  }
  close(pb)
  return(env)
}

pipeline.svmr.prediction <- function(env)
{
  util.info("Start SVMR prediction")

  # run the svm prediction

  pb <- txtProgressBar(0,nrow(env$prediction$predicted.metadata),style=3)
  
  for( i in 1:nrow(env$prediction$predicted.metadata) )
  {
    if(!is.null(env$prediction$svmr.model.list[[i]]))  
      env$prediction$predicted.metadata[i,] = e1071:::predict.svm(env$prediction$svmr.model.list[[i]], 
                                                                  t(env$prediction$prediction.indata[env$prediction$metagene.determinant.list[[i]],,drop=F]))
    setTxtProgressBar(pb, i)    
  }
  close(pb)
  return(env)
}
