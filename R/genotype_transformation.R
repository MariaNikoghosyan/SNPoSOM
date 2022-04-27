pipeline.genotype.transformation <- function(env)
{
  # check indata input parameters
  
  if (is.null(env$indata))
  {
    util.fatal("No indata supplied!")
    env$passedInputChecking <- FALSE
    return(env)
  }
  
  if (is(env$indata,"ExpressionSet"))
  {
    env$group.labels <- as.character(pData(env$indata)$group.labels)
    env$group.colors <- as.character(pData(env$indata)$group.colors)
    env$indata <- assayData(env$indata)$exprs
  }
  
  if (!is(env$indata,"matrix") && (is.null(dim(env$indata)) || dim(env$indata) < 1))
  {
    util.fatal("Invalid indata! Provide a two-dimensional numerical matrix.")
    env$passedInputChecking <- FALSE
    return(env)
  }
  
  #check for missing genotypes
  env$indata <- t(apply(env$indata, 1, function(x)
  {
    x[is.na(x)] <- "-"
    x[which(!check.genotype(x))] <- "-"
    return(x)
  }))
  

  #### global minor allele encoding   #### 
  if(env$preferences$genotype.transformation== 'global.minor.major.alleles')
  {
    biomart.table <- NULL
    
    try({
      mart <- useMart("ENSEMBL_MART_SNP", host=env$preferences$database.host)
      mart <- useDataset("hsapiens_snp", mart=mart)
      
      
      
      #query = c("refsnp_id","chr_name","ensembl_gene_stable_id")[ which( c("refsnp_id","chr_name","ensembl_gene_stable_id") %in% listAttributes(mart)[,1] ) ][1:2]
      suppressWarnings({  biomart.table <-
        getBM(c('refsnp_id', 'allele','minor_allele', 'chr_name'),
              "snp_filter",
              rownames(env$indata),
              mart, checkFilters=FALSE)  })
    }, silent=TRUE)
    
    if(is.null(biomart.table))
    {
      util.warn("No annotation for SNPs. Switching to minor.major.alleles genotype encoding")
      env$preferences$genotype.transformation <- "minor.major.alleles"
    }else
    {

      #define all multiallelic SNPs and remove from biomart.table
      ind <- c()
      for(i in 1:nrow(biomart.table))
      {
        if(lengths(regmatches(biomart.table$allele[i], gregexpr("/", biomart.table$allele[i]))) > 1)
        {
          ind <- c(ind,i)
        }
      }
      if(length(ind) > 1)
      {
        biomart.table <- biomart.table[-ind,]
      }
      # remove duplications
      biomart.table <- biomart.table[!duplicated(biomart.table$refsnp_id),]
      
      # separate two alleles into different columns
      biomart.table$allele_1 <- NA
      biomart.table$allele_2 <- NA
      
      for (i in 1:nrow(biomart.table)) 
      {
        k <- strsplit(biomart.table$allele[i], "/")[[1]]
        biomart.table$allele_1[i] <- k[1]
        biomart.table$allele_2[i] <- k[2]
      }
      
      biomart.table <- biomart.table[check.nucleotides(biomart.table$allele_1),]
      biomart.table <- biomart.table[check.nucleotides(biomart.table$allele_2),]
      biomart.table <- biomart.table[check.nucleotides(biomart.table$minor_allele),]

      #define major allele
      biomart.table$major_allele <- NA
      biomart.table$major_allele <- ifelse(biomart.table$minor_allele == biomart.table$allele_1,
                                           biomart.table$allele_2,  biomart.table$allele_1)

      # filter indata
      ind <- setdiff(rownames(env$indata), biomart.table$refsnp_id)
      if(length(ind) > 1)
      {
        util.warn(paste(length(ind), "SNPs do not have annotation: removed from data"))
      }
      env$indata <- env$indata[biomart.table$refsnp_id,]
      
      
      # define indata alleles

      indata_alleles <- apply(env$indata, 1,FUN =  function(x)
      {
        x <- unique(unlist(strsplit(as.character(x), split = "")))
        return(x)
      })
      
      # indata_alleles <- t(indata_alleles)
      # indata_alleles <- apply(indata_alleles, 1,unique)
      # 
      # indata_alleles <- t(indata_alleles)
      # indata_alleles <- setNames(split(indata_alleles, seq(nrow(indata_alleles))), rownames(indata_alleles))
      
      
      
      biomart.table$indata_minor_allele <- NA
      biomart.table$indata_major_allele <- NA
      

      for (i in 1:nrow(biomart.table))
      {
        #minor allele
        if(any(biomart.table$minor_allele[i] == indata_alleles[biomart.table$refsnp_id[i]][[1]]))
        {
          biomart.table$indata_minor_allele[i] <- biomart.table$minor_allele[i]
        }else
        {
          biomart.table$indata_minor_allele[i] <- pipeline.change.complementary.nucleotide(biomart.table$minor_allele[i])
        }
        #major allele
        if(any(biomart.table$major_allele[i] == indata_alleles[biomart.table$refsnp_id[i]][[1]]))
        {
          biomart.table$indata_major_allele[i] <- biomart.table$major_allele[i]
        }else
        {
          biomart.table$indata_major_allele[i] <- pipeline.change.complementary.nucleotide(biomart.table$major_allele[i])
        }
      }
      
      
      #biomart.table <- split(biomart.table[,c('minor_allele','major_allele', 'indata_minor_allele','indata_major_allele')], biomart.table$refsnp_id)
      
      for (i in 1:nrow(env$indata))
      {
        env$indata[i,] <- sapply(env$indata[i,], FUN = function(x)
        {
          x <- strsplit(as.character(x),split = "")[[1]]
          if(length(unique(x)) > 1)
          {
            x <- 1
            return(x)
            break()
          }else
          {
            if(unique(x) == biomart.table$minor_allele[which(biomart.table$refsnp_id %in% rownames(env$indata)[i])])
              
            {
              x <- 2
              return(x)
              break()
            }
            if(unique(x) == biomart.table$major_allele[which(biomart.table$refsnp_id %in% rownames(env$indata)[i])])
            {
              x <- 0
            }
            else
            {
              x <- 0
              return(x)
              break()
            }
          }

          
        })
      }
    }
    

    
  }
  
  #minor allele encoding
  if(preferences$genotype.transformation == 'minor.major.alleles') # calculate minor and major alleles
  {

    # define indata alleles
    indata_alleles <- apply(env$indata, 1,FUN =  function(x)
    {
      x <- unlist(strsplit(as.character(x), split = ""))
      return(x)
    })
    
    #remove multiallic SNPs
    multiallelic <- sapply(indata_alleles, function(x)
      {
      x <- x[check.nucleotides(x)]
      return(length(unique(x)))
    })
    
    # filter indata remove multiallic SNPs
    ind <- which(multiallelic > 2)
    if(length(ind) > 1)
    {
      util.warn(paste(length(ind), "SNPs are not diallelic: removed from data"))
      
    }
    env$indata <- env$indata[multiallelic <= 2,]

    
    #create minor.major.alleles
    env$minor.major.alleles <- data.frame(SNP_ID = rep("-",nrow(env$indata)),
                                          Minor.allele = rep("-",nrow(env$indata)),
                                          Minor.allele.frequency = rep(0,nrow(env$indata)),
                                          Major.allele = rep("-",nrow(env$indata)),
                                          Major.allele.frequency = rep(0,nrow(env$indata)))
    
    env$minor.major.alleles$SNP_ID <- rownames(env$indata)
    
    
    for (i in seq_along(indata_alleles))
    {
      nuc <- indata_alleles[[i]][check.nucleotides(indata_alleles[[i]])]

      if(length(nuc) > 0)
      {
        nuc.count <- as.vector(table(nuc))
        names(nuc.count) <- names(table(nuc))
        
        nuc.minor <- unname(ifelse(nuc.count[1] != nuc.count[2],
                                   yes = names(nuc.count)[which(nuc.count == min(nuc.count))],
                                   no = names(nuc.count)[1]))
        
        nuc.major <- unname(ifelse(nuc.count[1] != nuc.count[2],
                                   yes = names(nuc.count)[which(nuc.count != min(nuc.count))],
                                   no = names(nuc.count)[2]))
        
        #minor allele
        env$minor.major.alleles$Minor.allele[i] <- nuc.minor
        env$minor.major.alleles$Minor.allele.frequency[i] <- sum(indata_alleles[[i]] == nuc.minor)/length(indata_alleles[[i]])
        #major allele
        env$minor.major.alleles$Major.allele[i] <- nuc.major
        env$minor.major.alleles$Major.allele.frequency[i] <- sum(indata_alleles[[i]] == nuc.major)/length(indata_alleles[[i]])
      }else
      {
        
      }

    }
    
    
    for (i in 1:nrow(env$indata))
    {
      
      env$indata[i,] <- sapply(env$indata[i,], function(x)
      {
        x <- strsplit(as.character(x),split = "")[[1]]
        if(length(unique(x)) > 1)
        {
          x <- 1
          return(x)
          break()
        }
        if(unique(x) == env$minor.major.alleles$Minor.allele[i])
        {
          x <- 2
          return(x)
          break()
        }
        if(unique(x) == env$minor.major.alleles$Major.allele[i])
        {
          x <- 0
        }
        else
        {
          x <- 0
          return(x)
          break()
        }
        
      })
    }
    
  }

  #disease allele encoding
  
  if(preferences$genotype.transformation == 'disease.assocoated.alleles')
  {
    gwas <- read.table('data/gwas_catalogue.csv', sep = '\t', header = T, as.is = T)
    gwas <- gwas[,c("SNPS","RISK_ALLELE_starnd_1", 'neutral_allele_strand_1')]
    gwas <- gwas[!duplicated(gwas),]
    gwas <- gwas[which(gwas$SNPS %in% rownames(env$indata)),]
    
    ind <- setdiff(rownames(env$indata),gwas$SNPS)
    if(length(ind) > 1)
    {
      util.warn(paste(length(ind), "SNPs with no disease association are removed from data"))
      
    }

    env$indata <- env$indata[which(rownames(env$indata) %in% gwas$SNPS),]
    

    ## define indata alleles and corresponding disease associated alleles
    indata_alleles <- apply(env$indata, 1,FUN =  function(x)
    {
      x <- unlist(strsplit(as.character(x), split = ""))
      return(x)
    })
    
    gwas_alleles <- split(gwas[,c('RISK_ALLELE_starnd_1','neutral_allele_strand_1')], gwas$SNPS)
    
    
    
    env$disease.alleles <- data.frame(SNP_ID = vector("character"),
                                         SNP_uniq_ID = vector("character"),
                                         indata_allele_1= vector("character"),
                                         indata_allele_2 = vector("character"),
                                         disease_associated_allele= vector("character"), 
                                         neutral_allele = vector("character"),
                                         disease_associated_allele_in_indata= vector("character"), 
                                         neutral_allele_in_indata = vector("character"))
    
    for (i in 1:length(indata_alleles))
    {
      
      #i <- which(names(indata_alleles) == "rs10903129")
      
      snp <- as.data.frame(matrix(NA, nrow = nrow(gwas_alleles[names(indata_alleles[i])][[1]]), ncol = ncol(env$disease.alleles)))
      names(snp) <- colnames(env$disease.alleles)
      
      snp$SNP_ID <- names(indata_alleles)[i]
      if(nrow(snp) > 1)
      {
        snp$SNP_uniq_ID  <- paste(names(indata_alleles)[i], c("", paste("_", 1:(nrow(snp)-1), sep = "")), sep = "")
      }else
      {
        snp$SNP_uniq_ID  <- names(indata_alleles)[i]
      }

      
      snp$indata_allele_1 <- unique(indata_alleles[[i]])[1]
      snp$indata_allele_2 <- unique(indata_alleles[[i]])[2]
      
      snp$disease_associated_allele <-gwas_alleles[names(indata_alleles[i])][[1]][1]
      snp$neutral_allele <- gwas_alleles[names(indata_alleles[i])][[1]][2]
      
      snp$disease_associated_allele_in_indata <- sapply(unlist(snp$disease_associated_allele), function(x) 
      {
        if(pipeline.check.complementarity(x,unique(snp$indata_allele_1)))
        {
          return(unique(snp$indata_allele_1))
        }
        if(pipeline.check.complementarity(x,unique(snp$indata_allele_2)))
        {
          return(unique(snp$indata_allele_2))
        }})
      
      snp$neutral_allele_in_indata <- sapply(unlist(snp$neutral_allele), function(x) 
      {
        if(pipeline.check.complementarity(x,unique(snp$indata_allele_1)))
        {
          return(unique(snp$indata_allele_1))
        }
        if(pipeline.check.complementarity(x,unique(snp$indata_allele_2)))
        {
          return(unique(snp$indata_allele_2))
        }})
      env$disease.alleles <- rbind(env$disease.alleles, snp)
      
    }
    

    
    indata_numeric_genotypes <- as.data.frame(matrix(0, nrow = nrow(env$disease.alleles), ncol = ncol(env$indata),
                                                     dimnames = list( env$disease.alleles$SNP_uniq_ID,colnames(env$indata))))
    
    
    for (i in 1:nrow(indata_numeric_genotypes))
    {
      k <- env$disease.alleles[which(env$disease.alleles$SNP_uniq_ID %in% rownames(indata_numeric_genotypes)[i]),]
      indata_numeric_genotypes[i,] <- sapply(as.character(env$indata[which(rownames(env$indata) %in% k$SNP_ID),]), FUN = function(x)
      {
        x <- strsplit(as.character(x),split = "")[[1]]
        if(length(unique(x)) > 1)
        {
          x <- 1
          return(x)
          break()
        }
        if(unique(x) == k$disease_associated_allele_in_indata)
        {
          x <- 2
          return(x)
          break()
        }
        if(unique(x) == k$neutral_allele_in_indata)
        {
          x <- 0
        }
        else
        {
          x <- 0
          return(x)
          break()
        }
        
      })
    }
    
    
    env$indata <- indata_numeric_genotypes
    
  }
}

  
  
  

    
    
    
