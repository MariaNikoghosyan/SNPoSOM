# genotype encodding functions
check.genotype <- function(genotypes)
{

  nuc <- c("AA", "AT", "AC", "AG",
           "TA", "TT", "TC", "TG",
           "CA", "CT", "CC", "CG",
           "GA", "GT", "GC", "GG")
  nuc <- c(nuc, tolower(nuc))
  genotypes <- sapply(genotypes, function(x)
    {
    if(any(x == nuc))
    {
      return(TRUE)
    }else
    {
      return(FALSE)
    }
    
  })
  return(genotypes)
}

check.nucleotides <- function(genotypes)
{
  nuc <- c("A", "C", "T", "G")
  nuc <- c(nuc, tolower(nuc))
  genotypes <- sapply(genotypes, function(x)
  {
    if(any(x == nuc))
    {
      return(TRUE)
    }else
    {
      return(FALSE)
    }
    
  })
  return(genotypes)
}

pipeline.change.complementary.nucleotide <- function(allele)
{
  if(allele == 'A')
  {
    allele <- 'T'
    return(allele)
    break()
  }
  if(allele == 'T')
  {
    allele <- 'A'
    return(allele)
    break()
  }
  if(allele == 'G')
  {
    allele <- 'C'
    return(allele)
    break()
  }
  if(allele == 'C')
  {
    allele <- 'G'
    return(allele)
    break()
  }
}

pipeline.check.complementarity <- function(allele_1, allele_2)
{
  complementarity_list <- list(AT = c("A", "T"),
                               GC = c("G", "C"))
  
  if(any(allele_1 == complementarity_list$AT) & any(allele_2 == complementarity_list$AT)|
     any(allele_1 == complementarity_list$GC) & any(allele_2 == complementarity_list$GC))
  {
    return(T)
    break()
  }else
  {
    return(F)
  }
}

