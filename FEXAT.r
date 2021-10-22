# fexat 3.2
# Peter Kraft (pkraft@ucla.edu)
# 29 October 2002
#
# fexat (for Family Expression Association Test) takes as input:
#
# - x: a (number of genes) x (number of conditions) matrix of gene 
#   expression measurements
# - t: a (number of traits) x (number of conditions) matrix of trait 
#   measurements
# - fam: a (number of conditions) x 1 matrix of family (stratum) 
#   indentifiers
# 
# and returns a list with elements:
#
# - z.stat: a (number of genes) x (number of traits) matrix of 
#   statistics testing correlation between genes and traits
# - p.value: a (number of genes) x (number of traits) matrix of 
#   p-values associated with z.stat
# - est: a (number of genes) x (number of traits) matrix of 
#   stratum-adjusted correlation estimates
#
# The columns of x and t and rows of fam should correspond to the 
# same conditions. E.g., if the first column in x contains the gene
# expression measurements for Heinzl, then the first column of t
# should contain Heinzl's trait values and the first element of 
# fam should contain Heinzl's familiy identifier.
#
# - 29 Oct 2002 Rewritten to autmatically sort input matrices
# - 15 Oct 2002 Rewritten to accomodate missing phe,exp
# - 15 Oct 2002 Added correlation estimate


fexat <- function(x,t,fam,method=1)
{  
  fam <- as.matrix(fam)
  if (dim(fam)[2] != 1)  fam <- t(fam) #in case user did not read manual
  sort.fam <- order(t(fam))
  fam <- fam[sort.fam]
  x <- x[,sort.fam]
  t <- t[,sort.fam]
  
  fexat.1d1 <- function(exp,phe,fam)
  { 
    e.x <- rep(tapply(exp,fam,mean),tapply(exp,fam,length))
    e.t <- rep(tapply(phe,fam,mean),tapply(phe,fam,length))
    r.x <- exp-e.x
    r.t <- phe-e.t
    s <- tapply(t(phe*r.x),fam,sum)
    v <- sum(tapply(r.x^2,fam,sum)*tapply(r.t^2,fam,sum)/(tapply(phe,fam,length)-1)) #permutation test a la Cox and Hinkley (1974) p. 185
    #see also expression (13) in Mantel (1963) JASA 58:690-700
    d <- sum(sqrt(tapply(r.x^2,fam,sum)*tapply(r.t^2,fam,sum))) #stratified version of denominator for correlation estimate
    c(sum(s)/sqrt(v),sum(s)/d) 
  }
  
  fexat.1d2 <- function(exp,phe,fam)
  { 
    e.x <- rep(tapply(exp,fam,mean),tapply(exp,fam,length))
    r.x <- exp-e.x
    s <- tapply(phe*r.x,fam,sum)
    d <- sum(sqrt(tapply(r.x^2,fam,sum)*tapply(r.t^2,fam,sum))) #stratified version of denominator for correlation estimate
    c(sum(s)/sqrt(sum(s^2)),sum(s)/d) #original proposed fexat--uses empircal variance estimate, not recommended when number of families is small
  }
  
  z <- rep(666,NROW(x)*NROW(t))
  dim(z) <- c(NROW(x),NROW(t))
  e <- rep(666,NROW(x)*NROW(t))
  dim(e) <- c(NROW(x),NROW(t))
  p <- rep(666,NROW(x)*NROW(t))
  dim(p) <- c(NROW(x),NROW(t))
  
  if (method == 1)
  {
    for (i in c(1:NROW(t)))
    { 
      t1 <- t(t[i,])
      for (j in c(1:NROW(x)))
      {
        x1 <- t(x[j,])
        statest <- fexat.1d1(t1[complete.cases(t1,x1)],x1[complete.cases(t1,x1)],fam[complete.cases(t1,x1)])
        z[j,i] <- statest[1]
        p[j,i] <- 2*(1-pnorm(abs(z[j,i])))
        e[j,i] <- statest[2]
      }
    }
  }
  
  if (method == 2)
  {
    for (i in c(1:NROW(t)))
    { 
      t1 <- t[i,]
      for (j in c(1:NROW(x)))
      {
        x1 <- x[j,]
        z[j,i] <- fexat.1d2(t(x1),t(t1),fam)
        p[j,i] <- 2*(1-pnorm(abs(z[j,i])))
        e[j,i] <- statest[2]
      }
    }
  }
  
  out <- list(z.stat=z,p.value=p,est=e)
  out
}