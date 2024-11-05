This is a small example and all the source code related to Iuchi and Hamada, 2024.

```
library(tidyverse)
library(dyntoy)
library(edgeR)
library(spectral)

#Function difinition
LS <- function(g, time1, exp1, time2 = NULL, exp2 = NULL, min.freq = 0, max.freq = 1, length.freq = 101, m = "generalized", sim.num = 100){
  .quiet <- function(x){ 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }
  time1 <- time1[!is.infinite(time1)]
  exp1 <- exp1[,!is.infinite(time1)]
  if(is.null(time2)){
    seq1 <- exp1[which(rownames(exp1) == g),]
    ls1 <- spec.lomb(x = time1, y = seq1,
                     f = seq(min.freq, max.freq, length = length.freq), mode = m) %>% .quiet
    return(tibble(gene = g, p.dynamic = min(ls1$p), p.de = NA))
  }else{
    time2 <- time2[!is.infinite(time2)]
    exp2 <- exp2[,!is.infinite(time2)]
    seq1 <- exp1[which(rownames(exp1) == g),]
    seq2 <- exp2[which(rownames(exp2) == g),]
    ls1 <- spec.lomb(x = time1, y = seq1,
                     f = seq(min.freq, max.freq, length = length.freq), mode = m) %>% .quiet
    ls2 <- spec.lomb(x = time2, y = seq2,
                     f = seq(min.freq, max.freq, length = length.freq), mode = m) %>% .quiet
    
    null.hypothesis <- NULL
    for(i in seq(sim.num)){
      sim.ls1 <- spec.lomb(x = sample(time1, length(time1)),
                           y = seq1,
                           f = seq(min.freq, max.freq, length = length.freq), mode = m) %>% .quiet
      sim.ls2 <- spec.lomb(x = sample(time2, length(time2)),
                           y = seq2,
                           f = seq(min.freq, max.freq, length = length.freq), mode = m) %>% .quiet
      null.hypothesis <- c(null.hypothesis, dist(rbind(sim.ls1$A, sim.ls2$A), method = m))
    }
    p <- sum(dist(rbind(ls1$A, ls2$A)) <  null.hypothesis) / sim.num
    return(tibble(gene = g, p.dynamic = min(c(ls1$p, ls2$p)), p.de = p))
  }
}

#Generation of simulation data
set.seed(1)
num.cells <- 1000
dataset <- generate_dataset(model = model,
                                        num_cells = num.cells,
                                        num_features = num.features,
                                        differentially_expressed_rate = de.rate) %>%
      add_root(root_milestone_id = .$prior_information$start_milestones) %>%
      add_pseudotime()
expression <- cpm(t(as.matrix(dataset$counts)))
pseudotime <- dataset$pseudotime / max(dataset$pseudotime)
pt.df <- tibble(pt.method = "Ground_truth",
                cell.num = num.cells,
                cell = colnames(counts),
                pseudotime = dataset$pseudotime)
ls.res <- map_dfr(rownames(expression), LS, pseudotime, expression)
```
