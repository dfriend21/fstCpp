
gen_stats <- function(gen, pop = NULL){
  if(inherits(gen, "genind")){
    if(any(adegenet::ploidy(gen) != 2)){
      stop("This function currently only works with diploid data")
    }
    gen <- adegenet::genind2df(gen, oneColPerAll = TRUE)
    pop <- gen$pop
    gen <- as.data.frame(lapply(gen[,-1], function(x) suppressWarnings(as.integer(x))))
  }
  if(class(pop) != "factor") pop <- factor(pop)
  mat <- as.matrix(gen)
  
  stats <- wcCpp(mat, as.integer(pop))
  
  rownames(stats$pw_fst) <- levels(pop)[match(rownames(stats$pw_fst), 1:nlevels(pop))]
  colnames(stats$pw_fst) <- levels(pop)[match(colnames(stats$pw_fst), 1:nlevels(pop))]
  diag(stats$pw_fst) <- NA
  rownames(stats$pop_stats) <- levels(pop)[match(rownames(stats$pop_stats), 1:nlevels(pop))]
  
  return(stats)
}
