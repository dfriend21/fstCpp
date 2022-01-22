library(adegenet)
library(fstCpp)
library(data.table)
data(nancycats)
gen <- nancycats

gen_stats(nancycats)

format_genotype <- function(x){
  cols <- names(x)[grepl("allele", names(x))]
  
  alleles <- data.frame(x[,cols])
  alleles[alleles == -9] <- 999
  new_cols <- lapply(seq(2, ncol(alleles), by = 2), function(i){
    return(paste0(alleles[, i-1], "", alleles[, i]))
  })
  allele_df <- data.frame(do.call(cbind, new_cols))
  names(allele_df) <- paste0("L", 1:ncol(allele_df))
  final_df <- cbind(x[,c("tid", "pop")], allele_df)
  return(final_df)
}
load("/Users/dfriend/Documents/r_packages/genCpp/scratch/tort_gen.RData")
gen_dt <- fread("/Users/dfriend/OneDrive/Documents/Dissertation/Chapters/Ch2/data/genetics/genetics.csv")
gt <- format_genotype(data.frame(gen_dt))
gen <- df2genind(gt[,c(-1,-2)], pop = gen_dt$pop, ploidy=2, ncode=3, sep=NULL, NA.char = "999", ind.names = gen_dt$tid)
gen <- tort_gen
gen_stats(gen)


if(inherits(gen, "genind")){
  if(any(adegenet::ploidy(gen) != 2)){
    stop("This function currently only works with diploid data")
  }
  gen <- adegenet::genind2df(gen, oneColPerAll = TRUE)
  pop <- gen$pop
  names(gen) <- gsub("[.]", "_", names(gen))
  gen <- as.data.frame(lapply(gen[,-1], function(x) suppressWarnings(as.integer(x))))
}
if(class(pop) != "factor") pop <- factor(pop)
mat <- as.matrix(gen)

stats <- wcCpp(mat, as.integer(pop))

rownames(stats$pw_fst) <- levels(pop)[match(rownames(stats$pw_fst), 1:nlevels(pop))]
colnames(stats$pw_fst) <- levels(pop)[match(colnames(stats$pw_fst), 1:nlevels(pop))]
rownames(stats$pop_stats) <- levels(pop)[match(rownames(stats$pop_stats), 1:nlevels(pop))]

return(stats)