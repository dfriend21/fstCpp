#CREATE THE GENOTYPE FILE - combine columns of loci - ONLY DO THIS ONCE!
formatgenotype <- function(x){
  cols <- names(x)[grepl("allele", names(x))]
  
  alleles <- x[,cols]
  new_cols <- lapply(seq(2, ncol(alleles), by = 2), function(i){
    return(paste0(alleles[, i-1], alleles[, i]))
  })
  allele_df <- data.frame(do.call(cbind, new_cols))
  names(allele_df) <- paste0("L", 1:ncol(allele_df))
  final_df <- cbind(x[,c("tid", "zone")], allele_df)
  return(final_df)
}
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)

Rcpp::sourceCpp("~/Documents/clark_county_project/fstats_cpp.cpp")
Rcpp::sourceCpp("~/Documents/clark_county_project/fstats_cpp2.cpp")
df <- read.csv("/Users/dfriend/Documents/clark_county_project/data/sim_output/Oct_22/sampled/old/100_o.csv")

table(df$zone)
#df <- df[df$zone %in% c(1,2),]

pp <- df$zone
alls <- df[,grepl("allele", names(df))]
alls <- as.matrix(alls)

wc_cpp <- wcCpp(alls, pp)
wc_cpp2 <- wcCpp2(alls, pp)

wc_cpp
wc_cpp2$pop_Hs
wc_cpp2$loci_stats[,"Hs"]


###########################################
# check if results match other packages
###########################################

gt <- formatgenotype(df)
gen <- df2genind(gt[,c(-1,-2)], pop = df$zone, ploidy=2, ncode=3, sep=NULL, NA.char = "", ind.names = df$tid)
loc <- as.loci(gen)
hf <- genind2hierfstat(gen)


#---- overall stats

# Ho and Hs
hf_ov <- basic.stats(gen)$overall
ag_sm <- summary(gen)
ag_ov <- c(Ho = mean(ag_sm$Hobs), Hs = mean(ag_sm$Hexp))
hf_ov[1:2]
ag_ov
wc_cpp2$stats[1:2]

# F-stats
hf_wc <- wc(gen)
hf_wc
wc_cpp2$stats[3:5]

#---- loci stats

# number of alleles
ag_sm <- summary(gen)
ag_sm$loc.n.all
wc_cpp2$loci_stat[,"N_al"]

# Ho and Hs
hf_bs <- basic.stats(gen)$perloc
ag_sm <- summary(gen)
pp_h <- locus_table(gen)

hf_bs[,c("Ho", "Hs")]
cbind(ag_sm$Hobs, ag_sm$Hexp)
wc_cpp2$loci_stats[,c("Ho", "Hs")]

cbind(wc_cpp2$loci_stats[,"Hs"], pp_h[,"Hexp"])

# F-stats
pg_fst <- Fst(loc)
pg_fst
wc_cpp2$loci_stats[,c("Fit", "Fst", "Fis")]

#---- pw Fst

rn <- rownames(wc_cpp2$pop_stats)
m <- function(vals) match(names(vals), rn)

hf_pw <- pairwise.WCfst(hf)
hf_pw[m(hf_pw[,1]), m(hf_pw[1,])]
wc_cpp2$pw_fst

#---- population stats

rn <- rownames(wc_cpp2$pop_stats)
m <- function(vals) match(names(vals), rn)

# Ho
hf_ho <- hierfstat::Ho(gen)
n.pop <- seppop(gen) #https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2015-November/001379.html
ag_ho <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 

hf_ho[m(hf_ho)]
ag_ho[m(ag_ho)]
wc_cpp2$pop_stats[, "Ho"]

# Hs/He
hf_hs <- hierfstat::Hs(gen)
ag_hs <- adegenet::Hs(gen)

hf_hs[m(hf_hs)]
ag_hs[m(hf_hs)]
wc_cpp2$pop_stats[, "Hs"]




###########################################
# benchmark
###########################################

bench::mark(wcCpp(alls,pp), wcCpp2(alls,pp), check = FALSE)
bm <- bench::mark(wcCpp(alls, pp), pegas::Fst(loc) , hierfstat::wc(hf), hierfstat::pairwise.WCfst(hf), check = FALSE)
#bench::mark(wcCpp(alls, pp), pairwise.WCfst(hf), check = FALSE)
