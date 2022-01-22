library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
load("/Users/dfriend/Documents/r_packages/genCpp/scratch/tort_gen.RData")
gen <- tort_gen

data(nancycats)
gen <- nancycats

wc_cpp <- gen_stats(gen)
wc_cpp
gen[gen$pop == "P17"]
summary(gen[gen$pop == "P17"])
genind2df(gen[gen$pop == "P17"])
wc_cpp
###########################################
# check if results match other packages
###########################################
zero_est <- 1e-10

gen <- nancycats
loc <- as.loci(gen)
hf <- genind2hierfstat(gen)


#---- overall stats

# Ho and Hs
hf_ov <- basic.stats(gen)$overall
ag_sm <- summary(gen)
ag_ov <- c(Ho = mean(ag_sm$Hobs), Hs = mean(ag_sm$Hexp))
hf_ov[1:2]
ag_ov
wc_cpp$stats[1:2]
dif_ov <- ag_ov - wc_cpp$stats[1:2]
abs(dif_ov) < zero_est

# F-stats
hf_wc <- wc(gen)
hf_wc
wc_cpp$stats[4:5]
dif_wc <- c(hf_wc$FST, hf_wc$FIS) - wc_cpp$stats[4:5]
dif_wc
abs(wc_dif) < zero_est
#---- loci stats

# number of alleles
ag_sm <- summary(gen)
ag_sm$loc.n.all
wc_cpp$loci_stat[,"N_al"]
dif_n_al <- ag_sm$loc.n.all - wc_cpp$loci_stat[,"N_al"]
dif_n_al
abs(dif_n_al) < zero_est

# Ho and Hs
hf_bs <- basic.stats(gen)$perloc
ag_sm <- summary(gen)
pp_h <- locus_table(gen)

hf_bs[,c("Ho", "Hs")]
ag_h <- cbind(ag_sm$Hobs, ag_sm$Hexp)
ag_h
wc_cpp$loci_stats[,c("Ho", "Hs")]
dif_h <- ag_h - wc_cpp$loci_stats[,c("Ho", "Hs")]
abs(dif_h) < zero_est

cbind(wc_cpp$loci_stats[,"Hs"], pp_h[,"Hexp"])

# F-stats by locus
pg_fst <- Fst(loc)
pg_fst
wc_cpp$loci_stats[,c("Fit", "Fst", "Fis")]
dif_fst <- pg_fst - wc_cpp$loci_stats[,c("Fit", "Fst", "Fis")]
dif_fst
abs(dif_fst) < zero_est

#---- pw Fst
rn <- rownames(wc_cpp$pop_stats)
m <- function(vals) match(names(vals), rn)

hf_pw <- pairwise.WCfst(hf)
hf_pw[m(hf_pw[,1]), m(hf_pw[1,])]
wc_cpp$pw_fst
dif_pw <- hf_pw[m(hf_pw[,1]), m(hf_pw[1,])] - wc_cpp$pw_fst
dif_pw
abs(dif_pw) < zero_est

#---- population stats

rn <- rownames(wc_cpp$pop_stats)
m <- function(vals) match(names(vals), rn)

# Ho
hf_ho <- hierfstat::Ho(gen)
n.pop <- seppop(gen) #https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2015-November/001379.html
ag_ho <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 

hf_ho[m(hf_ho)]
ag_ho[m(ag_ho)]
wc_cpp$pop_stats[, "Ho"]
dif_ho <- ag_ho[m(ag_ho)] - wc_cpp$pop_stats[, "Ho"]
abs(dif_ho) < zero_est

# Hs/He
hf_hs <- hierfstat::Hs(gen)
ag_hs <- adegenet::Hs(gen)

hf_hs[m(hf_hs)]
ag_hs[m(hf_hs)]
wc_cpp$pop_stats[, "Hs"]
dif_hs <- ag_hs[m(hf_hs)] - wc_cpp$pop_stats[, "Hs"]
dif_hs
abs(dif_hs) < zero_est

###########################################
# benchmark
###########################################

bench::mark(wcCpp(alls,pp), wcCpp2(alls,pp), check = FALSE)
bm <- bench::mark(wcCpp(alls, pp), pegas::Fst(loc) , hierfstat::wc(hf), hierfstat::pairwise.WCfst(hf), check = FALSE)
#bench::mark(wcCpp(alls, pp), pairwise.WCfst(hf), check = FALSE)
