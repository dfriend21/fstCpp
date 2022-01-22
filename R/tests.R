
check_dif <- function(x, y, tol = 1.490116e-08, verbose = TRUE){
  dif <- x - y
  bool <- abs(dif) < tol
  all <- all(bool, na.rm = TRUE) & all(is.na(x) == is.na(y))
  if(verbose){
    cat("\n--------- VALUES ---------\n")
    cat("x:\n")
    print(x)
    cat("y:\n")
    print(y)
    cat("\n---- DIFS BY ELEMENT -----\n")
    print(dif)
    cat("\n---- EQUAL BY ELEMENT ----\n")
    print(bool)
    cat("\n------- ALL EQUAL: -------\n")
    cat(all, "\n")
  }
  return(all)
}

# compares against adegenet
test_Ho_Hs <- function(gen, ...){
  ag_sum <- adegenet::summary(gen)
  ag_ho_hs <- c(Ho = mean(ag_sum$Hobs), Hs = mean(ag_sum$Hexp))
  cpp_ho_hs <- gen_stats(gen)$stats[1:2]
  return(check_dif(ag_ho_hs, cpp_ho_hs, ...))
}

# compares against hierfstat
# it only gives Fst and Fis so those are the only two I can check
test_fstats <- function(gen, ...){
  hf_wc_l <- hierfstat::wc(gen)
  hf_wc <- c(Fst = hf_wc_l$FST, Fis = hf_wc_l$FIS)
  cpp_wc <- gen_stats(gen)$stats[4:5]
  return(check_dif(hf_wc, cpp_wc, ...))
}

# compares against adegenet
test_n_al <- function(gen, ...){
  ag_sum <- adegenet::summary(gen)
  ag_n_al <- ag_sum$loc.n.all
  cpp_n_al <- gen_stats(gen)$loci_stats[,"N_al"]
  return(check_dif(ag_n_al, cpp_n_al, ...))
}

# compares against adegenet
test_ho_hs_loci <- function(gen, ...){
  ag_sum <- adegenet::summary(gen)
  ag_h_loci <- cbind(ag_sum$Hobs, ag_sum$Hexp)
  cpp_h_loci <- gen_stats(gen)$loci_stats[,c("Ho", "Hs")]
  return(check_dif(ag_h_loci, cpp_h_loci, ...))
}

# compares against pegas
test_fstats_loci <- function(gen, ...){
  pg_fs <- pegas::Fst(genind2loci(gen))
  cpp_fs <- gen_stats(gen)$loci_stats[,c("Fit", "Fst", "Fis")]
  return(check_dif(pg_fs, cpp_fs, ...))
}

# compares against hierfstat
test_pw_fst <- function(gen, ...){
  hf_pw <- hierfstat::pairwise.WCfst(genind2hierfstat(gen))
  cpp_pw <- gen_stats(gen)$pw_fst
  ord_r <- match(rownames(cpp_pw), rownames(cpp_pw))
  ord_c <- match(colnames(cpp_pw), colnames(cpp_pw))
  diag(cpp_pw) <- NA
  return(check_dif(hf_pw[ord_r, ord_c], cpp_pw, ...))
}

# compares against adegenet
test_ho_hs_pop <- function(gen, ...){
  n.pop <- seppop(gen) #https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2015-November/001379.html
  ag_ho <- do.call("c", lapply(n.pop, function(x) mean(adegenet::summary(x)$Hobs, na.rm = TRUE))) 
  ag_hs <- adegenet::Hs(gen)
  ag_h <- cbind(ag_ho, ag_hs)
  cpp_h <- gen_stats(gen)$pop_stats
  return(check_dif(ag_h, cpp_h, ...))
}

test_benchmark <- function(gen){
  loc <- genind2loci(gen)
  hf <- genind2hierfstat(gen)
  bm <- bench::mark(gen_stats(gen), 
                    pegas::Fst(loc),
                    hierfstat::wc(hf),
                    hierfstat::pairwise.WCfst(hf),
                    check = FALSE)  
  return(bm)
}
