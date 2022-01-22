library(adegenet)
data(nancycats)

# data("dapcIllus")
data("eHGDP")
data("hybridtoy")
data("microbov")
data("nancycats")
data("rupica")
data("sim2pop")
data("spcaIllus")
data("swallowtails")

gen_stats(nancycats)
gen <- eHGDP # takes a long time
gen <- hybridtoy
gen <- microbov
gen <- nancycats
gen <- rupica #!!!!!!!!!!!!!!!!!!!
gen <- sim2pop
gen <- swallowtails

test_Ho_Hs(gen)
test_fstats(gen)
test_n_al(gen)
test_ho_hs_loci(gen)
test_fstats_loci(gen)
test_pw_fst(gen)
test_ho_hs_pop(gen)

test_benchmark(gen)
