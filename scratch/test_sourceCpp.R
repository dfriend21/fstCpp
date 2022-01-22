library(adegenet)
library(Rcpp)


sourceCpp("~/Documents/r_packages/fstCpp/src/fstats_cpp.cpp")

#---------

data(nancycats)
gen_nc <- nancycats

df_nc <- adegenet::genind2df(gen_nc, oneColPerAll = TRUE)
als_nc <- as.data.frame(lapply(df_nc[,-1], function(x) suppressWarnings(as.integer(x))))

mat_nc <- as.matrix(als_nc)
pop_nc <- as.integer(df_nc$pop)
pop_nc
#stats_nc <- wcCpp(mat_nc, pop_nc)
#stats_nc

bool_nc <- apply(mat_nc, 1, function(x) any(is.na(x)))
mat_nc2 <- mat_nc[!bool_nc,]
pop_nc2 <- pop_nc[!bool_nc]
dim(mat_nc2)
length(pop_nc2)
unique(pop_nc2)

stats_nc <- wcCpp(mat_nc2, pop_nc2)
#---------

load("/Users/dfriend/Documents/r_packages/fstCpp/scratch/tort_gen.RData")
gen_tr <- tort_gen

df_tr <- adegenet::genind2df(gen_tr, oneColPerAll = TRUE)
als_tr <- as.data.frame(lapply(df_tr[,-1], function(x) suppressWarnings(as.integer(x))))

mat_tr <- as.matrix(als_tr)
pop_tr <- as.integer(df_tr$pop)


unique(pop_tr)
stats_tr <- wcCpp(mat_tr, pop_tr)
stats_tr

#---------
head(mat_nc)
head(mat_tr)

apply(mat_nc, 2, unique)
apply(mat_tr, 2, unique)

nc_bool <- apply(mat_nc, 1, function(x) any(is.na(x)))
tr_bool <- apply(mat_tr, 1, function(x) any(is.na(x)))

sum(nc_bool)
sum(tr_bool)

mat_nc[nc_bool,]
mat_tr[tr_bool,]

class(mat_nc)
class(mat_tr)

dim(mat_nc)
dim(mat_tr)

pop_nc
pop_tr

length(pop_nc)
length(pop_tr)

class(pop_nc)
class(pop_tr)
