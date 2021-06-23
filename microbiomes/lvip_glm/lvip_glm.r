## get user-edited environmental variables input_prefix, output_prefix, and input file names
source(file.path(getwd(), 'lvip_glm.env'))
##

## set up output directory
dir.create(output_prefix)
outfile <- file(file.path(output_prefix, 'runlog.log'), open = 'wt')
sink(outfile, type = 'output', split = TRUE)
##

phy <- ''
counts <- ''
sample_dat <- ''


microbeTol <- exp(mean(log(microbeTree$edge.length)))
multiTree <- force.ultrametric(di2multi(microbeTree, tol=microbeTol))
finalMicrobeTree <- reorder(drop.tip(multiTree, multiTree$tip.label[!multiTree$tip.label %in% colnames(present)]), 'pruningwise')
finalMicrobeTree$edge.length <- finalMicrobeTree$edge.length / exp(mean(log(finalMicrobeTree$edge.length)))

## generate some summary numbers regarding microbes
microbeTips <- colnames(present)
NT <- length(microbeTips)
NI <- finalMicrobeTree$Nnode
NN <- NI + NT
sa <- rbind(finalMicrobeTree$edge,c(0,NT+1))[NN:1,]
time <- c(finalMicrobeTree$edge.length[1:length(finalMicrobeTree$tip.label)],1,finalMicrobeTree$edge.length[(length(finalMicrobeTree$tip.label)+1):length(finalMicrobeTree$edge.length)])

##
NSB <- length(sampleFactors)
NB_s <- ncol(modelMat)
##

standat <- list(NS           = NS,
               NI         = NI,
               NT         = NT,
               NB_s         = NB_s,
               NSB        = NSB,
               idx        = idx,
               count      = counts,
               modelMat     = modelMat,
               time       = time,
               self       = sa[,2],
               ancestor   = sa[,1],
               X_s          = modelMat)


save.image(file.path(output_prefix, 'setup.RData'))

cmdstanr::write_stan_json(standat, file.path(output_prefix, 'data.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], 'STANCFLAGS="--include-paths=', include_path, '" ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_commands[[engine]])
print(date())
system(sampling_commands[[engine]])

stan.fit <- cmdstanr::read_cmdstan_csv(file.path(output_prefix, paste0('samples_',engine,'.csv')),
                                       format = 'draws_array')

save.image(file.path(output_prefix, 'results_.RData'))
