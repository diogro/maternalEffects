source('./read_mouse_data.R')

library(MCMCglmm)
library(plyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(ggplot2)
library(doMC)

registerDoMC(50)

get_design = function(ind, locus, cromossome){
    gen_cols = paste0(c("A", "D", "I"), locus)
    mother = f3_pedigree[which(ind == f3_pedigree$ID), 'Dam']
    mother_gen = mouse_gen[[cromossome]][mouse_gen[[cromossome]]$ID == mother, gen_cols[-3]]
    names(mother_gen) = paste0("M", names(mother_gen))
    ind_gen = mouse_gen[[cromossome]][mouse_gen[[cromossome]]$ID == ind, gen_cols]
    genotypes = matrix(data.frame(ind_gen, mother_gen))
    if     (all(genotypes == c( 1, 0,  0,  1, 0)))  design = c( 0, -1/2,  0,  0, -1/2,  1,  1/4)
    else if(all(genotypes == c( 0, 1, -1,  1, 0)))  design = c( 0,  1/2,  0,  1, -1/2,  0, -1/4)
    else if(all(genotypes == c( 1, 0,  0,  0, 1)))  design = c( 1, -1/2,  0,  0,  1/2,  0, -1/4)
    else if(all(genotypes == c( 0, 1, -1,  0, 1)))  design = c( 0,  1/2, -1,  0,  1/2,  0,  1/4)
    else if(all(genotypes == c( 0, 1,  1,  0, 1)))  design = c( 0,  1/2,  1,  0,  1/2,  0,  1/4)
    else if(all(genotypes == c(-1, 0,  0,  0, 1)))  design = c(-1, -1/2,  0,  0,  1/2,  0, -1/4)
    else if(all(genotypes == c( 0, 1,  1, -1, 0)))  design = c( 0,  1/2,  0, -1, -1/2,  0, -1/4)
    else if(all(genotypes == c(-1, 0,  0, -1, 0)))  design = c( 0, -1/2,  0,  0, -1/2, -1,  1/4)
    else                                            design = rep(NA, 7)
    names(design) = paste(c("a_0", "d_0", "i_0", "a_m", "d_m", "c_m0", "dd_m0"), locus, sep = '_')
    return(design)
}

runCromossome <- function(cromossome){
    current_mouse_gen = mouse_phen_std$ID
    num_loci = (length(mouse_gen[[cromossome]])-1)/3
    for(locus in 1:num_loci){
        current_mouse_gen = cbind(current_mouse_gen,
                                  adply(mouse_phen_std$ID, 1, get_design, locus, cromossome)[,-1])
    }
    current_data = na.omit(merge(mouse_phen_std, current_mouse_gen, by.x = 'ID', by.y = 'current_mouse_gen'))

    num_traits = 7
    value = paste("cbind(",
                  paste(paste("grow",
                              paste(1:num_traits, 2:(num_traits+1), sep = ''),
                              sep = ''), collapse = ', '),
                  ")", sep = '')

    fixed.effects = "trait:SEX + trait:LSB + trait:LSW + trait:COHORT - 1"

    null_formula = paste(value, fixed.effects, sep = ' ~ ')

    runSingleLocusMCMCModel <- function(locus){
        genotype_formula = paste(null_formula,
                                 paste(paste('trait:',
                                             paste(c("a_0", "d_0", "i_0",
                                                     "a_m", "d_m",
                                                     "c_m0", "dd_m0"),
                                                   locus,
                                                   sep = '_'),
                                             sep = ''),
                                       collapse = ' + '),
                                 sep = ' + ')
        prior = list(R = list(V = diag(num_traits), n = 0.002),
                     G = list(G1 = list(V = diag(num_traits) * 0.02, n = 8)))
        mcmc.mouse.model = MCMCglmm(as.formula(genotype_formula),
                                    random = ~us(trait):FAMILY,
                                    data = current_data,
                                    rcov = ~us(trait):units,
                                    family = rep("gaussian", num_traits),
                                    prior = prior,
                                    verbose = FALSE)
        return(mcmc.mouse.model)
    }
    cromossome_model_list = alply(1:num_loci, 1, runSingleLocusMCMCModel, .parallel = TRUE)
    return(cromossome_model_list)
}
maternal_scan = llply(names(mouse_gen), runCromossome)
names(maternal_scan) = names(mouse_gen)
save(maternal_scan, file = "./Rdatas/maternalScan_MCMCglmm.Rdata")
#load("./Rdatas/maternalScan_MCMCglmm.Rdata")


