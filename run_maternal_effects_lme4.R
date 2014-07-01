source('./read_mouse_data.R')

library(lme4)
library(lmerTest)
library(ggplot2)
library(doMC)

registerDoMC(10)

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
        current_mouse_gen = cbind(current_mouse_gen, adply(mouse_phen_std$ID, 1, get_design, locus, cromossome)[,-1])
    }
    current_data = na.omit(merge(mouse_phen_std, current_mouse_gen, by.x = 'ID', by.y = 'current_mouse_gen'))
    melt_data = melt(current_data, id.vars = names(current_data)[c(1:6, 14:length(names(current_data)))])

    #ggplot(melt_data, aes(variable, value, color = SEX)) + geom_boxplot() + facet_wrap(~SEX)

    null_formula = "value ~ 1 + (0 + variable|FAMILY)"
    mouse_model_no_gen = lmer(as.formula(null_formula),
                              data = melt_data,
                              REML = FALSE)
    G_lme4 = VarCorr(mouse_model_no_gen)[[1]]
    attr(G_lme4,"correlation") = NULL

    runSingleLocusModel <- function(locus, cromossome){
        genotype_formula = paste(null_formula,
                                 paste(paste('variable*',  paste(c("a_0", "d_0", "i_0", "a_m", "d_m", "c_m0", "dd_m0"), locus, sep = '_'), sep = ''), collapse = ' + '),
                                 sep = ' + ')
        mouse_model = lmer(as.formula(genotype_formula),
                           data = melt_data,
                           REML = FALSE)
        test = anova(mouse_model_no_gen , mouse_model)
        print(test)
        return(list(model = mouse_model, anova = test, p.value = test$'Pr(>Chisq)'[2]))
    }
    cromossome_model_list = alply(1:num_loci, 1, runSingleLocusModel, cromossome, .parallel = TRUE)
    return(cromossome_model_list)
}
maternal_scan = llply(names(mouse_gen), runCromossome)
