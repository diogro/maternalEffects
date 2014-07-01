source('./read_mouse_data.R')

library(lme4)
library(ggplot2)

ind = mouse_phen_std$ID [1]
locus = 1
cromossome = 1
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
    else {print(genotypes);print(ind);print(mother);design = rep(NA, 7)}
    names(design) = c("a_0", "d_0", "i_0", "a_m", "d_m", "c_m0", "dd_m0")
    return(design)
}
get_design(ind, locus, cromossome)
current_mouse_gen = adply(mouse_phen_std$ID, 1, get_design, locus, cromossome)
design_matrix
na.omit(merge(mouse_phen_std, current_mouse_gen, by.x = 'ID', by.y = 'X1'))
