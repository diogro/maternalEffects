library(gdata)
library(plyr)
library(dplyr)
library(reshape2)

raw_mouse_phen = read.csv("./data/F3phenotypes_uncorrected.csv", as.is = T)
raw_mouse_phen = select(raw_mouse_phen, c(ID, FATPAD:PAIR))

raw_mouse_phen = mutate(raw_mouse_phen,
                          grow12 = WEEK2 - WEEK1,
                          grow23 = WEEK3 - WEEK2,
                          grow34 = WEEK4 - WEEK3,
                          grow45 = WEEK5 - WEEK4,
                          grow56 = WEEK6 - WEEK5,
                          grow67 = WEEK7 - WEEK6,
                          grow78 = WEEK8 - WEEK7)

raw_mouse_meta = read.csv("./data/F3Phenotypes_further corrected family data_corrected litter sizes.csv", as.is = T)
names(raw_mouse_meta) = gsub('SexAN', 'SEX', names(raw_mouse_meta))
raw_mouse_meta = select(raw_mouse_meta, ID:COHORT)

raw_mouse_meta = raw_mouse_meta[raw_mouse_meta$ID %in% raw_mouse_phen$ID,]
raw_mouse_phen = raw_mouse_phen[raw_mouse_phen$ID %in% raw_mouse_meta$ID,]

mouse_phen = data.frame(arrange(raw_mouse_meta, ID), arrange(raw_mouse_phen, ID)[,-1])
raw_mouse_gen = llply(paste0("./data/genotypes/chrom", 1:19, ".csv"), read.csv, as.is = TRUE)
names(raw_mouse_gen) = paste0("chrom", 1:19)
mouse_gen = raw_mouse_gen
mouse_phen = mouse_phen[mouse_phen$ID %in% mouse_gen[[1]]$ID,]

rm(list = ls(pattern='raw'))

mouse_phen = select(mouse_phen, ID, FAMILY, SEX, LSB, LSW, COHORT, grow12:grow78)
complete_rows = complete.cases(mouse_phen)
mouse_phen = mouse_phen[complete_rows,]

num_traits = 7
traits = c( "grow12", "grow23", "grow34", "grow45", "grow56", "grow67", "grow78")

m_data = melt(mouse_phen, id.vars = names(mouse_phen)[1:6])

null.formula = "value ~ 1 + variable * SEX + variable * LSB + variable * LSW + variable * COHORT"
mouse_no_fixed = lm(as.formula(null.formula), data = m_data)
m_data_std = m_data
m_data_std$value = residuals(mouse_no_fixed)
names(m_data_std)

exclude = c(dim(m_data_std)[2]-1, dim(m_data_std)[2])
cast_formula = paste(paste(names(m_data_std[,-exclude]), collapse = " + "),
                     'variable',
                     sep = " ~ ")
mouse_phen_std = dcast(m_data_std, as.formula(cast_formula))
