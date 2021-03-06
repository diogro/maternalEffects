library(gdata)
library(plyr)
library(dplyr)
library(reshape2)

raw_mouse_phen = read.csv("./data/F3phenotypes_uncorrected.csv", as.is = T)
raw_mouse_phen = select(raw_mouse_phen, c(ID, FATPAD:PAIR))

f3_pedigree = read.csv("./data/F3_pedigree.csv")
f3_pedigree = f3_pedigree[-c(1997:2000),]

raw_mouse_phen = mutate(raw_mouse_phen,
                          gr12 = WEEK2 - WEEK1,
                          gr23 = WEEK3 - WEEK2,
                          gr34 = WEEK4 - WEEK3,
                          gr45 = WEEK5 - WEEK4,
                          gr56 = WEEK6 - WEEK5,
                          gr67 = WEEK7 - WEEK6,
                          gr78 = WEEK8 - WEEK7)

raw_mouse_meta = read.csv("./data/F3Phenotypes_further corrected family data_corrected litter sizes.csv", as.is = T)
names(raw_mouse_meta) = gsub('SexAN', 'SEX', names(raw_mouse_meta))
raw_mouse_meta = select(raw_mouse_meta, ID:COHORT)

#selecting pups that where NOT cross-fostered
raw_mouse_meta = raw_mouse_meta[raw_mouse_meta$NURSE == raw_mouse_meta$Dam,]

raw_mouse_meta = raw_mouse_meta[raw_mouse_meta$ID %in% raw_mouse_phen$ID,]
raw_mouse_phen = raw_mouse_phen[raw_mouse_phen$ID %in% raw_mouse_meta$ID,]

mouse_phen = data.frame(arrange(raw_mouse_meta, ID), arrange(raw_mouse_phen, ID)[,-1])
raw_mouse_gen = llply(paste0("./data/genotypes/chrom", 1:19, ".csv"), read.csv, as.is = TRUE)
names(raw_mouse_gen) = paste0("chrom", 1:19)
mouse_gen = raw_mouse_gen
mouse_phen = mouse_phen[mouse_phen$ID %in% mouse_gen[[1]]$ID,]

mouse_phen = mouse_phen[mouse_phen$ID %in% f3_pedigree$ID,]

sum((na.omit(f3_pedigree[, 2]) %in% f3_pedigree[, 1]) == FALSE) > 0 & any(is.na(f3_pedigree[, 2]))
sum((na.omit(f3_pedigree[, 3]) %in% f3_pedigree[, 1]) == FALSE) > 0 & any(is.na(f3_pedigree[, 3]))

rm(list = ls(pattern='raw'))

mouse_phen = select(mouse_phen, ID, FAMILY, SEX, LSB, LSW, COHORT, gr12:gr78)
complete_rows = complete.cases(mouse_phen)
mouse_phen = mouse_phen[complete_rows,]

num_traits = 7
traits = c( "gr12", "gr23", "gr34", "gr45", "gr56", "gr67", "gr78")

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
