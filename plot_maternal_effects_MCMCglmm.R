source('./run_maternal_effects_MCMCglmm.R')

getEffects(maternal_scan_mcmc[[1]][[1]])
model_summary <- ldply(maternal_scan_mcmc, function(x) ldply(x, getEffects, .parallel = TRUE))
names(model_summary)[1:2] <- c('chrom', 'loci')
model_summary$chrom <- factor(model_summary$chrom, levels = paste0('chrom', 1:19))

x <- 'gr12'
plots_trait <- llply(unique(model_summary$trait),
                     function(x) {
                         ggplot(filter(model_summary, effect == 'd', trait == x), aes(loci, mean)) +
                             facet_wrap(~chrom, scale = 'free_x') +
                             geom_point() +
                             geom_errorbar(aes(ymin = lower, ymax = upper)) + theme_bw() + geom_hline(h = 0) + theme(title = element_text(x))
                     })
