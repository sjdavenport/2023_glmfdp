a = read.table('Cbeta.csv')
Cbeta = a$V1
b = read.table('orig_pvalues.csv')
orig_pvalues = b$V1

# remotes::install_github("pneuvial/sanssouci")
library(sanssouci)
alpha <- 0.1
lambda <- 0.22 # Obtained from the bootstrap in python

n_genes <- length(orig_pvalues)
thresholds <- sanssouci:::t_linear(lambda, 1:n_genes, n_genes)

pdf("volcano_plot.pdf", width=8, height=5)
sanssouci:::volcanoPlot.numeric(x = Cbeta, 
                                y = orig_pvalues,
                                pval = orig_pvalues,
                                thr = thresholds, 
                                p = 0.001,
                                r = 0.5, 
                                col = c("#33333333", "#202020", "#33333311"))
dev.off()
# not used:  convert -background white -flatten volcano_plot.pdf PNG64:volcano_plot.png
