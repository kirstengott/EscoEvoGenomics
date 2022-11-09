
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)

tree = "/Users/kirstengotting/Analysis/EscoEvoGenomicsArchive/timetree/ant_dates_mcmc/prune_2021_08_16/test/Pbruchi_MC1_SN_mcmctree_combr1-4_rename.tre"

phy <- readMCMCtree(tree, from.file = TRUE)



pdf(file = 'plots/MCMCtreeAnt.pdf')
MCMC.tree.plot(phy, analysis.type = "MCMCtree", cex.tips = 0.3,
               time.correction = 100, plot.type = "phylogram", lwd.bar = 2,
               scale.res = c("Epoch"), node.method = "bar", col.age = "#528DF7",
               no.margin = TRUE, label.offset = 1)
dev.off()


tree2 = "/Users/kirstengotting/Analysis/EscoEvoGenomicsArchive/timetree/esco_dates_mcmc/remcmctreeresults___finally/remcmctreeresults___finally/combined_r1-r6_FigTree.tre"
phy2 <- readMCMCtree(tree2, from.file = TRUE)


pdf(file = 'plots/MCMCtreeEscovopsis.pdf')
MCMC.tree.plot(phy2, analysis.type = "MCMCtree", 
                       cex.tips = 0.3,
                       time.correction = 100, 
                       plot.type = "phylogram", 
                       lwd.bar = 2,
                       scale.res = c("Epoch"), 
                       node.method = "bar", 
                       col.age = "#528DF7",
                       no.margin = TRUE, 
                       label.offset = 1)


dev.off()


MCMC.tree.plot(phy2, analysis.type = "MCMCtree", 
               cex.tips = 0.3,
               time.correction = 100, 
               plot.type = "phylogram", 
               lwd.bar = 2,
               scale.res = c("Epoch"), 
               node.method = "bar", 
               col.age = "#528DF7",
               no.margin = TRUE, 
               label.offset = 1)
