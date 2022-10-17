### Extracting phylogenetic eigenvectors
library(ape)

# get the phylogeny updated by Felix
tree <- read.tree("ANGIOtrees100.tre")
# tree1 <- tree[[1]] # get the #1 tree out of 100
# Phylogenetic tree with 329798 tips and 74544 internal nodes.

## reduce the species phylogeny to a genus level tree
### adapted from http://blog.phytools.org/2014/11/pruning-trees-to-one-member-per-genus.html
## get a list of all genera
tips<-tree[[1]]$tip.label
genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
head(genera)
length(genera)

## drop all tips but one of each genus
# be careful with the non-monophyletic genera 
genera.order<-sapply(strsplit(tips,"_"),function(x) x[1])
ii<-sapply(genera,function(x,y) match(x,y)[1],y=genera.order)
tree1d<-drop.tip(tree[[1]],setdiff(tree[[1]]$tip.label,tips[ii]))

# rename tips
tree1d$tip.label<-sapply(strsplit(tree1d$tip.label,"_"),function(x) x[1])
genus.tree <- tree1d

# extract eigenvectors using PVRdecomp
# library(PVR)
system.time(pvrd <- PVRdecomp(genus.tree, scale = TRUE)) # 6873.39s

# get the vectors and their names into a data frame
evec.df <- data.frame(pvrd@Eigen$vectors,row.names=pvrd@phylo$tip.label)
eval.df <- data.frame(pvrd@Eigen$values)

# write.csv(evec.df,"eigenvectors_0.1.1.csv") # large file
# write.csv(eval.df,"eigenvalues_0.1.1.csv")

# # read files
# evec.df <- read.csv("eigenvectors_0.1.1.csv")[,-1]
# eval.df <- read.csv("eigenvalues_0.1.1.csv")[,-1]

# use the Broken Stick method to select the eigenvectors
# Diniz-Filho, J. A. F., Sant'Ana, C. E. R. d., & Bini, L. M. (1998). An Eigenvector Method for Estimating Phylogenetic Inertia. Evolution, 52(5), 1247-1262. https://doi.org/10.2307/2411294 
# Broken stick model (according to Boccard et al 2011)

# ev <- pvrd@Eigen$values
ev <- eval.df$pvrd.Eigen.values
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n)
{
  bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p <- 100*bsm$p/n
bsm

# get a data frame with selected vectors
View(eval.df)
eval.df$pct.exp <- 100*ev/sum(ev)
eval.df$broken.stick <- sort(bsm$p, decreasing = T)
eval.df$diff <- eval.df$pct.exp - eval.df$broken.stick
sel.eval.df <- eval.df[eval.df$diff > 0,] 

# write.csv(eval.df,"eigenvalues_0.1.1.csv")
# write.csv(sel.eval.df,"eigenvalues_selected_BS_0.2.csv")

# variance explained by the selected vectors
sum(sel.eval.df$pct.exp) # 77% of the variance
nrow(sel.eval.df) # first 252 vectors selected

# subset the eigenvectors
sel.evec.df <- evec.df[,1:252]

# save csv with the selected 252 phylogenetic eigenvectors representing 77% of the variance associated with the phylogenetic tree
# write.csv(sel.evec.df,"eigenvectors_selected_BS_252_0.2.csv")

# plot % of variance explained by number of eigenvectors selected
library(ggplot2)

df.plot <- eval.df
df.plot$order <- 1:nrow(df.plot)
df.plot$pct.exp.cum <- cumsum(df.plot$pct.exp) # cumulative percentage explained by the eigenvectors

ggplot(df.plot[1:1000,], aes(x=order, y=pct.exp.cum))+
  geom_line(colour="red4", size=1) +
  geom_vline(xintercept = 252,colour="blue4",linetype = "dashed") +
  labs(x="Number of Eigenvectors Selected", y="Cumulative Explained Variance (%)") +
  geom_point(aes(x=252, y=77.39656),size = 4, shape = 21, fill = "white") +
  scale_y_continuous(breaks=seq(10,100,by=10)) +
  scale_x_continuous(breaks=c(1,seq(100,1000,by=100))) +
  theme_bw(base_size = 12)

ggsave("Eigen_CumExpVar.png", dpi=600, width=5, height=3.5)
# Caption: Cumulative variance of the angiosperm phylogeny captured by the phylogenetic eigenvectors. The vertical dashed blue line indicates the selection using the Broken Stick method (i.e., the first 252 vectors, explaining 77% of the variance).  


#############
## Note:
## The file "eigenvectors_selected_BS_252_0.2.csv" contains the 252 phylogenetic eigenvectors selected using the Broken Stick method. You import it as a data frame and merge it by genus (first column or row name) with the predictors' table. 
