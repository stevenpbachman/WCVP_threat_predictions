### Extracting phylogenetic eigenvectors
library(ape)

# get Felix's phylogeny

tree <- read.tree("ANGIOtrees100.tre")
# tree1 <- tree[[1]] # get the #1 tree out of 100
# Phylogenetic tree with 329798 tips and 74544 internal nodes.

## reduce to a genus tree
### adapted from http://blog.phytools.org/2014/11/pruning-trees-to-one-member-per-genus.html
## get a list of all genera
tips<-tree[[1]]$tip.label
genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
head(genera)
length(genera)

## now drop all tips but one of each genus
# be careful with the non-monophyletic genera 
genera.order<-sapply(strsplit(tips,"_"),function(x) x[1])
ii<-sapply(genera,function(x,y) match(x,y)[1],y=genera.order)
tree1d<-drop.tip(tree[[1]],setdiff(tree[[1]]$tip.label,tips[ii]))

# rename tips
tree1d$tip.label<-sapply(strsplit(tree1d$tip.label,"_"),function(x) x[1])
genus.tree <- tree1d

# extract eigenvectors using PVRdecomp
# library(PVR)
system.time(pvrd <- PVRdecomp(genus.tree, scale = TRUE)) # started 22:06 100822
ev <- pvrd@Eigen$values
evec <- pvrd@Eigen$vectors

# evec.df <- as.data.frame(evec[,])
# write.csv(evec.df,"eigenvectors_0.1.csv")
# 
# ev.df <- as.data.frame(ev)
# write.csv(ev.df,"eigenvalues_0.1.csv")

# Alternative and simpler way to get the eigenvectors (Swenson 2014)
#(only for Euclidian distance matrix)
# Extract the distance matrix from the phylogenetic tree
system.time(phylo.dist <- cophenetic(genus.tree))
system.time(pca <- princomp(phylo.dist))
evec.pca <- pca$scores

# ## try with a smaller tree first - estimate time
# # subsetting the tree
# drop <- genus.tree$tip.label[!genus.tree$tip.label %in% genera[1:4000]]
# genus.tree.sub <- drop.tip(genus.tree, as.character(drop))
# 
# system.time(pvrd <- PVRdecomp(genus.tree.sub, scale = TRUE) ) # 205.11s
# system.time(phylo.dist <- cophenetic(genus.tree.sub)) 
# system.time(pca <- princomp(phylo.dist))

# use the Broken Stick method to define the number of eigenvectors
# is this what Diniz-Filho recommend? check this
# Broken stick model (according to Boccard et al 2011)
# We will probably use more recent methods (this is a VERYpreliminary version)
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n)
{
  bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p <- 100*bsm$p/n
bsm

# Plot eigenvalues and % of variance for each axis
windows(title="PCA eigenvalues")
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")  	# average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, 
        main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), 
       pch=15, col=c("bisque",2), bty="n")


# Same plots using a single function: evplot()
# Plot eigenvalues and % of variance for each axis
evplot(ev)

# visualise how the taxa in our phylogeny load on these axes
library(adephylo)
library(phylobase)
# use the phylo4d() function to merge our PC axis data with the phylogenetic tree
obj4d <- phylo4d(genus.tree, evec[,1:10]) # observe that this example uses only the first two principal components

# use the table.phylo4d() function to plot the phylogeny with the "trait"
table.phylo4d(obj4d, cex.sym =1, cex.label = .5, ratio=.5,
              cex.legend=1, col=heat.colors(100),
              show.tip.label=F, show.node.label=F)

# table.phylo4d(obj4d, ratio=.5)
# bullseye(obj4d)


# Plot eigenvalues and % of variance for each axis
windows(title="PCA eigenvalues")
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")  	# average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev[1:100]/sum(ev[1:100]),bsm$p[n:(n-99)])), beside=TRUE, 
        main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), 
       pch=15, col=c("bisque",2), bty="n")

# cumulative percentage explained by the first components
s <- 100*ev/sum(ev)
sum(s[1:1000])
