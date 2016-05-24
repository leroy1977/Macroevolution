###### code for reconstructing ancestral states for morphological traits
### using models of maximum likelihood (Schluter 1997) and squared parsimony (removing branch length effects) and linear parsimony

cov.neotamias ## covariance matrices for each species
tamias.tree   ### phylogenetic tree for Neotamias
climate.means ### means of the climate variables for each species
morpho.mean.neotamias ## means of the morphological traits

tree.sq.parsimony <- compute.brlen(tamias.tree, 1) 


#########
#####################
##### for morphological traits
#####################

#### w.matrices pooling for sample sizes along the phylogeny
cov.neotamias ### covariance matrices estimated for each species. 
N_amostral.21 ## vector of sample sizes

require(evolqg)
w.cov <- PhyloW(tamias.tree, tip.data= cov.neotamias, tip.sample.size= N_amostral.21)

##########################
############################
morpho.mean.neotamias ### mean values for each of the 38 traits in neotamias

#### Maximum likelihood reconstruction

MLreconstruction<- function(tree, medias)
{
  MLrec<- matrix(NA, nrow= 19, ncol= ncol(medias))
  for ( i in 1: ncol(medias))
  {
    MLrec[,i] <- ace(medias[,i], tree, type= "continuous", method= "ML")$ace
    rownames(MLrec) <- names(ace(medias[,i],tree, type= "continuous", method= "ML")$ace)
  }
  MLrec <- rbind(medias, MLrec)
  return(MLrec)
}

ML.morpho.rec <- MLreconstruction(tamias.tree, morpho.mean.neotamias)

#### estimating delta z in each branch:
ML.morpho.rec <- rbind(ML.morpho.rec[order(rownames(ML.morpho.rec)[1:20]),], ML.morpho.rec[21:39,])


delta.z <-matrix(NA, nrow=38, ncol=38)
for(i in 1:38){
  delta.z[i,] <- (ML.morpho.rec[tamias.tree$edge[i,2],]- ML.morpho.rec[tamias.tree$edge[i,1],] )
}
 rownames(delta.z) <- paste(tamias.tree$edge[,2], tamias.tree$edge[,1], sep= "_")
 
############# let's standartize delta z by the mean value for each trait and species

delta.z.mi <- delta.z
for(i in 1:38){
  delta.z.mi[i,] <- delta.z[i,]/ ML.morpho.rec[tamias.tree$edge[i,1],]  
}


## lets calculate W-matrix standartized by the mean of each trait

w.cov <- w.cov[rownames(ML.morpho.rec)]
outer.means <- list()

for ( i in 1: 39){
  outer.means[[i]] <- ML.morpho.rec[i,] %o% ML.morpho.rec[i,]
}
outer.means[[1]]
names(outer.means) <- rownames(ML.morpho.rec)
########
names(outer.means)
names(w.cov)


Gmi <- w.cov
Gmi[[1]]
names(w.cov)

for(i in 1: length(w.cov))
{Gmi[[i]] <- w.cov[[i]]/outer.means[[i]]}

### Gmi noise controlled according to Marroig et al., 2012
require(evolqg)
Gmi.noise <- list()
for (i in 1: length(Gmi))
{
  Gmi.noise[[i]] <- ExtendMatrix(Gmi[[i]],ret.dim=5) ##### retaning 5 dimensions
}
names(Gmi.noise) <- names(Gmi)
#### Gmi.noise Gmi controlada para ruido!


### estimating the selection gradients
beta.mi <- delta.z.mi

for(i in 1: nrow(beta.mi)){
  beta.mi[i,]<-  delta.z.mi[i,] %*% solve(Gmi.noise[[tamias.tree$edge[i,1]]]$ExtMat)   
}

rownames(beta.mi) <- paste(tamias.tree$edge[,2], tamias.tree$edge[,1], sep= "_")

beta.mag <- apply(beta.mi, 1, Norm) ## magnitude of selection


############################

beta.mag

############################################
### SQUARED PARSIMONY RECONSTRUCTION- non weighted by branch length 

ML.morpho.rec.par <- MLreconstruction(tree.sq.parsimony, morpho.mean.neotamias)

#### estimating delta z in each branch:
ML.morpho.rec.par <- rbind(ML.morpho.rec.par[order(rownames(ML.morpho.rec.par)[1:20]),], ML.morpho.rec.par[21:39,])


delta.z.par <-matrix(NA, nrow=38, ncol=38)
for(i in 1:38){
  delta.z.par[i,] <- (ML.morpho.rec.par[tamias.tree$edge[i,2],]- ML.morpho.rec.par[tamias.tree$edge[i,1],] )
}
rownames(delta.z.par) <- paste(tamias.tree$edge[,2], tamias.tree$edge[,1], sep= "_")

############# let's standartize delta z by the mean value for each trait and species

delta.z.mi.par <- delta.z.par
for(i in 1:38){
  delta.z.mi.par[i,] <- delta.z.par[i,]/ ML.morpho.rec.par[tamias.tree$edge[i,1],]  
}


### estimating the selection gradients
beta.mi.par <- delta.z.mi.par

for(i in 1: nrow(beta.mi.par)){
  beta.mi.par[i,]<-  delta.z.mi.par[i,] %*% solve(Gmi.noise.par[[tamias.tree$edge[i,1]]]$ExtMat)   
}

rownames(beta.mi.par) <- paste(tamias.tree$edge[,2], tamias.tree$edge[,1], sep= "_")

beta.mag.par <- apply(beta.mi.par, 1, Norm) ## magnitude of selection

##############
##################
####### REconstructing using LINEAR PARSIMONY (estimated in Mesquite)

pl <- read.csv(file.choose(), sep=";")
pl
pl <- pl[order(rownames(pl)),]

#### estimating delta z in each branch:
parsimony <- rbind(ML.morpho.rec[order(rownames(ML.morpho.rec)[1:20]),], pl)

parsimony <- as.matrix(parsimony)


delta.z.pl <-matrix(NA, nrow=38, ncol=38)
for(i in 1:38){
  delta.z.pl[i,] <- c(parsimony[tamias.tree$edge[i,2],]- parsimony[tamias.tree$edge[i,1],])
}

rownames(delta.z.pl) <- paste(tamias.tree$edge[,2], tamias.tree$edge[,1], sep= "_")


############# let's standartize delta z by the mean value for each trait and species

delta.z.mi.pl <- delta.z.pl
for(i in 1:38){
  delta.z.mi.pl[i,] <- delta.z.pl[i,]/parsimony[tamias.tree$edge[i,1],]  
}

### estimating the selection gradients
beta.mi.pl <- delta.z.mi.pl

for(i in 1: nrow(beta.mi.pl)){
  beta.mi.pl[i,]<-  delta.z.mi.pl[i,] %*% solve(Gmi.noise[[tamias.tree$edge[i,1]]]$ExtMat)   
}

rownames(beta.mi.pl) <- paste(tamias.tree$edge[,2], tamias.tree$edge[,1], sep= "_")

beta.mag.pl <- apply(beta.mi.pl, 1, Norm) ## magnitude of selection


cor.test(beta.mag.par, beta.mag.pl)
###################



beta.mag[order(beta.mag)]
beta.mag.par[order(beta.mag.par)]

quartz()

pdf("figure_4_beta.pdf", width=3.5, height=4)
plotTree(tamias.tree, ftype="i", fsize=0.7)
edgelabels(text=round(beta.mag, 1), adj=c(0.5,-0.5), frame="none", cex=0.6)
dev.off()


pdf("figure_4_beta_par.pdf", width=3.5, height=4)
plotTree(tamias.tree, ftype="i", fsize=0.7)
edgelabels(text=round(beta.mag.par, 1), adj=c(0.5,-0.5), frame="none", cex=0.6)
dev.off()
