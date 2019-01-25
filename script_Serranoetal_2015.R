# Script to run multiple Auteur analyses and slice rates by time as presented in Serrano-Serrano et al. 2015 BMC Evol Biol
# August, 2012 
# Glenn Litsios - Martha Serrano-Serrano

library(ape)
library(picante)
library(geiger)
library(phangorn)
library(auteur)
library(phytools)


trees<-read.nexus(file="path_to_trees.tre")
qtraits<-read.csv(file="path_to_traits.csv", h=T)
#rownames of qtraits are species names

# Separating trees and trees into inputs for Auteur

dataset_creation<-function(trees, qtraits){

  if(class(trees)!="multiPhylo"){
		return("Tree is not a 'multiPhylo' class")}

    for (i in 1:length(trees)){
  
      tree<-trees[[i]] # sample a tree from multi.phylo object  
     
      for (k in 1:ncol(qtraits)) {               
	trait<-traits[,k]			# send it variable to a vector	
	names(trait)<-rownames(traits)
	name.check(tree, trait)									
	phylo<-match.phylo.data(tree,trait)  # relate that vector with each tree		
	save(phylo, file=paste("tree_",i,"_",colnames(traits)[k],".input",sep=""))	# a single qtrait-tree object is generated looping on qtraits and trees
      }
    }
  }

# Rates estimations with Auteur
  
 files=dir(pattern=".input") 

 for (i in 1:length(files)) {
  	
  load(files[i])
  r=paste("trait",i ,sep="") 
  lapply(1:2, function(x) rjmcmc.bm(phy=phylo$phy, dat=trait, ngen=100000, sample.freq=100, prob.mergesplit=0.1, simplestart=TRUE, prop.width=1, fileBase=paste(r,x,sep=".")))
  dirs=dir("./",pattern=paste("BM",r,sep="."))
  pool.rjmcmcsamples(base.dirs=dirs, lab=r)
  load(paste(paste(r,"combined.rjmcmc",sep="."),paste(r,"posteriorsamples.rda",sep="."),sep="/"))
  #to plot any single run analyses with standard Auteur output
  pdf(file=paste("trait",i, sep=""))
  shifts.plot(phy=phylo$phy, base.dir=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, edge.width=2, font=7, lab4ut="horizontal", x.lim=1.5)
  dev.off()
  save.image(file=paste(r,".space", sep="_"))

  }

  
# Summary of all estimations and slice by time into a single plot (also works for clade-specific time slicing)

files<- as.list(list.files(pattern=".space"))  # call all the trees
slices<-seq(0.01,1.01, 0.001)  # create time slices
i=1
filename<-files[[i]] 
load(file=filename)
vector_total<-phylo$phy$tip.label    #here a specific set of species or clade could be defined 

results_trait<-matrix(nrow=length(files),ncol=length(slices))
for (i in c(1:length(files))) {
	 print(i)
	 flush.console() 
	 filename<-files[[i]] 
	 load(file=filename)
	 rates=posteriorsamples$rates 				
	 burnin<-round(0.25*nrow(rates)) 
	 rates=rates[-c(1:burnin),]  
	 R.scl=apply(rates, 2, mean)
	 phy=phylo$phy                
	 rat<-treeSlice(phy,slices,rates=R.scl,vector_total) 
	 results_trait[i,]<-rat
    }

    # plot using the median value of the rate over all the trees in a point of time 
    rates_median<-apply(results_trait, 2, median,na.rm=T)
    trait_q25<-apply(results_trait, 2, function(x){
		         quantile(x, na.rm=T)[2]
			  })
    trait_q75<-apply(results_trait, 2, function(x){
                        quantile(x, na.rm=T)[4]
                         })

    plot(1,1, type="n", xlim=c(0,1), ylim=c(0,max(trait_q75), ylab="Trait Median Rates", xlab="Time")
    points(slices, rates_median, type="l", col="black", lwd=3)
    points(slices, trait_q25, type="l", col="gray", lwd=2, lty=2)
    points(slices, trait_q75, type="l", col="gray", lwd=2, lty=2)
    
  
# Run these two functions before execute the combination of plots. 
  
# Modified TreeSlice function 

treeSlice<-function(phy,slices,rates=rates,vector){
    if(class(phy)!="phylo") stop("tree should be object of class 'phylo'.")

    phy<-reorder(phy)   
    root<-length(phy$tip)+1
    node.height<-matrix(NA,nrow(phy$edge),3) 
        for(i in 1:nrow(phy$edge)){
        if(phy$edge[i,1]==root){
            node.height[i,1]<-0.0
            node.height[i,2]<-phy$edge.length[i]
        } else {
            node.height[i,1]<-node.height[match(phy$edge[i,1],phy$edge[,2]),2]
            node.height[i,2]<-node.height[i,1]+phy$edge.length[i]
               }
	}
    nodex<-1:length(node.height[,1])  
    node.height[,3]<-nodex   
    group<-match(vector, phylo$phy$tip.label) 		
    mrca<-oldest.mrca(phy, tips=group) 
    mrca<-append(mrca, Descendants(phy,mrca, type="all")) 
    new_edges<-which(phylo$phy$edge[,1]%in% mrca)     
    new_node.height<-node.height[(new_edges),]  
    
      ra<-vector("numeric", length=length(slices))
        for(i in 1:length(slices)){
            edges<-which(new_node.height[,2]>slices[i]&new_node.height[,1]<slices[i])
            edges<-new_node.height[edges,3]
	      if(length(edges)>0){
		ra[i]<-mean(rates[edges],na.rm=T)
		  }else{ra[i]<-NA}
				    }
        return(ra)
	}   
      
  

# Function for mrca from Phytools 


oldest.mrca<-function(tree,tips){
   H<-nodeHeights(tree)
   X<-mrca(tree)
   n<-length(tips)
   nodes<-height<-vector(); k<-1
   for(i in 1:(n-1)) for(j in (i+1):n){
      nodes[k]<-X[tips[i],tips[j]]
      height[k]<-H[match(nodes[k],tree$edge[,1]),1]
      k<-k+1
   }
   z<-match(min(height),height)
   return(nodes[z])
}
