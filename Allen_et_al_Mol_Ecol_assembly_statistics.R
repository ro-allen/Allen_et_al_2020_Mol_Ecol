# Packages
library(ggplot2)
library(phyloseq)
library(reshape2)
library(ape)
library(vegan)
library(vegan)
library(picante)
library(tidyverse)
library(viridis)

# Community assembly statistics --------------------------------------------------------

# filtering phyloseq object (gt.pro) to remove ASVs with zero reads. 
gt.pro.trim = filter_taxa(gt.pro, function(x) sum(x > 0.5) > (0.0000000000000001*length(x)), prune = TRUE) 

# Phylogenetic signal 
select.env = sample_data(gt.pro.trim)[,10:13] 
wa.asv = wascores(select.env, as.data.frame(otu_table(gt.pro.trim))) 

# euclidean distance between asv's
wa.dist = as.matrix(dist(wa.asv, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)) 
phy.dist = cophenetic.phylo(phy_tree(gt.pro.trim)) 

# mantel correlogram 
gt.correlog = mantel.correlog(wa.dist, phy.dist, r.type = "pearson") 
plot(gt.correlog, alpha = 0.01) 


# SES.MNTD 

# interspecific distance matrix to feed BMNTD alongside the community matrix 
otu.mat.2 = as.data.frame(t(otu_table(gt.pro.trim)))
gt.match = match.phylo.data(phy = phy_tree(gt.pro.trim), data = otu.mat.2) 
identical(colnames(t(gt.match$data)),colnames(cophenetic(gt.match$phy))) # TRUE
identical(colnames(t(gt.match$data)),rownames(cophenetic(gt.match$phy))) # TRUE

# actual BMNTD
gt.mntd = as.matrix(comdistnt(t(gt.match$data), cophenetic(gt.match$phy), abundance.weighted = TRUE)) # bmntd calculation

# BNTI  ---------------------------------------------------------

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(gt.match$data),ncol(gt.match$data),beta.reps));

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(gt.match$data),taxaShuffle(cophenetic(gt.match$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(gt.match$data),ncol=ncol(gt.match$data));

for (columns in 1:(ncol(gt.match$data)-1)) {
  for (rows in (columns+1):ncol(gt.match$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (gt.mntd[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
}

rownames(weighted.bNTI) = colnames(gt.match$data)
colnames(weighted.bNTI) = colnames(gt.match$data)

# histogram of results
hist.decomp = as.vector(weighted.bNTI)
hist.decomp = as.data.frame(hist.decomp)
colnames(hist.decomp) = c("Distribution")
head(hist.decomp)
hist.decomp$process[hist.decomp$Distribution>=-2] = "Selection Not Dominant"
hist.decomp$process[hist.decomp$Distribution>=2] = "Variable Selection"
hist.decomp$process[hist.decomp$Distribution<=-2] = "Homogenizing Selection"

process.histogram = ggplot(hist.decomp, aes(Distribution, fill = process)) + geom_histogram(bins = 100) + scale_fill_viridis(discrete = T, end = 0.8) + theme_bw() 

#end
# RCbray ------------------------------------------------------------------

raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_bray_curtis<-NULL
      for(i in 1:reps){
        
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); 
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; 
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); 
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; 
        rm('com1.samp.sp','com1.sp.counts');			
        
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); 
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; 
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); 
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; 
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); 
        
        null_bray_curtis[i] = vegdist(null.spXsite,method="bray");
        
      };
      
      obs.bray = vegdist(spXsite[c(null.one,null.two),],method="bray");
      
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); 
      
      print(c(null.one,null.two,date()));
      
    }; 
    
  }; 
  
  if(as.distance.matrix){ 
    results<-as.dist(results)
  }	
  
  return(results)
  
} 

# run rcBray 
gt.rc.bray = raup_crick_abundance(as.data.frame(otu_table(gt.pro.trim)), plot_names_in_col1 = FALSE, classic_metric = FALSE, split_ties = TRUE, reps = 999, set_all_species_equal = FALSE, as.distance.matrix = TRUE, report_similarity = FALSE)
