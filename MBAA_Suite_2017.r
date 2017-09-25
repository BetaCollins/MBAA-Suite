#-------------------------------------------------------------------------------
# Author: Chu Wenhan Collins (updated 30th June 2017)
#-------------------------------------------------------------------------------
# R source that contains the analysis functions of the new 16S analysis
# DO NOT ATTEMPT TO EDIT THIS SCRIPT, functions may not work !
#
# Annotation files are needed! make sure annnotation file is in the correct format
# as described in the main driver source script.
#
# Set the correct working directory containing the annnotation file!
#-------------------------------------------------------------------------------

library(vegan);library(ggplot2);library(reshape2);library(dplyr); library(qvalue); library(psych); library(RColorBrewer)
library(gplots); library(ggrepel)

setwd("R:/ID/ID1/Collins/MBAA/MBAA_Suite_2017")
param_settings=read.table(file="MBAA_Suite_param.txt",header=F,sep='=',stringsAsFactors=F)

#---------------------------------------------------------------------------------#
### Create Relative Abundance table with NA replaced by higher tax
#---------------------------------------------------------------------------------#
if(param_settings[which(param_settings[,1]=="create_RelAbunTable"),2]=="T"){
  setwd(param_settings[which(param_settings[,1]=="files_dir"),2])
  samplelist=param_settings[which(param_settings[,1]=="samplelist"),2]
  
  fname_suffix="raw-table"
  flist=list.files("./",pattern=fname_suffix)
  tmp=unlist(strsplit(flist,'_'))
  samplename=tmp[seq(from=1,to=length(tmp),by=2)]
  samplename=cbind(samplename,flist)
  sampleaff=read.table(file=samplelist,header=T,stringsAsFactors=F,sep='\t',check.names=F)
  sampleinfo=merge(samplename,sampleaff,by.x=1,by.y=1,sort=F)
  colnames(sampleinfo)[1:2]=c("samID","fname")
  
  family_rel_tb = data.frame(bac_name=NA);  genus_rel_tb = data.frame(bac_name=NA);  species_rel_tb = data.frame(bac_name=NA)
  
  for (k in 1:nrow(sampleinfo)){
    temp=read.table(as.character(sampleinfo$fname[k]),header=F,stringsAsFactors=F,sep='\t',check.names=F)
    
    #apply similarity threshold
    if(T){
      if(length(which(temp$V4<75))!=0){temp[which(temp$V4<75),c(6:11)]=rep(c("p:","c:","o:","f:","g:","s:"),each=length(which(temp$V4<75)))}
      if(length(which(temp$V4<86.5))!=0){temp[which(temp$V4<86.5),c(9:11)]=rep(c("f:","g:","s:"),each=length(which(temp$V4<86.5)))}
      if(length(which(temp$V4<94.5))!=0){temp[which(temp$V4<94.5),c(10:11)]=rep(c("g:","s:"),each=length(which(temp$V4<94.5)))}
      if(length(which(temp$V4<97))!=0){temp[which(temp$V4<97),c(11)]=rep(c("s:"),each=length(which(temp$V4<97)))}
    }
    
    #change species name
    change_idx = which(temp$V11!="s:")
    temp$V11[change_idx] = paste0("s:",gsub("g:","",temp$V10[change_idx]),"_",gsub("s:","",temp$V11[change_idx]))
    
    #family 
    if(T){
      bac_idx = t(apply(temp,1, FUN = function(x){match(c("p:","c:","o:","f:"),x[6:9])}))
      bac_idx = apply(bac_idx,1, FUN = function(x){length(which(is.na(x)))})
      temp_abun = temp[,c(3,2)]
      for(n in 1:nrow(temp)){temp_abun$V3[n] = temp[n,(bac_idx[n]+5)]}
      temp_abun = aggregate(V2~V3,data = temp_abun, sum)
      colnames(temp_abun) = c("bac_name",as.character(sampleinfo$samID[k]))
      family_rel_tb = merge(family_rel_tb, temp_abun, by.x = 1, by.y = 1, all = T)
    }
    #genus 
    if(T){
      bac_idx = t(apply(temp,1, FUN = function(x){match(c("p:","c:","o:","f:","g:"),x[6:10])}))
      bac_idx = apply(bac_idx,1, FUN = function(x){length(which(is.na(x)))})
      temp_abun = temp[,c(3,2)]
      for(n in 1:nrow(temp)){temp_abun$V3[n] = temp[n,(bac_idx[n]+5)]}
      temp_abun = aggregate(V2~V3,data = temp_abun, sum)
      colnames(temp_abun) = c("bac_name",as.character(sampleinfo$samID[k]))
      genus_rel_tb = merge(genus_rel_tb, temp_abun, by.x = 1, by.y = 1, all = T)
    }
    #species 
    if(T){
      bac_idx = t(apply(temp,1, FUN = function(x){match(c("p:","c:","o:","f:","g:","s:"),x[6:11])}))
      bac_idx = apply(bac_idx,1, FUN = function(x){length(which(is.na(x)))})
      temp_abun = temp[,c(3,2)]
      for(n in 1:nrow(temp)){temp_abun$V3[n] = temp[n,(bac_idx[n]+5)]}
      temp_abun = aggregate(V2~V3,data = temp_abun, sum)
      colnames(temp_abun) = c("bac_name",as.character(sampleinfo$samID[k]))
      species_rel_tb = merge(species_rel_tb, temp_abun, by.x = 1, by.y = 1, all = T)
    }
  }
  
  family_rel_tb = family_rel_tb[-nrow(family_rel_tb),];family_rel_tb[is.na(family_rel_tb)] = 0
  genus_rel_tb = genus_rel_tb[-nrow(genus_rel_tb),];genus_rel_tb[is.na(genus_rel_tb)] = 0
  species_rel_tb = species_rel_tb[-nrow(species_rel_tb),];species_rel_tb[is.na(species_rel_tb)] = 0
  
  family_rel_tb[,-1] = scale(family_rel_tb[,-1], center=FALSE, scale=colSums(family_rel_tb[,-1]))
  genus_rel_tb[,-1] = scale(genus_rel_tb[,-1], center=FALSE, scale=colSums(genus_rel_tb[,-1]))
  species_rel_tb[,-1] = scale(species_rel_tb[,-1], center=FALSE, scale=colSums(species_rel_tb[,-1]))
  
  setwd(param_settings[which(param_settings[,1]=="outdir"),2])
  write.csv(family_rel_tb, "family_RelAbun_NAreplaced.csv", row.names = F)
  write.csv(genus_rel_tb, "genus_RelAbun_NAreplaced.csv", row.names = F)
  write.csv(species_rel_tb, "species_RelAbun_NAreplaced.csv", row.names = F)
  write.csv(sampleinfo, "Metadata.csv", row.names = F)
}


#---------------------------------------------------------------------------------#
### Analysis Starts from Here   !!! 
#---------------------------------------------------------------------------------#
setwd(param_settings[which(param_settings[,1]=="outdir"),2])
level = param_settings[which(param_settings[,1]=="taxon_level"),2]
pri_grp = param_settings[which(param_settings[,1]=="pri_grp"),2]
color_match = as.matrix(read.csv(paste("R:/ID/ID1/Collins/MBAA/Color/all.csv",sep="")))
color_match[,1] = paste0(color_match[,3], color_match[,1])

metadata = read.table(file=param_settings[which(param_settings[,1]=="samplelist"),2],header=T,stringsAsFactors=F,sep='\t',row.names = 1)
relative_abun = read.csv(paste(level,"_RelAbun_NAreplaced.csv",sep = ""),header = T,row.names = 1,stringsAsFactors = F)
relative_abun = data.frame(metadata, t(relative_abun[,rownames(metadata)]))
tp_seq = param_settings[which(param_settings[,1]=="tp_seq"),2]
if(tp_seq!="NULL"){tp_seq = unlist(strsplit(tp_seq,split = ","))}

# Bar plot & Stacked Area plot
#-------------------------------------------------------------------------------------------------------
tmp_abun = data.frame(group = relative_abun[,pri_grp], relative_abun[,-c(1:ncol(metadata))])
tmp_abun = aggregate(. ~ group, tmp_abun, mean)
tmp = as.data.frame(t(tmp_abun[,-1]))
colnames(tmp) = tmp_abun[,"group"]
low_abun_idx = which(rowMeans(tmp)<0.001)
low_abun = colSums(tmp[low_abun_idx,])
tmp = tmp[-low_abun_idx,]
abun_rank = order(rowSums(tmp),decreasing=T)
tmp = tmp[abun_rank,]
sort_idx = c(grep("s.",rownames(tmp),fixed = T), grep("g.",rownames(tmp),fixed = T), grep("f.",rownames(tmp),fixed = T), 
             grep("o.",rownames(tmp),fixed = T), grep("c.",rownames(tmp),fixed = T), grep("p.",rownames(tmp),fixed = T), 
             grep("k.",rownames(tmp),fixed = T))
tmp_plot = as.matrix(rbind(tmp[sort_idx,],low_abun))
rownames(tmp_plot)[nrow(tmp_plot)] = "Low.Abundance"
#test = as.matrix(rownames(tmp_plot)[which(is.na(match(rownames(tmp_plot),color_match[,1])))])

# Bar plot
legcol = 1
if(nrow(tmp_plot)>60&nrow(tmp_plot)<=100){legcol = 2}else if(nrow(tmp_plot)>100){legcol = 3}
jpeg(filename=paste0("barplot_",level,"_",pri_grp,".jpeg"),width=1920,height=1280,quality=100,bg="white",pointsize=24)
barplot(tmp_plot,col=color_match[match(rownames(tmp_plot),color_match[,1]),2],cex.names=0.8,xlim=c(0,(ncol(tmp_plot)+legcol*20)),legend=rownames(tmp_plot),las=2,axis.lty=1,args.legend=list(ncol=legcol,cex=0.65,x = "topright"))
dev.off()


# PCoAplot
#-------------------------------------------------------------------------------------------------------
tmp_abun = relative_abun[,-c(1:ncol(metadata))]
if(length(which(colSums(tmp_abun)==0))!=0) tmp_abun = tmp_abun[,-which(colSums(tmp_abun)==0)]
tmp_abun_cap=capscale(tmp_abun ~ 1, distance="bray")
pcoa1 = round(summary(tmp_abun_cap)$cont$importance[2,1]*100,2)
pcoa2 = round(summary(tmp_abun_cap)$cont$importance[2,2]*100,2)

if(T){
  plot_df = data.frame(tmp_abun_cap$CA$u[,1:2], grp = relative_abun[,pri_grp], samID = rownames(relative_abun))
  p = ggplot(plot_df, aes(x=MDS1,y=MDS2,color=grp)) + geom_point() + 
    #stat_ellipse(type = "norm", linetype = 2) +
    geom_text_repel(aes(label = samID), size = 1.5,vjust=-1) +
    xlab(paste("MDS1  ",pcoa1,"%",sep = "")) + ylab(paste("MDS2  ",pcoa2,"%",sep = "")) + guides(color=guide_legend(title=pri_grp))
  ggsave(paste0("PCoAplot_bray_",pri_grp,".jpeg"), p, width = 6, height = 5)
}

#Constrained Ordinations
if(T){
  vf = envfit(tmp_abun_cap, tmp_abun, perm = 999)
  vf.sign.name = rownames(vf$vectors$arrows)[which(vf$vectors$pvals<0.01)]
  vf.sign = envfit(tmp_abun_cap, tmp_abun[,vf.sign.name], perm = 999, arrow.mul = 0.5)
  
  spp.scrs <- as.data.frame(scores(vf.sign, display = "vectors"))
  spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
  spp.scrs <- spp.scrs[rownames(vf.sign$vectors$arrows)[order(vf.sign$vectors$r,decreasing = T)[1:20]],] #top 10 vectors by r2
  
  plot_df = data.frame(tmp_abun_cap$CA$u[,1:2], grp = relative_abun[,pri_grp])
  mult = max(abs(plot_df[,1:2]))/max(abs(spp.scrs[,1:2]))-0.1
  p = ggplot(plot_df, aes(x=MDS1,y=MDS2)) + geom_point(aes(color=grp)) + guides(color=guide_legend(title=pri_grp)) +
    #stat_ellipse(type = "norm", linetype = 2) +
    xlim(-max(abs(c(plot_df[,1],mult*spp.scrs[,1])))-0.02,max(abs(c(plot_df[,1],mult*spp.scrs[,1])))+0.02) + 
    ylim(-max(abs(c(plot_df[,2],mult*spp.scrs[,2])))-0.02,max(abs(c(plot_df[,2],mult*spp.scrs[,2])))+0.02) +
    geom_segment(data = spp.scrs, aes(x = 0, xend = mult*MDS1, y = 0, yend = mult*MDS2),
                 arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
    geom_text_repel(data = spp.scrs, aes(x = mult*MDS1, y = mult*MDS2, label = Species),size = 2, alpha = 0.6) +
    xlab(paste("MDS1  ",pcoa1,"%",sep = "")) + ylab(paste("MDS2  ",pcoa2,"%",sep = ""))
  ggsave(paste0("PCoAplot_bray_arrow_",pri_grp,".jpeg"), p, width = 6, height = 5)
}

# PERMANOVA
perm_df = data.frame(grp = relative_abun[,pri_grp])
permanova = adonis(tmp_abun ~ grp, method = "bray", data = perm_df)
write.csv(as.matrix(permanova$aov.tab),paste0("PERMANOVA_result_",pri_grp,".csv"),row.names = T)


# Relative Abundance MW test
#-------------------------------------------------------------------------------------------------------
if(length(unique(relative_abun[,pri_grp]))==2){
  tmp_abun = relative_abun[,-c(1:ncol(metadata))]
  tmp_abun = tmp_abun[,order(colSums(tmp_abun),decreasing=T)]
  p_table = data.frame(bac_name = colnames(tmp_abun), pvalue = NA, adj.p.fdr = NA, qvalue = NA)
  grp_idx_1 = which(relative_abun[,pri_grp] == unique(relative_abun[,pri_grp])[1]); grp_idx_2 = which(relative_abun[,pri_grp] == unique(relative_abun[,pri_grp])[2])
  for(i in 1:nrow(p_table)){
    p_table$pvalue[i] = wilcox.test(x = tmp_abun[grp_idx_1,i], y = tmp_abun[grp_idx_2,i], p.adjust.method = "none")$p.value
  }
  p_table$adj.p.fdr = p.adjust(p_table$pvalue,method = "fdr")
  p_table$qvalue = qvalue(p = p_table$pvalue)$qvalues
  write.csv(p_table,paste(level,"_pvalue_MWtest_",pri_grp,".csv",sep=""),quote = F,row.names = F,na = "")
}

# Summary statistics
summary_stat = describeBy(tmp_abun, relative_abun[,pri_grp])
summary_stat = do.call("rbind",summary_stat)
summary_stat = data.frame(grp = gsub("\\..*","", rownames(summary_stat)), Bac_name = rep(colnames(tmp_abun),length(unique(relative_abun[,pri_grp]))), summary_stat[,-1])
write.csv(summary_stat,paste(level,"_summary_stat_",pri_grp,".csv",sep=""),quote = F,row.names = F,na = "")


### alpha diversity
#---------------------------------------------------------------------------------------------
tmp_abun = relative_abun[,-c(1:ncol(metadata))]
alpha_df = data.frame(grp = relative_abun[,pri_grp],shannon = diversity(tmp_abun, index = "shannon"),simpson = diversity(tmp_abun, index = "simpson"))
pval_shannon = pairwise.wilcox.test(x = alpha_df$shannon, g = factor(alpha_df$grp), p.adjust.method = "none" )$p.value
pval_simpson = pairwise.wilcox.test(x = alpha_df$simpson, g = factor(alpha_df$grp), p.adjust.method = "none" )$p.value
alpha_df = melt(alpha_df, measure.vars = c("shannon", "simpson"))
p = ggplot(data = alpha_df, aes(x=grp, y=value, fill=grp)) + geom_boxplot() + 
  xlab(paste0("Shannon p value: ",prettyNum(pval_shannon, digits=4, format="fg"),";  ","Simpson p value: ",prettyNum(pval_simpson, digits=4, format="fg"))) + 
  facet_wrap(~variable, scales = "free_y")
ggsave(paste0(level,"_alpha_diversity_",pri_grp,".jpeg"),p,width = 6,height = 4)


### Heatmap
#---------------------------------------------------------------------------------------------
tmp_abun = relative_abun[,-c(1:ncol(metadata))]
tmp_abun = tmp_abun[,order(colSums(tmp_abun),decreasing=T)]
mypalette=brewer.pal(9,"Set1")
tmp_abun = t(as.matrix(tmp_abun[,1:20]))

jpeg(filename=paste0(level," ",pri_grp," heatmap.jpeg"),width=1920,height=1280,quality=100,bg="white",pointsize=24)
par(mar=c(1, 1, 1, 1) + 0.1)
heatmap.2(tmp_abun,trace="none",density.info="none",ColSideColors=mypalette[unclass(factor(relative_abun[,pri_grp]))],margin=c(5,12))
legend("left",as.character(sort(unique(relative_abun[,pri_grp]))),fill=mypalette,title="Groups",box.col="white")
dev.off()






































