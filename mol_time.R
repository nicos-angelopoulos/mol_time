library(IRanges)
library("GenomicRanges")
library(mclust)
library(stringr)

#' Molecular timing from CN and SNV profiles
#'
#' During a copy number gain all SNVs already present on the involved allele will be duplicated as well.
#' Consequently, the VAF of all these substitutions will change from 50\% to  ~66\%, being present on two alleles of three. 
#' All mutations occurred anytime on the not duplicated allele and on one of duplicated alleles will have a VAF ~ 33\%.
#' Consequently, all subclonal substitutions will occur on one single allele with a \eqn{VAF < 33\%}.
#' 
#' This functions assigns each substitution to the 2 (pre gain) or 1 (post gain) alleles using the mclust R function. 
#' The input to the function is the corrected VAF profiles.
#' To avoid noise between the different clusters the function uses a threshold below which data points that belong to an alternative
#' cluster with probability larger than the threshold are excluded from the timing calculation. 
#' The timing calculation is based on the assumption that each patient has a constant mutation rate (r).
#' This function estimates the relative time of each chromosomal gain occurrence (“molecular time”; t)
#' using the model and formulas presented in the materials section of our paper.
#'
#' @author: Peter Campbell 
#' @author: Francesco Maura
#' @author: Nicos Angelopoulos
#' @maintainer: Nicos Angelopoulos    (this is not a proper tag)
#'
#' @param res.dir="../res-gth_0.2-hth_0.4-rsd_3"      
#'                             directory for output files (should not exist, else the script stops).
#'                             the directory name is constructed by bits of the parameters. The example shown above 
#'                             is produced if all parameters are the default ones
#'                             it is constructed by bits of the parameters as shown above if the parameters are the default ones
#'
#' @param samples=NA           list of samples for the analysis, else all samples in dataset are analysed
#'
#' @param data.dir="data"      inputs directory, should contain files: 
#'                             shared_clonal_mutations.txt, copy_number.txt and sample_normal_contamination.txt
#'
#' @param grey.thresh=0.2      data points that belong to an alternative cluster (see mclust() docs) with probability 
#'                             > grey.thresh are excluded from the timing calculations
#' 
#' @param hclust.thresh=0.4    level at which clusters are created from the hclust tree
#' 
#' @param seed=3               set random seed at the top of function
#' 
#' @param iter=1000            bootstrap operations
#'
#' @param warn=FALSE           when TRUE it prints messages about removed samples (due to missing CN profile)
#'
#' @param eps=FALSE            when TRUE also produces .eps versions of the _clusters plots
#'
#' @examples                   \dontrun{mol_time(samples="PD26419a")}
#'                             \dontrun{mol_time(samples=c("PD26400a","PD26401a"))}
#'                             \dontrun{mol_time(grey.thresh=0.1, data.dir="data", res.dir="res-thresh_0.1",warn=T)}
#'
#' @note      
#'

mol_time <- function(res.dir="",samples=NA,data.dir="data",grey.thresh=0.2,hclust.thresh=0.4,seed=3,iter=1000,warn=F,eps=F) {

  set.seed(seed)

  if (res.dir=="") {
      res.dir <- paste("../res",paste("gth",grey.thresh,sep="_"),paste("hth",hclust.thresh,sep="_"),paste("rsd",seed,sep="_"),sep="-")
  }
  print( sprintf("output dir: %s", res.dir ) )

  #### file containing all shared clonal mutations for each patients
  cave_all <- read.delim(file.path(data.dir,"shared_clonal_mutations.txt"),sep="\t", header=T, stringsAsFactors = F) 

  # cave_all <- read.delim(file.path(data.dir,"shared_clonal_mutations.txt"),sep="\t", header=T, stringsAsFactors = F) 
  # use the following if you are analysing a post germinal centre lymphoproliferative disorders (eg myeloma)
  cave_igh<- cave_all[cave_all$Chrom == 14 & cave_all$Pos >106000000 &  cave_all$Pos< 107400000 | 
                cave_all$Chrom == 22 & cave_all$Pos >22300000 &  cave_all$Pos< 23270000 | 
                cave_all$Chrom == 2 & cave_all$Pos >89880000 &  cave_all$Pos< 90300000 |
                cave_all$Chrom == 2 & cave_all$Pos >89100000 &  cave_all$Pos< 89700000,]

  code<- rownames(cave_igh)

  cave<- cave_all[!rownames(cave_all) %in% code,]

  # print( "cave_all" )
  # print( cave_all )

  #### file containing all copy_number data adapted 
  ascat_final<- read.delim(file.path(data.dir,"copy_number.txt"),sep="\t", header=T, stringsAsFactors = F) 

  ascat_final_cn<- ascat_final[ascat_final$class!="normal",]

  #### filter copy_number for: 2:0. 2:1, 3:1, 4:1 status
  batt<- ascat_final_cn[ascat_final_cn$class != "deletion" & ascat_final_cn$class !="Whole_chromosome_duplication" & ascat_final_cn$class !="LOH_3" &
                                   ascat_final_cn$class != "LOH_5" & ascat_final_cn$class != "LOH_2"& ascat_final_cn$class != "bi_del" &
                                   ascat_final_cn$class != "gain_6" & ascat_final_cn$class != "gain_7"& ascat_final_cn$class != "bi_del" & 
                                   ascat_final_cn$class != "LOH_4" & ascat_final_cn$class != "LOH_High" & ascat_final_cn$class != "LOH_8" 
                                 & ascat_final_cn$class != "LOH_9" & ascat_final_cn$class != "gain_4" & ascat_final_cn$class != "LOH_6" 
                                 & ascat_final_cn$class != "-" & ascat_final_cn$code_batt != "subclonal",]

  ascat_final_cn2<- subset(batt, sample %in% unique(cave$Sample))

  #### file with the % of contamination for each sample
  conc<- read.delim(file.path(data.dir,"sample_normal_contamination.txt"),sep="\t", header=T,stringsAsFactors = F) 

  if (is.na(samples[1])) {
      sample_list<- sort(unique(ascat_final_cn2$sample))
  } else {
      sample_list <- samples
  }

  # setwd("/home/na11/ac/wtsi/proj/myeloma/mol_time/res.thres.0.1.wG2")
  if (dir.exists(res.dir)) stop(paste("output directory exists: ",res.dir));  # make sure we are not over-riding results
  dir.create(res.dir)
  old.dir <- getwd()
  setwd(res.dir)

  #######  #######  #######  #######  #######  #######  #######  #######  #######
  #######  #######  #######  #######  #######  #######  #######  #######  #######
  #######  #######  #######  #######  #######  #######  #######  #######  #######
  #######
  ####### clusterization part
  #######
  #######  #######  #######  #######  #######  #######  #######  #######  #######
  #######  #######  #######  #######  #######  #######  #######  #######  #######
  #######  #######  #######  #######  #######  #######  #######  #######  #######

  for(j in (1:length(sample_list)))  
  {
    # j=11
    # set debug flag
    sample_code<- sample_list[j]
    print( sprintf("sample: %s", sample_code) )
    conc2<-conc[conc$SAMPLE == sample_code,]
    aberrant<- as.numeric(conc2$ABBR_CELL_FRAC) ### extract the purity of the selected sample
    cave_sample<- cave[cave$Sample == sample_code,]
    ascat_sample_spec_type<- ascat_final_cn2[ascat_final_cn2$sample == sample_code,]
    cave_sample$chr<- paste0("chr",cave_sample$Chrom)
    ascat_sample_spec_type$chr<- paste0("chr",ascat_sample_spec_type$Chrom)
    gr0 = with(cave_sample, GRanges(chr, IRanges(start=(Pos), end=(Pos))))
    gr1 = with(ascat_sample_spec_type, GRanges(chr, IRanges(start=start, end=end)))
    ranges2 <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
    ranges2 <- ranges2[with(ranges2, startB <= startA & endB >= endA),]
    ascat_brass_second<- ranges2[,c("seqnames","startA","startB","endB")] 
    colnames(ascat_brass_second)<- c("chr", "Pos","start","end")
    int2<- merge(ascat_sample_spec_type, ascat_brass_second, by=c("chr","start","end"))
    cave_ascat<- merge(int2, cave_sample, by=c("chr","Pos","Chrom"))
    if(nrow(cave_ascat)>0) {
        ##### create a code specific for each gain in order to filter out all gains with less than 50 SNV
        cave_ascat$ascat_code<-paste(cave_ascat$chr,cave_ascat$start, cave_ascat$end,cave_ascat$class,sep="-") 
        alfa_sample<-as.data.frame(table(cave_ascat$ascat_code))
        alfa_sample2<-alfa_sample[alfa_sample$Freq>=50,]
        gain_mut<-as.character(alfa_sample2$Var1)
        cave_ascat_filt <- subset(cave_ascat, ascat_code %in% gain_mut)
        } else { mol_time_warn(warn,"no cave_ascat %s",sample_code)
                 next  }     # goes to the next value of for loop
    if(nrow(cave_ascat_filt)>0) {
          all_code =NULL
        } else { mol_time_warn(warn,"no cave_ascat_filt %s",sample_code)
               next  }     # goes to the next value of for loop

    # setwd("/nfs/team78pc18/fm6/relapse_MM/molecular_time/plots")
    # setwd("/nfs/users/nfs_n/na11/ac/wtsi/proj/myeloma/mol_time")
    # setwd("/home/na11/ac/wtsi/proj/myeloma/mol_time/res.thres.0.1.wG2")
    DF<-NULL
    cyto<-list()
    mclust_mm<-list()
    chr_segment_list <-unique(cave_ascat_filt$ascat_code)

    # print( sprintf("for sample: %s", sample_code) )
    # print( sprintf("chr_segement_list: %i", length(chr_segment_list)) )
    for(i in (1:length(chr_segment_list)))
    { 
      cn_code<- chr_segment_list[i]
      chr1<- cave_ascat_filt[cave_ascat_filt$ascat_code ==chr_segment_list[i],]
      first<- chr1[,c("Sample", "Chrom","Pos","Ref","Alt","Gene", "PM.Tum","class")]
      first$Pos<- as.numeric(as.character(first$Pos))
      first$ccf<- (first$PM.Tum /aberrant)*100
      X<- first$ccf
      #mod5 = densityMclust(X)
      #plot(mod5, what = "density", type = "image", col = "dodgerblue3")  ##### add this if you want to see the density plot of the ccf
      myData <- data.frame(x=X)
      # get parameters for most optimal model. nicos: with at least 2 clusters

      if(unique(first$class) == "LOH") ###### single gain
         { myMclust <- Mclust(myData,G=2,verbose=FALSE)
            }
            else { myMclust <- Mclust(myData,G=2:9,verbose=FALSE) 
            }
      # add a column in myData CLUST with the cluster.
      myData$CLUST <- myMclust$classification
      colnames(myData)<-"ccf"
      first_def<- as.data.frame.matrix(cbind(first, myData[-1]))
      colnames(first_def)[ncol(first_def)]<- "code"
      grey.clust <- clusters_plot(sample_code, myMclust, first_def, cn_code, grey.thresh, eps=eps)
      first_def <- first_def[first_def$code != grey.clust,]

    ####### CN/allele assignment part
  
      number_clust<-sort(unique(first_def$code))
      median=NULL
   
      ##### estimate the median allelic fraction for each cluster 
    
    # print( sprintf("for cluster sample: %s", sample_code) )
    # # print( sprintf("for number_clust: %i", length(number_clust)) )
      for(h in (1:length(number_clust))) {
           cluster<- first_def[first_def$code == number_clust[h],]
           alfa<- median(cluster$ccf)
           median<- c(median, alfa)
      } # for h
      code_up<-which.max(median)
      a<- table(first_def$code)

      matrix<- as.data.frame.matrix(cbind(median, rownames(a), a))
      colnames(matrix)<-c("median","cluster_number","tot_mut")
      matrix$tot_mut<-as.numeric(as.character(matrix$tot_mut))
      matrix$median<-as.numeric(as.character(matrix$median))
      matrix<-matrix[order(matrix$median),]
    
      cyto[[i]] <- clust_assign(matrix,first_def,cn_code)

      mclust_mm[[i]]<-first_def


    }   # for i
  
    # merge results from all chromosomes  and create a data frame with all SNVs and clusters

    mclust_mm_final<-     as.data.frame.matrix(do.call("rbind", mclust_mm), stringsAsFactors = F) 
    #write.table(mclust_mm_final,sprintf("%s_with_cluster.txt",sample_code), quote=F, sep="\t", row.names=F, col.names = T)
  
    all_code <- as.data.frame.matrix(do.call("rbind", cyto), stringsAsFactors = F) # merge results from all chromosomes  
    all_code[is.na(all_code)]<-0

    out2 <- str_split_fixed((all_code$segment),'-',4)
    all_code$cyto <- out2[,4]
    all_code$chr<-out2[,1]
    all_code$chrom <- gsub("chr","",all_code$chr)
    all_code<- all_code[order(all_code$chrom),]
  
    cols = c(1:5,ncol(all_code))
    all_code[,cols] = apply(all_code[,cols], 2, function(x) as.numeric(as.character(x)))
  
    file_code<- paste0(sample_code,"_cluster_absolute_number_summary.txt")
    #write.table(all_code,file_code,sep="\t", quote=F)
  
  
    ####### calculation part
  
    ####### Peter's script part bootstrap to generate IC and dots: one function for each possible gain combination
    # all_code_final<- as.data.frame.matrix(rbind(gain_1,gain_2,gain_3,loh_1))
    # all_code_final <- mol_time_ic_gen(all_code)$pests
    
    # nicos: introducing _comb as variables that contain the point estimates ($pests) and bootstrap runs ($bootN)
    
    all_code_final_comb <- mol_time_ic_gen(all_code,iter=iter)
    all_code_final <- all_code_final_comb$pests
  
    cols = c(3:ncol(all_code_final))
    all_code_final[,cols] = apply(all_code_final[,cols], 2, function(x) as.numeric(as.character(x)))
    
    out <- str_split_fixed((all_code_final$segment),'-',4) 
    all_code_final$chr<- (out[,1])
    
    type <- str_split_fixed((out[,4]),'_',2) 
    all_code_final$type<- type[,1]
    all_code_final$chrom<-gsub("chr","",all_code_final$chr)
    all_code_final$chrom[all_code_final$chrom=="X"]<-23
    all_code_final$chrom<-as.numeric(as.character(all_code_final$chrom))
    all_code_final<-all_code_final[order(all_code_final$chrom),]
    all_code_final$chrom[all_code_final$chrom==23]<-"X"
  
    ###### summary of pre and post gain mutation loads for each chromosome
  
    file_code<- paste0(sample_code,"_cluster_summary.txt")
    write.table(all_code_final,file_code,sep="\t", quote=F)
  
    mol_time_plot(sample_code,all_code_final);

    #all_code_final_hc<- all_code_final # nicos: not used anywhere ?!
    #all_code_final_hc[is.na(all_code_final_hc)]<-0

    if(nrow(all_code_final)>3){

      second<-all_code_final[ !is.na(all_code_final$plot_dot2),]
      third<-all_code_final[ !is.na(all_code_final$plot_dot3),]
      
      k<- as.matrix(c(all_code_final$plot_dot1, second$plot_dot2, third$plot_dot3))
      rownames(k)<- c(all_code_final$chr, second$chr, third$chr)
      write.table(k,paste0(sample_code,"_h_clust_tbl.txt"),sep="\t", quote=F)
      all_code2<-  (dist(k))
      pdf(sprintf("%s_h_clust.pdf",sample_code), height=6, width=8)
      plot(hclust(all_code2))
      dev.off()
      write.table(all_code_final_comb$boots,paste0(sample_code,"_boots.txt"),sep="\t", quote=F)
      hcls <- character()

      for( i in (1:iter)){

        hmx <- as.matrix( all_code_final_comb$boots[[5+i]] )
        rownames(hmx) <- all_code_final_comb$boots[[5]]
        hds <- dist(hmx)
        hcl <- hclust(hds) 
        hct <- cutree(hcl, h=hclust.thresh)
        hct <- hct[order(names(hct))]   # align along lexico
                                     # we probably need to standarise cluster naming ?
        hcls[i] <- paste( hct, collapse="_" )
       }   # for i

      thcl <- table(hcls)

      max.len <- max( length(hct), length(thcl) )
      hnas <- rep(NA,max.len-length(thcl))

      edf <- data.frame( seg_cnvn=c(names(hct),rep(NA,max.len - length(hct))), hcls=c(names(thcl),hnas), hfreq=c(as.vector(thcl),hnas) )
      write.table(edf,paste0(sample_code,"_hfreqs.txt"),sep="\t", quote=F)

      nms_thcl <- names(thcl)
      cncl <- matrix(nrow=length(nms_thcl),ncol=length(as.integer(strsplit(nms_thcl[1],split="_")[[1]])) )
      # cncl <- matrix()

      for (i in 1:length(nms_thcl)) {
         cncl[i,] <- as.integer(strsplit(nms_thcl[i],split="_")[[1]])
      }

      colnames(cncl) <- names(hct)
      max.len.2 <- max( length(cncl), length(thcl) )
      hnas.2 <- rep(NA,max.len.2-length(thcl))
      hfdf <- data.frame(cncl, hfreq=c(as.vector(thcl),hnas.2))
      write.table(hfdf,paste0(sample_code,"_chrom_clust_freqs.txt"),sep="\t", quote=F)

   }  # IF THEN
 }  # for j
setwd(old.dir)
}  # function


clusters_plot <- function(sample_code, myMclust, first_def, cn_code, grey.thresh, eps=FALSE) {
    code<- unique(first_def$code)
    nof.clusts <- length(code)     #: nicos,    print(first_def$code); print(code);
    # here 17.09.04: resuffle the colours according to the mean ccf value of each cluster
    centroids <- numeric()
    for (t in 1:nof.clusts) {
          first_def_code <- first_def[first_def$code==code[t],]
          centroids[t] <- mean( first_def_code$ccf )
    }
    col<-c("brown3","forestgreen","cornflowerblue","blueviolet")[1:nof.clusts]
    col<-col[order(centroids)]
    col<-c(col,"grey50")
    max.clust <- max(code)
    grey.clust <- max.clust + 1

    par(mar=c(9,6,5,5), xpd=F)
    kk<- paste(sample_code,cn_code,sep="_")
    clust.prbs <- myMclust$z[,code]   # nicos: Mclust can have hidden clusters (remove them as they often are sub-clusts, with high co-prb
    ext.code <- c(code,grey.clust)   # fixme:

    if (is.matrix(clust.prbs) ) {   # fixme:
      for ( i in (1:nrow(clust.prbs)) ) { 
         row.i <- clust.prbs[i,]
         if (length(row.i[row.i>grey.thresh])>1 ) first_def$code[i] <- grey.clust;
      } 
    }

    pdf(sprintf("%s_clusters.pdf",kk), width=10, height=5) # was: "%s_evolution_plot_.pdf"
    clusters_plot_disp( nof.clusts, first_def, ext.code, col, cn_code )
    dev.off()

    if (eps) {
      setEPS()
      postscript(sprintf("%s_clusters.eps",kk),paper="special",family="NimbusSan") # also can try: ComputerModern
      clusters_plot_disp( nof.clusts, first_def, ext.code, col, cn_code )
      dev.off()
    }
    return(grey.clust)
}

clusters_plot_disp <- function( nof.clusts, first_def, ext.code, col, cn_code) {
    for ( z in (1:(nof.clusts+1)))  {
      first_def_col<- first_def[first_def$code == ext.code[z],]
      if (length(first_def_col$Pos) > 0) {
            plot(first_def_col$Pos, first_def_col$ccf, col=col[z],yaxt="n",xaxt="n", main=cn_code,ylab="", xlab="", ylim=c(0,150), pch=20,cex=1.5)
            par(new=TRUE)
         }
    }
    axis(2, at = seq(0,150, by=10),  col.axis = "black", labels=seq(0,150, by=10), lwd.ticks =1, cex.axis = 1.5, las=1)
}

# LOH
tetraploid.2.2.or.UPD.2.0.est <- function(ploidy.1, ploidy.2, bootstrap = TRUE, iter = 1000) {
    retv <- list()
    point.est <- ploidy.2 / (ploidy.2 + ploidy.1 / 2)
    if (bootstrap) {
      total.muts <- ploidy.1 + ploidy.2
      prob.2 <- ploidy.2 / total.muts
      resamps <- rbinom(n=iter, size=total.muts, p=prob.2)
      ests.bs <- resamps / (resamps + (total.muts - resamps)/2)
      ests.95 <- quantile(ests.bs, c(0.025,0.975))
      # return(c(point.est, ests.95))
      retv$pests <- c(point.est,ests.95)
      retv$noe   <- 1
      retv$boot1 <- ests.bs
    }
    else {retv$pests <- point.est}
    
    return(retv)
}
  
# 1 extra gain
triploid.2.1.est <- function(ploidy.1, ploidy.2, bootstrap = TRUE, iter=1000) {
    retv <- list()
    point.est <- ploidy.2 / (ploidy.2 + (ploidy.1 - ploidy.2)/3)
    if (bootstrap) {
      total.muts <- ploidy.1 + ploidy.2
      prob.2 <- ploidy.2 / total.muts
      resamps <- rbinom(n=iter, size=total.muts, p=prob.2)
      ests.bs <- resamps / (resamps + (total.muts - 2*resamps)/3)
      ests.95 <- quantile(ests.bs, c(0.025,0.975))
      # return(c(point.est, ests.95))
      retv$pests <- c(point.est,ests.95)
      retv$noe   <- 1
      retv$boot1 <- ests.bs
    }
    else {retv$pests <- point.est} # return(point.est)}

    return(retv)
}

# 2 extra gain
tetraploid.3.1.est <- function(ploidy.1, ploidy.2, ploidy.3, bootstrap=TRUE, iter=1000) {
    retv <- list()
    point.est.t1 <- ploidy.3 / (ploidy.3 + ploidy.2 + (ploidy.1 - 2*ploidy.2 - ploidy.3) / 4)
    point.est.t2 <- (ploidy.3 + ploidy.2) / (ploidy.3 + ploidy.2 + (ploidy.1 - 2*ploidy.2 - ploidy.3) / 4)
    if (bootstrap) {
      total.muts <- ploidy.1 + ploidy.2 + ploidy.3
      prob.2 <- ploidy.2 / total.muts
      prob.3 <- ploidy.3 / total.muts
      resamps <- rmultinom(n=iter, size=total.muts, p=c(1-prob.2-prob.3, prob.2, prob.3))
      ests.bs.t1 <- resamps[3,] / (resamps[3,] + resamps[2,] + (resamps[1,] - 2*resamps[2,] - resamps[3,])/4)
      ests.bs.t2 <- (resamps[3,] + resamps[2,]) / (resamps[3,] + resamps[2,] + (resamps[1,] - 2*resamps[2,] - resamps[3,])/4)
      ests.95.t1 <- quantile(ests.bs.t1, c(0.025,0.975))
      ests.95.t2 <- quantile(ests.bs.t2, c(0.025,0.975))
      retv$pests <- c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2)
      retv$noe   <- 2
      retv$boot1 <- ests.bs.t1
      retv$boot2 <- ests.bs.t2
      # return(c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2))
    }
    # else {return(c(point.est.t1, point.est.t2))}
    else {retv$pests <- c(point.est.t1, point.est.t2)}

    return(retv)
}

# 3 extra gain
tetraploid.4.1.est <- function(ploidy.1, ploidy.2, ploidy.3, ploidy.4,bootstrap=TRUE, iter=1000) {
    retv <- list()
    point.est.t1 <- (ploidy.4 * 5) / (ploidy.1 + ploidy.2 *2+ ploidy.3*3 +ploidy.4*4)
    point.est.t2 <- ((ploidy.4 +  ploidy.3)* 5) / (ploidy.1 + ploidy.2 *2+ ploidy.3*3 +ploidy.4*4)
    point.est.t3 <- ((ploidy.4 +  ploidy.3 + ploidy.2)* 5) / (ploidy.1 + ploidy.2 *2+ ploidy.3*3 +ploidy.4*4)
    if (bootstrap) {
      total.muts <- ploidy.1 + ploidy.2 + ploidy.3 + ploidy.4
      prob.2 <- ploidy.2 / total.muts
      prob.3 <- ploidy.3 / total.muts
      prob.4 <- ploidy.4 / total.muts
      resamps <- rmultinom(n=iter, size=total.muts, p=c(1-prob.2-prob.3- prob.4, prob.2, prob.3, prob.4))
      
      ests.bs.t1 <- ((resamps[4,]) * 5) / (resamps[1,] + resamps[2,] *2+ resamps[3,]*3 +resamps[4,]*4)
      ests.bs.t2 <- ((resamps[4,] +  resamps[3,])* 5) / (resamps[1,] + resamps[2,] *2+ resamps[3,]*3 +resamps[4,]*4)
      ests.bs.t3 <- ((resamps[4,] +  resamps[3,]+ resamps[2,])* 5) / (resamps[1,] + resamps[2,] *2+ resamps[3,]*3 +resamps[4,]*4)
      
      ests.95.t1 <- quantile(ests.bs.t1, c(0.025,0.975))
      ests.95.t2 <- quantile(ests.bs.t2, c(0.025,0.975))
      ests.95.t3 <- quantile(ests.bs.t3, c(0.025,0.975))
      
      retv$pests <- c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2,point.est.t3, ests.95.t3)
      retv$noe   <- 3
      retv$boot1 <- ests.bs.t1
      retv$boot2 <- ests.bs.t2
      retv$boot3 <- ests.bs.t3

      # return(c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2,point.est.t3, ests.95.t3))
    }
    # else {return(c(point.est.t1, point.est.t2,point.est.t3))}
    else {retv$pests <- c(point.est.t1,point.est.t2,point.est.t3)}

    return(retv)
}

mol_time_ic_gen <- function(all_code, iter=1000) {
###### generate molecular time and IC for gain
  pc<-list()
  bs<-list()  # nicos: bootstraps
  bs.code <- character()
  bs.segm <- character()
  bs.cnvn <- integer()
  bs.chrm <- character()
  bs.chvn <- character()
  j <- 1      # nicos: bootstraps index (accelerates faster than i)
  all_code_sgms <- all_code$segment
  ac_sgms_map <- list()
  for (i in 1:length(all_code_sgms)) {
      j = 1
      tchr = str_split_fixed(all_code_sgms[i],'-',4)[1]
      pt = paste0(tchr,letters[j])
      while (! is.null(ac_sgms_map[[pt]])) {
         j = j + 1  # let 's hope we don't run out of letters
         pt = paste0(tchr,letters[j])
         if (j>10) print(err2)
      }
      ac_sgms_map[[pt]] <- all_code_sgms[i]
  }
  all_code_gain<- all_code[all_code$cyto == "gain",]
  if(nrow(all_code_gain)>0)
  {
    for(i in (1:nrow(all_code_gain)))
    {
      CODE<-c(all_code_gain$cyto[i],all_code_gain$segment[i])
      est.res <- triploid.2.1.est(all_code_gain$cn1_clon[i],all_code_gain$cn2[i],iter=iter)
      pc_ic_row<- est.res$pests
      pc[[i]]<-c(CODE,pc_ic_row)

      for(k in(1:est.res$noe))
      { 
        bootN <- paste0("boot",k)
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,sep="_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    pc_plot_parameter_gain <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    gain_1<- as.data.frame.matrix(cbind(pc_plot_parameter_gain, matrix(, nrow = nrow(pc_plot_parameter_gain), ncol = 6)))
    colnames(gain_1)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {gain_1<-NULL}
  
  ###### generate molecular time and IC for gain_2

  pc<-list()
  all_code_gain_2<- all_code[all_code$cyto == "gain_2",]
  if(nrow(all_code_gain_2)>0)
  {
    for(i in (1:nrow(all_code_gain_2)))
    {
      CODE<-c(all_code_gain_2$cyto[i],all_code_gain_2$segment[i])
      est.res <- tetraploid.3.1.est(all_code_gain_2$cn1_clon[i],all_code_gain_2$cn2[i],all_code_gain_2$cn3[i],iter=iter)
      pc_ic_row<- est.res$pests
      pc[[i]]<-c(CODE,pc_ic_row)

      for(k in(1:est.res$noe))
      { 
        bootN <- paste("boot",k,sep="")
        # bs[[j]] <- c(CODE,est.res[bootN])
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,sep="_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    pc_plot_parameter_gain_2 <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    gain_2<- as.data.frame.matrix(cbind(pc_plot_parameter_gain_2, matrix(, nrow = nrow(pc_plot_parameter_gain_2), ncol = 3)))
    colnames(gain_2)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {gain_2<-NULL}
  ###### generate molecular time and IC for gain_3
  
  pc<-list()
  all_code_gain_3<- all_code[all_code$cyto == "gain_3",]
  if(nrow(all_code_gain_3)>0)
  {
    for(i in (1:nrow(all_code_gain_3)))
    {
      CODE<-c(all_code_gain_3$cyto[i],all_code_gain_3$segment[i])
      est.res   <- tetraploid.4.1.est(all_code_gain_3$cn1_clon[i],all_code_gain_3$cn2[i],all_code_gain_3$cn3[i],all_code_gain_3$cn4[i],iter=iter)
      pc_ic_row <- est.res$pests
      pc[[i]]   <- c(CODE,pc_ic_row)

      for(k in(1:est.res$noe))
      { 
        bootN <- paste("boot",k,sep="")
        # bs[[j]] <- c(CODE,est.res[bootN])
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,"_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    gain_3 <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    colnames(gain_3)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {gain_3<-NULL}
  
  ###### generate molecular time and IC for LOH
  pc<-list()
  all_code_gain_LOH<- all_code[all_code$cyto == "LOH",]
  if(nrow(all_code_gain_LOH)>0)
  {
    for(i in (1:nrow(all_code_gain_LOH)))
    {
      CODE<-c(all_code_gain_LOH$cyto[i],all_code_gain_LOH$segment[i])
      est.res   <- tetraploid.2.2.or.UPD.2.0.est(all_code_gain_LOH$cn1_clon[i],all_code_gain_LOH$cn2[i],iter=iter)
      pc_ic_row <- est.res$pests
      pc[[i]]<-c(CODE,pc_ic_row)

      for(k in(1:est.res$noe))
      { 
        bootN <- paste("boot",k,sep="")
        # bs[[j]] <- c(CODE,est.res[bootN])
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,sep="_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    pc_plot_parameter_loh <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    loh<- as.data.frame.matrix(cbind(pc_plot_parameter_loh))
    loh_1<- as.data.frame.matrix(cbind(pc_plot_parameter_loh, matrix(, nrow = nrow(pc_plot_parameter_loh), ncol = 6)))
    colnames(loh_1)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {loh_1<-NULL}
  
  #### create the final file
  
  # all_code_final<- as.data.frame.matrix(rbind(gain_1,gain_2,gain_3,loh_1))
  retv <- list()
  # rboots <- as.data.frame.matrix(do.call("rbind", bs), stringsAsFactors = F) # cast the list to a matrix
  # retv$boots <- as.data.frame.matrix(cbind(rboots))
  boots <- data.frame( code=bs.code, chrom=bs.chrm, segment=bs.segm, cnvn=bs.cnvn, chrom_cnvn=bs.chvn )
  for( k in (1:iter) ) {
   acc <- numeric()
   for(m in (1:(j-1))) {
      acc <- c(acc,bs[[m]][k])
   }
      boots <- cbind( boots, acc )
  }
  retv$boots <- boots
  retv$pests <- rbind(gain_1,gain_2,gain_3,loh_1)
  return(retv)
}

##### assign cluster to CN1, CN2, CN3, CN4 etc
clust_assign <- function(matrix, first_def, cn_code) {

    if(unique(first_def$class) == "gain") ###### single gain
    {
      if(nrow(matrix) == 2)
      {
        DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],0,0,as.numeric(sum(matrix$tot_mut)),cn_code)))
        colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
        ## cyto[[i]] <- DF
      }else
      { if(nrow(matrix) == 3)
      {
        if(matrix$median[2]<((matrix$median[1] +matrix$median[3])/2)) ##### avarage between the two extreme cluster to assign the one in the middle
        {
          DF <- t(as.matrix(c(matrix$tot_mut[1] + matrix$tot_mut[2],matrix$tot_mut[3],0,0,as.numeric(sum(matrix$tot_mut)),cn_code)))
          colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
          }
        else
          {
            DF <- t(as.matrix(c(matrix$tot_mut[1] , matrix$tot_mut[2] + matrix$tot_mut[3],0,0,as.numeric(sum(matrix$tot_mut)),cn_code)))
            colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
        }
        ## cyto[[i]] <- DF
      }}}
      else
      {
        if(unique(first_def$class) == "gain_2") ###### 2 extra gains
        { if(nrow(matrix) == 2)
        {
          DF <- t(as.matrix(c(matrix$tot_mut[1],0,matrix$tot_mut[2],0,as.numeric(sum(matrix$tot_mut)),cn_code)))
          colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
          ## cyto[[i]] <- DF
        }else
        {if(nrow(matrix) == 3)
        {
          DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],matrix$tot_mut[3],0,as.numeric(sum(matrix$tot_mut)),cn_code)))
          colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
          ## cyto[[i]] <- DF
        }}}
        else
        { if(unique(first_def$class) == "gain_3") ###### 3 extra gains
        { if(nrow(matrix) == 3)
          {
          if((matrix$median[2]>(matrix$median[3]/4))&(matrix$median[2]<(matrix$median[3]/2)))
          {DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],0,matrix$tot_mut[3],as.numeric(sum(matrix$tot_mut)),cn_code)))
          colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
          }else
            {if((matrix$median[2]>(matrix$median[3]/2))&(matrix$median[2]<(matrix$median[3]*3/4)))
              {
              DF <- t(as.matrix(c(matrix$tot_mut[1],0,matrix$tot_mut[2],matrix$tot_mut[3],as.numeric(sum(matrix$tot_mut)),cn_code)))
              colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
              }
              else {stop("gain_3.point1") }
            }
            ## cyto[[i]] <- DF
        }else
        {if(nrow(matrix) == 4)
        {
           DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],matrix$tot_mut[3],matrix$tot_mut[4],as.numeric(sum(matrix$tot_mut)),cn_code)))
           colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
           ## cyto[[i]] <- DF
        }}}else
        {
          if(unique(first_def$class) == "LOH")
          {if(nrow(matrix) == 2)
          {
          DF <- t(as.matrix(c(as.numeric(as.character(matrix$tot_mut[1])),as.numeric(as.character(matrix$tot_mut[2])),0,0,sum(as.numeric(as.character(matrix$tot_mut))),cn_code)))
          colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
          ## cyto[[i]] <- DF
          }else
          {if(nrow(matrix) == 1) # fixme:
            stop("less than 2 clusters in LOH")
            }}
          }
        }}

    return(DF)
}

mol_time_warn <- function(warn,mess,arg) {
   if (warn) {
      print(sprintf(mess,arg))
      }
   return(TRUE)
}

mol_time_plot <- function (sample_code,all_code_final) {
  pdf(sprintf("%s_mol_time_estimate.pdf",sample_code), height=6, width=13)
  par(mar=c(10,7,3,3), xpd=F)
  par(mfrow=c(1,1))
  name_all2<-paste(all_code_final$chr, all_code_final$type, sep=" ")
  
  plot(c(1:nrow(all_code_final)),all_code_final$plot_dot1, ylim=c(0,2), xlim=c(0,(length(all_code_final$plot_dot1)+1)), xaxt="n", ylab="", 
       xlab="", pch=16, col=("dodgerblue"), 
       las=2, cex.axis=2, cex.lab = 1.5, cex=1.5, main=sample_code, cex.main=2.5)
  axis(1, at = c(1:nrow(all_code_final)),  col.axis = "black", labels=name_all2 , lwd.ticks =1, cex.axis = 2, las=2)
  mtext("Molecular Time", side = 2, line = 4, cex=3)
  
  par(new=TRUE)
  for(i in (1:nrow(all_code_final))) {
    segments(i,all_code_final$IC_down_cn1[i], i,all_code_final$IC_up_cn1[i],lwd=1, lty = 2, col="grey50")
    segments(i-0.2,all_code_final$IC_down_cn1[i],i+0.2,all_code_final$IC_down_cn1[i],lwd=1, lty = 1, col="grey50")
    segments(i-0.2,all_code_final$IC_up_cn1[i],i+0.2,all_code_final$IC_up_cn1[i],lwd=1, lty = 1, col="grey50")
  }
  
  second<-all_code_final
  second[is.na(second)]<-0
  par(new=TRUE)
  plot(c(1:nrow(all_code_final)) +0.1, all_code_final$plot_dot2, ylim=c(0,2), xlim=c(0,(length(all_code_final$plot_dot1)+1)), xaxt="n", ylab="", xlab="", pch=16, col=("red"), 
       las=2, cex.axis=2, cex.lab = 2, cex=1.5)
  
  num_gain_2<-which(all_code_final$plot_dot2!=0)
  
  for(w in (num_gain_2)) {
    segments(w+0.1,second$IC_down_cn2[w], w+0.1, second$IC_up_cn2[w], lwd=1, lty = 2, col="coral1")
    segments(w-0.1,second$IC_down_cn2[w],w+0.3,second$IC_down_cn2[w],lwd=1, lty = 1, col="coral1")
    segments(w-0.1,second$IC_up_cn2[w],w+0.3,second$IC_up_cn2[w],lwd=1, lty = 1, col="coral1")
    par(new=TRUE)
  }
  
  third<-all_code_final
  third[is.na(third)]<-0
  par(new=TRUE)
  plot(c(1:nrow(all_code_final)) +0.2, all_code_final$plot_dot3, ylim=c(0,2), xlim=c(0,(length(all_code_final$plot_dot3)+1)), xaxt="n", ylab="", xlab="", pch=16, col=("forestgreen"), 
       las=2, cex.axis=2, cex.lab = 2, cex=1.5)
  
  num_gain_2<-which(all_code_final$plot_dot3!=0)
  
  for(w in (num_gain_2)) {
    segments(w +0.2,third$IC_down_cn3[w], w+0.2, third$IC_up_cn3[w], lwd=1, lty = 2, col="forestgreen")
    segments(w-0.2+0.2,third$IC_down_cn3[w],w+0.2+0.2,third$IC_down_cn3[w],lwd=1, lty = 1, col="forestgreen")
    segments(w-0.2+0.2,third$IC_up_cn3[w],w+0.2+0.2,third$IC_up_cn3[w],lwd=1, lty = 1, col="forestgreen")
    par(new=TRUE)
  }
  
  dev.off()
  return(TRUE)
}
