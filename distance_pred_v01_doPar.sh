#!/bin/bash
#SBATCH --job-name DIST8
#SBATCH --qos normal
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -o sim_pwr.logba3
#SBATCH -e sim_pwr.err3
#SBATCH --time=24:00:00
#SBATCH --array=1-9,11-25,27-60,62-158,160-232




module load R
R --quiet --no-save > /dev/null << EOF

source("/home/fyu/me_impute/code_server/rbf_functions.R")
source("/home/fyu/me_impute/code_server/functions.R")
library(doParallel)


df_beta_3 = read.table("/home/fyu/me_impute/code_server/df_beta_3.txt", header = TRUE, sep="\t", stringsAsFactors = F)

IR.annot.mapinfo = read.table("/home/fyu/me_impute/insulinResistance/A-GEOD-13534_comments.mapinfo.txt", header = TRUE, sep="\t", stringsAsFactors = F)
map = IR.annot.mapinfo\$Comment.MAPINFO[which(IR.annot.mapinfo\$Comment.CHR.==18)]
map=sort(map)#length=5922

clust.name="clust.3000.7.ave.bydbp5"

match.enough = get.match.enough(df_beta_3,clust.name, map,lowerb=0.01, upperb=1, min.num=5)
which.in.450= get.which.in.450(df_beta_3, clust.name="clust.3000.7.ave.bydbp5",map = map)
cluster.WGBS=df_beta_3[,"clust.3000.7.ave.bydbp5"]
cluster.WGBS.unique=unique(cluster.WGBS)

out = read.table("/home/fyu/me_impute/code_server/flexmix_LSc_misc/out_crCrossCheck.txt", header = TRUE, sep="\t", stringsAsFactors = F)
colnames(out)[c(1,4,7)] <- c("cluster", "gamma.cor" ,"gamma.rmse")

df_beta_4 = df_beta_3[which(df_beta_3\$clust.3000.7.ave.bydbp5 %in% match.enough),]

DIST = 100

if(FALSE){
for(a in 1:25){


    f_O_1 <- paste("/home/fyu/me_impute/code_server/distance/cp3/bb_", DIST, "_", a, ".txt", sep="")


slct.unknown = as.vector(as.matrix(read.table("/home/fyu/me_impute/code_server/distance/slct.txt",sep="\t", header = FALSE)[sim,]))
df_beta_4$slct=0
df_beta_4$slct[slct.unknown]<-1



  aa = match.enough[a]
  df = df_beta_4[which(df_beta_4[,clust.name]==aa),]
  T = which(df$slct==1)

  library(randomForest)

no_cores <- 10
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores) 
  
  
  
  if(length(T)>0){
    
	foreach( t = T) %dopar%   {
	
      
      close = which(abs(df$start-df$start[t] )<DIST)
      
	   if(length(close)==0){
	    betas0 = df[,3:22]
		xx1 = df[,2]
	  }else{
	     betas0 = df[-close,3:22]
	     xx1 = df[-close,2]
	  }
	  	  
      betas1 = df[t,3:22]
      xxpred =  df[t,2]
      num.center=max(5, 0.1*length(xx1)) 
      xx3=seq(first(xx1),last(xx1), floor((last(xx1)-first(xx1)) /num.center ) )
      
      mx=median(xx1)
      X1=xx1-mx
      X3=xx3-mx
      Xpred = xxpred-mx
      gamma0=out[which(out$cluster==aa), "gamma.rmse"]
      
      
      betas.pred = matrix(data=NA, ncol=22, nrow=1)
      betas.pred[,1]=aa
      betas.pred[,2]=xxpred
      
	  if(FALSE){
	  for(j in 1:ncol(betas1)){
	  print(c(aa,t,j))
        phi01 = make.phi(X = X1,  gamma0, mus=X3)
        dataset=data.frame(Y = betas0[,j],phi01[,-1])
        
        phi.pred = make.phi(X = Xpred, gamma0, mus=X3)
        if(nrow(phi.pred)>1) {
          dataset.pred = data.frame(Y = betas1[,j], phi.pred[,-1])
        }else{     
          dataset.pred = data.frame(Y = betas1[,j],t(phi.pred[,-1]) )  
        }
        
        
        if(length(which(is.na(dataset$Y)))>0  ){
          dataset=dataset[-which(is.na(dataset$Y)),]
        }
        
        
        if(nrow(dataset)<ncol(dataset)){
          betas.pred[,j+2] = NA
        }else{
          
          robustM=function(dataset){
            tryCatch({
              flexmix(Y ~ ., data = dataset, k = 2)
            },  error=function(e){
              cat("ERROR :",conditionMessage(e), "\n"); NULL
            })
          }
          
          
          m1 <- robustM(dataset=dataset)
          if(!is.null(m1)){
            yy1 <- predict(m1, newdata = dataset.pred)
            cc <- clusters(m1, newdata = dataset.pred)
            betas.pred[,j+2] <- cluster.assigned.pred(yy1, cc)
          }
          
          
        }
      }
      
      
	  }
      
      write.table(betas.pred, f_O_1, col.names = F, row.names = F, sep="\t", append = TRUE)
      
  
   }
   
stopCluster(cl)
	
	}
	
	
	
	
	
    

}


}






for(sim in 8:10){

    f_O_1 <- paste("/home/fyu/me_impute/code_server/distance/cp8/bb_", DIST, "_", sim,"_", $SLURM_ARRAY_TASK_ID, ".txt", sep="")


slct.unknown = as.vector(as.matrix(read.table("/home/fyu/me_impute/code_server/distance/slct.txt",sep="\t", header = FALSE)[sim,]))
df_beta_4\$slct=0
df_beta_4\$slct[slct.unknown]<-1



  aa = match.enough[$SLURM_ARRAY_TASK_ID]
  df = df_beta_4[which(df_beta_4[,clust.name]==aa),]
  T = which(df\$slct==1)

  library(randomForest)

no_cores <- 10
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores) 
  
  
  
  if(length(T)>0){
    
	foreach( t = T) %dopar%   {
	
      
      close = which(abs(df\$start-df\$start[t] )<DIST)
      
	   if(length(close)==0){
	    betas0 = df[,3:22]
		xx1 = df[,2]
	  }else{
	     betas0 = df[-close,3:22]
	     xx1 = df[-close,2]
	  }
	  	  
      betas1 = df[t,3:22]
      xxpred =  df[t,2]
      num.center=max(5, 0.1*length(xx1)) 
      xx3=seq(first(xx1),last(xx1), floor((last(xx1)-first(xx1)) /num.center ) )
      
      mx=median(xx1)
      X1=xx1-mx
      X3=xx3-mx
      Xpred = xxpred-mx
      gamma0=out[which(out\$cluster==aa), "gamma.rmse"]
      
      
      betas.pred = matrix(data=NA, ncol=22, nrow=1)
      betas.pred[,1]=aa
      betas.pred[,2]=xxpred
      
	  if(TRUE){
	  for(j in 1:ncol(betas1)){
	  print(c(aa,t,j))
        phi01 = make.phi(X = X1,  gamma0, mus=X3)
        dataset=data.frame(Y = betas0[,j],phi01[,-1])
        
        phi.pred = make.phi(X = Xpred, gamma0, mus=X3)
        if(nrow(phi.pred)>1) {
          dataset.pred = data.frame(Y = betas1[,j], phi.pred[,-1])
        }else{     
          dataset.pred = data.frame(Y = betas1[,j],t(phi.pred[,-1]) )  
        }
        
        
        if(length(which(is.na(dataset\$Y)))>0  ){
          dataset=dataset[-which(is.na(dataset\$Y)),]
        }
        
        
        if(nrow(dataset)<ncol(dataset)){
          betas.pred[,j+2] = NA
        }else{
          
          robustM=function(dataset){
            tryCatch({
              flexmix(Y ~ ., data = dataset, k = 2)
            },  error=function(e){
              cat("ERROR :",conditionMessage(e), "\n"); NULL
            })
          }
          
          
          m1 <- robustM(dataset=dataset)
          if(!is.null(m1)){
            yy1 <- predict(m1, newdata = dataset.pred)
            cc <- clusters(m1, newdata = dataset.pred)
            betas.pred[,j+2] <- cluster.assigned.pred(yy1, cc)
          }
          
          
        }
      }
      
      
	  }
      
      write.table(betas.pred, f_O_1, col.names = F, row.names = F, sep="\t", append = TRUE)
      
  
   }
   
stopCluster(cl)
	
	}
	
	
	
	
	
    
}





EOF
