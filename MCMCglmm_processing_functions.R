#30/06/20
#Mads F. Schou
#Functions for data processing of MCMCglmm

#1. ReportRandomVarianceMCMC: Summarizes variance from random effects
#2. ReportFixedMCMC: Summarizes effect sizes and signifcance of fixed effects
#3. ReportCorrelationsMCMC: Estimate correlations using covariances and variances
#4. Settings for workbook
#5. Custom contrasts fixed
#6. Custom contrasts random
#7. R2 marginal (not finished!)

#######################################
###--- 1. ReportRandomVarianceMCMC
#######################################


ReportRandomVarianceMCMC = function(x, roundto = 3){
  
  library(reshape);library(MCMCglmm)
  
  MyRan <- data.frame(matrix(unlist(lapply(colnames(x$VCV), FUN = function(y){cbind(posterior.mode(x$VCV[,y]), HPDinterval(x$VCV[,y]))})), ncol = 3, byrow = T))
  MyRan$Effect <- colnames(x$VCV)
  MyRan[,c(1,2,3)] <- round(MyRan[,c(1,2,3)],10) #rounding as sometimes covars are not identical 
  MyRan <- MyRan[!duplicated(MyRan[,c(1,2,3)]) & !duplicated(MyRan[,c(1,2,3)], fromLast = TRUE),] #Remove covariances by deleting rows appearing twice
  MyRan[,c("X1","X2","X3")] <- round(MyRan[,c("X1","X2","X3")],roundto)
  
  #Split by levels and adjust variance names
  if(length(unlist(strsplit(MyRan$Effect, split = "\\."))) == nrow(MyRan)){
    MyRan[c("Random Effects: Variances","Level")] <- MyRan$Effect #If not levels (no slopes), just repeat names across columns
  }else{
    MySplit <- colsplit(MyRan$Effect, split = "\\.", names = c("RandomEffect","Level"))
    MyRan$Level <- as.character(MySplit[, ncol(MySplit)])
    MyRan$RandomEffect <- apply(MyRan,1, FUN = function(x) gsub(paste(".",x[5],sep = ""),"", x[4])) #Remove level part
    #Adjust variance names
    
    for(i in 1:nrow(MyRan)){
      
      MySplit <- colsplit(MyRan$RandomEffect[i], split = "\\:", names = c("out"))
      if(ncol(MySplit) == 1){
        MyRan[i,c("Random Effects: Variances")] <- as.character(MyRan$RandomEffect[i])
        }
      if(ncol(MySplit) == 2){
        MyRan[i,c("Random Effects: Variances")] <- as.character(colsplit(MyRan$RandomEffect[i], split = "\\:", names = c("out","Random Effects: Variances"))[,2])
        }
      if(ncol(MySplit) == 4){
        MyRan[i,c("Random Effects: Variances")] <- as.character(apply(colsplit(MyRan$RandomEffect[i], split = "\\:", names = c("out","Random Effects: Variances"))[,c(1,2)], 1, function(x) paste(x, collapse = ":")))
      }
      if(ncol(MySplit) == 6){
        MyRan[i,c("Random Effects: Variances")] <- as.character(apply(colsplit(MyRan$RandomEffect[i], split = "\\:", names = c("out","Random Effects: Variances"))[,c(1,2,3)], 1, function(x) paste(x, collapse = ":")))
      }
    }
  }
  
  MyRan[,c("Random Effects: Variances")] <-  gsub("units","residuals",MyRan[,c("Random Effects: Variances")])
  MyRan[,c("Level")] <-  gsub("units","residuals",MyRan[,c("Level")])
  
  MyRan["Posterior Mode (CI)"] <- paste(MyRan$X1, " (", MyRan$X2, ",", MyRan$X3, ")", sep = "")
  return(MyRan[,c("Random Effects: Variances","Posterior Mode (CI)","Level")])
}


#Usage: ranmodeOut <- ReportRandomVarianceMCMC(mcmc2)

#######################################
###--- 2. ReportFixedMCMC
#######################################

ReportFixedMCMC = function(x, remove = NULL){
  MySol <- data.frame(summary(x)$solutions)
  colnames(MySol) <- c("PostMean","CIlwr","CIupr","Sampling","pMCMC")
  MySol["Fixed Effects"] <- row.names(MySol)
  MySol <- MySol[!MySol$`Fixed Effects` %in% remove,] #Remove unwanted effects
  row.names(MySol) <- NULL
  MySol[,c("Fixed Effects")] <-  gsub("\\(Intercept\\)","Intercept",MySol[,c("Fixed Effects")])
  
  MySol[,c("PostMean","CIlwr","CIupr")] <- round(MySol[,c("PostMean","CIlwr","CIupr")],2)
  MySol[,c("pMCMC")] <- round(MySol[,c("pMCMC")],3)
  MySol[,c("Sampling")] <- round(MySol[,c("Sampling")])
  MySol["Posterior Mean (CI)"] <- paste(MySol$PostMean, " (", MySol$CIlwr, ",", MySol$CIupr, ")", sep = "")
  return(MySol[, c("Fixed Effects","Posterior Mean (CI)","pMCMC")])
}

#Usage: fixefecOut <- ReportFixedMCMC(mcmc2)

#######################################
###--- 3. ReportCorrelationsMCMC
#######################################

ReportCorrelationsMCMC = function(x, roundto = 2){
  #Get vars
  Vars <- ReportRandomVarianceMCMC(x)[c("Random Effects: Variances","Level")]
  colnames(Vars)[1] <- "Var"
  #"units" was automatically renamed to "residuals"; we reverse it
  Vars$Level[Vars$Level %in% "residuals"] <- "units"
  #Vars <- Vars[Vars$Level != "units",]
  
  #Get covars within each level
  Covars <- as.numeric()
  for(mylevel in levels(factor(Vars$Level))){
    tempCovars <- expand.grid(Vars[Vars$Level %in% mylevel ,c(1,1)])
    tempCovars$Level <- mylevel
    Covars <- rbind(Covars,tempCovars)
  }
  Covars <- Covars[Covars$Var != Covars$Var.1,]
  Covars <- Covars[!duplicated(Covars),]
  Covars$CovarNames <- paste(paste(Covars$Var,Covars$Var.1,sep = ":"),Covars$Level, sep = ".") #Gives two of each covar, just as in summary$VCV
  Covars$VarNames1 <- paste(paste(Covars$Var,Covars$Var,sep = ":"),Covars$Level, sep = ".")
  Covars$VarNames2 <- paste(paste(Covars$Var.1,Covars$Var.1,sep = ":"),Covars$Level, sep = ".")
  #Estimate correlations
  corrs <- NULL
  corrs=matrix(nrow = nrow(x$VCV), ncol = 0)
  for(i in 1:nrow(Covars)){
    covar.post <- x$VCV[,colnames(x$VCV) %in% Covars$CovarNames[i]]
    var1.post <- x$VCV[,colnames(x$VCV) %in% Covars$VarNames1[i]]
    var2.post <- x$VCV[,colnames(x$VCV) %in% Covars$VarNames2[i]]
    tmpcor<-covar.post/sqrt(var1.post*var2.post)
    
    if(!all(is.na(tmpcor[]))){ #If not empty: covar estimated by model
      corrs<-cbind(corrs,tmpcor)
      colnames(corrs)[ncol(corrs)] <- Covars$CovarNames[i]
    }
  }
  corrs<-as.mcmc(corrs)
  
  #Cor Summaries
  cor1=paste(round(posterior.mode(corrs),roundto)," (",round(HPDinterval(corrs)[,1],roundto), ",",round(HPDinterval(corrs)[,2],roundto),")",sep="")
  ncors<-ifelse(is.null(dim(corrs)), 1,dim(corrs)[2])
  nits<-ifelse(is.null(dim(corrs)),length(corrs), dim(corrs)[1])
  if(ncors >1){
    pCor=pmax(0.5/nits, pmin(colSums(corrs[,1:ncors, drop = FALSE] > 0)/nits, 1 - colSums(corrs[, 1:ncors, drop = FALSE] > 0)/nits))*2
  } else  {
    pCor=pmax(0.5/nits, pmin(sum(corrs[,drop = FALSE] > 0)/nits, 1 - sum(corrs[, drop = FALSE] > 0)/nits))*2
  }
  randomCorr<-data.frame("Random Effects: Correlations"=colnames(corrs),"Posterior Mode (CI)"=cor1,"pMCMC"=round(pCor,3), check.names=FALSE)
  #Remove duplicates (each corr will apear twice as they do so in the VCV)
  randomCorr <- randomCorr[!duplicated(randomCorr$`Posterior Mode (CI)`),]
  randomCorr[,c("Random Effects: Correlations")] <-  gsub("units","residuals",randomCorr[,c("Random Effects: Correlations")])
  return(randomCorr)
}

#Usage: CorrOut <- ReportCorrelationsMCMC(mcmc2)

#######################################
###--- 4. Settings for workbook
#######################################

OpenWorkBookMCMCglmm <- function(title, sheetname){
  #Table headers
  hs1=createStyle(fgFill = "grey80", halign = "LEFT", textDecoration = "bold",border = "TopBottom")
  assign("hs1", hs1, envir = .GlobalEnv) #add to env outside function
  hs2=createStyle(halign = "LEFT",border = "TopBottom",textDecoration = "bold")
  assign("hs2", hs2, envir = .GlobalEnv) #add to env outside function
  #table title
  header=data.frame(col1=c(""),col2=c(""),col3=c(""))
  colnames(header)<-c(title,"","")
  assign("workbook", createWorkbook(), envir = .GlobalEnv) #add to env outside function
  addWorksheet(workbook, sheetname)
  #Bold pMCMC values less than 0.05
  conditionalFormatting(workbook, sheetname, cols=3, rows=1:10000, rule="<0.05", style = createStyle(textDecoration="bold"))
  #Write header into notebook
  writeData(workbook, sheetname, header, startCol = 1, startRow = 1,headerStyle = hs1)
}


#######################################
###--- 5. Custom contrasts fixed
#######################################

FixedEffectContrasts = function(model,Eff1,Eff2){
  Diff <- model$Sol[,Eff1]-model$Sol[,Eff2]
  Name <- paste(Eff1, " vs ",Eff2, sep = " ")
  pDiff=pmax(0.5/length(Diff), pmin(sum(Diff[,drop = FALSE] > 0)/length(Diff), 1 - sum(Diff[, drop = FALSE] > 0)/length(Diff)))*2
  ModeCI=paste(round(posterior.mode(Diff),2)," (",round(HPDinterval(Diff)[,1],2), ",",round(HPDinterval(Diff)[,2],2),")",sep="")
  out <- data.frame('Effects'=Name,"Estimate"=ModeCI, "pMCMC"=round(pDiff,3))
  colnames(out) = c("Fixed Effect Contrasts","Posterior Mode (CI)","pMCMC")
  return(out)
}


#######################################
###--- 6. Custom contrasts random
#######################################

RandomEffectContrasts = function(model,Eff1,Eff2){
  Diff <- model$VCV[,Eff1]-model$VCV[,Eff2]
  Name <- paste(Eff1, " vs ",Eff2, sep = " ")
  pDiff=pmax(0.5/length(Diff), pmin(sum(Diff[,drop = FALSE] > 0)/length(Diff), 1 - sum(Diff[, drop = FALSE] > 0)/length(Diff)))*2
  ModeCI=paste(round(posterior.mode(Diff),3)," (",round(HPDinterval(Diff)[,1],3), ",",round(HPDinterval(Diff)[,2],3),")",sep="")
  out <- data.frame('Effects'=Name,"Estimate"=ModeCI, "pMCMC"=round(pDiff,3))
  colnames(out) = c("Variance Contrasts","Posterior Mode (CI)","pMCMC")
  return(out)
}

#######################################
###--- 7. R2 marginal (not finished!)
#######################################

#This is the raw version I got from Charlie
#I need to create the distvar variable, which is probably something like this:

# link_var[link_var=="gaussian"]<-0
# link_var[link_var=="poisson"]<-0
# link_var[link_var=="logit"]<-pi^2/3
# link_var[link_var=="probit"]<-1

#R2 Marginal function: Code from Sanchez-Tojar et al 2018 eLife modified into a function for calculation R2marginal re Nakagawa and Schielzeth 2013. distVar: logit = pi^2/3, log = log(1/exp(intercept)+1)
R2marginal<-function(model,distVar=0) {
  model$VCV<-model$VCV[, colnames(model$VCV) !="sqrt(mev):sqrt(mev).meta"]#remove inverse weights
  vmVarF<-numeric(dim(model$Sol)[1])
  for(i in 1:dim(model$Sol)[1]){
    Var<-var(as.vector(model$Sol[i,c(1:model$Fixed$nfl)] %*% t(model$X)))
    vmVarF[i]<-Var}
  R2m<-as.mcmc(100*(vmVarF/(vmVarF+rowSums(model$VCV)+distVar)))
  R2m<-paste(round(mean(R2m),2)," (",round(HPDinterval(R2m)[,1],2), " , ",round(HPDinterval(R2m)[,2],2),")",sep="")
  return(R2m)
}

