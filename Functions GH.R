# Survival in group 3 pulmonary hypertension
#
# Code relates to the manuscript: Dawes et al
# "The use of type 5 phosphodiesterase inhibitors is associated with improved survival in patients with group 3 pulmonary hypertension"
#




      # Load packages
      options(warn=-1)
      
      suppressPackageStartupMessages({
        library(Amelia)
        library(blme)
        library(BSDA)
        library(DescTools)
        library(epitools)
        library(lubridate)
        library(mice)
        library(miceadds)
        library(ordinal)
        library(psfmi)
        library(reshape2)
        library(utils)
        library(readxl)
        library(stats)
        library(survival)
        library(survminer)
        library(grDevices)
        library(RColorBrewer)
        library(graphics)
        library(ggplot2)
        library(stats)
        
        #library(mitools)
        #library(mix)
        #library(glmnet)
        #library(Hmisc)
        #library(GGally)
        #library(car)
        #library(nlme)
        #library(splines)
        #library(spBayesSurv)
        #library(coda)
        #library(tidyverse)
        #library(flexsurv)
        
        
        
      })


# Define Functions



#categorise<- function(cat.vec, a)
#          {
#            a.categorised<- rep(0, length(a))
#            a.categorised[which(is.na(a)==TRUE)]<- NA
#            for (i in 1:length(cat.vec)) {
#              x<- which(a<=cat.vec[i])
#              y<- which(a.categorised==0)
#              a.categorised[intersect(x,y)]<- i}
#            return(a.categorised)
#          }




pooling.by.Rubins.rules<- function(matrix.mean, matrix.sd, imputations)
{

  Vw.mat<- matrix.sd^2
  Vb.mat<- matrix.mean
  Vw<- apply(Vw.mat,2,mean)
  Vb<- apply(Vb.mat,2,var)
  
    # Degrees of freedom correction
    DoFc<- 1 + (1/imputations)
    
    # Final variance and SD
    V = Vw + Vb * DoFc
    SD = sqrt(V)
    
    # Final coefficient means and SD
    coef.df<- data.frame(mean=apply(Vb.mat,2,mean), SD=SD)
    coef.df$LowerCI<- coef.df$mean - 1.96*coef.df$SD
    coef.df$UpperCI<- coef.df$mean + 1.96*coef.df$SD
    coef.df$PerCentChance<- 100 * pnorm(abs(coef.df[,1]) / coef.df[,2])
    
    return(coef.df)
  }



median.cat<- function(x)
  {
  n<- nrow(x)
  if (n==2) {x.sorted<- x[,order(x[1,], x[2,])]}
  if (n==3) {x.sorted<- x[,order(x[1,], x[2,], x[3,])]}
  if (n==4) {x.sorted<- x[,order(x[1,], x[2,], x[3,], x[4,])]}
  
  p<- ncol(x)
  return(x.sorted[,round(p,0)])
  }
  


# t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
# {
#   if( equal.variance==FALSE ) 
#   {
#     se <- sqrt( (s1^2/n1) + (s2^2/n2) )
#     # welch-satterthwaite df
#     df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
#   } else
#   {
#     # pooled standard deviation, scaled by the sample sizes
#     se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
#     df <- n1+n2-2
#   }      
#   t <- (m1-m2-m0)/se 
#   dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
#   names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
#   return(dat) 
# }


# make.final.fmla<- function(fmla.final.fxd.effects, data, selectedVars, fxd.effects.factors)
#           {
#             initial.y<- fmla.final.fxd.effects
#             add<- NULL
#             
#             # Fixed effects
#                   for (i in 1:length(colnames(data)[selectedVars])) {
#                     if ((colnames(data)[selectedVars][i] %in% fxd.effects.factors) == TRUE) {
#                       add<- paste("as.factor(",colnames(data)[selectedVars][i],")",sep="")} else {
#                         add<- paste(colnames(data)[selectedVars][i],sep="")}
#                     if (fmla.final.fxd.effects!= initial.y) {add<- paste(" + ",add)}
#                     fmla.final.fxd.effects<- paste(fmla.final.fxd.effects, add, sep="")}
#             
#                   fmla.final.fxd.effects<- as.formula(fmla.final.fxd.effects)
#             
#             # Random effects
#                   fmla.final.rdm.effects<- "random = ~ tstart | IDn + "
#                   initial.y<- fmla.final.rdm.effects
#                   
#                   for (i in 1:length(colnames(data)[selectedVars])) {
#                     if ((colnames(data)[selectedVars][i] %in% fxd.effects.factors) == TRUE) {
#                       add<- paste("as.factor(",colnames(data)[selectedVars][i],") | IDn",sep="")} else {
#                         add<- paste(colnames(data)[selectedVars][i]," | IDn",sep="")}
#                     if (fmla.final.rdm.effects!= initial.y) {add<- paste(" + ",add)}
#                     fmla.final.rdm.effects<- paste(fmla.final.rdm.effects, add, sep="")}
#                   
#                   fmla.final.rdm.effects<- as.formula(fmla.final.rdm.effects)
#             
#             # Print results
#                   
#                   cat("\n \n Final formula \n")
#                   
#                   message<- paste(colnames(data)[selectedVars], collapse = " + ")
#                   line<- paste(paste(rep("-", nchar(message)+1),collapse=""),"\n")
#                   
#                   cat(line)
#                   cat("\n", message," \n")
#                   cat(line)
#             
#             # Update log of which covariates are important
#                   
#             return(list(FE = fmla.final.fxd.effects, RE = fmla.final.rdm.effects))
#           }
# 



make.contiguous.dates<- function(e, drug.start, drug.stop, drug.cols, dose.equivalents, drug.class.for.each.drug)
          {
            # Identify the diagnosis and censoring dates as 'end buffers' to the drug data
            diagnosis.date<- ymd(t(e$dateDIAGNOSIS))
            censoring.date<- ymd(t(e$Date.censoring))
            
            start.dates<- ymd(t(e[drug.start + (drug.cols - 1)]))
            stop.dates<- ymd(t(e[drug.stop + (drug.cols - 1)]))
            if (any(difftime(stop.dates, start.dates) < 0)) {cat("\n Drug stop dates for",d2[i,]$URN,"are BEFORE the given start dates")}
            
            o<- order(start.dates)
            start.dates<- start.dates[o]
            stop.dates<- stop.dates[o]
            
            # If the drug history doesn't START at the diagnosis date, then: 
            if (diagnosis.date < head(start.dates,1)) {
              # ...add the period between diagnosis and first drug date onto the start and stop dates
                    stop.dates<- c(start.dates[1], stop.dates)
                    start.dates<- c(diagnosis.date, start.dates)
              # ...add a dose of 'zero' to the beginning of the drug history
                    dose.equivalents<- c(0, dose.equivalents)
              # ... and add a drug class of 'Nil' to the list of drug classes
                    drug.class.for.each.drug<- c("Nil",drug.class.for.each.drug)}
            
            # If the drug history doesn't END at the censoring date, then: 
            if (censoring.date > tail(stop.dates,1)) {
              # ...add the period between last drug date and censoring on to the start and stop dates
                    start.dates<- c(start.dates, tail(stop.dates,1))
                    stop.dates<- c(stop.dates, censoring.date)
              # ...add a dose of 'zero' to the end of the drug history
                    dose.equivalents<- c(dose.equivalents, 0)
              # ... and add a drug class of 'Nil' to the end of the list of drug classes
                    drug.class.for.each.drug<- c(drug.class.for.each.drug,"Nil")}
            
            return(list(start.dates = start.dates, stop.dates = stop.dates, dose.equivalents = dose.equivalents,
                        diagnosis.date = diagnosis.date, censoring.date = censoring.date, DateOfBirth = e$DoB, drug.class.for.each.drug = drug.class.for.each.drug))
          }


# MeanAndSD<- function(column, Ix.at.diagnosis, imputations)
#           {
#             w.col<- which(colnames(Ix.at.diagnosis)==column)
#             
#             hist(Ix.at.diagnosis[,w.col], main=column)
#             
#             imputed.mean<- sapply(1:imputations, function(x) {mean(Ix.at.diagnosis[which(Ix.at.diagnosis$Imp==x),w.col])})
#             imputed.summary<- sapply(1:imputations, function(x) {summary(Ix.at.diagnosis[which(Ix.at.diagnosis$Imp==x),w.col])})
#             imputed.sd<- sapply(1:imputations, function(x) {sd(Ix.at.diagnosis[which(Ix.at.diagnosis$Imp==x),w.col])})
#             
#             Mean = mean(imputed.mean)  
#             StdDev = mean(imputed.sd)
#             Summary = apply(imputed.summary,1,mean)
#             RawData = Ix.at.diagnosis[,w.col]
#             
#             return(list(Mean = Mean, SD = StdDev, Summary = Summary, RawData))
#           }
# 



drug.equivalents<- function(a,b,c,d) {return(round(100*b*1000 / d$`DOSE/DAY_MCG`[match(paste(a,c,sep=""), d$ABBREVIATION)],0))}


predictFEV1<- function(v)
            {
              if (v[1]=="0") {FEV1<- 0.0395*v[2]*100 - 0.025*v[3] - 2.6}
              if (v[1]=="1") {FEV1<- 0.043*v[2]*100 - 0.029*v[3] - 2.49}
              return(unlist(FEV1))
            }


predictFVC<- function(v)
            {
              if (v[1]=="0") {FVC<- 0.0443*v[2]*100 - 0.026*v[3] - 2.89}
              if (v[1]=="1") {FVC<- 0.0576*v[2]*100 - 0.026*v[3] - 4.34}
              return(unlist(FVC))
             }


eGFR<- function(creat, age, gender, black, method)
            {
              g<- b<- 1
              age<- as.numeric(age)
              
              if (method=="MDRD") { # Abbreviated MDRD equation
                    if (gender=="0") {g<- 0.742} 
                    if (black==TRUE) {b<- 1.210}
                    eGFR = 186 * (creat/88.4)^(-1.154) * (age^(-0.203)) * g * b}
                    
              if (method=="CKD-EPI") { #CKD-EPI Creatinine Equation
                    creat<- creat / 88.39998
                    K<- 0.9
                    alpha<- -0.411
                    
                    if (gender=="0") {g<- 0.7; alpha<- -0.329}
                    if (black==TRUE) {b<- 1.159}
                    
                    part.one<- sapply(creat/K, function(x) {min(x,1)})^alpha
                    part.two<- sapply(creat/K, function(x) {max(x,1)})
                    
                    eGFR = 141 * part.one * part.two * 0.993^age * b}
                    
              return(eGFR)
              
            }



LOCFandB<- function(y)
            {
              if (("ID" %in% colnames(y))==FALSE) {cat("\n No clear ID column - please rename one column name as 'ID'")}
              unique.IDs<- unique(y$ID)
              fields<- match(setdiff(colnames(y), "ID"), colnames(y))
              cat("\n Last observation carried forward and back...")
              for (i in 1:length(unique.IDs)){
                  for (j in fields){
                    rows<- which(y$ID==unique.IDs[i])
                    locf.data<- unlist(na.locf(na.locf(y[rows,j], na.rm=FALSE), fromLast=TRUE, na.rm=FALSE))
                    if (length(locf.data)>0) {y[rows,j]<- locf.data}
                    }
              }
              cat(".. completed")
              return (y)
            }


# LOCFandBplusMICE<- function(y, m, maxit, method)
#           {
#             y<- LOCFandB(y)  
#             cat("\n Multiple imputation with chained equations in progress...")
#             
#             colnames(y)<- make.names(colnames(y))
#             # Identify the columns which are NOT 'ID'
#                 numeric.column.names<- match(setdiff(colnames(y),"ID"), colnames(y))
#             # Convert all non-ID columns to numeric data
#                 for (i in numeric.column.names) {y[,i]<- as.numeric(unlist(y[,i]))}
#             # Impute
#                 test.imp<- mice(y[,numeric.column.names], m=m, maxit=maxit, method=method, print=F)
#             # Long version, convert all 1/2 columns to 0/1, and stick back to the IDs
#                 results<- complete(test.imp, 'long')
#                 results<- zero.and.one.cols(results,3)
#                 #results<- cbind(ID = rep(y$ID,m),results)
#                 
#             # Clear up any remaining zeros or NAs    
#                 for (i in numeric.column.names)
#                 {
#                 zeros.or.NAs<- union(which(y[,i]==0), which(is.na(y[,i])==TRUE))
#                 not.zeros.or.NAs<- setdiff(1:nrow(y), zeros.or.NAs)
#                 y[zeros.or.NAs,i]<- mean(y[not.zeros.or.NAs,i])
#                 }
#                 
#             cat("..completed")
#             return(results)
#           }
# 
# 
# 
# JustMICE<- function(y, m, maxit, method)
#           {
#                 cat("\n Multiple imputation with chained equations in progress...")
#                 colnames(y)<- make.names(colnames(y))
#                 numeric.cols<- na.omit(match(c("Outcome","tstop","tstart","futime","time","imp","HR","SBP","FAC","TAPSE","S","TRvel","RVSP","RA.pressure",
#                                                "RA.areaI","RA.volumeI","PAAT","LVEF.Simpson","Edecel","E.A","EE.","LA.areaI","LA.volumeI",
#                                                "FEV1","FEV1perc","FVC","FVCperc","TLcoc","TLcoperc","Kco","Kcoperc","mPAP","PAP","CO","PCWP"), colnames(y)))
#                 
#                 factor.cols<- na.omit(match(c("Gender","SR","RV.dilat","Long.dysfunc","Rad.dysfunc","RA.dilat","Effusion","LA.dilat"), colnames(y)))
#                 
#             # Impute
#                 test.imp<- mice(y, m=m, maxit=maxit, method=method, print=F)
#             # Long version, convert all 1/2 columns to 0/1, and stick back to the IDs
#                 results<- complete(test.imp, 'long')
#                 results<- zero.and.one.cols(results,3)
#             # Clear up any remaining zeros or NAs    
#                 
#                 
#                 for (i in numeric.cols){
#                   zeros.or.NAs<- union(which(y[,i]==0), which(is.na(y[,i])==TRUE))
#                   not.zeros.or.NAs<- setdiff(1:nrow(y), zeros.or.NAs)
#                   y[zeros.or.NAs,i]<- mean(y[not.zeros.or.NAs,i])}
#             
#                 cat("..completed")
#                 return(results)
#           }
# 

JustAMELIA<- function(y, m, maxit, max.emburn, logs)
          {
            cat("\n Amelia: Multiple imputation with chained equations in progress...")
            
            # Find the bounds of the current data and use this as a bound for imputation
                  bds<- apply(y, 2, range, na.rm=TRUE)
                  bds<- cbind(1:ncol(y), matrix(as.numeric(bds),ncol=2, byrow=T))
                  rownames(bds)<- colnames(y)
                  bds[match(logs, rownames(bds)),2:3]<- log(bds[match(logs, rownames(bds)),2:3])
                  
            # Make the IDs numeric
                  y$IDn<- match(y$ID, unique(y$ID))
            
            # Impute
                  test.imp<- amelia(y, m=m, maxit=maxit, empri = 0.1, emburn = c(1,max.emburn), idvars=c("ID","id","Diagnosis","Diagnosis_sub","AgeAtDiag"), #parallel = "multicore", ncpus = 4,  
                                     p2s=2, polytime = 1, intercs = FALSE, ts = "tstart", cs = "IDn", bounds = bds,
                                     noms = c("Gender","SR","RVdys","RV.dilat","Long.dysfunc","Rad.dysfunc","RA.dilat","Effusion","LA.dilat","Treated","TreatedPDE5"),
                                     logs = logs,
                                     ords=c("FC","DoseFctrPDE5","EmPHasis10"))
                              
            
            cat("..completed")
            
            return(test.imp)
          }


# draw.KM<- function(KMmodel, pval) 
#           {
#             p<- ggsurvplot(KMmodel, pval=pval, pval.method = TRUE, surv.median.line = "hv", conf.int = TRUE, risk.table = TRUE, ggtheme = theme_bw())
#             
#             surv_median <- as.vector(summary(KMmodel)$table[, "median"])
#             df <- data.frame(x1 = surv_median, x2 = surv_median, y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
#             
#             #p$plot <- p$plot + 
#             #  geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
#             #               linetype = "dashed", size = 0.5)+ # horizontal segment
#             #  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
#             #               linetype = "dashed", size = 0.5) # vertical segments
#             
#             p$plot<- p$plot + ggplot2::annotate("text", x = surv_median[1], y = 0.01, # x and y coordinates of the text
#                               label = paste(round(surv_median[1],2),"years"), size = 4)
#             
#             p$plot<- p$plot + ggplot2::annotate("text", x = surv_median[2], y = 0.01, # x and y coordinates of the text
#                               label = paste(round(surv_median[2],2),"years"), size = 4)
#             
#             
#             print(p)
#           }


#zero.and.one.cols<- function(x,start) {for (i in start:ncol(x)) {if (identical(as.numeric(as.character(names(table(x[,i])))),c(1,2))==TRUE) {x[,i]<- as.numeric(as.character(x[,i])) - 1}}; return(x)}


# Function to remove duplicate investigations
remove.duplicates<- function(test.data, test.dates, test.demographics, test.present)
          {
            date.and.ID<- paste(test.demographics$ID,"_",test.dates$DateOfScan,sep="")
            duplicate.Ix<- which(duplicated(date.and.ID)==TRUE)
            
            if (length(duplicate.Ix>0)) { test.data<- test.data[-duplicate.Ix,]
                                          test.dates<- test.dates[-duplicate.Ix,]
                                          test.demographics<- test.demographics[-duplicate.Ix,]
                                          test.present<- test.present[-duplicate.Ix,]
                                          cat("\n ",length(duplicate.Ix),"duplicates removed")
                                          cat("\n", date.and.ID[duplicate.Ix])
                                        } else {cat("\n No duplicate test results found")}
            
            return(list(test.data, test.dates, test.demographics, test.present))
          }




# LOCFandBplusMICEtest<- function(echo.data.minus, echo.data.i1, errors.cum, errors.mean, i, m, maxit, factors, na.locations, echo.data, nbr_missing)
#           {
#             # Imputation stage 1: Use LOCF and LOCB within a particular subject if there are pre-existing data points
#             
#             echo.data.i1<- LOCFandB(echo.data.minus, echo.data.i1)          
#             
#             # Imputation stage 2: MICE
#             colnames(echo.data.i1)<- make.names(colnames(echo.data.i1))
#             for (i in 2:ncol(echo.data.i1)) {echo.data.i1[,i]<- as.numeric(echo.data.i1[,i])}
#             echo.imp<- mice(echo.data.i1[,-1], m=m, maxit=maxit, method="pmm", seed=100, verbose=F)
#             echo.data.i2<- complete(echo.imp, 'long')
#             
#             # Aggregate numerical results
#             echo.data.i3<- cbind(echo.data.minus[,1], aggregate(echo.data.i2, by = list(echo.data.i2$.id),FUN= mean)[,-c(1:3)])
#             
#             # Compare imputed results with gold standard
#             #      errors.this.iteration<- as.numeric(echo.data[as.matrix(na.locations)]) - as.numeric(echo.data.i3[as.matrix(na.locations)])
#             #      for (i in 1:nbr_missing) {errors.cum[2, na.locations[i,2]] <- errors.cum[2,na.locations[i,2]] + abs(errors.this.iteration[i]);
#             #      errors.cum[1, na.locations[i,2]] <- errors.cum[1,na.locations[i,2]] + 1}
#             # Calculate mean error for this iteration
#             #      errors.mean.it<- 100*(errors.cum[2,-1] / errors.cum[1,-1]) / colMeans(data.matrix(echo.data[,-1]), na.rm=T)
#             #      errors.mean[maxit.counter, 3:(length(errors.mean.it)+2)]<- errors.mean.it 
#             # Compare imputed results with gold standard
#             results<- compare.imputed.with.gold(factors, na.locations,echo.data, echo.data.i3, nbr_missing, errors.cum, errors.mean)
#             return(results)
#           }
# 




change.data.type<- function(m)
          {
            numeric.cols<- na.omit(match(c("Height","Weight","HR","SBP","BMI","BSA","Outcome","tstop","tstart","futime","time","imp","HR","SBP","FAC","TAPSE","S","TRvel","RVSP","RA.pressure",
                                   "RA.areaI","RA.volumeI","PAAT","LVEF.Simpson","Edecel","E.A","EE.","LA.areaI","LA.volumeI",
                                   "FEV1","FEV1perc","FVC","FVCperc","TLcoc","TLcoperc","Kco","Kcoperc","mPAP","PAP","CO","PCWP"), colnames(m)))
            factor.cols<- na.omit(match(c("Diagnosis","Gender","imp","SR","RV.dilat","Long.dysfunc","Rad.dysfunc","RA.dilat","Effusion","LA.dilat"), colnames(m)))
            for (i in numeric.cols) {m[,i]<- as.numeric(as.character(m[,i]))}
            for (i in factor.cols) {m[,i]<- as.factor(unlist(m[,i]))}
            return(m)
          }



make.tdc.matrix<- function(data1, test.complete, test.column.name)
          {
            data2<- NULL
            column<- which(colnames(test.complete)==test.column.name)
            for (i in 1:nrow(data1))
            {
              scans.for.one.ID<- which(test.complete$ID %in% data1$ID[i])
              scan.dates.for.one.ID<- test.complete[scans.for.one.ID,column]
              
              futime<- as.Date(data1$DateOfCensoring[i]) - scan.dates.for.one.ID[which.min(scan.dates.for.one.ID)]
              times<- scan.dates.for.one.ID - scan.dates.for.one.ID[which.min(scan.dates.for.one.ID)]
              
              t<- length(times)
              o<- order(times)
              
              data2<- rbind(data2, cbind.data.frame(test.complete[scans.for.one.ID,], time=as.numeric(times) / 365, futime=rep(futime / 365,t))[o,])                                     
            }
            return (data2)
          }


# sampler<- function(t,nmax)
#           {
#             cat("\n Sampling..")
#             if (nmax==0) {n<- 1:nrow(t)} else
#             {
#               j=1; n<- NULL
#                   while (length(n)<nmax)
#                   {
#                     cat("\n ",j,"..",length(n))
#                     for (i in 1:185) {
#                       IDns<- which(t$IDn==i)
#                       n<- c(n, which(t$IDimp == unique(t[IDns,]$IDimp)[j]))}
#                       j<- j + 1
#                   }
#             }
#               cat("done.")
#               return(n)
#           }
#   

extract.class.columns.and.dates<- function(d2,class)
          {
            test.present<- NULL
            # a. Identify which columns the PFT data and the PFT dates are stored in
                test.cols<- unlist(sapply(paste(class,1:9,sep=""), function(x) {grep(x, colnames(d2))}))
                test.date.cols<- sort(unique(test.cols[grep(paste("date",class,sep=""), colnames(d2[,test.cols]))]))
                
            # b. Build the matrix storing the PFT date row and columns in    
                for (i in 1:length(test.date.cols)) {
                  r<- which(is.na(d2[,test.date.cols[i]])==FALSE)
                  #cat("\n number:",i," frequency:", length(r))
                  c<- rep(test.date.cols[i], length(r))
                  test.present<- rbind(test.present, cbind(r,c))
                }
                
                cat("\n Number of tests in this class available:", nrow(test.present))
            return(list(test.present, test.cols))
          }



extract.class.data.and.demographics<- function(test.present, test.cols, d2, ncols)
          {
            # Add progress bar
                pb <- txtProgressBar(min = 0, max = nrow(test.present), style = 3)
                
            # Make a dataframe to hold all the data pulled from the spreadsheet
                test.data<- data.frame(matrix(0, nrow=nrow(test.present), ncol=ncols))
                column.names.E1<- colnames(d2)[test.cols[1]+(2:(ncols+1))]
                dot.positions<- unlist(lapply(strsplit(column.names.E1, ''), function(x) max(which(x == '.'))))
                column.names<- substr(column.names.E1, 1, dot.positions-1)
                colnames(test.data)<- make.names(column.names)
                
            # Make dataframe to hold all the demographic data pulled from the spreadsheet
                test.demographics<- data.frame(ID = vector(length=nrow(test.present), mode="character"),
                                               Diagnosis = vector(length=nrow(test.present), mode="character"),
                                               
                                               Gender = vector(length=nrow(test.present), mode="integer"),
                                               AgeAtTest = vector(length=nrow(test.present), mode="integer"),
                                               Height = vector(length=nrow(test.present), mode="integer"),
                                               Outcome = rep(0, nrow(test.present)),
                                               stringsAsFactors = FALSE)
            # Make dataframe for dates
                test.dates<- data.frame(
                        DateOfScan = as.Date(rep("1900-01-01", length=nrow(test.present), format="%y/%m/%d")),
                        DateOfBirth = as.Date(rep("1900-01-01", length=nrow(test.present), format="%y/%m/%d")),
                        DateOfCensoring = as.Date(rep("1900-01-01", length=nrow(test.present), format="%y/%m/%d")))
                
                dateE.columns<- grep("date.Echo",colnames(d2))[-1]
                
                for (i in 1:nrow(test.present))
                        {
                        c<- (test.present[i,2]+2):(test.present[i,2]+ncols+1)
                        test.data[i,]<- as.numeric(unlist(c(d2[test.present[i,1],c])))
                        
                        test.dates$DateOfScan[i]<- ymd(unlist(c(d2[test.present[i,1], c[1]-2])[[1]]))
                        test.dates$DateOfBirth[i]<- ymd(d2[test.present[i,1],]$DoB)
                        test.dates$DateOfCensoring[i]<- ymd(as.Date(unlist(d2[test.present[i,1],]$Date.censoring), origin="1900-01-01"))
                        test.demographics$Gender[i]<- d2[test.present[i,1],"Gender"]
                        test.demographics$ID[i]<- unlist(d2[test.present[i,1],"URN"])
                        test.demographics$Diagnosis[i]<- d2[test.present[i,1],"Diagnosis_sub"]
                        test.demographics$Outcome[i]<- d2[test.present[i,1],"Dead"]
                        test.demographics$AgeAtTest[i]<- (test.dates$DateOfScan[i] - test.dates$DateOfBirth[i]) / 365
                        
                        nearest.date<- order(abs(ymd(as.Date(unlist(d2[test.present[i,1],dateE.columns]), origin="1900-01-01")) - ymd(as.Date(test.dates$DateOfScan[i], origin="1900-01-01"))))
                        test.demographics$Height[i]<- na.omit(unlist(d2[test.present[i,1], dateE.columns[nearest.date]+2]))[1]
                        test.demographics$Weight[i]<- na.omit(unlist(d2[test.present[i,1], dateE.columns[nearest.date]+3]))[1]
                        
                        setTxtProgressBar(pb, i)
                        }
                
  
              # f. Stick the test data to the test demographics and reformat
                    factors<- 1
                    for (i in 1:ncol(test.data)) {if (length(unique(test.data[,i]))<=6) {factors<- c(factors, i); test.data[,i]<- as.factor(unlist(test.data[,i]))}}
                    test.demographics$Gender<- as.numeric(test.demographics$Gender)
                    test.demographics$AgeAtTest<- (test.dates$DateOfScan - test.dates$DateOfBirth) / 365
                    
              # Check for missing height values and substitute in the population mean
                    missing.height<- which(is.na(test.demographics$Height)==TRUE)
                    missing.weight<- which(is.na(test.demographics$Weight)==TRUE)
                    
                    if (length(missing.height)>0) {
                        cat("\n Warning: ",length(missing.height)," subjects (",round(100*length(missing.height)/nrow(test.demographics),0),"%) had no height readings - left as NAs for imputation", sep="")}
                    #  test.demographics$Height[missing.height]<- mean(na.omit(test.demographics$Height))}
                    
              # Re-order the test data by (a) imputed ID and (b) chronological order of tests in preparation for LOCF/B
                          o<- order(test.demographics$ID, test.demographics$AgeAtTest)
                          test.data<- data.frame(test.data[o,])
                          colnames(test.data)<- make.names(column.names)
                          test.dates<- test.dates[o,]
                          test.demographics<- test.demographics[o,]
                    
                          
              return(list(test.data, test.dates, test.demographics))
}


# compare.demographics<- function(field, y, fig)
# {
#   col<- which(colnames(y)==field)
#   
#   All.mean<- round(mean(unlist(y[,col])),fig)
#   All.sd<- round(sd(unlist(y[,col])),fig)
#   
#   alive<- which(y$Outcome==0)
#   dead<- which(y$Outcome==1)
#   
#   Alive.mean<- round(mean(unlist(y[alive,col])),fig)
#   Alive.sd<- round(sd(unlist(y[alive,col])),fig)
#   Dead.mean<- round(mean(unlist(y[dead,col])),fig)
#   Dead.sd<- round(sd(unlist(y[dead,col])),fig)
#   
#   a.vs.d.ttest<- round(t.test(y[alive,col], y[dead,col])$p.value,4)
#   #prop.test(matrix(c(Alive.mean, 1-Alive.mean, Dead.mean, 1-Dead.mean), nrow=2))
#   
#   tstart.col<- which(colnames(y)=="tstart")
#   tstop.col<- which(colnames(y)=="tstop")
#   outcome.col<- which(colnames(y)=="Outcome")
#   
#   
#   fit<- lapply(1:imputations, function(x) {coxph(Surv(y[,tstart.col], y[,tstop.col], y[,outcome.col]==1) ~ y[,col], data=subset(Ix.at.diagnosis, Imp==x))})
#   cox.model<- summary(pool(fit))
#   
#   hr<- round(exp(cox.model[2]),2)
#   ci.lower<- round(exp(cox.model[2] - 1.96*cox.model[3]),2)
#   ci.upper<- round(exp(cox.model[2] + 1.96*cox.model[3]),2)
#   hr.pvalue<- round(cox.model[6],4)
#   
#   return(list(All.mean = All.mean, All.sd = All.sd, Alive.mean = Alive.mean, Alive.sd = Alive.sd, Dead.mean = Dead.mean,
#               Dead.sd = Dead.sd, p.value.ttest = a.vs.d.ttest, hr = hr, ci.lower = ci.lower, ci.upper = ci.upper,
#               hr.pvalue = hr.pvalue))
# }
# 


# interpret.CPH.interactions<- function(r)
#         {
#           # Identify the regression terms and the interaction terms
#           interactions<- grep(":",rownames(r))
#           non.interactions<- setdiff(1:nrow(r),interactions)
#           # Interpret the regression terms
#           answer<- paste("Higher",rownames(r)[non.interactions], c("decreases the risk of death","increases the risk of death")[(1.5 + sign(r[non.interactions,1]) / 2)])
#           cat(answer, sep="\n")
#           # Interpret the interaction terms
#           cat("\n -------------------------------")
#           if (length(interactions)>0) 
#             {
#             for (i in 1:length(interactions)){
#               f<- strsplit(rownames(r)[interactions],":")[[i]]
#               cat("\n ", c(paste("If",f[1],"is higher then the effect of",f[2],"is less"), paste("If",f[1],"is higher then the effect of",f[2],"is also greater"))[(1.5 + sign(r[interactions,1]) / 2)])}
#             } else {cat("\n No interaction terms found.")}
#         }
# 


# image.scale <- function(z, zlim, col = heat.colors(12),
#                         breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
#               if(!missing(breaks)){
#                 if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
#               }
#               if(missing(breaks) & !missing(zlim)){
#                 breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
#               }
#               if(missing(breaks) & missing(zlim)){
#                 zlim <- range(z, na.rm=TRUE)
#                 zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
#                 zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
#                 breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
#               }
#               poly <- vector(mode="list", length(col))
#               for(i in seq(poly)){
#                 poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
#               }
#               xaxt <- ifelse(horiz, "s", "n")
#               yaxt <- ifelse(horiz, "n", "s")
#               if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
#               if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
#               if(missing(xlim)) xlim=XLIM
#               if(missing(ylim)) ylim=YLIM
#               plot(1,1,type="n",ylim=ylim, xlim=xlim, axes = F, xlab="", ylab="")  
#               axis(1, at=seq(XLIM[1]+0.5, XLIM[2]-0.5,1), labels=seq(XLIM[1], XLIM[2]-1, 1))
#               for(i in seq(poly)){
#                 if(horiz){
#                   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
#                 }
#                 if(!horiz){
#                   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
#                 }
#               }
#             }



# Extra functions which aren't needed at the moment

#draw1<- function(errors.mean.it, errors.mean, maxit.counter,new,col,xoffset)
#          {
#            if (new==TRUE) {barplot(errors.mean.it, horiz=T, names.arg=names(errors.mean.it), col=col, cex.names=0.4, las=2)}
#            if (new==FALSE) {barplot(errors.mean.it, horiz=T, names.arg=names(errors.mean.it), col=col, cex.names=0.4, las=2, add=TRUE)}
#          }


#draw2<- function(errors.mean.it, errors.mean, maxit.counter,new,col,xoffset, yoffset)
#           {
#            if (new==TRUE) {plot(0,0,xlim=c(1,20), ylim=c(1,25), type='n', xlab="mtry",ylab="Ntree / Multivariate")}
#            if (maxit.counter==1) {labels<- mean(na.omit(errors.mean[1,2:ncol(errors.mean)]))} else {labels<- apply(errors.mean[1:maxit.counter,2:ncol(errors.mean)],1,mean, na.rm=T)}
#            text(errors.mean[1:maxit.counter,1]+xoffset, errors.mean[1:maxit.counter,2]+yoffset, labels=round(labels,1), cex=0.5)
#           }
# 
# 
# MICE<- function(echo.data.minus, echo.data.i1, errors.cum, errors.mean, m, maxit, method, factors, na.locations, echo.data, nbr_missing)
# {
#   colnames(echo.data.minus)<- make.names(colnames(echo.data.minus))
#   for (i in 2:ncol(echo.data.minus)) {echo.data.minus[,i]<- as.numeric(echo.data.minus[,i])}
#   if (method=="pmm") {echo.imp<- mice(echo.data.minus[,-1], m=m, maxit=maxit, method="pmm", seed=100, verbose=F)}
#   if (method=="rf") {echo.imp<- mice(echo.data.minus[,-1], ntree=m*10, method="rf", seed=100, verbose=F)}
#   
#   echo.data.i2<- complete(echo.imp, 'long')
#   
#   # Aggregate numerical results
#   echo.data.i3<- cbind(ID = echo.data.minus[,1], aggregate(echo.data.i2, by = list(echo.data.i2$.id),FUN= mean)[,-c(1:3)])
#   
#   # Compare imputed results with gold standard
#   results<- compare.imputed.with.gold(factors, na.locations,echo.data, echo.data.i3, nbr_missing, errors.cum, errors.mean)
#   
#   return(results)
# }
# 
# 
# HMISC<- function(echo.data.minus, echo.data.i1, errors.cum, nk)
# {
#   colnames(echo.data.minus)<- make.names(colnames(echo.data.minus))
#   for (i in 2:ncol(echo.data.minus)) {echo.data.minus[,i]<- as.numeric(echo.data.minus[,i])}
#   f<- as.formula(paste("~",paste(colnames(echo.data.minus)[-1], collapse=" + "),sep=""))
#   imputed<- aregImpute(f, echo.data.minus, n.impute=10, nk=0, type="pmm")
#   for (i in 1:(ncol(echo.data.minus)-1))
#   {
#     imputed.values<- apply(imputed$imputed[[i]],1,Mode)
#     for (j in 1:length(imputed.values)) {echo.data.i1[which(is.na(echo.data.minus[,i+1])),i+1] <- as.factor(imputed.values[j]-1)}
#   }
#   
#   # Compare imputed results with gold standard
#   results<- compare.imputed.with.gold(factors, na.locations, echo.data, echo.data.i1, nbr_missing, errors.cum, errors.mean)
#   errors.mean.it<- results[[1]]
#   errors.mean<- results[[2]]
#   errors.cum<- results[[3]]
#   return(list(errors.mean.it, errors.mean, errors.cum))
# }
# 
# 
# compare.imputed.with.gold<- function(factors, na.locations,echo.data, echo.data.i3, nbr_missing, errors.cum, errors.mean)
# {
#   na.factors<- (na.locations[,2])[na.locations[,2] %in% factors[factors>=3]] # Vector of rows in na.locations which have put NAs in factor fields
#   # Calculate errors.this.iteration for all fields where NAs have been inserted (continuous and factors)
#   errors.this.iteration<- abs(as.numeric(echo.data[as.matrix(na.locations)]) - as.numeric(echo.data.i3[as.matrix(na.locations)]))
#   # Calculate the errors.this.iteration for the fields with factors in
#   errors.this.iteration[na.factors]<- abs(as.numeric(echo.data[as.matrix(na.locations[na.factors,])]) - (as.numeric(echo.data.i3[as.matrix(na.locations[na.factors,])])-1))
#   
#   for (i in 1:nbr_missing) {
#     errors.cum[2, na.locations[i,2]] <- errors.cum[2,na.locations[i,2]] + errors.this.iteration[i]
#     errors.cum[1, na.locations[i,2]] <- errors.cum[1,na.locations[i,2]] + 1}
#   
#   #errors.cum.2<- errors.cum.5<- matrix(0, nrow=2, ncol=ncol(echo.data), dimnames=list(c("n","Error"),colnames(echo.data)))
#   
#   # Calculate mean error for this iteration
#   options(warn=-1)
#   errors.mean.it<- 100*(errors.cum[2,] / errors.cum[1,]) / colMeans(data.matrix(echo.data), na.rm=T)
#   errors.mean.it[factors[-1]]<- errors.cum[2,factors[-1]] / errors.cum[1,factors[-1]]
#   errors.mean[maxit.counter, 3:(length(errors.mean.it)+1)]<- errors.mean.it[-1]
#   options(warn=0)
#   return (list(errors.mean.it, errors.mean, errors.cum))
# } 

