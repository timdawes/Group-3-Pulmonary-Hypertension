# Phosphodiesterase 5 inhibitor treatment and survival in interstitial lung disease pulmonary hypertension: a Bayesian retrospective observational cohort study
#
# Timothy JW Dawes [1], Colm McCabe [1,2], Konstantinos Dimopoulos [1,2,3], Iain Stewart [1], Simon Bax [2], 
# Carl Harries [2], Chinthaka Samaranayake [1], Aleksander Kempny [1,2,3], Philip L Molyneaux [1,4],
# Samuel Seitler [2], Thomas Semple [5], Wei Li [1,6], Peter George [1,4], Vasilis Kouranos [1,4],
# Felix Chua [1,4], Elisabetta A Renzoni [1,4], Maria Kokosi [1,4], Gisli Jenkins [1,4], Athol U Wells [1,4],
# S John Wort [1,2]*, Laura C Price [1,2,3]*
#
# [1] National Heart and Lung Institute, Imperial College London, UK
# [2] National Pulmonary Hypertension Service
# [3] Adult Congenital Heart Disease Service, Royal Brompton Hospital, London, UK
# [4] Department of Interstitial Lung Disease, Royal Brompton Hospital, London, UK
# [5] Department of Radiology, Royal Brompton Hospital, London
# [6] Department of Echocardiography, Royal Brompton Hospital, London, UK.
#
# *Joint senior authors
#
# Correspondence details: Laura C Price. laura.price@rbht.nhs.uk


   # Load packages
      options(warn=-1)
      
      suppressPackageStartupMessages({
        library(readxl)
        library(survival)
        library(Amelia)
        library(mitools)
        library(epitools)
        library(mix)
        library(glmnet)
        library(Hmisc)
        library(lubridate)
        library(GGally)
        library(survminer)
        library(car)
        library(reshape2)
        library(ordinal)
        library(MCMCglmm)
        library(splines)
        library(spBayesSurv)
        library(coda)
        library(ggplot2)
        library(tidyverse)
        library(flexsurv)
        library(RColorBrewer)
        library(scales)
        library(MatchIt)
        library(JointAI)
        
        
        # Pooling of imputed dataset for Table 1
        library(BSDA)
        library(miceadds)
        library(psfmi)
        library(DescTools)
        library(MASS)
        
        # Arranging composite ggplot figures
        library(egg)
        library(gtable)
        library(grid)
        
        # Imputation accuracy plot ("overimpute")
        library(png)
        library(ggthemes)
        
        library(rpsychi) # ind.twoway.second function
        library(km.ci)
        
        # Propensity Score Matching
        library(MatchIt)
        library(gridExtra)
        
        
      })

   


categorise<- function(cat.vec, a)
          {
            a.categorised<- rep(0, length(a))
            a.categorised[which(is.na(a)==TRUE)]<- NA
            for (i in 1:length(cat.vec)) {
              x<- which(a<=cat.vec[i])
              y<- which(a.categorised==0)
              a.categorised[intersect(x,y)]<- i}
            return(a.categorised)
          }


CPI<- function(DLCOperc, FVCperc, FEV1perc) {return(91.0 - 0.65*DLCOperc - 0.53*FVCperc + 0.34*FEV1perc)}


percent.nas<- function(m, fields.match, i)
          {
          return(length(c(
            which(is.na(m[,na.omit(fields.match)[i]])==TRUE),
            which(as.character(m[,i])=="NA"))) / (nrow(m)/100))
          }


pooling.by.Rubins.rules<- function(matrix.mean, matrix.sd, imputations, ncoefs.survival.coefficients)
{

  Vw.mat<- sapply(1:ncoefs.survival.coefficients, function(y) {sapply(1:imputations, function(x) {matrix.sd[x,y]^2})})
  colnames(Vw.mat)<- rownames(matrix.sd[[1]])
  rownames(Vw.mat)<- paste("I",1:imputations,sep="")
  
  Vb.mat<- sapply(1:ncoefs.survival.coefficients, function(y) {sapply(1:imputations, function(x) {matrix.mean[x,y]})})
  colnames(Vb.mat)<- rownames(matrix.sd[[1]])
  rownames(Vb.mat)<- paste("I",1:imputations,sep="")
  
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


# Calculate differences in median survival for each variable (i.e. an effect size for each variable in the model) + credible intervals (lower/upper)

survival.effect<- function(res.bayes,ncoefs.survival.coefficients,coef.df)
{
  res<- res.bayes[[1]][[1]]
  a<- matrix(0, nrow=ncoefs.survival.coefficients, ncol=1000)
  for (i in 1:ncoefs.survival.coefficients) {a[i,]<- rnorm(1000,mean=coef.df[[1]]$mean[i], sd=coef.df[[1]]$SD[i])}  
  res$beta<- a
  par(mfrow=c(2,5));
  wide=0.01;
  tgrid = seq(1e-10,40,wide);
  ngrid = length(tgrid);
  median.surv0<- median.surv1<- lower.surv0<- lower.surv1<- upper.surv0<- upper.surv1<- rep(0,10)
  names(median.surv0)<- names(median.surv1)<- names(lower.surv0)<- names(lower.surv1)<- names(upper.surv0)<- names(upper.surv1)<- c("DoseEqPDE5","DoseEqERA","AgeAtDiag","TAPSE","Gender","FVC","PVR","BNP","SubgroupOne","SubgroupTwo")
  
  for (i in 1:10)
  {
    cat(i)
    # Set up new dataframe of covariates (all zeros)
    newdata = data.frame(DoseEqPDE5=c(0,0), DoseEqERA=c(0,0), AgeAtDiag=c(0,0), TAPSE=c(0,0), Gender=c(0,0), FVC=c(0,0), PVR=c(0,0), BNP=c(0,0), SubgroupOne=c(0,0), SubgroupTwo=c(0,0), tstart=0)
    
    # Change one value sequentially to look at change in median survival
    newdata[2,i]<- 1
    
    p<- plot(res, xnewdata=newdata, tgrid=tgrid, PLOT=TRUE) # Shat = survival, hhat = hazard, fhat = density
    res.pooled.df<- data.frame(mean = 100 * c(p$Shat[,1], p$Shat[,2]),
                               upper = 100 * c(p$Shatup[,1], p$Shatup[,2]),
                               lower = 100 * c(p$Shatlow[,1], p$Shatlow[,2]),
                               curve = c(rep("a",ngrid), rep("b",ngrid)),
                               colour = c(rep("blue",ngrid), rep("red", ngrid)),
                               interval = rep(tgrid, times=2))
    
    median.surv0[i]<- round(tgrid[which.min(abs(p$Shat[,1]-0.5))],2)
    median.surv1[i]<- round(tgrid[which.min(abs(p$Shat[,2]-0.5))],2)
    
    lower.surv0[i]<- round(tgrid[which.min(abs(p$Shatlow[,1]-0.5))],2)
    lower.surv1[i]<- round(tgrid[which.min(abs(p$Shatlow[,2]-0.5))],2)
    
    upper.surv0[i]<- round(tgrid[which.min(abs(p$Shatup[,1]-0.5))],2)
    upper.surv1[i]<- round(tgrid[which.min(abs(p$Shatup[,2]-0.5))],2)
    
  }
  
  # Difference in survival
  median.surv.change<- median.surv1 - median.surv0
  lower.surv.change<- lower.surv1 - lower.surv0
  upper.surv.change<- upper.surv1 - upper.surv0
  
  return(list(median.surv.change=median.surv.change, lower.surv.change=lower.surv.change, upper.surv.change=upper.surv.change))
  
}


surv.median<- function(S, q,B) {sapply(q, function(x) {median(bootkm(S, x, B, pr=FALSE))})}


median.cat<- function(x)
  {
  n<- nrow(x)
  if (n==2) {x.sorted<- x[,order(x[1,], x[2,])]}
  if (n==3) {x.sorted<- x[,order(x[1,], x[2,], x[3,])]}
  if (n==4) {x.sorted<- x[,order(x[1,], x[2,], x[3,], x[4,])]}
  
  p<- ncol(x) / 2
  return(x.sorted[,round(p,0)])
  }
  


t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}


one.way.anova<- function (fields,W1,x,IDgroup)
        {
  
        group.means<- as.numeric(W1[x,(1:length(IDgroup))*2-1])
        group.ns<- unlist(lapply(IDgroup, function(x) {length(x)}))
        overall.mean<- group.means[1]
        variance<- as.numeric(W1[x,(1:length(IDgroup))*2])^2
        
        SSB<- sum(group.ns[-1] * (overall.mean - group.means[-1]) ^ 2)
        SSE<- sum(variance[-1] * (group.ns[-1]-1))
        
        k<- length(IDgroup)-1
        n<- group.ns[1]
        
        df.treatment<- k-1
        df.error<- n-k
        df.total<- n-1
        
        MS.treatment<- SSB/df.treatment
        MS.error<- SSE/df.error
        
        f.stat<- MS.treatment / MS.error
        
        pf(f.stat, df.treatment, df.error, lower.tail=FALSE)
        }


make.final.fmla<- function(fmla.final.fxd.effects, data, selectedVars, fxd.effects.factors)
          {
            initial.y<- fmla.final.fxd.effects
            add<- NULL
            
            # Fixed effects
                  for (i in 1:length(colnames(data)[selectedVars])) {
                    if ((colnames(data)[selectedVars][i] %in% fxd.effects.factors) == TRUE) {
                      add<- paste("as.factor(",colnames(data)[selectedVars][i],")",sep="")} else {
                        add<- paste(colnames(data)[selectedVars][i],sep="")}
                    if (fmla.final.fxd.effects!= initial.y) {add<- paste(" + ",add)}
                    fmla.final.fxd.effects<- paste(fmla.final.fxd.effects, add, sep="")}
            
                  fmla.final.fxd.effects<- as.formula(fmla.final.fxd.effects)
            
            # Random effects
                  fmla.final.rdm.effects<- "random = ~ tstart | IDn + "
                  initial.y<- fmla.final.rdm.effects
                  
                  for (i in 1:length(colnames(data)[selectedVars])) {
                    if ((colnames(data)[selectedVars][i] %in% fxd.effects.factors) == TRUE) {
                      add<- paste("as.factor(",colnames(data)[selectedVars][i],") | IDn",sep="")} else {
                        add<- paste(colnames(data)[selectedVars][i]," | IDn",sep="")}
                    if (fmla.final.rdm.effects!= initial.y) {add<- paste(" + ",add)}
                    fmla.final.rdm.effects<- paste(fmla.final.rdm.effects, add, sep="")}
                  
                  fmla.final.rdm.effects<- as.formula(fmla.final.rdm.effects)
            
            # Print results
                  
                  cat("\n \n Final formula \n")
                  
                  message<- paste(colnames(data)[selectedVars], collapse = " + ")
                  line<- paste(paste(rep("-", nchar(message)+1),collapse=""),"\n")
                  
                  cat(line)
                  cat("\n", message," \n")
                  cat(line)
            
            # Update log of which covariates are important
                  
            return(list(FE = fmla.final.fxd.effects, RE = fmla.final.rdm.effects))
          }




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


MeanAndSD<- function(column, Ix.at.diagnosis, imputations)
          {
            w.col<- which(colnames(Ix.at.diagnosis)==column)
            
            hist(Ix.at.diagnosis[,w.col], main=column)
            
            imputed.mean<- sapply(1:imputations, function(x) {mean(Ix.at.diagnosis[which(Ix.at.diagnosis$Imp==x),w.col])})
            imputed.summary<- sapply(1:imputations, function(x) {summary(Ix.at.diagnosis[which(Ix.at.diagnosis$Imp==x),w.col])})
            imputed.sd<- sapply(1:imputations, function(x) {sd(Ix.at.diagnosis[which(Ix.at.diagnosis$Imp==x),w.col])})
            
            Mean = mean(imputed.mean)  
            StdDev = mean(imputed.sd)
            Summary = apply(imputed.summary,1,mean)
            RawData = Ix.at.diagnosis[,w.col]
            
            return(list(Mean = Mean, SD = StdDev, Summary = Summary, RawData))
          }




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


LOCFandBplusMICE<- function(y, m, maxit, method)
          {
            y<- LOCFandB(y)  
            cat("\n Multiple imputation with chained equations in progress...")
            
            colnames(y)<- make.names(colnames(y))
            # Identify the columns which are NOT 'ID'
                numeric.column.names<- match(setdiff(colnames(y),"ID"), colnames(y))
            # Convert all non-ID columns to numeric data
                for (i in numeric.column.names) {y[,i]<- as.numeric(unlist(y[,i]))}
            # Impute
                test.imp<- mice(y[,numeric.column.names], m=m, maxit=maxit, method=method, print=F)
            # Long version, convert all 1/2 columns to 0/1, and stick back to the IDs
                results<- complete(test.imp, 'long')
                results<- zero.and.one.cols(results,3)
                #results<- cbind(ID = rep(y$ID,m),results)
                
            # Clear up any remaining zeros or NAs    
                for (i in numeric.column.names)
                {
                zeros.or.NAs<- union(which(y[,i]==0), which(is.na(y[,i])==TRUE))
                not.zeros.or.NAs<- setdiff(1:nrow(y), zeros.or.NAs)
                y[zeros.or.NAs,i]<- mean(y[not.zeros.or.NAs,i])
                }
                
            cat("..completed")
            return(results)
          }






JustAMELIA<- function(y, m, maxit, max.emburn, logs, noms, ords, idvars, all.as.numeric)
          {
  
            cat("\n Amelia: Multiple imputation with chained equations in progress...")
            
            # Identify the priors
                  # Number of subjects (n) and fields (p)
                        n<- length(unique(y$ID))
                  
                  # Find field types in which a distribution can sensibly be estimated
                        col.type<- rep(0, ncol(y))
                        for (i in 1:ncol(y)) {col.type[i]<- typeof(y[,i])}
                        cols.integer<- which(col.type=="integer")
                        cols.double<- which(col.type=="double")
                        cols.to.exclude<- c(unlist(sapply(c("Date","ID","id","Diagnosis","Gender","Outcome","futime","tstart","tstop","Dose","Treated"), function(x) {grep(x, colnames(y))})), match(c(noms, idvars), colnames(y)))
                        cols.to.exclude<- sort(unique(cols.to.exclude))
      
                        cols.priors<- sort(setdiff(c(cols.integer, cols.double), cols.to.exclude))
                        p<- length(cols.priors)
                        
                  # Define matrices for the mean/sd of each missing point, the ID/field name of each missing point
                        prior.distr<- matrix(0, nrow=p, ncol=4)
                        prior.distr.names<- cbind(rep(sort(unique(y$ID)), p), rep(colnames(y)[cols.priors], each=n))
                        
                        
                  # Find the mean and sd for each missing value using data from  1.The subject + specific field or 2. The field overall if (1) not possible
                  
                        for (i in 1:p)
                                    {
                                    field.col<- cols.priors[i]
                                    field.type.double<- typeof(y[,field.col])
                                    if (field.type.double == "double") {prior.distr[i,]<- c(0, field.col, mean(na.omit(y[,field.col])), sd(na.omit(y[,field.col])))}
                                    if (field.type.double != "double") {prior.distr[i,]<- c(0, field.col, as.numeric(c(names(sort(table(y[, field.col]), decreasing=T))[1])), 0.1)}
                                    }
                        
                  # Set bounds: Find the bounds of the current data and use this as a bound for imputation
                                bds<- apply(y, 2, range, na.rm=TRUE)
                                bds<- cbind(1:ncol(y), matrix(as.numeric(bds),ncol=2, byrow=T))
                                rownames(bds)<- colnames(y)
                                #match.rows<- na.omit(match(logs, rownames(bds)))
                                #bds[match.rows,2:3]<- log(bds[match.rows,2:3])
                                
            # Make the IDs numeric
                  y$IDn<- match(y$ID, unique(y$ID))
            
            # Impute
                  if (all.as.numeric == TRUE)
                    
                    {
                    for (i in 1:ncol(y)) {y[,i]<- as.numeric(y[,i])}
                    test.imp<- amelia(y, m=m, maxit=maxit, autopri = 1, emburn = c(1,max.emburn), idvars = idvars, #parallel = "multicore", ncpus = 4,  
                                      p2s=2, polytime = 1, intercs = FALSE, ts = "tstart", cs = "IDn", bounds = bds, priors = prior.distr)
                    
                    } else {
                    test.imp<- amelia(y, m=m, maxit=maxit, autopri = 1, emburn = c(1,max.emburn), idvars = idvars, #parallel = "multicore", ncpus = 4,  
                                      p2s=2, polytime = 1, intercs = FALSE, ts = "tstart", cs = "IDn", bounds = bds, priors = prior.distr,
                                      noms = noms, ords = ords)}
                  
                              
            cat("..completed")
            
            return(test.imp)
          }







draw.KM<- function(KMmodel, pval) 
          {
            p<- ggsurvplot(KMmodel, pval=pval, pval.method = TRUE, surv.median.line = "hv", conf.int = TRUE, risk.table = TRUE, ggtheme = theme_bw())
            
            surv_median <- as.vector(summary(KMmodel)$table[, "median"])
            df <- data.frame(x1 = surv_median, x2 = surv_median, y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
            
            #p$plot <- p$plot + 
            #  geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
            #               linetype = "dashed", size = 0.5)+ # horizontal segment
            #  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
            #               linetype = "dashed", size = 0.5) # vertical segments
            
            p$plot<- p$plot + ggplot2::annotate("text", x = surv_median[1], y = 0.01, # x and y coordinates of the text
                              label = paste(round(surv_median[1],2),"years"), size = 4)
            
            p$plot<- p$plot + ggplot2::annotate("text", x = surv_median[2], y = 0.01, # x and y coordinates of the text
                              label = paste(round(surv_median[2],2),"years"), size = 4)
            
            
            print(p)
          }


zero.and.one.cols<- function(x,start) {for (i in start:ncol(x)) {if (identical(as.numeric(as.character(names(table(x[,i])))),c(1,2))==TRUE) {x[,i]<- as.numeric(as.character(x[,i])) - 1}}; return(x)}


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


sampler<- function(t,nmax)
          {
            cat("\n Sampling..")
            if (nmax==0) {n<- 1:nrow(t)} else
            {
              j=1; n<- NULL
                  while (length(n)<nmax)
                  {
                    cat("\n ",j,"..",length(n))
                    for (i in 1:185) {
                      IDns<- which(t$IDn==i)
                      n<- c(n, which(t$IDimp == unique(t[IDns,]$IDimp)[j]))}
                      j<- j + 1
                  }
            }
              cat("done.")
              return(n)
          }
  

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
                
                dateE.columns<- grep("date.Echo",colnames(d2))
                
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





compare.demographics<- function(field, y, fig)
{
  col<- which(colnames(y)==field)
  
  All.mean<- round(mean(unlist(y[,col])),fig)
  All.sd<- round(sd(unlist(y[,col])),fig)
  
  alive<- which(y$Outcome==0)
  dead<- which(y$Outcome==1)
  
  Alive.mean<- round(mean(unlist(y[alive,col])),fig)
  Alive.sd<- round(sd(unlist(y[alive,col])),fig)
  Dead.mean<- round(mean(unlist(y[dead,col])),fig)
  Dead.sd<- round(sd(unlist(y[dead,col])),fig)
  
  a.vs.d.ttest<- round(t.test(y[alive,col], y[dead,col])$p.value,4)
  #prop.test(matrix(c(Alive.mean, 1-Alive.mean, Dead.mean, 1-Dead.mean), nrow=2))
  
  tstart.col<- which(colnames(y)=="tstart")
  tstop.col<- which(colnames(y)=="tstop")
  outcome.col<- which(colnames(y)=="Outcome")
  
  
  fit<- lapply(1:imputations, function(x) {coxph(Surv(y[,tstart.col], y[,tstop.col], y[,outcome.col]==1) ~ y[,col], data=subset(Ix.at.diagnosis, Imp==x))})
  cox.model<- summary(pool(fit))
  
  hr<- round(exp(cox.model[2]),2)
  ci.lower<- round(exp(cox.model[2] - 1.96*cox.model[3]),2)
  ci.upper<- round(exp(cox.model[2] + 1.96*cox.model[3]),2)
  hr.pvalue<- round(cox.model[6],4)
  
  return(list(All.mean = All.mean, All.sd = All.sd, Alive.mean = Alive.mean, Alive.sd = Alive.sd, Dead.mean = Dead.mean,
              Dead.sd = Dead.sd, p.value.ttest = a.vs.d.ttest, hr = hr, ci.lower = ci.lower, ci.upper = ci.upper,
              hr.pvalue = hr.pvalue))
}




image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
              if(!missing(breaks)){
                if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
              }
              if(missing(breaks) & !missing(zlim)){
                breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
              }
              if(missing(breaks) & missing(zlim)){
                zlim <- range(z, na.rm=TRUE)
                zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
                zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
                breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
              }
              poly <- vector(mode="list", length(col))
              for(i in seq(poly)){
                poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
              }
              xaxt <- ifelse(horiz, "s", "n")
              yaxt <- ifelse(horiz, "n", "s")
              if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
              if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
              if(missing(xlim)) xlim=XLIM
              if(missing(ylim)) ylim=YLIM
              plot(1,1,type="n",ylim=ylim, xlim=xlim, axes = F, xlab="", ylab="")  
              axis(1, at=seq(XLIM[1]+0.5, XLIM[2]-0.5,1), labels=seq(XLIM[1], XLIM[2]-1, 1))
              for(i in seq(poly)){
                if(horiz){
                  polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
                }
                if(!horiz){
                  polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
                }
              }
            }
