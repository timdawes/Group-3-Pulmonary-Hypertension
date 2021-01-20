# Group 3 PH Project
# Tim Dawes September 2020



# ILD survival data analysis


# Setup
          imputations<- 100
          maxit<- 3
          
          # Which variables need to be log transformed to be normal?
              logs = c("RVSP","RA.pressure","RA.areaI","RA.volumeI","Edecel","E.A","EE.","LA.areaI","LA.volumeI","FEV1",
                       "FEV1perc","FVC","FVCperc","TLcoc","TLcoperc","PVR","BNP","Ur","Creat","Bili","ALP","ALT","Fe","TFSat")
              

# Load up libraries, functions and data

          working.dir<- c("~/Dropbox/ILD-PH")
          setwd(working.dir)
          source("Code/Functions.R")
          source("Code/Col_types.R")
          source("Code/Setup.R") 
          load("Data/SetUp.RData")


# START OF CODE
          
          # Set up colours
                display.brewer.all(colorblindFriendly=TRUE)
                display.brewer.pal(20,"Paired")
                cols<- brewer.pal(20,"Paired")[c(1:10)]
                barplot(1:10, col=cols)
                
          # Data formatted for tables in write-up        
          
                # Make a matrix of data near the time of diagnosis (Ix.at.diagnosis)
                      Ix.at.diagnosis<- NULL
                      for (j in 1:nrow(d2)) {
                        cat(j)
                        a<- which(results.olps.i.amelia$imputation[[1]]$ID == d2[j,]$URN)
                        b<- which.min(abs(results.olps.i.amelia$imputations[[1]][a,]$DateOfScan - ymd(d2[j,]$dateDIAGNOSIS)))
                        for (k in 1:imputations) {Ix.at.diagnosis<- rbind(Ix.at.diagnosis, cbind(d2[j,]$URN, Imp = k, results.olps.i.amelia$imputations[[k]][a[b],]))}}

                      # Add FEV1/FVC
                            Ix.at.diagnosis$'FEV1FVC'<- 100 * Ix.at.diagnosis$FEV1 / Ix.at.diagnosis$FVC
                
                # Make a matrix of diagnoses      
                      fields<- c("AgeAtDiag","Gender","BSA","FEV1","FVC","FEV1FVC","TLcoc","TLcoperc","Kco","Kcoperc","mPAP","PCWP","CO","PVR","SR","RV.dilat","FAC","TAPSE","S","TRvel","RA.pressure","RA.areaI","PAAT","LVEF.Simpson","Edecel","E.A","EE.","LA.areaI","FC","distExTest","EmPHasis10","BNP")
                      fields.columns<- match(fields, colnames(Ix.at.diagnosis))
                      dead.IDs<- d2$URN[which(d2$Dead==1)]
                      alive.IDs<- d2$URN[which(d2$Dead==0)]
                      
                      
                      # Rows of dead/alive patients
                          rows.dead<- which(Ix.at.diagnosis$ID %in% dead.IDs)
                          rows.alive<- which(Ix.at.diagnosis$ID %in% alive.IDs)
                          rows.all<- 1:nrow(Ix.at.diagnosis)
                          
                          W1<- matrix(0, nrow=length(fields), ncol=11, dimnames=list(fields, c("All","SD","Alive","SD","Dead","SD","P","HR","Lower95CI","Upper95CI","P")))
                          W1.categories<- list()
                          W1.quantiles<- matrix(0, nrow=length(fields), ncol=9, dimnames=list(fields, c("All_LQ","All_median","All_UQ","Alive_LQ","Alive_median","Alive_UQ","Dead_LQ","Dead_median","Dead_UQ")))
                          
                          for (k in 1:3)
                          {
                            cat("\n",k,":",sep="")
                            if (k==1) {rows.in.use<- rows.all}
                            if (k==2) {rows.in.use<- rows.alive}
                            if (k==3) {rows.in.use<- rows.dead}
                            matrix.mean<- matrix.sd<- matrix.hr<- matrix.hr.sd<- matrix.lowerQTL<- matrix.median<- matrix.upperQTL<- matrix(0, nrow = imputations, ncol = length(fields), dimnames = list(paste("I",1:100,sep=""), fields))
                            
                                for (i in 1:length(fields))
                                        {
                                          cat(i) 
                                  
                                          rows.imp<- which(Ix.at.diagnosis$Imp == j)
                                          summary.data<- Ix.at.diagnosis[intersect(rows.in.use, rows.imp),]
                                          if ((fields[i] %in% logs)==TRUE) {summary.data[,fields.columns[i]]<- log(summary.data[,fields.columns[i]])}
                                          
                                                       if (class(Ix.at.diagnosis[,fields[i]]) %in% c("numeric","difftime"))
                                                                    {
                                                                        for (j in 1:imputations)
                                                                        {
                                                                          # Store the means and sd for each imputation for this field
                                                                              matrix.mean[j,i]<- mean(summary.data[,fields.columns[i]])
                                                                              matrix.sd[j,i]<- sd(summary.data[,fields.columns[i]])
                                                                              
                                                                              matrix.lowerQTL[j,i]<- quantile(summary.data[,fields.columns[i]], probs=0.25)
                                                                              matrix.median[j,i]<- quantile(summary.data[,fields.columns[i]], probs=0.5)
                                                                              matrix.upperQTL[j,i]<- quantile(summary.data[,fields.columns[i]], probs=0.75)
                                                                              
                                                                              
                                                                          # Store the HR and sd for each imputation for this field
                                                                              matrix.hr[j,i]<- summary(coxph(Surv(summary.data$futime, summary.data$Outcome) ~ scale(summary.data[,fields[i]])))$conf.int[2]
                                                                              matrix.hr.sd[j,i]<- sqrt(coxph(Surv(summary.data$futime, summary.data$Outcome) ~ scale(summary.data[,fields[i]]))$var)
                                                                        }
                                                                    }
                                            
                                                        if (class(Ix.at.diagnosis[,fields[i]])=="factor")
                                                                    {
                                                                      matrix.cat<- NULL
                                                                      counter<- which(c("Gender","SR","RV.dilat","FC") == fields[i]) + (k-1)*4
                                                                      
                                                                        for (j in 1:imputations)
                                                                        {
                                                                          # Store the table of categorical data for this imputation only for this field
                                                                                matrix.cat<- cbind(matrix.cat, unclass(table(as.numeric(as.character(Ix.at.diagnosis[intersect(rows.in.use, rows.imp), fields.columns[i]])))))
                                                                          # Store the HR and SD for each imputation for this categorical field
                                                                                matrix.hr[j,i]<- summary(coxph(Surv(summary.data$futime, summary.data$Outcome) ~ as.numeric(summary.data[,fields[i]])))$conf.int[2]
                                                                                matrix.hr.sd[j,i]<- sqrt(coxph(Surv(summary.data$futime, summary.data$Outcome) ~ as.numeric(summary.data[,fields[i]]))$var)
                                                                        }
                                                                      
                                                                          # Store the table of categorical data generated above (matrix.cat) for each imputation so that it can be analysed below
                                                                                colnames(matrix.cat)<- paste("I",1:imputations, sep="")
                                                                                W1.categories[[counter]]<- matrix.cat
                                                                                names(W1.categories)[[counter]]<- rep(c("Gender","SR","RV.dilat","FC"),3)[counter]
                                                                      }          
                                        } # Fields loop
                            
                            
                                          # Store the mean, SD, HR and CI data in the final matrix for continuous variables
                                              start<- 1 + (k-1)*2
                                              stop<- start + 1
                                              W1[, start:stop]<- data.matrix(pooling.by.Rubins.rules(matrix.mean, matrix.sd, imputations)[,1:2])
                                              
                                              start<- 1 + (k-1)*3
                                              stop<- start + 2
                                              W1.quantiles[,start:stop]<- cbind(apply(matrix.lowerQTL,2,mean), apply(matrix.median,2,mean), apply(matrix.upperQTL,2,mean))
                                  } # k loop
                          
                          
                          

                            # Columns 7-11 for continuous variables
                                          
                                          # Column 7: P-values for comparing alive and dead patients
                                                p.vals<- sapply(1:length(fields), function(x) {tsum.test(as.numeric(W1[x,3]), as.numeric(W1[x,4]), nrow(d2) - sum(d2$Dead), as.numeric(W1[x,5]), as.numeric(W1[x,6]), sum(d2$Dead), 0, alternative = "two.sided", var.equal=FALSE)$p.value})
                                                W1[,7]<- p.vals
                                                
                                          # Column 8: Hazard ratio for survival  
                                                hazards<- pooling.by.Rubins.rules(matrix.hr, matrix.hr.sd, imputations)
                                                W1[, 8]<- hazards[,1]
                                                
                                          # Columns 9,10: 95% Confidence intervals
                                                W1[, 9] <- hazards[,3]
                                                W1[, 10] <- hazards[,4]
                                                
                                          # Column 11: p-value
                                                W1[, 11] <- pnorm(q = abs(hazards[,1] - 1), sd = hazards[,2], lower.tail = FALSE)
                                          
                                          
                            # Columns 1,3,5 & 7-11 for categorical variables
                                                
                                          for (i in 1:4)
                                          {
                                            # Select the right data for each categorical variable
                                                all.cats <- W1.categories[[i]]
                                                alive.cats<- W1.categories[[i+4]]
                                                dead.cats<- W1.categories[[i+8]]
                                                
                                                row.W1<- match(names(W1.categories)[[i]], rownames(W1))
                                                
                                                # Find the median values (equivalent to the means for the continuous variables)
                                                    all.mode<- median.cat(all.cats)
                                                    alive.mode<- median.cat(alive.cats)
                                                    dead.mode<- median.cat(dead.cats)
                                                
                                                # Columns 1,3,5
                                                    # Reverse the listings of gender, SR and RV.dilat so that the final results appears as "Yes/No"
                                                          if (i!=4) {all.mode<- rev(all.mode); alive.mode<- rev(alive.mode); dead.mode<- rev(dead.mode)}
                                                    
                                                    W1[row.W1, 1]<- paste(all.mode, collapse="/")
                                                    W1[row.W1, 3]<- paste(alive.mode, collapse="/")
                                                    W1[row.W1, 5]<- paste(dead.mode, collapse="/")
                                                    
                                                    W1[row.W1, 2]<- W1[row.W1, 4]<- W1[row.W1, 6]<- "-"
                                                     
                                                    
                                                # Column 7 
                                                    
                                                      if (i<4) {
                                                                # For age, gender, RV.dilat: do a chi-squared on each imputation and pool the results    
                                                                    chisq.imp<- df.imp<- rep(0, imputations)
                                                                    
                                                                    for (j in 1:imputations)
                                                                    {
                                                                      data.chisq<- cbind(alive.cats[,j], dead.cats[,j])
                                                                      chisq.imp[j]<- chisq.test(data.chisq)$statistic
                                                                      df.imp[j]<- chisq.test(data.chisq)$parameter
                                                                    }
                                                                    
                                                                    W1[row.W1,7]<- micombine.chisquare(chisq.imp, df=2, display=FALSE)[2]
                                                                }
                                                       
                                                      if (i==4) {   
                                                                # For FC: do a Cochrane Armitage test on each imputation and pool the results
                                                                    chisq.imp<- df.imp<- rep(0, imputations)
                                                            
                                                              
                                                                    for (j in 1:imputations)
                                                                    {
                                                                      data.chisq<- cbind(alive.cats[,j], dead.cats[,j])
                                                                      chisq.imp[j]<- CochranArmitageTest(data.chisq)$statistic ^ 2
                                                                      df.imp[j]<- CochranArmitageTest(data.chisq)$parameter
                                                                    }
                                                      
                                                                    W1[row.W1,7]<- micombine.chisquare(chisq.imp, df=3, display=TRUE)[2]
                                                                  }          
                                                # Columns 8-11 (hazards, CIs and p-values)
                                                    
                                                    # FC treated as a continuous variable, otherwise hazard ratios and CIs can't be produced
                                                      predictor.variable<- names(W1.categories)[i]
                                                      fmla<- as.formula(paste("Surv(futime, Outcome) ~ ", predictor.variable, sep=""))
                                                      pool_coxr<- psfmi_coxr(formula = fmla, data = Ix.at.diagnosis, predictors = predictor.variable, nimp = imputations, impvar = "Imp", method = "D1", p.crit = 1)
                                                      hazards<- pool_coxr$RR_model$`Step 1 - no variables removed -`
                                                      W1[row.W1,8:11]<- unlist(hazards[c(7:9,6)])
                                                      
                                                }
                                       
                                                
                                    # Finally, convert the logged fields (l) into confidence intervals       
                                                
                                                for (k in 1:3)
                                                {
                                                  start<- 1 + (k-1)*2
                                                  stop<- start + 1
                                                  
                                                  start2<- 1 + (k-1)*3
                                                  stop2<- start2 + 2
                                                  
                                                  sub.in.IQR.4.these.fields<- which(fields %in% logs)
                                                  
                                                    for (l in sub.in.IQR.4.these.fields) 
                                                    {
                                                      lowerCI<- round(exp(W1.quantiles[l,start2]),2)
                                                      medianCI<- round(exp(W1.quantiles[l,start2+1]),2)
                                                      upperCI<- round(exp(W1.quantiles[l,stop2]),2)
                                                      
                                                      W1[l,stop]<- paste(lowerCI,upperCI, sep="-")
                                                      W1[l,start]<- exp(as.numeric(as.character(W1[l,start])))
                                                    }
                                                }
                                                
                       
                             
                                              
                                              
                
                            # Save the output
                                  write.csv(W1, file="Data/Table1.csv", col.names = T, row.names = T)
                              
                                  
                      
                      
                # Identify the majority diagnosis for each patient from the imputed datasets (x100) - Supplementary Table 1
                      # Choose a field of interest and its column number
                            field<- "RV.dilat"
                            col.field<- which(colnames(Ix.at.diagnosis)==field)
                            
                      # Choose a subject group of interest, check it's the length you expect and form a vector of the same length
                            ID.group<- 1:184
                            #ID.group<- which(d2$ERA==0)
                            length(ID.group)
                            majority.diagnosis<- rep(0, length(ID.group))
                      
                      # Find all the imputations for each subject and choose the majority diagnosis
                            for (i in 1:length(ID.group)){
                              IDn.i<- which(Ix.at.diagnosis$IDn == ID.group[i])
                              majority.diagnosis[i]<- names(sort(table(Ix.at.diagnosis[IDn.i,col.field]), decreasing=T)[1])}
                            
                      # Tabulate the results
                            table(majority.diagnosis)
                      
                            
                            
                # W2: Describe the PFTs in groups 3.2 and 3.1
                
                      fields<- c("FVCperc","TLcoperc","Kcoperc", "FEV1perc","FEV1FVC","TLcoperc","Kcoperc")
                      patient.group<- c(rep(3.2,3), rep(3.1,4))
                      
                      W2<- matrix(0, nrow=length(fields), ncol=3, dimnames=list(paste(patient.group,"_",fields,sep=""), c("mean","sd","comb")))
                      fit.pooled<- list()
                      
                      for (i in 1:length(fields))
                      {
                        fmla<- as.formula(paste("as.numeric(",fields[i],")~TreatedPDE5"))
                        fit<- lapply(1:imputations, function(x) {lm(fmla, data=subset(Ix.at.diagnosis, Imp==x & Diagnosis==patient.group[i]))})
                        fit.pooled[[i]]<- pool(fit)
                        summ<- summary(fit.pooled[[i]])
                        W2[i,1:2]<- c(format(round(summ[1,2],1),nsmall=1), format(round(summ[1,3]*sqrt(sum(Ix.at.diagnosis$Imp==1 & Ix.at.diagnosis$Diagnosis==patient.group[i])),1),nsmall=1))
                      }
                      
                      W2[,3]<- paste(W2[,1],"±",W2[,2],sep="")
          
                 
                # W2.1: Describe the RHC data in the whole group
          
                      fields<- c("mPAP","PVR","CO","TAPSE","S")
                      
                      W2.1<- matrix(0, nrow=length(fields), ncol=3, dimnames=list(fields, c("mean","sd","comb")))
                      
                      for (i in 1:length(fields))
                      {
                        fmla<- as.formula(paste("as.numeric(",fields[i],")~TreatedPDE5"))
                        fit<- lapply(1:imputations, function(x) {lm(fmla, data=subset(Ix.at.diagnosis, Imp==x))})
                        fit.pooled[[i]]<- pool(fit)
                        summ<- summary(fit.pooled[[i]])
                        W2.1[i,1:2]<- c(round(summ[1,2],1), round(summ[1,3]*sqrt(sum(Ix.at.diagnosis$Imp==1)),1))
                      }
                      
                      W2.1[,3]<- paste(W2.1[,1],"±",W2.1[,2],sep="")
                
                # How many have PH vs 'severe PH' as defined by Werner Seeger's paper (mPAP>=35 or mPAP>=25 and CI<2.0 L/min/m2)
                      length(which(Ix.at.diagnosis$mPAP>=35))
                      severePH<- matrix(0, nrow = nrow(d2), ncol=1)
                      for (i in 1:nrow(d2))
                      {
                        IDs<- which(Ix.at.diagnosis$ID==d2$URN[i])
                        severePH[i]<- (mean(Ix.at.diagnosis$CO[IDs] / Ix.at.diagnosis$BSA[IDs]) < 2.0) | (mean(Ix.at.diagnosis$mPAP[IDs]) >= 35)
                      }
                      
                      table(severePH)
                      
                      
                # Data for PDE5+ vs PDE5- comparison table
                      fields<- c("AgeAtTest","distExTest","FEV1perc","FVCperc","TLcoperc","Kcoperc","mPAP","PVR","CO","TAPSE","S","BNP")
                      W3<- matrix(0, nrow=length(fields), ncol=3, dimnames=list(fields, c("PDE+","PDE-","p")))
                      for (i in 1:length(fields))
                      {
                        fmla<- as.formula(paste("as.numeric(",fields[i],")~TreatedPDE5"))
                        fit<- lapply(1:imputations, function(x) {lm(fmla, data=subset(Ix.at.diagnosis, Imp==x))})
                        fit.pooled[[i]]<- pool(fit)
                        summ<- summary(fit.pooled[[i]])
                        W3[i,]<- c(round(sum(summ[,2]),1), round(summ[1,2],1), round(summ[2,6],2))
                      }
                      
                # Data for ERA+ vs ERA- comparison table
                      fields<- c("AgeAtDiag","BSA","FEV1","FVC","FEV1FVC","TLcoc","Kco","mPAP","PCWP","CO","PVR","FAC","TAPSE","S","TRvel","RA.pressure","RA.areaI","PAAT","LVEF.Simpson",
                                 "Edecel","E.A","EE.","LA.areaI","distExTest","EmPHasis10","BNP")
                      W4<- matrix(0, nrow=length(fields), ncol=5, dimnames=list(fields, c("ERA+","SD+","ERA-","SD-","p")))
                      
                      for (i in 1:length(fields))
                      {
                        # Fit with linear model directly
                            fmla.lin<- as.formula(paste("as.numeric(",fields[i],")~TreatedERA"))
                            fit.lin<- lapply(1:imputations, function(x) {lm(fmla.lin, data=subset(Ix.at.diagnosis, Imp==x))})
                            
                        # Fit with linear model with log transform    
                            fmla.log<- as.formula(paste("log(as.numeric(",fields[i],"))~TreatedERA"))
                            fit.log<- lapply(1:imputations, function(x) {lm(fmla.log, data=subset(Ix.at.diagnosis, Imp==x))})
                        
                        # Calculate both the variance and the IQR 
                            field.column<- match(fields[i], colnames(Ix.at.diagnosis))
                            vars<- apply(matrix(unlist(sapply(c(0,1), function (t) {lapply(1:imputations, function(x) {var(subset(Ix.at.diagnosis, Imp==x & TreatedERA==t)[,field.column])})})),ncol=2,byrow=F),2,mean)
                            IQRs<- apply(matrix(unlist(sapply(c(0,1), function (t) {lapply(1:imputations, function(x) {quantile(subset(Ix.at.diagnosis, Imp==x & TreatedERA==t)[,field.column], probs=c(0.25, 0.75))})})),ncol=4,byrow=T),2,median)
                            
                        # Pool the models  
                            summ.log<- summary(pool(fit.log))
                            summ.lin<- summary(pool(fit.lin))
                            
                        # Work out if the field is skewed (i.e. in the 'logs' list) and use the correct data accordingly
                            if ((fields[i] %in% logs) == FALSE)
                              {W4[i,]<- round(c(sum(summ.lin[,2]), sqrt(vars[1]), summ.lin[1,2], sqrt(vars[2]), summ.lin[2,6]),2)} else 
                              {W4[i,]<- c(round(exp(summ.log[1,2])*exp(summ.log[2,2]),2), paste(round(IQRs[1:2],2), collapse="-"), round(exp(summ.log[1,2]),2), paste(round(IQRs[3:4],2), collapse="-"), round(summ.log[2,6],2))}
                      }
                      
                
                
                # Data for describing the length and peak doses of treatment prescribing
                      
                      # Median treatment duration
                            w1.col<- which(colnames(d2)=="Drug1Duration")
                            w2.col<- which(colnames(d2)=="Drug16Duration")
                            w3.col<- which(colnames(d2)=="DrugName1")
                            w4.col<- which(colnames(d2)=="DrugName16")
                            w5.col<- which(colnames(d2)=="DrugDose1")
                            w6.col<- which(colnames(d2)=="DrugDose16")
                            w7.col<- which(colnames(d2)=="date.DrugStart1")
                            w8.col<- which(colnames(d2)=="date.DrugStart16")
                            
                            
                            
          
          
                      # Choose which treatment group you're interested in 
                            #Rx.given.IDs<- unique(union(which(d2$PDE==1), which(d2$ERA==1))) # On PDE5is AND/OR an ERA
                            #Rx.given.IDs<- intersect(which(d2$PDE==1), which(d2$ERA==1)) # On a PDE5i AND an ERA
                            Rx.given.IDs<- which(d2$PDE==1) # On a PDE5i
                            #Rx.given.IDs<- which(d2$ERA==1) # On an ERA
                            
                      # Find median and IQR dose duration for the patients on the drugs (Rx.given.IDs) for the periods of the drugs in question (match statement below)
                            Rx.length<- NULL
                            for (i in 1:length(Rx.given.IDs)) 
                                  {
                                    drugs.in.group<- which(is.na(match(d2[Rx.given.IDs[i],w3.col:w4.col], c(1,2,3,4)))==FALSE)
                                    Rx.length<- c(Rx.length, na.omit(c(data.matrix(d2[Rx.given.IDs[i],(w1.col:w2.col)[drugs.in.group]]))))
                                  }
                            summary(Rx.length)
                            
                      # How many patients are on a particular drug?
                            drug.ID.number<- 1
                            sum(apply(data.matrix(d2[,w3.col:w4.col]),1,function(x) {drug.ID.number %in% x}))
                            
                      # How many patients achieve maximum dose?
                            dose.achieved<- NULL
                            drug.ID.number<- 3
                            
                            for (i in 1:length(Rx.given.IDs)) 
                                  {
                                    drugs.in.group<- which(is.na(match(d2[Rx.given.IDs[i],w3.col:w4.col], c(drug.ID.number)))==FALSE)
                                    if (length(drugs.in.group)>0) {dose.achieved<- c(dose.achieved, max(na.omit(c(data.matrix(d2[Rx.given.IDs[i],(w5.col:w6.col)[drugs.in.group]])))))}
                                  }
                            summary(dose.achieved)
                            
                      # How many patients were up-titrated?
                            up.titrated<- down.titrated<- 0
                            counter<- drug.ID.number<- 1
                            up.titration.df<- data.frame(IDn=rep(0,24), Date = rep(as.Date("1900-01-01"), 24))
                            
                            
                            for (i in 1:length(Rx.given.IDs)) 
                                  {
                            
                              # Identify drug dosing of the drug of interest (e.g. 1=sild) in the patient of interest
                                    drugs.in.group<- which(is.na(match(d2[Rx.given.IDs[i],w3.col:w4.col], c(drug.ID.number)))==FALSE)
                              # Find the drug doses for the drug of interest in the patient of interest (g = a vector of drug doses)
                                    g<- na.omit(c(data.matrix(d2[Rx.given.IDs[i],(w5.col:w6.col)[drugs.in.group]])))
                              # Find how many times the drug doses are up-/down-titrated by:
                                    # 1. frame-shifting the vector 'g' along by one
                                    # 2. Substracting the frame-shifted version from the original vector
                                    # 3. Removing the first and last values as they're artefactual from the process
                                    # 4. Finding negative (down-titration) or positive (up-titration) values
                                    # 5. Summing up the number of positive/negative values from the vector
                                    
                                    if (length(g)>1) 
                                      {                       #    5   4      2    1           3  
                                      up.titrated<- up.titrated + sum(sign((g - c(0,g))[-c(1,length(g)+1)]) > 0)
                                      down.titrated<- down.titrated + sum(sign((g - c(0,g))[-c(1,length(g)+1)]) < 0)
                                      up.titration.df$IDn[counter]<- Rx.given.IDs[i]
                                      up.titration.df$Date[counter]<- ymd(c(d2[Rx.given.IDs[i], w7.col])[[1]])
                                      counter<- counter + 1
                                      }
                                    
                                  }
                            
                      
          
                      
                        # Compare patients taking PDE5i with patients not taking PDE5i
                            d2.PDE<- which(d2$PDE==1)
                            d2.notPDE<- which(d2$PDE==0)
                            
                            T1<- data.frame(Group=c("Demographic","Functional","Spirometry","","","","RHC","","","Echo","","Bloods"),
                                            Measure=c("Age (years)","6MWD (m)","FEV1","FVC","TLCO","KCO","mPAP","PVR","CO (L/min)","TAPSE (cm)","S' (cm/s)","BNP"),
                                            PDE5iplus=as.character(format(round(W3[,1],1), nsmall=1)),
                                            PDE5iminus=as.character(format(round(W3[,2],1), nsmall=1)),
                                            pvalue=as.character(format(round(W3[,3],2), nsmall=2)), stringsAsFactors=FALSE)
                            
                            T1$pvalue[which(T1$pvalue=="0.00")]<- c("<0.01")
                            T1$pvalue[which(T1$pvalue=="0.00")]<- 0.01
                            
                            colnames(T1)<- c("Group","Measure","PDE5i+","PDE5i-","p-value")
                            rownames(T1)<- NULL
                            
                            
                        # Compare patients taking ERA with patients not taking ERAs
                            d2.ERA<- which(d2$ERA==1)
                            d2.notERA<- which(d2$ERA==0)
                            
                            T1b<- data.frame(Group=c("Demographic","Functional","Spirometry","","","","RHC","","","Echo","","Bloods"),
                                            Measure=c("Age (years)","6MWD (m)","FEV1","FVC","TLCO","KCO","mPAP","PVR","CO (L/min)","TAPSE (cm)","S' (cm/s)","BNP"),
                                            ERAplus=as.character(format(round(W4[,1],1), nsmall=1)),
                                            ERAplussd=as.character(format(round(W4[,2],1), nsmall=1)),
                                            ERAminus=as.character(format(round(W4[,3],1), nsmall=1)),
                                            ERAminussd=as.character(format(round(W4[,4],1), nsmall=1)),
                                            pvalue=as.character(format(round(W4[,5],2), nsmall=2)), stringsAsFactors=FALSE)
                            
                            T1b$pvalue[which(T1b$pvalue=="0.00")]<- c("<0.01")
                            T1b$pvalue[which(T1b$pvalue=="0.00")]<- 0.01
                            
                            colnames(T1b)<- c("Group","Measure","ERA+","ERA+sd","ERA-","ERA-sd","p-value")
                            rownames(T1b)<- NULL
                            T1b                  
                            
                           
          # Odds ratio of being gender vs treatment (i.e. women vs PDE5 association)
                          # Contingency table:  Not treated with PDE5i        Treated with PDE5i
                          # Male                  80                          31
                          # Female                36                          37
                            
                        contingency.table.gender.PDE5i<- matrix(c(80,31,36,37), nrow=2, byrow=T, dimnames=list(c("Male","Female"),c("-","+")))
                        contingency.table.gender.ERA<- matrix(c(100,11,64,9), nrow=2, byrow=T, dimnames=list(c("Male","Female"),c("-","+")))
                        contingency.table.SR.ERA<- matrix(c(100,11,64,9), nrow=2, byrow=T, dimnames=list(c("Male","Female"),c("-","+")))
                        contingency.table<- matrix(0, nrow=4, ncol=2, dimnames=list(c("FC1","FC2","FC3","FC4"),c("-","+")))  
                          
                            for (i in 1:nrow(d2))
                            {
                              IDs<- which(Ix.at.diagnosis$IDn==i)
                              criterion<- median(as.numeric(as.character(Ix.at.diagnosis$FC[IDs])))
                              drug<- median(Ix.at.diagnosis$TreatedPDE5[IDs])
                              contingency.table[(criterion),(drug+1)]<- contingency.table[(criterion),(drug+1)] + 1 
                            }
                          
                        oddsratio(contingency.table)
                        
             
          # Treatment by diagnostic subgroup
                        subgroups.on.PDE5<- matrix(0, nrow=2, ncol=6, dimnames=list(c("PDE5-","PDE5+"), names(table(d2$Diagnosis_sub))))
                        subgroups.on.PDE5[1,c(1:3,6)]<- table(d2$Diagnosis_sub[which(d2$PDE==0)])
                        subgroups.on.PDE5[2,1:5]<- table(d2$Diagnosis_sub[which(d2$PDE==1)])
                        round(100*subgroups.on.PDE5 / apply(subgroups.on.PDE5,1,sum))
                        chisq.test(subgroups.on.PDE5)
                        
                        subgroups.on.ERA<- matrix(0, nrow=2, ncol=6, dimnames=list(c("ERA-","ERA+"), names(table(d2$Diagnosis_sub))))
                        subgroups.on.ERA[1,]<- table(d2$Diagnosis_sub[which(d2$ERA==0)])
                        subgroups.on.ERA[2,c(1,2,5)]<- table(d2$Diagnosis_sub[which(d2$ERA==1)])
                        round(100*subgroups.on.ERA / apply(subgroups.on.ERA,1,sum))
                        
                        chisq.test(subgroups.on.ERA)
                        
                        
          # Follow-up times of those who died
                        t<- table(ceiling(d2$FUTimeYears[which(d2$Dead==1)]))
                        round((100 * t) / sum(d2$Dead))
                        
          # KM plot by PDE5+ / PDE5-  
                        # NB printed at 5 x 8 inch pdfs, and then zoomed to 37% in Word
                        
                        # Choose a specific subgroup?
                              COPD.IDn<- which(d2$Diagnosis_sub=="3.1")
                              IPF.IDn<- which(d2$Diagnosis_sub=="3.2")
                              Other.IDn<- setdiff(1:nrow(d2), union(COPD.IDn, IPF.IDn))
                              All.IDn<- 1:nrow(d2)
                              
                              
                        # Choose a subgroup and a label to remind you!
                              IDn<- COPD.IDn
                              ID.label<- "COPD"
                              
                        km<- with(d2[IDn,], Surv(FUTimeYears, Dead))
                        km_PDE5<- survfit(Surv(FUTimeYears, Dead) ~ PDE, data = d2[IDn,])
                        p<- ggsurvplot(km_PDE5, data = d2[IDn,], size = 1.5,                 # change line size
                          palette = cols[c(2,6)],# custom color palettes
                          conf.int = TRUE,          # Add confidence interval
                          conf.int.alpha = 0.4,
                          pval = TRUE,              # Add p-value
                          risk.table = TRUE,        # Add risk table
                          risk.table.col = "strata",# Risk table color by groups
                          legend.labs = c("PDE5-", "PDE5+"), 
                          censor = FALSE,
                          surv.median.line = "hv",
                          xlim = c(0,10),
                          break.x.by = 1,
                          break.y.by = 0.1,
                          xlab = "Time (years)",
                          #fun = "cumhaz",
                          legend = "none", # Change legend labels
                          ggtheme = theme_classic()      # Change ggplot2 theme
                        )
                        p
                        surv_median(km_PDE5)
                        
                        df<- data.frame(surv=summary(km_PDE5)$upper[1:16], time=summary(km_PDE5)$time[1:16])
                        fit<- lm(log(surv) ~ time, data=df)
                        seq(1,10,0.1)[which.min(exp(predict(fit, data.frame(time=seq(1,10,0.1)))) - 0.5)]
                        
                        # Export as 6 x 6" .pdf for Figure 2 resize to 50% in Word
                        # Export as 5 x 8" .pdf for Figure 5 resize to 35% in Word
                        
                        
                  # KM plot by ERA+ / ERA-
                        km<- with(d2[IDn,], Surv(FUTimeYears, Dead))
                        km_ERA<- survfit(Surv(FUTimeYears, Dead) ~ ERA, data = d2[IDn,])
                        p<- ggsurvplot(km_ERA, data = d2[IDn,], size = 1.5,                 # change line size
                                       palette = cols[c(4,10)],# custom color palettes
                                       conf.int = TRUE,          # Add confidence interval
                                       conf.int.alpha = 0.4,
                                       pval = TRUE,              # Add p-value
                                       risk.table = TRUE,        # Add risk table
                                       risk.table.col = "strata",# Risk table color by groups
                                       legend.labs = c("ERA-", "ERA+"), 
                                       censor = FALSE,
                                       surv.median.line = "hv",
                                       xlim = c(0,10),
                                       break.x.by = 1,
                                       break.y.by = 0.1,
                                       xlab = "Time (years)",
                                       #fun = "cumhaz",
                                       legend = "none", # Change legend labels
                                       ggtheme = theme_classic()      # Change ggplot2 theme
                        )
                        surv_median(km_ERA)
                        p
                        
                        
          # Survival by subtype of group 3 PH
                        km<- with(d2, Surv(FUTimeYears, Dead))
                        d2$Diagnosis_sub_grouped<- d2$Diagnosis_sub
                        d2$Diagnosis_sub_grouped[which(d2$Diagnosis_sub_grouped %in% c("3.3","3.4","3.5","3.7"))]<- "Other"
                        
                        km_group<- survfit(Surv(FUTimeYears, Dead) ~ Diagnosis_sub_grouped, data = d2)
                        p<- ggsurvplot(km_group, data = d2, size = 1.5,                 # change line size
                          palette = cols[c(10,6,4)],
                          conf.int = TRUE,          # Add confidence interval
                          conf.int.alpha = 0.3,
                          pval = TRUE,              # Add p-value
                          risk.table = TRUE,        # Add risk table
                          risk.table.col = "strata",# Risk table color by groups
                          legend.labs = c("Group 3.1", "Group 3.2","Groups 3.3-3.7 "), 
                          censor = FALSE,
                          surv.median.line = "hv",
                          xlim = c(0,10),
                          break.x.by = 1,
                          break.y.by = 0.1,
                          xlab = "Time (years)",
                          #fun = "cumhaz",
                          legend = "none", # Change legend labels
                          ggtheme = theme_classic()      # Change ggplot2 theme
                        )
          
                        surv_median(km_group)
                                      p
                                                                                                                       
          # Association of PDE5 +/- with Functional Markers
                        # List of markers of interest (ns) and their labels( ns.labels) and a count of how many their are (n)
                              ns<- c("distExTest","FC","EmPHasis10","TAPSE","RA.pressure","CO","PVR","mPAP","BNP","Kcoperc") 
                              ns.labels<- c("6MWD","FC","emPHasis-10","TAPSE","RAP","CO","PVR","mPAP","BNP","Kco")
                              
                              n<- length(ns)
                              timeframe<- "after" # "after" or "all"
                        
                        # Construct formulae to analyses these markers
                              fmla.mm<- paste(rep(ns,times=2), c("~ as.factor(Gender) + AgeAtTest + "), rep(c("DoseEqPDE5","DoseEqERA"),each=n), c("+ (1 | IDn)"))
                              #fmla.mm<- paste(rep(ns,times=2), c("~ "), rep(c("DoseEqPDE5","DoseEqERA"),each=n), c("+ (1 | IDn)"))
                        
                        # Make matrices to store the analysis results
                              matrix.mean<- matrix.sd<- matrix(0, nrow=imputations, ncol=length(fmla.mm), dimnames=list(paste("I",1:imputations,sep=""), paste(rep(ns,time=2), rep(c("_PDE5","_ERA"),each=n), sep="")))
                        
                        # Find the corresponding columns from the data
                              cols.mm<- match(ns, colnames(results.olps.i.amelia$imputations[[1]]))
                        
                        # Choose a specific subgroup of patients?
                            COPD.IDn<- which(d2$Diagnosis_sub!="3.2")
                            IPF.IDn<- which(d2$Diagnosis_sub=="3.2")
                            All.IDn<- 1:nrow(d2)
                            
                            IDn<- COPD.IDn
                            ID.label<- "COPD"
                            
                        # Loop through each formula analysing it
                              for (j in 1:length(fmla.mm))
                                {
                                if (j==1) {cat("\n Association between treatment and functional markers in this subgroup:  ",ID.label, " (n=",length(IDn),")  patients.\n Only data acquired at time-points ",timeframe," up-titration of doses will be used \n", sep="")} else {cat(", ")}
                                cat(c(ns,ns)[j])
                                
                                  for (x in 1:imputations) 
                                          {
                                          IDn.long<- which(results.olps.i.amelia$imputations[[1]]$IDn %in% IDn)
                                          sc.data<- results.olps.i.amelia$imputations[[x]][IDn.long,cols.mm]
                                          
                                                for (i in 1:ncol(sc.data)) {sc.data[,i]<- as.numeric(as.character(sc.data[,i]))}
                                                sc.data<- scale(sc.data)
                                                sc.data<- cbind(sc.data, results.olps.i.amelia$imputations[[x]][IDn.long, match(c("Gender","AgeAtTest","DoseEqPDE5","DoseEqERA","DateOfScan","IDn"), colnames(results.olps.i.amelia$imputations[[x]]))])
                                                
                                          # Remove data taken after an up-titration
                                                data.after.up.titration<- NULL
                                                for (y in 1:nrow(up.titration.df)) {data.after.up.titration<- c(data.after.up.titration, intersect(which(sc.data$IDn == up.titration.df$IDn[y]), which(sc.data$DateOfScan >= up.titration.df$Date[y])))}
                                                if (timeframe=="before") {sc.data<- sc.data[-data.after.up.titration,]}
                                                if (timeframe=="after") {sc.data<- sc.data[data.after.up.titration,]}
                                                
                                                which.col<- match(c(ns,ns)[j], colnames(sc.data))
                                                sc.data[,which.col]<- sc.data[,which.col] + rnorm(nrow(sc.data),0,0.00001)
                                              
                                                fit_me<- blmer(as.formula(fmla.mm[j]), data = sc.data, cov.prior = gamma)
                                          
                                          # Identify the coefficients associated with the drug effect and store the beta and se
                                                drug.effect.row<- grep("DoseEq", rownames(summary(fit_me)$coefficients))
                                                matrix.mean[x,j]<- summary(fit_me)$coefficients[drug.effect.row,1]
                                                matrix.sd[x,j]<- summary(fit_me)$coefficients[drug.effect.row,2]
                                          }
                                }
                        
                        
                        # Pool the results
                              results.mm<- pooling.by.Rubins.rules(matrix.mean, matrix.sd, imputations)
                              results.mm$Drug<- rep(c("PDE5i","ERA"),each=n)
                              results.mm$Outcome<- rep(ns.labels, times=2)
                              
                              results.after.up.titration<- results.mm
                              #results.before.up.titration<- results.mm
                              
                              #results.mm.BA<- rbind(results.before.up.titration, results.after.up.titration)
                              results.mm.BA.IPF<- rbind(results.before.up.titration, results.after.up.titration)
                              
                              
                              
                        # Draw Forest plot
                        
                              results.for.Forest.plot<- results.mm.BA.IPF
                              
                              # Upper and lowers
                                  o<- c(1,7,10,5,4,6,3,9,8,2)
                                  #o<- order(results.for.Forest.plot$mean[1:n])
                                  o<- rep(o,4) + rep(seq(0,30,10),each=n)
                                  #o<- 1:40
                                  lower<- results.for.Forest.plot$LowerCI[o]
                                  mean<- results.for.Forest.plot$mean[o]
                                  upper<-  results.for.Forest.plot$UpperCI[o]
                                  pvals<- results.for.Forest.plot$PerCentChance[o]
                                  pvals[which(sign(results.for.Forest.plot$mean[o])<0)]<- 100-pvals[which(sign(results.for.Forest.plot$mean[o])<0)]
                                  labels<- results.for.Forest.plot$Outcome[o]
                                  
                                  # x-limits
                                  #xmax<- ceiling(max(abs(min(lower)), abs(max(upper))))
                                  xmax<- 1
                                  xmin<- -xmax
                                  
                                  oor.min<- which(lower<xmin)
                                  lower[oor.min]<- xmin
                                  oor.max<- which(upper>xmax)
                                  upper[oor.max]<- xmax
                                  
                                  
                                  par(mfrow=c(1,1))
                                  par(mar=c(4, 10, 2, 20), xpd=TRUE)
                                  plot(1:n,1:n, type="n", xlim=c(xmin,xmax), ylim=c(1,n), yaxt="n", xaxt="n",xlab="Beta", cex.lab=1.5, font.lab=2, ylab="", bty="n", line=2)
                                  
                              # y axis
                                  
                                  #ys<- seq(1, n,length.out=n)
                                  #ys<- seq(1,5,length.out=n)
                                  ys<- seq(1,5,length.out=n)
                                  ys<- c(ys, ys + 5)
                                  ys<- sort(ys, decreasing=T)
                                  ys.bars<- c(ys + 0.05, ys - 0.05)
                                  ys.numbers<- c(ys + 0.08, ys - 0.08)
                                  
                                  
                                  #ys<- c(ys-0.1, ys+0.3)
                                  
                              # Add ticks to x- and y-axes
                                  axis(1, seq(xmin,xmax,0.5), labels = FALSE, lwd=3)
                                  #axis(2, ys, labels = FALSE, lwd=3)
                              
                              # Add labels to axes      
                                  text(rep(xmin-0.2,n), ys[seq(1,length(ys),1)], adj = 1, labels = unique(labels), xpd = TRUE, cex = 1.5, font=2)
                                  text(seq(xmin,xmax,1), rep(0.4, (xmax*2)+1), srt = 0, adj = 0.5, labels = sprintf("%1.1f",seq(xmin,xmax,1)), xpd = TRUE, cex = 1.5, font=2)
                              
                              # Add segments  
                                  cols.bars<- rep(cols[c(6,10,5,9)], each=n)
                                  segments(lower, ys.bars, upper, ys.bars, lwd=4, col=cols.bars)
                                  
                                  tips<- setdiff(1:length(lower), oor.min)
                                  segments(lower[tips], ys.bars[tips]-0.02,lower[tips],ys.bars[tips]+0.02, lwd=4, col=cols.bars[tips])
                                  
                                  tips<- setdiff(1:length(upper), oor.max)
                                  segments(upper[tips],ys.bars[tips]-0.02,upper[tips],ys.bars[tips]+0.02, lwd=4, col=cols.bars[tips])
                                  
                                  for (i in oor.min) {segments(lower[i]+0.05, ys.bars[i]-0.05, xmin, ys.bars[i], lwd=4, col=cols.bars[i]);
                                                      segments(lower[i]+0.05, ys.bars[i]+0.05, xmin, ys.bars[i], lwd=4, col=cols.bars[i])}
                                  
                                  for (i in oor.max) {segments(upper[i]-0.05, ys.bars[i]+0.05, xmax, ys.bars[i], lwd=4, col=cols.bars[i]);
                                                      segments(upper[i]-0.05, ys.bars[i]-0.05, xmax, ys.bars[i], lwd=4, col=cols.bars[i])}
                                  
                                  points(mean, ys.bars, lwd=3, cex=0.8, col=cols.bars)
                                  segments(0, 0.5, 0, 5.2, lty=3, lwd=3)
                                  segments(0, 5.8, 0, 10.2, lty=3, lwd=3)
                                  
                                  
                              # Text to the right of the plot
                                  text(rep(xmax+0.2,23), c(ys.numbers,n+.6), labels=c(sprintf("%1.2f",mean),"Beta"),cex=0.9, font=1)
                                  text(rep(xmax+0.5,41), c(ys.numbers,n+.3), labels=c(sprintf("%1.2f", results.for.Forest.plot$LowerCI[o]), "Lower"),cex=0.9, font=1)
                                  text(rep(xmax+0.7,41), c(ys.numbers,n+.285), labels=c(sprintf("%1.2f", results.for.Forest.plot$UpperCI[o]), "Upper"),cex=0.9, font=1)
                                  text(xmax+0.6, max(ys.numbers)+.515, labels="Credible Interval",cex=0.9, font=1)
                                  
                                  text(rep(xmax+1,23), c(ys.numbers,n+.6), labels=c(paste(sprintf("%1.1f", pvals)),"P"),cex=0.9, font=1)
                                  text(0, 10.4, labels="PDE5i", cex=1.5, font=2, col=cols[6])
                                  text(0, 5.4, labels="ERA", cex=1.5, font=2, col=cols[10])
                                  
                                  
                              
                              # Export as 12" x 11" .pdf portrait
                              
                        
                        
          # Find orders of treatment regimens
                        DrugOrder<- matrix(0, nrow=4, ncol=4, dimnames=list(DrugConversion$DRUG[1:4], c("1st","2nd","3rd","4th")))
                        
                        for (i in 1:nrow(d2))
                        {
                          a<- unique(na.omit(unlist(d2[i,w3.col:w4.col])))
                          DrugOrderRow<- match(a, DrugConversion$INFOFLEX_CODE)
                          for (j in 1:length(a)) {DrugOrder[DrugOrderRow[j], j]<- DrugOrder[DrugOrderRow[j], j] + 1}
                          
                        }
          
                   
        
          DrugOrder.melted<- melt(DrugOrder)
          colnames(DrugOrder.melted)<- c("Drug","Order","Percent")
          DrugOrder.melted$Percent<- DrugOrder.melted$Percent
          cols<- c("white", colorRampPalette(c("white","lightblue","blue","darkblue"))(100)[10:100])
          
          ggplot(data = DrugOrder.melted, aes(x=Drug, y=Order, fill=Percent)) + 
            geom_tile(color = "white") +
            scale_fill_gradientn(colours=cols, limit = c(0,100), space = "Lab", name="Percent") +
            theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
          # Survival analysis
                            
                      # Functional class groups by PDE5i+ and PDE5i-
                            df<- data.frame(AD = as.factor(c(rep("Alive",4), rep("Dead",4))), Count=c(0,21,117,27,0,2,13,5))
                            model<- clm(AD ~ Count, data=df)         
                            anova(model, type = "II")
                            
                      # Survival prediction by functional class
                            fit<- lapply(1:imputations, function(x) {coxph(Surv(tstart, tstop, Outcome==1) ~ as.numeric(FC), data=subset(Ix.at.diagnosis, Imp==x))})
                            cox.model<- summary(pool(fit))
                            exp(cox.model[2])
                            exp(cox.model[2] + 1.96*cox.model[3])
                            exp(cox.model[2] - 1.96*cox.model[3])
                            
          
                            
                      # Fit CPH model for multiple sets of imputed data (mitools package, imputationList & MIcombine functions)
                            
                            # Demonstrate that CPH model doesn't fit these data adequately (i.e. needs time-dependent covariates)
                                fit.coxph<- with(results.olps.i.amelia$imputations[[1]],
                                                 coxph(Surv(tstart, tstop, Outcome==1) ~ AgeAtDiag + TAPSE + TreatedPDE5))
                                cox.zph(fit.coxph) # GLOBAL p<<0.001    
                            
                          
                      # TIMECOX model: frequentist, time-dependent coefficients, time-dependent covariates
                                
                                    imputations<- 1
                                    fit<- lapply(1:imputations, function(x) {
                                      cat(x);
                                      with(subset(results.olps.i.amelia$imputations[[x]]),
                                           timecox(Surv(tstart, tstop, Outcome==1) ~ const(Gender) +
                                                     const(AgeAtDiag) + TAPSE + DoseEqPDE5, n.sim=50,
                                                   id = IDn, degree=1))})
                                
                                    par(mfrow=c(2,2))
                                    plot(fit[[1]])
                                    summary(fit[[1]])
                                    coef(fit[[1]])
                                    
                                    #https://arxiv.org/pdf/1509.03253.pdf
                                    # Obtain estimates and variances for the PH covariates
                                        covarPH = fit[[1]]$gamma
                                        covarPH.var = fit[[1]]$var.gamma
                                        
                                    # Extract estimates and variances for the CUMULATIVE NON-proportional hazards (NPH) covariates
                                    #     These are all in the 'cum' value of the timecox function
                                    # Col1  = timepoints at which the parameters are estimated
                                    # Col2  = cumulative estimate for alpha0
                                    # Col3+ = cumulative estimate for the beta(t)
                                        #
                                        # HAZARD function:
                                        #
                                        # Lambda(t) = Y(t) . Lambda0(t) . exp{t[X].Beta(t) + t[Z].Gamma}
                                        #     Lambda(t) = Cumulative Hazard function at time t
                                        #     Y(t) = scaling factor
                                        #     Lambda0(t) = Baseline hazard function (time-varying)
                                        #     t[X] = transpose of X, a vector of covariate values relating to the time varying variables
                                        #     Beta(t) = matrix of time-varying coefficients for each covariate
                                        #     t[Z] = transpose of Z, a vector of covariate values relating to the time-constant variables (i.e. CPH)
                                        #     Gamma = vector of (time-constant) coefficients for each covariate
                                        
                                        
                                        # Find the baseline hazard rate: Lambda0(t)
                                              # Convert cumulative log baseline hazard -> within-bin log baseline hazard (Alpha0t)
                                                  predicted.timepts<- seq(0,2,0.01)
                                                  decumulated.ests = CsmoothB(fit[[1]]$cum, predicted.timepts, b = 1)
                                                  Alpha0t<- decumulated.ests[,2]
                                                  
                                              # Convert log baseline hazard rate -> baseline hazard rate
                                                  Lambda0t = exp(Alpha0t)
                                                  
                                              # Baseline hazard rate (lambda0(t))
                                                  par(mfrow=c(1,1)); plot(predicted.timepts, Lambda0t, type='l', lwd=4, col="red")
                                                  
                                        # Get estimates and variances for the cumulative NPH covariates
                                                  covarNPH = decumulated.ests[,-c(1:2)]
                                                  
                                                  
                                              # Multiply by the covariates for that patient (Xt * B(t))
                                                  NPH.covariates<- matrix(c(54,1),nrow=1,ncol=2)
                                                  
                                                      swe<- apply(sweep(covarNPH, 2, NPH.covariates, '*'),1,sum) # Age=54, Rx=1
                                                      swe.smoothed<- lowess(swe, delta=0.1)$y
                                                      covarNPH.Lambda<- exp(swe.smoothed)
                                                  
                                                  plot(covarNPH.Lambda)
                                                      
                                        # Cumulative log of the PH covariate effects
                                            ZtGamma<- unlist(c(1.5 * covarPH)) # (Zt * Gamma) TAPSE=1.5
                                        
                                            
                                        # Put together to make the hazard function
                                            LambdaT = Lambda0t * exp(covarNPH.Lambda + ZtGamma)
                                            S<- Area<- rep(0, length(LambdaT))
                                            for (i in 1:length(LambdaT))
                                            {
                                                Area[i]<- trapz(predicted.timepts[1:i], LambdaT[1:i])
                                                S[i]<- exp(-Area[i])
                                            }
                                            
                                            plot(1:length(S), S, type='l', ylim=c(0,1))
                                            
                                         
                                            
                                            
                                            
                                               
                                        # Smooth
                                            Surv.Func.smoothed<- matrix(0, nrow=nrow(Surv.Func), ncol=3)
                                            for (i in 1:ncol(Surv.Func)) {
                                              Surv.Func.smoothed[,i]<- lowess(Surv.Func[,i], delta=0.1)$y}
                                            
                                            
                                        # Scale-up to initial survival of 1
                                            Yt = 1/(1-Surv.Func.smoothed[1,1])
                                            
                                        # Find median survival time  
                                            fit[[1]]$cum[,1][which.min(abs(0.5 - Yt*(1-Surv.Func.smoothed[,2])))]
                                            
                                        # Make into a nice dataframe
                                            
                                            df<- data.frame(time=fit[[1]]$cum[,1],
                                                            Survival = Surv.Func.smoothed[,1],
                                                            Upper = Yt * (1 - Surv.Func.smoothed[,2]),
                                                            Lower = Yt * (1 - Surv.Func.smoothed[,3]))
                                            plot(fit[[1]]$cum[,1], Surv.Func.smoothed[,1], type='l', lwd=4, col="red", xlim=c(0,8))
                                            points(fit[[1]]$cum[,1], Yt * (1-Surv.Func.smoothed[,2]), type='l', lwd=4, col="blue")
                                            points(fit[[1]]$cum[,1], Yt * (1-Surv.Func.smoothed[,3]), type='l', lwd=4, col="green")
                                            p<- ggplot(data=df, aes(x=time, y=Survival)) + geom_line() + geom_smooth()
                                            p<- p + geom_ribbon(aes(ymin=df$Lower, ymax=df$Upper), linetype=2, alpha=0.1)
                                            p                
                            
                                            
                            
                
                                            
                                            
                                            
                                            
                                            
                                                  
                            survival.results<- summary(fit.p)[,1:4]
                            survival.results$p.values<- pnorm(-abs(summary(fit.p)$results / summary(fit.p)$se))
                            r<- round(survival.results,3)
                            
                            interpret.CPH.interactions(r)
                            
                            
                            
                            # Draw some survival curves to show the above results
                                  cox.combined<- coxph(Surv(tstart, tstop, Outcome == 1) ~ tt(TreatedPDE5) +tt(TAPSE), data=results.olps.i.amelia$imputations[[1]])
                                  cox.combined$coefficients<- fit.p$coefficients
                                  cox.combined$var<- fit.p$variance
                                  
                                  new.data<- data.frame(expand.grid(TreatedPDE5 = c(0,1), TAPSE = c(0.5 * mean(Ix.at.diagnosis$TAPSE),mean(Ix.at.diagnosis$TAPSE))))
                                  s<- survfit(cox.combined, newdata=new.data)
                                  
                            # Fit CPH with covariates
                                  cox.combined<- coxph(Surv(tstart, tstop, Outcome == 1) ~ AgeAtDiag, data=results.olps.i.amelia$imputations[[1]])
                            
                      # Add a p-value based on a normal distribution
                            coxph.pvalue<- pnorm(fit.p$coefficients / sqrt(fit.p$variance), mean=0, sd=1)
                            if (coxph.pvalue<0.001) {coxph.pvalue="<0.001"} else {coxph.pvalue=paste("=",coxph.pvalue,sep="")}
                            
                      # Find the median survival and 95% CIs for PDE5i+ and PDE5i- groups in each imputed set
                            surv_median(survfit(fit[[1]], newdata = data.frame(TreatedPDE5 = c(1,0))))
                            median.survivals.list<- lapply(1:imputations, function(x) {surv_median(survfit(fit[[x]], newdata = data.frame(TreatedPDE5 = c(1,0))))})
                            median.survivals.matrix<- matrix(round(as.numeric(unlist(median.survivals.list)),2),nrow=8, dimnames=list(c("Untreated","Treated","Treated_Median","Untreated_Median","Treated_Lower_95CI","Untreated_Lower_95CI","Treated_Upper_95CI","Untreated_Upper_95CI"), paste("Imp",1:imputations,sep="")))
                            median.survivals<- apply(median.survivals.matrix,1,mean)
                            
                      # Make a table of the results
                            T1<- rbind(T1, c("Outcome", "Survival (years)", median.survivals, coxph.pvalue))
                  
          
          
          
          # Save the data
                #setwd("~/Dropbox/ILD-PH/ILD PDE5 Paper")        
                #save(list=ls(all=TRUE), file="SetUp.RData")         
          
                                
                              
                              # https://rpkgs.datanovia.com/survminer/
                              # http://www.drizopoulos.com/courses/emc/ep03_%20survival%20analysis%20in%20r%20companion
                              
                              pdf(file="Figures/SupFigure2.pdf")
                              p<- ggsurvplot(
                                s,
                                data = results.olps.i.amelia$imputations[[1]],
                                size = 1,                 # change line size
                                palette =
                                  c("darkred", "darkblue","red","blue"),# custom color palettes
                                conf.int = TRUE,          # Add confidence interval
                                pval = TRUE,              # Add p-value
                                risk.table = FALSE,        # Add risk table
                                risk.table.col = "strata",# Risk table color by groups
                                legend.labs = c("PDE5i- / low TAPSE", "PDE5i+ / low TAPSE","PDE5i- / high TAPSE","PDE5i+ / high TAPSE"), 
                                censor = FALSE,
                                surv.median.line = "hv",
                                xlim = c(0,10),
                                break.x.by = 1,
                                xlab = "Time (years)",
                                #fun = "cumhaz",
                                legend = "top", # Change legend labels
                                ggtheme = theme_classic()      # Change ggplot2 theme
                              )
                              p
                              # Find median survival time
                                  size<- 5
                                  for (i in 1:ncol(s$surv))
                                    {
                                      median<- format(round(s$time[which.min(abs(0.5 - s$surv[,i]))],2),nsmall=2)
                                      y.height<- 0.1*match(summary(s)$table[,7], sort(summary(s)$table[,7]))[i]
                                      p$plot<- p$plot + ggplot2::annotate("text", 
                                                             x = as.numeric(median)-0.1, y = 0.5 - y.height, # x and y coordinates of the text
                                                             label = as.character(median), size = size, )
                              
                                      p$plot<- p$plot + ggplot2::annotate("text", 
                                                              x = as.numeric(median)-0.1, y = 0.47 - y.height, # x and y coordinates of the text
                                                              label = "years", size = size)
                                    }
                              
                              
                              p
                              dev.off()
                                 
              
   

