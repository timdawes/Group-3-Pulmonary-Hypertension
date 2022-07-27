# Pulmonary vasodilator treatment and survival in group 3 pulmonary hypertension: a Bayesian observational cohort study
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
#
#
# Survival analysis code



# Take the non-imputed dataset
      
      d5<- results.olps.non.i
      d5$IDn<- match(d5$ID, unique(d5$ID))
      d5$Diagnosis_IPF<- as.factor(d5$ID %in% d2$URN[which(d2$Diagnosis_IPF==1)])
      

# For fields with zero values, find a patient who matches nearest by age and gender
# and take one value of this field

      
      for (j in c("AgeAtDiag","BNP","Gender","PVR","TAPSE","TLcoperc","distExTest"))
      {
        col<- match(j, colnames(d5))
        
            for (i in unique(d5$IDn))
              {
              rows<- which(d5$IDn==i)
              d5$label<- 0
              d5$label[rows]<- 1
              
                
                    if (all(is.na(d5[rows,col]))) {
                      cat("\n",unique(d5$ID[rows]),j)
                      form<- as.formula(c("label ~ AgeAtDiag + Gender"))
                      d6<- d5[c(rows[1], which(is.na(d5[,col])==FALSE)),]
                      m<- matchit(form, method="nearest", data=d6)$match.matrix[,1]
                      d5[rows[1],col]<- d6[match(m, rownames(d6)),col]
                    }
                
              d5$Outcome[rows]<- max(d5$Outcome[rows])
              }
      }

      
# Re-format outcome and age
      d5$Outcome<- as.factor(d5$Outcome)
      d5$AgeAtDiag<- as.numeric(d5$AgeAtDiag)
      d5$Gender<- as.factor(d5$Gender)
      
    n.chains<- 3
    n.iter<- 2000

    
# Scale these columns selectively
        data.sc<- d5
        scaling.data<- rep(0, 18)
        names(scaling.data)<-  c("DoseEqPDE5_Mean","DoseEqPDE5_sd", "AgeAtDiag_Mean","AgeAtDiag_sd", "TAPSE_Mean","TAPSE_sd", "FVCperc_Mean","FVCperc_sd", 
                                    "PVR_Mean","PVR_sd", "BNP_Mean","BNP_sd", "TLcoperc_Mean","TLcoperc_sd","BMI_Mean","BMI_sd","6MWD_Mean","6MWD_sd")
        cols.sc<- match(c("DoseEqPDE5","AgeAtDiag", "TAPSE","FVCperc","PVR","BNP","TLcoperc","BMI","distExTest"), colnames(data.sc))
        data.sc$BNP<- log(data.sc$BNP)
        data.sc$PVR<- log(data.sc$PVR)
        data.sc$BMI<- log(data.sc$BMI)
    
    for (j in 1:length(cols.sc)){
      from<- 2*(j-1) + 1
      to<- from + 1
      scaling.data[from:to]<- c(mean(na.omit(data.sc[,cols.sc[j]])), sd(na.omit(data.sc[,cols.sc[j]])))
      data.sc[,cols.sc[j]]<- (data.sc[,cols.sc[j]] - scaling.data[from]) / scaling.data[to]}

# Run joint models for survival
    
    prior.list<- data.frame(mu=rep(c(-2,0,2),each=5), tau=rep(10^seq(-3,1,1),3))
    sensitivity.analysis<- TRUE
    n.iter<- 200
    
        if (sensitivity.analysis==TRUE) {
                    coxAI.sens<- list()
                    for (model in 1:nrow(prior.list))
                    {
                      cat(model,"...",sep="")
                      
                      hyp<- default_hyperpars()
                      hyp$norm["mu_reg_norm"]<- prior.list[model,"mu"]
                      hyp$norm["tau_reg_norm"]<- prior.list[model,"tau"]
                      
                      coxAI.sens[[model]]<- coxph_imp(Surv(futime, Outcome) ~ DoseEqPDE5*TAPSE + AgeAtDiag + BNP + FVCperc + Gender + PVR + TLcoperc + distExTest + (1 | ID),
                                                data = data.sc, 
                                                df_basehaz = 6,
                                                n.chains = n.chains,
                                                n.adapt = 200,
                                                thin = 1,
                                                monitor_params = c(imps = TRUE),
                                                timevar = 'tstart',
                                                n.iter = n.iter,
                                                hyperpars = hyp)
                      
                      if (model>1) {boxplot(t(sapply(1:model, function(x) {exp(summary(coxAI.sens[[x]])[[6]]$'Surv(futime, Outcome)'$regcoef[,1])})))}
                    }
          
        } else {
          
                  
                      coxAI.RVdys<- coxph_imp(Surv(futime, Outcome) ~ DoseEqPDE5*TAPSE + AgeAtDiag + BNP + FVCperc + Gender + PVR + TLcoperc + distExTest + (1 | ID),
                                        data = data.sc, 
                                        df_basehaz = 6,
                                        n.chains = n.chains,
                                        n.adapt = 200,
                                        thin = 1,
                                        monitor_params = c(imps = TRUE),
                                        timevar = 'tstart',
                                        n.iter = 1000)
                      
                  
                  cbind(round(exp(summary(coxAI.RVdys)[[6]]$'Surv(futime, Outcome)'$regcoef[,c(1,3,4)]),2), p = round(summary(coxAI.RVdys)[[6]]$'Surv(futime, Outcome)'$regcoef[,5],3))
                  
                  
          
                  }
                
    
    # Frequentist Cox analysis for the Supplementary Data
    
      
      imputations<- 100
      plot(0,0,xlim=c(1,imputations), ylim=c(0,1), type='n')
      coxFreq.RVdys<- list()
      
      for (k in 1:imputations) {
        
        Ix.at.diagnosis<- NULL
        cat(k)
        
        for (j in 1:128) {
        a<- which(results.olps.i.amelia$imputation[[1]]$ID == d2[j,]$URN)
        b<- order(abs(results.olps.i.amelia$imputations[[1]][a,]$DateOfScan - ymd(d2[j,]$dateDIAGNOSIS)), decreasing=FALSE)[1]
        Ix.at.diagnosis<- rbind(Ix.at.diagnosis, cbind(d2[j,]$URN, futime = d2[j,]$FUTimeYears, Outcome = d2[j,]$Outcome, Imp = k, results.olps.i.amelia$imputations[[k]][a[b],]))
        }
        
        Ix.at.diagnosis.sc<- cbind(Ix.at.diagnosis[,c("TreatedPDE5","Gender","futime","Outcome")], scale(Ix.at.diagnosis[,c("AgeAtDiag","BNP","distExTest","FVCperc","PVR","TAPSE","TLcoperc")]))
        
        coxFreq.RVdys[[k]]<- coxph(Surv(futime, Outcome) ~ TreatedPDE5 + AgeAtDiag + BNP + FVCperc + Gender + PVR + TAPSE + TLcoperc + distExTest, data = Ix.at.diagnosis.sc)
        points(k, median(sapply(1:k, function(x) {summary(coxFreq.RVdys[[x]])$coefficients[1,]})[5,]), col="blue", pch=19)
      }
      
      
      
    matrix.mean<- t(sapply(1:imputations, function(x) {summary(coxFreq.RVdys[[x]])$coefficients[,1]}))
    matrix.se<- t(sapply(1:imputations, function(x) {summary(coxFreq.RVdys[[x]])$coefficients[,3]}))
    dim(matrix.mean)
    dim(matrix.se)
    
    ncoefs.survival.coefficients<- 9
    names<- rownames(summary(coxFreq.RVdys[[1]])$coefficients)
    
    results.pooling.coxPH<- pooling.by.Rubins.rules(matrix.mean, matrix.se, imputations, ncoefs.survival.coefficients)

    rownames(results.pooling.coxPH)<- names
    

    
