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

    # Choose a specific subgroup?
                                  cats<- names(table(d2$Diagnosis_cat4))
                                  km.cols<- brewer.pal(10,"Paired")[c(2,4,6,8,10)]
                                  
                                  
                                  for (i in 1:length(cats))
                                  {
                                        IDn<- which(d2$Diagnosis_cat4 != cats[i])
                                        
                                        km_PDE5<- survfit(Surv(FUTimeYears, Outcome) ~ PDE, data = d2[IDn,])
                                        p2<- ggsurvplot(km_PDE5, data = d2[IDn,], size = 1.5,                 # change line size
                                                        palette = km.cols,
                                                        conf.int = TRUE,          # Add confidence interval
                                                        conf.int.alpha = 0.2,
                                                        risk.table = TRUE,        # Add risk table
                                                        risk.table.col = "strata",# Risk table color by groups
                                                        legend.labs = c("PDE5i-","PDE5i+"), 
                                                        surv.median.line = "hv",
                                                        xlim = c(0,5),
                                                        break.x.by = 1,
                                                        break.y.by = 0.1,
                                                        xlab = "Time (years)",
                                                        legend = "none", # Change legend labels
                                                        gg_theme = theme_classic())
                                                        
                                        surv_pvalue(km_PDE5)
                                        km.ci(km_PDE5, method="peto")
                                                        
                                                        
                                        ) 
                                  
                                        title<- paste(cats[i],".pdf", sep="")
                                        
                                        pdf(file=title, width=6, height=6, onefile = FALSE)
                                        p2
                                        dev.off()
                                        
                                  }      
                                  
                         
