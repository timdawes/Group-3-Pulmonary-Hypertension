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


   # Choose appropriate colours
                        km.cols<- c(brewer.pal(9,"Reds")[8], brewer.pal(9,"Blues")[8],
                                    brewer.pal(9,"Greens")[8], brewer.pal(9, "Purples")[8],
                                    brewer.pal(9,"Greys")[8])
                        
                        diagnosis.cat<- d2$Diagnosis_cat4
                        diagnosis_cat<- model.matrix(~ diagnosis.cat)
                        diagnosis_cat[,1]<- rep(0, nrow(diagnosis_cat))
                        SubgroupOne<- which(apply(diagnosis_cat,1,sum)==0)
                        diagnosis_cat[SubgroupOne,1]<- 1
                        colnames(diagnosis_cat)<- c("Undifferentiated","Other","NSIP","IPF","HP")
                        d2.plus<- cbind(d2, diagnosis_cat)
                        
                        
                        km_Diag<- survfit(Surv(FUTimeYears, Outcome) ~ HP + IPF + NSIP + Other + Undifferentiated, data = d2.plus)

                        p<- ggsurvplot(km_Diag, data = d2.plus, size = 1.5,                # change line size
                          palette = km.cols,
                          conf.int = TRUE,          # Add confidence interval
                          conf.int.alpha = 0.2,
                          pval = TRUE,
                          risk.table = TRUE,        # Add risk table
                          risk.table.col = "strata",# Risk table color by groups
                          legend.labs = c("HP","IPF","NSIP","Undifferentiated","Other"), 
                          xlim = c(0,5),
                          surv.median.line = "hv",
                          break.x.by = 1,
                          break.y.by = 0.1,
                          xlab = "Time (years)",
                          legend = "none",
                          censor = TRUE,
                          gg_theme = theme_classic()
                          )
                        
                      
                        pdf("SFigure5.pdf", width = 8, height = 8, onefile = FALSE)
                        p
                        dev.off()            
                   
