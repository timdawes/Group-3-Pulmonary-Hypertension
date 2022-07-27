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


tiff("SFigure4.tiff", width = 5, height = 7, units='in', res=300, compression='lzw')


km_subtype<- with(d2, Surv(FUTimeYears, Dead))
d2$Diagnosis_sub_grouped<- d2$Diagnosis_sub
d2$Diagnosis_sub_grouped[which(d2$Diagnosis_sub_grouped %in% c("3.3","3.4","3.5","3.7"))]<- "Other"
d2$Diagnosis_sub_grouped<- factor(d2$Diagnosis_sub_grouped, levels = c("3.1", "3.2", "Other"))
table(d2$Diagnosis_sub_grouped)


km_group<- survfit(Surv(FUTimeYears, Dead) ~ Diagnosis_sub_grouped, data = d2)
p<- ggsurvplot(km_group, data = d2, size = 2,                # change line size
               palette = c(cols[c(4,6,2)]),#, add.alpha(cols[3],0.99)),
               conf.int = TRUE,          # Add confidence interval
               conf.int.alpha = 0.4,
               pval = TRUE,              # Add p-value
               pval.size = 6,
               pval.coord = c(0,0.1),
               risk.table = TRUE,        # Add risk table
               risk.table.col = "strata",# Risk table color by groups
               legend.labs = rev(c("Groups 3.3-3.7","Group 3.2","Group 3.1")), 
               censor = FALSE,
               xlim = c(0,10),
               break.x.by = 1,
               break.y.by = 0.2,
               xlab = "Time (years)",
               legend = c(0.85,0.9),
               legend.title = "",
               font.x = c(18, "plain"),
               font.y = c(18, "plain"),
               font.tickslab = c(14, "plain"),
               risk.table.height = 0.35,
               risk.table.y.text = F,
               risk.table.font.x = c(32, "plain"))
p

dev.off()                                  


surv_median(km_group)

# 1,2 and 3 year survival by subgroup
      surv_pvalue(survfit(Surv(FUTimeYears, Dead) ~ Diagnosis_sub_grouped, data = d2))
      summary(survfit(Surv(FUTimeYears, Dead) ~ 1, data = d2), times=c(1,2,3,4,5))

