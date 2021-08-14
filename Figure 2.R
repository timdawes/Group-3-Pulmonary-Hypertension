# Pulmonary vasodilator treatment and survival in group 3 pulmonary hypertension: an observational study
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
# Copyright Tim Dawes, August 2021
#

# Figure 2 (left panel): KM plot by PDE5+ / PDE5-  

# Choose a subgroup
      IDn<- 1:nrow(d2)

# Draw the plot
      km<- with(d2[IDn,], Surv(FUTimeYears, Dead))
      km_PDE5<- survfit(Surv(FUTimeYears, Dead) ~ PDE, data = d2[IDn,])
      p2<- ggsurvplot(km_PDE5, data = d2[IDn,], size = 1.5,                 # change line size
                      palette = cols[c(2,6)],# custom color palettes
                      conf.int = TRUE,          # Add confidence interval
                      conf.int.alpha = 0.2,
                      pval = TRUE,              # Add p-value
                      risk.table = TRUE,        # Add risk table
                      risk.table.col = "strata",# Risk table color by groups
                      censor = FALSE,
                      surv.median.line = "hv",
                      xlim = c(0,10),
                      break.x.by = 1,
                      break.y.by = 0.1,
                      xlab = "Time (years)",
                      legend = "none", # Change legend labels
                      ggtheme = theme_classic()      # Change ggplot2 theme
            )
      dev.off()
      p2

# Output the plot
      pdf(file="Figure3 All PDE5i Rx.pdf", width=8, height=5, onefile = FALSE)
      print(
        p2,
        surv.plot.height = 1,
        risk.table.height = .25,
        ncensor.plot.height = NULL,
        newpage = TRUE,
      )
      dev.off()

      # Find Kaplan Meier confidence intervals
            with(d2[IDn,], round(surv.median(Surv(FUTimeYears, Dead), c(0.025, 0.5, 0.0975), 500),2))
            km.ci(km_PDE5, method="peto")

            df<- data.frame(surv=summary(km_PDE5)$upper[1:16], time=summary(km_PDE5)$time[1:16])
            fit<- lm(log(surv) ~ time, data=df)
            seq(1,10,0.1)[which.min(exp(predict(fit, data.frame(time=seq(1,10,0.1)))) - 0.5)]


            
            
# Figure 2 (centre panel): KM plot by ERA+ / ERA-
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
                     legend = "none", # Change legend labels
                     ggtheme = theme_classic()      # Change ggplot2 theme
      )
      surv_median(km_ERA)
      
      
# Output the plot
      pdf(file="Figure3 All ERA Rx.pdf", width=8, height=5, onefile = FALSE)
      print(
        p,
        surv.plot.height = 1,
        risk.table.height = .25,
        ncensor.plot.height = NULL,
        newpage = TRUE,
      )
      dev.off()
      
      
# Find Kaplan Meier confidence intervals
      km.ci(km_ERA, method="peto")


      
      

# Figure 2 (right panel): KM plot by number of therapies
      km<- with(d2[IDn,], Surv(FUTimeYears, Dead))
      km_NoTherapies<- survfit(Surv(FUTimeYears, Dead) ~ NoTherapies, data = d2[IDn,])
      p<- ggsurvplot(km_NoTherapies, data = d2[IDn,], size = 1.5,                 # change line size
                     palette = brewer.pal(11,"Blues")[c(4,7,9)],# custom color palettes
                     conf.int = TRUE,          # Add confidence interval
                     conf.int.alpha = 0.4,
                     pval = TRUE,              # Add p-value
                     risk.table = TRUE,        # Add risk table
                     risk.table.col = "strata",# Risk table color by groups
                     legend.labs = c("None", "One","Two"), 
                     censor = FALSE,
                     surv.median.line = "hv",
                     xlim = c(0,10),
                     break.x.by = 1,
                     break.y.by = 0.1,
                     xlab = "Time (years)",
                     legend = "none", # Change legend labels
                     ggtheme = theme_classic()      # Change ggplot2 theme
      )
      
      round(surv.median(Surv(d2[NoTherapies == 1,]$FUTimeYears, d2[NoTherapies == 1,]$Dead), c(0.025, 0.5, 0.0975), 500),2)
      surv_median(km_NoTherapies)
      

      
# Output the plot
      pdf(file="Figure2 All NoTherapies.pdf", width=5, height=5, onefile = FALSE)
      print(
        p,
        surv.plot.height = 1,
        risk.table.height = .25,
        ncensor.plot.height = NULL,
        newpage = TRUE,
      )
      dev.off()      
