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
#



# Figure 3

# BSE Criteria for RV dysfunction
        d2$RVdys<- unlist(sapply(1:nrow(d2), function(x) {d2[x, which(colnames(d2) == paste("RVdys.Echo", d2$Nearest.Echo.To.Dx.Number[x],sep=""))]}))

           # RV dysfunction by FAC criterion
                    d2$FAClow<- 0
                    d2$FAClow[which(d2$Gender==0)]<- unlist(sapply(which(d2$Gender==0), function(x) {d2[x, which(colnames(d2) == paste("FAC.Echo", d2$Nearest.Echo.To.PDE.Number[x],sep=""))] < 35}))
                    d2$FAClow[which(d2$Gender==1)]<- unlist(sapply(which(d2$Gender==1), function(x) {d2[x, which(colnames(d2) == paste("FAC.Echo", d2$Nearest.Echo.To.PDE.Number[x],sep=""))] < 30}))

           # RV dysfunction by TAPSE criterion
                    d2$TAPSElow<- sapply(1:nrow(d2), function(x) {d2[x, which(colnames(d2) == paste("TAPSE.Echo", d2$Nearest.Echo.To.PDE.Number[x],sep=""))] < 1.7})
                    d2$TAPSElow[is.na(d2$TAPSElow) == TRUE]<- FALSE

           # RV dysfunction by S' criterion
                    d2$Slow<- sapply(1:nrow(d2), function(x) {d2[x, which(colnames(d2) == paste("S.Echo", d2$Nearest.Echo.To.PDE.Number[x],sep=""))] < 9})
                    d2$Slow[is.na(d2$Slow) == TRUE]<- FALSE


           # Any of the above
                    RVdys.IDn<- which(d2$RVdys==TRUE | d2$TAPSElow==TRUE | d2$FAClow==TRUE | d2$Slow==TRUE)
                    d2$RVdys[RVdys.IDn]<- 1
                                      
           # ..or those without
                    noRVdys.IDn<- setdiff(1:nrow(d2), RVdys.IDn)
                    d2$RVdys[noRVdys.IDn]<- 0
                                      
                                      
                                      
                                      IDn<- All.IDn
                                      IDn<- RVdys.IDn
                                      IDn<- noRVdys.IDn
                                      
                        # Choose appropriate colours
                                      if (identical(IDn, All.IDn)) {km.cols<- c(brewer.pal(9,"Greys")[c(8)], brewer.pal(9,"Blues")[c(8)])}
                                      if (identical(IDn, RVdys.IDn)) {km.cols<- c(brewer.pal(9,"Greys")[c(8)], brewer.pal(9,"Reds")[c(8)])}
                                      if (identical(IDn, noRVdys.IDn)) {km.cols<- c(brewer.pal(9,"Greys")[c(8)], c(brewer.pal(9,"Greens")[c(8)]))}
                    
            # Plot the Kaplan Meier plot
                        km_PDE5<- survfit(Surv(FUTimeYears, Outcome) ~ PDE, data = d2[IDn,])
                        p2<- ggsurvplot(km_PDE5, data = d2[IDn,], size = 1.5,                 # change line size
                                        palette = km.cols,# custom color palettes
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
                                        gg_theme = theme_classic()
                        ) 
                        p2
                        
                        pdf(file="Figure2 Rx.pdf", width=6, height=6, onefile = FALSE)
                        #pdf(file="Figure3 noRVDys.pdf", width=6, height=6, onefile = FALSE)
                        
                        p2
                        dev.off()
                        
                        surv_pvalue(km_PDE5)
                        km.ci(km_PDE5, method="peto")
                        
                        
