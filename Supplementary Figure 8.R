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


# Define a new variable called 'TreatedAndImproved' which will define whether patients were treated and improved
# Takes 1 of 4 values:
# 1. Insufficient Data to make a decision (default)
# 2. Patient wasn't treated with PDE5i
# 3. Patient was treated with PDE5i but didn't improve
# 4. Patient was treated with PDE5i and improved as judged by BNP, PVR, TAPSE, 6MWD

                d2$TreatedAndImproved<- "InsufficientData"
                d2$TreatedAndImproved[which(d2$PDE==0)]<- "NotTreated"

# Define improvement.marker, which is a new matrix which will store the details of all the cases which showed improvement
                improvement.marker<- NULL

# Define FC1 and FC2 which are two new variables which will store the Functional Class variables before and after 'improvement' as defined above
                d2$FC1<- d2$FC2<- NA
                
                
                
    # Loop through each formula analysing it
                  
                  
        for (i in 1:length(fmla.mm))
        {
                if (i==1) {cat("\n Association between treatment and functional markers: ")}
                cat(" ",ns[i])
              
                # Untreated patients
                      fit_me_noRx[[i]]<- lme_imp(fmla.mm.noRx[[i]],
                                             data = data.sc.func2[-treated.ids,],
                                             random = ~ 1 | IDn,
                                             family = gaussian(link = "identity"),
                                             n.chains = 3, n.adapt = 500, n.iter = 2000, thin = 1)
                      
                # Treated patients before treatment starts
                      fit_me_preRx[[i]]<- lme_imp(fmla.mm.noRx[[i]],
                                             data = data.sc.func2[t.beforeRx,],
                                             random = ~ 1 | IDn,
                                             family = gaussian(link = "identity"),
                                             n.chains = 3, n.adapt = 500, n.iter = 2000, thin = 1)
                      
                # Treated patients after treatment starts
                      fit_me_postRx[[i]]<- lme_imp(fmla.mm[[i]],
                                             data = data.sc.func2[t.afterRx,],
                                             random = ~ 1 | IDn,
                                             family = gaussian(link = "identity"),
                                             n.chains = 3, n.adapt = 500, n.iter = 2000, thin = 1)
                      
                # Find subjects whose trajectory improves after treatment
                      
                          unique.treated.ids<- unique(data.sc.func2$ID[treated.ids])
                          
                              for (id in unique.treated.ids)
                              {
                                dataset<- d5[t.beforeRx,]
                                y1<- na.omit(dataset[which(dataset$ID==id),ns[i]])
                                
                                dataset<- d5[t.afterRx,]
                                y2<- na.omit(dataset[which(dataset$ID==id),ns[i]])
                              
                                change.y2<- NA
                                if (length(y1)!=0 && length(y2)!=0) {change.y2<- head(na.omit(y2),1) - tail(na.omit(y1),1)}
                                
                                if (ns[i] == "distExTest") {d2$ChangeIndistExTest[match(id, d2$URN)]<- change.y2}
                                if (ns[i] == "TAPSE") {d2$ChangeInTAPSE[match(id, d2$URN)]<- change.y2}
                                if (ns[i] == "FC") 
                                  {
                                  if (length(y1)!=0) {d2$FC1[which(d2$URN == id)]<- tail(na.omit(y1),1)}
                                  if (length(y2)!=0) {d2$FC2[which(d2$URN == id)]<- head(na.omit(y2),1)}
                                  }
                                          
                                  
                                if (ns[i] == "BNP") {d2$ChangeInBNP[match(id, d2$URN)]<- change.y2}
                                
                                
                                
                                 # Binary idea of whether a parameter improves with treatment
                                       better<- NA
                                       
                                       if (is.na(change.y2)==FALSE) {
                                           if (change.y2 < 0 && ns[i]=="BNP") {better<- TRUE; improvement.marker<- rbind(improvement.marker, c(id, change.y2, d2[which(d2$URN == id),c("Rxfutime","Outcome")], "BNP"))}
                                           if (change.y2 < 0 && ns[i]=="PVR") {better<- TRUE; improvement.marker<- rbind(improvement.marker, c(id,change.y2, d2[which(d2$URN == id), c("Rxfutime","Outcome")], "PVR"))}
                                           if (change.y2 > 0 && ns[i]=="TAPSE") {better<- TRUE; improvement.marker<- rbind(improvement.marker, c(id,change.y2, d2[which(d2$URN == id), c("Rxfutime","Outcome")], "TAPSE"))}
                                           if (change.y2 > 0 && ns[i]=="distExTest") {better<- TRUE; improvement.marker<- rbind(improvement.marker, c(id, change.y2, d2[which(d2$URN == id), c("Rxfutime","Outcome")], "distExTest"))}
                                         
                                         
                                           if (change.y2 >= 0 && ns[i]=="BNP") {better<- FALSE}
                                           if (change.y2 >= 0 && ns[i]=="PVR") {better<- FALSE}
                                           if (change.y2 <= 0 && ns[i]=="TAPSE") {better<- FALSE}
                                           if (change.y2 <= 0 && ns[i]=="distExTest") {better<- FALSE}
                                         
                                         }
                                
                                       if (is.na(better)==FALSE) {if (better == FALSE) {d2[which(d2$URN == id),]$TreatedAndImproved<- "TreatedAndDidNotImprove"}}
                                
                                       
                                     
                                             
                                rm(dataset); rm(y1); rm(y2); rm(x1); rm(x2); rm(better)
                                
                                
                                
                              }
                      }
                  
                
                # Supplementary Figure 8: Responders and Non-responders
                
                      d2$TreatedAndImproved[match(unlist(improvement.marker[,1]), d2$URN)]<- "TreatedAndImproved"
                      
                      Rx<- which(d2$TreatedAndImproved %in% c("NotTreated", "TreatedAndDidNotImprove","TreatedAndImproved"))
                      table(d2$TreatedAndImproved)
                      
                      km_PDE5<- survfit(Surv(Rxfutime, Outcome) ~ TreatedAndImproved, data = d2[Rx,])
                      km.cols<- c(brewer.pal(9,"Greens")[c(8)], brewer.pal(9,"Reds")[8], brewer.pal(9,"Blues")[8])
                      
                      
                      titles<- c("Not Treated","Treated, no improvement", "Treated, improvement")
                      
                      p2<- ggsurvplot(km_PDE5, data = d2[Rx,], size = 1.5,                 # change line size
                                      palette = km.cols,
                                      conf.int = TRUE,          # Add confidence interval
                                      conf.int.alpha = 0.4,
                                      pval = FALSE,              # Add p-value
                                      risk.table = TRUE,        # Add risk table
                                      risk.table.col = "strata",# Risk table color by groups
                                      legend.labs = titles,
                                      censor = FALSE,
                                      surv.median.line = "hv",
                                      xlim = c(0,5),
                                      break.x.by = 1,
                                      break.y.by = 0.1,
                                      xlab = "Time (years)",
                                      legend = "none", # Change legend labels
                                      ggtheme = theme_classic()      # Change ggplot2 theme
                      )
                      p2
                      surv_pvalue(km_PDE5)
                      km.ci(km_PDE5, method="peto")
                      
                      pdf("SFigure8.pdf", width = 8, height = 8)
                      print(
                        p2,
                        surv.plot.height = 1,
                        risk.table.height = .20,
                        ncensor.plot.height = NULL,
                        newpage = FALSE,
                      )
                      dev.off()
                    


                # Difference in Functional Classes before/after treatment in patients who DID respond
                      before.DID<- d2$FC1[match(unique(unlist(improvement.marker[,1])), d2$URN)]
                      after.DID<- d2$FC2[match(unique(unlist(improvement.marker[,1])), d2$URN)]
                      length(before.DID)
                      
                # Difference in Functional Classes before/after treatment in patients who DID NOT respond
                      before.DIDNOT<- d2$FC1[which(d2$URN %in% setdiff(unique.treated.ids, unlist(improvement.marker[,1])))]
                      after.DIDNOT<- d2$FC2[which(d2$URN %in% setdiff(unique.treated.ids, unlist(improvement.marker[,1])))]
                      length(before.DIDNOT)      
                      
                      
                      time<- c(rep(0, length(before.DID)), rep(1, length(after.DID)),
                               rep(0, length(before.DIDNOT)), rep(1, length(after.DIDNOT)))
                      length(time)
                      
                      Rx<- c(rep(1, each = (length(before.DID) + length(after.DID))), rep(0, (length(before.DIDNOT) + length(after.DIDNOT))))
                      length(Rx)
                      
                      df<- data.frame(index = c(rep(1:length(before.DID),2), rep(1:length(before.DIDNOT), 2)),
                                      variable = rep("FC", 2 * (length(before.DID) + length(before.DIDNOT))),
                                      FC = c(before.DID, after.DID, before.DIDNOT, after.DIDNOT),
                                      time = time,
                                      Rx = Rx)
                      
                      fit.FC<- lm(FC ~ time*Rx, data=df[which(df$FC!=0),])
                      
                      summary(fit.FC)
                      confint(fit.FC)
                      par(mfrow=c(1,1), mar = c(4,4,4,4))
                      barplot(FC ~ time + Rx, data = df)
                      median(before.DID)
                      
                      # Stacked + percent
                      df<- data.frame(Proportion = c(table(factor(before.DID, 1:4)), table(factor(after.DID, 1:4)), table(factor(before.DIDNOT, 1:4)), table(factor(after.DIDNOT, 1:4))),
                                      FC = factor(rep(1:4, times = 4)),
                                      Time = factor(rep(rep(c("Before Treatment","After Treatment"), each=4),2), levels = c("Before Treatment","After Treatment")),
                                      Rx = rep(c("Responders","Non-Responders"), each = 8))
                      
                      df<- df[-which(df$FC == 1),]
                      names(df)[names(df)=="FC"]<- "Functional Class"
                      
                      pdf(file="Figure S9.pdf", width=6, height=6)
                      
                      ggplot(df, aes(fill=FC, y=Proportion, x=Time)) + 
                        geom_bar(position="fill", stat="identity") + facet_grid( ~ Rx) + theme_bw() +
                        scale_fill_brewer(c(brewer.pal(9,"Reds")[c(8)], brewer.pal(9,"Blues")[c(8)], brewer.pal(9, "Greens")[c(8)])) +
                        theme(legend.position="bottom") +
                        scale_fill_discrete(name = "Functional Class") +
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
                        
                      
                      dev.off()
                
    # Fig

