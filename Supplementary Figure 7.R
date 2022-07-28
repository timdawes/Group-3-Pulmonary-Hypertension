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


# Take the non-imputed dataset
      
      d5<- results.olps.non.i
      d5$IDn<- match(d5$ID, unique(d5$ID))
      d5$Diagnosis_IPF<- as.factor(d5$ID %in% d2$URN[which(d2$Diagnosis_IPF==1)])
      

# For fields with zero values, find a patient who matches nearest by age and gender and take one value of this field

      
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
    n.iter<- 2000
    
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

      }


# Print output to scatterplot

      pdf(file="Sensitivity Analysis.pdf", width=13, height=8)

      par(mfrow=c(1,1), mar=c(6, 11, 4, 35), xpd=TRUE)
      plot(1:n.survival,1:n.survival, type="n", xlim=c(xmin,xmax), ylim=c(1,y.height), yaxt="n", xaxt="n",xlab="BETA", cex.lab=1.5, font.lab=2, ylab="", bty="n", line=3)
          
    # Add ticks to x- and y-axes
          axis(1, seq(xmin,xmax,0.5), labels = FALSE, lwd=3)
    
    # Add labels to x-axis      
          text(seq(xmin,xmax,length.out=3), rep(0.65, 3), srt = 0, adj = 0.5, labels = sprintf("%+1.1f",seq(xmin,xmax,length.out=3)), xpd = TRUE, cex = 1.5, font=2)
    
    # Y-AXIS
    # Define the y-coordinates
          ys.bars<- seq(1,y.height,length.out=n.survival)
          z<- length(na.omit(coefs.survival[,1]))
          y.counter<- (1:z)
          
    # Add segments  
          marker.size<- 0.05
          
        
                    n.priors<- nrow(prior.list)
                    priors.to.use<- order(prior.list$mu)
                    marker.size<- 0.02
                    lwd.priors<- 3
                    
                    cols<- c(colorRampPalette(c("white",brewer.pal(9,"Greens")[c(8)]))(6)[-1],
                             colorRampPalette(c("white",brewer.pal(9,"Reds")[c(8)]))(6)[-1],
                             colorRampPalette(c("white",brewer.pal(9,"Blues")[c(8)]))(6)[-1])
                    
                    
                             
                        for (i in priors.to.use)
                          { 
                            # Bars
                                cols.bars<- rep(cols[i], each=n.survival)
                                ys.bars.priors<- outer(ys.bars, (seq(-0.5,0.5,length.out=n.priors) * (ys.bars[2] - ys.bars[1]) * 0.5)[i], '+')
                                
                                coefs.priors<- summary(coxAI.sens[[i]])[[6]]$'Surv(futime, Outcome)'$regcoef[o,c("Mean","2.5%","97.5%","tail-prob.")]
                                lower<- coefs.priors[,"2.5%"]
                                mean<- coefs.priors[,"Mean"]
                                upper<-  coefs.priors[,"97.5%"]
                                
                                segments(lower, ys.bars.priors, upper, ys.bars.priors, lwd=lwd.priors, col=cols[i])
                            
                            # LEFT end markers
                                segments(lower, ys.bars.priors-marker.size, lower, ys.bars.priors+marker.size, lwd=lwd.priors, col=cols[i])
                                
                            # RIGHT end markers
                                segments(upper, ys.bars.priors-marker.size, upper, ys.bars.priors+marker.size, lwd=lwd.priors, col=cols[i])
                                
                            # OOR arrows
                                #marker.size<- 0.05
                                #for (i in oor.min) {segments(lower[i]+marker.size, ys.bars.priors[i]+marker.size, xmin, ys.bars.priors[i], lwd=lwd, col=cols[i]);
                                #  segments(lower[i]+marker.size, ys.bars.priors[i]-marker.size, xmin, ys.bars[i], lwd=lwd, col=cols[i])}
                                
                                #for (i in oor.max) {segments(upper[i]-0.05, ys.bars.priors[i]+marker.size, xmax, ys.bars.priors[i], lwd=lwd, col=cols[i]);
                                #  segments(upper[i]-0.05, ys.bars.priors[i]-marker.size, xmax, ys.bars[i], lwd=lwd, col=cols[i])}
                                
                            # Add point spheres
                            points(mean, ys.bars.priors, cex=1, col=cols[i], pch=20)
                            
                            # Posterior samples
                            for (j in 1:n.chains) # was 'imputations but the .pdf is massive
                            {
                              xs<- coxAI.sens[[i]]$MCMC[[j]][tail(1:n.iter, n.iter/5),o]
                              xs[xs>xmax]<- NA
                              xs[xs<xmin]<- NA
                              random.y<- runif(nrow(xs),0,0.001) - 0.00005
                              
                              for (k in 1:z)
                              {
                                ys<- rep(ys.bars.priors[k],nrow(xs)) + random.y
                                points(xs[,k], ys, col=alpha(cols[i], 0.2), pch=20, cex=0.4)
                              }
                            }
                        }
                    
                    for (k in 1:5) {
                      a<- 0.15+(k/8)
                      segments(graph.inc,a,graph.inc+col.inc*1.5,a,lwd=20, col=cols[priors.to.use[k]], lend=2);
                      segments(graph.inc+col.inc*2,a,graph.inc+col.inc*3.5,a,lwd=20, col=cols[priors.to.use[k]+5], lend=2);
                      segments(graph.inc+col.inc*4,a,graph.inc+col.inc*5.5,a,lwd=20, col=cols[priors.to.use[k]+10], lend=2);
                      #segments(graph.inc+col.inc*6,a,graph.inc+col.inc*7.5,a,lwd=20, col=cols[priors.to.use+15], lend=2)
                    }
                    
                    for (k in 1:5) {
                      a<- 0.15+(k/8)
                      text(graph.inc+col.inc*0.75,a,labels = paste("tau=",prior.list[priors.to.use[k],1], "  mu=",prior.list[priors.to.use[k],2],sep=""), adj=0.5)
                      text(graph.inc+col.inc*2.75,a,labels = paste("tau=",prior.list[priors.to.use[k+5],1], "  mu=",prior.list[priors.to.use[k+5],2],sep=""), adj=0.5)
                      text(graph.inc+col.inc*4.75,a,labels = paste("tau=",prior.list[priors.to.use[k+10],1], "  mu=",prior.list[priors.to.use[k+10],2],sep=""), adj=0.5);
                      #text(graph.inc+col.inc*6.75,a,labels = paste("tau=",prior.list[priors.to.use+15,1], "  mu=",prior.list[priors.to.use[k+15],2],sep=""), adj=0.5)
                    }
                       
            
    # Add vertical lines at 'zero'
          segments(0, .8, 0, y.height+0.2, lty=2, lwd=4)
    
    # Add the text around the plot
          text.cex<- 1.4
          lab.cex<- 1.4
          col.inc<- 0.8
          graph.inc<- xmax + 0.5
          top.header.inc<- 0.52
          
    # Better/Worse Survival
          #text(rep(xmax,n.survival)+0.1, ys.bars, adj = 0, labels = ns.labels.survival.worse[o], xpd = TRUE, cex = text.cex, font=2)
          text(rep(xmin,n.survival)-0.1, ys.bars, adj = 1, labels = ns.labels.survival.better[o], xpd = TRUE, cex = text.cex, font=2)
          
    # Left of graph
          text(xmin-0.1, y.height+top.header.inc, labels="BETTER",cex=text.cex, font=2, adj=1)
          text(xmin-0.1, y.height+0.35, labels="SURVIVAL",cex=text.cex, font=2, adj=1)
                   
    # Beta
          text(rep(graph.inc,20), c(ys.bars,y.height+top.header.inc), labels=c(sprintf("%+1.2f",mean),"BETA"),xpd = TRUE, cex=text.cex, font=2)
    
    # Credible intervals for beta
          text(rep(graph.inc+col.inc,length(ys.bars)+1), c(ys.bars,y.height+0.35), labels=c(sprintf("%+1.2f", c(coefs.survival[,2])), "LOWER"),cex=lab.cex, font=2, adj=0.5)
          text(rep(graph.inc+col.inc*2,length(ys.bars)+1), c(ys.bars,y.height+0.35), labels=c(sprintf("%+1.2f", c(coefs.survival[,3])), "UPPER"),cex=lab.cex, font=2, adj=0.5)
          text(graph.inc+col.inc*1.5, max(ys.bars)+top.header.inc, labels="CI",cex=lab.cex, font=2, adj=0.5)
    
    # Effect size & Credible Intervals
          rounding.digits<- rep(1,n.survival)
          effect.size<- sapply(mean.cols, function(x) {sprintf("%1.2f", exp(mean)[x])})
          text(rep(graph.inc+col.inc*3,12), c(ys.bars,y.height+0.35, y.height+top.header.inc), labels=c(effect.size,"RATIO","HAZARD"), cex=lab.cex, font=2, adj=0.5)
          
          
    # Credible Intervals
          effect.size.lower<- sapply(mean.cols, function(x) {sprintf("%1.2f", exp(lower)[x])})
          text(rep(graph.inc+col.inc*4,11), c(ys.bars,y.height+0.35), labels=c(effect.size.lower, "LOWER"),cex=lab.cex, font=2, adj=0.5)
          effect.size.upper<- sapply(mean.cols, function(x) {sprintf("%1.2f", exp(upper)[x])})
          text(rep(graph.inc+col.inc*5,11), c(ys.bars,y.height+0.35), labels=c(effect.size.upper, "UPPER"),cex=lab.cex, font=2, adj=0.5)
          
          text(graph.inc+col.inc*4.5, max(ys.bars)+top.header.inc, labels="CI",cex=lab.cex, font=2, adj=0.5)
          text(rep(graph.inc+col.inc*6,12), c(ys.bars,y.height+0.35, y.height+top.header.inc), labels=c(pvals,"","Pr(>0)"),cex=lab.cex, font=2, adj=0.5)
          
          dev.off()
    
    
