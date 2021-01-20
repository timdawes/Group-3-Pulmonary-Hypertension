# Group 3 PH Project
# Tim Dawes September 2020
# Bayesian CPH method

# Survreg to check distributions

      t<- results.olps.i.amelia$imputations[[1]]
      surv.obj<- with(t, Surv(tstart, tstop, Outcome == 1))





      
# So now fit the Bayesian Cox model with covariate transformation as decided by the above
        # MCMC parameters
      
      imputations<- 100
      res.bayes<- list()
      nburn=1000; nsave=1000; nskip=10
      mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000)
      
      # Setting priors
      
          # Fit a gamma distribution to the FU times and extract the two coefficients which define it
                estimates<- fitdistr(na.omit(d$FUTimeYears[-which(d$FUTimeYears<=0)]), "gamma")$estimate
          # Define the priors using the gamma coefficients and theta0 from the baseline survival function
                prior = list(maxL=15, a0=estimates[1], b0=estimates[2],
                             theta0 = fit_lognormal$icoef)
                
                scaling.data<- matrix(0, nrow=imputations, ncol=6, 
                                      dimnames=list(paste("I",1:imputations,sep=""), 
                                                    c("DoseEqPDE5_Mean","DoseEqPDE5_sd",
                                                      "AgeAtDiag_Mean","AgeAtDiag_sd",
                                                      "TAPSE_Mean","TAPSE_sd")))
                
      # Choose a specific subgroup?
            COPD.IDn<- which(d2$Diagnosis_sub=="3.1")
            IPF.IDn<- which(d2$Diagnosis_sub=="3.2")
            Other.IDn<- setdiff(1:nrow(d2), union(COPD.IDn, IPF.IDn))
            All.IDn<- 1:nrow(d2)
            
      # Choose a subgroup and a label to remind you!
            IDn<- COPD.IDn
            ID.label<- "COPD"
  
            
                   
      for (i in 1:imputations)
              {
                        cat("\n \n ",i)
                  # Pick out the repeat observations of the patients in the subset of patients of interest
                              which.rows<- which(results.olps.i.amelia$imputations[[1]]$IDn %in% IDn)
                        # Screen text to confirm which cohort of patients you're currently working on
                              cat("\n Diagnoses in cohort = ", names(table(results.olps.i.amelia$imputations[[1]]$Diagnosis_sub[which.rows])), "\n \n")
                        # Select out the subset of patients
                              which.cols<- match(c("DoseEqERA","AgeAtDiag", "TAPSE","Gender","tstop","tstart","Outcome","IDn","ALT"), colnames(results.olps.i.amelia$imputations[[i]]))
                              data.sc<- results.olps.i.amelia$imputations[[i]][which.rows,which.cols]
                              
                        
                  # Scale and centre the data in the columns of interest
                        # Pick out the columns in the analysis, so they can be scaled
                              cols.sc<- match(c("DoseEqERA","AgeAtDiag", "TAPSE"), colnames(data.sc))
                              
                        # Scale these columns selectively
                              for (j in 1:length(cols.sc)){
                                    from<- 2*(j-1) + 1
                                    to<- from + 1
                                    scaling.data[i,from:to]<- c(mean(data.sc[,cols.sc[j]]), sd(data.sc[,cols.sc[j]]))
                                    data.sc[,cols.sc[j]]<- (data.sc[,cols.sc[j]] - scaling.data[i,from]) / scaling.data[i,to]}
                                    #data.sc<- data.sc + rnorm(nrow(data.sc)*ncol(data.sc), 0, .1)
                              
                              
                  # Bayesian survival regression
                              options(warn=-1)
                              repeat {
                                      t<- tryCatch(survregbayes(Surv(tstart, tstop, Outcome == 1) ~ DoseEqERA + AgeAtDiag + as.factor(Gender) + TAPSE,
                                                                    data = data.sc,
                                                                    prior = prior, mcmc=mcmc, survmodel = "AFT", # AFT because KM curves diverge and then converge
                                                                    dist = "weibull", subject.num = IDn))
                                      
                                      fail<- FALSE
                                      if(("try-error" %in% class(t))==TRUE) {fail<- TRUE}
                                      if (abs(summary(t)[[4]][1,1])>5) {fail<- TRUE}
                                      if (fail == FALSE) {res.bayes[[i]]<- t; break}
                              
                                      cat("\n Repeating... \n")
                              }
                              options(warn=0)
                              # Plot the data so far as reality check
                              matrix.mean<- lapply(1:i, function(x) {summary(res.bayes[[x]])[[4]]})
                              matrix.mean.plot<- matrix(unlist(matrix.mean), ncol=20, byrow=T)
                              plot(1:imputations, seq(-2,2,length.out=imputations), type='n', ylab="Coefficient")
                              for (j in 1:i) {points(rep(j, 4), matrix.mean.plot[j,1:4], type='p', pch=19, col=c("blue","red","green","black"), cex=0.2)}
                              
                              
                              
                              
                }
        
            
                
            # Pool coefficients   
                matrix.mean<- lapply(1:imputations, function(x) {summary(res.bayes[[x]])[[4]]})
                matrix.sd<- lapply(1:imputations, function(x) {summary(res.bayes[[x]])[[4]]})
                coef.df<- pooling.by.Rubins.rules(matrix.mean, matrix.sd, imputations)                  
        
                
                        
            # Plot survival curves
                  res<- res.bayes[[100]]
                  a<- matrix(0, nrow=4, ncol=1000)
                  for (i in 1:4) {a[i,]<- rnorm(1000,mean=coef.df$mean[i], sd=coef.df$SD[i])}  
                  res$beta<- a
                  par(mfrow=c(1,1));
                  wide=0.01;
                  tgrid = seq(1e-10,10,wide);
                  ngrid = length(tgrid);
                  
                  newdata = data.frame(DoseEqERA=c(0,1), AgeAtDiag=c(0,0), TAPSE=c(0,0), Gender=c(0,0), tstart=0)
                  matching.cols.1<- na.omit(match(paste(colnames(newdata),"_Mean",sep=""), colnames(scaling.data)))
                  matching.cols.2<- match(colnames(scaling.data)[matching.cols.1], paste(colnames(newdata),"_Mean",sep=""))
                  
                  
                  newdata.c<- sweep(newdata[,matching.cols.2], 2, apply(scaling.data[,matching.cols.1],2,mean), "-") # Newdata centred
                  newdata.cs<- cbind(sweep(newdata.c, 2, apply(scaling.data[,matching.cols.1+1],2,mean), "/"), newdata[,4:5]) # Newdata centred and scaled
                  
                  
                  p<- plot(res, xnewdata=newdata, tgrid=tgrid, PLOT=TRUE) # Shat = survival, hhat = hazard, fhat = density
                  res.pooled.df<- data.frame(mean = 100 * c(p$Shat[,1], p$Shat[,2]),
                                             upper = 100 * c(p$Shatup[,1], p$Shatup[,2]),
                                             lower = 100 * c(p$Shatlow[,1], p$Shatlow[,2]),
                                             curve = c(rep("a",ngrid), rep("b",ngrid)),
                                             colour = c(rep("blue",ngrid), rep("red", ngrid)),
                                             interval = rep(tgrid, times=2))
                  
                  median.surv.no.Rx<- round(tgrid[which.min(abs(p$Shat[,1]-0.5))],2)
                  median.surv.Rx<- round(tgrid[which.min(abs(p$Shat[,2]-0.5))],2)
                  
                  res.pooled.df %>%
                    ggplot(aes(interval, mean, fill=curve)) + 
                    geom_ribbon(aes(ymin = lower,
                                    ymax = upper,
                                    group = curve),
                                fill = "grey", alpha=0.4) +
                    geom_line(aes(color=colour), size=2) +
                    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                       legend.position = "none", panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
                    xlab("Time (years)") + ylab("Survival (%)") +
                    scale_x_continuous(breaks=seq(0,10,2), labels=seq(0,10,2)) +
                    scale_y_continuous(breaks=seq(0,100,10), labels=seq(0,100,10)) +
                    geom_segment(aes(x=0, y=50, xend=median.surv.Rx, yend=50), linetype="dashed", size=1) +
                    geom_segment(aes(x=median.surv.no.Rx, y=0, xend = median.surv.no.Rx, yend=50), size=1, linetype="dashed") +
                    geom_segment(aes(x=median.surv.Rx, y=0, xend = median.surv.Rx, yend=50), size=1, linetype="dashed") +
                    geom_text(x=median.surv.no.Rx, y=-2, label=paste(as.character(median.surv.no.Rx),"years")) +
                    geom_text(x=median.surv.Rx-0.5, y=-2, label=paste(as.character(median.surv.Rx),"years"))
                    
                  
# Calculate difference in median survival for different age and TAPSE and display it
                  
        # Decide the age and TAPSE ranges
            age.seq<- seq(20,80,5)
            TAPSE.seq<- seq(0.5,2.5,0.25)
        
        # Make a long format matrix to store the results
            l<- length(TAPSE.seq)*length(age.seq)
            median.surv.mats<- matrix(0, nrow=l, ncol=2)

        # Make a wide format matrix to store all this                  
            e<- expand.grid(age.seq, TAPSE.seq)
            
            # Calculate fitted values for women
                  newdata.women = data.frame(DoseEqPDE5=rep(c(0,1), each=l),
                                 AgeAtDiag=rep(e$Var1, 2), 
                                 Gender=as.factor(rep(0, l*2)),
                                 TAPSE=rep(e$Var2, 2),
                                 tstart=0)
                  p.women<- plot(res, xnewdata=newdata.women, tgrid=tgrid, PLOT=FALSE) # Shat = survival, hhat = hazard, fhat = density
                  for (j in 1:l) {median.surv.mats[j,]<- c(tgrid[which.min(abs(p.women$Shat[,j]-0.5))], tgrid[which.min(abs(p.women$Shat[,j+l]-0.5))])}
                  
            # Calculate the fitted values for men
                  newdata.men = data.frame(DoseEqPDE5=rep(c(0,1), each=l),
                                             AgeAtDiag=rep(e$Var1, 2), 
                                             Gender=as.factor(rep(1, l*2)),
                                             TAPSE=rep(e$Var2, 2),
                                             tstart=0)
                  p.men<- plot(res, xnewdata=newdata.men, tgrid=tgrid, PLOT=FALSE) # Shat = survival, hhat = hazard, fhat = density
                  for (j in 1:l) {median.surv.mats[j,]<- c(tgrid[which.min(abs(p.men$Shat[,j]-0.5))], tgrid[which.min(abs(p.men$Shat[,j+l]-0.5))])}
                  
        # Create vectors containing coordinates for interpolation
            interp_coords <- expand_grid(age = seq(20, 80,0.2), TAPSE = seq(0.5, 2, 0.02))
                  
        # Interpolate using the interp() from the pracma library
            mat<- matrix(c(median.surv.mats[,2] - median.surv.mats[,1]), nrow=length(TAPSE.seq), byrow=T, dimnames=list(TAPSE.seq, age.seq))
            interp_vals <- pracma::interp2(x = age.seq, y = TAPSE.seq, Z = mat, xp = interp_coords$age, yp = interp_coords$TAPSE)
                  
        # Bind in new df, convert to matrix and plot
                  cols<- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
                  df_2 <- cbind(interp_coords, interp_vals)
                  mat2 <- reshape2::acast(df_2, age~TAPSE, value.var="interp_vals" )
        
        # Plot
                  plot.new()
                  par(mar=c(4,4,1,1))
                  layout(matrix(c(1,2), nrow=2), heights=c(3,1))
                  image(mat2, axes=F, col = cols(1000), add=F)
                  age.seq.for.axis<- seq(min(age.seq), max(age.seq), 10)
                  TAPSE.seq.for.axis<- seq(min(TAPSE.seq), max(TAPSE.seq), length.out=5)
                  
                  axis(1, at=seq(0,1, length.out=length(age.seq.for.axis)), labels=age.seq.for.axis, line=0.5)
                  axis(2, at=seq(0,1, length.out=5), labels=seq(0.5, 2.5, 0.5), line=0.5)
                  title(xlab="Age (years)", ylab="TAPSE (cm)")
                  
                  par(mar=c(5,4,2,1))
                  scale<- seq(1, 8, length.out=8)
                  image.scale(mat2, col=cols(length(scale)-1), breaks=scale, horiz=TRUE)
                  title(xlab="Difference in survival (years)")
                       
        # Save as .pdf 8" x 6"     
              
              
