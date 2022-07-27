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
# Bayesian mixed linear model code
# Copyright Tim Dawes, October 2021


# Links to Figure 5 (Association of PDE5 with echocardiographic, haemodynamic and functional markers)

  # List of markers of interest (ns) and their labels( ns.labels) and a count of how many their are (n)
                              ns<- c("distExTest","FC","EmPHasis10","TAPSE","HR","CO","PVR","mPAP","BNP","Kcoperc","SaO2","SvO2") 
                              ns.labels<- c("Six-minute walk","Functional class","emPHasis-10","Tricuspid annulus planar","Heart rate (bpm)","Cardiac output (L/min)","Pulmonary vascular","Mean pulmonary artery","B-type natriuretic","Kco (/s)","Arterial saturation (%)","Venous saturation (%)")
                              ns.labels2<- c("distance (m)","","","systolic excursion (cm)","","","resistance (WU)","pressure (mmHg)","peptide (ng/L)","","","")
                              n<- length(ns)
                              
                        # Construct formulae to analyses these markers
                              fmla.mm<- paste(ns, c("~ "), c("DoseEqPDE5 +  RxTime + DoseEqPDE5*RxTime + (1 | IDn)"))
                              fmla.mm.MCMC.fixed<- paste(ns, c("~ "), c("DoseEqPDE5 +  RxTime + DoseEqPDE5*RxTime"))
                              
                        # Make matrices to store the analysis results
                              func.cols<- c("Intercept","DoseEqPDE5","RxTime","DoseEqPDE5:RxTime")
                              func.ncols<- length(func.cols)
                              
                              func.colnames<- paste(rep(ns,each=12), "_", rep(func.cols, times=n),sep="")
                              func.rownames<- paste("I",1:imputations,sep="")
                              
                              scaling.colnames<- paste(rep(ns,each=2), c("_Mean","_sd"),sep="")
                              scaling.rownames<- ns
                              fit_me2<- list()
                              scaling.data<- matrix(0, nrow=length(fmla.mm), ncol=length(scaling.colnames), dimnames=list(scaling.rownames, scaling.colnames))
                              
                        # Choose a specific subgroup of patients?
                              
                            
                        # Loop through each formula analysing it
                              for (i in 1:length(fmla.mm))
                                {
                                if (i==1) {cat("\n Association between treatment and functional markers: ")}
                                cat(" ",c(ns,ns)[i])
                                fit_me2[[i]]<- list()
                                
                                for (x in 1:imputations)
                                {
                                  if ((x/10)==round(x/10)) {cat(".")}
                                            # Find the right rows (IDs) and columns (fields) and extract these data only
                                                  cols.mm<- match(c("DoseEqPDE5","AgeAtDiag", "tstart","tstop","ID", "IDn", ns), colnames(results.olps.i.amelia$imputations[[1]]))
                                                  rows.mm<- which(results.olps.i.amelia$imputations[[1]]$IDn %in% IDn)
                                                  data.sc.func<- results.olps.i.amelia$imputations[[x]][rows.mm,cols.mm]
                                                        
                                                  # All columns with 2 or fewer members = factors, >2 or more members = continuous
                                                        these.are.factors<- which(sapply(apply(data.sc.func,2,table),length)<=2)
                                                        these.are.not.factors<- setdiff(1:ncol(data.sc.func), these.are.factors)
                                                        for (k in these.are.factors) {data.sc.func[,k]<- as.factor(data.sc.func[,k])} 
                                                        for (k in these.are.not.factors) {data.sc.func[,k]<- as.numeric(data.sc.func[,k])} 
                                                        
                                                 # Align the timings so that the start of treatment = time zero for all patients       
                                                        data.sc.func2<- NULL
                                                        for (k in unique(data.sc.func$ID))
                                                        {
                                                          ages<- data.sc.func[which(data.sc.func$ID==k),]$tstart
                                                          sc.data.per.pt<- data.sc.func[which(data.sc.func$ID==k)[order(ages)],]
                                                          first.Rx.date<- min(na.omit(sc.data.per.pt[which(sc.data.per.pt$DoseEqPDE5!=0)[1],]$tstart))
                                                          if (first.Rx.date==Inf) {first.Rx.date<- min(sc.data.per.pt$tstart)}
                                                          sc.data.per.pt$RxTime<- sc.data.per.pt$tstart - first.Rx.date
                                                          data.sc.func2<- rbind(data.sc.func2,sc.data.per.pt)
                                                        }
                                                  
                                                      
                                                  # Scale and centre the data in the columns of interest
                                                        cols.sc<- match(ns, colnames(data.sc.func2))
                                                  
                                                        # Scale these columns selectively
                                                            for (j in 1:length(cols.sc))
                                                                {
                                                                  from<- 2*(j-1) + 1
                                                                  to<- from + 1
                                                                  scaling.data[i,from:to]<- c(mean(data.sc.func[,cols.sc[j]], na.rm=T), sd(data.sc.func[,cols.sc[j]], na.rm=T))
                                                                  data.sc.func2[,cols.sc[j]]<- (data.sc.func2[,cols.sc[j]] - scaling.data[i,from]) / scaling.data[i,to]}
                                                 
                                                 V1<- list(V=c(0.4, 2.2, 0.1, 0.0, 1.3, -1.1, 0.2)[match(ns[i], c("FC","distExTest","EmPHasis10","mPAP","CO","PVR","RAP"))], nu=0.002)
                                                 V2<- V3<- V4<- V5<- V6<- V7<- V8<- V9<- list(V=1, nu=0.002)
                                                               
                                            # Fit the model        
                                                            fit_me2[[i]][[x]]<- MCMCglmm(fixed = as.formula(fmla.mm.MCMC.fixed[i]),
                                                                                         random = ~ ID,
                                                                                         data = data.sc.func2,
                                                                                         #prior = list(R = list(V1=V1, V2=V2, V3=V3, V4=V4, V5=V5, V6=V6, V7=V7, V8=V8, V9=V9),
                                                                                        #              G = list(G1 = list(V = 1, nu = 0.002))),
                                                                                         nitt=2000, burnin=1000, thin=10, verbose=FALSE)
                                                            
                                            }
                                  }
                              
                        
                        
                        # Pool the results
                              
                              # Store the mean and se coefficients
                                          matrix.mean<- matrix.lowerCI<- matrix.upperCI<- matrix.sd<- coef.df.mm<- list()
                                          imputations<- 100

                                          for (model in 1:length(fmla.mm))
                                          { 
                                            cat(model)
                                            matrix.mean[[model]]<- t(sapply(1:imputations, function(x) {summary(fit_me2[[model]][[x]])$solutions[,1]}))
                                            matrix.lowerCI[[model]]<- t(sapply(1:imputations, function(x) {summary(fit_me2[[model]][[x]])$solutions[,2]}))
                                            matrix.upperCI[[model]]<- t(sapply(1:imputations, function(x) {summary(fit_me2[[model]][[x]])$solutions[,3]}))
                                            matrix.sd[[model]]<- (matrix.upperCI[[model]] - matrix.lowerCI[[model]]) / (2*1.96)
                                            coef.df.mm[[model]]<- pooling.by.Rubins.rules(matrix.mean[[model]], matrix.sd[[model]], imputations, ncol(matrix.mean[[model]]))
                                            rownames(coef.df.mm[[model]])<- colnames(matrix.mean[[model]])
                                          }
                              
                              
