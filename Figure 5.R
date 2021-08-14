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
#
# Survival analysis code
# Copyright Tim Dawes, August 2021



                  
# Figure 5: Predicted difference in survival for women and men
# Calculate difference in median survival for different age and TAPSE and display it
                  
        # Decide the age and TAPSE ranges
            age.seq<- TAPSE.seq<- seq(-3,3,1)
            
        # Make a long format matrix to store the results
            l<- length(TAPSE.seq)*length(age.seq)
            median.surv.mats.women<- median.surv.mats.men<- matrix(0, nrow=l, ncol=2)

        # Make a wide format matrix to store all this                  
            e<- expand.grid(age.seq, TAPSE.seq)
            
            # Calculate fitted values for women
                  newdata.women.cs = data.frame(DoseEqPDE5=rep(c(0,1), each=l), DoseEqERA=rep(c(0,0), each=l), mPAP=c(0,0), 
                                             AgeAtDiag=rep(e$Var1, 2), PVR=c(0,0), FVC=c(0,0), TAPSE=rep(e$Var2, 2))
                  
                  # Centre and scale the data
                      matching.cols.1<- na.omit(match(paste(colnames(newdata.women.cs),"_Mean",sep=""), colnames(scaling.data)))
                      matching.cols.2<- match(colnames(scaling.data)[matching.cols.1], paste(colnames(newdata.women.cs),"_Mean",sep=""))
                      newdata.women.c<- sweep(newdata.women.cs, 2, apply(scaling.data[,matching.cols.1+1],2,mean), "*")
                      newdata.women<- sweep(newdata.women.c, 2, apply(scaling.data[,matching.cols.1],2,mean), "+")
                      newdata.women.cs<- cbind(newdata.women.cs, SubgroupOne=0, SubgroupTwo=1, Gender=1, tstart=0)
                      
                      
                  p.women<- plot(res, xnewdata=newdata.women.cs, tgrid=tgrid, PLOT=TRUE) # Shat = survival, hhat = hazard, fhat = density
                  for (j in 1:l) {median.surv.mats.women[j,]<- c(tgrid[which.min(abs(p.women$Shat[,j]-0.5))], tgrid[which.min(abs(p.women$Shat[,j+l]-0.5))])}
                  
            # Calculate the fitted values for men
                  newdata.men.cs = data.frame(DoseEqPDE5=rep(c(0,1), each=l), DoseEqERA=rep(c(0,0), each=l), mPAP=c(0,0),
                                          AgeAtDiag=rep(e$Var1, 2), PVR=c(0,0), FVC=c(0,0), TAPSE=rep(e$Var2, 2))
                                           
                  # Centre and scale the data
                      matching.cols.1<- na.omit(match(paste(colnames(newdata.men.cs),"_Mean",sep=""), colnames(scaling.data))) # Find the columns in newdata.men which have centring data
                      matching.cols.2<- match(colnames(scaling.data)[matching.cols.1], paste(colnames(newdata.men.cs),"_Mean",sep="")) # Find the columns in the scaling data with the necessary mean and sd
                      newdata.men.c<- sweep(newdata.men.cs, 2, apply(scaling.data[,matching.cols.1+1],2,mean), "*") 
                      newdata.men<- sweep(newdata.men.c, 2, apply(scaling.data[,matching.cols.1],2,mean), "+")
                      newdata.men.cs<- cbind(newdata.men.cs, SubgroupOne=0, SubgroupTwo=1, Gender=0, tstart=0)
                      
                  p.men<- plot(res, xnewdata=newdata.men.cs, tgrid=tgrid, PLOT=TRUE) # Shat = survival, hhat = hazard, fhat = density
                  for (j in 1:l) {median.surv.mats.men[j,]<- c(tgrid[which.min(abs(p.men$Shat[,j]-0.5))], tgrid[which.min(abs(p.men$Shat[,j+l]-0.5))])}
                  
        # Create vectors containing coordinates for interpolation
                  
            #interp_coords <- expand_grid(age = seq(20,80,.1), TAPSE = seq(0.5, 2, 0.1))
            interp_coords <- expand_grid(age = (seq(30,80,length.out=100) - mean(scaling.data[,5])) / mean(scaling.data[,6]),
                                         TAPSE = (seq(0.5,3,length.out=100) - mean(scaling.data[,7])) / mean(scaling.data[,8]))
            
                  
        # Interpolate using the interp() from the pracma library for women
            cols<- rev(blue2green2red(1000))
            
            mat.women<- matrix(c(median.surv.mats.women[,2] - median.surv.mats.women[,1]), nrow=length(TAPSE.seq), byrow=T, dimnames=list(TAPSE.seq, age.seq))
            interp_vals <- interp2(x = age.seq, y = TAPSE.seq, Z = mat.women, xp = interp_coords$age, yp = interp_coords$TAPSE)
            df_2 <- cbind(interp_coords, interp_vals)
            mat2.women <- reshape2::acast(df_2, age~TAPSE, value.var="interp_vals" )
            
        # Interpolate using the interp() from the pracma library for men
            mat.men<- matrix(c(median.surv.mats.men[,2] - median.surv.mats.men[,1]), nrow=length(TAPSE.seq), byrow=T, dimnames=list(TAPSE.seq, age.seq))
            interp_vals <- interp2(x = age.seq, y = TAPSE.seq, Z = mat.men, xp = interp_coords$age, yp = interp_coords$TAPSE)
            df_2 <- cbind(interp_coords, interp_vals)
            mat2.men <- reshape2::acast(df_2, age~TAPSE, value.var="interp_vals" )
            
            
        # Plot for women
                  layout(matrix(c(1,2), ncol=2, byrow=T))
                  z.range<- c(0.3,5.9)
                  plot.new()
                  image(mat2.women, zlim=z.range, axes=F, col = cols, add=T)
                  contour(mat2.women, add=T, levels=seq(0,5,0.5), lwd=4, labcex=2)
                  
                  age.seq.for.axis<- seq(30,80,10)
                  TAPSE.seq.for.axis<- seq(0.5,3,length.out=6)
                  
                  axis(1, at=seq(0,1, length.out=length(age.seq.for.axis)), labels=age.seq.for.axis, line=-1, lwd=4, cex.axis=2)
                  axis(2, at=seq(0,1, length.out=6), labels=TAPSE.seq.for.axis, line=-0.9, lwd=4, cex.axis=2, las=2)
                  title(xlab="Age (years)", ylab="TAPSE (cm)", line=+2.5, cex.lab=2)
                  
        
        # Plot for men
                  plot.new()
                  image(mat2.men, zlim=z.range, axes=F, col = cols, add=T)
                  contour(mat2.men, add=T, levels=seq(0,5,0.5), lwd=4, labcex=2)
                  
                  age.seq.for.axis<- seq(30, 80, 10)
                  TAPSE.seq.for.axis<- seq(0.5, 3, length.out=6)
                  
                  axis(1, at=seq(0,1, length.out=length(age.seq.for.axis)), labels=age.seq.for.axis, line=-1,lwd=4, cex.axis=2)
                  axis(2, at=seq(0,1, length.out=6), labels=TAPSE.seq.for.axis, line=-0.9, lwd=4, cex.axis=2, las=2)
                  title(xlab="Age (years)", ylab="TAPSE (cm)", line=+2.5, cex.lab=2)
                  
                  plot.new()
                  image.scale(z=mat2.men)
                  image.scale(z=mat2.men, zlim=z.range, col=rev(blue2green2red(10000)), breaks=seq(floor(z.range[1]), ceiling(z.range[2]),length=10001), axis.pos=4, add.axis=TRUE)
                  dev.off()
                  
                  
                       
                    
                    image.scale <- function(z, zlim, col = heat.colors(12),
                                            breaks, axis.pos=1, add.axis=TRUE, ...){
                      if(!missing(breaks)){
                        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
                      }
                      if(missing(breaks) & !missing(zlim)){
                        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
                      }
                      if(missing(breaks) & missing(zlim)){
                        zlim <- range(z, na.rm=TRUE)
                        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
                        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
                        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
                      }
                      poly <- vector(mode="list", length(col))
                      for(i in seq(poly)){
                        poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
                      }
                      if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
                      if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
                      plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")  
                      for(i in seq(poly)){
                        if(axis.pos %in% c(1,3)){
                          polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
                        }
                        if(axis.pos %in% c(2,4)){
                          polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
                        }
                      }
                      
                      axis(at=seq(floor(z.range[1]), ceiling(z.range[2]), 1), labels=TRUE, cex.axis=2, cex.ticks=2, side=axis.pos, las=2)
                    }
                  
                    
                    
            