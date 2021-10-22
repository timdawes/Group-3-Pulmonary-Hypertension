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
# Figure 6 code: Bayesian mixed linear model scatterplot
# Copyright Tim Dawes, October 2021



      ns.labels.mm<- c("PDE5i treatment","ERA treatment","Age","Gender","TAPSE","FVC","PVR","Subgroup 3.1","Subgroup 3.2")#,"PDE5i*Subgroup 3.1","PDE5i*Subgroup 3.2")
      baseline<- c(sapply(c("^RxTime$","DoseEqPDE5:RxTime","DoseEqERA:RxTime"), function(x) grep(x, rownames(coef.df.mm[[1]]))))

      n.mm<- length(fmla.mm)
      lower<- mean<- upper<- pvals<- rep(0, length=n.mm*3)


      for (model in 1:n.mm)
      {
        counter<- ((model-1)*3 + 1) : ((model-1)*3+3)
        lower[counter]<- coef.df.mm[[model]]$LowerCI[baseline]
        mean[counter]<- coef.df.mm[[model]]$mean[baseline]
        upper[counter]<-  coef.df.mm[[model]]$UpperCI[baseline]
        pvals[counter]<- coef.df.mm[[model]]$PerCentChance[baseline]
      }

      pvals[which(mean<0)]<- 100-pvals[which(mean<0)]
      pvals<- pvals/100
      pvals<- sprintf("%1.2f", pvals)
      pvals[pvals=="1.00"]<- ">0.99"
      pvals[pvals=="0.00"]<- "<0.01"


# Y height
      y.height<- 6.2

# X-AXIS
      xmax<- 0.3
      xmin<- -xmax
      lwd<- 3

# Change any coefficients which are out-of-range (oor) to the x-limit value
      oor.min<- which(lower<xmin)
      lower[oor.min]<- xmin
      oor.max<- which(upper>xmax)
      upper[oor.max]<- xmax


# Set-up the diagram
    tiff("Figure6.tiff", width = 14, height = 7, units='in', res=300, compression='lzw')
    lab.cex<- 0.7
    par(mfrow=c(1,1), mar=c(2, 9.5, 1.2, 20), xpd=TRUE)
    plot(1:n.mm,1:n.mm, type="n", xlim=c(xmin,xmax), ylim=c(1,y.height), yaxt="n", xaxt="n",
         xlab="BETA", cex.lab=0.7, font.lab=2, ylab="", bty="n", line=1)

# Add ticks to x- and y-axes
    axis(1, seq(xmin,xmax,0.1), labels = FALSE, lwd=3)

# Add labels to x-axis      
    text(seq(xmin,xmax,length.out=7), rep(0.65, (xmax*2)+1), srt = 0, adj = 0.5, labels = sprintf("%1.1f",seq(xmin,xmax,length.out=7)), xpd = TRUE, cex = lab.cex, font=2)

# Add titles above the top
    text(xmax/2, max(ys.numbers)+0.23, labels="INCREASES WITH TIME",cex=lab.cex, font=2, adj=0.5)
    text(xmin/2, max(ys.numbers)+0.23, labels="DECREASES WITH TIME",cex=lab.cex, font=2, adj=0.5)


# Y-AXIS
    # Define the y-coordinates
    ys<- seq(1.2,y.height,length.out=n.mm)
    ys.numbers<- ys.bars<- c(ys - 0.1, ys, ys + 0.1)
    ys.bars<- sort(ys.bars)

# Add labels to y-axis
    text(rep(xmin,n)-0.02, ys, adj = 1, labels = ns.labels, xpd = TRUE, cex = lab.cex, font=2)
    text(rep(xmin,n)-0.02, ys-0.15, adj = 1, labels = ns.labels2, xpd = TRUE, cex = lab.cex, font=2)


    for (k in 1:12)
    {
      counter<- c(k, k+12, k+24)
      
      for (j in 1:imputations)
      {
        xs<- fit_me2[[k]][[j]]$Sol[,baseline]
        xs[xs<xmin]<- xmin
        xs[xs>xmax]<- xmax
        random.y<- runif(nrow(xs),0,0.02) - 0.01
        
        
        ys<- ys.numbers[counter]
        o2<- order(abs(random.y), decreasing=FALSE)
        o3<- order(abs(xs[,1]-mean[k]), decreasing=TRUE)
        points(xs[o3,1],rep(ys[1],nrow(xs)) + random.y[o2] + rnorm(10,0,0.02) - 0.01, cex=0.5, col=alpha(cols[3], 0.05), pch=19)
        
        o2<- order(abs(random.y), decreasing=FALSE)
        o3<- order(abs(xs[,2]-mean[k+12]), decreasing=TRUE)
        points(xs[o3,2],rep(ys[2],nrow(xs)) + random.y[o2] + rnorm(10,0,0.02) - 0.01, cex=0.5, col=alpha(cols[5], 0.05), pch=19)
        
        o2<- order(abs(random.y), decreasing=FALSE)
        o3<- order(abs(xs[,3]-mean[k+24]), decreasing=TRUE)
        points(xs[,3],rep(ys[3],nrow(xs)) + random.y[o2] + rnorm(10,0,0.02) - 0.01, cex=0.5, col=alpha(cols[9], 0.05), pch=19)
        
        
      }
    }


# Add segments  
    # Bars
    cols.bars<- rep(cols[c(4,6,10)], times=n)
    segments(lower, ys.bars, upper, ys.bars, lwd=lwd, col=cols.bars)

# LEFT end markers
    tips<- setdiff(1:length(lower), oor.min)
    segments(lower[tips], ys.bars[tips]-0.02,lower[tips],ys.bars[tips]+0.02, lwd=lwd, col=cols.bars[tips])

# RIGHT end markers
    tips<- setdiff(1:length(upper), oor.max)
    segments(upper[tips],ys.bars[tips]-0.02,upper[tips],ys.bars[tips]+0.02, lwd=lwd, col=cols.bars[tips])



# Beta & Credible Intervals
    col.inc<- 0.05
    s<- xmax+0.01

    text(rep(s,20), c(ys.bars,max(ys.numbers)+0.23), labels=c(sprintf("%1.2f",mean),"BETA"),cex=lab.cex, font=2, adj=0.5)
    text(rep(s+col.inc,length(ys.bars)+1), c(ys.bars,y.height+.21), labels=c(sprintf("%1.2f", lower), "LOWER"),cex=lab.cex, font=2, adj=0.5)
    text(rep(s+col.inc*2,length(ys.bars)+1), c(ys.bars,y.height+.21), labels=c(sprintf("%1.2f", upper), "UPPER"),cex=lab.cex, font=2, adj=0.5)
    text(s+col.inc*1.5, max(ys.numbers)+0.23, labels="CI",cex=lab.cex, font=2, adj=0.5)

# Effect size
    cols.sd<- rep(grep("_sd",colnames(scaling.data)), each=3)
    rows.means<- rep(match(ns, rownames(scaling.data)), each=3)
    rounding.digits<- rep(c(1,2,1,2,1,2,2,1,0,1,1,1), each=3)


# Effect size & Credible Intervals
    effect.size<- sapply(1:length(rounding.digits), function(x) {sprintf(paste("%1.",rounding.digits[x],"f",sep=""), (scaling.data[1,cols.sd] * mean)[x])})
    text(rep(s+col.inc*3,21), c(ys.bars,y.height+0.21, max(ys.numbers)+0.23), labels=c(effect.size, "SIZE","EFFECT"), cex=lab.cex, font=2, adj=0.5)

# Credible Intervals
    effect.size.lower<- sapply(1:length(rounding.digits), function(x) {sprintf(paste("%1.",rounding.digits[x],"f",sep=""), (scaling.data[1,cols.sd] * lower)[x])})
    text(rep(s+col.inc*4,length(ys.bars)+1), c(ys.bars,y.height+.21), labels=c(effect.size.lower, "LOWER"),cex=lab.cex, font=2, adj=0.5)
    effect.size.upper<- sapply(1:length(rounding.digits), function(x) {sprintf(paste("%1.",rounding.digits[x],"f",sep=""), (scaling.data[1,cols.sd] * upper)[x])})
    text(rep(s+col.inc*5,length(ys.bars)+1), c(ys.bars,y.height+.20), labels=c(effect.size.upper, "UPPER"),cex=lab.cex, font=2, adj=0.5)

    text(s+col.inc*4.5, max(ys.numbers)+0.23, labels="CI",cex=lab.cex, font=2, adj=0.5)
    text(rep(s+col.inc*6,length(ys.bars)+1), c(ys.bars,max(ys.numbers)+0.23), labels=c(pvals,"Pr(>0)"),cex=lab.cex, font=2, adj=0.5)

    segments(s+col.inc, 0.5, s+col.inc*5, 0.5, lwd=lwd*5, col=cols.bars[1], lend=2)
    text(s+col.inc*3, 0.5, labels = "No treatment", cex = lab.cex, font=2, adj=0.5, col="black")
    segments(s+col.inc, 0.65, s+col.inc*5, 0.65, lwd=lwd*5, col=cols.bars[2], lend=2)
    text(s+col.inc*3, 0.65, labels = "Treated with PDE5i", cex = lab.cex, font=2, adj=0.5, col="black")
    segments(s+col.inc, 0.8, s+col.inc*5, 0.8, lwd=lwd*5, col=cols.bars[3], lend=2)
    text(s+col.inc*3, 0.8, labels = "Treated with ERA", cex = lab.cex, font=2, adj=0.5, col="black")


# Add point spheres
    points(mean, sort(ys.bars), lwd=3, cex=1.2, col=cols.bars, pch=19)

# Add vertical lines at 'zero'
      segments(0, .8, 0, y.height+0.2, lty=2, lwd=2)


      
dev.off()


