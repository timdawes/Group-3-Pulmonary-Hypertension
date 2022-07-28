# Phosphodiesterase 5 inhibitor treatment and survival in interstitial lung disease pulmonary hypertension: a Bayesian retrospective observational cohort study

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
# Supplementary Figure 4: Is the diagnosis associated with number of RHCs?
                
# Find the URNs for all the right heart catheters
    y<- table(d2$URN[match(rhc.complete.non.i$ID, d2$URN)])
    x<- d2$Diagnosis_sub[match(names(y), d2$URN)]

# Group them into COPD, ILD or other
    x[which(x %in% c("3.3","3.4","3.5","3.7"))]<- "Other"
    x[which(x %in% c("3.2"))]<- "ILD"
    x[which(x %in% c("3.1"))]<- "COPD"

# GLM 
    x.mat<- model.matrix(~ x)[,-1]
    poi.fit<- glm(y ~ x.mat, family = quasipoisson(link="log"))
    summary(poi.fit)

# Reformat the data for a barplot  
    barplot.data<- sort(unique(x))
    t<- matrix(0, nrow=3, ncol=3, dimnames=list(c("One","Two","Three"),barplot.data))
    for (i in 1:3) {t[,i]<- table(factor(y[which(x==barplot.data[i])], levels=(1:3)))}

# Draw the barplot
    cols<- brewer.pal(11,"Blues")[c(5,7,9)]
    t.long<- cbind(gather(data.frame(t), Diagnosis, Frequency), RHC=factor(rep(c("One","Two","Three"),3), levels=c("One","Two","Three")))
    p1<- ggplot(data=data.frame(t.long), aes(x=Diagnosis, y=Frequency, fill=RHC, label=Frequency)) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          geom_col(position = position_dodge(w = -0.5)) +
            scale_fill_manual(values = rep(cols,2)) +
            theme(axis.line = element_line(colour = 'black', size = 1.5),
                  axis.text = element_text(size=16),
                  axis.title=element_text(size=16,face="bold")) +
            geom_text(aes(y=Frequency+2), size=6, position=position_dodge(w = -0.5)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p1   




  # Is PDE5i treatment associated with number of RHCs? (GLM with binomial link function)
                    x<- table(d2$URN[match(rhc.complete.non.i$ID, d2$URN)])
                    y<- d2$PDE[match(names(x), d2$URN)]
                    
                    bi.fit<- glm(y ~ x, family = binomial(link="logit"))
                    summary(bi.fit)
                    
                    d<- sort(unique(y))
                    t2<- matrix(0, nrow=3, ncol=2, dimnames=list(c("One","Two","Three"),c("No","Yes")))
                    for (i in 1:2)
                    {
                      z<- d[i]
                      s<- table(factor(x[which(y==z)], levels=(1:3)))
                      t2[,i]<- s
                    }
                    
                    cols<- brewer.pal(11,"Reds")[c(5,7,9)]
                    t.long2<- cbind(gather(data.frame(t2), Treatment, Frequency), RHC=factor(rep(c("One","Two","Three"),2), levels=c("One","Two","Three")),
                                   col=rep(cols,2))
                    p2<- ggplot(data=data.frame(t.long2), aes(x=Treatment, y=Frequency, fill=RHC, label=Frequency)) +
                        geom_col(position = position_dodge(width = +0.3)) +
                        scale_fill_manual(values = rep(cols,2)) +
                        theme(axis.line = element_line(colour = 'black', size = 1.5),
                              axis.text = element_text(size=16),
                              axis.title=element_text(size=16,face="bold")) +
                        geom_text(aes(y=Frequency+2), size=6, position=position_dodge(width= +0.3)) +
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                      xlab("Treated with phophodiesterase 5 inhibitors")
            
                    p3<- ggarrange(p1, p2, labels = c("A", "B"), ncol = 1, nrow = 2)
                    
                    pdf(file="TreatmentVsRHCs.pdf", width=6, height=6, onefile = FALSE)
                    p3
                    dev.off()
                    
          
