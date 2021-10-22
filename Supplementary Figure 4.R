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
