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
# Supplementary Figure 2: Graphs of imputation accuracy by field

n.draws<- 2
oi.plot<- oi<- list()
overimpute.vec<- NULL

for (i in 1:length(df$field)) {overimpute.vec<- c(overimpute.vec, is.numeric(results.olps.i.amelia$imputations[[1]][,match(df$field[i], colnames(results.olps.i.amelia$imputations[[1]]))]))}
mae<- mae.pc<- rep(0, length=length(overimpute.vec))

for (i in 1:nrow(df))
{
  cat("\n",i,"...", df$field[i])
  
  # Generate some overimputed values
  oi[[i]]<- overimpute(results.olps.i.amelia.as.numeric, df$field[i], draws=n.draws)
  oi.df<- data.frame(x=oi[[i]]$orig, y=oi[[i]]$mean.overimputed, fill=oi[[i]]$prcntmiss)
  
  plot.title.size<- 28
  plot.subtitle.size<- 26
  axis.title.size<- 28
  
  
  if (overimpute.vec[i]==FALSE)
  {
    # Plot heatmap for discrete variables (e.g. RV dilatation, SR)    
    
    tab<- table(round(oi.df[,1:2]))
    tab.pc<- as.numeric(format(round(100 * tab / sum(tab),1), nsmall=1))
    oi.df3<- data.frame(Observed = c("+","-","+", "-"), Imputed = c("+", "+", "-", "-"), Frequency = c(tab.pc)[c(4,3,2,1)])
    conf.mat<- confusionMatrix(tab)
    
    oi.plot[[i]]<- ggplot(data = oi.df3, aes(x=Observed, y=Imputed, fill=Frequency)) +
      geom_tile() +
      geom_text(aes(label=sprintf("%0.1f", round(Frequency, digits = 2))), size=24) +
      scale_fill_gradient(limits=c(0,max(tab.pc)), low="white", high=cols[2]) +
      theme(legend.position = "none", panel.background = element_blank()) +
      xlab(df$x.title[i]) +
      ylab(df$y.title[i]) +
      ggtitle(as.character(paste("Accuracy: ", format(round(100*conf.mat$byClass[11],digits=1),nsmall=1), "%",sep=""))) +
      theme(axis.text=element_text(family="Arial", size=24),
            axis.title = element_text(size=axis.title.size),
            plot.title = element_text(size=plot.title.size, face="bold", hjust = 0.5),
            plot.subtitle = element_text(size=plot.subtitle.size, hjust = 0.5),
            axis.line = element_line(colour = "black", size=2))
  } else 
  {
    # Plot scatterplot for continuous variables
    
    mae[i]<- mean(abs(oi[[i]]$mean.overimputed - oi[[i]]$orig))
    mae.pc[i]<- 100 * mae[i] / mean(oi[[i]]$orig)
    cat("mae:",round(mae[i],1), "   mae(%):", round(mae.pc[i],1))
    # df for the background points (s reduces the number for smaller file sizes)
    #s<- sample(1:(nrow(oi$overimps)*ncol(oi$overimps)), 5000)
    oi.df2<- data.frame(x=rep(oi[[i]]$orig, n.draws), y=unlist(c(oi[[i]]$overimps)),
                        fill=rep(oi[[i]]$prcntmiss, n.draws))#[s[1:10000],]
    
    range<- range(c(oi.df$x, oi.df2$x, oi.df$y, oi.df2$y))
    range.min<- floor(range[1])
    range.max<- ceiling(range[2])
    span<- (range.max - range.min)
    step.options<- c(1,2,5,10,50,100,200,500,1000)
    step.size<- step.options[which.min(abs(log(span / step.options) - log(5)))]
    if (step.size>range.min) {range.min<- max(0, range.min - step.size)}  
    
    oi.plot[[i]] = ggplot(oi.df2, aes(x=x, y=y)) +
      geom_point(aes(x=x, y=y, fill=fill), alpha=0.05, size=0.2, shape=19, colour=cols[1],
                 position=position_jitter(width=0.1, height=0.1)) +
      geom_point(data = oi.df, aes(x=x, y=y, fill=fill), alpha=1, size=2, shape=19, colour=cols[2],
                 position=position_jitter(width=0.1, height=0.1)) +
      geom_abline(intercept=0, slope=1, colour=cols[6], size=2) +
      theme(legend.position = "none", panel.background = element_blank()) +
      xlab(df$x.title[i]) +
      ylab(df$y.title[i]) +
      #labs(subtitle = as.character(paste(format(round(mae.pc[i],digits=1),nsmall=1), "%", sep=""))) +
      ggtitle(df$fielddisplay[i]) +
      labs(subtitle = as.character(paste(format(round(mae.pc[i],digits=1),nsmall=1), "%", sep=""))) +
      theme(axis.text=element_text(family="Arial", size=24),
            axis.title = element_text(size=axis.title.size),
            plot.title = element_text(size=plot.title.size, face="bold", hjust = 0.5),
            plot.subtitle = element_text(size=plot.subtitle.size, hjust = 0.5),
            axis.line = element_line(colour = "black", size=2)) +
      scale_x_continuous(limits=c(range.min, range.max), breaks=seq(range.min, range.max, step.size)) +
      scale_y_continuous(limits=c(range.min, range.max), breaks=seq(range.min, range.max, step.size))
  }   
  
  filename.tiff<- paste(df$field[i],".tiff",sep="")
  tiff(file=filename.tiff, height=500,width=500)
  grid.arrange(oi.plot[[i]], ncol=1)
  dev.off()
  
} 

# Save the plot  







filenames<- list.files(pattern=".tiff")
margin = theme(plot.margin = unit(c(0,0,0,0), "in"))

cat.plots<- match(c("SR.tiff","RV.dilat.tiff"), filenames)
cont.plots<- setdiff(1:length(filenames), cat.plots)
all.plots<- c(cont.plots,cat.plots)

plots<- list()
for(j in 1:length(all.plots)) {plots[[j]]<-rasterGrob(readTIFF(filenames[all.plots[j]]),
                                                      width = unit(1.8, "in"),
                                                      height = unit(1.8, "in"))}


tiff(file="All.tiff", height=9.05, width=7.09, units="in", res=300)
lay<- rbind(1:3, 4:6, 7:9, 10:12, 13:15)

grid.arrange(grobs = plots, ncol = 3, margin = margin, layout_matrix = lay)
dev.off()

plots<- list()
for(j in 1:length(cat.plots)) {plots[[j]]<-rasterGrob(readTIFF(filenames[cat.plots[j]]),
                                                      width = unit(4, "in"),
                                                      height = unit(4, "in"))}


tiff(file="All_IvsO_categorical.tiff", height=4, width=8, units="in", res=300)
lay<- matrix(c(1,2),nrow=1)
grid.arrange(grobs = plots, ncol = 2, margin = margin, layout_matrix = lay)
dev.off()



plots<- list()
for(j in 1:length(cont.plots)) {plots[[j]]<-rasterGrob(readTIFF(filenames[cont.plots[j]]),
                                                       width = unit(4, "in"),
                                                       height = unit(4, "in"))}
tiff(file="All_IvsO_continuous.pdf", height=10, width=8, units="in", res=300)
lay<- rbind(1:3, 4:6, 7:9, 10:12, c(NA,13,NA))
grid.arrange(grobs = plots, ncol = 3, margin = margin, layout_matrix = lay)
dev.off()


plots<- list()
cont.and.cat<- c(cont.plots, cat.plots)
for(j in 1:length(filenames)) {plots[[j]]<-rasterGrob(readPNG(filenames[cont.and.cat[j]]),
                                                      width = unit(2, "in"),
                                                      height = unit(2, "in"))}
pdf(file="All_Ivs.pdf", height=11, width=8)
lay<- rbind(1:3, 4:6, 7:9, 10:12, 13:15)
grid.arrange(grobs = plots, ncol = 3, margin = margin, layout_matrix = lay)
dev.off()


