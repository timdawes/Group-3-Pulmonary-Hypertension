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
# Copyright Tim Dawes, October 2021
#
# Supplementary Figure 1: Barchart of missingness of data

fields<- c("Outcome","futime","DoseEqPDE5","DoseEqERA",
           "distExTest","mPAP","Kcoperc","CO","PVR","SaO2","SvO2","EmPHasis10",
           "BNP","FC","AgeAtDiag","Gender","BSA","FEV1perc","FVCperc","TLcoperc","Kcoperc","mPAP","CO","PVR",
           "SR","RV.dilat","TAPSE","TRvel","LVEF.Simpson","FC",
           "EmPHasis10","BNP")

length(fields)

units<- c("","(years)","","",
          "(m)","(mmHg)","(%)","(L/min)","(WU)","(%)","(%)","",
          "(ng/L)","","(years)","(M)", "(kg/m2)", "(L)","(L)","%","%","(mmHg)","(L/min)","(WU)",
          "", "","(cm)","(m/s)","(%)","",
          "","(ng/ml)")

fielddisplay<- c("Survival","Follow-up Time", "Phosphodiesterase 5 Inhibitor Treatment", "Endothelin Receptor Antagonist Treatment",
                 "Six Minute Walk Distance", "Mean Pulmonary Artery Pressure", "CO Transfer Coefficient (%)","Cardiac Output", "Pulmonary Vascular Resistance", "Arterial saturations","Venous saturations","emPHasis-10",
                 "B-type Natriuretic Peptide", "WHO Functional Class", "Age","Gender","Body Surface Area","Forced Expiratory Volume (%)","Forced Vital Capacity (%)","Diffusing Capacity for CO (%)","CO Transfer Coefficient (%)","Mean Pulmonary Artery Pressure","Cardiac Output","Pulmonary Vascular Resistance",
                 "Sinus Rhythm","Right Ventricular Dilatation","TAPSE","Tricuspid Regurgitation Velocity","Left Ventricular Ejection Fraction","WHO Functional Class",
                 "emPHasis-10","B-type Natriuretic Peptide")


nas<- rep(NA, length(fields))
names(nas)<- fields


missing.data<- list(echo.complete.non.i, rhc.complete.non.i, pft.complete.non.i, 
                    ex.complete.non.i, Bloods.complete.non.i, FC.complete.non.i,
                    EmPHasis10.complete.non.i, Rx.data)

      # Screen echo data
      for (j in 1:8)
      {
        
        fields.match<- match(fields, colnames(missing.data[[j]]))
        col.nas<- which(is.na(fields.match)==FALSE)
        col.fields<- na.omit(fields.match)
        
        for (i in 1:length(col.nas)){nas[col.nas[i]]<- percent.nas(missing.data[[j]], fields.match, i)}
        cat("\n Missing fields:",sum(is.na(nas)))
      }


      main<- 4
      fu<- 10
      baseline<- 18

# Reformat the data on missing fields for a chart
      df<- data.frame(field = fields, units = units, nas = nas, fielddisplay = fielddisplay, percentcomplete = 100 - nas,
                      group = c(rep("main",main),
                                rep("fu",fu),
                                rep("baseline", baseline)))

# Run this lot if you want repeat fields between survival/baseline/follow-up categories
      dups<- duplicated(df$field)
      df$field[dups]<- paste(df$field[dups],"'",sep="")
      df<- arrange(df, group, percentcomplete, field)
      df$field<- factor(df$field, levels = df$field)

# Run this lot if you don't want repeat fields between survival/baseline/follow-up categories
      df<- df[match(unique(df$field), df$field),]

# Percentage imputation in baseline data
      r<- 1:baseline
      100-mean(df[r,]$percentcomplete)

      p1<- ggplot(df[r,], aes(fill=group, y=percentcomplete, x=field)) + geom_bar(position="dodge", stat="identity", width = 0.9, fill = cols[2]) + coord_flip() +
        theme(axis.text.x = element_text(face = "bold", colour = "black", size = 10)) +
        theme(axis.text.y = element_text(face = "bold", colour = "black", size = 10)) +
        scale_x_discrete(labels = df$fielddisplay[r]) +
        labs(y = "Complete (%)", x = "") + 
        theme(panel.background = element_blank()) + theme(legend.position = "none")


# Percentage imputation in follow-up data
      r<- (baseline+1):(baseline+fu)
      100-mean(df[r,]$percentcomplete)

      p2<- ggplot(df[r,], aes(fill=group, y=percentcomplete, x=field)) + geom_bar(position="dodge", stat="identity", width = 0.9, fill = cols[4]) + coord_flip() +
        theme(axis.text.x =  element_text(face = "bold", colour = "black", size = 10)) +
        theme(axis.text.y = element_text(face = "bold", colour = "black", size = 10)) +
        scale_x_discrete(labels = df$fielddisplay[r]) +
        labs(y = "Complete (%)", x = "") + 
        theme(panel.background = element_blank()) + theme(legend.position = "none")


# Percentage imputation in main data
      r<-(baseline+fu+1):(baseline+fu+main)
      100-mean(df[r,]$percentcomplete)

      p3<- ggplot(df[r,], aes(fill=group, y=percentcomplete, x=field)) + geom_bar(position="dodge", stat="identity", width = 0.9, fill = cols[6]) + coord_flip() +
        theme(axis.text.y = element_text(face = "bold", colour = "black", size = 10)) +
        theme(axis.text.x = element_text(face = "bold", colour = "black", size = 10)) +
        scale_x_discrete(labels = df$fielddisplay[r]) +
        labs(y = "Complete (%)", x = "") + 
        theme(panel.background = element_blank()) + theme(legend.position = "none")

# Combine plots
      g1 <- ggplotGrob(p1)
      g2 <- ggplotGrob(p2)
      g3 <- ggplotGrob(p3)


      fg1 <- gtable_frame(g1, width = unit(1, "null"), height = unit(9, "null"), debug = FALSE)
      fg2 <- gtable_frame(g2, width = unit(1, "null"), height = unit(4.5, "null"), debug = FALSE)
      fg3 <- gtable_frame(g3, width = unit(1, "null"), height = unit(2.2, "null"), debug = FALSE)

      fg123 <- gtable_frame(gtable_rbind(fg3, fg2, fg1), width = unit(1, "null"), height = unit(1, "null"))
      grid.newpage()


# Output plot
      setEPS()
      postscript(file="SFigure1.eps", width=7.08, height=8.85)
      
      grid.draw(fg123)
      dev.off()
