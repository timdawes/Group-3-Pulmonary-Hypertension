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

# Input data

      # Change to data subfolder and read in Excel files
            setwd(working.dir)
      
      # Load a vector of the expected data types (col_types) and then load the data
            
            col_types = as.character(read.csv("Data/col_types.csv")[,1])
            d<- read_excel("Data/Data.xlsx", col_types=col_types)
            DrugConversion<- read_excel("Data/DrugConversion.xlsx")
            d2<- d[which(d$Include.3.PH.PFT.Echo=="Yes"),]
            
            cat("done")
            

# Treatment data
            
            cat("\n Loading treatment data...")
            
            r<- range(grep("^DrugDose",colnames(d2)))
                dose.start<- r[1]
                dose.stop<- r[2]
                
            r<- range(grep("^DrugName",colnames(d2)))
                drug.name.start<- r[1]
                drug.name.stop<- r[2]
            
                drug.start<- which(colnames(d2)=="date.DrugStart1")
                drug.stop<- which(colnames(d2)=="date.DrugStop1")
                
            
            
            
      # Loop through patients sequentially building 'data2' type dataframe with treatment data
              
              no.rows<- 432
          
          # Treatment dates
              Rx.dates<- data.frame(DateOfBirth = as.Date(rep("1900-01-01", length=no.rows, format="%y/%m/%d")),
                                    DateOfDiagnosis = as.Date(rep("1900-01-01", length=no.rows, format="%y/%m/%d")),
                                    DateOfScan = as.Date(rep("1900-01-01", length=no.rows, format="%y/%m/%d")),
                                    DateOfCensoring = as.Date(rep("1900-01-01", length=no.rows, format="%y/%m/%d")))
           
          # Treatment demographics
              Rx.demographics<- vector(length=no.rows, mode="character")
              
          # Treatment data
              Rx.data<- data.frame(futime = rep(0, length=no.rows),
                                   DoseEqPDE5 = rep(0, length=no.rows),
                                   DoseEqERA = rep(0, length=no.rows),
                                   DoseEqPG = rep(0, length=no.rows),
                                   DoseEqCCB = rep(0, length=no.rows),
                                   DoseEqGCS = rep(0, length=no.rows),
                                   DoseEqSLP = rep(0, length=no.rows),
                                   DoseEqLTOT = rep(0, length=no.rows),
                                   DoseEqAny = rep(0, length=no.rows),
                                   Outcome = rep(0, no.rows),
                                   Treated = rep(FALSE, no.rows),
                                   TreatedPDE5 = rep(FALSE, no.rows),
                                   TreatedERA = rep(FALSE, no.rows))
              
          # Row counter for adding rows
              counter<- 1
          
          # Treatment families
              drug.groups<- c(unique(DrugConversion$GROUP), "Any") # Plus 'any treatment'
              
          # Loop through the cohort adding rows for each treatment (often >1 per patient) or for no treatment (one row)
              for (i in 1:nrow(d2))
                  {
                      drug.cols<- which(is.na(d2[i,dose.start:(dose.start+15)])==FALSE)
                      
                      # Add rows to the Rx dataset if this patient is on treatment
                      if (length(drug.cols)>0)
                            {
                              # Find the right row of the Drug Conversion.xlsx file for each drug this patient has taken 
                                    rows.in.DrugConversion<- match(d2[i,drug.name.start + drug.cols - 1], DrugConversion$INFOFLEX_CODE)
                              # Find the max dose for each drug the patient has taken from Drug Conversion.xlsx 
                                    max.dose.for.each.drug<- DrugConversion$`MAX_DOSE/DAY_(MG)`[rows.in.DrugConversion]
                              # Find the drug class for each drug the patient has taken from Drug Conversion.xlsx
                                    drug.class.for.each.drug<- DrugConversion$`GROUP`[rows.in.DrugConversion]
                              # Express the doses taken as a fraction of the max dose
                                    dose.equivalents<- unlist(d2[i,dose.start + drug.cols - 1] / max.dose.for.each.drug)
                              # Make a continuous timeline for this patient, starting at diagnosis, through all treatment titrations, finishing at censoring
                                    new.drug.dates<- make.contiguous.dates(d2[i,], drug.start, drug.stop, drug.cols, dose.equivalents, drug.class.for.each.drug)
                              
                                    unique.dates<- sort(unique(c(new.drug.dates$start.dates, new.drug.dates$stop.dates)))
                                    Rx.matrix<- matrix(0, nrow=length(unique.dates), ncol=length(drug.groups), dimnames=list(as.character(unique.dates), c(drug.groups)))
                                    
                                    for (j in 1:length(drug.groups))
                                          {
                                            # Find the drugs this patient has received which fall in each class
                                                  doses.by.class<- which(new.drug.dates$drug.class.for.each.drug %in% drug.groups[j])
                                                  
                                                  if (length(doses.by.class)>0){
                                                        r.start<- match(new.drug.dates$start.dates[doses.by.class], ymd(rownames(Rx.matrix)))
                                                        r.end<- match(new.drug.dates$stop.dates[doses.by.class], ymd(rownames(Rx.matrix)))
                                                        
                                                  # Add the dose equivalences to the appropriate columns of the treatment matrix (Rx.matrix)
                                                        for (k in 1:length(r.start)) {Rx.matrix[r.start[k]:(r.end[k]-1), j]<- new.drug.dates$dose.equivalents[doses.by.class[k]]}}
                                          }
                                    
                              x<- counter:(counter+length(unique.dates)-1)
                              n<- length(x)
                              
                          # Add n rows to 'demographics'
                              Rx.demographics[x]<- d2$URN[i]
                              
                          # Add n rows to 'dates'
                              Rx.dates$DateOfBirth[x]<- new.drug.dates$DateOfBirth
                              Rx.dates$DateOfDiagnosis[x]<- new.drug.dates$diagnosis.date
                              Rx.dates$DateOfScan[x]<- rownames(Rx.matrix)
                              Rx.dates$DateOfCensoring[x]<- new.drug.dates$censoring.date
                              
                              
                          # Add n rows to 'data'
                              Rx.data$DoseEqPDE5[x]<- Rx.matrix[,1]
                              Rx.data$DoseEqERA[x]<- Rx.matrix[,2]
                              Rx.data$DoseEqPG[x]<- Rx.matrix[,3]
                              Rx.data$DoseEqCCB[x]<- Rx.matrix[,4]
                              Rx.data$DoseEqGCS[x]<- Rx.matrix[,5]
                              Rx.data$DoseEqSLP[x]<- Rx.matrix[,6]
                              Rx.data$DoseEqLTOT[x]<- Rx.matrix[,7]
                              Rx.data$DoseEqAny[x]<- apply(Rx.matrix,1,sum)
                              
                              
                              Rx.data$futime[x]<- rep(d2$Date.censoring[i] - d2$dateDIAGNOSIS[i], n) / 365.25
                              Rx.data$Outcome[x]<- rep(d2$Censored[i], n)
                              Rx.data$Treated[x]<- TRUE
                              Rx.data$TreatedPDE5[x]<- any(Rx.matrix[,1]!=0)
                              Rx.data$TreatedERA[x]<- any(Rx.matrix[,2]!=0)
                              
                           
                          # Add 'n' to the row counter
                              counter<- counter + n
                          
                      } else {
                        
                        # Add a single row to the Rx dataset if the patient had no treatment
                              x<- counter
                             
                              # Add n rows to 'demographics'
                                  Rx.demographics[x]<- d2$URN[i]
                              
                              # Add n rows to 'dates'
                                  Rx.dates$DateOfBirth[x]<- ymd(d2$DoB[i])
                                  Rx.dates$DateOfDiagnosis[x]<- ymd(d2$dateDIAGNOSIS[i])
                                  Rx.dates$DateOfScan[x]<- ymd(d2$dateDIAGNOSIS[i])
                                  Rx.dates$DateOfCensoring[x]<- ymd(d2[i,]$Date.censoring)
                              
                              # Add n rows to 'data'
                                  Rx.data$DoseEqPDE5[x]<- Rx.data$DoseEqERA[x]<- Rx.data$DoseEqPG[x]<- Rx.data$DoseEqCCB[x]<- Rx.data$DoseEqGCS[x]<- Rx.data$DoseEqSLP[x]<- Rx.data$DoseEqLTOT[x]<- Rx.data$DoseEqAny[x]<- 0
                                  Rx.data$futime[x]<- (d2$Date.censoring[i] - d2$dateDIAGNOSIS[i]) / 365.25
                                  Rx.data$Outcome[x]<- d2$Censored[i]
                                  Rx.data$Treated[x]<- FALSE
                                  
                              # Add 'n' to the row counter     
                                  counter<- counter + 1
                      }
                      
              }
          
          
          Rx.data<- cbind(ID = Rx.demographics, Rx.dates, Rx.data)
          
          Rx.data$AgeAtTest<- (ymd(Rx.data$DateOfScan) - ymd(Rx.data$DateOfBirth)) / 365.25
          Rx.data$AgeAtDiag<- (ymd(Rx.data$DateOfDiagnosis) - ymd(Rx.data$DateOfBirth)) / 365.25
          Rx.data$DoseFctrPDE5<- match(Rx.data$DoseEqPDE5, sort(unique(Rx.data$DoseEqPDE5)))
          Rx.data$DoseFctrERA<- match(Rx.data$DoseEqERA, sort(unique(Rx.data$DoseEqERA)))
          Rx.data$DoseFctrAny<- match(Rx.data$DoseEqAny, sort(unique(Rx.data$DoseEqAny)))
          
          
          # Time-dependent covariates    
          
                
          # Make a dataframe which has one line per patient/imputation
                d2.demo<- data.frame(d2[,c("DoB","Date.censoring","URN","Diagnosis","Gender","Dead","FUTimeYears")]) #temp = data1.joint
                colnames(d2.demo)<- c("DateOfBirth","DateOfCensoring","ID","Diagnosis","Gender","Outcome","futime")
                d2.demo$Outcome<- as.numeric(d2.demo$Outcome)
                futime<- Outcome<- rep(0, nrow(d2.demo))
                data1<- tmerge(d2.demo, d2.demo, id=ID, Outcome = event(futime, Outcome))
                
                
          results.Rx.olps<- tmerge(data1, make.tdc.matrix(data1, Rx.data, "DateOfScan"), id=ID, AgeAtDiag = tdc(time, AgeAtDiag), DoseEqPDE5 = tdc(time, DoseEqPDE5), DoseEqERA = tdc(time, DoseEqERA), DoseEqPG = tdc(time, DoseEqPG),
                                   DoseEqCCB = tdc(time, DoseEqCCB), DoseEqGCS = tdc(time, DoseEqGCS), DoseEqSLP = tdc(time, DoseEqSLP), DoseEqLTOT = tdc(time, DoseEqLTOT), Treated = tdc(time, Treated),
                                   TreatedPDE5 = tdc(time, TreatedPDE5), DoseFctrPDE5 = tdc(time, DoseFctrPDE5),
                                   TreatedERA = tdc(time, TreatedERA), DoseFctrERA = tdc(time, DoseFctrERA),
                                   DoseFctrAny = tdc(time, DoseFctrAny))
          
          results.Rx.olps<- change.data.type(results.Rx.olps)
          
          
          
          
          
          cat("...done")
          
          

          
          
          
# Echo data
          
          cat("\n Loading echo data...")
          
          
          results<- extract.class.columns.and.dates(d2,".Echo")
                echo.present<- results[[1]]
                echo.cols<- results[[2]]
                
          results<- extract.class.data.and.demographics(echo.present, echo.cols, d2, 29)
                echo.data<- results[[1]]
                echo.dates<- results[[2]]
                echo.demographics<- results[[3]]
                
          # Remove duplicate investigations
                t<- remove.duplicates(echo.data, echo.dates, echo.demographics, echo.present)
                  echo.data<- t[[1]]
                  echo.dates<- t[[2]]
                  echo.demographics<- t[[3]]
                  echo.present<- t[[4]]
                  
                   
          
          
          
          
          
          # Imputation
          
                echo.complete.non.i<- cbind(echo.dates, echo.demographics, echo.data)
                
          # Time-dependent covariates    
          
                # Make a dataframe which has one line per patient/imputation
                      d2.demo<- data.frame(d2[,c("DoB","Date.censoring","URN","Diagnosis","Diagnosis_sub","Gender","Dead","FUTimeYears")]) #temp = data1.joint
                      colnames(d2.demo)<- c("DateOfBirth","DateOfCensoring","ID","Diagnosis","Diagnosis_sub","Gender","Outcome","futime")
                      d2.demo$Outcome<- as.numeric(d2.demo$Outcome)
                      futime<- Outcome<- rep(0, nrow(d2.demo))
      
                # Make a second dataframe which has one line per study (data2)
                      data1<- tmerge(d2.demo, d2.demo, id=ID, Outcome = event(futime, Outcome))
                      results.echo.olps<- tmerge(data1, make.tdc.matrix(d2.demo, echo.complete.i, "DateOfScan"), id=ID, RVdys = tdc(time, RVdys), HR = tdc(time, HR), SBP = tdc(time, SBP), SR = tdc(time, SR),
                                                 RV.dilat = tdc(time, RV.dilat), Long.dysfunc = tdc(time, Long.dysfunc), Rad.dysfunc = tdc(time, Rad.dysfunc), RVdys = tdc(time, RVdys), FAC = tdc(time, FAC),
                                                 TAPSE = tdc(time, TAPSE), S = tdc(time, S), TRvel = tdc(time, TRvel), RVSP = tdc(time, RVSP), RA.pressure = tdc(time, RA.pressure), RA.dilat = tdc(time, RA.dilat),
                                                 RA.areaI = tdc(time, RA.areaI), RA.volumeI = tdc(time, RA.volumeI), PAAT = tdc(time, PAAT), Effusion = tdc(time, Effusion), LVEF.Simpson = tdc(time, LVEF.Simpson),
                                                 Edecel = tdc(time, Edecel), E.A = tdc(time, E.A), EE. = tdc(time, EE.), LA.dilat = tdc(time, LA.dilat), LA.areaI = tdc(time, LA.areaI), LA.volumeI = tdc(time, LA.volumeI))
                      
          
                # Convert numbers-to-numbers and factors-to-factors    
                      results.echo.olps<- change.data.type(results.echo.olps)

                      
          cat("...done")

          
# 2. Pulmonary Function Tests
          
          cat("\n Loading pulmonary function tests data...")
          
          results<- extract.class.columns.and.dates(d2,".PFT")
                pft.present<- results[[1]]
                pft.cols<- results[[2]]
          
          results<- extract.class.data.and.demographics(pft.present, pft.cols, d2, 8)
                pft.data<- results[[1]]
                pft.dates<- results[[2]]
                pft.demographics<- results[[3]]
          
          # Remove duplicate investigations
                t<- remove.duplicates(pft.data, pft.dates, pft.demographics, pft.present)
                      pft.data<- t[[1]]
                      pft.dates<- t[[2]]
                      pft.demographics<- t[[3]]
                      pft.present<- t[[4]]
                      
          
          
          # Complete absolute vs percentage values
          absolute.pft.cols<- match(c("FEV1","FVC"),colnames(pft.data))
          
          for (i in absolute.pft.cols) {
             na.absolute<- which(is.na(pft.data[,i])==TRUE)
             na.perc<- which(is.na(pft.data[,i+1])==TRUE)
             percs.missing<- setdiff(na.perc, na.absolute)
             cat(paste("\n Completing missing percentages...",colnames(pft.data)[i],sep=""))
             
             for (j in percs.missing)
             {
                pft.AHW<- pft.demographics[which(pft.demographics$ID %in% pft.demographics$ID[j]), c("Gender","Height","AgeAtTest")]
                if (nrow(pft.AHW)>0)
                {
                   nearest<- which.min(abs(pft.AHW$AgeAtTest-pft.demographics$AgeAtTest[j])) 
                   if (i==absolute.pft.cols[1]) {Predicted.FEV1<- predictFEV1(pft.AHW[nearest,]); pft.data[j,i+1]<- 100 * pft.data[j,i] / Predicted.FEV1}
                   if (i==absolute.pft.cols[2]) {Predicted.FVC<- predictFVC(pft.AHW[nearest,]); pft.data[j,i+1]<- 100 * pft.data[j,i] / Predicted.FVC}
                }
             }
             
             cat(paste("\n Completing missing absolute values...",colnames(pft.data)[i],sep=""))
             abs.missing<- setdiff(na.absolute, na.perc)  
             
             for (j in abs.missing)
             {
                pft.AHW<- pft.demographics[which(pft.demographics$ID %in% pft.demographics$ID[j]), c("Gender","Height","AgeAtTest")]
                if (nrow(pft.AHW)>0)
                {
                   nearest<- which.min(abs(pft.AHW$AgeAtTest-pft.demographics$AgeAtTest[j])) 
                   if (i==absolute.pft.cols[1]) {Predicted.FEV1<- predictFEV1(pft.AHW[nearest,]); pft.data[j,i]<- Predicted.FEV1 * (pft.data[j,i+1] / 100)}
                   if (i==absolute.pft.cols[2]) {Predicted.FVC<- predictFVC(pft.AHW[nearest,]); pft.data[j,i+1]<- Predicted.FVC * (pft.data[j,i+1] / 100)}
                }
             }
          }
          
          
          
          # Imputation
          
                pft.complete.non.i<- cbind(pft.dates, pft.demographics, pft.data)
                
          
          # Time-dependent covariates    
          
                # Make a second dataframe (data2) which has one line per echo
                      data1<- tmerge(d2.demo, d2.demo, id=ID, Outcome = event(futime, Outcome))
                      data2<- make.tdc.matrix(d2.demo, pft.complete.i,"DateOfScan")
                      
                      data1<- change.data.type(data1)
                      data2<- change.data.type(data2)
                      
                      results.pft.olps<- tmerge(data1, data2, id=ID, FEV1 = tdc(time, FEV1), FEV1perc = tdc(time, FEV1perc), FVC = tdc(time, FVC), FVCperc = tdc(time, FVCperc),
                                                TLcoc = tdc(time, TLcoc),  TLcoperc = tdc(time, TLcoperc), Kco = tdc(time, Kco), Kcoperc = tdc(time, Kcoperc))
                      
          cat("...done")
          
          
          
          
# Right Heart Catheter data
          
          cat("\n Loading RHC data...")
          
          # Define some variables
                rhc.present<- NULL
                rhc.date<- list()
                factors<- 1
                
          # a. Identify which columns the rhc data and the rhc dates are stored in
                rhc.cols<- unlist(c(sapply(paste(".RHC",1:4,sep=""), function(x) {grep(x, colnames(d2))})))
                rhc.date.cols<- sort(unique(rhc.cols[grep("Date.RHC", colnames(d2[,rhc.cols]))]))
                
          # b. Build the matrix storing the rhc date row and columns in    
                for (i in 1:4) {
                   r<- which(is.na(d2[,rhc.date.cols[i]])==FALSE)
                   #cat("\n rhc number:",i," frequency:", length(r))
                   c<- rep(rhc.date.cols[i], length(r))
                   rhc.present<- rbind(rhc.present, cbind(r,c))
                }
                
                cat("\n Number of rhcs available:", nrow(rhc.present))
                
          
          # c. Build the matrix storing the PFT data in
                results<- extract.class.data.and.demographics(rhc.present, rhc.cols, d2, 4)
                      rhc.data<- results[[1]]
                      rhc.dates<- results[[2]]
                      rhc.demographics<- results[[3]]
                      
            rhc.complete.non.i<- cbind(rhc.dates, rhc.demographics, rhc.data)
          
          # Time-dependent covariates    
                results.rhc.olps<- tmerge(data1, make.tdc.matrix(d2.demo, rhc.complete.i, "DateOfScan"), id=ID, mPAP = tdc(time, mPAP), PVR = tdc(time, PVR), CO = tdc(time, CO), PCWP = tdc(time, PCWP))
          
          # Fit COX model
                results.rhc.olps$PVR.bin<- results.rhc.olps$PVR>mean(results.rhc.olps$PVR)
         
          cat("...done")
          
          
          
          
# Exercise data
          
          cat("\n\n Loading 6MWT data...")
          
          
          results<- extract.class.columns.and.dates(d2,".ExTest")
                ex.present<- results[[1]]
                ex.cols<- results[[2]]
                
          results<- extract.class.data.and.demographics(ex.present, ex.cols, d2, 1)
                ex.data<- data.frame(distExTest = results[[1]])
                ex.dates<- results[[2]]
                ex.demographics<- results[[3]]
                
          
          # Remove duplicate investigations
                t<- remove.duplicates(ex.data, ex.dates, ex.demographics, ex.present)
                ex.data<- t[[1]]
                ex.dates<- t[[2]]
                ex.demographics<- t[[3]]
                ex.present<- t[[4]]
                
          
              ex.complete.non.i<- cbind(ex.dates, ex.demographics, ex.data)
          
                results.ex.olps$distExTest[which(is.na(results.ex.olps$distExTest)==TRUE)]<- mean(ex.data$distExTest)
                
          cat("...done")
          
          
          
          
          
          
          
# Quality of life scores
          
          cat("\n\n Loading EmPHasis10 data...")
          
          results<- extract.class.columns.and.dates(d2,".QoLScore")
                EmPHasis10.present<- results[[1]]
                EmPHasis10.cols<- results[[2]]
                
          results<- extract.class.data.and.demographics(EmPHasis10.present, EmPHasis10.cols, d2, 1)
                EmPHasis10.data<- results[[1]]
                EmPHasis10.dates<- results[[2]]
                EmPHasis10.demographics<- results[[3]]
                
          
          # Remove duplicate investigations
                t<- remove.duplicates(EmPHasis10.data, EmPHasis10.dates, EmPHasis10.demographics, EmPHasis10.present)
                EmPHasis10.data<- t[[1]]
                EmPHasis10.dates<- t[[2]]
                EmPHasis10.demographics<- t[[3]]
                EmPHasis10.present<- t[[4]]
                
             EmPHasis10.complete.non.i<- cbind(EmPHasis10.dates, EmPHasis10.demographics, EmPHasis10.data)
                
          # Make a second dataframe which has one line per study (data2)
                results.EmPHasis10.olps<- tmerge(data1, make.tdc.matrix(d2.demo, EmPHasis10.complete.i, "DateOfScan"), id=ID, EmPHasis10 = tdc(time, EmPHasis10))
                results.EmPHasis10.olps$EmPHasis[which(is.na(results.EmPHasis10.olps$EmPHasis10)==TRUE)]<- mean(na.omit(EmPHasis10.data$EmPHasis10))
                
          cat("...done")
          
          
          
# Functional class scores
          
          cat("\n\n Loading Functional Class data...")
          
          results<- extract.class.columns.and.dates(d2,".FC")
                FC.present<- results[[1]]
                FC.cols<- results[[2]]
                
          results<- extract.class.data.and.demographics(FC.present, FC.cols, d2, 1)
                FC.data<- results[[1]]
                FC.dates<- results[[2]]
                FC.demographics<- results[[3]]
                
          # Remove duplicate investigations
                t<- remove.duplicates(FC.data, FC.dates, FC.demographics, FC.present)
                FC.data<- t[[1]]
                FC.dates<- t[[2]]
                FC.demographics<- t[[3]]
                FC.present<- t[[4]]
                
              FC.complete.non.i<- cbind(FC.dates, FC.demographics, FC.data)
                
          
          # Make a second dataframe which has one line per study (data2)
                results.FC.olps<- tmerge(data1, make.tdc.matrix(d2.demo, FC.complete.i, "DateOfScan"), id=ID, FC = tdc(time, FC))
                results.FC.olps$FC[which(is.na(results.FC.olps$FC)==TRUE)]<- mean(na.omit(FC.data$FC))
                
          cat("...done")
          
          
          
          
          
          
# Bloods including BNP
          
          cat("\n Loading Bloods data...")
          
          
          Bloods.all<- read_excel("Data/Bloods.xlsx", sheet=2)
          
          in.cohort<- which((Bloods.all$PATID %in% d2$URN)==TRUE)
          no.in.d2<- na.omit(match(Bloods.all$PATID, d2$URN))
          
          Bloods.in.cohort<- Bloods.all[in.cohort, c('PATID','SampleCollectDateTime','DOB', 'BNP_ng/L', 'Sodium', 'Potassium', 'Urea','Creatinine','GFR','Bilirubin','ALP','ALT',
                                                     'Protein_Tot','Albumin','Corr_calcium', 'Iron',
                                                     'Transferrin_saturation','Transferrin','TIBC','Free_T4')]
          
          
          
          # Format for LOCFandB function (coming up)
                Bloods.in.cohort.for.imputation<- Bloods.in.cohort[order(Bloods.in.cohort$PATID, Bloods.in.cohort$SampleCollectDateTime),-c(2,3)]
                Bloods.names<- c("ID","BNP","Na","K","Ur","Creat","eGFR","Bili","ALP","ALT","Prot","Alb","CCa","Fe","TFSat","TF","TIBC","T4")
                colnames(Bloods.in.cohort.for.imputation)<- Bloods.names
 
         # Make BNP dates
                DateOfScan = as.Date(Bloods.in.cohort$SampleCollectDateTime, format="%Y-%m-%d")
                DateOfBirth = as.Date(Bloods.in.cohort$DOB, format="%Y-%m-%d")
                DateOfCensoring = as.Date(d2$Date.censoring[no.in.d2], format="%Y-%m-%d")
                Bloods.dates<- data.frame(DateOfScan=DateOfScan, DateOfBirth=DateOfBirth, DateOfCensoring=DateOfCensoring)
                
          # Make BNP demographics
                dateE.columns<- grep("date.Echo",colnames(d2))[-1]
                Height<- rep(0, length(no.in.d2))
                
          # Add progress bar
                pb <- txtProgressBar(min = 0, max = length(no.in.d2), style = 2)
          
          for (i in 1:length(no.in.d2))
          {
             nearest.date<- order(abs(ymd(as.Date(unlist(d2[no.in.d2[i],dateE.columns]), origin="1900-01-01")) - ymd(as.Date(Bloods.dates$DateOfScan[i], origin="1900-01-01"))))
             Height[i]<- na.omit(unlist(d2[no.in.d2[i], dateE.columns[nearest.date]+2]))[1]
             setTxtProgressBar(pb, i)
          } 
          
          
          Bloods.demographics<- do.call(rbind, lapply(1:imputations, function(x) {
             data.frame(ID = d2$URN[no.in.d2], Diagnosis = d2$Diagnosis[no.in.d2], Gender = d2$Gender[no.in.d2],
                        AgeAtTest = (ymd(as.Date(DateOfScan)) - ymd(as.Date(d2$DoB[no.in.d2], format="%y-%m-%d"))) / 365.25,
                        Height = Height, Outcome = d2$Censored[no.in.d2])}))
          
          Bloods.complete.non.i<- cbind(Bloods.dates, Bloods.demographics, Bloods.in.cohort.for.imputation) 
          Bloods.complete.non.i$eGFR<- with(Bloods.complete.non.i, eGFR(creat = Creat, age = AgeAtTest, gender = Gender, black = FALSE, method = "CKD-EPI"))
          
          # Time dependent covariates
                # Make a second dataframe which has one line per study (data2)
                results.Bloods.olps<- tmerge(data1, make.tdc.matrix(d2.demo, Bloods.complete.i, "DateOfScan"), id=ID, BNP = tdc(time, BNP), Na = tdc(time, Na), K = tdc(time, K),
                                          Ur = tdc(time, Ur), Creat = tdc(time, Creat), eGFR = tdc(time, eGFR), Bili = tdc(time, Bili), ALP = tdc(time, ALP), ALT = tdc(time, ALT),
                                          Prot = tdc(time, Prot), Alb = tdc(time, Alb), CCa = tdc(time, CCa), Fe = tdc(time, Fe), TFSat = tdc(time, TFSat),
                                          TF = tdc(time, TF), TIBC = tdc(time, TIBC), T4 = tdc(time, T4))
                
                results.Bloods.olps$BNP[which(is.na(results.Bloods.olps$BNP)==TRUE)]<- mean(Bloods.in.cohort$`BNP_ng/L`)
                results.Bloods.olps$Creat[which(is.na(results.Bloods.olps$Creat)==TRUE)]<- mean(Bloods.in.cohort$`Creatinine`)
                results.Bloods.olps$BNP.bin<- results.Bloods.olps$BNP>mean(results.Bloods.olps$BNP)
                
          
          cat("\n ...done")
          
                    
          

 
# Add all the NON-IMPUTED covariates together
          
          cat("\n Combining the non-imputed versions of the fields...")
          
          
          results.olps.non.i<- tmerge(data1, make.tdc.matrix(d2.demo, echo.complete.non.i, "DateOfScan"), DateOfScan = tdc(time, DateOfScan), id=ID, Diagnosis = tdc(time, Diagnosis), AgeAtTest = tdc(time, AgeAtTest),
                                      Height = tdc(time, Height), Weight = tdc(time, Weight), BMI = tdc(time, BMI), BSA = tdc(time, BSA), HR = tdc(time, HR), SBP = tdc(time, SBP), SR = tdc(time, SR),
                                      RV.dilat = tdc(time, RV.dilat), Long.dysfunc = tdc(time, Long.dysfunc), Rad.dysfunc = tdc(time, Rad.dysfunc), RVdys = tdc(time, RVdys), FAC = tdc(time, FAC),
                                      TAPSE = tdc(time, TAPSE), S = tdc(time, S), TRvel = tdc(time, TRvel), RVSP = tdc(time, RVSP), RA.pressure = tdc(time, RA.pressure), RA.dilat = tdc(time, RA.dilat),
                                      RA.areaI = tdc(time, RA.areaI), RA.volumeI = tdc(time, RA.volumeI), PAAT = tdc(time, PAAT), Effusion = tdc(time, Effusion), LVEF.Simpson = tdc(time, LVEF.Simpson),
                                      Edecel = tdc(time, Edecel), E.A = tdc(time, E.A), EE. = tdc(time, EE.), LA.dilat = tdc(time, LA.dilat), LA.areaI = tdc(time, LA.areaI), LA.volumeI = tdc(time, LA.volumeI))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, pft.complete.non.i,"DateOfScan"), id=ID, AgeAtTest = tdc(time, AgeAtTest), FEV1 = tdc(time, FEV1), FEV1perc = tdc(time, FEV1perc),
                                      FVC = tdc(time, FVC), FVCperc = tdc(time, FVCperc), TLcoc = tdc(time, TLcoc), TLcoperc = tdc(time, TLcoperc), Kco = tdc(time, Kco), Kcoperc = tdc(time, Kcoperc))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, rhc.complete.non.i, "DateOfScan"), id=ID, AgeAtTest = tdc(time, AgeAtTest), PVR = tdc(time, PVR), mPAP = tdc(time, mPAP), PCWP = tdc(time, PCWP), CO = tdc(time, CO))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, ex.complete.non.i, "DateOfScan"), id=ID, AgeAtTest = tdc(time, AgeAtTest), distExTest = tdc(time, distExTest))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, EmPHasis10.complete.non.i, "DateOfScan"), id=ID, AgeAtTest = tdc(time, AgeAtTest), EmPHasis10 = tdc(time, EmPHasis10))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, FC.complete.non.i, "DateOfScan"), id=ID, AgeAtTest = tdc(time, AgeAtTest), FC = tdc(time, FC))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, Rx.data, "DateOfScan"), id=ID, AgeAtDiag = tdc(time, AgeAtDiag), AgeAtTest = tdc(time, AgeAtTest), 
                                      DoseEqPDE5 = tdc(time, DoseEqPDE5), DoseFctrPDE5 = tdc(time, DoseFctrPDE5), TreatedPDE5 = tdc(time, TreatedPDE5),
                                      DoseEqERA = tdc(time, DoseEqERA), DoseFctrERA = tdc(time, DoseFctrERA), TreatedERA = tdc(time, TreatedERA),
                                      Treated = tdc(time, Treated))
          
          results.olps.non.i<- tmerge(results.olps.non.i, make.tdc.matrix(d2.demo, Bloods.complete.non.i, "DateOfScan"), id=ID, AgeAtTest = tdc(time, AgeAtTest), BNP = tdc(time, BNP), Na = tdc(time, Na), K = tdc(time, K),
                                      Ur = tdc(time, Ur), Creat = tdc(time, Creat), eGFR = tdc(time, eGFR), Bili = tdc(time, Bili), ALP = tdc(time, ALP), ALT = tdc(time, ALT),
                                      Prot = tdc(time, Prot), Alb = tdc(time, Alb), CCa = tdc(time, CCa), Fe = tdc(time, Fe), TFSat = tdc(time, TFSat),
                                      TF = tdc(time, TF), TIBC = tdc(time, TIBC), T4 = tdc(time, T4))
          
          results.olps.non.i<- change.data.type(results.olps.non.i)
          results.olps.i.amelia<- JustAMELIA(results.olps.non.i, imputations, maxit, 500, logs)
          
          
          cat("\n done.")
          
          
  
          
   
