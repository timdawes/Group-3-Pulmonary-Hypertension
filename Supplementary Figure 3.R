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
# Supplementary Figure 3: Sankey Plot for treatment options

library(dplyr)
library(networkD3)
library(tidyr)

# Add a 'group' column to the nodes data frame:
my_color <- 'd3.scaleOrdinal() .domain(["", "Treatment","No Treatment", "Sildenafil", "Tadalafil", "Ambrisentan","Bosentan", "Epoprostenol", "Macitentan"]) .range(["#69b3a2", "#69b3a2", "#69b3a2", "#69b3a2", "#69b3a2","steelblue", "steelblue", "steelblue", "steelblue"])'

# Give a color for each group:
links <- data.frame(
  source2=c("",            "",
            "Treatment",   "Treatment", "Treatment",
            "1st line Sildenafil", "1st line Sildenafil", "1st line Sildenafil", "1st line Sildenafil", "1st line Sildenafil",
            "2nd line Tadalafil","1st line Tadalafil", "2nd line Ambrisentan",
            "3rd line Ambrisentan"), 
  target2=c("No Treatment","Treatment",
            "1st line Sildenafil", "1st line Tadalafil", "1st line Bosentan",
            "2nd line Tadalafil", "2nd line Ambrisentan","2nd line Bosentan", "2nd line Epoprostenol","2nd line Macitentan",
            "3rd line Ambrisentan", "2nd line Bosentan", "3rd line Bosentan",
            "4th line Bosentan"),
  value=c(102,82,71,1,10,6,1,5,1,1,2,1,1,1))



# Node data frame
      nodes <- data.frame(
        name=c(as.character(links$source2), 
               as.character(links$target2)) %>% unique())

# Connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source2, nodes$name)-1 
links$IDtarget <- match(links$target2, nodes$name)-1

# Make the Network
      p <- sankeyNetwork(Links = links, Nodes = nodes,
                         Source = "IDsource", Target = "IDtarget",
                         Value = "value", NodeID = "name", nodeWidth = 30,
                         sinksRight=FALSE, colourScale=my_color)
      p

(links$value) / .82
      