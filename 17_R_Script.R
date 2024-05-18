##### PCoA for Prokaryotic Communities ####

#Set up the folder

getwd()
setwd("C:/Users/LuisM/OneDrive - NIOO/Desktop/data")

#Importing the file 

library(readxl)
Bac_OTU_Table <- read_excel("01_prokaryotic_communities.xlsx")

#View(Bac_OTU_Table)

#Omitting missing samples
Bac_OTU_Table<-na.omit(Bac_OTU_Table)

#Row name
row.names(Bac_OTU_Table)<-Bac_OTU_Table$Sample_ID

### Active Selection

library(dplyr)

data_active_filt = Bac_OTU_Table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,8:1772]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_bac_active <- ggplot(data = dataset_active,
       
       aes(x = PCoA_active1,
           
           y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+

  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
 
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_active,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                         y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)")) +

  theme(legend.position = "none") + 
  ggtitle("Bacteria-Archaea - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot_bac_active

# Multivariate Adonis test
set.seed(123)
Adonis_Active_1 <-  adonis2(data_active_filt[,8:1772] ~ data_active_filt$Treatment + data_active_filt$Area + data_active_filt$Square + data_active_filt$Soil_class, permutations = 999, method = "bray")
print(Adonis_Active_1)


#Pairwise comparison
library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:1772],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')


### Assisted selection
library(dplyr)
data_assisted_filt = Bac_OTU_Table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:1772]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_bac_assisted <- ggplot(data = dataset_assisted,
       
       aes(x = PCoA_assisted1,
           
           y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_assisted,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Bacteria-Archaea - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))


#Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:1772] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square + data_assisted_filt$Soil_class, permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,8:1722],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

##### PCoA for Fungi Communities ####

#Importing the file 

library(readxl)
ITS_table <- read_excel("02_fungi_communities.xlsx")
#View(ITS_table)

#Omitting missing samples
ITS_table<-na.omit(ITS_table)

#Row name
row.names(ITS_table)<-ITS_table$Sample_ID

# Active selection

library(dplyr)

data_active_filt = ITS_table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,8:1391]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_fungi_active <- ggplot(data = dataset_active,
                            
                            aes(x = PCoA_active1,
                                
                                y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_active,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  
  theme(legend.position = "none") + 
  ggtitle("Fungi - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

#Multivariate Adonis

set.seed(123)
Adonis_Assisted_1 <-  adonis2(data_active_filt[,8:1391] ~ data_active_filt$Treatment  + data_active_filt$Area + data_active_filt$Square + data_active_filt$Soil_class, permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:1391],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

# Assisted selection

library(dplyr)

data_assisted_filt = ITS_table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:1391]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_fungi_assisted <- ggplot(data = dataset_assisted,
                              
                              aes(x = PCoA_assisted1,
                                  
                                  y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=16)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_assisted,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                            y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  
  theme(legend.position = "none") + 
  ggtitle("Fungi - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))


# Multivariate Adonis

set.seed(123)
Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:1391] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square + data_assisted_filt$Soil_class , permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:1391],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

#### PCoA for Protist Communities ####

#Importing the file 

library(readxl)
ITS_table <- read_excel("03_protist_communities.xlsx")
head(ITS_table)

#Omitting missing samples
ITS_table<-na.omit(ITS_table)

#Row name
row.names(ITS_table)<-ITS_table$Sample_ID

# Active selection

library(dplyr)

data_active_filt = ITS_table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,8:569]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_protist_active <- ggplot(data = dataset_active,
                              
                              aes(x = PCoA_active1,
                                  
                                  y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_active,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  
  theme(legend.position = "none") + 
  ggtitle("Protist - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot_protist_active

# Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_active_filt[,8:569] ~ data_active_filt$Treatment + data_active_filt$Area + data_active_filt$Square + data_active_filt$Soil_class, permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:569],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

### Assisted selection

library(dplyr)

data_assisted_filt = ITS_table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:569]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_protist_assisted <- ggplot(data = dataset_assisted,
                                
                                aes(x = PCoA_assisted1,
                                    
                                    y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_assisted,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                            y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  
  theme(legend.position = "none") + 
  ggtitle("Protist - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot_protist_assisted

#Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:569] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square + data_assisted_filt$Soil_class, permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,8:569],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')


#### PCoA for Prokaryotic Potential Functions (FAPROTAX) ####

#Importing the file 
library(readxl)
Bac_OTU_Table <- read_excel("04_prokaryotic_potential_functions_faprotax.xlsx")

#View(Bac_OTU_Table)

#Omitting missing samples
Bac_OTU_Table<-na.omit(Bac_OTU_Table)

#Row name
row.names(Bac_OTU_Table)<-Bac_OTU_Table$Sample_ID

# Active test

library(dplyr)

data_active_filt = Bac_OTU_Table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,8:67]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_bac_active_func <- ggplot(data = dataset_active,
                               
                               aes(x = PCoA_active1,
                                   
                                   y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_active,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Potential Functions (Bac-Arc) - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Multivariate Adonis test
set.seed(123)
Adonis_Active_1 <-  adonis2(data_active_filt[,8:67] ~ data_active_filt$Treatment + data_active_filt$Area + data_active_filt$Square + data_active_filt$Soil_class, permutations = 999, strata = data_active_filt$Area)
print(Adonis_Active_1)


#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:67],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

# Assisted selection
library(dplyr)

data_assisted_filt = Bac_OTU_Table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:67]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_bac_assisted_func <- ggplot(data = dataset_assisted,
                                 
                                 aes(x = PCoA_assisted1,
                                     
                                     y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_assisted,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                            y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Potential Functions (Bac-Arc) - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

##Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:67] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square + data_assisted_filt$Soil_class, permutations = 999, strata = data_assisted_filt$Area)
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,8:67],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

#### PCoA for Fungi Potential Functions (FUNGuild) ####

#Importing the file 

library(readxl)
Bac_OTU_Table <- read_excel("05_fungi_potential_functions_funguild.xlsx")
View(Bac_OTU_Table)

#Omitting missing samples
Bac_OTU_Table<-na.omit(Bac_OTU_Table)

#Row name
row.names(Bac_OTU_Table)<-Bac_OTU_Table$Sample_ID

# Active selection

library(dplyr)

data_active_filt = Bac_OTU_Table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,8:131]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_fungi_active_func <- ggplot(data = dataset_active,
                                 
                                 aes(x = PCoA_active1,
                                     
                                     y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_active,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Potential Functions (Fungi) - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

#Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Active_1 <-  adonis2(data_active_filt[,8:131] ~ data_active_filt$Treatment + data_active_filt$Area + data_active_filt$Square + data_active_filt$Soil_class, permutations = 999, strata = data_active_filt$Area)
print(Adonis_Active_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:131],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

# Assisted selection

library(dplyr)

data_assisted_filt = Bac_OTU_Table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:131]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_fungi_assisted_func <- ggplot(data = dataset_assisted,
                                   
                                   aes(x = PCoA_assisted1,
                                       
                                       y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_assisted,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                            y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Potential Functions (Fungi) - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))


# Multivariate Adonis test

set.seed(123)
Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:131] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square+ data_assisted_filt$Soil_class, permutations = 999, strata = data_assisted_filt$Area, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,8:131],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

#### PCoA for Protist Potential Functions (In house database) ####

#Importing the file 

library(readxl)
Bac_OTU_Table <- read_excel("06_protist_potential_functions_inhouse_database.xlsx")

#View(Bac_OTU_Table)

#Omitting missing samples
Bac_OTU_Table<-na.omit(Bac_OTU_Table)

#Row name
row.names(Bac_OTU_Table)<-Bac_OTU_Table$Sample_ID

# Active selection

library(dplyr)

data_active_filt = Bac_OTU_Table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,8:14]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_prot_func_active <- ggplot(data = dataset_active,
                                
                                aes(x = PCoA_active1,
                                    
                                    y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_active,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Potential Functions (Protist) - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Multivariate Adonis test

set.seed(123)
Adonis_Active_1 <-  adonis2(data_active_filt[,8:14] ~ data_active_filt$Treatment + data_active_filt$Area + data_active_filt$Square+ data_active_filt$Soil_class, permutations = 999, strata = data_active_filt$Area, method = "bray")
print(Adonis_Active_1)

#Pairwise comparison
library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,8:14],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')


### Assisted selection

library(dplyr)

data_assisted_filt = Bac_OTU_Table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:14]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_prot_func_assisted <- ggplot(data = dataset_assisted,
                                  
                                  aes(x = PCoA_assisted1,
                                      
                                      y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbarh(dataframe_assisted,
                 
                 mapping = aes(x = PCoA1m, y = PCoA2m,
                               
                               xmin = PCoA1m - se1,
                               
                               xmax = PCoA1m + se1),
                 
                 height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                            y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Potential Functions (Protist) - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:14] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square + data_assisted_filt$Soil_class, permutations = 999)
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,8:14],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m='fdr')

#### PCoA for Soil Characteristics ####


#Importing the file 

library(readxl)
ITS_table <- read_excel("07_soil_characteristics.xlsx")
View(ITS_table)

#Omitting missing samples
ITS_table<-na.omit(ITS_table)

#Row name
row.names(ITS_table)<-ITS_table$Sample_ID


### Active selection

set.seed(123)
library(dplyr)

data_active_filt = ITS_table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,c(8:30)]

#Creating a distance matrix
set.seed(123)
library(vegan)

dis_active <- vegdist(data_active,method='gower')

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)
#View(PCOA_active$vectors)


PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:7]

head(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_soil_active <- ggplot(data = dataset_active,
                           
                           aes(x = PCoA_active1,
                               
                               y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  # geom_errorbarh(dataframe_active,
  #                
  #                mapping = aes(x = PCoA1m, y = PCoA2m,
  #                              
  #                              xmin = PCoA1m - se1,
  #                              
  #                              xmax = PCoA1m + se1),
  #                
  #                height = 0, alpha= 0.3)  +
  # 
geom_point(data = dataframe_active, aes(x = PCoA1m,
                                        y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Soil Characteristics - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot_soil_active

# Multivariate Adonis test

set.seed(123)
Adonis_Assisted_1 <-  adonis2(data_active_filt[,c(8:30)] ~ data_active_filt$Treatment + data_active_filt$Area + data_active_filt$Square + data_active_filt$Soil_class, permutations = 999, method = 'gower')
print(Adonis_Assisted_1)


#Pairwise comparison
set.seed(123)
library(pairwiseAdonis)
pairwise.adonis(data_active_filt[,c(8:30)],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='gower',p.adjust.m='fdr')

# Assisted selection 

library(dplyr)

data_assisted_filt = ITS_table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,8:30]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="gower")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

head(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_soil_assisted <- ggplot(data = dataset_assisted,
                             
                             aes(x = PCoA_assisted1,
                                 
                                 y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  # geom_errorbarh(dataframe_assisted,
  #                
  #                mapping = aes(x = PCoA1m, y = PCoA2m,
  #                              
  #                              xmin = PCoA1m - se1,
  #                              
  #                              xmax = PCoA1m + se1),
  #                
  #                height = 0, alpha= 0.3)  +
  
geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Soil Characteristics - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot_soil_assisted

# Multivariate Adonis test
set.seed(123)
Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,8:30] ~ data_assisted_filt$Treatment + data_assisted_filt$Area + data_assisted_filt$Square + data_assisted_filt$Soil_class , permutations = 999, method = "gower")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,8:30],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='gower',p.adjust.m='fdr')


#### PCoA for Plant Communities ####

#Importing the file 

library(readxl)
ITS_table <- read_excel("08_plant_communities.xlsx")
View(ITS_table)

#Omitting missing samples
ITS_table<-na.omit(ITS_table)

# Active Selection

library(dplyr)

data_active_filt = ITS_table %>% filter(Active == "1")
dim(data_active_filt)
data_active = data_active_filt[,5:143]

#Creating a distance matrix

library(vegan)

dis_active <- vegdist(data_active,method="bray")

library(ape)

PCOA_active <- pcoa(dis_active, correction = "none", rn = NULL)
biplot.pcoa(PCOA_active, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_active)

#View(PCOA_active$values)

#View(PCOA_active$vectors)

PCOA1<- PCOA_active$values$Relative_eig[1]*100

PCOA2<- PCOA_active$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_active$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_active',1:2)

treatment<- data_active_filt[,1:4]

#View(treatment)

dataset_active<- cbind(treatment, sampledata)

head(dataset_active)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_active <-
  
  dataset_active %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_active1),
                   
                   PCoA1sd = sd(PCoA_active1),
                   
                   PCoA2m = mean(PCoA_active2),
                   
                   PCoA2sd = sd(PCoA_active2),
                   
                   se1 =  standard_error(PCoA_active1),
                   
                   se2 =  standard_error(PCoA_active2))

head(dataframe_active)


#Plotting the PCoA

dataset_active$Treatment <- as.factor(dataset_active$Treatment)

plot_plant_active <- ggplot(data = dataset_active,
                            
                            aes(x = PCoA_active1,
                                
                                y = PCoA_active2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  scale_color_manual(values=c("#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbar(dataframe_active,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              xmin = PCoA1m - se1,
                              
                              xmax = PCoA1m + se1),
                
                height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_active, aes(x = PCoA1m,
                                          y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_active$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_active$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Plant - Active Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot_plant_active

# Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_active_filt[,5:143] ~ data_active_filt$Treatment, permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)

pairwise.adonis(data_active_filt[,5:143],factors=data_active_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m = "bonferroni", perm = 999)


### Assisted Selection

library(dplyr)

data_assisted_filt = ITS_table %>% filter(Assisted == "1")
dim(data_assisted_filt)
data_assisted = data_assisted_filt[,5:143]

#Creating a distance matrix

library(vegan)

dis_assisted <- vegdist(data_assisted,method="bray")

library(ape)

PCOA_assisted <- pcoa(dis_assisted, correction = "none", rn = NULL)
biplot.pcoa(PCOA_assisted, plot.axes = c(1,2))


#Extracting the axes and creating a data frame with the axes 1 and 2
library(phyloseq)

plot_scree(PCOA_assisted)

#View(PCOA_assisted$values)

#View(PCOA_assisted$vectors)

PCOA1<- PCOA_assisted$values$Relative_eig[1]*100

PCOA2<- PCOA_assisted$values$Relative_eig[2]*100

sampledata<- as.data.frame(PCOA_assisted$vectors[,1:2])

colnames(sampledata)<- paste0('PCoA_assisted',1:2)

treatment<- data_assisted_filt[,1:7]

#View(treatment)

dataset_assisted<- cbind(treatment, sampledata)

head(dataset_assisted)

#PCoA plot

library(ggpubr)
library(tidyverse)
standard_error <- function(x) sd(x) / sqrt(length(x))

dataframe_assisted <-
  
  dataset_assisted %>%
  
  group_by(Treatment) %>%
  
  dplyr::summarise(count = n(),
                   
                   PCoA1m = mean(PCoA_assisted1),
                   
                   PCoA1sd = sd(PCoA_assisted1),
                   
                   PCoA2m = mean(PCoA_assisted2),
                   
                   PCoA2sd = sd(PCoA_assisted2),
                   
                   se1 =  standard_error(PCoA_assisted1),
                   
                   se2 =  standard_error(PCoA_assisted2))

head(dataframe_assisted)


#Plotting the PCoA

dataset_assisted$Treatment <- as.factor(dataset_assisted$Treatment)

plot_plant_assisted<- ggplot(data = dataset_assisted,
                             
                             aes(x = PCoA_assisted1,
                                 
                                 y = PCoA_assisted2)) +
  
  stat_ellipse(linetype=2, aes(group=Treatment),level = 0.95)+
  
  stat_ellipse( aes(linetype=Treatment, fill =Treatment, color = Treatment), geom = "polygon",level = 0.95,alpha = 0.1)+ #,type='t'
  
  #scale_fill_jco()+
  
  scale_color_manual(values=c("#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  scale_fill_manual(values=c("#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  
  theme(text = element_text(size=12)) +
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),
        
        panel.background = element_blank())+
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              ymin = PCoA2m - se2,
                              
                              ymax = PCoA2m + se2),
                
                width = 0, alpha= 0.3) +
  
  geom_errorbar(dataframe_assisted,
                
                mapping = aes(x = PCoA1m, y = PCoA2m,
                              
                              xmin = PCoA1m - se1,
                              
                              xmax = PCoA1m + se1),
                
                height = 0, alpha= 0.3)  +
  
  geom_point(data = dataframe_assisted, aes(x = PCoA1m,
                                            y = PCoA2m, color =Treatment),size = 5)+
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  theme(legend.title = element_blank())+
  
  labs(x = paste0("PCoA 1 (", round(PCOA_assisted$values$Relative_eig[1] * 100, 2), "%)"),
       
       y = paste0("PCoA 2 (", round(PCOA_assisted$values$Relative_eig[2] * 100, 2), "%)"))+
  theme(legend.position = "none") + 
  ggtitle("Plant - Assisted Restoration")+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Multivariate Adonis test

set.seed(123)

###With the composition table

Adonis_Assisted_1 <-  adonis2(data_assisted_filt[,5:143] ~ data_assisted_filt$Treatment, permutations = 999, method = "bray")
print(Adonis_Assisted_1)

#Pairwise comparison

library(pairwiseAdonis)
pairwise.adonis(data_assisted_filt[,5:143],factors=data_assisted_filt$Treatment,sim.function='vegdist', sim.method='bray',p.adjust.m = "bonferroni", perm = 999)


#### Diversity from the plant and soil communities - Plots ####

library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(hrbrthemes)
library(viridis)
library(ggeasy)
library(ggplot2)
library(readr)

#Plant Data

#### Importing the data

library(readr)
plant_diversity <- read_delim("09_plant_diversity.csv", 
                              delim = ",", escape_double = FALSE, trim_ws = TRUE)
View(plant_diversity)


#Filtering Active (- assisted)

Active_plants = plant_diversity %>% filter(Active == "1")

#Filtering Assisted (- active)

Assisted_plants = plant_diversity %>% filter(Assisted == "1")

#Active - Plot - Shannon Plants
box_plot_Active_shannon_plants<-ggplot(Active_plants, aes(x=Treatments, y=Shannon_H, color=Treatments)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200",  "#006A00"))+
  #scale_color_manual(values=c("#ED7D31","#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200",  "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Plant Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Active Restoration")+
  ylab("Shannon Index")

plot(box_plot_Active_shannon_plants)


#Assisted - Plot - Shannon Plants
box_plot_Assisted_shannon_plants<-ggplot(Assisted_plants, aes(x=Treatments, y=Shannon_H, color=Treatments)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  #scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Plant Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Assisted Restoration")+
  ylab("Shannon Index")

plot(box_plot_Assisted_shannon_plants)

#Importing microbiome diversity data

tudo <- read_delim("10_soil_diversity.csv", delim = ",", 
                   escape_double = FALSE, trim_ws = TRUE)
View(tudo)

#Filtering 0-10 cm

my_data = tudo %>% filter(Depth == "0_10")

#Filtering Active (- Active)

Active = my_data %>% filter(Active == "1")

Active_Bacteria = Active %>% filter(Barcode == "a_16S rRNA")

Active_Fungi = Active %>% filter(Barcode == "b_ITS")

Active_Protist = Active %>% filter(Barcode == "c_18S rRNA")


#Filtering Assisted (- Assisted)

Assisted = my_data %>% filter(Assisted == "2")

Assisted_Bacteria = Assisted %>% filter(Barcode == "a_16S rRNA")

Assisted_Fungi = Assisted %>% filter(Barcode == "b_ITS")

Assisted_Protist = Assisted %>% filter(Barcode == "c_18S rRNA")

#Active selection

#Bacteria-Archaea-Prokarya
box_plot_Active_shannon_bac<-ggplot(Active_Bacteria, aes(x=Treatment, y=Shannon_H, color=Treatment)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200",  "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Prokaryotic Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Active Restoration")+
  ylab("Shannon Index")

plot(box_plot_Active_shannon_bac)

#Fungi
box_plot_Active_shannon_fungi<-ggplot(Active_Fungi, aes(x=Treatment, y=Shannon_H, color=Treatment)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200",  "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Fungi Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Active Restoration")+
  ylab("Shannon Index")

plot(box_plot_Active_shannon_fungi)

#Protist
box_plot_Active_shannon_protist<-ggplot(Active_Protist, aes(x=Treatment, y=Shannon_H, color=Treatment)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200",  "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Protist Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Active Restoration")+
  ylab("Shannon Index")

plot(box_plot_Active_shannon_protist)

#Assisted selection

#Bacteria-Archaea-Prokarya
box_plot_Assisted_shannon_bac<-ggplot(Assisted_Bacteria, aes(x=Treatment, y=Shannon_H, color=Treatment)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Prokaryotic Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Assisted Restoration")+
  ylab("Shannon Index")

plot(box_plot_Assisted_shannon_bac)

#Fungi
box_plot_Assisted_shannon_fungi<-ggplot(Assisted_Fungi, aes(x=Treatment, y=Shannon_H, color=Treatment)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Fungi Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Assisted Restoration")+
  ylab("Shannon Index")

plot(box_plot_Assisted_shannon_fungi)

#Protist
box_plot_Assisted_shannon_protist<-ggplot(Assisted_Protist, aes(x=Treatment, y=Shannon_H, color=Treatment)) +
  geom_violin(size= 1.0)+
  theme_gray()+
  geom_jitter(shape=16, position=position_jitter(0.1), size = 0.4)+
  scale_color_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  ggtitle("Protist Diversity")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+ #element_text(angle = 45, hjust=1)
  xlab("Assisted Restoration")+
  ylab("Shannon Index")

plot(box_plot_Assisted_shannon_protist)


#### Diversity from the plant and soil communities - Statistics ####

#Activating the package

library(agricolae)
library(dplyr)
library(readr)

#Importing the data_set using the readr - .csv

Tudo_2 <- read_delim("11_statistics_diversity.csv")
View(Tudo_2)
my_data = Tudo_2

# #Filtering the depth

#Filtering the barcode
bacteria = my_data %>% filter(Barcode == "a_16S rRNA")
fungi = my_data %>% filter(Barcode =="b_ITS")
protist = my_data %>% filter(Barcode =="c_18S rRNA")

#Filtering Restorations

active_bac = bacteria %>% filter(Active == "1")
active_fungi = fungi %>% filter(Active == "1")
active_protist = protist %>% filter(Active == "1")

assisted_bac = bacteria %>% filter(Assisted == "2")
assisted_fungi = fungi %>% filter(Assisted == "2")
assisted_protist = protist %>% filter(Assisted == "2")

#Activating packages
set.seed(123)
library(lme4)

#Active Prokaryotic Diversity
active_bac_model_div <- lmer(Shannon_H ~ Treatment + (1 | Square) + (1 | Soil_class), data=active_bac)
summary(active_bac_model_div)  

# summary() of mixed model object
anova(active_bac_model_div)

library(jtools)
summ(active_bac_model_div)
library(lmerTest)
ranova(active_bac_model_div)

#Post-hoc test
library("multcomp")
stat <- glht(active_bac_model_div, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(stat)
tuk.cld <- cld(stat ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#Active Fungi Diversity
active_fungi_model_div <- lmer(Shannon_H ~ Treatment + (1 | Square) + (1 | Soil_class), data=active_fungi)

summary(active_fungi_model_div)  

# summary() of mixed model object
anova(active_fungi_model_div)

library(jtools)
summ(active_fungi_model_div)

library(lmerTest)
ranova(active_fungi_model_div)

#Post-hoc test
library("multcomp")
stat <- glht(active_fungi_model_div, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(stat)
tuk.cld <- cld(stat ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#Active Protist Diversity

active_protist_model_div <- lmer(Shannon_H ~ Treatment + (1 | Square) + (1 | Soil_class), data=active_protist)

summary(active_protist_model_div)  

# summary() of mixed model object
anova(active_protist_model_div)

library(jtools)
summ(active_protist_model_div)

library(lmerTest)
ranova(active_protist_model_div)

#Post-hoc test
library("multcomp")
stat <- glht(active_protist_model_div, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(stat)
tuk.cld <- cld(stat ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)



#Assisted Prokaryotic Diversity

assisted_bac_model_div <- lmer(Shannon_H ~ Treatment + (1 | Square) + (1 | Soil_class), data=assisted_bac)

summary(assisted_bac_model_div)  

# summary() of mixed model object
anova(assisted_bac_model_div)

library(jtools)
summ(assisted_bac_model_div)

library(lmerTest)
ranova(assisted_bac_model_div)

#Post-hoc test
library("multcomp")
stat <- glht(assisted_bac_model_div, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(stat)
tuk.cld <- cld(stat ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)


#Assisted Fungi Diversity

assisted_fungi_model_div <- lmer(Shannon_H ~ Treatment + (1 | Square) + (1 | Soil_class), data=assisted_fungi)

summary(assisted_fungi_model_div)  

# summary() of mixed model object
anova(assisted_fungi_model_div)
library(jtools)
summ(assisted_fungi_model_div)

library(lmerTest)
ranova(assisted_fungi_model_div)

#Post-hoc test
library("multcomp")
stat <- glht(assisted_fungi_model_div, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(stat)
tuk.cld <- cld(stat ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#Assisted Protist Diversity

assisted_protist_model_div <- lmer(Shannon_H ~ Treatment + (1 | Square) + (1 | Soil_class), data=assisted_protist)

summary(assisted_protist_model_div)  

# summary() of mixed model object
anova(assisted_protist_model_div)

library(jtools)
summ(assisted_protist_model_div)

library(lmerTest)
ranova(assisted_protist_model_div)

#Post-hoc test
library("multcomp")
stat <- glht(assisted_protist_model_div, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(stat)
tuk.cld <- cld(stat ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)


#Plant Ddiversity

library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(hrbrthemes)
library(viridis)
library(ggeasy)
library(ggplot2)
library(readr)

#### Importing the data

library(readr)
plant_diversity <- read_csv("09_plant_diversity.csv")
View(plant_diversity)


#Filtering Active (- assisted)

Active_plants = plant_diversity %>% filter(Active == "1")

#Filtering Assisted (- active)

Assisted_plants = plant_diversity %>% filter(Assisted == "1")


# Plant Active Diversity

# Make sure Treatments is a factor
Active_plants$Treatments <- as.factor(Active_plants$Treatments)

# Fit the ANOVA model
active_protist_model_div <- aov(Shannon_H ~ Treatments, data=Active_plants)

# Display the summary of the ANOVA model
summary(active_protist_model_div)

# Fit the ANOVA model
active_protist_model_div <- lm(Shannon_H ~ Treatments, data=Active_plants)

# Display the summary of the ANOVA model
summary(active_protist_model_div)


# Post-hoc test
library("multcomp")
# Specify that Treatments is a factor
Active_plants$Treatments <- as.factor(Active_plants$Treatments)

# Perform post-hoc testing with Bonferroni adjustment
tuk.cld <- cld( 
  glht(active_protist_model_div, linfct = mcp(Treatments = "Tukey")), 
  alpha = 0.05
)
tuk.cld

plot(tuk.cld)


# Plant Assisted Diversity

# Make sure Treatments is a factor
Assisted_plants$Treatments <- as.factor(Assisted_plants$Treatments)

# Fit the ANOVA model
Assisted_protist_model_div <- aov(Shannon_H ~ Treatments, data=Assisted_plants)

# Display the summary of the ANOVA model
summary(Assisted_protist_model_div)

Assisted_protist_model_div <- lm(Shannon_H ~ Treatments, data=Assisted_plants)

summary(Assisted_protist_model_div)

# Post-hoc test
library("multcomp")
# Specify that Treatments is a factor
Assisted_plants$Treatments <- as.factor(Assisted_plants$Treatments)

# Perform post-hoc testing with Bonferroni adjustment
tuk.cld <- cld( 
  glht(Assisted_protist_model_div, linfct = mcp(Treatments = "Tukey")), 
  alpha = 0.05
)
tuk.cld

plot(tuk.cld)

#### C stocks from plant can C stocks and isotopes from soil ####

#Packages

library(ggplot2)
library(readr)
library(dplyr)
library(readxl)
library(reshape2)
library(viridis)
library(dplyr)

#Import the excel file
my_data <- read_excel("12_carbon_data_soil.xlsx")
View(my_data)

#Filtering Assisted (- assisted)
c_table_filt_active <- my_data %>% filter(Active == "1")
dim(c_table_filt_active)
c_table_active = c_table_filt_active[,6:9]
c_data_active = c_table_active[,-c(4)]

#Filtering Assisted (- active)
c_table_filt_assisted <- my_data %>% filter(Assisted == "1")
dim(c_table_filt_assisted)
c_table_assisted = c_table_filt_assisted[,6:9]
c_data_assisted = c_table_assisted[,-c(4)]

#Set up the treatment
c_data_active$Treatment <- as.factor(c_data_active$Treatment)
c_data_assisted$Treatment <- as.factor(c_data_assisted$Treatment)

#Functions to calculate the metrics to error bar

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# Setup the treatments

pcm_active = melt(c_data_active, id = c("Treatment"))

c_data_active_2 <- data_summary(pcm_active, varname="value", 
                                groupnames=c("Treatment", "variable"))

pcm_assisted = melt(c_data_assisted, id = c("Treatment"))

c_data_assisted_2 <- data_summary(pcm_assisted, varname="value", 
                                  groupnames=c("Treatment", "variable"))

# Adding the factors to the new table

c_data_active_2$Treatment=as.factor(c_data_active_2$Treatment)
head(c_data_active_2)

c_data_assisted_2$Treatment=as.factor(c_data_assisted_2$Treatment)
head(c_data_assisted_2)

soil_carbon_active_total <- data_summary(pcm_active, varname="value", 
                                         groupnames=c("Treatment"))
soil_carbon_assisted_total <- data_summary(pcm_assisted, varname="value", 
                                           groupnames=c("Treatment"))


#Carbon stocks together
#Adapted from https://stackoverflow.com/questions/39723606/add-error-bars-to-stacked-bar-plot-in-ggplot2-r-solved

#Active
# dev.off()

c_data_active_2$y_pos = NA
c_data_active_2$y_pos[c_data_active_2$variable =="C_C3"] = c_data_active_2$value[c_data_active_2$variable =="C_C3"]
c_data_active_2$y_pos[c_data_active_2$variable =="C_C4"] = c_data_active_2$value[c_data_active_2$variable =="C_C3"]+
  c_data_active_2$value[c_data_active_2$variable == "C_C4"]

c_data_active_2$C_type = NA
c_data_active_2$C_type[c_data_active_2$variable =="C_C3"] = "C3"
c_data_active_2$C_type[c_data_active_2$variable =="C_C4"] = "C4"

#Plotting active 

#Plotting C - Active Stack

c_active_stack <- ggplot(c_data_active_2, aes(fill=Treatment, y=value, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position="stack") +
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=y_pos-sd, ymax=y_pos+sd), width=0.2, position = "identity")+
  labs(title="Soil C stocks", x="", y = "Mg C ha-1")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(c_active_stack)


#assisted

c_data_assisted_2$y_pos = NA
c_data_assisted_2$y_pos[c_data_assisted_2$variable =="C_C3"] = c_data_assisted_2$value[c_data_assisted_2$variable =="C_C3"]
c_data_assisted_2$y_pos[c_data_assisted_2$variable =="C_C4"] = c_data_assisted_2$value[c_data_assisted_2$variable =="C_C3"]+
  c_data_assisted_2$value[c_data_assisted_2$variable == "C_C4"]

c_data_assisted_2$C_type = NA
c_data_assisted_2$C_type[c_data_assisted_2$variable =="C_C3"] = "C3"
c_data_assisted_2$C_type[c_data_assisted_2$variable =="C_C4"] = "C4"

#Plotting assisted 

#Plotting C - assisted Stack

c_assisted_stack <- ggplot(c_data_assisted_2, aes(fill=Treatment, y=value, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position="stack") +
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=y_pos-sd, ymax=y_pos+sd), width=0.2, position = "identity")+
  labs(title="Soil C stocks", x="", y = "Mg C ha-1")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(c_assisted_stack)

# Soil Carbon Isotopes
set.seed(123)
library(lme4)
library(lmerTest)
library(jtools)
library("multcomp")

#For Active data
#C3

#Model

model_c_soil_c3_active<- lmer(C_C3 ~ Treatment + (1 | Square) + (1| Soi_class), data=c_table_filt_active)
summary(model_c_soil_c3_active)
print(model_c_soil_c3_active, correlation = TRUE)

# summary() of mixed model object
anova(model_c_soil_c3_active, ddf="Satterthwaite")
summ(model_c_soil_c3_active)
ranova(model_c_soil_c3_active)

#Post-hoc test

post_hoc_test <- glht(model_c_soil_c3_active, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#C4

#Model

model_c_soil_C4_active<- lmer(C_C4 ~ Treatment  + (1 | Square) + (1| Soi_class), data=c_table_filt_active)
summary(model_c_soil_C4_active)
print(model_c_soil_C4_active, correlation = TRUE)

# summary() of mixed model object
anova(model_c_soil_C4_active, ddf="Satterthwaite")
summ(model_c_soil_C4_active)
ranova(model_c_soil_C4_active)

#Post-hoc test

post_hoc_test <- glht(model_c_soil_C4_active, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#C-Total

#Model
model_c_soil_total_active<- lmer(C_Total ~ Treatment + (1| Square) + (1| Soi_class), data=c_table_filt_active)
summary(model_c_soil_total_active)
print(model_c_soil_total_active, correlation = TRUE)

# summary() of mixed model object
anova(model_c_soil_total_active, ddf="Satterthwaite")
summ(model_c_soil_total_active)
ranova(model_c_soil_total_active)

#Post-hoc test

post_hoc_test <- glht(model_c_soil_total_active, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#For Assisted dataset

#C3

#Model

model_c_soil_c3_assisted<- lmer(C_C3 ~ Treatment + (1 | Square) + (1| Soi_class), data=c_table_filt_assisted)
summary(model_c_soil_c3_assisted)
print(model_c_soil_c3_assisted, correlation = TRUE)

# summary() of mixed model object
anova(model_c_soil_c3_assisted, ddf="Satterthwaite")
summ(model_c_soil_c3_assisted)
ranova(model_c_soil_c3_assisted)

#Post-hoc test

post_hoc_test <- glht(model_c_soil_c3_assisted, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#C4

#Model

model_c_soil_C4_assisted<- lmer(C_C4 ~ Treatment + (1 | Square) + (1| Soi_class), data=c_table_filt_assisted)
summary(model_c_soil_C4_assisted)
print(model_c_soil_C4_assisted, correlation = TRUE)

# summary() of mixed model object
anova(model_c_soil_C4_assisted, ddf="Satterthwaite")
summ(model_c_soil_C4_assisted)
ranova(model_c_soil_C4_assisted)

#Post-hoc test

post_hoc_test <- glht(model_c_soil_C4_assisted, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#C-Total

#Model
model_c_soil_total_assisted<- lmer(C_Total ~ Treatment + (1 | Square) + (1| Soi_class), data=c_table_filt_assisted)
summary(model_c_soil_total_assisted)
print(model_c_soil_total_assisted, correlation = TRUE)

# summary() of mixed model object
anova(model_c_soil_total_assisted, ddf="Satterthwaite")
# anova(model_c_soil_total_assisted,ddf="lme4")
# anova(model_c_soil_total_assisted,ddf="Kenward-Roger")
summ(model_c_soil_total_assisted)
ranova(model_c_soil_total_assisted)

#Post-hoc test

post_hoc_test <- glht(model_c_soil_total_assisted, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#Plant Carbon
#Importing the data

library(readxl)
plant_c <- read_excel("13_carbon_data_plant.xlsx")
View(plant_c)


#Filtering plant -  Active

c_plant_active <- plant_c %>% filter(Active == "1")
#Average and standard deviation to the plot
c_plant_active_summary <- data_summary(c_plant_active, varname="Plant_C",groupnames=c("Treatment"))
#Adding the factor
c_plant_active_summary$Treatment=as.factor(c_plant_active_summary$Treatment)

#Active Barplot for plants

plants_c_active<- ggplot(c_plant_active_summary , aes(fill=Treatment, y=Plant_C, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_fill_manual(values=c("#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=Plant_C-sd, ymax=Plant_C+sd), width=.2,
                position=position_dodge(.9))+
  labs(title="Plant C stocks", x="", y = "Mg C ha-1")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(plants_c_active)


#Filtering plant -  Assisted

c_plant_assisted <- plant_c %>% filter(Assisted == "1")
#Average and standard deviation to the plot
c_plant_assisted_summary <- data_summary(c_plant_assisted, varname="Plant_C",groupnames=c("Treatment"))
#Adding the factor
c_plant_assisted_summary$Treatment=as.factor(c_plant_assisted_summary$Treatment)

#Assisted Barplot for plants

plants_c_assisted<- ggplot(c_plant_assisted_summary , aes(fill=Treatment, y=Plant_C, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_fill_manual(values=c("#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=Plant_C-sd, ymax=Plant_C+sd), width=.2,
                position=position_dodge(.9))+
  labs(title="Plant C stocks", x="", y = "Mg C ha-1")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(plants_c_assisted)


##Plant Carbon Statistics


#Importing the data

library(readxl)
plant_c <- read_excel("13_carbon_data_plant.xlsx")
View(plant_c)

#Filtering plant -  Active
c_plant_active <- plant_c %>% filter(Active == "1")

#Stat
set.seed(123)
model_c_plant_active <- aov(Plant_C ~ Treatment, data=c_plant_active)
model_c_plant_active_2 <- lm(Plant_C ~ Treatment, data=c_plant_active)
summary(model_c_plant_active)
anova(model_c_plant_active)
summary(model_c_plant_active_2) # extraction the R2

#emmeans
library(multcompView)
library(lsmeans)
library(emmeans)
inter.test2 <- emmeans(model_c_plant_active, "Treatment", type = "response", adjusted = "bonferroni")
multcomp::cld(inter.test2, Letter="abcdefg")


#Filtering plant -  assisted

c_plant_assisted <- plant_c %>% filter(Assisted == "1")

#Stat
set.seed(123)
model_c_plant_assisted <- aov(Plant_C ~ Treatment, data=c_plant_assisted)
model_c_plant_assisted_2 <- lm(Plant_C ~ Treatment, data=c_plant_assisted)
summary(model_c_plant_assisted)
anova(model_c_plant_assisted)
summary(model_c_plant_assisted_2) # extraction the R2

#emmeans
library(multcompView)
library(lsmeans)
library(emmeans)
inter.test2 <- emmeans(model_c_plant_assisted, "Treatment", type = "response", adjusted = "bonferroni")
multcomp::cld(inter.test2, Letter="abcdefg")

#### Greenhouse gas plots ####

#Packages

library(ggplot2)
library(readr)
library(dplyr)
library(customerChurnUU )

#Import the .csv table
my_data <- read_delim("14_greenhouse_gases_data.csv", delim = ",", 
                      escape_double = FALSE, na = "NA", trim_ws = TRUE)
View(my_data)


#Filtering Assisted (- assisted)
c_table_filt_active <- my_data %>% filter(Active == "1")
dim(c_table_filt_active)
c_table_active = c_table_filt_active[,5:10]
c_data_active = c_table_active[,-c(3,4,5)]

#Filtering Assisted (- active)
c_table_filt_assisted <- my_data %>% filter(Assisted == "1")
dim(c_table_filt_assisted)
c_table_assisted = c_table_filt_assisted[,5:10]
c_data_assisted = c_table_assisted[,-c(3,4,5)]

#Set up the treatment
c_data_active$Treatment <- as.factor(c_data_active$Treatment)
c_data_assisted$Treatment <- as.factor(c_data_assisted$Treatment)


#CH4 - Methane

#Active

c_table_active = c_table_filt_active[,5:10]
#only 0-10 cm layer
ghg_table_active_filt <- c_table_active %>% filter(Depth == "0-10") 
#Removed columns not useful for the analysis
CH4_active <- ghg_table_active_filt[,-c(1,4,5,6)] 
#Set up the treatment
CH4_active$Treatment <- as.factor(CH4_active$Treatment) 
#Average and standard deviation to the plot
CH4_active_2 <- data_summary(CH4_active, varname="CH4",groupnames=c("Treatment"))
# Convert dose to a factor variable
CH4_active_2$Treatment=as.factor(CH4_active_2$Treatment)
head(CH4_active_2)

#Plotting CH4 - Active
CH4_active_plot<- ggplot(CH4_active_2, aes(fill=Treatment, y=CH4, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=CH4-sd, ymax=CH4+sd), width=.2,
                position=position_dodge(.9))+
  labs(title="CH4 emission/comsuption", x="", y = "(µg C/m2/h)")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(CH4_active_plot)

#Assisted

c_table_assisted = c_table_filt_assisted[,5:10]
#only 0-10 cm layer
ghg_table_assisted_filt <- c_table_assisted %>% filter(Depth == "0-10") 
#Removed columns not useful for the analysis
CH4_assisted <- ghg_table_assisted_filt[,-c(1,4,5,6)] 
#Set up the treatment
CH4_assisted$Treatment <- as.factor(CH4_assisted$Treatment) 
#Average and standard deviation to the plot
CH4_assisted_2 <- data_summary(CH4_assisted, varname="CH4",groupnames=c("Treatment"))
# Convert dose to a factor variable
CH4_assisted_2$Treatment=as.factor(CH4_assisted_2$Treatment)
head(CH4_assisted_2)

#Plotting CH4 - Active
CH4_assisted_plot<- ggplot(CH4_assisted_2, aes(fill=Treatment, y=CH4, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=CH4-sd, ymax=CH4+sd), width=.2,
                position=position_dodge(.9))+
  labs(title="CH4 emission/comsuption", x="", y = "(µg C/m2/h)")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(CH4_assisted_plot)


#CO2 - Dioxide of carbon

#Active

c_table_active = c_table_filt_active[,5:10]
#only 0-10 cm layer
ghg_table_active_filt <- c_table_active %>% filter(Depth == "0-10") 
#Removed columns not useful for the analysis
CO2_active <- ghg_table_active_filt[,-c(1,3,5,6)] 
#Set up the treatment
CO2_active$Treatment <- as.factor(CO2_active$Treatment) 
#Average and standard deviation to the plot
CO2_active_2 <- data_summary(CO2_active, varname="CO2",groupnames=c("Treatment"))
# Convert dose to a factor variable
CO2_active_2$Treatment=as.factor(CO2_active_2$Treatment)
head(CO2_active_2)

#Plotting CO2 - Active
CO2_active_plot<- ggplot(CO2_active_2, aes(fill=Treatment, y=CO2, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_fill_manual(values=c("#ED7D31", "#9BCDFF", "#3383CB", "#005AB4", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=CO2-sd, ymax=CO2+sd), width=.2,
                position=position_dodge(.9))+
  labs(title="CO2 emission", x="", y = "(mg C/m2/h)")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(CO2_active_plot)


#Assisted

c_table_assisted = c_table_filt_assisted[,5:10]
#only 0-10 cm layer
ghg_table_assisted_filt <- c_table_assisted %>% filter(Depth == "0-10") 
#Removed columns not useful for the analysis
CO2_assisted <- ghg_table_assisted_filt[,-c(1,3,5,6)] 
#Set up the treatment
CO2_assisted$Treatment <- as.factor(CO2_assisted$Treatment) 
#Average and standard deviation to the plot
CO2_assisted_2 <- data_summary(CO2_assisted, varname="CO2",groupnames=c("Treatment"))
# Convert dose to a factor variable
CO2_assisted_2$Treatment=as.factor(CO2_assisted_2$Treatment)
head(CO2_assisted_2)

#Plotting CO2 - assisted
CO2_assisted_plot<- ggplot(CO2_assisted_2, aes(fill=Treatment, y=CO2, x=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_fill_manual(values=c("#fec44f", "#CFC9FF","#A375FF", "#7731EE", "#52C088", "#00A200", "#006A00"))+
  geom_errorbar(aes(ymin=CO2-sd, ymax=CO2+sd), width=.2,
                position=position_dodge(.9))+
  labs(title="CO2 emission", x="", y = "(mg C/m2/h)")+
  theme(axis.text.x = element_blank(), legend.position = "none")

plot(CO2_assisted_plot)


#### Greenhouse gas statistics ####

#Importing the data set for statistics

my_data_stat <- read_delim("15_greenhouse_gases_statistics.csv", delim = ",", 
                           escape_double = FALSE, na = "NA", trim_ws = TRUE)
View(my_data_stat)

#Filtering Active
active_data_stat <- my_data_stat%>% filter(Active == "1")

#filtering Assisted
assisted_data_stat <- my_data_stat%>% filter(Assisted == "1")

#Package
set.seed(123)
library(lme4)
library(lmerTest)
library(jtools)
library("multcomp")

#Active data 

#Filtering the 0-10 cm layer
active_0_10_layer <- active_data_stat %>% filter(Depth == "0-10") 

#CH4

#Model
model_active_CH4 <- lmer(CH4 ~ Treatment + (1 | Square) + (1 | Soil_class), data=active_0_10_layer)
summary(model_active_CH4)          
print(model_active_CH4, correlation = TRUE)

# summary() of mixed model object
anova(model_active_CH4, ddf="Satterthwaite")
summ(model_active_CH4)
ranova(model_active_CH4)

#Post-hoc test

post_hoc_test <- glht(model_active_CH4, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#CO2

#Model
model_active_CO2 <- lmer(CO2 ~ Treatment + (1 | Square) + (1 | Soil_class), data=active_0_10_layer)
summary(model_active_CO2)          
print(model_active_CO2, correlation = TRUE)

# summary() of mixed model object
anova(model_active_CO2, ddf="Satterthwaite")
summ(model_active_CO2)
ranova(model_active_CO2)

#Post-hoc test

post_hoc_test <- glht(model_active_CO2, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#Assisted data

#Filtering the 0-10 cm layer
assisted_0_10_layer <- assisted_data_stat %>% filter(Depth == "0-10")

#CH4

#Model
model_assisted_CH4 <- lmer(CH4 ~ Treatment + (1 | Square) + (1 | Soil_class), data=assisted_0_10_layer)
summary(model_assisted_CH4)
print(model_assisted_CH4, correlation = TRUE)

# summary() of mixed model object
anova(model_assisted_CH4, ddf="Satterthwaite")
summ(model_assisted_CH4)
ranova(model_assisted_CH4)

#Post-hoc test

post_hoc_test <- glht(model_assisted_CH4, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)

#CO2
#Model
model_assisted_CO2 <- lmer(CO2 ~ Treatment + (1 | Square) + (1 | Soil_class), data=assisted_0_10_layer)
summary(model_assisted_CO2)
print(model_assisted_CO2, correlation = TRUE)

# summary() of mixed model object
anova(model_assisted_CO2, ddf="Satterthwaite")
summ(model_assisted_CO2)
ranova(model_assisted_CO2)

#Post-hoc test

post_hoc_test <- glht(model_assisted_CO2, linfct = mcp(Treatment = "Tukey"), adjusted = "bonferroni")
summary(post_hoc_test)
tuk.cld <- cld(post_hoc_test ,level = 0.05, Letter = "abcdefg")
tuk.cld

#Plotting with the letters
plot(tuk.cld)
