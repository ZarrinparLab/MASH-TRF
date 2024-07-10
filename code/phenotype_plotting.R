setwd("/mnt/zarrinpar/scratch/sfloresr/NASH_KF")

library(tidyverse)
library(data.table)
library(ggvenn)
library(gplots)
library(ggpubfigs)
library(viridis)
library(RColorBrewer)
library(BioVenn) 
library(VennDiagram)
library("qiime2R")
library(ggpubr)

#plotting the phenotypic data from Kelly
#############################################################
#food intake

food_NASH <- read.csv("files_from_KF/STAM_TRF/Mouse_and_food_weights/STAM_TRF_mouse_food_weights_NASH_separated_mice_grouped.csv",  header=TRUE)
food_NASH$condition[is.na(food_NASH$condition)] <- "NA" #rename NA group

food_NASH <- group_by(food_NASH, week, condition) #group by week & condition 

food_NASH_mean<- summarise(food_NASH,avg_cumulative_kcal_per_mouse= mean(avg_cumulative_kcal_per_mouse)) #average weight per group per week
food_NASH_SEM<- summarise(food_NASH, SEM= sd(avg_cumulative_kcal_per_mouse)/sqrt(length(avg_cumulative_kcal_per_mouse))) #SEM per group per week
food_NASH_summarized <- merge(food_NASH_mean, food_NASH_SEM) %>% # merge mean weight & SEM tables
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

#NASH cohort cumulative food intake
NASH_food <- ggplot(food_NASH_summarized, aes(x= week, y= avg_cumulative_kcal_per_mouse))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=avg_cumulative_kcal_per_mouse-SEM, ymax=avg_cumulative_kcal_per_mouse+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition), colour="black",pch=21, size=2) +
  labs(title="Average Cumulative Food\nIntake vs Time") +
  ylab("Food Intake (kcal/mouse)")+
  xlab("Week (post wean)")+
  theme_classic()+
  scale_x_continuous(breaks = seq(1, 7, by = 1))+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top")

ggsave("phenotype/cumul_food_intake.pdf", plot=NASH_food,height=3.2, width=3)


# NASH food consumption stats
for (i in 1:7){
  nash_food_t_test <- pairwise.wilcox.test(x= subset(food_NASH, week==i)$avg_cumulative_kcal_per_mouse, 
                                           g=subset(food_NASH, week==i)$condition,
                                           p.adjust.method = "fdr")
  print(i)
  print(nash_food_t_test)
}

#############################################################
#weight

weight_NASH <- read.csv("files_from_KF/STAM_TRF/Mouse_and_food_weights/STAM_TRF_mouse_weights_NASH.csv",  header=TRUE)
weight_NASH$condition[is.na(weight_NASH$condition)] <- "NA" #rename NA group

weight_NASH <- group_by(weight_NASH, week, condition) #group by week & condition 

weight_NASH_mean<- summarise(weight_NASH,weight_g= mean(weight_g)) #average weight per group per week
weight_NASH_SEM<- summarise(weight_NASH, SEM= sd(weight_g)/sqrt(length(weight_g))) #SEM per group per week
weight_NASH_summarized <- merge(weight_NASH_mean, weight_NASH_SEM) %>% # merge mean weight & SEM tables
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

#NASH cohort weights 
NASH_weights_plot <- ggplot(weight_NASH_summarized, aes(x= week, y= weight_g))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=weight_g-SEM, ymax=weight_g+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition), colour="black",pch=21, size=2) +
  labs(title="Average Mouse\nWeight vs Time") +
  ylab("Weight (g)")+
  xlab("Week (after wean)")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 6, by = 1))+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") 

ggsave("phenotype/avg_mouse_weight.pdf", plot=NASH_weights_plot,height=3.2, width=3)

for (i in 0:6){
  nash_weight_t_test <- pairwise.wilcox.test(x= subset(weight_NASH, week==i)$weight_g, 
                                             g=subset(weight_NASH, week==i)$condition,
                                             p.adjust.method = "fdr")
  print(i)
  print(nash_weight_t_test)
}

#############################################################
#insulin

df_insulin<- data.frame(read.csv("files_from_KF/STAM_TRF/ELISA_Results/stam_trf_insulin_measurements.csv", header=TRUE))%>%
  filter(cohort=="NASH")%>%
  mutate(Condition=ifelse(Condition=="NA_","NA",Condition))%>%
  mutate(Condition=factor(Condition,levels=c("NA","FA","FT")),
         zt=factor(zt,levels=c("ZT1","ZT13")))

#without control samples
insulin <- ggplot(df_insulin, aes(x=Condition, y=Concentration_ng_ml, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  labs(title="Serum Insulin Concentration, 12wk") +
  ylab("Serum Insulin Concentration (ng/ml)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 

ggsave("phenotype/serum_insulin_conc.pdf", plot=insulin,height=4, width=4)

##by ZT
insulin <- ggplot(df_insulin, aes(x=Condition, y=Concentration_ng_ml, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  facet_grid(~zt)+
  labs(title="Serum Insulin Concentration, 12wk") +
  ylab("Serum Insulin Concentration (ng/ml)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 

ggsave("phenotype/serum_insulin_conc_byZT.pdf", plot=insulin,height=4, width=6)


#Mann Whitney ZT1 vs ZT13 per condition
W <- c()
pval <- c()
Conditions <- c("NA", "FA", "FT")

for (i in 1:length(Conditions)){
  df_loop <- df_insulin%>% filter(Condition==Conditions[i])
  res <- wilcox.test(data=df_loop, Concentration_ng_ml ~ zt)
  W <- c(W, res$statistic)
  pval <- c(pval, res$p.value)
}
mann_whitney_insulin_zt1_V_zt13 <- data.frame('Condition'=Conditions, 'W stat'= W, "p value"=pval)
mann_whitney_insulin_zt1_V_zt13 #diff in FT for ZT1 vs ZT13

#Kruskal Wallis FA vs FT vs NA
chisq <- c()
df <- c()
pval <- c()
comparison <- c("FAvFT", "NAvFT","NAvFA")
Conditions <- c("NA", "FA", "FT")

for (i in 1:length(comparison)){
  df_loop <- subset(df_insulin, Condition!=Conditions[i])
  res <- kruskal.test(data=df_loop, Concentration_ng_ml ~ Condition)
  chisq <- c(chisq, res$statistic)
  df <- c(df, res$parameter)
  pval <- c(pval, res$p.value)
}
kruskal_wallis_insulin_conditions <- data.frame('Comparison'=comparison, 'chi_squared'= chisq, 'df'=df, "p value"=pval)
kruskal_wallis_insulin_conditions #no diff b/w conditions when combining ZT1 ans ZT13

pairwise.wilcox.test(df_insulin$Concentration_ng_ml, df_insulin$Condition,
                     p.adjust.method = "fdr")

# NA   FA  
# FA 0.34 -   
#   FT 0.56 0.34

#Kruskal Wallis FA vs FT vs NA, split by ZT
chisq <- c()
df <- c()
pval <- c()
comparison <- c("ZT1_FAvFT", "ZT1_NAvFT","ZT1_NAvFA","ZT13_FAvFT", "ZT13_NAvFT","ZT13_NAvFA")
Conditions <- c("NA", "FA", "FT","NA", "FA", "FT")
ZT <- c("ZT1","ZT1","ZT1","ZT13","ZT13","ZT13")

for (i in 1:length(comparison)){
  df_loop <- subset(df_insulin, zt==ZT[i] & Condition!=Conditions[i])
  res <- kruskal.test(data=df_loop, Concentration_ng_ml ~ Condition)
  chisq <- c(chisq, res$statistic)
  df <- c(df, res$parameter)
  pval <- c(pval, res$p.value)
}
kruskal_wallis_insulin_conditions_zt <- data.frame('Comparison'=comparison, 'chi_squared'= chisq, 'df'=df, "p value"=pval)
kruskal_wallis_insulin_conditions_zt

#############################################################
#glucose

#Serum Glucose Results 
df_glucose <- data.frame(read.csv("files_from_KF/STAM_TRF/ELISA_Results/stam_trf_glucose_measurements.csv", header=TRUE))%>%
  filter(cohort=="NASH") %>%
  mutate(Condition=ifelse(Condition=="NA_","NA",Condition))%>%
  mutate(Condition=factor(Condition,levels=c("NA","FA","FT")),
         zt=factor(zt,levels=c("ZT1","ZT13")))

#without control samples
glucose <- ggplot(df_glucose, aes(x=Condition, y=Real.Concentration.mg_dL, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  labs(title="Serum Glucose Across Concentration, 12 wk") +
  ylab("Serum Glucose Concentration (mg/dL)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggsave("phenotype/serum_glucose_conc.pdf", plot=glucose,height=4, width=4)

##split by ZT
glucose <- ggplot(df_glucose, aes(x=Condition, y=Real.Concentration.mg_dL, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  facet_grid(~zt)+
  labs(title="Serum Glucose Across Concentration, 12 wk") +
  ylab("Serum Glucose Concentration (mg/dL)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggsave("phenotype/serum_glucose_conc_byZT.pdf", plot=glucose,height=4, width=6)

#Mann Whitney ZT1 vs ZT13 per condition
W <- c()
pval <- c()
Conditions <- c("NA", "FA", "FT")

for (i in 1:length(Conditions)){
  df_loop <- df_glucose%>% filter(Condition==Conditions[i])
  res <- wilcox.test(data=df_loop, Real.Concentration.mg_dL ~ zt)
  W <- c(W, res$statistic)
  pval <- c(pval, res$p.value)
}
mann_whitney_glucose_zt1_V_zt13 <- data.frame('Condition'=Conditions, 'W stat'= W, "p value"=pval)
mann_whitney_glucose_zt1_V_zt13 #not ZT1 vs ZT13 diff for any condition

#Kruskal Wallis FA vs FT vs NA
chisq <- c()
df <- c()
pval <- c()
comparison <- c("FAvFT", "NAvFT","NAvFA")
Conditions <- c("NA", "FA", "FT")

for (i in 1:length(comparison)){
  df_loop <- subset(df_glucose, Condition!=Conditions[i])
  res <- kruskal.test(data=df_loop, Real.Concentration.mg_dL ~ Condition)
  chisq <- c(chisq, res$statistic)
  df <- c(df, res$parameter)
  pval <- c(pval, res$p.value)
}
kruskal_wallis_glucose_conditions <- data.frame('Comparison'=comparison, 'chi_squared'= chisq, 'df'=df, "p value"=pval)
kruskal_wallis_glucose_conditions #FA is sig diff from FA and FT

pairwise.wilcox.test(df_glucose$Real.Concentration.mg_dL,df_glucose$Condition,p.adjust.method = "fdr")

# NA     FA    
# FA 0.0003 -     
#   FT 0.0003 0.1782

#Kruskal Wallis FA vs FT vs NA, split by ZT
chisq <- c()
df <- c()
pval <- c()
comparison <- c("ZT1_FAvFT", "ZT1_NAvFT","ZT1_NAvFA","ZT13_FAvFT", "ZT13_NAvFT","ZT13_NAvFA")
Conditions <- c("NA", "FA", "FT","NA", "FA", "FT")
ZT <- c("ZT1","ZT1","ZT1","ZT13","ZT13","ZT13")

for (i in 1:length(comparison)){
  df_loop <- subset(df_glucose, zt==ZT[i] & Condition!=Conditions[i])
  res <- kruskal.test(data=df_loop, Real.Concentration.mg_dL ~ Condition)
  chisq <- c(chisq, res$statistic)
  df <- c(df, res$parameter)
  pval <- c(pval, res$p.value)
}
kruskal_wallis_glucose_conditions_zt <- data.frame('Comparison'=comparison, 'chi_squared'= chisq, 'df'=df, "p value"=pval)
kruskal_wallis_glucose_conditions_zt #regarless of ZT time NA is sig diff from FA and NA
#############################################################
#NASH score

df <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         NASH_category=ifelse(NASH_category=="N/A","not_applicable",NASH_category),
         steatosis_grade=as.character(steatosis_grade),
         lobular_inflammation=as.character(lobular_inflammation),
         hepatocyte_ballooning=as.character(hepatocyte_ballooning),
         portal_inflammation=as.character(portal_inflammation),
         nucleomegallyanisonucleosis=as.character(nucleomegallyanisonucleosis),
         fibrosis_stage=as.character(fibrosis_stage),
         NASH_score=as.character(NASH_score),
         glycogenated_nuclei=factor(glycogenated_nuclei, level=c('absent', 'rare', 'few', 'many')))%>%
  mutate(NASH_category=factor(NASH_category,level=c("Non_NASH","NASH","not_applicable")),
         condition=factor(condition,level=c("NA","FA","FT")))
  

# X = grade, Y = freq of each condition
relative_freq <- function(df_column){
  my_plot <- ggplot(data=df, aes(x=df_column, fill=condition))+
    geom_bar(position="fill")+
    #theme_pubr()+
    scale_fill_manual(values = friendly_pal("ito_seven")) +
    scale_x_discrete(expand = c(0, 0.6)) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("relative frequency")
  
  my_plot
}

p<-relative_freq(df$steatosis_grade) + xlab("steatosis_grade")
ggsave("phenotype/NASH_steatosis_grade_bycond.pdf", plot=p,height=4, width=4)

# X= condition, Y= % at each grade
myplot <- function(column){
  this_plot <- ggplot(data=df, aes(x=condition, fill=column))+
    geom_bar(position="fill", color="black")+
    ylab("Relative Frequency")+
    scale_x_discrete(expand = c(0, 0.6)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())#remove grid lines
  
  this_plot
}

#Graphing. Each separately so that title can be added to legend 
p<-myplot(df$NASH_score) + guides(fill=guide_legend(title="NAS Score")) + 
  scale_fill_manual(values=c("#648FFF","#785EF0","#DC267F","#FFB000","#FE6100"))
ggsave("phenotype/NASH_cond_byNASscore.pdf", plot=p,height=3, width=3)

p<-myplot(df$steatosis_grade) + guides(fill=guide_legend(title="Steatosis Grade"))+ 
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave("phenotype/NASH_cond_bysteatosisgrd.pdf", plot=p,height=3, width=4)

p<-myplot(df$lobular_inflammation) + guides(fill=guide_legend(title="Lobular Inflammation"))+
  scale_fill_manual(values=c("#648FFF","#DC267F"))
ggsave("phenotype/NASH_cond_bylobimfl.pdf", plot=p,height=3, width=4.5)

p<-myplot(df$hepatocyte_ballooning) + guides(fill=guide_legend(title="Hepatocyte Ballooning"))+ 
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave("phenotype/NASH_cond_byhepballon.pdf", plot=p,height=3, width=4.5)

p<-myplot(df$lipogranuloma) + guides(fill=guide_legend(title="Lipogranuloma"))+
  scale_fill_manual(values=c("#648FFF","#DC267F"))
ggsave("phenotype/NASH_cond_bylipgranul.pdf", plot=p,height=3, width=4)

p<-myplot(df$portal_inflammation) + guides(fill=guide_legend(title="Portal Inflammation"))+
  scale_fill_manual(values=c("#648FFF","#DC267F"))
ggsave("phenotype/NASH_cond_byportalinflm.pdf", plot=p,height=3, width=4.2)

p<-myplot(df$glycogenated_nuclei) + guides(fill=guide_legend(title="Glycogenated Nuclei"))+ 
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave("phenotype/NASH_cond_byglyconucl.pdf", plot=p,height=3, width=4.5)

p<-myplot(df$nucleomegallyanisonucleosis) + guides(fill=guide_legend(title="Nucleomegaly"))+ 
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave("phenotype/NASH_cond_bynuclmega.pdf", plot=p,height=3, width=4)

p<-myplot(df$fibrosis_stage) + guides(fill=guide_legend(title="Fibrosis Stage"))+ 
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave("phenotype/NASH_cond_fibrosisstage.pdf", plot=p,height=3, width=3)


#statistics, is there a diff among all three conditions as they relate to NASH scoring
fisher.test(table(df$NASH_score, df$condition)) #p-value = 0.004539
fisher.test(table(df$steatosis_grade, df$condition)) #p-value = 0.00005287
fisher.test(table(df$lobular_inflammation, df$condition)) #p-value = 1
fisher.test(table(df$hepatocyte_ballooning, df$condition)) #p-value = 0.881
fisher.test(table(df$lipogranuloma, df$condition)) #p-value = 0.3319
#fisher.test(table(df$portal_inflammation, df$condition)) #doesnt work, only one grade
fisher.test(table(df$glycogenated_nuclei, df$condition)) #p-value = 0.08653
fisher.test(table(df$nucleomegallyanisonucleosis, df$condition)) #p-value = 0.1073
fisher.test(table(df$fibrosis_stage, df$condition)) #p-value = 0.003068

df_sub<-df%>%filter(condition!="NA")
#statistics, is there a diff among FA vs. FT as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.5717
fisher.test(table(df_sub$steatosis_grade, df_sub$condition)) #p-value = 1
fisher.test(table(df_sub$lobular_inflammation, df_sub$condition)) #p-value = 1
fisher.test(table(df_sub$hepatocyte_ballooning, df_sub$condition)) #p-value = 0.6802
fisher.test(table(df_sub$lipogranuloma, df_sub$condition)) #p-value = 1
#fisher.test(table(df_sub$portal_inflammation, df_sub$condition)) #doesnt work, only one grade
fisher.test(table(df_sub$glycogenated_nuclei, df_sub$condition)) #p-value = 0.2138
fisher.test(table(df_sub$nucleomegallyanisonucleosis, df_sub$condition)) #p-value = 0.06865
fisher.test(table(df_sub$fibrosis_stage, df_sub$condition)) #p-value = 0.8278

df_sub<-df%>%filter(condition!="FT")
#statistics, is there a diff among NA vs. Fa as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.004748
#fisher.test(table(df_sub$steatosis_grade, df_sub$condition)) #doesnt work
#fisher.test(table(df_sub$lobular_inflammation, df_sub$condition)) #doesnt work
fisher.test(table(df_sub$hepatocyte_ballooning, df_sub$condition)) #p-value = 1
fisher.test(table(df_sub$lipogranuloma, df_sub$condition)) #p-value = 0.4783
#fisher.test(table(df_sub$portal_inflammation, df_sub$condition)) #doesnt work, only one grade
fisher.test(table(df_sub$glycogenated_nuclei, df_sub$condition)) #p-value = 0.6668
fisher.test(table(df_sub$nucleomegallyanisonucleosis, df_sub$condition)) #p-value = 0.4003
fisher.test(table(df_sub$fibrosis_stage, df_sub$condition)) #p-value = 0.01373

df_sub<-df%>%filter(condition!="FA")
#statistics, is there a diff among NA vs. FT as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.007352
fisher.test(table(df_sub$steatosis_grade, df_sub$condition)) #p-value = 0.0000673
fisher.test(table(df_sub$lobular_inflammation, df_sub$condition)) #p-value = 1
fisher.test(table(df_sub$hepatocyte_ballooning, df_sub$condition)) #p-value = 0.6802
fisher.test(table(df_sub$lipogranuloma, df_sub$condition)) #p-value = 0.2174
#fisher.test(table(df_sub$portal_inflammation, df_sub$condition)) #doesnt work, only one grade
fisher.test(table(df_sub$glycogenated_nuclei, df_sub$condition)) #p-value = 0.03913
fisher.test(table(df_sub$nucleomegallyanisonucleosis, df_sub$condition)) #p-value = 0.5901
fisher.test(table(df_sub$fibrosis_stage, df_sub$condition)) #p-value = 0.001346

#############################################################
#metabolic cages

md <- read.csv("files_from_KF/STAM_TRF/Metabolic_Cage_Analysis/stam_trf_3_3_21/data/stam_trf_3_3_21_metadata.csv", 
                 colClasses = c(Animal="factor"),strip.white = TRUE)%>%
  mutate(Group=ifelse(is.na(Group),"NA",Group))

# df_old <- read.csv("files_from_KF/STAM_TRF/Metabolic_Cage_Analysis/stam_trf_3_3_21/data/stam_trf_3_3_21_1hr_macro_timeseries.csv", na.strings = '.') %>%
# filter(!is.na(Animal))%>%
#   mutate(DateTime=as.POSIXct(DateTime, tz = "GMT"))
df <- read.csv("phenotype/stam_trf_run3_3_3_21_fogelson_macro_OneClickMacroV.2.53.2-slice1hr.mac_1_slide1.csv", na.strings = '.') %>%
  filter(!is.na(Animal))%>%
  mutate(DateTime=as.POSIXct(DateTime, tz = "GMT"))

#Add Water and Food Consumption Column (per 1 hr bin), and cumulative food consumption 
for (i in 1:nrow(df)){
  if (i-1 == 0){
    df$WaterCons[i] <- 0
    df$FoodCons[i] <- 0
    df$food_cumulative[i] <- 0
  } else if (df$Animal[i] != df$Animal[i-1]){
    df$WaterCons[i] <- 0
    df$FoodCons[i] <- 0
    df$food_cumulative[i] <- 0
  }else {
    water_diff <- df$WaterInA_M[i] - df$WaterInA_M[i-1]
    food_diff <- df$FoodInA_M[i] - df$FoodInA_M[i-1]
    food_cumulative <- df$FoodInA_M[i] + df$FoodInA_M[i-1]
    df$WaterCons[i] <- water_diff
    df$FoodCons[i] <- food_diff
    df$food_cumulative[i] <- food_cumulative
  }
}

#Isolate the first 48 hours (need to drop data from acclimation period)
acc_period <- df$DateTime[c(1:48)]
#Drop data from acclimation period (dropping 672 rows = 14 mice * 48 timepoints)
df_cleaned <- df[!(df$DateTime %in% acc_period),]
#Isolate time after 3 NA mice died (need to drop)
df_cleaned_only_NA <- subset(df_cleaned, df_cleaned$Animal %in% c("3",'11','15'))
#na_died <- unique(df_cleaned$DateTime[(df_cleaned$DateTime > as.POSIXct('2021-03-09 00:08:00', tz="GMT"))]) #NA 9bR, Animal 3, NA 11R, Animal 11, NA 12L, Animal 15
na_died <- unique(df_cleaned$DateTime[(df_cleaned$DateTime > as.POSIXct('2021-03-09 00:37:33', tz="GMT"))])
df_cleaned_na_died <- df_cleaned_only_NA[(df_cleaned_only_NA$DateTime %in% na_died),]
#Drop data after 9bR, 11R, and 12L died (only for NA group)
df_dropped <- anti_join(df_cleaned,df_cleaned_na_died, by = c("DateTime","Animal")) #dropped 312 rows (new df its 390)
#Drop ALL data after 9bR, 11R, and 12L died (for ALL groups)
df_cleanest <- df_cleaned[!(df_cleaned$DateTime %in% na_died),]

#zt_time_values <- data.frame(DateTime_hour = c(20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19), zt_hour = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),stringsAsFactors=FALSE)
zt_time_values <- data.frame(DateTime_hour = c(21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), zt_hour = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),stringsAsFactors=FALSE)


df_cleanest$MotorActivity <- df_cleanest$AllMeters_R - df_cleanest$PedMeters_R #adds column title "MotorActivity" and populates with values specified on right side of the equation
df1_cleanest <- mutate(df_cleanest, DateTime=hour(DateTime)) #create a new df with DateTime in proper formation. Mutate function adds new variable to df while preserving the existing variable. Hour function takes data from POSIXct object, and assigns numeric value to the hour component. 
df2_cleanest <- group_by(df1_cleanest, DateTime, Animal) #create new df where DateTime and Animal are grouped
#df3 <- summarize(df2, vo2m=mean(VO2_M)) 
df3_cleanest <- summarise_all(df2_cleanest, mean, na.rm = TRUE) #create new df and average all data from a single DateTime across all days of the metabolic cage run. There are about 7 days of data for each animal (i.e. 7 replicate DateTimes). summarise_all function averages data across all days  for each hour slice for each animal
df3_cleanest <- as.data.frame(df3_cleanest) #drops all non-summarized columns

for (i in 1:nrow(df3_cleanest)) {
  for (z in 1:nrow(zt_time_values)) {
    if (df3_cleanest$DateTime[i] == zt_time_values$DateTime_hour[z]) {
      df3_cleanest$zt[i] <- zt_time_values$zt_hour[z]
    }
  }
}

df_group_cleanest<- merge(x=df3_cleanest, y=md, by="Animal") #merge df and md by Animal
df_group_cleanest$Animal <- as.factor(df_group_cleanest$Animal) #factor animal column (prev. numerical)

#reorder dataframe by zt and Animal
df_group_cleanest <- 
  df_group_cleanest[order(df_group_cleanest$Animal,
                          df_group_cleanest$zt), ]

#Add Kcal to FoodCons (per pin)
for (i in 1:nrow(df_group_cleanest)) {
  if (df_group_cleanest$Group[i]== 'FT'){
    df_group_cleanest$kcalCons[i] <- df_group_cleanest$FoodCons[i] * 5.21
  }
  else if (df_group_cleanest$Group[i]== 'FA'){
    df_group_cleanest$kcalCons[i] <- df_group_cleanest$FoodCons[i] * 5.21
  }
  else if (df_group_cleanest$Group[i]== 'NA'){
    df_group_cleanest$kcalCons[i] <- df_group_cleanest$FoodCons[i] * 3.1
  }
}


#Avg cumulative kcal intake 
for (i in 1:nrow(df_group_cleanest)){
  if (i-1 == 0){
    df_group_cleanest$cumulative_Kcal[i] <-
      df_group_cleanest$kcalCons[i]
  } else if (df_group_cleanest$Animal[i] != df_group_cleanest$Animal[i-1]){
    df_group_cleanest$cumulative_Kcal[i] <-
      df_group_cleanest$kcalCons[i]
  }else {
    cumulative_Kcal<- df_group_cleanest$cumulative_Kcal[i-1] +
      df_group_cleanest$kcalCons[i]
    df_group_cleanest$cumulative_Kcal[i] <- cumulative_Kcal
  }
}

#Remove FT condition, comment out if you want to plot all 3 
#df_group_cleanest <- subset(df_group_cleanest, Group != "FT")

df_group_cleanest<-df_group_cleanest %>%
group_by(Animal, Group) %>%
mutate(tot_kcal_cons= sum(kcalCons))%>%
mutate(percent_cumulative_Kcal=cumulative_Kcal/tot_kcal_cons)

#Avg cumulative food intake 
# avg_cumfood_summary<-df_group_cleanest %>%
#   group_by(Animal, Group) %>%
#   summarise(tot_kcal_cons= sum(kcalCons))

#Percent daily kcal intake 
# for (i in 1:nrow(df_group_cleanest)) {
#   if (df_group_cleanest$Animal[i]== 1){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/11.27
#   }else if (df_group_cleanest$Animal[i]== 2){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/15.37056
#   }else if (df_group_cleanest$Animal[i]== 3){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/18.36547
#   }else if (df_group_cleanest$Animal[i]== 4){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/6.5
#   }else if (df_group_cleanest$Animal[i]== 5){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/10.76
#   }else if (df_group_cleanest$Animal[i]== 7){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/12.84069
#   }else if (df_group_cleanest$Animal[i]== 8){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/11.35
#   }else if (df_group_cleanest$Animal[i]== 9){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/12.39
#   }else if (df_group_cleanest$Animal[i]== 11){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/24.47451
#   }else if (df_group_cleanest$Animal[i]== 12){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/11.04
#   } else if (df_group_cleanest$Animal[i]== 13){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/9.13
#   } else if (df_group_cleanest$Animal[i]== 14){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/26.13016
#   } else if (df_group_cleanest$Animal[i]== 15){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/22.54457
#   }else if (df_group_cleanest$Animal[i]== 16){
#     df_group_cleanest$percent_cumulative_Kcal[i] <-
#       df_group_cleanest$cumulative_Kcal[i]/17
#   }
# }

df_group_cleanest<-df_group_cleanest%>%mutate(Group=factor(Group,levels = c("NA","FA","FT")))%>%
  dplyr::rename(condition=Group)

RER_M_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=RER_M, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0.5,1)) +
  labs(title="Mean respiratory exchange ratio (RER)") +
  theme_classic()+#theme_pubr()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_mRER_v2.pdf", plot=RER_M_SEM_cleanest,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_mRER_L4h.pdf", plot=RER_M_SEM_cleanest,height=3.5, width=3.5)


MotorActivity_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=MotorActivity, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title='Motor Activity (m)') + theme_classic()+#theme_pubr()+
  theme(legend.position="top", plot.title = element_text(size = 12), axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_motoractivity_v2.pdf", plot=MotorActivity_SEM_cleanest,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_motoractivity_L4h.pdf", plot=MotorActivity_SEM_cleanest,height=3.5, width=3.5)

WaterCons_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=WaterCons, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme_classic()+#theme_pubr()+
  labs(title='Mass of water consumed (g)') + theme(legend.position="top",plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_masswaterconsum_v2.pdf", plot=WaterCons_SEM_cleanest,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_masswaterconsum_L4h.pdf", plot=WaterCons_SEM_cleanest,height=3.5, width=3.5)


#Mean Energy Expenditure 
kcal_hr_M_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=kcal_hr_M, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title="Mean energy expenditure") +
  ylab("Kcal/hr")+
  xlab("ZT") + 
  theme_classic()+#theme_pubr()+
  theme(legend.position="top",plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_energyexp_v2.pdf", plot=kcal_hr_M_SEM_cleanest,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_energyexp_L4h.pdf", plot=kcal_hr_M_SEM_cleanest,height=3.5, width=3.5)

#Food consumed, grams 
FoodCons_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=FoodCons, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') + 
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,1.5)) + theme_classic()+ #theme_pubr()+
  labs(title='Mass of food consumed (g)') + theme(legend.position="top",plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_foodconsumed_g_v2.pdf", plot=FoodCons_SEM_cleanest,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_foodconsumed_g_L4h.pdf", plot=FoodCons_SEM_cleanest,height=3.5, width=3.5)


#Kcal consumed
KcalCons_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=kcalCons, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,4)) +
  labs(title="Food Consumption")+
  ylab("Avg Kcal Consumed") +
  xlab("ZT")+
  theme_classic()+#theme_pubr()+
  theme(legend.position="top",plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_foodconsumed_kcal_v2.pdf", plot=KcalCons_SEM_cleanest ,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_foodconsumed_kcal_L4h.pdf", plot=KcalCons_SEM_cleanest ,height=3.5, width=3.5)

#Percent cumulative Kcal
percent_cumulative_KcalCons_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=percent_cumulative_Kcal, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,1)) +
  labs(title="% of Daily Food Consumption")+
  ylab("% Kcal")+
  xlab("ZT")+
  theme_classic()+#theme_pubr()+
  theme(legend.position="top",plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

#ggsave("phenotype/metab_cage_cumulkcalconsumed_v2.pdf", plot=percent_cumulative_KcalCons_SEM_cleanest ,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_cumulkcalconsumed_L4h.pdf", plot=percent_cumulative_KcalCons_SEM_cleanest ,height=3.5, width=3.5)


#sleep
sleep_SEM_cleanest <- ggplot(data=df_group_cleanest, aes(x=zt, y=Sleep_pct_M, colour=condition, fill=condition)) +
  #geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf),
  geom_rect(aes(xmin = 4, xmax = 23, ymin = -Inf, ymax = Inf),
            color="transparent", fill = 'gray88') +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title="Sleep")+
  ylab("Sleep_pct_M")+
  xlab("ZT")+
  theme_classic()+#theme_pubr()+
  theme(legend.position="top",plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))

ggsave("phenotype/metab_cage_sleep_v2.pdf", plot=sleep_SEM_cleanest ,height=3.5, width=3.5)
ggsave("phenotype/metab_cage_sleep_L4h.pdf", plot=sleep_SEM_cleanest ,height=3.5, width=3.5)




# SEM_figure_cleanest <- ggarrange(RER_M_SEM_cleanest,MotorActivity_SEM_cleanest,kcal_hr_M_SEM_cleanest, ncol = 3, nrow = 1)
# ggsave("phenotype/metab_cage_energyplot.pdf", plot=SEM_figure_cleanest ,height=3.5, width=10)
# 
# SEM_food_water_cleanest <- ggarrange(WaterCons_SEM_cleanest, FoodCons_SEM_cleanest, KcalCons_SEM_cleanest, percent_cumulative_KcalCons_SEM_cleanest, ncol = 4, nrow = 1)
# ggsave("phenotype/metab_cage_foodwaterplot.pdf", plot=SEM_food_water_cleanest ,height=3.5, width=12)
