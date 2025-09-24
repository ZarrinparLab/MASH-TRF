setwd("~/Notebooks/sfloresr/MASH-TRF/JBR_2025")

library(tidyverse)
library(data.table)
library(ggpubfigs)
library(ggpubr)
library(agricolae)

#############################################################
#inputs
food_intake<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/STAM_TRF_mouse_food_weights_NASH_separated_mice_grouped.csv"
weights<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/STAM_TRF_mouse_weights_NASH.csv"
glucoselev<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/stam_trf_glucose_measurements.csv"
insulinlev<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/stam_trf_insulin_measurements.csv"
histology_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/nas_scoring_analysis_nash.csv"
tg_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/mash_liver_tg.txt"
results<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/results/"

#############################################################
#functions

myplot <- function(column){ # X= condition, Y= % at each grade
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

#############################################################
#NASH cohort cumulative food intake--Fig. 1b

food_NASH <- read.csv(food_intake,  header=TRUE)
food_NASH$condition[is.na(food_NASH$condition)] <- "NA" #rename NA group

food_NASH <- group_by(food_NASH, week, condition) #group by week & condition 

food_NASH_mean<- summarise(food_NASH,avg_cumulative_kcal_per_mouse= mean(avg_cumulative_kcal_per_mouse)) #average weight per group per week
food_NASH_SEM<- summarise(food_NASH, SEM= sd(avg_cumulative_kcal_per_mouse)/sqrt(length(avg_cumulative_kcal_per_mouse))) #SEM per group per week
food_NASH_summarized <- merge(food_NASH_mean, food_NASH_SEM) %>% # merge mean weight & SEM tables
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

NASH_food <- ggplot(food_NASH_summarized, aes(x= week, y= avg_cumulative_kcal_per_mouse))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=avg_cumulative_kcal_per_mouse-SEM, ymax=avg_cumulative_kcal_per_mouse+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition), colour="black",pch=21, size=2) +
  labs(title="Average Cumulative Food\nIntake vs Time") +
  ylab("Food Intake (kcal/mouse)")+
  xlab("Week (post intervention)")+
  theme_classic()+
  scale_x_continuous(breaks = seq(1, 7, by = 1))+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top")

ggsave(paste0(results,"cumul_food_intake.pdf"), plot=NASH_food,height=3.2, width=3)

#ANOVA + LSD test
food_NASH<-food_NASH%>%mutate(week=as.factor(week))
foodint_m <- aov(avg_cumulative_kcal_per_mouse ~condition * week, data = food_NASH)
summary(foodint_m)

LSD_Test<- (LSD.test(foodint_m, c("condition", "week")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","week"))
write.table(LSD_Test,paste0(results,"SFR24_1027_foodintvweek_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#############################################################
#NASH cohort weights--Fig. 1c

weight_NASH <- read.csv(weights,  header=TRUE)
weight_NASH$condition[is.na(weight_NASH$condition)] <- "NA" #rename NA group

weight_NASH <- group_by(weight_NASH, week, condition) #group by week & condition 

weight_NASH_mean<- summarise(weight_NASH,weight_g= mean(weight_g)) #average weight per group per week
weight_NASH_SEM<- summarise(weight_NASH, SEM= sd(weight_g)/sqrt(length(weight_g))) #SEM per group per week
weight_NASH_summarized <- merge(weight_NASH_mean, weight_NASH_SEM) %>% # merge mean weight & SEM tables
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

NASH_weights_plot <- ggplot(weight_NASH_summarized, aes(x= week, y= weight_g))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=weight_g-SEM, ymax=weight_g+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition), colour="black",pch=21, size=2) +
  labs(title="Average Mouse\nWeight vs Time") +
  ylab("Weight (g)")+
  xlab("Week (after intervention)")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 6, by = 1))+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") 

ggsave(paste0(results,"avg_mouse_weight.pdf"), plot=NASH_weights_plot,height=3.2, width=3)

#ANOVA + LSD test
weight_NASH<-weight_NASH%>%mutate(week=as.factor(week))
weight_m <- aov(weight_g ~condition * week, data = weight_NASH)
summary(weight_m)

LSD_Test<- (LSD.test(weight_m, c("condition", "week")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","week"))
write.table(LSD_Test,paste0(results,"SFR24_1027_weightvweek_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#############################################################
#Serum Glucose Results--Fig. 1d

df_glucose <- data.frame(read.csv(glucoselev, header=TRUE))%>%
  filter(cohort=="NASH") %>%
  mutate(Condition=ifelse(Condition=="NA_","NA",Condition))%>%
  mutate(Condition=factor(Condition,levels=c("NA","FA","FT")),
         zt=factor(zt,levels=c("ZT1","ZT13")))

glucose <- ggplot(df_glucose, aes(x=Condition, y=Real.Concentration.mg_dL, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  labs(title="Serum Glucose Across Concentration, Wk 7") +
  ylab("Serum Glucose Concentration (mg/dL)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggsave(paste0(results,"serum_glucose_conc.pdf"), plot=glucose,height=4, width=4)

#ANOVA + LSD test
glucose_m <- aov(Real.Concentration.mg_dL ~Condition, data = df_glucose)
summary(glucose_m)

LSD_Test<- (LSD.test(glucose_m, c("Condition")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Condition"))
write.table(LSD_Test,paste0(results,"SFR24_1027_glucose_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#############################################################
#Serum Glucose Results by ZT--Fig. S1a

glucose <- ggplot(df_glucose, aes(x=Condition, y=Real.Concentration.mg_dL, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  facet_grid(~zt)+
  labs(title="Serum Glucose Across Concentration, Wk 7") +
  ylab("Serum Glucose Concentration (mg/dL)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggsave(paste0(results,"serum_glucose_conc_byZT.pdf"), plot=glucose,height=4, width=6)

#ANOVA + LSD test
glucose_m <- aov(Real.Concentration.mg_dL ~Condition * zt, data = df_glucose)
summary(glucose_m)

LSD_Test<- (LSD.test(glucose_m, c("Condition", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Condition","zt"))
write.table(LSD_Test,paste0(results, "SFR24_1027_glucosevzt_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#############################################################
#Serum Insulin Results -- Fig. 1e

df_insulin<- data.frame(read.csv(insulinlev, header=TRUE))%>%
  filter(cohort=="NASH")%>%
  mutate(Condition=ifelse(Condition=="NA_","NA",Condition))%>%
  mutate(Condition=factor(Condition,levels=c("NA","FA","FT")),
         zt=factor(zt,levels=c("ZT1","ZT13")))

insulin <- ggplot(df_insulin, aes(x=Condition, y=Concentration_ng_ml, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  labs(title="Serum Insulin Concentration, Wk 7") +
  ylab("Serum Insulin Concentration (ng/ml)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 

ggsave(paste0(results,"serum_insulin_conc.pdf"), plot=insulin,height=4, width=4)

#ANOVA + LSD test
insulin_m <- aov(Concentration_ng_ml ~Condition, data = df_insulin)
summary(insulin_m)

LSD_Test<- (LSD.test(insulin_m, c("Condition")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Condition"))
write.table(LSD_Test,paste0(results,"SFR24_1027_insulin_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#############################################################
#Serum Insulin Results by ZT -- Fig. s1b

insulin <- ggplot(df_insulin, aes(x=Condition, y=Concentration_ng_ml, group=Condition, label=Condition,fill=Condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  facet_grid(~zt)+
  labs(title="Serum Insulin Concentration, Wk 7") +
  ylab("Serum Insulin Concentration (ng/ml)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 

ggsave(paste0(results,"serum_insulin_conc_byZT.pdf"), plot=insulin,height=4, width=6)

#ANOVA + LSD test
insulin_m <- aov(Concentration_ng_ml ~Condition * zt, data = df_insulin)
summary(insulin_m)

LSD_Test<- (LSD.test(insulin_m, c("Condition", "zt")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("Condition","zt"))
write.table(LSD_Test,paste0(results,"SFR24_1027_insulinvzt_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#############################################################
#NASH scores by condition-Fig. 2a-c

df <- fread(histology_df, header=TRUE)%>%
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
  
p<-myplot(df$NASH_score) + guides(fill=guide_legend(title="NAS Score")) + #Fig. 2a
  scale_fill_manual(values=c("#648FFF","#785EF0","#DC267F","#FFB000","#FE6100"))
ggsave(paste0(results,"NASH_cond_byNASscore.pdf"), plot=p,height=2.5, width=3.5)

p<-myplot(df$fibrosis_stage) + guides(fill=guide_legend(title="Fibrosis Stage"))+ #Fig. 2b
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave(paste0(results,"NASH_cond_fibrosisstage.pdf"), plot=p,height=3, width=3)

p<-myplot(df$steatosis_grade) + guides(fill=guide_legend(title="Steatosis Grade"))+ #Fig. 2c
  scale_fill_manual(values=c("#648FFF","#FFB000","#DC267F"))
ggsave(paste0(results,"NASH_cond_bysteatosisgrd.pdf"), plot=p,height=3, width=3)

df_sub<-df%>%filter(condition!="NA")
#statistics, is there a diff among FA vs. FT as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.5717
fisher.test(table(df_sub$steatosis_grade, df_sub$condition)) #p-value = 1
fisher.test(table(df_sub$fibrosis_stage, df_sub$condition)) #p-value = 0.8278

df_sub<-df%>%filter(condition!="FT")
#statistics, is there a diff among NA vs. Fa as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.004748
fisher.test(table(df_sub$steatosis_grade, df_sub$condition)) #doesnt work
fisher.test(table(df_sub$fibrosis_stage, df_sub$condition)) #p-value = 0.01373

df_sub<-df%>%filter(condition!="FA")
#statistics, is there a diff among NA vs. FT as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.007352
fisher.test(table(df_sub$steatosis_grade, df_sub$condition)) #p-value = 0.0000673
fisher.test(table(df_sub$fibrosis_stage, df_sub$condition)) #p-value = 0.001346

#############################################################
#Liver Triglycerides Results -- Fig. 2d

df_trigly <- fread(tg_df,na.strings = "")%>%
  dplyr::rename(conc_mgg=`Liver TG concentration (mg/g)`)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         ZT=factor(ZT,levels=c("ZT1","ZT13")))

trigly <- ggplot(df_trigly, aes(x=condition, y=conc_mgg, group=condition, label=condition,fill=condition ))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  labs(title="Liver Triglyceride Concentration (Wk 7)") +
  ylab("Liver Triglyceride Concentration (mg/g liver)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") 
ggsave(paste0(results,"liver_triglyceride_conc.pdf"), plot=trigly,height=3.5, width=3.5)

#ANOVA + LSD test
trigly_m <- aov(conc_mgg ~condition, data = df_trigly)
summary(trigly_m)

LSD_Test<- (LSD.test(trigly_m, c("condition")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition"))
write.table(LSD_Test,paste0(results,"SFR24_1027_triglycerides_pvals_anovaLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
