setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubfigs)
library(ggpubr)
library(ggbeeswarm)
library(agricolae)

#############################################################
#inputs
NASH_histology<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/nas_scoring_analysis_nash.csv"
liver_tgs<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/mash_liver_tg.txt"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/phenotype_replot_JBR2025/"
#############################################################
#NASH scores by condition-Fig. S1B (replotted from JBR 2025)

df <- fread(NASH_histology, header=TRUE)%>%
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
  

relative_freq <- function(df_column){
  my_plot <- ggplot(data=df, aes(fill=df_column, x=fct_rev(condition)))+
    geom_bar(position = "fill")+
    theme_pubr()+
    scale_fill_brewer(palette = "RdYlBu", direction = -1)+
    scale_x_discrete(expand = c(0, 0.6)) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("relative frequency")+coord_flip()
  
  my_plot
}

p<-relative_freq(df$NASH_score) + xlab("condition")
ggsave(paste0(results_dir,"NASH_cond_byNASscore_v2.pdf"), plot=p,height=3, width=4.5)

df_sub<-df%>%filter(condition!="NA")
#statistics, is there a diff among FA vs. FT as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.5717

df_sub<-df%>%filter(condition!="FT")
#statistics, is there a diff among NA vs. Fa as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.004748

df_sub<-df%>%filter(condition!="FA")
#statistics, is there a diff among NA vs. FT as they relate to NASH scoring
fisher.test(table(df_sub$NASH_score, df_sub$condition)) #p-value = 0.007352

#############################################################
#Liver Triglycerides Results -- Fig. S1C (replotted from JBR 2025)

df_trigly <- fread(liver_tgs,na.strings = "")%>%
  dplyr::rename(conc_mgg=`Liver TG concentration (mg/g)`)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         ZT=factor(ZT,levels=c("ZT1","ZT13")))

ggplot(df_trigly, aes(x = condition, y = conc_mgg, colour = condition)) +
  geom_beeswarm(size = 2, cex = 3) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.5, fatten = 2, color = "black", 
               aes(ymin = ..y.., ymax = ..y..)) +
  theme_bw(base_size = 14) +
  labs(x = "Condition", y = "Liver TG (mg/g)") +
  theme_pubr() +
  scale_colour_manual(values = friendly_pal("ito_seven")) +
  labs(title = "Liver Triglyceride Concentration (Wk 7)") +
  ylab("Liver Triglyceride Concentration (mg/g liver)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

ggsave(paste0(results_dir,"liver_triglyceride_conc_v2.pdf"), height=3.5, width=3.5)

#ANOVA + LSD test
trigly_m <- aov(conc_mgg ~condition, data = df_trigly)
summary(trigly_m)

LSD_Test<- (LSD.test(trigly_m, c("condition")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition"))