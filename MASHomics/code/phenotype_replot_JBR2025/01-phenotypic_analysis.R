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
metag_cage<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/first48h_metabcagetable.txt"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/phenotype_replot_JBR2025/"
#############################################################
#functions

overZTplt<-function(df,yval){
  p<-ggplot(data=df, aes(x=zt_hour, y=yval, colour=condition, fill=condition)) +
    geom_rect(aes(xmin = 13, xmax = 21, ymin = -Inf, ymax = Inf),
              color="transparent", fill = 'gray88') +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, aes(colour=condition)) +
    scale_fill_manual(values = friendly_pal("ito_seven")) +
    scale_color_manual(values = friendly_pal("ito_seven")) +
    scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
    theme_classic()+
    theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  return(p)
}
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

#############################################################
#metabolic cage results by condition-Fig. S3A (replotted from JBR 2025)

f48h<-fread(metag_cage)

#Food Consumption (Kcal, left)
KcalCons_SEM_cleanest <- overZTplt(f48h,f48h$kcalCons)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,4)) +
  labs(title="Food Consumption")+
  ylab("Avg Kcal Consumed") +
  xlab("ZT")
ggsave(paste0(results_dir,"SFR24_0216_metab_cage_foodconsumed_kcal_f48h_new.pdf"), plot=KcalCons_SEM_cleanest ,height=3.5, width=3.5)

kcalCons_mod <- aov(kcalCons~condition * zt_hour, data = f48h)
summary(kcalCons_mod)

LSD_Test_kcalCons<- (LSD.test(kcalCons_mod, c("condition", "zt_hour")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_hour"))
write.table(LSD_Test_kcalCons,paste0(results_dir,"SFR24_1027_metab_cage_foodconsumed_kcal_bytmpt_48h_anovaFisherLSD_new.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#Motor Activity (middle)
MotorActivity_SEM_cleanest <- overZTplt(f48h,f48h$MotorActivity)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title='Motor Activity (m)',y="MotorActivity")
ggsave(paste0(results_dir,"SFR24_0216_metab_cage_motoractivity_48h_new.pdf"), plot=MotorActivity_SEM_cleanest,height=3.5, width=3.5)

mamod <- aov(MotorActivity ~condition * zt_hour, data = f48h)
summary(mamod)

LSD_Test_ma<- (LSD.test(mamod, c("condition", "zt_hour")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_hour"))
write.table(LSD_Test_ma,paste0(results_dir,"SFR24_1027_metab_cage_motoract_bytmpt_48h_anovaFisherLSD_new.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

# % Sleep (right)
sleep_SEM_cleanest <-overZTplt(f48h,f48h$Sleep_pct_M)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title="Sleep")+
  ylab("Sleep_pct_M")+
  xlab("ZT")
ggsave(paste0(results_dir,"SFR24_0216_metab_cage_sleep_48h_new.pdf"), plot=sleep_SEM_cleanest ,height=3.5, width=3.5)

sleepmod <- aov(Sleep_pct_M~ condition * zt_hour, data = f48h)
summary(sleepmod)

LSD_Test_sleep<- (LSD.test(sleepmod, c("condition", "zt_hour")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_hour"))
write.table(LSD_Test_sleep,paste0(results_dir,"SFR24_1027_metab_cage_sleep_bytmpt_48h_anovaFisherLSD_new.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

