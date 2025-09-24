setwd("~/Notebooks/sfloresr/MASH-TRF/JBR_2025")

library(tidyverse)
library(data.table)
library(ggpubfigs)
library(ggpubr)
library(car)
library(agricolae)

###########################################################
#inputs
metadata_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/stam_trf_3_3_21_metadata.csv"
metabdata_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/stam_trf_run3_3_3_21_macro_OneClickMacroV.2.53.2-slice1hr.mac_1_fromsable.csv"
metabdata48h_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/first48h_metabcagetable.txt"
metabdata48hmean_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/first48h_metabcagetable_mean.txt"
results<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/results/"
###########################################################
#functions

overZTplt<-function(df,yval){
  p<-ggplot(data=df, aes(x=zt_hour, y=yval, colour=condition, fill=condition)) +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
    scale_fill_manual(values = friendly_pal("ito_seven")) +
    scale_color_manual(values = friendly_pal("ito_seven")) +
    scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
    theme_classic()+
    theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  return(p)
}

###########################################################
#read and process metabolic cage data

md <- read.csv(metadata_df, colClasses = c(Animal="factor"),strip.white = TRUE)%>%
  mutate(Group=ifelse(is.na(Group),"NA",Group))%>%
  arrange(Group)

df<- read.csv(metabdata_df, na.strings = '.',
               colClasses = c(Animal="factor")) %>%
  filter(!is.na(Animal))%>%
  mutate(DateTime=as.POSIXct(DateTime, tz = "GMT"))

zt_time_values <- data.frame(DateTime_hour = c(21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), zt_hour = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),stringsAsFactors=FALSE)

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

df_wmd<-df%>%left_join(.,md,by="Animal")%>%
  filter(!(DateTime > as.POSIXct('2021-03-09 00:37:33', tz="GMT") & (Animal %in% c("3",'11','15'))))%>%
  mutate(MotorActivity=AllMeters_R - PedMeters_R,
         DateTime_hour=hour(DateTime),
         condition=factor(Group,levels=c("NA","FA","FT")))%>%
  left_join(.,zt_time_values,by="DateTime_hour")%>%
  mutate(phase=ifelse(zt_hour>-1 & zt_hour<5, "light","dark"))

#keep just the first 48 hours
f48h <- df_wmd %>%
  arrange(Animal, DateTime) %>%
  group_by(Animal) %>%
  slice_head(n = 48) %>%
  ungroup()%>%
  mutate(kcalCons=ifelse(Group=="NA",FoodCons*3.1,FoodCons*5.21))%>%
  group_by(Animal,Group) %>%
  arrange(zt_hour)%>%
  mutate(cumulative_Kcal=cumsum(kcalCons))%>%
  mutate(tot_kcal_cons= sum(kcalCons))%>%
  mutate(percent_cumulative_Kcal=cumulative_Kcal/tot_kcal_cons)%>%
  mutate(zt_time_binned = case_when(
    between(zt_hour, 0, 5) ~ "5",
    between(zt_hour, 4, 9) ~ "9",
    between(zt_hour, 8, 13) ~ "13",
    between(zt_hour, 12, 17) ~ "17",
    between(zt_hour, 16, 21) ~ "21",
    TRUE ~ "1" ),
    phase=ifelse(zt_hour<13, "light","dark"))%>%
  mutate(zt_time_binned=factor(zt_time_binned,levels=c("1","5","9","13","17","21")),
         phase=factor(phase,levels=c("light","dark")))%>%
  filter(!(Animal==11 & zt_hour<17 & zt_hour>14))

write.table(f48h,"data/first48h_metabcagetable.txt",sep = "\t",row.names = FALSE,quote=FALSE)

f48h_ZT<-f48h%>%
  group_by(Animal,zt_time_binned,condition)%>%
  summarise(mn_RER=mean(RER_M),sd_RER = sd(RER_M),
            mn_MotorAct=mean(MotorActivity),sd_MotorAct=sd(MotorActivity),
            mn_WaterCons=mean(WaterCons),sd_WaterCons=sd(WaterCons),
            mn_kcal_hr_M=mean(kcal_hr_M),sd_kcal_hr_M=sd(kcal_hr_M),
            mn_FoodCons=mean(FoodCons),sd_FoodCons=sd(FoodCons),
            mn_kcalCons=mean(kcalCons),sd_kcalCons=sd(kcalCons),
            mn_Sleep_pct_M=mean(Sleep_pct_M),sd_Sleep_pct_M=sd(Sleep_pct_M))

write.table(f48h_ZT,"data/first48h_metabcagetable_mean.txt",sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#load data for plotting metab cage results-- Fig. 2e-i

f48h<-fread(metabdata48h_df)%>%mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))

f48h_ZT<-fread(metabdata48hmean_df)%>%mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))

###########################################################
#Motor Activity--Fig. 2e

MotorActivity_SEM_cleanest <- overZTplt(f48h,f48h$MotorActivity)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,2.3)) +
  labs(title='Motor Activity (m)',y="MotorActivity")
ggsave(paste0(results,"SFR24_0216_metab_cage_motoractivity_48h.pdf"), plot=MotorActivity_SEM_cleanest,height=3.5, width=3.5)

#ANOVA + LSD test
mamod <- aov(mn_MotorAct ~condition * zt_time_binned, data = f48h_ZT)
summary(mamod)

LSD_Test_ma<- (LSD.test(mamod, c("condition", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_time_binned"))
write.table(LSD_Test_ma,paste0(results,"SFR24_1027_metab_cage_MotorActivity_bytmpt_48h_anovaFisherLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#sleep--Fig. 2f

sleep_SEM_cleanest <-overZTplt(f48h,f48h$Sleep_pct_M)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,100)) +
  labs(title="Sleep")+
  ylab("Sleep_pct_M")+
  xlab("ZT")
ggsave(paste0(results,"SFR24_0216_metab_cage_sleep_48h.pdf"), plot=sleep_SEM_cleanest ,height=3.5, width=3.5)

#ANOVA + LSD test
sleepmod <- aov(Sleep_pct_M~ condition * zt_time_binned, data = f48h)
summary(sleepmod)

LSD_Test_sleep<- (LSD.test(sleepmod, c("condition", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_time_binned"))
write.table(LSD_Test_sleep,paste0(results,"SFR24_1027_metab_cage_sleep_bytmpt_48h_anovaFisherLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#RER--Fig. 2g

RER_M_SEM_cleanest <- overZTplt(f48h,f48h$RER_M)+
  scale_y_continuous(expand = c(0, 0), limits=c(0.5,1.0)) +
  labs(title="Mean respiratory exchange ratio (RER)",y="RER")
ggsave(paste0(results,"SFR25_0407_metab_cage_mRER_48h.pdf"), plot=RER_M_SEM_cleanest,height=3.5, width=3.5)

#ANOVA + LSD test
rermod <- aov(mn_RER~ condition * zt_time_binned, data = f48h_ZT)
summary(rermod)

LSD_Test_rer<- (LSD.test(rermod, c("condition", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_time_binned"))
write.table(LSD_Test_rer,paste0(results,"SFR24_1027_metab_cage_mRER_bytmpt_48h_anovaFisherLSD.tx"),sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#Mean Energy Expenditure--Fig. 2h

kcal_hr_M_SEM_cleanest <- overZTplt(f48h,f48h$kcal_hr_M)+
  scale_y_continuous(expand = c(0, 0), limits=c(0.25,0.75)) +
  labs(title="Mean energy expenditure") +
  ylab("Kcal/hr")+
  xlab("ZT")
ggsave(paste0(results,"SFR24_0216_metab_cage_energyexp_48h.pdf"), plot=kcal_hr_M_SEM_cleanest,height=3.5, width=3.5)

#ANOVA + LSD test
kcal_hr_mod <- aov(mn_kcal_hr_M~condition * zt_time_binned, data = f48h_ZT)
summary(kcal_hr_mod)

LSD_Test_kcal_hr<- (LSD.test(kcal_hr_mod, c("condition", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_time_binned"))
write.table(LSD_Test_kcal_hr,paste0(results,"SFR24_1027_metab_cage_energyexp_bytmpt_48h_anovaFisherLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#Food consumed, Kcal--Fig. 2i

KcalCons_SEM_cleanest <- overZTplt(f48h,f48h$kcalCons)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,4)) +
  labs(title="Food Consumption")+
  ylab("Avg Kcal Consumed") +
  xlab("ZT")
ggsave(paste0(results,"SFR24_0216_metab_cage_foodconsumed_kcal_f48h.pdf"), plot=KcalCons_SEM_cleanest ,height=3.5, width=3.5)

#ANOVA + LSD test
kcalCons_mod <- aov(mn_kcalCons~condition * zt_time_binned, data = f48h_ZT)
summary(kcalCons_mod)

LSD_Test_kcalCons<- (LSD.test(kcalCons_mod, c("condition", "zt_time_binned")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("condition","zt_time_binned"))
write.table(LSD_Test_kcalCons,paste0(results,"SFR24_1027_metab_cage_foodconsumed_kcal_bytmpt_48h_anovaFisherLSD.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
