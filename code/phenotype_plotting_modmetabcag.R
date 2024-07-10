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

###########################################################
#function

overZTplt<-function(df,yval){
  p<-ggplot(data=df, aes(x=zt_hour, y=yval, colour=condition, fill=condition)) +
    geom_rect(aes(xmin = 13, xmax = 21, ymin = -Inf, ymax = Inf),
              color="transparent", fill = 'gray88') +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
    scale_fill_manual(values = friendly_pal("ito_seven")) +
    scale_color_manual(values = friendly_pal("ito_seven")) +
    scale_x_continuous(breaks = seq(1, 23, by = 4 ), expand = c(0, 0)) +
    theme_classic()+
    theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  return(p)
}

overbinnedZT<-function(df,yval){
  p<-ggplot(data=df, aes(x=zt_time_binned, y=yval, colour=condition, fill=condition)) +
    geom_boxplot(color = "black",alpha=0.3)+
    geom_point(shape=21, color = "black", position=position_jitterdodge(),size = 1)+
    scale_fill_manual(values = friendly_pal("ito_seven")) +
    scale_color_manual(values = friendly_pal("ito_seven")) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic()+
    theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  return(p)
}

overbinnedLD<-function(df,yval){
  p<-ggplot(data=df, aes(x=condition, y=yval, colour=condition, fill=condition)) +
    geom_boxplot(color = "black",alpha=0.3)+
    geom_point(shape=21, color = "black", position=position_jitterdodge())+
    facet_grid(~phase)+
    scale_fill_manual(values = friendly_pal("ito_seven")) +
    scale_color_manual(values = friendly_pal("ito_seven")) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic()+
    theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  return(p)
}

run_ttest<-function(data,yval){
  if(yval=="RER"){
    x<-pairwise.wilcox.test(data$mn_RER, data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]]) 
  }
  else if(yval=="MotorActivity"){
    x<-pairwise.wilcox.test(data$mn_MotorAct,data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]])  
  }
  else if(yval=="WaterCons"){
    x<-pairwise.wilcox.test(data$mn_WaterCons,data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]])   
  }
  else if(yval=="kcal_hr_M"){
    x<-pairwise.wilcox.test(data$mn_kcal_hr_M,data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]])   
  }
  else if(yval=="FoodCons"){
    x<-pairwise.wilcox.test(data$mn_FoodCons,data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]])   
  }
  else if(yval=="kcalCons"){
    x<-pairwise.wilcox.test(data$mn_kcalCons,data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]])   
  }
  else{
    x<-pairwise.wilcox.test(data$mn_Sleep_pct_M,data$condition, p.adjust.method="fdr")
    df <- data.frame(NAvFA= x$p.value[[1]], NAvFT=x$p.value[[2]],FAvFT=x$p.value[[4]])  
  }
  
  return(df)
}


get_pvals<-function(df,time_compar,yval){
  if(time_compar=="ZT"){
    df_pval<-df %>% 
      group_by(zt_time_binned) %>% 
      nest()%>%
      mutate(pval = map2(data,zt_time_binned,  ~ run_ttest(.x,yval)))%>%
      dplyr::select(pval)%>%
      unnest()
  }
  else{
    df_pval<-df %>% 
      group_by(phase) %>% 
      nest()%>%
      mutate(pval = map2(data,phase,  ~ run_ttest(.x,yval)))%>%
      dplyr::select(pval)%>%
      unnest()
  }

  return(df_pval)
}
###########################################################
#metabolic cages

md <- read.csv("files_from_KF/STAM_TRF/Metabolic_Cage_Analysis/stam_trf_3_3_21/data/stam_trf_3_3_21_metadata.csv", 
                 colClasses = c(Animal="factor"),strip.white = TRUE)%>%
  mutate(Group=ifelse(is.na(Group),"NA",Group))%>%
  arrange(Group)

df<- read.csv("phenotype/stam_trf_run3_3_3_21_fogelson_macro_OneClickMacroV.2.53.2-slice1hr.mac_1_fromsable_fromsable.csv", na.strings = '.',
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

# create two vectors for start and end of the dark cycle
RER <- ggplot(data=df_wmd, aes(x=DateTime, y=RER_M, colour=condition, fill=condition)) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Mean respiratory exchange ratio (RER)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_RER_DateTime.pdf", plot=RER,height=3.5, width=6)
ggsave("phenotype/SFR24_0216_metab_cage_RER_DateTime_byAnimal.pdf", plot=RER,height=10, width=10)

xmins <- (df_wmd %>% filter(zt_hour=="4" & Animal==5))$DateTime
xmaxs <-c((df_wmd %>% filter(zt_hour=="23" & Animal==5))$DateTime,"2021-03-14 10:37:33 GMT")

RER <- ggplot(data=df_wmd, aes(x=DateTime, y=RER_M, colour=condition, fill=condition)) +
  annotate(geom = "rect",
           xmin = xmins,
           xmax = xmaxs,
           ymin = -Inf,
           fill="lightgrey",
           ymax = +Inf,
           alpha = 0.5) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Mean respiratory exchange ratio (RER)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_RER_DateTime_LD_byanimal.pdf", plot=RER,height=8, width=10)

Vo2 <- ggplot(data=df_wmd, aes(x=DateTime, y=VO2_M, colour=condition, fill=condition)) +
  annotate(geom = "rect",
           xmin = xmins,
           xmax = xmaxs,
           ymin = -Inf,
           fill="lightgrey",
           ymax = +Inf,
           alpha = 0.5) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="VO2") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_VO2_DateTime_LD_byanimal.pdf", plot=Vo2,height=8, width=10)

Vco2 <- ggplot(data=df_wmd, aes(x=DateTime, y=VCO2_M, colour=condition, fill=condition)) +
  annotate(geom = "rect",
           xmin = xmins,
           xmax = xmaxs,
           ymin = -Inf,
           fill="lightgrey",
           ymax = +Inf,
           alpha = 0.5) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="VCO2") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_VCO2_DateTime_LD_byanimal.pdf", plot=Vco2,height=8, width=10)


foodconsump <- ggplot() +
  geom_rect(
    aes(xmin=xmins, xmax=xmaxs, ymin=-Inf, ymax=Inf),
    fill= rep("lightgrey", 11),
    alpha=0.5) +
  geom_line(df_wmd, mapping=aes(x=DateTime, y=log10(FoodCons+1),colour=condition)) +
  #stat_summary(fun = mean, geom = "line") +
  #stat_summary(fun.data = mean_se, geom = "ribbon", alpha=0.25, aes(colour=NULL)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Food Consumption") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_foodcons_DateTime.pdf", plot=foodconsump,height=3.5, width=6)
ggsave("phenotype/SFR24_0216_metab_cage_foodcons_DateTime_12h.pdf", plot=foodconsump,height=3.5, width=6)

waterconsump <- ggplot() +
  geom_rect(
    aes(xmin=xmins, xmax=xmaxs, ymin=-Inf, ymax=Inf),
    fill= rep("lightgrey", 11),
    alpha=0.5) +
  geom_line(df_wmd, mapping=aes(x=DateTime, y=WaterCons,colour=condition)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Water Consumption") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/metab_cage_watercons_DateTime.pdf", plot=waterconsump,height=3.5, width=6)
ggsave("phenotype/metab_cage_watercons_DateTime_12h.pdf", plot=waterconsump,height=3.5, width=6)
###########################################################

#just plot the first 48 hours
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

f48h_ZT<-f48h%>%
  group_by(Animal,zt_time_binned,condition)%>%
  summarise(mn_RER=mean(RER_M),mn_MotorAct=mean(MotorActivity),mn_WaterCons=mean(WaterCons),
            mn_kcal_hr_M=mean(kcal_hr_M),mn_FoodCons=mean(FoodCons),mn_kcalCons=mean(kcalCons),
            mn_Sleep_pct_M=mean(Sleep_pct_M))

f48h_mm<-f48h%>%
  group_by(Animal,phase,condition)%>%
  summarise(mn_RER=mean(RER_M),mn_MotorAct=mean(MotorActivity),mn_WaterCons=mean(WaterCons),
            mn_kcal_hr_M=mean(kcal_hr_M),mn_FoodCons=mean(FoodCons),mn_kcalCons=mean(kcalCons),
            mn_Sleep_pct_M=mean(Sleep_pct_M))


xmins <- (f48h %>% filter(zt_hour=="12" & Animal==5))$DateTime
xmaxs <-(f48h %>% filter(zt_hour=="23" & Animal==5))$DateTime

RER <- ggplot(data=f48h, aes(x=DateTime, y=RER_M, colour=condition, fill=condition)) +
  annotate(geom = "rect",
           xmin = xmins,
           xmax = xmaxs,
           ymin = -Inf,
           fill="lightgrey",
           ymax = +Inf,
           alpha = 0.5) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="Mean respiratory exchange ratio (RER)") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_mRER_48h_byanimal.pdf", plot=RER,height=8, width=10)

VO2 <- ggplot(data=f48h, aes(x=DateTime, y=VO2_M, colour=condition, fill=condition)) +
  annotate(geom = "rect",
           xmin = xmins,
           xmax = xmaxs,
           ymin = -Inf,
           fill="lightgrey",
           ymax = +Inf,
           alpha = 0.5) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="VO2") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_VO2_48h_byanimal.pdf", plot=VO2,height=8, width=10)

VCO2 <- ggplot(data=f48h, aes(x=DateTime, y=VCO2_M, colour=condition, fill=condition)) +
  annotate(geom = "rect",
           xmin = xmins,
           xmax = xmaxs,
           ymin = -Inf,
           fill="lightgrey",
           ymax = +Inf,
           alpha = 0.5) +
  geom_line()+
  facet_wrap(~Animal, scales="free") +
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  labs(title="VO2") +
  theme_classic()+
  theme(legend.position = "top", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("phenotype/SFR24_0216_metab_cage_VCO2_48h_byanimal.pdf", plot=VCO2,height=8, width=10)


RER_M_SEM_cleanest <- overZTplt(f48h,f48h$RER_M)+
  scale_y_continuous(expand = c(0, 0), limits=c(0.5,1)) +
  labs(title="Mean respiratory exchange ratio (RER)",y="RER")
ggsave("phenotype/SFR24_0216_metab_cage_mRER_48h.pdf", plot=RER_M_SEM_cleanest,height=3.5, width=3.5)

RER_M_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_RER)+
  scale_y_continuous(expand = c(0, 0), limits=c(0.45,1)) +
  labs(title="Mean respiratory exchange ratio (RER)", x="ZT time",y="RER")
ggsave("phenotype/SFR24_0301_metab_cage_mRER_binned_48h.pdf", plot=RER_M_SEM_boxplot,height=3.5, width=6)

RER_pval <-get_pvals(f48h_ZT,"ZT","RER") 
write.table(RER_pval,"phenotype/SFR24_0301_metab_cage_mRER_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

RER_M_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_RER)+
  scale_y_continuous(expand = c(0, 0), limits=c(0.5,0.85)) +
  labs(title="Mean respiratory exchange ratio (RER)", x="ZT time",y="RER")
ggsave("phenotype/SFR24_0301_metab_cage_mRER_LD_48h.pdf", plot=RER_M_SEM_boxplotLD,height=3.5, width=4)

RER_pval2 <-get_pvals(f48h_mm,"phase","RER") 
write.table(RER_pval2,"phenotype/SFR24_0301_metab_cage_mRER_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

MotorActivity_SEM_cleanest <- overZTplt(f48h,f48h$MotorActivity)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title='Motor Activity (m)',y="MotorActivity")
ggsave("phenotype/SFR24_0216_metab_cage_motoractivity_48h.pdf", plot=MotorActivity_SEM_cleanest,height=3.5, width=3.5)

pairwise.wilcox.test(f48h_light$RER_M_mn,f48h_light$condition,p.adjust.method = "fdr")
# NA  FA 
# FA 0.2 -  
#   FT 0.2 0.2
pairwise.wilcox.test(f48h_dark$RER_M_mn,f48h_dark$condition,p.adjust.method = "fdr")
# NA   FA  
# FA 0.15 -   
#   FT 0.15 0.38

#get signif of Light vs. dark 

f48h_light<-f48h%>%filter(zt_hour<13)%>%
  group_by(Mouse_ID,condition)%>%
  summarise(MotorActivity_mn=mean(MotorActivity),
            kcal_hr_M_mn=mean(kcal_hr_M),
            kcalCons_mn=mean(kcalCons),
            Sleep_pct_M_mn=mean(Sleep_pct_M),
            RER_M_mn=mean(RER_M))

pairwise.wilcox.test(f48h_light$MotorActivity_mn,f48h_light$condition,p.adjust.method = "fdr")
# NA   FA  
# FA 0.29 -   
#   FT 0.50 0.29
kruskal.test(f48h_light$MotorActivity_mn,f48h_light$condition)

f48h_dark<-f48h%>%filter(zt_hour>12)%>%
  group_by(Mouse_ID,condition)%>%
  summarise(MotorActivity_mn=mean(MotorActivity),
            kcal_hr_M_mn=mean(kcal_hr_M),
            kcalCons_mn=mean(kcalCons),
            Sleep_pct_M_mn=mean(Sleep_pct_M),
            RER_M_mn=mean(RER_M))

pairwise.wilcox.test(f48h_dark$MotorActivity_mn,f48h_dark$condition,p.adjust.method = "fdr")
# NA   FA  
# FA 0.92 -   
#   FT 0.92 0.92
kruskal.test(f48h_dark$MotorActivity_mn,f48h_dark$condition)

MotorActivity_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_MotorAct)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title='Motor Activity (m)',y="MotorActivity")
ggsave("phenotype/SFR24_0301_metab_cage_motoractivity_binned_48h.pdf", plot=MotorActivity_SEM_boxplot,height=3.5, width=6)

MotorActivity_pval <-get_pvals(f48h_ZT,"ZT","MotorActivity") 
write.table(MotorActivity_pval,"phenotype/SFR24_0301_metab_cage_motoractivity_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

MotorActivity_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_MotorAct)+
  labs(title='Motor Activity (m)',y="MotorActivity")
ggsave("phenotype/SFR24_0301_metab_cage_motoractivity_binnedLD_48h.pdf", plot=MotorActivity_SEM_boxplotLD,height=3.5, width=4)

MotorActivity_pval2 <-get_pvals(f48h_mm,"phase","MotorActivity") 
write.table(MotorActivity_pval2,"phenotype/SFR24_0301_metab_cage_motoractivity_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

WaterCons_SEM_cleanest <- overZTplt(f48h,f48h$WaterCons)+
  scale_y_continuous(expand = c(0, 0))+
  labs(title='Mass of water consumed (g)', y="WaterCons")
ggsave("phenotype/SFR24_0216_metab_cage_masswaterconsum_48h.pdf", plot=WaterCons_SEM_cleanest,height=3.5, width=3.5)

WaterCons_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_WaterCons)+
  labs(title='Mass of water consumed (g)',y="WaterCons")
ggsave("phenotype/SFR24_0301_metab_cage_masswaterconsum_binned_48h.pdf", plot=WaterCons_SEM_boxplot,height=3.5, width=6)

WaterCons_pval <-get_pvals(f48h_ZT,"ZT","WaterCons") 
write.table(WaterCons_pval,"phenotype/SFR24_0301_metab_cage_masswaterconsum_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

WaterCons_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_WaterCons)+
  labs(title='Mass of water consumed (g)',y="WaterCons")
ggsave("phenotype/SFR24_0301_metab_cage_masswaterconsum_binnedLD_48h.pdf", plot=WaterCons_SEM_boxplotLD,height=3.5, width=4)

WaterCons_pval2 <-get_pvals(f48h_mm,"phase","WaterCons") 
write.table(WaterCons_pval2,"phenotype/SFR24_0301_metab_cage_masswaterconsum_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#Mean Energy Expenditure 
kcal_hr_M_SEM_cleanest <- overZTplt(f48h,f48h$kcal_hr_M)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title="Mean energy expenditure") +
  ylab("Kcal/hr")+
  xlab("ZT")

pairwise.wilcox.test(f48h_light$kcal_hr_M_mn,f48h_light$condition, p.adjust.method = "fdr")
# NA  FA 
# FA 0.5 -  
#   FT 0.5 0.5
pairwise.wilcox.test(f48h_dark$kcal_hr_M_mn,f48h_dark$condition, p.adjust.method = "fdr")
# NA    FA   
# FA 0.600 -    
#   FT 0.630 0.073

ggsave("phenotype/SFR24_0216_metab_cage_energyexp_48h.pdf", plot=kcal_hr_M_SEM_cleanest,height=3.5, width=3.5)

kcal_hr_M_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_kcal_hr_M)+
  labs(title="Mean energy expenditure") +
  ylab("Kcal/hr")+
  xlab("ZT")
ggsave("phenotype/SFR24_0301_metab_cage_energyexp_binned_48h.pdf", plot=kcal_hr_M_SEM_boxplot,height=3.5, width=6)

kcal_hr_M_pval <-get_pvals(f48h_ZT,"ZT","kcal_hr_M") 
write.table(kcal_hr_M_pval,"phenotype/SFR24_0301_metab_cage_energyexp_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

kcal_hr_M_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_kcal_hr_M)+
  labs(title="Mean energy expenditure") +
  ylab("Kcal/hr")+
  xlab("ZT")
ggsave("phenotype/SFR24_0301_metab_cage_energyexp_binnedLD_48h.pdf", plot=kcal_hr_M_SEM_boxplotLD,height=3.5, width=4)

kcal_hr_M_pval2 <-get_pvals(f48h_mm,"phase","kcal_hr_M") 
write.table(kcal_hr_M_pval2,"phenotype/SFR24_0301_metab_cage_energyexp_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#Food consumed, grams 
FoodCons_SEM_cleanest <- overZTplt(f48h,f48h$FoodCons)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title='Mass of food consumed (g)', y="FoodCons")
ggsave("phenotype/SFR24_0216_metab_cage_foodconsumed_g_48h.pdf", plot=FoodCons_SEM_cleanest,height=3.5, width=3.5)

FoodCons_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_FoodCons)+
  labs(title='Mass of food consumed (g)',y="FoodCons")
ggsave("phenotype/SFR24_0301_metab_cage_foodconsumed_g_binned_48h.pdf", plot=FoodCons_SEM_boxplot,height=3.5, width=6)

FoodCons_pval <-get_pvals(f48h_ZT,"ZT","FoodCons") 
write.table(FoodCons_pval,"phenotype/SFR24_0301_metab_cage_foodconsumed_g_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

FoodCons_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_FoodCons)+
  labs(title='Mass of food consumed (g)',y="FoodCons")
ggsave("phenotype/SFR24_0301_metab_cage_foodconsumed_g_binnedLD_48h.pdf", plot=FoodCons_SEM_boxplotLD,height=3.5, width=4)

FoodCons_pval2 <-get_pvals(f48h_mm,"phase","FoodCons") 
write.table(FoodCons_pval2,"phenotype/SFR24_0301_metab_cage_foodconsumed_g_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#Kcal consumed
KcalCons_SEM_cleanest <- overZTplt(f48h,f48h$kcalCons)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,4)) +
  labs(title="Food Consumption")+
  ylab("Avg Kcal Consumed") +
  xlab("ZT")
ggsave("phenotype/SFR24_0216_metab_cage_foodconsumed_kcal_f48h.pdf", plot=KcalCons_SEM_cleanest ,height=3.5, width=3.5)

KcalCons_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_kcalCons)+
  labs(title='Food Consumption',y="Kcal Consumed")
ggsave("phenotype/SFR24_0301_metab_cage_foodconsumed_kcal_binned_48h.pdf", plot=KcalCons_SEM_boxplot,height=3.5, width=6)

kcalCons_pval <-get_pvals(f48h_ZT,"ZT","kcalCons") 
write.table(kcalCons_pval,"phenotype/SFR24_0301_metab_cage_foodconsumed_kcal_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

KcalCons_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_kcalCons)+
  labs(title='Food Consumption',y="Kcal Consumed")
ggsave("phenotype/SFR24_0301_metab_cage_foodconsumed_kcal_binnedLD_48h.pdf", plot=KcalCons_SEM_boxplotLD,height=3.5, width=4)

kcalCons_pval2 <-get_pvals(f48h_mm,"phase","kcalCons") 
write.table(kcalCons_pval2,"phenotype/SFR24_0301_metab_cage_foodconsumed_kcal_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

f48h_mm_rm11<-f48h_mm%>%filter(Animal!=11)

KcalCons_SEM_boxplotLD <- overbinnedLD(f48h_mm_rm11,f48h_mm_rm11$mn_kcalCons)+
  labs(title='Food Consumption',y="Kcal Consumed")
ggsave("phenotype/SFR24_0301_metab_cage_foodconsumed_kcal_binnedLD_48h_rm11.pdf", plot=KcalCons_SEM_boxplotLD,height=3.5, width=4)

kcalCons_pval2 <-get_pvals(f48h_mm_rm11,"phase","kcalCons") 
write.table(kcalCons_pval2,"phenotype/SFR24_0301_metab_cage_foodconsumed_kcal_binnedLD_48h_pvalFDR_rm11.txt",sep = "\t",row.names = FALSE,quote=FALSE)

x<-pairwise.t.test(f48h_mm_rm11$mn_kcalCons,f48h_mm_rm11$condition, p.adjust.method="fdr")

pairwise.wilcox.test(f48h_light$kcalCons_mn,f48h_light$condition, p.adjust.method = "fdr")
# NA    FA   
# FA 1.000 -    
#   FT 0.018 0.018 #dont count FT
pairwise.wilcox.test(f48h_dark$kcalCons_mn,f48h_dark$condition, p.adjust.method = "fdr")
# NA   FA  
# FA 0.15 -   
#   FT 0.15 0.38

#Percent cumulative Kcal
percent_cumulative_KcalCons_SEM_cleanest <- overZTplt(f48h,f48h$percent_cumulative_Kcal)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1)) +
  labs(title="% of Daily Food Consumption")+
  ylab("% Kcal")+
  xlab("ZT")
ggsave("phenotype/SFR24_0216_metab_cage_cumulkcalconsumed_48h.pdf", plot=percent_cumulative_KcalCons_SEM_cleanest ,height=3.5, width=3.5)


#sleep
sleep_SEM_cleanest <-overZTplt(f48h,f48h$Sleep_pct_M)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(title="Sleep")+
  ylab("Sleep_pct_M")+
  xlab("ZT")
ggsave("phenotype/SFR24_0216_metab_cage_sleep_48h.pdf", plot=sleep_SEM_cleanest ,height=3.5, width=3.5)

sleep_SEM_boxplot <- overbinnedZT(f48h_ZT,f48h_ZT$mn_Sleep_pct_M)+
  labs(title='Sleep',y="Sleep_pct_M")
ggsave("phenotype/SFR24_0301_metab_cage_sleep_binned_48h.pdf", plot=sleep_SEM_boxplot,height=3.5, width=6)

sleep_SEM_pval <-get_pvals(f48h_ZT,"ZT","mn_Sleep_pct_M") 
write.table(sleep_SEM_pval,"phenotype/SFR24_0301_metab_cage_sleep_binned_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sleep_SEM_boxplotLD <- overbinnedLD(f48h_mm,f48h_mm$mn_Sleep_pct_M)+
  labs(title='Sleep',y="Sleep_pct_M")
ggsave("phenotype/SFR24_0301_metab_cage_sleep_binnedLD_48h.pdf", plot=sleep_SEM_boxplotLD,height=3.5, width=4)

sleep_SEM_pval2 <-get_pvals(f48h_mm,"phase","mn_Sleep_pct_M") 
write.table(sleep_SEM_pval2,"phenotype/SFR24_0301_metab_cage_sleep_binnedLD_48h_pvalFDR.txt",sep = "\t",row.names = FALSE,quote=FALSE)

pairwise.wilcox.test(f48h_light$Sleep_pct_M_mn,f48h_light$condition,p.adjust.method = "fdr")
# NA   FA  
# FA 0.30 -   
#   FT 0.63 0.25
pairwise.wilcox.test(f48h_dark$Sleep_pct_M_mn,f48h_dark$condition,p.adjust.method = "fdr")
# NA   FA  
# FA 0.92 -   
#   FT 0.92 0.92

# x<-f48h%>%filter(Animal==11 & zt_hour>9)
# ggplot(data=x, aes(x=zt_hour, y=kcalCons)) +
#   geom_line()+
#   geom_point()
