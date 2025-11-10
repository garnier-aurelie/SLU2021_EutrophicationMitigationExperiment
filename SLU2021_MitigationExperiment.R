Project_data_folder = "~/Desktop/_RESEARCH/PROJECTS/2021_SLU_Garnier/data/" ## create folder at the beginning
Project_figure_folder = "~/Desktop/_RESEARCH/PROJECTS/2021_SLU_Garnier/figures/" ## create folder at the beginning


### Load the libraries used for the analyses
library(plyr); library(tidyr)
library(ggplot2); library(ggpubr); library(forcats)
library(stats); library(rstatix) ; library(car); library(MASS)
library(lme4); library(performance) ; library(lmerTest); library(lmtest)
library(vegan)


### Load the experimental design 
expt_design2021 <- read.csv("expt_design.csv", header=T, stringsAsFactors = T)
expt_design2021$Mitigation_action_name <- factor(expt_design2021$Mitigation_action_name, 
                                                 levels=c("no mitigation", "fish removed", "nutrient stopped", "fish removed and nutrient stopped"))

color_mitigation <- c("#8DE60A","#1D2EAC","#FBB80B","#D30965") ## colors used in the figures for the mitigation actions


### Creation of the "dd_mesocosm" gathering all variables for the analysis for each tank
dd_mesocosm <- expt_design2021


################################################### 
###################################################  Light attenuation coefficient and Napieran coefficient
################################################### 

### >> Light attenuation coefficient from PAR measurements
raw_PAR <- read.csv(file=paste0(Project_data_folder,"raw_PARmeasurements.csv"), header=T, sep="\t", dec=",", stringsAsFactors = T)
raw_PAR$Date <- as.Date(raw_PAR$Date, format="%d/%m/%Y")

### calculation 
dd_Kd_coeff <- plyr::ddply(raw_PAR, .(Mesocosm,Date,Weather,In_the_shade), function(df){
  slope=summary(lm(log(df$PAR)~df$Depth_m))$coefficients[2]
  intercept=summary(lm(log(df$PAR)~df$Depth_m))$coefficients[1]
  return(c(abs(slope),intercept,slope))
}) 
names(dd_Kd_coeff)[c(5:7)] <- c("Kd","intercept","slope") 

## removing the first date (21/07/2021) - before the addition of the Sera Blackwater conditioner
mean_Kd = plyr::ddply(subset(dd_Kd_coeff, Date!="2021-07-21"), .(Mesocosm), summarise,
                      mean_Kd = mean(Kd, na.rm=T))

## add information to "dd_mesocosm"
dd_mesocosm <- merge(dd_mesocosm, mean_Kd, by="Mesocosm")

### >> Proportion of "In_the_shade = yes" 
dd_shade = plyr::ddply(subset(raw_PAR,Depth_m==0.1 & Date!="2021-07-21"), .(Mesocosm), 
                       function(df){
                         prop_shade = nrow(subset(df, In_the_shade=="yes")) / nrow(subset(df, In_the_shade=="yes" | In_the_shade=="no"));
                         return(round(prop_shade,3))})
names(dd_shade)[2]="prop_shade"

## add information to "dd_mesocosm"
dd_mesocosm <- merge(dd_mesocosm, dd_shade, by="Mesocosm")

##### Figure Suppl 3C --- influence of weather and shade on the PAR at 10cm depth
raw_PAR <- merge(expt_design2021[,c("Mesocosm","Light_treatment")], raw_PAR, by="Mesocosm")
ggplot2::ggplot(subset(raw_PAR, Depth_m==0.1), aes(x=Weather, y=PAR, fill=In_the_shade)) + geom_boxplot() +
  facet_wrap(. ~ Light_treatment)+
  labs(y="PAR at 10 cm below the surface")+
  theme_classic()

##### Figure Suppl 3D --- influence of the proportion of shade on the mean Kd coefficient
ggplot2::ggplot(dd_mesocosm, aes(x=Light_treatment, y=mean_Kd)) +
  geom_jitter(size=3,width=0.15, height=0, aes(fill=prop_shade), shape=21) +
  labs(x="Light treatment", title ="mean (over time) Kd per tank", y="Kd", fill="Proportion \nof shade")+
  scale_fill_gradient(low="white",high="black")+
  theme_classic()


  
### >> Absorbance measured at 420nm
raw_absorbance <- read.csv(file=paste0(Project_data_folder,"raw_absorbance.csv"))
raw_absorbance$Date <- as.Date(raw_absorbance$Date, format="%d/%m/%Y")

# Estimation of the blank for each method (5cm or 1cm cuvette)
blanck_1cm_cuvette <- mean(subset(raw_absorbance, Method=="1cm_cuvette" & ID=="MQ")$Abs_420)
blanck_5cm_cuvette <- mean(subset(raw_absorbance, Method=="5cm_flow_through_cuvette" & ID=="MQ")$Abs_420)

# Five samples (from the 6/08/2021) have two measurements, taking the mean of the measurements:
raw_absorbance2 <- plyr::ddply(.data=raw_absorbance, .(Date, ID, Method), 
                               summarize,
                               mean_abs = mean(Abs_420, na.rm=T))
## removal of the blanck
raw_absorbance2$Abs_420 <- ifelse(raw_absorbance2$Method=="1cm_cuvette",
                                  raw_absorbance2$mean_abs - blanck_1cm_cuvette, 
                                  raw_absorbance2$mean_abs - blanck_5cm_cuvette) 

## estimation of the Napieran absorption coefficient (cf. Renée's paper and Hu. Mueller-Karger and Zepp 2022 Limno&Oceano 47(4))
raw_absorbance2$Napieran_Coef <- ifelse(raw_absorbance2$Method=="1cm_cuvette",
                                        (raw_absorbance2$Abs_420*log(10))/0.01, 
                                        (raw_absorbance2$Abs_420*log(10))/0.05) 

## Quality check
raw_absorbance2 <- merge(expt_design2021, raw_absorbance2, by.x="Mesocosm", by.y="ID")
ggplot2::ggplot(subset(raw_absorbance2), aes(x=Date, y=Abs_420, group=Mesocosm)) + geom_point()+geom_line() +geom_vline(xintercept=as.Date(c("2021-07-22","2021-08-15")), colour="red") +
  geom_point(data=subset(raw_absorbance2, Date=="2021-08-25"), aes(x=Date, y=Abs_420), color="red")+
  facet_wrap(.~Light_treatment)
## For the samples on 25th August, we could not filter directly the samples before frozen. 
## It appears, from the figure above, that the samples were not filtered either before the spectro analysis.
## We remove these samples.

## Additionally, two measures on the 25th July (tank #5 - clear) and on the 27th July (tank #15 - dark) are outliers.
## We remove the data.

raw_absorbance3a <- raw_absorbance2[-which(raw_absorbance2$Date=="2021-08-25"),]
raw_absorbance3b <- raw_absorbance3a[-which(raw_absorbance3a$Date=="2021-07-25" & raw_absorbance3a$Mesocosm=="5"),]
raw_absorbance3 <- raw_absorbance3b[-which(raw_absorbance3b$Date=="2021-07-27" & raw_absorbance3b$Mesocosm=="15"),]
rm("raw_absorbance3b", "raw_absorbance3a")

ggplot2::ggplot(subset(raw_absorbance3), aes(x=Date, y=Abs_420, group=Mesocosm)) + geom_point()+geom_line() +geom_vline(xintercept=as.Date(c("2021-07-22","2021-08-15")), colour="red") +
  facet_wrap(.~Light_treatment)

# get the mean Napieran coefficient per mesocosm (without the data before the addition of water conditioner)
meanNapCoef = plyr::ddply(subset(raw_absorbance3, Date>"2021-07-21"), .(Mesocosm), summarize,
                          mean_NapieranCoef = mean(Napieran_Coef, na.rm=T))

## add information to "dd_mesocosm"
dd_mesocosm <- merge(dd_mesocosm, meanNapCoef, by="Mesocosm")

################################################### 
###################################################  Abiotic parameters measured by Aquaread
################################################### 

raw_abiotic <- read.csv(paste0(Project_data_folder,"raw_Aquaread.csv"), header=TRUE)
raw_abiotic$Date = as.Date(raw_abiotic$Date, format="%d/%m/%Y")

## Add the day of the experiment, which started (perturbation added) on the 22/07/2021.
raw_abiotic$Day_expt <- raw_abiotic$Date - as.Date("2021-07-22")

abiotic_expt <- merge(expt_design2021, raw_abiotic, by="Mesocosm")

## For the temperature, we compare the air temperature measured in Örskär, an island close to the experimental site. 
# the SMHI's station "Örskär" is located Latitude: 60.5256, longitude: 18.3729, altitude: 7m. 
# data can be downloaded from the following link:
# https://www.smhi.se/data/hitta-data-for-en-plats/ladda-ner-vaderobservationer/airtemperatureMean24h/108320
temperature_Orskar <- read.csv(file=paste0(Project_data_folder,"temperature_Orskar_SE.csv"))
temperature_Orskar$Date = as.Date(temperature_Orskar$Date, format="%d/%m/%Y")
temperature_Orskar$Day_expt <- temperature_Orskar$Date - as.Date("2021-07-22")

# information to add on the figure
anno <- data.frame(x=c(0,23), y=c(-Inf, -Inf), Light_treatment=c("clear","clear"),labs=c("Jul-22", "Aug-14"))

fig_temperature = ggplot(abiotic_expt, aes(x=Day_expt, y=Temp..C.)) + 
  geom_vline(xintercept=c(0.5, 23.5), colour="red", linetype=2)+
  geom_text(data=anno, aes(x=x, y=y, label=labs), angle=90, hjust=-0.2, vjust=0)+
 geom_point(color="grey50")+
  geom_line(aes(group=Mesocosm),color="grey50")+
  geom_point(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=Temp..C., group=Mesocosm, color=Mitigation_action_name))+
  geom_line(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=Temp..C., group=Mesocosm, color=Mitigation_action_name))+
  geom_line(data=subset(temperature_Orskar, Date<"2021-09-05"), aes(x=Day_expt, y=AverageDay_Air_temperature_degreeC), color="black", linewidth=1.2)+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  facet_wrap(.~Light_treatment, ncol=3)+
  labs(y="Temperature (°C)", x="Day")+
  theme(legend.position = "bottom")

fig_salinity = ggplot(abiotic_expt, aes(x=Day_expt, y=SAL..PSU.)) + 
  geom_vline(xintercept=c(0.5, 23.5), colour="red", linetype=2)+
  geom_point(color="grey50")+
  geom_line(aes(group=Mesocosm),color="grey50")+
  geom_point(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=SAL..PSU., group=Mesocosm, color=Mitigation_action_name))+
  geom_line(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=SAL..PSU., group=Mesocosm, color=Mitigation_action_name))+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  facet_wrap(.~Light_treatment, ncol=3)+
  labs(y="Salinity (PSU)", x="Day")+
  theme(legend.position = "bottom")

fig_DOsat = ggplot(abiotic_expt, aes(x=Day_expt, y=DO....Sat.)) + 
  geom_vline(xintercept=c(0.5, 23.5), colour="red", linetype=2)+
  geom_point(color="grey50")+
  geom_line(aes(group=Mesocosm),color="grey50")+
  geom_point(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=DO....Sat., group=Mesocosm, color=Mitigation_action_name))+
  geom_line(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=DO....Sat., group=Mesocosm, color=Mitigation_action_name))+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  facet_wrap(.~Light_treatment, ncol=3)+
  labs(y="Dissolved Oxygen (%)", x="Day")+
  theme(legend.position="bottom")

fig_pH = ggplot(abiotic_expt, aes(x=Day_expt, y=pH)) + 
  geom_vline(xintercept=c(0.5, 23.5), colour="red", linetype=2)+
  geom_point(color="grey50")+
  geom_line(aes(group=Mesocosm),color="grey50")+
  geom_point(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=pH, group=Mesocosm, color=Mitigation_action_name))+
  geom_line(data=subset(abiotic_expt, Date > "2021-08-12"), aes(x=Day_expt, y=pH, group=Mesocosm, color=Mitigation_action_name))+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  facet_wrap(.~Light_treatment, ncol=3)+
  labs(y="pH", x="Day")+
  theme(legend.position = "bottom")

ggpubr::ggarrange(fig_temperature, fig_salinity, fig_DOsat, fig_pH,
          nrow=4, ncol=1, align = "hv", labels="AUTO",
          common.legend = T, legend = "bottom")
ggplot2::ggsave(paste0(Project_figure_folder,"fig_abiotic.png"), units="cm", width=20, height=30, dpi=300)

################### Table Suppl 1
### Linear regression on abiotic parameters (salinity, DO and pH) at the end of the eutrophication phase

abiotic_endEutrophication <- subset(abiotic_expt, Date=="2021-08-09")
abiotic_endEutrophication <- merge(abiotic_endEutrophication, meanNapCoef, by="Mesocosm")


## Salinity
ggplot(abiotic_endEutrophication, aes(x=mean_NapieranCoef, y=SAL..PSU.)) + geom_point() + 
  geom_smooth(method="lm")

lm_Salinity_endEutrophication <- lm(SAL..PSU. ~ mean_NapieranCoef, data=abiotic_endEutrophication)
stats::shapiro.test(residuals(lm_Salinity_endEutrophication)) 
plot(lm_Salinity_endEutrophication)
summary(lm_Salinity_endEutrophication)

## Dissolved oxygen (% saturation)
ggplot(abiotic_endEutrophication, aes(x=mean_NapieranCoef, y=DO....Sat.)) + geom_point() + 
  geom_smooth(method="lm")

lm_DO_endEutrophication <- lm(DO....Sat. ~ mean_NapieranCoef, data=abiotic_endEutrophication)
stats::shapiro.test(residuals(lm_DO_endEutrophication)) 
plot(lm_DO_endEutrophication)
summary(lm_DO_endEutrophication)

## pH
ggplot(abiotic_endEutrophication, aes(x=mean_NapieranCoef, y=pH)) + geom_point() + 
  geom_smooth(method="lm")

lm_pH_endEutrophication <- lm(pH ~ mean_NapieranCoef, data=abiotic_endEutrophication)
stats::shapiro.test(residuals(lm_pH_endEutrophication)) 
plot(lm_pH_endEutrophication)
summary(lm_pH_endEutrophication)

################### Table Suppl 2
### Linear mixed effect modal on abiotic parameters (salinity, DO and pH) during the mitigation phase

abiotic_Mitigation = subset(abiotic_expt, Date >= "2021-08-13" & Date <="2021-09-02")
abiotic_Mitigation <- merge(abiotic_Mitigation, meanNapCoef, by="Mesocosm")


## salinity
hist(abiotic_Mitigation$SAL..PSU.)
lmer_Salinity_P2 <- lme4::lmer(SAL..PSU. ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                            data=abiotic_Mitigation,
                            REML=TRUE)

performance::check_model(lmer_Salinity_P2,
                         check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_Salinity_P2, check = "reqq")

mod_Salinity_Mitigation=lmerTest::lmer(SAL..PSU. ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                                 data=abiotic_Mitigation,
                                 REML=TRUE); anova(mod_Salinity_Mitigation)


## DO
hist(abiotic_Mitigation$DO....Sat.)
lmer_DOsat_P2 <- lme4::lmer(DO....Sat. ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                            data=abiotic_Mitigation,
                            REML=TRUE)

performance::check_model(lmer_DOsat_P2,
                         check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_DOsat_P2, check = "reqq")

mod_DO_Mitigation=lmerTest::lmer(DO....Sat. ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                   data=abiotic_Mitigation,
                   REML=TRUE) ; anova(mod_DO_Mitigation)

## pH
hist(abiotic_Mitigation$pH)
lmer_pH_P2 <- lme4::lmer(pH ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                            data=abiotic_Mitigation,
                            REML=TRUE)

performance::check_model(lmer_pH_P2,
                         check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_pH_P2, check = "reqq")

mod_pH_Mitigation=lmerTest::lmer(pH ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                                 data=abiotic_Mitigation,
                                 REML=TRUE); anova(mod_pH_Mitigation)

################################################### 
###################################################  Log-Response ratio for pelagic chlorophyll-a
################################################### 

## Pelagic chlorophyll-a concentrations
raw_chl_pel <- read.csv(file=paste0(Project_data_folder,"raw_chlorophyll_pelagic.csv"))
raw_chl_pel$Date <- as.Date(raw_chl_pel$Date, format="%d/%m/%Y")

## Estimate the chl-a concentration from measurements
# calibration curve from UMF data : concentration = 1.2221 + 0.1182*fluorescence
raw_chl_pel$conc_Pelagic_chla <- (1.2221 + (0.1182*raw_chl_pel$Mean_fluorescence)) * (raw_chl_pel$Volume_ethanol_ml/ raw_chl_pel$Volume_filtered_ml) * raw_chl_pel$Dilution 


# For the estimation of chl-a at the end of the eutrophication phase, we took the average of the last two measurements (from 31st July to 9th August 2021)
# The end of the mitigation phase was the 2nd September 2021

## end of Eutrophication phase (P1)
chl_pel_endP1 = plyr::ddply(subset(raw_chl_pel, Date<="2021-08-09" & Date>="2021-07-31"), 
                                  .(Mesocosm), summarize,
                                  conc_Pelagic_chla_P1 = mean(conc_Pelagic_chla, na.rm=T))
## end of Mitigation phase (P2)
chl_pel_endP2 = subset(raw_chl_pel, Date=="2021-09-02")[,c("Mesocosm","conc_Pelagic_chla")]
names(chl_pel_endP2)[2] = "conc_Pelagic_chla_P2"

## Get the log response ratio
LRR_chl_pel <- merge(chl_pel_endP1, chl_pel_endP2, by="Mesocosm")
LRR_chl_pel$LRR_chl_pel <-log(LRR_chl_pel$conc_Pelagic_chla_P2/LRR_chl_pel$conc_Pelagic_chla_P1)

dd_mesocosm <- merge(dd_mesocosm, LRR_chl_pel, by="Mesocosm")

fig_LRR_chlPel_grad = ggplot(dd_mesocosm, aes(x=mean_NapieranCoef, y=LRR_chl_pel, color=Mitigation_action_name, fill=Mitigation_action_name)) + 
  geom_point()+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  geom_hline(yintercept = 0, linetype=2)+
  geom_smooth(method="lm", alpha=.1)+
  ylab("LRR pelagic chl-a") +
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "bottom"); fig_LRR_chlPel_grad

## Ancova

## linearity assumption
ggpubr::ggscatter(
  dd_mesocosm, x = "mean_NapieranCoef", y = "LRR_chl_pel",
  facet.by  = c("Nutrient_stop", "Fish_removed"), 
  short.panel.labs = FALSE
)+
  ggplot2::stat_smooth(method = "loess", span = 0.9)

## homogeneity of regression slopes
dd_mesocosm %>%
  tidyr::unite(col = "group", Nutrient_stop, Fish_removed) %>%
  rstatix::anova_test(LRR_chl_pel ~ group*mean_NapieranCoef)

## normality of residuals
model <- lm(LRR_chl_pel ~ mean_NapieranCoef * Nutrient_stop * Fish_removed, data=dd_mesocosm)
model.metrics <- augment(model)  # Remove details
head(model.metrics,3)
shapiro_test(model.metrics$.resid)

## homogeneity of variances
rstatix::levene_test(.resid ~ Nutrient_stop * Fish_removed, data = model.metrics)

## outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

res.aov <- dd_mesocosm %>% 
  rstatix::anova_test(LRR_chl_pel ~ mean_NapieranCoef * Nutrient_stop * Fish_removed)
rstatix::get_anova_table(res.aov)


################################################### 
###################################################  Log-Response ratio for benthic chlorophyll-a
################################################### 


## Benthic chlorophyll-a concentrations
raw_chl_ben <- read.csv(file=paste0(Project_data_folder,"raw_chlorophyll_benthic.csv"))
raw_chl_ben$Date <- as.Date(raw_chl_ben$Date, format="%d/%m/%Y")

# estimate the chl-a concentration from measurements:
raw_chl_ben$conc_ugL_in_EtOH <- (1.2221 + (0.1182*raw_chl_ben$Mean_fluorescence))*raw_chl_ben$Dilution
# estimate the chl-a quantity (mg/m2):
## 1. total amount (µg) of chl-a extracted from the subsample in the 30ml of Ethanol
raw_chl_ben$amount_ug_in_subsample <- raw_chl_ben$conc_ugL_in_EtOH * raw_chl_ben$Volume_ethanol_ml * 10^-3 
## 2. extrapolate the quantity to the entire core sampled
raw_chl_ben$amount_ug_in_Total_core <- raw_chl_ben$amount_ug_in_subsample * (raw_chl_ben$Dryweight_total_sediment_g / round(raw_chl_ben$Dryweight_subsample_g, digits=3))
## 3. adjust to the area of the core (2.9cm diameter)
raw_chl_ben$mgm2_Benthic_chla <- raw_chl_ben$amount_ug_in_Total_core * 10^-3 / (pi*0.0145*0.0145) 

# For the end of the eutrophication phase, we spread the sampling during the 10th, 11th and 12th August
# The end of the mitigation phase was the 3rd and 4th September 2021

## end of Eutrophication phase (P1)
chl_ben_endP1 = subset(raw_chl_ben, Date>="2021-08-10" & Date <="2021-08-12")[,c("Mesocosm","mgm2_Benthic_chla")]
names(chl_ben_endP1)[2] = "mgm2_Benthic_chla_P1"
## end of Mitigation phase (P2)
chl_ben_endP2 = subset(raw_chl_ben, Date>="2021-09-03")[,c("Mesocosm","mgm2_Benthic_chla")]
names(chl_ben_endP2)[2] = "mgm2_Benthic_chla_P2"

## Get the log response ratio
LRR_chl_ben <- merge(chl_ben_endP1, chl_ben_endP2, by="Mesocosm")
LRR_chl_ben$LRR_chl_ben <-log(LRR_chl_ben$mgm2_Benthic_chla_P2/LRR_chl_ben$mgm2_Benthic_chla_P1)

dd_mesocosm <- merge(dd_mesocosm, LRR_chl_ben, by="Mesocosm")

fig_LRR_chlBen_grad = ggplot(dd_mesocosm, aes(x=mean_NapieranCoef, y=LRR_chl_ben, color=Mitigation_action_name, fill=Mitigation_action_name)) + 
  geom_point()+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  geom_hline(yintercept = 0, linetype=2)+
  geom_smooth(method="lm", alpha=.1)+
  ylab("LRR benthic chl-a") +
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "bottom"); fig_LRR_chlBen_grad

## Ancova

## linearity assumption
ggpubr::ggscatter(
  dd_mesocosm, x = "mean_NapieranCoef", y = "LRR_chl_ben",
  facet.by  = c("Nutrient_stop", "Fish_removed"), 
  short.panel.labs = FALSE
)+
  ggplot2::stat_smooth(method = "loess", span = 0.9)

## homogeneity of regression slopes
dd_mesocosm %>%
  tidyr::unite(col = "group", Nutrient_stop, Fish_removed) %>%
  rstatix::anova_test(LRR_chl_ben ~ group*mean_NapieranCoef)

## normality of residuals
model <- lm(LRR_chl_ben ~ mean_NapieranCoef * Nutrient_stop * Fish_removed, data=dd_mesocosm)
model.metrics <- augment(model)  # Remove details
head(model.metrics,3)
rstatix::shapiro_test(model.metrics$.resid)

## homogeneity of variances
rstatix::levene_test(.resid ~ Nutrient_stop * Fish_removed, data = model.metrics)

## outliers
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

res.aov <- dd_mesocosm %>% 
  rstatix::anova_test(LRR_chl_ben ~ mean_NapieranCoef * Nutrient_stop * Fish_removed)
rstatix::get_anova_table(res.aov)


# FIGURE 2 - CHL-a
ggpubr::ggarrange(fig_LRR_chlPel_grad, fig_LRR_chlBen_grad,
          nrow=1, ncol=2, align = "hv", labels="AUTO",
          common.legend = T, legend = "bottom")
ggplot2::ggsave(paste0(Project_figure_folder,"Figure2_Chla.png"), units="cm", width=20, height=10, dpi=300)
ggplot2::ggsave(paste0(Project_figure_folder,"Figure2_Chla.tiff"), units="cm", width=20, height=10, dpi=300, compression = 'lzw')


################################################### 
###################################################  Nutrients - Total Nitrogen and Total Phosphorus
################################################### 

raw_nutrients <- read.csv(file=paste0(Project_data_folder,"raw_nutrients.csv"), header=TRUE)

# some times, two concentrations were estimated per sample. We took the average. 
mean_nutrients <- plyr::ddply(raw_nutrients, .(Mesocosm, Date, Variable), summarize,
                        conc = mean(Concentration, rm.na=TRUE))
nutrients <- tidyr::spread(mean_nutrients, key=Variable, value = conc)
names(nutrients)[c(3,4)] <- c("Total_Nitrogen", "Total_Phosphorus")
nutrients$Date <- as.Date(nutrients$Date, format="%d/%m/%Y")
nutrients$Days <- as.numeric(nutrients$Date - as.Date("2021-07-22"))

nutrients_NapCoef <- merge(nutrients, meanNapCoef, by="Mesocosm") ## lost information regarding nutrients at Forsmark
nutrients_NapCoef <- merge(expt_design2021, nutrients_NapCoef, by="Mesocosm")
nutrients_NapCoef$Mesocosm <- as.factor(nutrients_NapCoef$Mesocosm)

Mean_SD_nutrients <- plyr::ddply(nutrients_NapCoef, .(Date, Days, Light_treatment, Mitigation_action_name), summarize,
                                 meanN = mean(Total_Nitrogen, na.rm=T), 
                                 sdN = sd(Total_Nitrogen, na.rm=T),
                                 meanP = mean(Total_Phosphorus, na.rm=T), 
                                 sdP = sd(Total_Phosphorus, na.rm=T))

initial_TotNitrogen = mean(subset(nutrients, Date=="2021-07-21")$Total_Nitrogen, rm.na=T) ## mean of total nitrogen before the addition of nutrients
initial_TotPhosphorus = mean(subset(nutrients, Date=="2021-07-21")$Total_Phosphorus, rm.na=T) ## mean of total phosphorus before the addition of nutrients

### Fig. 3A - Total nitrogen over time during the Mitigation phase (P2)

fig_totN_Light_ttt <-  ggplot(subset(Mean_SD_nutrients, Date>="2021-08-09"), 
                              aes(x=Days, y=meanN, group=Mitigation_action_name, 
                                  color=Mitigation_action_name, fill=Mitigation_action_name))+
  geom_ribbon(aes(x=Days, ymin=meanN-sdN, ymax=meanN+sdN), alpha=.3, color=NA)+
  geom_line(linewidth=1)+
  facet_wrap(.~Light_treatment)+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  labs(y="Total Nitrogen (µg/L)")+
  geom_hline(yintercept=initial_TotNitrogen, color="black", linetype=2)+
  geom_point(data=subset(Mean_SD_nutrients, Date=="2021-08-09"), fill="white", shape=21)+
  theme_classic()+
  theme(legend.position="bottom");fig_totN_Light_ttt

### linear mixed models for Total Nitrogen over time during the Mitigation phase (P2) -- Table 2
lmer_TotN_P2 <- lme4::lmer(log(Total_Nitrogen) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                        data=subset(nutrients_NapCoef, Date>="2021-08-09"),
                        REML=TRUE)

performance::check_model(lmer_TotN_P2, 
            check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_TotN_P2, check = "reqq")
stats::shapiro.test(residuals(lmer_TotN_P2))
plot(lmer_TotN_P2)

mod=lmerTest::lmer(log(Total_Nitrogen) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                   data=subset(nutrients_NapCoef, Date>="2021-08-09"),
                   REML=TRUE); anova(mod)


### Fig. 3B - Total phosphorus over time during the Mitigation phase (P2)

fig_totP_Light_ttt <-  ggplot(subset(Mean_SD_nutrients, Date>="2021-08-09"), 
                              aes(x=Days, y=meanP, group=Mitigation_action_name, 
                                  color=Mitigation_action_name, fill=Mitigation_action_name))+
  geom_ribbon(aes(x=Days, ymin=meanP-sdP, ymax=meanP+sdP), alpha=.3, color=NA)+
  geom_line(linewidth=1)+
  facet_wrap(.~Light_treatment)+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  labs(y="Total Phosphorus (µg/L)")+
  geom_hline(yintercept=initial_TotPhosphorus, color="black", linetype=2)+
  geom_point(data=subset(Mean_SD_nutrients, Date=="2021-08-09"), fill="white", shape=21)+
  theme_classic()+
  theme(legend.position="bottom");fig_totP_Light_ttt

### linear mixed models for Total Phosphorus over time during the Mitigation phase (P2) -- Table 2
lmer_TotP_P2 <- lme4::lmer(Total_Phosphorus ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                           data=subset(nutrients_NapCoef, Date>="2021-08-09"),
                           REML=TRUE)

performance::check_model(lmer_TotP_P2, 
            check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_TotP_P2, check = "reqq")
stats::shapiro.test(residuals(lmer_TotP_P2))
plot(lmer_TotP_P2)

mod=lmerTest::lmer(Total_Phosphorus ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) * as.factor(Date) + (1|Mesocosm) ,
                   data=subset(nutrients_NapCoef, Date>="2021-08-09"),
                   REML=TRUE); anova(mod)

### Figure 3. 
ggpubr::ggarrange(fig_totN_Light_ttt, 
          fig_totP_Light_ttt, 
          nrow=2, ncol=1, align = "hv", labels="AUTO",
          common.legend = T, legend = "bottom")
ggplot2::ggsave(paste0(Project_figure_folder,"Figure3_Nutrients.png"), units="cm", width=20, height=15, dpi=300)



### Fig. S5A - Total nitrogen over time during the Eutrophication phase (P1)

Nitrogen_P1 <- ggplot(subset(nutrients_NapCoef, Date >"2021-07-21"&Date<="2021-08-09"), aes(x=Days, y=Total_Nitrogen, group=Mesocosm, color=mean_NapieranCoef))+
  geom_line(aes(linewidth=mean_NapieranCoef))+
  scale_linewidth(range = c(0.3, 2))+
  ylim(200,1500)+
  guides(linewidth=FALSE)+
  scale_color_gradient(low = "grey", high = "coral4",
                       name = expression(Absorption~coefficient~(m^{-1})))+
  labs(y="Total Nitrogen (µg/L)")+
  geom_hline(yintercept=initial_TotNitrogen, color="black", linetype=2)+
  theme_classic()+
  theme(legend.position="bottom") ;Nitrogen_P1

### Fig S5B - Total nitrogen at the end of Eutrophication phase
Nitrogen_endP1 <- ggplot(subset(nutrients_NapCoef, Date=="2021-08-09"), aes(x=mean_NapieranCoef, y=Total_Nitrogen)) + 
  ylim(200,1500)+
  geom_point()+ geom_smooth(method="lm", color="black")+
  labs(y="Total Nitrogen (µg/L) \n(end eutrophication phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_classic() ;Nitrogen_endP1

lm_Nitrogen_endEutroph <- lm(log(Total_Nitrogen) ~ mean_NapieranCoef, 
                               data=subset(nutrients_NapCoef, Date=="2021-08-09"))
stats::shapiro.test(residuals(lm_Nitrogen_endEutroph))
summary(lm_Nitrogen_endEutroph)

### Fig. S5C - Total phosphorus over time during the Eutrophication phase (P1)

Phosphorus_P1 <- ggplot(subset(nutrients_NapCoef, Date >"2021-07-21"&Date<="2021-08-09"), aes(x=Days, y=Total_Phosphorus, group=Mesocosm, color=mean_NapieranCoef))+
  geom_line(aes(linewidth=mean_NapieranCoef))+
  scale_linewidth(range = c(0.3, 2))+
  ylim(0,200)+
  guides(linewidth=FALSE)+
  scale_color_gradient(low = "grey", high = "coral4",
                       name = expression(Absorption~coefficient~(m^{-1})))+
  labs(y="Total Phosphorus (µg/L)")+
  geom_hline(yintercept=initial_TotPhosphorus, color="black", linetype=2)+
  theme_classic()+
  theme(legend.position="bottom") ;Phosphorus_P1

### Fig S5D - Total phosphorus at the end of Eutrophication phase
Phosphorus_endP1 <- ggplot(subset(nutrients_NapCoef, Date=="2021-08-09"), aes(x=mean_NapieranCoef, y=Total_Phosphorus)) + 
  ylim(0,200)+
  geom_point()+ geom_smooth(method="lm", color="black")+
  labs(y="Total Phosphorus (µg/L) \n(end eutrophication phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_classic() ;Phosphorus_endP1

lm_Phosphorus_endEutroph <- lm(log(Total_Phosphorus) ~ mean_NapieranCoef, 
                               data=subset(nutrients_NapCoef, Date=="2021-08-09"))
stats::shapiro.test(residuals(lm_Phosphorus_endEutroph))
summary(lm_Phosphorus_endEutroph)

### Table S3 - Mixed linear models of total nitrogen and phosphorus concentrations over time during the eutrophication phases. 
############ Nitrogen
lmer_TotN_P1 <- lme4::lmer(log(Total_Nitrogen) ~ mean_NapieranCoef * as.factor(Date) + (1|Mesocosm) ,
                           data=subset(nutrients_NapCoef, Date<="2021-08-09"& Date>"2021-07-21"),
                           REML=TRUE)
performance::check_model(lmer_TotN_P1, 
            check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_TotN_P1, check = "reqq")

mod_N=lmerTest::lmer(log(Total_Nitrogen) ~ mean_NapieranCoef * as.factor(Date) + (1|Mesocosm) ,
                     data=subset(nutrients_NapCoef, Date<="2021-08-09"& Date>"2021-07-21"),
                     REML=TRUE) ; anova(mod_N)

############ Phosphorus
lmer_TotP_P1 <- lme4::lmer(log(Total_Phosphorus) ~ mean_NapieranCoef * as.factor(Date) + (1|Mesocosm) ,
                           data=subset(nutrients_NapCoef, Date<="2021-08-09"& Date>"2021-07-21"),
                           REML=TRUE)
performance::check_model(lmer_TotP_P1,
                         check = c("linearity", "homogeneity", "qq", "outliers"))
performance::check_model(lmer_TotP_P1, check = "reqq")

mod_P=lmerTest::lmer(log(Total_Phosphorus) ~ mean_NapieranCoef * as.factor(Date) + (1|Mesocosm) ,
                   data=subset(nutrients_NapCoef, Date<="2021-08-09"& Date>"2021-07-21"),
                   REML=TRUE) ; anova(mod_P)



################################################### 
###################################################  Zooplankton
################################################### 

raw_zooplankton_data <- read.csv(file=paste0(Project_data_folder,"raw_zooplankton_biomass.csv"))
str(raw_zooplankton_data)

raw_zooplankton_data$Date <- as.Date(raw_zooplankton_data$Date, format="%d/%m/%Y")
raw_zooplankton_data <- merge(expt_design2021, raw_zooplankton_data, by="Mesocosm")

## Estimation of total zooplankton biomass (in µg/L)
# 1. The rotifers Keratella were counted twice (i.e., 'Keratella' as total all genera and "K.quadrata, etc" for each genus)
# To have a better estimation of total biomass, we used the biomass estimated from each genus 
# Indeed, the Keratella's biomass is based on mean_size of all Keratella 

# 2. We removed as well the Chironomidae (not zooplankton)

total_biomass <- plyr::ddply(.data=subset(raw_zooplankton_data, Taxa!="Keratella" & Taxa!="Chironomid"), .(Mesocosm,Date), summarize,
                             total_biomass = sum(Sample_biomass_µg.L, na.rm=T))

total_biomassP2a <- subset(total_biomass, Date=="2021-09-02")
total_biomassP2 <- total_biomassP2a[,-2]
names(total_biomassP2)[2] <- "Total_ZP_biomassP2"
dd_mesocosm <- merge(dd_mesocosm, total_biomassP2, by="Mesocosm")

total_biomassP1a <- subset(total_biomass, Date=="2021-08-09")
total_biomassP1 <- total_biomassP1a[,-2]
names(total_biomassP1)[2] <- "Total_ZP_biomassP1"
dd_mesocosm <- merge(dd_mesocosm, total_biomassP1, by="Mesocosm")

# >> Biomass per group
group_biomass1 <- plyr::ddply(.data=subset(raw_zooplankton_data, Taxa!="Keratella" & Taxa!="Chironomid"), .(Mesocosm,Date,Group), summarize,
                              biomass = sum(Sample_biomass_µg.L, na.rm=T))

group_biomass <- stats::reshape(group_biomass1, idvar = c("Mesocosm","Date"), v.names = "biomass", sep="_", timevar = "Group", direction = "wide")

group_biomassP2a <- subset(group_biomass, Date=="2021-09-02")
group_biomassP2b <- group_biomassP2a[,-2]
names(group_biomassP2b) <- c("Mesocosm", "biomass_Copepod_P2", "biomass_Rotifer_P2", "biomass_Cladoceran_P2")
group_biomassP2 <- group_biomassP2b

dd_mesocosm <- merge(dd_mesocosm, group_biomassP2, by="Mesocosm")


## Figure 4A - Total zooplankton biomass ~ absorption coefficient 
fig_ZPbiom <- ggplot(dd_mesocosm, aes(x=mean_NapieranCoef, y=(Total_ZP_biomassP2), color=Mitigation_action_name, fill=Mitigation_action_name)) +
  geom_point()+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  geom_smooth(method="lm", alpha=.1)+
  ylab("Total zooplankton biomass (µg/L) \n(end of mitigation phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "bottom"); fig_ZPbiom

## Model // gradient 
model_ZPbiomass <- lm(log(Total_ZP_biomassP2) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed),
                      data=dd_mesocosm)
##Check assumptions:
plot(model_ZPbiomass)
hist(model_ZPbiomass$residuals)
stats::shapiro.test(residuals(model_ZPbiomass)) # check normality (model residuals) - OK if p.value > 0.05
lmtest::bptest(log(Total_ZP_biomassP2) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) , data=dd_mesocosm) # check homoscedasticity (Breusch-Pagan test) - OK if p.value > 0.05

## get the results (Table 3)
car::Anova(model_ZPbiomass, type=2)
summary(model_ZPbiomass)


## Figure 4B - Copepods
fig_Copepod <- ggplot(dd_mesocosm, aes(x=mean_NapieranCoef, y=biomass_Copepod_P2, color=Mitigation_action_name, fill=Mitigation_action_name)) +
  geom_point()+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  geom_smooth(method="lm", alpha=.1)+
  ylab("Copepod biomass (µg/L) \n(end of mitigation phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "bottom") ; fig_Copepod


model_Copepod <- lm(sqrt(biomass_Copepod_P2) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed),
                         data=dd_mesocosm)
plot(model_Copepod)
hist(model_Copepod$residuals)
stats::shapiro.test(residuals(model_Copepod))
lmtest::bptest(sqrt(biomass_Copepod_P2) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) , data=dd_mesocosm)

car::Anova(model_Copepod, type=2)
summary(model_Copepod)

## Figure 4C - Rotifers
fig_Rotifer <- ggplot(dd_mesocosm, aes(x=mean_NapieranCoef, y=biomass_Rotifer_P2, color=Mitigation_action_name, fill=Mitigation_action_name)) +
  geom_point()+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  geom_smooth(method="lm", alpha=.1)+
  ylab("Rotifer biomass (µg/L) \n(end of mitigation phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "bottom") ; fig_Rotifer


model_Rotifer <- lm(log(biomass_Rotifer_P2) ~ mean_NapieranCoef * as.factor(Nutrient_stop)* as.factor(Fish_removed),
                         data=dd_mesocosm)
plot(model_Rotifer)
hist(model_Rotifer$residuals)
stats::shapiro.test(residuals(model_Rotifer))
lmtest::bptest(log(biomass_Rotifer_P2) ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) , data=dd_mesocosm)

car::Anova(model_Rotifer, type=2)
summary(model_Rotifer)

## Cladoceran -- not enough data to run an analysis
ggplot(dd_mesocosm, aes(x=mean_NapieranCoef, y=biomass_Cladoceran_P2, color=Mitigation_action_name, fill=Mitigation_action_name)) +
  geom_point()+
  scale_color_manual(values=color_mitigation, name="Mitigation actions:")+
  scale_fill_manual(values=color_mitigation, name="Mitigation actions:")+
  theme_classic()+
  geom_smooth(method="lm", alpha=.1)+
  ylab("Cladoceran biomass (µg/L) \n(end of mitigation phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "bottom") 


#### Figure 4D - Zooplankton composition at the end of the Mitigation phase (P2)

spp.Keratella <- c("K.quadrata", "K.cruciformis", "K.cochlearis", "K.elongated")
taxa.ordered <- c("Bosmina", "Chydorus", "Polyphemus", "polyartha",
                  "Nauplii", "Cyclopoid", "Calanoid", 
                  "Keratella", "Brachionus","Euchlanis", "Asplanchna", "Notholca", "unidentified", "unidentified_2")
taxa.col <- c("#6a44bd", "#b5a2de",
              "#1198bd", "#60bad3","#b0dde9", 
              "#BD0000","#c82b2b","#d35555","#de8080","#e9aaaa","#f4d5d5")

## remove the family Keratella (because there is different species identified) and the family Chironomidae (benthic organism).
zp_data <- subset(raw_zooplankton_data, !(Taxa %in% spp.Keratella) & Taxa!="Chironomid")
zp_data$Taxa <- factor(zp_data$Taxa, levels=taxa.ordered, ordered = TRUE)

## order the mesocosm by Napieran coefficient >>> for figure 4C
order_Napieran <- dd_mesocosm[,c("Mesocosm","mean_NapieranCoef")]
order_Napieran$ord_all <- rank(order_Napieran$mean_NapieranCoef)
zp_data2 <- merge(zp_data, order_Napieran, by="Mesocosm")

fig_ZPcomm <- ggplot2::ggplot(subset(zp_data2, Date=="2021-09-02"), 
                              aes(x=as.factor(ord_all), y=Sample_biomass_µg.L, fill=Taxa)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=taxa.col)+
  facet_grid(Mitigation_action_name~.,scales = "free_x", drop=F, labeller = label_wrap_gen(multi_line = TRUE))+
  ylab("Zooplankton biomass per taxa (µg/L)")+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.grid = element_blank()) ;fig_ZPcomm


### permutation test (adonis) CONTINUOUS

## keep only the mesocosm ID, the group, the taxa and the biomass
## look at the last date of the experiment "2021-09-02"
ZP_P2 <- subset(zp_data, Date=="2021-09-02")[,c("Mesocosm","Taxa","Sample_biomass_µg.L")]; names(ZP_P2) <- c("Mesocosm","Taxa","Biomass")

## transform the data set to have the species as a column
t_ZP_P2 <- stats::reshape(ZP_P2, idvar = c("Mesocosm"), v.names = "Biomass", sep="_", timevar = "Taxa", direction = "wide", new.row.names = unique(ZP_P2$Mesocosm))

# # remove information on mesocosm
t_sppZP_P2 <- t_ZP_P2[,-1]
row.names(t_sppZP_P2)

## add zero to missing values
t_sppZP_P2 <- replace(t_sppZP_P2, is.na(t_sppZP_P2), 0)

# calculate the dissimilarity (Bray-Curtis distance)
sppZP_P2.dist <- vegan::vegdist(t_sppZP_P2, method="bray")

vegan::adonis2(sppZP_P2.dist ~ mean_NapieranCoef*as.factor(Nutrient_stop)*as.factor(Fish_removed), data=dd_mesocosm[order(dd_mesocosm$Mesocosm),], 
               permutations = 9999, add=TRUE, by="terms")


######## FIGURE 4 - 

fig_ZPbiom_2 <- fig_ZPbiom + guides(color=guide_legend(nrow=4, byrow=TRUE))
fig_Copepod_2 <- fig_Copepod + guides(color=guide_legend(nrow=4, byrow=TRUE))
fig_Rotifer_2 <- fig_Rotifer + guides(color=guide_legend(nrow=4, byrow=TRUE))

ggpubr::ggarrange(ggarrange(fig_ZPbiom_2, fig_Copepod_2, fig_Rotifer_2, ncol=1, nrow=3, heights=c(2,1.5,1.5),
                    labels=c("A","B","C"),
                    common.legend = T, legend = "bottom"),
          fig_ZPcomm, labels=c("","D"),
          nrow = 1, ncol=2, widths=c(3,5))
ggplot2::ggsave(paste0(Project_figure_folder,"fig_ZPcomm.png"), units="cm", width=25, height=25, dpi=300)
ggplot2::ggsave(paste0(Project_figure_folder,"fig_ZPcomm.tiff"), units="cm", width=25, height=25, dpi=300)


######### FIGURE 5

str(t_ZP_P2)
dd_ratio <- data.frame(Mesocosm = t_ZP_P2$Mesocosm,
                       diff_ratio_CycloKerat = (t_ZP_P2$Biomass_Cyclopoid-t_ZP_P2$Biomass_Keratella)/(t_ZP_P2$Biomass_Cyclopoid+t_ZP_P2$Biomass_Keratella))
dd_mesocosm <- merge(dd_mesocosm, dd_ratio, by="Mesocosm")

ggplot(dd_mesocosm, 
       aes(x=LRR_chl_pel, y=Total_ZP_biomassP2, color=Mitigation_action_name)) + 
  geom_smooth(method="lm", alpha=.1, se=F)+
  geom_point(aes(fill=diff_ratio_CycloKerat), shape=21,size=4, stroke=1.5)+
  scale_color_manual(name = "Mitigation actions:", values=color_mitigation)+
  scale_fill_stepsn(colors=c("#BD0000","white","#60bad3"), name="Zooplankton dominance: ", breaks=c(-.5, .5), labels=c("<-0.5 Rotifers (Keratella)",">0.5 Copepods (Cyclopoid)"))+
  theme_classic()+
  labs(y="Final (end phase 2) \ntotal zoopplankton biomass (µg/L)", x="LRR Chl-a pelagic")+
  theme(legend.position = "right")
ggplot2::ggsave(paste0(Project_figure_folder,"Figure5_ZP_Chl.png"), units="cm", width=20, height=15, dpi=300)


################################################### 
###################################################  End of Eutrophication phase
################################################### 

######### Figure S6 - Pelagic chl-a, zooplankton biomass and composition at the end of the Eutrophication phase (P1)
### S6a) Pelagic chl-a ~ Napieran coef
fig_pelP1 <- ggplot(dd_mesocosm,
                    aes(x=mean_NapieranCoef, y=conc_Pelagic_chla_P1)) + 
  geom_point()+ geom_smooth(method="lm", color="black")+
  labs(y="Pelagic chl-a (µg/L) \n(end eutrophication phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_bw() ; fig_pelP1

m3 <- lm(conc_Pelagic_chla_P1 ~ mean_NapieranCoef, data=dd_mesocosm)
stats::shapiro.test(residuals(m3)) ; summary(m3)

### S6b) Zooplankton biomass ~ Pelagic chl-a and Napieran coef
fig_ZP_pelP1 <- ggplot(dd_mesocosm, 
                       aes(x=conc_Pelagic_chla_P1, y=Total_ZP_biomassP1, fill=mean_NapieranCoef)) + 
  geom_smooth(method="lm", se=T, color="black")+
  geom_point(shape=21, size=3)+
  scale_fill_gradient(low = "grey", high = "coral4",
                      name = expression(Absorption~coefficient~(m^{-1})))+
  theme_classic()+
  labs(y="Total zoopplankton biomass (µg/L) \n(end eutrophication phase)", x="Chl-a pelagic (µg/L) \n(end eutrophication phase)")+
  theme(legend.position = "bottom") ; fig_ZP_pelP1

### S6c) Zooplankton composition at the end of Eutrophication phase (P1)
fig_ZP_P1 <- ggplot2::ggplot(subset(zp_data2, Date=="2021-08-09"), aes(x=ord_all,y=Sample_biomass_µg.L,fill=Taxa)) +
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=taxa.col)+
  ylab("Zooplankton biomass per taxa (µg/L) -- end Eutrophication phase")+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) ; fig_ZP_P1


### overall figure S6
ggpubr::ggarrange(ggarrange(fig_pelP1, fig_ZP_pelP1, ncol=2, labels="AUTO"), 
          fig_ZP_P1, nrow=2, labels=c("","C"), heights = c(2,3))
ggplot2::ggsave(paste0(Project_figure_folder,"FigureS6_initial_ZP_chla.png"), units="cm", width=30, height=30, dpi=300)


################################################### 
###################################################  Chironomid
################################################### 

raw_chironomid <- read.csv(file=paste0(Project_data_folder,"raw_chironomids.csv"), sep=";")
dd_mesocosm <- merge(dd_mesocosm, raw_chironomid, by="Mesocosm")


fig_chiroNapieran = ggplot(dd_mesocosm, 
                           aes(x=mean_NapieranCoef, y=Nb_chironomids, color=Mitigation_action_name)) + 
  geom_smooth(method="lm", alpha=.1,se=F)+
  geom_point(shape=19,size=2)+
  scale_color_manual(name = "Mitigation actions:", values=color_mitigation)+
  theme_classic()+
  labs(y="Number of chironomids \n (end mitigation phase)")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme(legend.position = "right") ; fig_chiroNapieran


### Negative binomial
glm_chiro = MASS::glm.nb(Nb_chironomids ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed), data = dd_mesocosm)
summary(glm_chiro)

#overdispersion
40.677/28

1 - stats::pchisq(glm_chiro$deviance, glm_chiro$df.residual) # if not significant, we don't need to worry with overdispersion
#>>> seems to be OK

# model selection based on this distribution
m1 = MASS::glm.nb(Nb_chironomids ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed), data = dd_mesocosm)
m2 = MASS::glm.nb(Nb_chironomids ~ mean_NapieranCoef * as.factor(Nutrient_stop) * as.factor(Fish_removed) - as.factor(Nutrient_stop):as.factor(Fish_removed), data = dd_mesocosm)
m3 = MASS::glm.nb(Nb_chironomids ~ mean_NapieranCoef + as.factor(Nutrient_stop) + as.factor(Fish_removed), data = dd_mesocosm)

stats::AIC(m1, m2, m3)
## best model with m3 without interaction between the explanatory variables

GLMnb.chiro = MASS::glm.nb(Nb_chironomids ~ mean_NapieranCoef + as.factor(Nutrient_stop) + as.factor(Fish_removed), data = dd_mesocosm)
summary(GLMnb.chiro) 
# overdispersion
40.593/32
1 - stats::pchisq(GLMnb.chiro$deviance, GLMnb.chiro$df.residual) # if not significant, we don't need to worry with overdispersion
#>>> seems to be OK


##### Link between benthic chl-a increase/decrease (i.e., LRR benthic chl-a) and Nb of chironomids

fig_chiroChlaBen = ggplot(dd_mesocosm, 
                          aes(x=LRR_chl_ben, y=Nb_chironomids, color=Mitigation_action_name)) + 
  geom_smooth(method="lm", alpha=.1, se=F)+
  geom_point(shape=19,size=2)+
  scale_color_manual(name = "Mitigation actions:", values=color_mitigation)+
  #scale_fill_stepsn(colors=c("#BD0000","white","#60bad3"), name="Zooplankton dominance: ", breaks=c(-.5, .5))+
  theme_classic()+
  labs(y="", x="LRR Chl-a benthic")+
  theme(legend.position = "right") ; fig_chiroChlaBen

ggpubr::ggarrange(fig_chiroNapieran, fig_chiroChlaBen, nrow=1, ncol = 2, 
          labels = "AUTO", align="hv",
          common.legend = T, legend = "right")
ggplot2::ggsave(paste0(Project_figure_folder,"FigureS7_chiro_Chl-ben.png"), units="cm", width=22, height=10, dpi=300)

## model with Mitigation actions
model_LRR_ChlBen_chiro <- lm(LRR_chl_ben ~ Nb_chironomids * as.factor(Nutrient_stop) * as.factor(Fish_removed), data=dd_mesocosm)

plot(model_LRR_ChlBen_chiro)
hist(model_LRR_ChlBen_chiro$residuals)
stats::shapiro.test(residuals(model_LRR_ChlBen_chiro))
lmtest::bptest(LRR_chl_ben ~ Nb_chironomids * as.factor(Nutrient_stop) * as.factor(Fish_removed) , data=dd_mesocosm)

car::Anova(model_LRR_ChlBen_chiro, type=2)
summary(model_LRR_ChlBen_chiro)

## Model without Mitigation actions
model_LRR_ChlBen_chiro2 <- lm(LRR_chl_ben ~ Nb_chironomids, data=dd_mesocosm)

plot(model_LRR_ChlBen_chiro2)
hist(model_LRR_ChlBen_chiro2$residuals)
stats::shapiro.test(residuals(model_LRR_ChlBen_chiro2))
lmtest::bptest(LRR_chl_ben ~ Nb_chironomids , data=dd_mesocosm)

car::Anova(model_LRR_ChlBen_chiro2, type=2)
summary(model_LRR_ChlBen_chiro2)

## correlation
stats::cor.test(dd_mesocosm$LRR_chl_ben, dd_mesocosm$Nb_chironomids)


################################################### 
###################################################  Fish (length, weight, Fulton index)
################################################### 

## Measurements of fish initial population
fish_population <- read.csv(file=paste0(Project_data_folder,"raw_fish_population_measurement.csv"))
fish_population$Date <- as.Date(fish_population$Date, format = "%d/%m/%Y")
fish_population$Fulton_index <- 100 * (fish_population$wet_weight_g/((fish_population$total_length_mm/10)^3))


## Fish growth during the experiment
fish_tank <- read.csv(file=paste0(Project_data_folder,"raw_fish_tank_measurement.csv"))

fish_tank$Date <- as.Date(fish_tank$Date, format="%d/%m/%Y")
fish_tank2 <- merge(expt_design2021, fish_tank, by="Mesocosm")
fish_tank2 <- merge(fish_tank2, dd_mesocosm[,c("Mesocosm","mean_NapieranCoef")])

fish_tank2$Fulton_index <- 100 * (fish_tank2$wet_weight_g/((fish_tank2$total_length_mm/10)^3))

## FIGURE S8 Fish measurements

length_startP1 <- ggplot(subset(fish_population, Date=="2021-07-20" & Size_selection=="yes"), aes(y=standard_length_mm)) +
  ylim(8,35)+ 
  labs(y="standard length (mm)")+
  geom_boxplot() + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()); length_startP1

length_endP1 <- ggplot(subset(fish_tank2, PhaseID_missing=="phase1"), aes(x=mean_NapieranCoef, y=standard_length_mm)) + 
  geom_point() + 
  ylim(8,35)+ 
  geom_smooth(method="lm", color="black")+
  labs(y="", title="Eutrophication phase")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_classic() ; length_endP1

length_startP2 <- ggplot(subset(fish_population, Date=="2021-08-11" & Size_selection=="yes"), aes(y=standard_length_mm)) +
  labs(y="")+
  ylim(8,35)+ 
  geom_boxplot()+ 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()); length_startP2

length_endP2 <- ggplot2::ggplot(subset(fish_tank2, PhaseID=="phase2"), aes(x=mean_NapieranCoef, y=standard_length_mm, color=Mitigation_action_name, fill=Mitigation_action_name)) + 
  scale_color_manual(values=c(color_mitigation[1],color_mitigation[3]))+
  scale_fill_manual(values=c(color_mitigation[1],color_mitigation[3]))+
  geom_smooth(method="lm")+
  ylim(8,35)+ 
  labs(y="", title="Mitigation phase")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  geom_point() + 
  geom_point(data=subset(fish_tank2, PhaseID_missing=="phase1_bis"),aes(x=mean_NapieranCoef, y=standard_length_mm), color="red", position = position_jitter(width=0.05, height=0), shape=8, show.legend = FALSE) +
  theme_classic() ; length_endP2

### Weight
weight_startP1 <- ggplot(subset(fish_population, Date=="2021-07-20" & Size_selection=="yes"), aes(y=wet_weight_g)) +
  ylim(0.01,0.4)+
  labs(y="wet weight (g)")+
  geom_boxplot() + 
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()); weight_startP1

weight_endP1 <- ggplot(subset(fish_tank2, PhaseID_missing=="phase1"), aes(x=mean_NapieranCoef, y=wet_weight_g)) + 
  geom_point() +
  ylim(0.01,0.4)+
  geom_smooth(method="lm", color="black")+labs(y="")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_classic() ; weight_endP1

weight_startP2 <- ggplot(subset(fish_population, Date=="2021-08-11" & Size_selection=="yes"), aes(y=wet_weight_g)) +
  labs(y="")+
  ylim(0.01,0.4)+
  geom_boxplot()+ 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()); weight_startP2

weight_endP2 <- ggplot2::ggplot(subset(fish_tank2, PhaseID=="phase2"), aes(x=mean_NapieranCoef, y=wet_weight_g, color=Mitigation_action_name, fill=Mitigation_action_name)) + 
  scale_color_manual(values=c(color_mitigation[1],color_mitigation[3]))+
  scale_fill_manual(values=c(color_mitigation[1],color_mitigation[3]))+
  geom_smooth(method="lm")+
  ylim(0.01,0.4)+ 
  labs(y="")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  geom_point() + theme_classic() + 
  geom_point(data=subset(fish_tank2, PhaseID_missing=="phase1_bis"),aes(x=mean_NapieranCoef, y=wet_weight_g), color="red", position = position_jitter(width=0.05, height=0), shape=8, show.legend = FALSE) ; weight_endP2


### Fulton index
Fulton_startP1 <- ggplot(subset(fish_population, Date=="2021-07-20" & Size_selection=="yes"), aes(y=Fulton_index)) +
  ylim(0.4,1.2)+
  labs(y="Fulton index")+
  geom_boxplot() + 
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()); Fulton_startP1

Fulton_endP1 <- ggplot(subset(fish_tank2, PhaseID_missing=="phase1"), aes(x=mean_NapieranCoef, y=Fulton_index)) + 
  geom_point() +
  ylim(0.4,1.2) +
  geom_smooth(method="lm", color="black")+
  labs(y="")+
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_classic() ; Fulton_endP1

Fulton_startP2 <- ggplot(subset(fish_population, Date=="2021-08-11" & Size_selection=="yes"), aes(y=Fulton_index)) +labs(y="")+
  ylim(0.4,1.2)+
  geom_boxplot()+ 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()); Fulton_startP2

Fulton_endP2 <- ggplot2::ggplot(subset(fish_tank2, PhaseID=="phase2"), aes(x=mean_NapieranCoef, y=Fulton_index, color=Mitigation_action_name, fill=Mitigation_action_name)) + 
  geom_point() + 
  scale_color_manual(values=c(color_mitigation[1],color_mitigation[3]))+
  scale_fill_manual(values=c(color_mitigation[1],color_mitigation[3]))+
  geom_smooth(method="lm")+
  ylim(0.4,1.2)+ 
  labs(y="")+ 
  geom_point(data=subset(fish_tank2, PhaseID_missing=="phase1_bis"),aes(x=mean_NapieranCoef, y=Fulton_index), color="red", position = position_jitter(width=0.05, height=0), shape=8, show.legend=FALSE) + 
  xlab(expression(Absorption~coefficient~(m^{-1})))+
  theme_classic() ; Fulton_endP2



## FIGURE S8 - fish measurements

ggpubr::ggarrange(length_startP1, length_endP1, length_startP2, length_endP2,
          weight_startP1, weight_endP1, weight_startP2, weight_endP2,
          Fulton_startP1, Fulton_endP1, Fulton_startP2, Fulton_endP2,
          nrow=3,ncol=4, widths=c(1,3,1,3), legend = "bottom",common.legend = T,align="hv", labels = "AUTO")

ggplot2::ggsave(paste0(Project_figure_folder,"fig_fish.png"), units="cm", width=20, height=30, dpi=300)








