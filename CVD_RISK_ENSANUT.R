
## 10-Year Fatal Cardiovascular Risk Estimator
## Recalibrated for Mexico (Globorisk MEX 2020, fatal events)

## Authors of original analysis:
##   Neftali Eduardo Antonio Villa
##   Omar Yaxmehen Bello-Chavolla
## 

#### Dataset: Working Directory ####

setwd("/Users/nefoantonio/Library/CloudStorage/GoogleDrive-nantonio@facmed.unam.mx/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/0_Finalized projects/CVD risk scores/Table_Risks")

#### Dataset: Library ####

library(survival);library(haven);library(survey);library("jtools")
library("remotes");library("readr");library("dplyr");library("xgboost")
library("caret");library(purrr);library(dplyr);library(readr)
library(tidyr);library("na.tools");library("preventr");library("radiant.data")
library("Hmisc");library("gtsummary");library("flextable")
library("treemapify");library(risks)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(officer)
library(globorisk)

#### Dataset: Loading Datasets#####

#ENSANUT MERGED

base<-read_csv("base_CMKS_final.csv")

#### Dataset: Filter Dataset ####

base.2 <- base 

table(base.2$Year)

#Filter 1: Age >40
base.2<-base.2%>%
  dplyr::filter(edad>=40)

table(base.2$Year)
nrow(base.2)

#Filter 2: Complete weighted and Prevent calculation

base.2<-base.2%>%
  dplyr::filter(!is.na(ponde_bioq))%>%
  dplyr::filter(!is.na(PAS))%>%
  dplyr::filter(!is.na(ct))

table(base.2$Year)
nrow(base.2)

#Filter 4: Extreme values
base.2<-base.2%>%
  dplyr::mutate(PAS.meanval=mean(PAS,na.rm=T),PAS.stdev=sd(PAS,na.rm=T),
                ct.meanval=mean(ct,na.rm=T),ct.stdev=sd(ct,na.rm = T))


base.2<-base.2%>%
  dplyr::filter(abs((PAS-PAS.meanval)/PAS.stdev)<4)%>%
  dplyr::filter(abs((ct-ct.meanval)/ct.stdev)<4)

table(base.2$Year)
nrow(base.2)


#### Dataset: Recoding Dataset ####

#Region
base.2$region<-NULL
base.2$region[base.2$entidad %in% c(2,3,5,8,10,19,26,28)]<-1
base.2$region[base.2$entidad %in% c(1,6,11,14,16,17,18,22,24,25,32)]<-2
base.2$region[base.2$entidad %in% c(9,15)]<-3
base.2$region[base.2$entidad %in% c(4,7,12,13,20,21,29,23,27,30,31)]<-4
base.2$region<-factor(base.2$region,labels = c("North-Region","Central-Region","Metropolitan-Area", "South-Region"))

#Subregion
base.2$SUBREGION[base.2$entidad %in% c(2,3,18,25,26)]<-1
base.2$SUBREGION[base.2$entidad %in% c(5,8,19,28)]<-2
base.2$SUBREGION[base.2$entidad %in% c(6,14,16)]<-3
base.2$SUBREGION[base.2$entidad %in% c(1,10,11,22,24,32)]<-4
base.2$SUBREGION[base.2$entidad %in% c(13,29,30)]<-5
base.2$SUBREGION[base.2$entidad %in% c(9)]<-6
base.2$SUBREGION[base.2$entidad %in% c(15)]<-7
base.2$SUBREGION[base.2$entidad %in% c(12,17,20,21)]<-8
base.2$SUBREGION[base.2$entidad %in% c(4,7,27,23,31)]<-9

#Age groups
base.2$age_group <- NULL
base.2$age_group [base.2$edad<45] <- 1 # 30 - 44 years
base.2$age_group [base.2$edad>=45 & base.2$edad<65] <- 2 # 44 - 64 years 
base.2$age_group [base.2$edad>=65] <- 3 # 65 - 79 years

# Categorizacion BMI
base.2$IMC_CAT<-NULL
base.2$IMC_CAT[base.2$IMC<18.5]<-0
base.2$IMC_CAT[base.2$IMC>=18.5 & base.2$IMC<25]<-1
base.2$IMC_CAT[base.2$IMC>=25 & base.2$IMC<30]<-2
base.2$IMC_CAT[base.2$IMC>=30 & base.2$IMC<40]<-3
base.2$IMC_CAT[base.2$IMC>=40]<-4

#Categorizacion BMI2
base.2$IMC_CAT2<-NA
base.2$IMC_CAT2[base.2$IMC<18.5]<-0
base.2$IMC_CAT2[base.2$IMC>=18.5 & base.2$IMC<25]<-1
base.2$IMC_CAT2[base.2$IMC>=25 & base.2$IMC<30]<-2
base.2$IMC_CAT2[base.2$IMC>=30]<-3

base.2$Normopeso <- NULL
base.2$Normopeso[base.2$IMC>=18.5 & base.2$IMC<25] <- 1
base.2$Normopeso[base.2$IMC>25] <- 0

base.2$Sobrepeso <- NULL
base.2$Sobrepeso[base.2$IMC>=25 & base.2$IMC<30] <- 1
base.2$Sobrepeso[base.2$IMC<25] <- 0
base.2$Sobrepeso[base.2$IMC>=30] <- 0

base.2$Obesidad <- NULL
base.2$Obesidad[base.2$IMC>=30] <- 1
base.2$Obesidad[base.2$IMC<30] <- 0

# Metabolic Syndrome

base.2$FEMALE <- NULL
base.2$FEMALE[base.2$sexo==0] <- 1
base.2$FEMALE[base.2$sexo==1] <- 0

base.2$HDL_ALT <- ifelse(
  base.2$FEMALE==1 & base.2$hdl<50,1,
  ifelse(base.2$FEMALE==0&base.2$hdl<40,1,0))

base.2$TG_ALT.SM<- NULL
base.2$TG_ALT.SM [base.2$trig>=150] <- 1
base.2$TG_ALT.SM [base.2$trig<150] <- 0

base.2$CT_ALT <- NULL
base.2$CT_ALT[base.2$ct>=200] <- 1
base.2$CT_ALT[base.2$ct<200] <- 0

base.2$Obesidad_C <- ifelse(
  base.2$FEMALE==1 & base.2$cintura>=88,1,
  ifelse(base.2$FEMALE==0 & base.2$cintura>=102,1,0))

base.2$PA_ALT <- NULL
base.2$PA_ALT[base.2$antihta == 1 |
                base.2$PAS >= 130| base.2$PAD >= 85]<-1
base.2$PA_ALT<-na.tools::na.replace(base.2$PA_ALT,0)

base.2$GLC_ALT <- NULL
base.2$GLC_ALT[base.2$glc>=110|
                 base.2$hipoglc == 1 ] <- 1
base.2$GLC_ALT <- na.tools::na.replace(base.2$GLC_ALT,0)

base.2$Criterios_SM <- with(base.2, Obesidad_C + PA_ALT + GLC_ALT + TG_ALT.SM + HDL_ALT)
base.2$SM <- ifelse(base.2$Criterios_SM>=3,1,0)

base.2 <- base.2 %>%
  mutate(
    WHtR        = cintura / talla)

#Altered tryglicerides for Cardiovascular-kidney-metabolic syndrome

base.2$tg_alt <- NULL
base.2$tg_alt [base.2$trig>=135] <- 1
base.2$tg_alt [base.2$trig<135] <- 0

# Prediabetes
base.2$prediabetes <- NULL
base.2$prediabetes[base.2$hba1c>=5.7 & base.2$hba1c <6.5|
                     base.2$glc>=100 & base.2$glc<126] <- 1
base.2$prediabetes<-na.tools::na.replace(base.2$prediabetes,0)

#Diabetes
base.2$DM<-NULL
base.2$DM[base.2$dxdm==1 | 
            base.2$hipoglc == 1 |
            base.2$glc >=126 |
            base.2$hba1c >= 6.5]<-1
base.2$DM<-na.tools::na.replace(base.2$DM,0) 

# Hipertension
base.2$HAS<-NULL
base.2$HAS[base.2$dxhta==1 | 
             base.2$antihta == 1 |
             base.2$PAS >= 140 | base.2$PAD >= 90]<-1
base.2$HAS<-na.tools::na.replace(base.2$HAS,0)

#ERC
base.2$erc<-na.tools::na.replace(base.2$erc,0)

#Prevent estimation
#smoking dummy
base.2$smoking_dummy [base.2$smoking == 2] <- 1
base.2$smoking_dummy [base.2$smoking == 0] <- 0
base.2$smoking_dummy [base.2$smoking == 1] <- 0
base.2$smoking_dummy<-na.tools::na.replace(base.2$smoking_dummy,0)

#antihta dummy
base.2$antihta_dummy <- base.2$antihta
base.2$antihta_dummy<-na.tools::na.replace(base.2$antihta_dummy,0)

#dummy estatinas
base.2$estatinas_dummy <- base.2$estatinas
base.2$estatinas_dummy<-na.tools::na.replace(base.2$estatinas_dummy,0)

#Dummy Age
base.2$edad_dummy<-NULL
base.2$edad_dummy<-base.2$edad
base.2$edad_dummy[base.2$edad_dummy<=40]<-40
base.2$edad_dummy[base.2$edad_dummy>=80]<-79

#Colesterol Dummy
base.2$colesterol_dummy<- base.2$ct
base.2$colesterol_dummy[base.2$colesterol_dummy<=130]<-130
base.2$colesterol_dummy[base.2$colesterol_dummy>=320]<-320

#Colestero HDL Dummy 2 (Prevent)
base.2$hdl_dummy<-base.2$hdl
base.2$hdl_dummy[base.2$hdl_dummy<=20]<-20
base.2$hdl_dummy[base.2$hdl_dummy>=100]<-100

#TFG Dummy 2 (Prevent)
base.2$eGFR_dummy<-base.2$eGFR
base.2$eGFR_dummy[base.2$eGFR_dummy<=15] <- 15
base.2$eGFR_dummy[base.2$eGFR_dummy>=140] <- 140

#BMI Dummy  (Prevent)
base.2$imc_dummy<-base.2$IMC
base.2$imc_dummy[base.2$imc_dummy<=18.5] <- 18.5
base.2$imc_dummy[base.2$imc_dummy>=39.9]<-39.9

#SBP Dummy  (Prevent)
base.2$PAS_dummy<-base.2$PAS
base.2$PAS_dummy[base.2$PAS_dummy<=90]<-90
base.2$PAS_dummy[base.2$PAS_dummy>=180]<-180

#Gender Dummy
base.2$gender<-NULL
base.2$gender[base.2$sexo==1]<-"male"
base.2$gender[base.2$sexo==0]<-"female"

#Globorisk - Recalibrated Function 

globorisk_recalib_fast <- function(sex, age, sbp, tcol_mmol, diab, smoke, validate = TRUE) {
  
  if (validate) {
    # Verificación mínima y vectorizada
    vars <- list(sex = sex, age = age, sbp = sbp,
                 tcol_mmol = tcol_mmol, diab = diab, smoke = smoke)
    
    no_numericas <- names(vars)[!vapply(vars, is.numeric, logical(1))]
    if (length(no_numericas) > 0) {
      stop(sprintf(
        "Estas variables deben ser numéricas: %s",
        paste(no_numericas, collapse = ", ")
      ))
    }
    
    if (anyNA(sex) || anyNA(age) || anyNA(sbp) || anyNA(tcol_mmol) || anyNA(diab) || anyNA(smoke)) {
      stop("No se permiten valores NA en las variables de entrada.")
    }
    
    if (any(sex != 0 & sex != 1)) {
      stop("sexo debe ser 0 = hombre y 1 = mujer")
    }
    
    if (any(diab != 0 & diab != 1)) {
      stop("diabetes debe ser 0 = no y 1 = sí")
    }
    
    if (any(smoke != 0 & smoke != 1)) {
      stop("tabaquismo debe ser 0 = no y 1 = sí")
    }
    
    n <- length(sex)
    if (!all(lengths(vars) == n)) {
      stop("Todas las variables deben tener la misma longitud.")
    }
  }
  
  # Riesgo base Globorisk
  score <- globorisk(
    sex = sex,
    age = age,
    sbp = sbp,
    tc = tcol_mmol,
    dm = diab,
    smk = smoke,
    iso = "MEX",
    year = 2020,
    version = "fatal",
    type = "risk"
  )
  
  # Parámetros por sexo: 0 = hombre, 1 = mujer
  alpha <- c(-0.8639137, -0.881977)
  beta  <- c( 0.7745622,  1.002298)
  
  surv_model <- 1 - score
  log_term <- log(-log(surv_model))
  
  # sex + 1 sirve para indexar: 0->1 (hombre), 1->2 (mujer)
  log_recalib <- alpha[sex + 1] + beta[sex + 1] * log_term
  surv_recalib <- exp(-exp(log_recalib))
  rsk_recalib <- 1 - surv_recalib
  
  round(rsk_recalib * 100, 0)
}

base.2$sexo_rec<-NULL
base.2$sexo_rec[base.2$sexo==0]<-1
base.2$sexo_rec[base.2$sexo==1]<-0

base.2$ct_mmol<-base.2$ct*0.02586

base.2<-base.2 %>%
  mutate(
    globorisk_recalib = globorisk_recalib_fast(
      sex = sexo_rec,
      age = edad_dummy,
      sbp = PAS,
      tcol_mmol = ct_mmol,
      diab = DM,
      smoke = smoking_dummy))


base.2$globorisk_risk_cat<-NULL
base.2$globorisk_risk_cat[base.2$globorisk_recalib<5]<-1
base.2$globorisk_risk_cat[base.2$globorisk_recalib>=5 & base.2$globorisk_recalib<=10]<-2
base.2$globorisk_risk_cat[base.2$globorisk_recalib>=10 & base.2$globorisk_recalib<=20]<-3
base.2$globorisk_risk_cat[base.2$globorisk_recalib>=20]<-4

#Globorisk Low
base.2$globorisk_low<-NULL
base.2$globorisk_low[base.2$globorisk_risk_cat==1]<-1
base.2$globorisk_low[base.2$globorisk_risk_cat!=1]<-0

#Globorisk Intermediate 
base.2$globorisk_intermediate<-NULL
base.2$globorisk_intermediate[base.2$globorisk_risk_cat==2]<-1
base.2$globorisk_intermediate[base.2$globorisk_risk_cat!=2]<-0

#Globorisk High
base.2$globorisk_high<-NULL
base.2$globorisk_high[base.2$globorisk_risk_cat==3]<-1
base.2$globorisk_high[base.2$globorisk_risk_cat!=3]<-0

#Globorisk Very-High
base.2$globorisk_very_high<-NULL
base.2$globorisk_very_high[base.2$globorisk_risk_cat==4]<-1
base.2$globorisk_very_high[base.2$globorisk_risk_cat!=4]<-0


base.2$globorisk_original<-globorisk(
  sex = base.2$sexo_rec,
  age = base.2$edad_dummy,
  sbp = base.2$PAS,
  tc = base.2$ct_mmol,
  dm = base.2$DM,
  smk = base.2$smoking_dummy,
  iso = rep("MEX", nrow(base.2)),
  year = "2020",
  version = "lab",
  type = "risk"
)*100


base.2$globorisk_risk_original_cat<-NULL
base.2$globorisk_risk_original_cat[base.2$globorisk_original<5]<-1
base.2$globorisk_risk_original_cat[base.2$globorisk_original>=5 & base.2$globorisk_original<=10]<-2
base.2$globorisk_risk_original_cat[base.2$globorisk_original>=10 & base.2$globorisk_original<=20]<-3
base.2$globorisk_risk_original_cat[base.2$globorisk_original>=20]<-4

#globorisk_original_original Low
base.2$globorisk_original_low<-NULL
base.2$globorisk_original_low[base.2$globorisk_risk_original_cat==1]<-1
base.2$globorisk_original_low[base.2$globorisk_risk_original_cat!=1]<-0

#globorisk_original Intermediate 
base.2$globorisk_original_intermediate<-NULL
base.2$globorisk_original_intermediate[base.2$globorisk_risk_original_cat==2]<-1
base.2$globorisk_original_intermediate[base.2$globorisk_risk_original_cat!=2]<-0

#globorisk_original High
base.2$globorisk_original_high<-NULL
base.2$globorisk_original_high[base.2$globorisk_risk_original_cat==3]<-1
base.2$globorisk_original_high[base.2$globorisk_risk_original_cat!=3]<-0

#globorisk_original Very-High
base.2$globorisk_original_very_high<-NULL
base.2$globorisk_original_very_high[base.2$globorisk_risk_original_cat==4]<-1
base.2$globorisk_original_very_high[base.2$globorisk_risk_original_cat!=4]<-0

#### Dataset: Survey Objects - ENSANUT 2024####

options(survey.lonely.psu = "adjust")
options(survey.adjust.domain.lonely = TRUE)

base_2024 <- base.2%>%filter(Year==2024)

base_2024 <- base_2024 %>% mutate(
  "pond.ant"=ponde_bioq, "estr.ant"=varianza,
  "sampleid"=paste(entidad, as.numeric(estr.ant))) %>%
  filter(!is.na(pond.ant),!is.na(estr.ant))

wts <- with(base_2024, tapply(pond.ant, sampleid, sum)) %>%
  data.frame() %>% `names<-`("wt") %>% rownames_to_column("ids")

t1<-as.data.frame(table(base_2024$sampleid))

wts2<-data.frame(
  ids=wts$ids, sumwt=wts$wt, jn=t1$Freq)
base_2024 <- merge(base_2024, wts2, by.x="sampleid", by.y="ids", all.x=T)
base_2024$swts <- base_2024$pond.ant*(base_2024$jn/base_2024$sumwt)

ensanut_2024_survey <- svydesign(data=base_2024, id=~id_persona, strata=~varianza, weights=~ponde_bioq, nest=TRUE) #SD1

### All prevalence ###

prev.2024.all.low<-svyby(~globorisk_low, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_low*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Low") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

prev.2024.all.intermediate<-svyby(~globorisk_intermediate, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_intermediate*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Intermediate") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

prev.2024.all.high<-svyby(~globorisk_high, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_high*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="High") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

prev.2024.all.very_high<-svyby(~globorisk_very_high, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_very_high*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Very-High") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

dplyr::bind_rows(
  prev.2024.all.low,
  prev.2024.all.intermediate,
  prev.2024.all.high,
  prev.2024.all.very_high) %>%
  mutate(
    Year = as.integer(Year),
    cluster = factor(
      cluster,
      levels = c("Low", "Intermediate", "High", "Very-High"),
      ordered = TRUE)) %>%
  mutate(
    label = sprintf("%.1f (%.1f–%.1f)", prop, lIC95, uIC95)) %>%
  ggplot(aes(area = prop, fill = cluster, label = label)) +
  geom_treemap(layout="squarified") +
  geom_treemap_text(colour="black", place="centre", size=16) +
  theme_void() +
  labs(fill="CVD Risk",
       title = "CVD Risk in Mexico in 2024")+
  scale_fill_manual(
    values = c(
      "Low" = "#4CAF50",
      "Intermediate" = "#FFEB3B",
      "High" = "#FF9800",
      "Very-High" =  "#F44336"))+
  theme_nice()+
  xlab("")+
  ylab("")





prev.2024.all.low_og<-svyby(~globorisk_original_low, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_original_low*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Low") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

prev.2024.all.intermediate_og<-svyby(~globorisk_original_intermediate, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_original_intermediate*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Intermediate") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

prev.2024.all.high_og<-svyby(~globorisk_original_high, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_original_high*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="High") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

prev.2024.all.very_high_og<-svyby(~globorisk_original_very_high, by=~Year, design=ensanut_2024_survey, svymean, na.rm=T) %>%
  mutate(prop=round(globorisk_original_very_high*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="Very-High") %>%
  dplyr::select(Year, cluster, prop, lIC95, uIC95)

dplyr::bind_rows(
  prev.2024.all.low_og,
  prev.2024.all.intermediate_og,
  prev.2024.all.high_og,
  prev.2024.all.very_high_og)

