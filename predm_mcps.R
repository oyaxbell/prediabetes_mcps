# Prediabetes and risk of all-cause and cause-specific mortality: a prospective study of 114,062 adults in Mexico City
# Data Analysis: Omar Yaxmehen Bello-Chavolla & Carlos A. Fermin-Martinez
# Latest version of Analysis February 2025
# For any question regarding analysis contact:
# Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#### Library------------------------ ####
pacman::p_load(
  haven, tidyverse, data.table, ggpubr, survival, flexsurv, coxme, coxphw,
  flextable, gtsummary, officer, pROC, timeROC, survivalROC, lmtest, nortest,
  jtools, viridis, forestplot, forestmodel, forestploterm, PropCIs, Epi, popEpi,
  caret, glmnet, bestNormalize, blandr, BlandAltmanLeh, corrplot, rms, dummy,
  ggalluvial, ggsankey, ggbreak, ggedit, ggimage, ggsci, ggplotify, gtools,
  survminer, magrittr, foreach, doParallel, patchwork, lubridate, naniar,
  circlize, shadowtext, khroma, FactoMineR, fpc, NbClust, cluster, UpSetR,
  TimeVTree, simPH, factoextra, gridExtra,  nhanesA, wesanderson, chorddiag)

#### Data set management------------ ####
#Official working directories:
setwd("~/Mi unidad (obello@facmed.unam.mx)/Mexico City Cohort Study")
#setwd("/Users/carlosfermin/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/Mexico City Cohort Study")

### Load dataset ###
mcps<-read.csv("Data/2021-004 MCPS BASELINE.csv") %>%
  left_join(read.csv("Data/2022-012 MCPS MORTALITY.csv"), by="PATID") %>%
  left_join(read.csv("Data/2022-012 MCPS RESURVEY.csv"), by="PATID") %>% 
  left_join(read.csv("Data/2022-012_OmarBelloChavolla_NMR_Second_Release-Corrected/2022-012 MCPS BASE_NMR_DATA_RECALIB-corrected.csv"), by="PATID")

### Recode dataset ###
#Diabetes
mcps$diabetes <- with(mcps, ifelse((
  BASE_HBA1C>=6.5|BASE_DIABETES==1|DRUG_D1==1|DRUG_D2==1|DRUG_D3==1|DRUG_D4==1), 1, 0))
#Prediabetes
mcps$predm<- ifelse(mcps$BASE_HBA1C>=5.7, 1,0) #ADA
mcps$predm_ice<- ifelse(mcps$BASE_HBA1C>=6.0,1,0) #IEC
mcps$predm[mcps$diabetes==1] <- 0; mcps$predm_ice[mcps$diabetes==1] <- 0
mcps$pred_ADA <- mcps$predm; mcps$pred_IEC <- mcps$predm_ice
#HbA1C categories
mcps$hba1c_categories <-  mcps$BASE_HBA1C %>% cut(c(-Inf, 5.7, 6, Inf), right = F) %>%
  factor(labels = c("<5.7%", "5.7-5.9%", "6.0-6.4%")); mcps$hba1c_cat<- mcps$hba1c_categories
#Follow-up
mcps <- mcps %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps$age_risk<-mcps$AGE+mcps$PERSON_YEARS

#Deaths
#All cause
mcps$D000_p<-ifelse(mcps$D000==1 & mcps$age_risk<75,1,0); mcps$D000_o<-ifelse(mcps$D000==1 & mcps$age_risk>=75,1,0)
mcps$D000_p5<-ifelse(mcps$D000_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D000_o5<-ifelse(mcps$D000_o==1 & mcps$PERSON_YEARS>5,1,0)
#Cardiac
mcps$D003_p<-ifelse(mcps$D003==1 & mcps$age_risk<75,1,0); mcps$D003_o<-ifelse(mcps$D003==1 & mcps$age_risk>=75,1,0)
mcps$D003_p5<-ifelse(mcps$D003_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D003_o5<-ifelse(mcps$D003_o==1 & mcps$PERSON_YEARS>5,1,0)
#Cerebrovascular
mcps$D008_p<-ifelse(mcps$D008==1 & mcps$age_risk<75,1,0); mcps$D008_o<-ifelse(mcps$D008==1 & mcps$age_risk>=75,1,0)
mcps$D008_p5<-ifelse(mcps$D008_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D008_o5<-ifelse(mcps$D008_o==1 & mcps$PERSON_YEARS>5,1,0)
#Other vascular
mcps$D014_p<-ifelse(mcps$D014==1 & mcps$age_risk<75,1,0); mcps$D014_o<-ifelse(mcps$D014==1 & mcps$age_risk>=75,1,0)
mcps$D014_p5<-ifelse(mcps$D014_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D014_o5<-ifelse(mcps$D014_o==1 & mcps$PERSON_YEARS>5,1,0)
#Renal
mcps$D019_p<-ifelse(mcps$D019==1 & mcps$age_risk<75,1,0); mcps$D019_o<-ifelse(mcps$D019==1 & mcps$age_risk>=75,1,0)
mcps$D019_p5<-ifelse(mcps$D019_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D019_o5<-ifelse(mcps$D019_o==1 & mcps$PERSON_YEARS>5,1,0)
#Acute diabetic
mcps$D023_p<-ifelse(mcps$D023==1 & mcps$age_risk<75,1,0); mcps$D023_o<-ifelse(mcps$D023==1 & mcps$age_risk>=75,1,0)
mcps$D023_p5<-ifelse(mcps$D023_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D023_o5<-ifelse(mcps$D023_o==1 & mcps$PERSON_YEARS>5,1,0)
#Hepatobiliary
mcps$D022_p<-ifelse(mcps$D022==1 & mcps$age_risk<75,1,0); mcps$D022_o<-ifelse(mcps$D022==1 & mcps$age_risk>=75,1,0)
mcps$D022_p5<-ifelse(mcps$D022_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D022_o5<-ifelse(mcps$D022_o==1 & mcps$PERSON_YEARS>5,1,0)
#Other GI
mcps$D045_p<-ifelse(mcps$D045==1 & mcps$age_risk<75,1,0); mcps$D045_o<-ifelse(mcps$D045==1 & mcps$age_risk>=75,1,0)
mcps$D045_p5<-ifelse(mcps$D045_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D045_o5<-ifelse(mcps$D045_o==1 & mcps$PERSON_YEARS>5,1,0)
#Neoplastic
mcps$D044_p<-ifelse(mcps$D044==1 & mcps$age_risk<75,1,0); mcps$D044_o<-ifelse(mcps$D044==1 & mcps$age_risk>=75,1,0)
mcps$D044_p5<-ifelse(mcps$D044_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D044_o5<-ifelse(mcps$D044_o==1 & mcps$PERSON_YEARS>5,1,0)
#Respiratory
mcps$D040_p<-ifelse(mcps$D040==1 & mcps$age_risk<75,1,0); mcps$D040_o<-ifelse(mcps$D040==1 & mcps$age_risk>=75,1,0)
mcps$D040_p5<-ifelse(mcps$D040_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$D040_o5<-ifelse(mcps$D040_o==1 & mcps$PERSON_YEARS>5,1,0)
#External, ill-defined, other
mcps$ext_ill <- with(mcps, ifelse((D000+D003+D008+D014+D019+D023+D022+D045+D044+D040)==1, 1, 0))
mcps$ext_ill_p<-ifelse(mcps$ext_ill==1 & mcps$age_risk<75,1,0); mcps$ext_ill_o<-ifelse(mcps$ext_ill==1 & mcps$age_risk>=75,1,0)
mcps$ext_ill_p5<-ifelse(mcps$ext_ill_p==1 & mcps$PERSON_YEARS>5,1,0); mcps$ext_ill_o5<-ifelse(mcps$ext_ill_o==1 & mcps$PERSON_YEARS>5,1,0)

#Sex
mcps$Sex<-ifelse( mcps$MALE==1, 0, 1)
mcps$female<-ifelse( mcps$MALE==0, 1, 0)
#Age
mcps$edad60<-ifelse( mcps$AGE>=60, 1, 0)
mcps$edad60<-factor( mcps$edad60, labels = c("<60 years", "≥60 years"))

#Alcohol Intake
mcps$alcohol<- mcps$ALCGP;  mcps$alcohol[ mcps$ALCGP %in% c(3:5)]<-3
mcps$alcohol<-factor( mcps$alcohol,labels = c("Never","Former","Current"))
#Smoking
mcps$smoke<- mcps$SMOKEGP;  mcps$smoke[ mcps$SMOKEGP %in% c(3:5)]<-3
mcps$smoke<-factor( mcps$smoke,labels = c("Never", "Former","Current"))
#Physical Activity
mcps$PHYSGP<-factor( mcps$PHYSGP,labels = c("None", "≤2 times a week","≥3 times a week"))
mcps$phys<-ifelse(mcps$PHYSGP=="None",0,1)

#Others
mcps$BMI<- mcps$WEIGHT/(( mcps$HEIGHT/100)^2)
mcps$ICE<- mcps$WAISTC/( mcps$HEIGHT)
mcps$R_HXDIAB[is.na( mcps$R_HXDIAB)]<-0
mcps$ASCVD<-ifelse( mcps$D003==1 |  mcps$D008==1, 1, 0)
mcps$college<-ifelse(mcps$EDU_LEVEL==1,1,0)

#Labels
setattr( mcps$AGE, "label", "Age, (Years)")
setattr( mcps$female, "label", "Female sex, (%)")
setattr( mcps$COYOACAN, "label", "Residence in Coyoacan, (%)")
setattr( mcps$EDU_LEVEL, "label", "Educationnal Attainments, (%)")
setattr( mcps$INCOME, "label", "Income, (pesos/month)")
setattr( mcps$smoke, "label", "Smoking, (%)")
setattr( mcps$alcohol, "label", "Alcohol intake, (%)")
setattr( mcps$PHYSGP, "label", "Leisure-time physical activity, (%)")
setattr( mcps$phys, "label", "Regular leisure-time physical activity, (%)")
setattr( mcps$BASE_DIABETES, "label", "Diabetes, (%)")
setattr( mcps$BASE_HYPERTENSION, "label", "Hypertension, (%)")
setattr( mcps$BASE_STROKE, "label", "Stroke, (%)")
setattr( mcps$BASE_CKD, "label", "CKD, (%)")
setattr( mcps$BASE_HEARTATTACK, "label", "Myocardial infarction, (%)")
setattr( mcps$BMI, "label", "Body Mass Index, (kg/m2)")
setattr( mcps$ICE, "label", "Waist-to-Height Ratio")
setattr( mcps$HIPC, "label", "Hip Circunference, (cm)")
setattr( mcps$SBP1, "label", "Systolic Blood Pressure, (mmHg)")
setattr( mcps$DBP1, "label", "Diastolic Blood Pressure, (mmHg)")
setattr( mcps$BASE_HBA1C, "label", "HbA1c, (%)")
setattr( mcps$predm, "label", "HbA1c ≥5.7%")
setattr( mcps$predm_ice, "label", "HbA1c ≥6.0%")
setattr( mcps$edad60, "label", "Age ≥60")
setattr( mcps$college, "label", "University/college educated (%)")
setattr( mcps$Clinical_LDL_C, "label", "LDL cholesterol, mmol/L")
setattr( mcps$HDL_C, "label", "HDL cholesterol, mmol/L")
setattr( mcps$Total_TG, "label", "Triglycerides, mmol/L")
setattr( mcps$ApoA1, "label", "Apolipoprotein A1, mmol/L")
setattr( mcps$ApoB, "label", "Apolipoprotein B, mmol/L")

### Filter dataset ###
#Filter 1: diagnosed or undiagnosed diabetes
mcps1 <- mcps %>% filter(
  BASE_DIABETES==0&DRUG_D1==0&DRUG_D2==0&DRUG_D3==0&DRUG_D4==0) %>%
  filter(BASE_HBA1C<6.5|is.na(BASE_HBA1C))
#Filter 2: missing data on HbA1c, covariates and mortality
mcps2 <- mcps1 %>% filter(
  !is.na(D000), !is.na(BASE_HBA1C), 
  !is.na(MALE), !is.na(AGE),!is.na(COYOACAN), !is.na(EDU_LEVEL),
  !is.na(SMOKEGP), !is.na(ALCGP), !is.na(PHYSGP), !is.na(BMI), !is.na(WAISTC),
  !is.na(SBP1), !is.na(DBP1), !is.na(Clinical_LDL_C), !is.na(HDL_C),
  !is.na(Total_TG), !is.na(ApoA1), !is.na(ApoB))
#Filter 3: comorbidities (IHD, stroke, cancer, cirrhosis, COPD, CKD)
mcps3 <- mcps2 %>%
  filter(BASE_HEARTATTACK==0 &
           BASE_ANGINA==0 &
           BASE_STROKE==0 &
           BASE_BREASTCANCER==0 &
           BASE_LUNGCANCER==0 &
           BASE_STOMCANCER==0 &
           BASE_ORALCANCER==0 &
           BASE_PROSTATECANCER==0 &
           BASE_CERVCANCER==0 &
           BASE_OTHCANCER==0 &
           BASE_EMPHYSEMA==0 &
           BASE_CKD==0 &
           BASE_CIRR==0)
#Filter 4: aged >85 years old
mcps4 <- mcps3 %>% filter(AGE<=84)

### Population flowchart ###
nrow(mcps) #Total 
nrow(mcps)-nrow(mcps1) #Filter 1: Diabetes
nrow(mcps1) #Diabetes-free
nrow(mcps1)-nrow(mcps2) #Filter 2: Missing data
nrow(mcps2) #Complete data
nrow(mcps2)-nrow(mcps3) #Filter 3: Comorbidities
nrow(mcps3) #Without comorbidities
nrow(mcps3)-nrow(mcps4) #Filter 4: Aged >84 years
nrow(mcps4) #35-84 years old

### Manuscript overview ###
(nrow(mcps) - nrow(mcps4)) #Total number of participants removed
((nrow(mcps) - nrow(mcps4))/nrow(mcps)*100) %>% round(1) #Percentage of participants removed
nrow(mcps)-nrow(mcps1); ((nrow(mcps)-nrow(mcps1))/nrow(mcps)*100) %>% round(1) #Diabetes
nrow(mcps1)-nrow(mcps2); ((nrow(mcps1)-nrow(mcps2))/nrow(mcps)*100) %>% round(1) #Missing data
nrow(mcps2)-nrow(mcps3); ((nrow(mcps2)-nrow(mcps3))/nrow(mcps)*100) %>% round(1) #Comorbidities
nrow(mcps3)-nrow(mcps4); ((nrow(mcps3)-nrow(mcps4))/nrow(mcps)*100) %>% round(1) #Aged >95 years
nrow(mcps)-nrow(mcps) #Duplicated
nrow(mcps4) #Final sample size
(mcps4$AGE>=75) %>% table #Aged ≥75 total
((mcps4$AGE>=75) %>% table %>% prop.table*100) %>% round(1) #Aged ≥75 percentage
(with(mcps4, age_risk>=75) %>% table) #Aged ≥75 at the end of follow-up
(with(mcps4, age_risk>=75) %>% table) %>% prop.table() %>% round(3)*100 #Percentage
mcps4$hba1c_categories %>% table #A1c groups total
(mcps4$hba1c_categories %>% table %>% prop.table*100) %>% round(1) #A1c groups pecentage

mcps_main <- mcps4 %>% filter(AGE<75); nrow(mcps_main)
(mcps_main$pred_ADA %>% table)[2]; (mcps_main$pred_IEC %>% table)[2] #Prediabetes definitions
((mcps_main$pred_ADA %>% table %>% prop.table)[2]*100) %>% round(1)
((mcps_main$pred_IEC %>% table %>% prop.table)[2]*100) %>% round(1)
remove(mcps_main)

### Rename datasets ###
mcps_fin <- mcps4 #MAIN ANALYSIS
mcps_sens <- mcps2 #SENSITIVITY ANALYSIS (with comorbidities)
remove(mcps1, mcps2, mcps3)

### Recodify ###
mcps_sens$n_comorb <- mcps_sens %>% dplyr::select(
  BASE_HEARTATTACK, BASE_ANGINA, BASE_STROKE, BASE_BREASTCANCER,
  BASE_LUNGCANCER, BASE_STOMCANCER, BASE_ORALCANCER, BASE_CERVCANCER,
  BASE_PROSTATECANCER, BASE_OTHCANCER, BASE_EMPHYSEMA, BASE_CKD,
  BASE_CIRR, BASE_PAD,BASE_ASTHMA, BASE_PEP) %>% apply(1, sum)

### Death comparison ###
### Follow-up
quantile(mcps_fin$PERSON_YEARS) %>% round(1)
#Without comorbidities
death.distribution <- data.frame(
  "All"=mcps_fin  %>% summarise(
    "All"=sum(D000),"Cardiac"=sum(D003),"Cerebrovascular"=sum(D008),"Other_vascular"=sum(D014),"Renal"=sum(D019),
    "Cardiovascular"=sum(D003)+sum(D008)+sum(D014),
    "Acute_diabetic"=sum(D023), "Hepatobiliary"=sum(D022), "Other_GI"=sum(D045),
    "Neoplastic"=sum(D040),"Respiratory"=sum(D044), "External_ill"=sum(ext_ill)) %>% t,
  "Premature"=mcps_fin %>% summarise(
    "All"=sum(D000_p),"Cardiac"=sum(D003_p),"Cerebrovascular"=sum(D008_p),"Other_vascular"=sum(D014_p),"Renal"=sum(D019_p),
    "Cardiovascular"=sum(D003_p)+sum(D008_p)+sum(D014_p),
    "Acute_diabetic"=sum(D023_p),"Hepatobiliary"=sum(D022_p), "Other_GI"=sum(D045_p),
    "Neoplastic"=sum(D040_p),"Respiratory"=sum(D044_p), "External_ill"=sum(ext_ill_p)) %>% t,
  "Old"=mcps_fin %>% summarise(
    "All"=sum(D000_o),"Cardiac"=sum(D003_o),"Cerebrovascular"=sum(D008_o),"Other_vascular"=sum(D014_o),"Renal"=sum(D019_o),
    "Cardiovascular"=sum(D003_o)+sum(D008_o)+sum(D014_o),
    "Acute_diabetic"=sum(D023_o),"Hepatobiliary"=sum(D022_o), "Other_GI"=sum(D045_o),
    "Neoplastic"=sum(D040_o),"Respiratory"=sum(D044_o), "External_ill"=sum(ext_ill_o)) %>% t); death.distribution

(death.distribution[,2]*100/death.distribution[1,2]) %>% round(1) #<75
(death.distribution[,3]*100/death.distribution[1,3]) %>% round(1) #≥75
9260-5333
6495-2568

#Follow-up
mcps_fin$PERSON_YEARS %>% sum
mcps_fin$PERSON_YEARS %>% summary %>% round(1)
#Number of deaths with and without prediabetes
rbind( #ADA
  (mcps_fin %>% group_by(predm) %>% summarise("D"=sum(D000),"Prop"=sum(D000)/sum(!predm)))[1,],
  (mcps_fin %>% group_by(predm) %>% summarise("D"=sum(D000),"Prop"=sum(D000)/sum(predm)))[2,])
rbind( #IEC
  (mcps_fin %>% group_by(predm_ice) %>% summarise("D"=sum(D000),"Prop"=sum(D000)/sum(!predm_ice)))[1,],
  (mcps_fin %>% group_by(predm_ice) %>% summarise("D"=sum(D000),"Prop"=sum(D000)/sum(predm_ice)))[2,])

#With comorbidities
data.frame(
  "All"=mcps_sens %>% summarise(
    "All"=sum(D000),"Cardiac"=sum(D003),"Cerebrovascular"=sum(D008),"Other_vascular"=sum(D014),"Renal"=sum(D019),
    "Acute_diabetic"=sum(D023), "Hepatobiliary"=sum(D022), "Other_GI"=sum(D045),
    "Neoplastic"=sum(D040),"Respiratory"=sum(D044), "External_ill"=sum(ext_ill)) %>% t,
  "Premature"=mcps_sens %>% summarise(
    "All"=sum(D000_p5),"Cardiac"=sum(D003_p5),"Cerebrovascular"=sum(D008_p5),"Other_vascular"=sum(D014_p5),"Renal"=sum(D019_p5),
    "Acute_diabetic"=sum(D023_p5),"Hepatobiliary"=sum(D022_p5), "Other_GI"=sum(D045_p5),
    "Neoplastic"=sum(D040_p5),"Respiratory"=sum(D044_p5), "External_ill"=sum(ext_ill_p5)) %>% t,
  "Old"=mcps_sens %>% summarise(
    "All"=sum(D000_o5),"Cardiac"=sum(D003_o5),"Cerebrovascular"=sum(D008_o5),"Other_vascular"=sum(D014_o5),"Renal"=sum(D019_o5),
    "Acute_diabetic"=sum(D023_o5),"Hepatobiliary"=sum(D022_o5), "Other_GI"=sum(D045_o5),
    "Neoplastic"=sum(D040_o5),"Respiratory"=sum(D044_o5), "External_ill"=sum(ext_ill_o5)) %>% t)

#At-risk participants
(with(mcps_fin, age_risk<75) %>% table)[2] #At-risk for premature
(with(mcps_fin, age_risk>=75) %>% table)[2] #At risk for older
(with(mcps_sens, age_risk<75) %>% table)[2] #At-risk for premature with comorb

#Deaths occur disproportionally higher in men for cardiac and hepatobiliary causes
(mcps_fin %>% group_by(MALE) %>% summarise(
  A1=sum(D000_p)/sum(PERSON_YEARS), A2=sum(D003_p)/sum(PERSON_YEARS), A3=sum(D008_p)/sum(PERSON_YEARS),
  A4=sum(D014_p)/sum(PERSON_YEARS), A5=sum(D019_p)/sum(PERSON_YEARS), A6=sum(D023_p)/sum(PERSON_YEARS),
  A7=sum(D022_p)/sum(PERSON_YEARS), A8=sum(D045_p)/sum(PERSON_YEARS), A9=sum(D044_p)/sum(PERSON_YEARS),
  A10=sum(D040_p)/sum(PERSON_YEARS), A11=sum(ext_ill_p)/sum(PERSON_YEARS))*100)->GGG
GGG[2,]/GGG[1,]

#### Table  1: Study population----- ####
mcps_fin$icc<-mcps_fin$WAISTC/mcps_fin$HIPC
setattr(mcps_fin$icc, "label", "Waist-Hip Ratio")

tab1 <- mcps_fin %>% filter(AGE<75) %>%
  dplyr::select(hba1c_categories,AGE,BASE_HBA1C,female,
                COYOACAN,college,smoke,alcohol,phys,
                BMI,ICE,icc,SBP1,DBP1,
                Clinical_LDL_C,HDL_C,Total_TG,ApoA1,ApoB)%>%
  tbl_summary(by = hba1c_categories,missing = "ifany",
              statistic = list(all_continuous() ~ "{mean} ({sd})"))%>%
  bold_labels() %>% add_overall() %>% add_p() %>%  
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(dplyr::mutate,
    label = ifelse(label == "N missing (% missing)", "Unknown", label))%>%
  as_flex_table() %>% align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = tab1, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Proyectos/Prediabetes/tabla1.docx")

(86518*100/nrow(mcps_fin)) %>% round(1)
(23901*100/nrow(mcps_fin)) %>% round(1)
(8074*100/nrow(mcps_fin)) %>% round(1)

mcps_fin$pred_ADA %>% table %>% prop.table()
mcps_fin$pred_IEC %>% table %>% prop.table()

#### SuppTab1: Study pop by age----- ####

mcps_fin$age_decade <- cut(
  x = mcps_fin$AGE, breaks = c(-Inf,45,55,65,75,Inf), right=F,
  labels = c("35-44", "45-54", "55-64", "65-74", "≥75"))
setattr(mcps_fin$hba1c_categories, "label", "HbA1c category")

s_tab1 <- mcps_fin %>%
  dplyr::select(age_decade,AGE,BASE_HBA1C,hba1c_categories,female,
                COYOACAN,college,smoke,alcohol,phys,
                BMI,ICE,icc,SBP1,DBP1,
                Clinical_LDL_C,HDL_C,Total_TG,ApoA1,ApoB)%>%
  tbl_summary(by = age_decade, missing = "ifany",
              statistic = list(all_continuous() ~ "{mean} ({sd})"))%>%
  bold_labels() %>% add_overall() %>% add_p() %>%  
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(dplyr::mutate, label = ifelse(
    label == "N missing (% missing)", "Unknown", label))%>%
  as_flex_table() %>% align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = s_tab1, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Proyectos/Prediabetes/tabla_supp1.docx")


#### SuppMethods: ICD-10 codes------ ####
ICD10_ <- function(x){
  (as.data.frame(table(x$ICD10_UNDERLYING))%>% `names<-`(c("Code", "Num")) %>%
     transmute("A"=paste0(Code, " (", Num, ")")))$A %>% paste(collapse=", ")}

set_flextable_defaults(
  font.size = 9, font.family = "Arial", table.layout="autofit", split=F,
  table_align="center"); supp_methods <- data.frame(paste0(c(
    "Cardiac", "Cerebrovascular", "Other vascular", "Renal", "Acute diabetic",
    "Neoplastic", "Respiratory", "Gastrointestinal (including hepatobiliary)",
    "External, ill-defined, or others"), " deaths, N=", mcps_fin %>% summarise(
      sum(D003),sum(D008),sum(D014),sum(D019),sum(D023),sum(D040),sum(D044),
      sum(D022+D045),sum(ext_ill)) %>% `names<-`(NULL)),
    c(ICD10_(filter(mcps_fin, D003==1)), ICD10_(filter(mcps_fin, D008==1)),
      ICD10_(filter(mcps_fin, D014==1)), ICD10_(filter(mcps_fin, D019==1)),
      ICD10_(filter(mcps_fin, D023==1)), ICD10_(filter(mcps_fin, (D022+D045)==1)),
      ICD10_(filter(mcps_fin, D040==1)), ICD10_(filter(mcps_fin, D044==1)),
      ICD10_(filter(mcps_fin, ext_ill==1)))) %>% `colnames<-`(c(
        "Cause-specific mortality, Number of deaths", "ICD-10 codes (N)")) %>%
  flextable() %>% align(align="center", part="all") %>%
  align(align="justify", j=2) %>% italic(part="all") %>%
  bold(part="header") %>% autofit()

read_docx() %>% body_add_flextable(value=supp_methods) %>%
  print(target = "Proyectos/Prediabetes/SuppMethods.docx")

#### Figure 1: Prediabetes prevalence---- ####
mcps_fin %>% filter(!is.na(BASE_HBA1C)) %>% nrow
#ADA
prev.all.ADA_male <- round(mean(mcps_fin$pred_ADA[mcps_fin$MALE==1],na.rm=T)*100,1)
prev.all.ADA_female <- round(mean(mcps_fin$pred_ADA[mcps_fin$MALE==0],na.rm=T)*100,1)

prev_ADA<- mcps_fin %>%  mutate(category=cut(
  AGE, breaks=c(-Inf, 40, 50, 55, 60, 65, 70, 75, Inf),
  labels=c("35-39","40-49","50-54", "55-59", "60-64", "65-69", "70-74", "75-84"))) %>%
  group_by(category, MALE) %>% summarise(n=n(), predm=sum(predm, na.rm = T)) %>%
  mutate(prev=(predm/n)*100) %>% mutate(Sex=factor(MALE, labels=c("Female", "Male")))
temp<-NULL;L<-NULL;S<-NULL; for(i in 1:nrow(prev_ADA)) {
  temp[[i]]<-exactci(prev_ADA$predm[[i]],prev_ADA$n[[i]],conf.level=0.95)
  L[i]<-temp[[i]]$conf.int[1]; S[i]<-temp[[i]]$conf.int[2]}
prev_ADA$lower<-L*100; prev_ADA$upper<-S*100
prevalence_ADA<-prev_ADA %>%
  ggplot(aes(y=prev, x=category, col=Sex, group=Sex))+
  geom_point(size=2)+ geom_line(size=1)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+
  labs(x="Age groups", y="Prevalence (%)")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A"))+
  theme_bw() + ggtitle("ADA-defined prediabetes") +
  scale_y_continuous(limits = c(0,50.1))+
  geom_hline(yintercept = prev.all.ADA_female, linetype="dashed", color=c("#3B9AB2"))+
  geom_hline(yintercept = prev.all.ADA_male, linetype="dashed", color=c("#EBCC2A"))+
  theme(text=element_text(size=15),
        plot.title = element_text(face = "bold", hjust=0.5, size=15.5));prevalence_ADA

#IEC
prev.all.IEC_male <- round(mean(mcps_fin$predm_ice[mcps_fin$MALE==1],na.rm=T)*100,1)
prev.all.IEC_female <- round(mean(mcps_fin$predm_ice[mcps_fin$MALE==0],na.rm=T)*100,1)

prev_iec<- mcps_fin %>% mutate(category=cut(
  AGE, breaks=c(-Inf, 40, 45, 50, 55, 60, 65, 70, 75,Inf),
  labels=c("35-39","40-44","45-49","50-54", "55-59", "60-64", "65-69", "70-74", "75-84"))) %>%
  group_by(category, MALE) %>% summarise(n=n(), predm=sum(predm_ice, na.rm = T)) %>%
  mutate(prev=(predm/n)*100) %>% mutate(Sex=factor(MALE, labels=c("Female", "Male"))) 
temp<-NULL;L<-NULL;S<-NULL; for(i in 1:nrow(prev_iec)) {
  temp[[i]]<-exactci(prev_iec$predm[[i]],prev_iec$n[[i]],conf.level=0.95)
  L[i]<-temp[[i]]$conf.int[1]; S[i]<-temp[[i]]$conf.int[2]}
prev_iec$lower<-L*100; prev_iec$upper<-S*100
prevalence_iec<-prev_iec%>%
  ggplot(aes(y=prev, x=category, col=Sex, group=Sex))+
  geom_point(size=2)+ geom_line(size=1)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+
  labs(x="Age groups", y="Prevalence (%)")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A"))+
  theme_bw() + ggtitle("IEC-defined prediabetes") +
  scale_y_continuous(limits = c(0,50.1))+
  geom_hline(yintercept = prev.all.IEC_female, linetype="dashed", color=c("#3B9AB2"))+
  geom_hline(yintercept = prev.all.IEC_male, linetype="dashed", color=c("#EBCC2A"))+
  theme(text=element_text(size=15),
        plot.title = element_text(face = "bold", hjust=0.5, size=15.5));prevalence_iec

fig1<-ggarrange(prevalence_ADA, prevalence_iec, labels = c("A", "B"),
                font.label = list(size = 16), common.legend = T, legend = "bottom")
ggsave(fig1, file="Proyectos/Prediabetes/Figure1.tif", bg="transparent",
         width=25*1.3, height=15*1.3, units=c("cm"), dpi=300, limitsize = FALSE)

#### SuppFig2: Prediabetes risk by sex--- ####
# Lexis expansion
remove(mcps_lexis, mcps_fin2)
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

#Deaths at 35-74 years old
(mcps_fin2 %>% filter(MALE==0 & between(AGE, 35, 74)) %>% 
  select(premature,hba1c_cat) %>% table)[2,] -> labels_ev3
m0.4A_female <- coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(age.cat)+cluster(id), data=mcps_fin2 %>% filter(MALE==0))
m0.4A_female_float <- float(object = m0.4A_female); dfAf <- data.frame(
  "A1c"=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%"), "mean"=m0.4A_female_float$coef,
  "Lower"=m0.4A_female_float$coef-qnorm(.975)*(m0.4A_female_float$var %>% sqrt), "Upper"=m0.4A_female_float$coef+qnorm(.975)*(m0.4A_female_float$var %>% sqrt)) %>%
  mutate("HR"=exp(mean), "Lower"=exp(Lower), "Upper"=exp(Upper), "A1c"=ordered(A1c,levels=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%")), "psize"=1/(mean+1))  %>%
  `rownames<-`(NULL) %>% arrange(A1c) %>% mutate("N"=labels_ev3, "HR2"=sprintf("%#.2f",HR), "Sex"="Females")

(mcps_fin2 %>% filter(MALE==1 & between(AGE, 35, 74)) %>% 
    select(premature,hba1c_cat) %>% table)[2,] -> labels_ev4
m0.4A_male <- coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(age.cat)+cluster(id), data=mcps_fin2 %>% filter(MALE==1))
m0.4A_male_float <- float(object = m0.4A_male); dfAm <- data.frame(
  "A1c"=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%"), "mean"=m0.4A_male_float$coef,
  "Lower"=m0.4A_male_float$coef-qnorm(.975)*(m0.4A_male_float$var %>% sqrt), "Upper"=m0.4A_male_float$coef+qnorm(.975)*(m0.4A_male_float$var %>% sqrt)) %>%
  mutate("HR"=exp(mean), "Lower"=exp(Lower), "Upper"=exp(Upper), "A1c"=ordered(A1c,labels=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%")), "psize"=1/(mean+1))  %>%
  `rownames<-`(NULL) %>% arrange(A1c) %>% mutate("N"=labels_ev4, "HR2"=sprintf("%#.2f",HR), "Sex"="Males")

rbind(dfAf, dfAm) %>%
  ggplot(aes(x=A1c,color=A1c,y=HR,ymin=Lower,ymax=Upper,size=psize)) + geom_hline(yintercept = 1, linetype=2) +
  geom_pointrange(shape=15) + facet_wrap(~Sex, nrow=1, ncol=3) + coord_cartesian(ylim=c((0.85),(1.75))) +
  scale_color_manual(values = wes_palette("Zissou1", 3)) + theme_bw() + scale_size_continuous(guide="none", range = c(0.7,1.3)) +
  labs(color="HbA1c", y="Mortality Rate Ratio (95%CI)", x="\nHbA1c category", title=NULL) +
  scale_y_continuous(trans = scales::log_trans()) +facet_wrap(~Sex)+
  ggtitle("Deaths ocurring at ages 35-74 years")+
  theme(legend.position="top", text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text.x = element_text(face = "italic", size=11),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust=0.5, size=13.5)) +
  geom_text(aes(x=A1c,y=HR,label=HR2), color="black", nudge_x=0.2, nudge_y=0.05, size=3, fontface="bold.italic") +
  geom_text(aes(x=A1c,y=Lower,label=N), color="black", nudge_x=0, nudge_y=-0.025, size=2.75) -> fs2b.f; fs2b.f

#Deaths at ≥75 years old
(mcps_fin2 %>% filter(MALE==0 & age_risk>=75) %>% 
    select(old.death,hba1c_cat) %>% table)[2,] -> labels_ev5
m0.4A_female <- coxph(Surv(time_at_entry, time_at_exit, old.death)~hba1c_cat+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(age.cat)+cluster(id), data=mcps_fin2 %>% filter(MALE==0))
m0.4A_female_float <- float(object = m0.4A_female); dfAf_o <- data.frame(
  "A1c"=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%"), "mean"=m0.4A_female_float$coef,
  "Lower"=m0.4A_female_float$coef-qnorm(.975)*(m0.4A_female_float$var %>% sqrt), "Upper"=m0.4A_female_float$coef+qnorm(.975)*(m0.4A_female_float$var %>% sqrt)) %>%
  mutate("HR"=exp(mean), "Lower"=exp(Lower), "Upper"=exp(Upper), "A1c"=ordered(A1c,labels=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%")), "psize"=1/(mean+1))  %>%
  `rownames<-`(NULL) %>% arrange(A1c) %>% mutate("N"=labels_ev5, "HR2"=sprintf("%#.2f",HR), "Sex"="Females")

(mcps_fin2 %>% filter(MALE==1 & age_risk>=75) %>% 
    select(old.death,hba1c_cat) %>% table)[2,] -> labels_ev6
m0.4A_male <- coxph(Surv(time_at_entry, time_at_exit, old.death)~hba1c_cat+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(age.cat)+cluster(id), data=mcps_fin2 %>% filter(MALE==1))
m0.4A_male_float <- float(object = m0.4A_male); dfAm_o <- data.frame(
  "A1c"=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%"), "mean"=m0.4A_male_float$coef,
  "Lower"=m0.4A_male_float$coef-qnorm(.975)*(m0.4A_male_float$var %>% sqrt), "Upper"=m0.4A_male_float$coef+qnorm(.975)*(m0.4A_male_float$var %>% sqrt)) %>%
  mutate("HR"=exp(mean), "Lower"=exp(Lower), "Upper"=exp(Upper), "A1c"=ordered(A1c,labels=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%")), "psize"=1/(mean+1))  %>%
  `rownames<-`(NULL) %>% arrange(A1c) %>% mutate("N"=labels_ev6, "HR2"=sprintf("%#.2f",HR), "Sex"="Males")

rbind(dfAf_o, dfAm_o) %>%
  ggplot(aes(x=A1c,color=A1c,y=HR,ymin=Lower,ymax=Upper,size=psize)) + geom_hline(yintercept = 1, linetype=2) +
  geom_pointrange(shape=15) + facet_wrap(~Age, nrow=1, ncol=3) + coord_cartesian(ylim=c((0.85),(1.75))) +
  scale_color_manual(values = wes_palette("Zissou1", 3)) + theme_bw() + scale_size_continuous(guide="none", range = c(0.7,1.3)) +
  labs(color="HbA1c", y="Mortality Rate Ratio (95%CI)", x="\nHbA1c category", title=NULL) +
  scale_y_continuous(trans = scales::log_trans()) +facet_wrap(~Sex)+
  ggtitle("Deaths ocurring at ages 75-84 years")+
  theme(legend.position="top", text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text.x = element_text(face = "italic", size=11),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust=0.5, size=13.5)) +
  geom_text(aes(x=A1c,y=HR,label=HR2), color="black", nudge_x=0.2, nudge_y=0.05, size=3, fontface="bold.italic") +
  geom_text(aes(x=A1c,y=Lower,label=N), color="black", nudge_x=0, nudge_y=-0.025, size=2.75) -> fs2b.f_o; fs2b.f_o

#### Build figure 2
figs2a <- ggarrange(fs2b.f,fs2b.f_o, nrow=2, ncol=1, common.legend = T,
                    labels=c("A","B"), font.label = list(size = 17))
ggsave(figs2a, file="Proyectos/Prediabetes/FSupp2.tiff", bg="white", width=15*1.3, height=15*1.5, units=c("cm"), dpi=600, limitsize = FALSE)
remove(mcps_lexis, mcps_fin2)


#### Figure 2: (A) Death rates per cause- ####
mcps_fin$age_group <- mcps_fin$AGE %>%
  cut(c(-Inf, 45,55,65,75,85, Inf), labels = 1:6, right = F)

table(mcps_fin$age_group, mcps_fin$D000)
table(mcps_fin$age_group, mcps_fin$D000_o)
table(mcps_fin$age_group, mcps_fin$D000_p)

get.death.rate <- function(x){
  a <- x %>% summarise(cases=n(), time=sum(PERSON_YEARS))
  time <- sum(a$time); cases <- a$cases[2]
  lambda<-(cases/time)*1000; lambda.se<-(sqrt(cases/time^2))*1000
  LI<-lambda-(qnorm(.975)*lambda.se); LS<-lambda+(qnorm(.975)*lambda.se)
  cases.L<-LI*time/1000; cases.U<-LS*time/1000
  data.frame(lambda, cases, cases.L, cases.U, time)} #Crude death rates
un.ag.st <- function(x){
  a <- data.frame(); for(i in 1:length(table(mcps_fin$age_group))){
    b <- x %>% filter(age_group==i) %>% get.death.rate()
    a <- rbind(a, b)}; c <- a %>% cbind(
      "age"=mcps_fin$age_group %>% levels); c} #Uniformly age and sex standardized
un.ag.st.CI <- function(x){
  x %>% na.omit %>% mutate("w"=1/length(age)) %>%
    mutate("var"=(((w)^2*lambda)/time)) %>% summarise(
      lambda=mean(lambda), var=sum(var), lambda.se=sqrt(sum(var)),
      cases=sum(cases), cases.L=sum(cases.L), cases.U=sum(cases.U)) %>% 
    mutate("Lower"=lambda+(cases.L-cases)*sqrt(var/cases),
           "Upper"=lambda+(cases.U-cases)*sqrt(var/cases))} #Standardized confidence intervals

#------------------ PREMATURE DEATHS ------------------#
## All-cause mortality ##
#HBA1C <5.7
rates1A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D000_p) %>% un.ag.st
rates1A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D000_p) %>% un.ag.st
rates1A <- rbind(rates1A.M, rates1A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates1B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D000_p) %>% un.ag.st
rates1B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D000_p) %>% un.ag.st
rates1B <- rbind(rates1B.M, rates1B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates1C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D000_p) %>% un.ag.st
rates1C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D000_p) %>% un.ag.st
rates1C <- rbind(rates1C.M, rates1C.F) %>% un.ag.st.CI
rates1 <- rbind(rates1A,rates1B,rates1C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Overall")

## Cardiac ## D003
#HBA1C <5.7
rates2A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D003_p) %>% un.ag.st
rates2A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D003_p) %>% un.ag.st
rates2A <- rbind(rates2A.M, rates2A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates2B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D003_p) %>% un.ag.st
rates2B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D003_p) %>% un.ag.st
rates2B <- rbind(rates2B.M, rates2B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates2C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D003_p) %>% un.ag.st
rates2C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D003_p) %>% un.ag.st
rates2C <- rbind(rates2C.M, rates2C.F) %>% un.ag.st.CI
rates2 <- rbind(rates2A,rates2B,rates2C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Cardiac")

## Cerebrovascular ## D008
#HBA1C <5.7
rates3A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D008_p) %>% un.ag.st
rates3A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D008_p) %>% un.ag.st
rates3A <- rbind(rates3A.M, rates3A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates3B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D008_p) %>% un.ag.st
rates3B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D008_p) %>% un.ag.st
rates3B <- rbind(rates3B.M, rates3B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates3C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D008_p) %>% un.ag.st
rates3C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D008_p) %>% un.ag.st
rates3C <- rbind(rates3C.M, rates3C.F) %>% un.ag.st.CI
rates3 <- rbind(rates3A,rates3B,rates3C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Cerebrovascular")

## Other vascular ## D014
#HBA1C <5.7
rates4A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D014_p) %>% un.ag.st
rates4A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D014_p) %>% un.ag.st
rates4A <- rbind(rates4A.M, rates4A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates4B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D014_p) %>% un.ag.st
rates4B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D014_p) %>% un.ag.st
rates4B <- rbind(rates4B.M, rates4B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates4C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D014_p) %>% un.ag.st
rates4C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D014_p) %>% un.ag.st
rates4C <- rbind(rates4C.M, rates4C.F) %>% un.ag.st.CI
rates4 <- rbind(rates4A,rates4B,rates4C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Other vascular")

## Renal ## D019
#HBA1C <5.7
rates5A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D019_p) %>% un.ag.st
rates5A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D019_p) %>% un.ag.st
rates5A <- rbind(rates5A.M, rates5A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates5B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D019_p) %>% un.ag.st
rates5B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D019_p) %>% un.ag.st
rates5B <- rbind(rates5B.M, rates5B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates5C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D019_p) %>% un.ag.st
rates5C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D019_p) %>% un.ag.st
rates5C <- rbind(rates5C.M, rates5C.F) %>% un.ag.st.CI
rates5 <- rbind(rates5A,rates5B,rates5C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Renal")

## Acute diabetic ## D023
#HBA1C <5.7
rates6A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D023_p) %>% un.ag.st
rates6A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D023_p) %>% un.ag.st
rates6A <- rbind(rates6A.M, rates6A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates6B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D023_p) %>% un.ag.st
rates6B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D023_p) %>% un.ag.st
rates6B <- rbind(rates6B.M, rates6B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates6C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D023_p) %>% un.ag.st
rates6C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D023_p) %>% un.ag.st
rates6C <- rbind(rates6C.M, rates6C.F) %>% un.ag.st.CI
rates6 <- rbind(rates6A,rates6B,rates6C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Acute diabetic")

## Other gastrointestinal ## D045
#HBA1C <5.7
rates8A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D045_p) %>% un.ag.st
rates8A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D045_p) %>% un.ag.st
rates8A <- rbind(rates8A.M, rates8A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates8B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D045_p) %>% un.ag.st



rates8B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D045_p) %>% un.ag.st
rates8B <- rbind(rates8B.M, rates8B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates8C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D045_p) %>% un.ag.st
rates8C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D045_p) %>% un.ag.st
rates8C <- rbind(rates8C.M, rates8C.F) %>% un.ag.st.CI
rates8 <- rbind(rates8A,rates8B,rates8C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Other\ngastrointestinal")

## Respiratory ## D044
#HBA1C <5.7
rates9A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D044_p) %>% un.ag.st
rates9A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D044_p) %>% un.ag.st


rates9A <- rbind(rates9A.M, rates9A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates9B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D044_p) %>% un.ag.st
rates9B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D044_p) %>% un.ag.st
rates9B <- rbind(rates9B.M, rates9B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates9C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D044_p) %>% un.ag.st
rates9C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D044_p) %>% un.ag.st
rates9C <- rbind(rates9C.M, rates9C.F) %>% un.ag.st.CI
rates9 <- rbind(rates9A,rates9B,rates9C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Respiratory")

## Neoplastic ## D040
#HBA1C <5.7
rates10A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D040_p) %>% un.ag.st
rates10A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D040_p) %>% un.ag.st
rates10A <- rbind(rates10A.M, rates10A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates10B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D040_p) %>% un.ag.st
rates10B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D040_p) %>% un.ag.st
rates10B <- rbind(rates10B.M, rates10B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates10C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D040_p) %>% un.ag.st
rates10C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D040_p) %>% un.ag.st
rates10C <- rbind(rates10C.M, rates10C.F) %>% un.ag.st.CI
rates10 <- rbind(rates10A,rates10B,rates10C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Neoplastic")

## Other ## ext_ill
#HBA1C <5.7
rates11A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(ext_ill_p) %>% un.ag.st
rates11A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(ext_ill_p) %>% un.ag.st
rates11A <- rbind(rates11A.M, rates11A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates11B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>%group_by(ext_ill_p) %>% un.ag.st
rates11B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(ext_ill_p) %>% un.ag.st
rates11B <- rbind(rates11B.M, rates11B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates11C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(ext_ill_p) %>% un.ag.st
rates11C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(ext_ill_p) %>% un.ag.st
rates11C <- rbind(rates11C.M, rates11C.F) %>% un.ag.st.CI
rates11 <- rbind(rates11A,rates11B,rates11C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Other")

rbind(rates1, rates2, rates3, rates5,
      rates6, rates8, rates9, rates10) %>% arrange(desc(lambda)) %>%
  mutate(pred=reorder(pred, -lambda, decreasing=T), cause=reorder(cause, -lambda, decreasing=T)) %>% 
  ggplot(aes(y=lambda, fill=pred, x=cause)) + geom_bar(stat="identity", position='dodge', color="black", width=0.8) +
  theme_bw() + scale_fill_manual(values = wes_palette("Zissou1", 3)) +
  labs(fill="HbA1c", y="Deaths per 1,000 person-years (uniformly age- and sex-standardized)") +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(.8)) +
  coord_flip() + scale_y_break(c(1.5,3), scales = 0.25, ticklabels = c(3,4)) + xlab(NULL) +
  theme(legend.position="top", text=element_text(size=15), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust=0.5, size=13.5))  -> f2a.1


#------------------ DEATHS AT OLDER AGES ------------------#
## All-cause mortality ##
#HBA1C <5.7
rates1A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D000_o) %>% un.ag.st
rates1A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D000_o) %>% un.ag.st
rates1A <- rbind(rates1A.M, rates1A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates1B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D000_o) %>% un.ag.st
rates1B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D000_o) %>% un.ag.st
rates1B <- rbind(rates1B.M, rates1B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates1C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D000_o) %>% un.ag.st
rates1C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D000_o) %>% un.ag.st
rates1C <- rbind(rates1C.M, rates1C.F) %>% un.ag.st.CI
rates1o <- rbind(rates1A,rates1B,rates1C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Overall")

## Cardiac ## D003
#HBA1C <5.7
rates2A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D003_o) %>% un.ag.st
rates2A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D003_o) %>% un.ag.st
rates2A <- rbind(rates2A.M, rates2A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates2B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D003_o) %>% un.ag.st
rates2B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D003_o) %>% un.ag.st
rates2B <- rbind(rates2B.M, rates2B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates2C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D003_o) %>% un.ag.st
rates2C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D003_o) %>% un.ag.st
rates2C <- rbind(rates2C.M, rates2C.F) %>% un.ag.st.CI
rates2o <- rbind(rates2A,rates2B,rates2C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Cardiac")

## Cerebrovascular ## D008
#HBA1C <5.7
rates3A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D008_o) %>% un.ag.st
rates3A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D008_o) %>% un.ag.st
rates3A <- rbind(rates3A.M, rates3A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates3B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D008_o) %>% un.ag.st
rates3B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D008_o) %>% un.ag.st
rates3B <- rbind(rates3B.M, rates3B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates3C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D008_o) %>% un.ag.st
rates3C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D008_o) %>% un.ag.st
rates3C <- rbind(rates3C.M, rates3C.F) %>% un.ag.st.CI
rates3o <- rbind(rates3A,rates3B,rates3C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Cerebrovascular")

## Other vascular ## D014
#HBA1C <5.7
rates4A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D014_o) %>% un.ag.st
rates4A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D014_o) %>% un.ag.st
rates4A <- rbind(rates4A.M, rates4A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates4B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D014_o) %>% un.ag.st
rates4B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D014_o) %>% un.ag.st
rates4B <- rbind(rates4B.M, rates4B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates4C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D014_o) %>% un.ag.st
rates4C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D014_o) %>% un.ag.st
rates4C <- rbind(rates4C.M, rates4C.F) %>% un.ag.st.CI
rates4o <- rbind(rates4A,rates4B,rates4C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Other vascular")

## Renal ## D019
#HBA1C <5.7
rates5A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D019_o) %>% un.ag.st
rates5A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D019_o) %>% un.ag.st
rates5A <- rbind(rates5A.M, rates5A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates5B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D019_o) %>% un.ag.st
rates5B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D019_o) %>% un.ag.st
rates5B <- rbind(rates5B.M, rates5B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates5C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D019_o) %>% un.ag.st
rates5C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D019_o) %>% un.ag.st
rates5C <- rbind(rates5C.M, rates5C.F) %>% un.ag.st.CI
rates5o <- rbind(rates5A,rates5B,rates5C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Renal")

## Acute diabetic ## D023
#HBA1C <5.7
rates6A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D023_o) %>% un.ag.st
rates6A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D023_o) %>% un.ag.st
rates6A <- rbind(rates6A.M, rates6A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates6B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D023_o) %>% un.ag.st
rates6B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D023_o) %>% un.ag.st
rates6B <- rbind(rates6B.M, rates6B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates6C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D023_o) %>% un.ag.st
rates6C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D023_o) %>% un.ag.st
rates6C <- rbind(rates6C.M, rates6C.F) %>% un.ag.st.CI
rates6o <- rbind(rates6A,rates6B,rates6C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Acute diabetic")

## Other gastrointestinal ## D045
#HBA1C <5.7
rates8A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D045_o) %>% un.ag.st
rates8A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D045_o) %>% un.ag.st
rates8A <- rbind(rates8A.M, rates8A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates8B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D045_o) %>% un.ag.st
rates8B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D045_o) %>% un.ag.st
rates8B <- rbind(rates8B.M, rates8B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates8C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D045_o) %>% un.ag.st
rates8C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D045_o) %>% un.ag.st
rates8C <- rbind(rates8C.M, rates8C.F) %>% un.ag.st.CI
rates8o <- rbind(rates8A,rates8B,rates8C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Other\ngastrointestinal")

## Respiratory ## D044
#HBA1C <5.7
rates9A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D044_o) %>% un.ag.st
rates9A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D044_o) %>% un.ag.st
rates9A <- rbind(rates9A.M, rates9A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates9B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D044_o) %>% un.ag.st
rates9B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D044_o) %>% un.ag.st
rates9B <- rbind(rates9B.M, rates9B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates9C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D044_o) %>% un.ag.st
rates9C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D044_o) %>% un.ag.st
rates9C <- rbind(rates9C.M, rates9C.F) %>% un.ag.st.CI
rates9o <- rbind(rates9A,rates9B,rates9C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Respiratory")

## Neoplastic ## D040
#HBA1C <5.7
rates10A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(D040_o) %>% un.ag.st
rates10A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(D040_o) %>% un.ag.st
rates10A <- rbind(rates10A.M, rates10A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates10B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(D040_o) %>% un.ag.st
rates10B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(D040_o) %>% un.ag.st
rates10B <- rbind(rates10B.M, rates10B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates10C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(D040_o) %>% un.ag.st
rates10C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(D040_o) %>% un.ag.st
rates10C <- rbind(rates10C.M, rates10C.F) %>% un.ag.st.CI
rates10o <- rbind(rates10A,rates10B,rates10C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Neoplastic")

## Other ## ext_ill
#HBA1C <5.7
rates11A.M <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==1) %>% group_by(ext_ill_o) %>% un.ag.st
rates11A.F <- mcps_fin %>% filter(hba1c_cat=="<5.7%", MALE==0) %>% group_by(ext_ill_o) %>% un.ag.st
rates11A <- rbind(rates11A.M, rates11A.F) %>% un.ag.st.CI
#HBA1C 5.7-5.9
rates11B.M <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==1) %>% group_by(ext_ill_o) %>% un.ag.st
rates11B.F <- mcps_fin %>% filter(hba1c_cat=="5.7-5.9%", MALE==0) %>% group_by(ext_ill_o) %>% un.ag.st
rates11B <- rbind(rates11B.M, rates11B.F) %>% un.ag.st.CI
#HBA1C ≥6.0
rates11C.M <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==1) %>% group_by(ext_ill_o) %>% un.ag.st
rates11C.F <- mcps_fin %>% filter(hba1c_cat=="6.0-6.4%", MALE==0) %>% group_by(ext_ill_o) %>% un.ag.st
rates11C <- rbind(rates11C.M, rates11C.F) %>% un.ag.st.CI
rates11o <- rbind(rates11A,rates11B,rates11C) %>% as_tibble %>% mutate(
  pred=c("<5.7%","≥5.7 to <6.0%","≥6.0 to <6.5%"), cause="Other")

rbind(rates1o, rates2o, rates3o, rates5o,
      rates6o, rates8o, rates9o, rates10o) %>% arrange(desc(lambda)) %>%
  mutate(pred=reorder(pred, -lambda, decreasing=T), cause=reorder(cause, -lambda, decreasing=T)) %>% 
  ggplot(aes(y=lambda, fill=pred, x=cause)) + geom_bar(stat="identity", position='dodge', color="black", width=0.8) +
  theme_bw() + scale_fill_manual(values = wes_palette("Zissou1", 3)) +
  labs(fill="HbA1c", y="Deaths per 1,000 person-years (uniformly age- and sex-standardized)") +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(.8)) +
  coord_flip() + scale_y_break(c(7,15), scales = 0.25, ticklabels = c(15,18)) + xlab(NULL) +
  theme(legend.position="top", text=element_text(size=15), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust=0.5, size=13.5)) -> f2a.2

####
rbind(rates1, rates2, rates3, rates4, rates5, rates6, rates8, rates9, rates10, rates11) %>% group_by(cause) %>%
  mutate(RR=lambda/lambda[pred=="<5.7%"], RR.se=sqrt((1/cases[pred=="<5.7%"])+(1/cases)),
         RR.LL=RR-qnorm(.975)*RR.se,RR.UL=RR+qnorm(.975)*RR.se)
rbind(rates1o, rates2o, rates3o, rates4o, rates5o, rates6o, rates8o, rates9o, rates10o, rates11o) %>% group_by(cause) %>%
  mutate(RR=lambda/lambda[pred=="<5.7%"], RR.se=sqrt((1/cases[pred=="<5.7%"])+(1/cases)),
         RR.LL=RR-qnorm(.975)*RR.se,RR.UL=RR+qnorm(.975)*RR.se)


#### Figure 2: (B) Age*HbA1c (+SuppFig3)- ####
# Lexis expansion
remove(mcps_lexis, mcps_fin2)
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

# Baseline a1c and age categories
mcps_fin2 <- mcps_fin2 %>% mutate(
  "a1c_int"=cut(BASE_HBA1C, c(-Inf, 5.7, 6, Inf), right=F) %>% factor(labels=c("<5.7%","5.7-5.9%","6.0-6.4%")),
  "age_int"=cut(AGE, c(-Inf, 50, 64, Inf), right=F) %>% factor(labels=c("35-50","50-64","≥65")),
  "age_int2"=cut(AGE, c(-Inf, 65, 69, Inf), right=F) %>% factor(labels=c("<65","65-69","≥70")),
  "a1c_age"=paste0(a1c_int,", ",age_int) %>% factor(levels=c("<5.7%, 35-50","<5.7%, 50-64","<5.7%, ≥65",
                                                             "5.7-5.9%, 35-50","5.7-5.9%, 50-64","5.7-5.9%, ≥65",
                                                             "6.0-6.4%, 35-50","6.0-6.4%, 50-64","6.0-6.4%, ≥65")),
  "a1c_age2"=paste0(a1c_int,", ",age_int2) %>% factor(levels=c("<5.7%, <65","<5.7%, 65-69","<5.7%, ≥70",
                                                               "5.7-5.9%, <65","5.7-5.9%, 65-69","5.7-5.9%, ≥70", 
                                                               "6.0-6.4%, <65","6.0-6.4%, 65-69","6.0-6.4%, ≥70")))

#Deaths at 35-74 years old
(mcps_fin2 %>% filter(between(AGE,35,74)) %>% transmute(
  "Out"=premature, "A1c"=a1c_age, MALE, EDU_LEVEL) %>% na.omit %>%
    select(Out,A1c) %>% table)[2,] -> labels_ev1
m0.4A <- coxph(Surv(time_at_entry, time_at_exit, premature)~a1c_age+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.4A_float <- float(object = m0.4A); dfA <- data.frame(
  "A1c"=(m0.4A_float$coef%>%names) %>% str_extract("([<≥]?\\d)([.])(\\d)([-]\\d[.]\\d)?([%])"),
  "Age"=(m0.4A_float$coef%>%names) %>% str_extract("([<≥]?\\d{2})([-]\\d{2})?"), "mean"=m0.4A_float$coef,
  "Lower"=m0.4A_float$coef-qnorm(.975)*(m0.4A_float$var %>% sqrt), "Upper"=m0.4A_float$coef+qnorm(.975)*(m0.4A_float$var %>% sqrt)) %>%
  mutate("Age"=ordered(Age,levels=c("35-50","50-64","≥65"), labels=c("Baseline age <50 years","≥50 to <65 years","≥65 years"))) %>%
  mutate("HR"=exp(mean), "Lower"=exp(Lower), "Upper"=exp(Upper), "A1c"=ordered(A1c,levels=c("<5.7%","5.7-5.9%","6.0-6.4%"), labels=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%")), "psize"=1/(mean+1))  %>%
  `rownames<-`(NULL) %>% arrange(A1c,Age) %>% mutate("N"=labels_ev1, "HR2"=sprintf("%#.2f",HR)); dfA %>%
  ggplot(aes(x=A1c,color=A1c,y=HR,ymin=Lower,ymax=Upper,size=psize)) + geom_hline(yintercept = 1, linetype=2) +
  geom_pointrange(shape=15) + facet_wrap(~Age, nrow=1, ncol=3) + coord_cartesian(ylim=c((0.85),(2.15))) +
  scale_color_manual(values = wes_palette("Zissou1", 3)) + theme_bw() + scale_size_continuous(guide="none", range = c(0.7,1.3)) +
  labs(color="HbA1c", y="Mortality Rate Ratio (95%CI)", x="\nAge-HbA1c category", title=NULL) +
  scale_y_continuous(trans = scales::log_trans()) +
  theme(legend.position="top", text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text.x = element_text(face = "italic", size=11),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust=0.5, size=13.5)) +
  geom_text(aes(x=A1c,y=HR,label=HR2), color="black", nudge_x=0.2, nudge_y=0.05, size=3, fontface="bold.italic") +
  geom_text(aes(x=A1c,y=Lower,label=N), color="black", nudge_x=0, nudge_y=-0.025, size=2.75) -> f2b.1; f2b.1

#Deaths at ≥75 years old
(mcps_fin2 %>% filter(age_risk>=75) %>% transmute(
  "Out"=old.death, "A1c"=a1c_age2, MALE, EDU_LEVEL) %>% na.omit %>%
    select(Out,A1c) %>% table)[2,] -> labels_ev2
m0.4B <- coxph(Surv(time_at_entry, time_at_exit, old.death)~a1c_age2+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.4B_float <- float(object = m0.4B); dfB <- data.frame(
  "A1c"=(m0.4B_float$coef%>%names) %>% str_extract("([<≥]?\\d)([.])(\\d)([-]\\d[.]\\d)?([%])"),
  "Age"=(m0.4B_float$coef%>%names) %>% str_extract("([<≥]?\\d{2})([-]\\d{2})?"), "mean"=m0.4B_float$coef,
  "Lower"=m0.4B_float$coef-qnorm(.975)*(m0.4B_float$var %>% sqrt), "Upper"=m0.4B_float$coef+qnorm(.975)*(m0.4B_float$var %>% sqrt)) %>%
  mutate("Age"=ordered(Age,levels=c("<65","65-69","≥70"), labels=c("Baseline age <65 years","≥65 to <70 years","≥70 years"))) %>%
  mutate("HR"=exp(mean), "Lower"=exp(Lower), "Upper"=exp(Upper), "A1c"=ordered(A1c,levels=c("<5.7%","5.7-5.9%","6.0-6.4%"), labels=c("<5.7%","≥5.7% to <6.0%","≥6.0% to <6.5%")), "psize"=1/(mean+1))  %>%
  `rownames<-`(NULL) %>% arrange(A1c,Age) %>% mutate("N"=labels_ev2, "HR2"=sprintf("%.2f",HR)); dfB %>%
  ggplot(aes(x=A1c,color=A1c,y=HR,ymin=Lower,ymax=Upper,size=psize)) + geom_hline(yintercept = 1, linetype=2) +
  geom_pointrange(shape=15) + facet_wrap(~Age, nrow=1, ncol=3) + coord_cartesian(ylim=c(0.85,2.15)) +
  scale_color_manual(values = wes_palette("Zissou1", 3)) + theme_bw() + scale_size_continuous(guide="none", range = c(0.7,1.3))+
  labs(color="HbA1c", y="Mortality Rate Ratio (95%CI)", x="\nAge-HbA1c category", title=NULL) +
  scale_y_continuous(trans = scales::log_trans()) +
  theme(legend.position="top", text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text.x = element_text(face = "italic", size=11),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust=0.5, size=13.5)) +
  geom_text(aes(x=A1c,y=HR,label=HR2), color="black", nudge_x=0.2, nudge_y=0.05, size=3, fontface="bold.italic") +
  geom_text(aes(x=A1c,y=Lower,label=N), color="black", nudge_x=0, nudge_y=-0.025, size=2.75) -> f2b.2; f2b.2

#### Build figure 2
fig2 <- ggarrange(print(f2a.1), f2b.1, nrow=1, ncol=2, labels=c("A","B"), font.label = list(size = 17)) 
ggsave(fig2, file="Proyectos/Prediabetes/Figure2.tif", bg="white", width=30*1.3, height=15*1.5, units=c("cm"), dpi=300, limitsize = FALSE)

#### Build figure 2
supp2 <- ggarrange(print(f2a.2), f2b.2, nrow=1, ncol=2, labels=c("A","B"), font.label = list(size = 17)) %>%
  annotate_figure(top=text_grob("Deaths ocurring at ages 75-84 years", face="bold.italic", size=18, hjust=0.5))
ggsave(supp2, file="Proyectos/Prediabetes/FSupp3.png", bg="white", width=30*1.3, height=15*1.5, units=c("cm"), dpi=600, limitsize = FALSE)


#### RevDCare: A1c*Age*Sex + splines----- ####
### TEST TRIPLE INTERACTION
# Lexis expansion
remove(mcps_lexis, mcps_fin2)
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)


#DC revision: Quantify HbA1c*Age and HbA1c*Sex interactions
#Premature deaths
coxph( #HbA1c*Sex
  Surv(time_at_entry,time_at_exit,premature)~BASE_HBA1C*MALE+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+ factor(ALCGP)+COYOACAN+
    strata(age.cat), data=mcps_fin2) -> test.int_1
coxph( #HbA1c*AGE
  Surv(time_at_entry,time_at_exit,premature)~BASE_HBA1C*AGE+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+ factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE), data=mcps_fin2) -> test.int_2

#Older deaths
coxph( #HbA1c*Sex
  Surv(time_at_entry,time_at_exit,old.death)~BASE_HBA1C*MALE+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+ factor(ALCGP)+COYOACAN+
    strata(age.cat), data=mcps_fin2) -> test.int_4
coxph( #HbA1c*AGE
  Surv(time_at_entry,time_at_exit,old.death)~BASE_HBA1C*AGE+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+ factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE), data=mcps_fin2) -> test.int_5

#Premature deaths
(test.int_1 %>% summary)$coefficients[c(1:2,15),] #A1c*Sex
(test.int_2 %>% summary)$coefficients[c(1:2,15),] #A1c*Age
#Older deaths
(test.int_4 %>% summary)$coefficients[c(1:2,15),] #A1c*Sex
(test.int_5 %>% summary)$coefficients[c(1:2,15),] #A1c*Age



### SPLINES ###
### Lexis expansion ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
test_poly3<-coxph(
  Surv(time_at_entry, time_at_exit, premature)~poly(BASE_HBA1C,3)+
    BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+
    COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
test_rcs3<-coxph(
  Surv(time_at_entry, time_at_exit, premature)~rcs(BASE_HBA1C,3)+
    BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+
    COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
test_rcs4<-coxph(
  Surv(time_at_entry, time_at_exit, premature)~rcs(BASE_HBA1C,4)+
    BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+
    COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
test_rcs5<-coxph(
  Surv(time_at_entry, time_at_exit, premature)~rcs(BASE_HBA1C,5)+
    BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+
    COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)

test_poly3 %>% summary
test_rcs3 %>% summary
test_rcs4 %>% summary
test_rcs5 %>% summary


#### Tab2: Deaths at 35-74y -------------- ####
### Lexis expansion ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG

mcps_fin2<- mcps_fin2 %>% filter(Clinical_LDL_C>0, Total_TG>0, HDL_C>0)

#Sex, municipality
m0<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat  + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA   + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC   + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle
m0_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat  + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA   + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC   + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle, adiposity
m0_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C  +BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m0.5_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat   +BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m2_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA    +BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m3_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC    +BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C + BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat  + BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA   + BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC   + BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)

### Table quick ###
quicktab <- function(a,b,c,d,x){
  rbind(
    (a %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[2,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (c %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (d %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f")) %>%
    as.data.frame %>% `names<-`(LETTERS[1:3]) %>%
    transmute(A=paste0(A, " (", B, ", ", C, ")")) %>% `names<-`(x)}

qt1 <- data.frame(
  quicktab(m0.5, m1, m2, m3, "Baseline"),
  quicktab(m0.5_adj1, m1_adj1, m2_adj1, m3_adj1, "Lifestyle"),
  quicktab(m0.5_adj2, m1_adj2, m2_adj2, m3_adj2, "Adiposity"),
  quicktab(m0.5_adj3, m1_adj3, m2_adj3, m3_adj3, "Adiposity + BP/Lipids")
  ); rbind(qt1[1,], "1.00", qt1[2:5,]) %>%
  `rownames<-`(c("HbA1c (%)", "<5.7%", "5.7-5.9%", "6.0-6.4%",
                 "ADA (5.7-6.4%)", "IEC (6.0-6.4%)")) %>% View


### Regression table ###
#Continuous HBA1C pero 1.0% increase
#models <- list(m0,m0_adj1,m0_adj2)
#fit_model0 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), label="BASE_HBA1C"~"HbA1c (%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model0)
#row0 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Continuous HBA1C per 0.5% increase
#models <- list(m0.5,m0.5_adj1,m0.5_adj2)
#fit_model0 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), label="BASE_HBA1C"~"HbA1c (%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model0)
#row0.5 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#HBA1C categories
#models <- list(m1,m1_adj1,m1_adj2)
#fit_model1 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("hba1c_cat"), label="hba1c_cat"~"HbA1c levels") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model1)
#row1 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Prediabetes ADA
#models <- list(m2,m2_adj1,m2_adj2)
#fit_model2 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("pred_ADA"), label="pred_ADA"~"ADA (≥5.7%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model2)
#row2 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Prediabetes IEC
#models <- list(m3,m3_adj1,m3_adj2)
#fit_model3 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("pred_IEC"), label="pred_IEC"~"IEC (≥6.0%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model3)
#row3 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Save table
#tab2 <-tbl_stack(list(row0, row1, row2, row3), group_header = c("Continuous", "Categorical", "Prediabetes definition", ""))%>% as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
#doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%  body_end_section_landscape() %>%
#  print(target = "Proyectos/Prediabetes/tabla2.5.docx"); remove(row0, row1, row2, row3, results, mcps_fin2)




#### SuppTab3: Participants w/comorb---------- ####
### Recodify ###
mcps_sens$n_comorb <- mcps_sens %>% dplyr::select(
  BASE_HEARTATTACK, BASE_ANGINA, BASE_STROKE, BASE_BREASTCANCER,
  BASE_LUNGCANCER, BASE_STOMCANCER, BASE_ORALCANCER, BASE_CERVCANCER,
  BASE_PROSTATECANCER, BASE_OTHCANCER, BASE_EMPHYSEMA, BASE_CKD, BASE_CIRR,
  BASE_PAD,BASE_ASTHMA, BASE_PEP) %>% apply(1, sum)
(mcps_sens %>% filter(AGE<75)) %>% nrow

### Lexis expansion ###
mcps_sens <- mcps_sens %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_sens)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C + COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat  + COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA   + COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC   + COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle
m0_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat  + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA   + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC   + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle, BP, lipids
m0_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat  + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA   + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC   + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m0.5_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat   +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m2_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_ADA    +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m3_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~pred_IEC    +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+n_comorb+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)

### Table quick ###
quicktab <- function(a,b,c,d,x){
  rbind(
    (a %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[2,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (c %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (d %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f")) %>%
    as.data.frame %>% `names<-`(LETTERS[1:3]) %>%
    transmute(A=paste0(A, " (", B, ", ", C, ")")) %>% `names<-`(x)}

qt1 <- data.frame(
  quicktab(m0.5, m1, m2, m3, "Baseline"),
  quicktab(m0.5_adj1, m1_adj1, m2_adj1, m3_adj1, "Lifestyle"),
  quicktab(m0.5_adj2, m1_adj2, m2_adj2, m3_adj2, "Adiposity"),
  quicktab(m0.5_adj3, m1_adj3, m2_adj3, m3_adj3, "Adiposity + BP/Lipids")
); rbind(qt1[1,], "1.00", qt1[2:5,]) %>%
  `rownames<-`(c("HbA1c (%)", "<5.7%", "5.7-5.9%", "6.0-6.4%",
                 "ADA (5.7-6.4%)", "IEC (6.0-6.4%)")) %>% View

### Regression table ###
#Continuous HBA1C per 1.0% increase
#models <- list(m0,m0_adj1,m0_adj2)
#fit_model0 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), label="BASE_HBA1C"~"HbA1c (%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model0)
#row0 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Continuous HBA1C per 0.5% increase
#models <- list(m0.5,m0.5_adj1,m0.5_adj2)
#fit_model0 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), label="BASE_HBA1C"~"HbA1c (%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model0)
#row0.5 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#HBA1C categories
#models <- list(m1,m1_adj1,m1_adj2)
#fit_model1 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("hba1c_cat"), label="hba1c_cat"~"HbA1c levels") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#registerDoParallel(cores = detectCores()); results <- foreach(model = models) %dopar% {fit_model1(model)}
#row1 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Prediabetes ADA
#models <- list(m2,m2_adj1,m2_adj2)
#fit_model2 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("pred_ADA"), label="pred_ADA"~"ADA (≥5.7%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#registerDoParallel(cores = detectCores()); results <- foreach(model = models) %dopar% {fit_model2(model)}
#row2 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Prediabetes IEC
#models <- list(m3,m3_adj1,m3_adj2)
#fit_model3 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("pred_IEC"), label="pred_IEC"~"IEC (≥6.0%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#registerDoParallel(cores = detectCores()); results <- foreach(model = models) %dopar% {fit_model3(model)}
#row3 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Save table
#tab2 <-tbl_stack(list(row0, row0.5, row1, row2, row3), group_header = c("Continuous", "Categorical", "Prediabetes definition", ""))%>% as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
#doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%  body_end_section_landscape() %>%
#  print(target = "Proyectos/Prediabetes/tabla_supp2.docx"); remove(row0, row1, row2, row3, results, mcps_fin2)


#### SuppTab4: Deaths at ≥75y----------------- ####
mcps_fin %>% filter(age_risk>=75) %>% nrow
mcps_fin %>% filter(age_risk>=75) %>% select(D000) %>% table
mcps_fin %>% filter(age_risk<75) %>% nrow
mcps_fin %>% filter(age_risk<75) %>% select(D000) %>% table

### Lexis expansion ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~hba1c_cat  + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_ADA   + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_IEC   + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle
m0_adj1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5_adj1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~hba1c_cat  + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2_adj1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_ADA   + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3_adj1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_IEC   + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle, BP, lipids
m0_adj2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m0.5_adj2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~hba1c_cat  + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m2_adj2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_ADA   + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
m3_adj2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_IEC   + BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0_adj3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m0.5_adj3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
m1_adj3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~hba1c_cat   +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m2_adj3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_ADA    +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)
m3_adj3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~pred_IEC    +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2)

### Table quick ###
quicktab <- function(a,b,c,d,x){
  rbind(
    (a %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[2,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (c %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (d %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f")) %>%
    as.data.frame %>% `names<-`(LETTERS[1:3]) %>%
    transmute(A=paste0(A, " (", B, ", ", C, ")")) %>% `names<-`(x)}

qt1 <- data.frame(
  quicktab(m0.5, m1, m2, m3, "Baseline"),
  quicktab(m0.5_adj1, m1_adj1, m2_adj1, m3_adj1, "Lifestyle"),
  quicktab(m0.5_adj2, m1_adj2, m2_adj2, m3_adj2, "Adiposity"),
  quicktab(m0.5_adj3, m1_adj3, m2_adj3, m3_adj3, "Adiposity + BP/Lipids")
); rbind(qt1[1,], "1.00", qt1[2:5,]) %>%
  `rownames<-`(c("HbA1c (%)", "<5.7%", "5.7-5.9%", "6.0-6.4%",
                 "ADA (5.7-6.4%)", "IEC (6.0-6.4%)")) %>% View

### Regression table ###
#Continuous HBA1C per 1.0% increase
#models <- list(m0,m0_adj1,m0_adj2)
#fit_model0 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), label="BASE_HBA1C"~"HbA1c (%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model0)
#row0 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Continuous HBA1C per 0.5% increase
#models <- list(m0.5,m0.5_adj1,m0.5_adj2)
#fit_model0 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), label="BASE_HBA1C"~"HbA1c (%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model0)
#row0.5 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#HBA1C categories
#models <- list(m1,m1_adj1,m1_adj2)
#fit_model1 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("hba1c_cat"), label="hba1c_cat"~"HbA1c levels") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model1)
#row1 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Prediabetes ADA
#models <- list(m2,m2_adj1,m2_adj2)
#fit_model2 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("pred_ADA"), label="pred_ADA"~"ADA (≥5.7%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model2)
#row2 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Prediabetes IEC
#models <- list(m3,m3_adj1,m3_adj2)
#fit_model3 <- function(model, group) {
#  model %>% tbl_regression(exponentiate = TRUE, include = c("pred_IEC"), label="pred_IEC"~"IEC (≥6.0%)") %>%
#    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
#results <- lapply(models,fit_model3)
#row3 <- tbl_merge(results, tab_spanner = c("Model 1", "Model 2","Model 3"))
#Save table
#tab2 <-tbl_stack(list(row0, row0.5, row1, row2, row3), group_header = c("Continuous", "Categorical", "Prediabetes definition", ""))%>% as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
#doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%  body_end_section_landscape() %>%
#  print(target = "Proyectos/Prediabetes/tabla_supp1.docx"); remove(row0, row1, row2, row3, results, mcps_fin2)

#### SuppTab6: sample size and deaths table--- ####
#DC revision: Sample size and deaths for each stratum
# Lexis expansion
remove(mcps_lexis, mcps_fin2)
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

# Baseline a1c and age categories
mcps_fin2 <- mcps_fin2 %>% mutate(
  "a1c_int"=cut(BASE_HBA1C, c(-Inf, 5.7, 6, Inf), right=F) %>% factor(labels=c("<5.7%","5.7-5.9%","6.0-6.4%")),
  "age_int"=cut(AGE, c(-Inf, 50, 64, Inf), right=F) %>% factor(labels=c("35-50","50-64","≥65")),
  "age_int2"=cut(AGE, c(-Inf, 65, 69, Inf), right=F) %>% factor(labels=c("<65","65-69","≥70")),
  "a1c_age"=paste0(a1c_int,", ",age_int) %>% factor(levels=c("<5.7%, 35-50","<5.7%, 50-64","<5.7%, ≥65",
                                                             "5.7-5.9%, 35-50","5.7-5.9%, 50-64","5.7-5.9%, ≥65",
                                                             "6.0-6.4%, 35-50","6.0-6.4%, 50-64","6.0-6.4%, ≥65")),
  "a1c_age2"=paste0(a1c_int,", ",age_int2) %>% factor(levels=c("<5.7%, <65","<5.7%, 65-69","<5.7%, ≥70",
                                                               "5.7-5.9%, <65","5.7-5.9%, 65-69","5.7-5.9%, ≥70", 
                                                               "6.0-6.4%, <65","6.0-6.4%, 65-69","6.0-6.4%, ≥70")))


##-- Baseline age --##
#Deaths at 35-74
mcps_tab <- mcps_fin2 %>% filter(between(AGE,35,74)) %>% transmute(
  "Out"=premature, "A1c"=a1c_age, MALE, EDU_LEVEL, id); table(filter(
    mcps_tab, !duplicated(id)) %>% select(A1c)) %>% cbind((
      mcps_tab %>% na.omit %>% select(Out,A1c) %>% table)[2,]) %>%
  data.frame %>% `names<-`(c("A","B")) %>% #apply(2,sum)
  transmute("Overall"=paste0(format(A,big.mark=",")," (", B,")")) -> st4.A
#Deaths at 75-84
mcps_tab <- mcps_fin2 %>% filter(age_risk>=75) %>% transmute(
  "Out"=old.death, "A1c"=a1c_age2, MALE, EDU_LEVEL, id); table(filter(
    mcps_tab, !duplicated(id)) %>% select(A1c)) %>% cbind((
      mcps_tab %>% na.omit %>% select(Out,A1c) %>% table)[2,]) %>%
  data.frame %>% `names<-`(c("A","B")) %>% #apply(2,sum)
  transmute("Overall"=paste0(format(A,big.mark=",")," (", B,")")) -> st4.B

st4.A %>% View
st4.B %>% View

##-- Sex --##
#Female, deaths at 35-74 years
mcps_tab <- mcps_fin2 %>% filter(between(AGE,35,74),MALE==0) %>% transmute(
  "Out"=premature, "A1c"=a1c_int, MALE, EDU_LEVEL, id); table(filter(
    mcps_tab, !duplicated(id)) %>% select(A1c)) %>% cbind((
      mcps_tab %>% na.omit %>% select(Out,A1c) %>% table)[2,]) %>%
  data.frame %>% `names<-`(c("A","B")) %>% 
  transmute("Female"=paste0(format(A,big.mark=",")," (", B,")")) -> st4.C
#Male, deaths at 35-74 years
mcps_tab <- mcps_fin2 %>% filter(between(AGE,35,74),MALE==1) %>% transmute(
  "Out"=premature, "A1c"=a1c_int, MALE, EDU_LEVEL, id); table(filter(
    mcps_tab, !duplicated(id)) %>% select(A1c)) %>% cbind((
      mcps_tab %>% na.omit %>% select(Out,A1c) %>% table)[2,]) %>%
  data.frame %>% `names<-`(c("A","B")) %>% 
  transmute("Male"=paste0(format(A,big.mark=",")," (", B,")")) -> st4.D
c(st4.C[1,], st4.D[1,], st4.C[2,], st4.D[2,], st4.C[3,], st4.D[3,]) %>%
  t %>% t %>% View

#Female, deaths at 75-84 years
mcps_tab <- mcps_fin2 %>% filter(age_risk>=75,MALE==0) %>% transmute(
  "Out"=old.death, "A1c"=a1c_int, MALE, EDU_LEVEL, id); table(filter(
    mcps_tab, !duplicated(id)) %>% select(A1c)) %>% cbind((
      mcps_tab %>% na.omit %>% select(Out,A1c) %>% table)[2,]) %>%
  data.frame %>% `names<-`(c("A","B")) %>% 
  transmute("Female"=paste0(format(A,big.mark=",")," (", B,")")) -> st4.E
#Male, deaths at 75-84 years
mcps_tab <- mcps_fin2 %>% filter(age_risk>=75,MALE==1) %>% transmute(
  "Out"=old.death, "A1c"=a1c_int, MALE, EDU_LEVEL, id); table(filter(
    mcps_tab, !duplicated(id)) %>% select(A1c)) %>% cbind((
      mcps_tab %>% na.omit %>% select(Out,A1c) %>% table)[2,]) %>%
  data.frame %>% `names<-`(c("A","B")) %>% 
  transmute("Male"=paste0(format(A,big.mark=",")," (", B,")")) -> st4.F
c(st4.E[1,], st4.F[1,], st4.E[2,], st4.F[2,], st4.E[3,], st4.F[3,]) %>%
  t %>% t %>% View

remove(mcps_lexis, mcps_fin2)


#### Rev JCEM: A1c*Adipose markers (SuppFig4)- ####
# Normal vs High BMI and Waist stratifications
mcps_fin$bmi_cat <- with(mcps_fin, case_when(
  BMI<25~"Normal", BMI>=30~"High"))
mcps_fin$waist_cat <- with(mcps_fin, case_when(
  (WAISTC<80&MALE==0)~"Normal", (WAISTC<90&MALE==1)~"Normal",
  (WAISTC>=80&MALE==0)~"High", (WAISTC>=90&MALE==1)~"High"))

#ADA - BMI
mcps_fin %>% select(bmi_cat, pred_ADA) %>% table
mcps_fin %>% select(bmi_cat, pred_ADA) %>% table %>% prop.table(1) %>% round(3)*100

#IEC - BMI
mcps_fin %>% select(bmi_cat, pred_IEC) %>% table
mcps_fin %>% select(bmi_cat, pred_IEC) %>% table %>% prop.table(1) %>% round(3)*100


#ADA - Waist
mcps_fin %>% select(waist_cat, pred_ADA) %>% table
mcps_fin %>% select(waist_cat, pred_ADA) %>% table %>% prop.table(1) %>% round(3)*100

#IEC - Waist
mcps_fin %>% select(waist_cat, pred_IEC) %>% table
mcps_fin %>% select(waist_cat, pred_IEC) %>% table %>% prop.table(1) %>% round(3)*100



# Lexis expansion
remove(mcps_lexis, mcps_fin2)
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

#JCEM revision: Check for A1c*BMI and A1c*WC interactions
coxph( #A1C*BMI PREMATURE
  Surv(time_at_entry,time_at_exit,premature)~BASE_HBA1C*BMI+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> inter1
coxph( #A1C*WC PREMATURE
  Surv(time_at_entry,time_at_exit,premature)~BASE_HBA1C*WAISTC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> inter2
coxph( #A1C*BMI OLDER
  Surv(time_at_entry,time_at_exit,old.death)~BASE_HBA1C*BMI+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> inter3
coxph( #A1C*WC OLDER
  Surv(time_at_entry,time_at_exit,old.death)~BASE_HBA1C*WAISTC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> inter4


summary(inter1)$coefficients[15,5] %>% round(4) -> ann1
summary(inter2)$coefficients[15,5] %>% round(4) -> ann2
summary(inter3)$coefficients[15,5] %>% round(4) -> ann3
summary(inter4)$coefficients[15,5] %>% round(4) -> ann4



#JCEM revision: Stratify by normal vs high BMI/WC
coxph( #ADA Normal BMI
  Surv(time_at_entry,time_at_exit,premature)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> cox1
coxph( #IEC Normal BMI
  Surv(time_at_entry,time_at_exit,premature)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> cox2
coxph( #ADA High BMI
  Surv(time_at_entry,time_at_exit,premature)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="High")) -> cox3
coxph( #IEC High BMI
  Surv(time_at_entry,time_at_exit,premature)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="High")) -> cox4

coxph( #ADA Normal WC
  Surv(time_at_entry,time_at_exit,premature)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="Normal")) -> cox5
coxph( #IEC Normal WC
  Surv(time_at_entry,time_at_exit,premature)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="Normal")) -> cox6
coxph( #ADA High WC
  Surv(time_at_entry,time_at_exit,premature)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="High")) -> cox7
coxph( #IEC High WC
  Surv(time_at_entry,time_at_exit,premature)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="High")) -> cox8

coxph( #ADA Normal BMI
  Surv(time_at_entry,time_at_exit,old.death)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> cox9
coxph( #IEC Normal BMI
  Surv(time_at_entry,time_at_exit,old.death)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="Normal")) -> cox10
coxph( #ADA High BMI
  Surv(time_at_entry,time_at_exit,old.death)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="High")) -> cox11
coxph( #IEC High BMI
  Surv(time_at_entry,time_at_exit,old.death)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(bmi_cat=="High")) -> cox12

coxph( #ADA Normal WC
  Surv(time_at_entry,time_at_exit,old.death)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="Normal")) -> cox13
coxph( #IEC Normal WC
  Surv(time_at_entry,time_at_exit,old.death)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="Normal")) -> cox14
coxph( #ADA High WC
  Surv(time_at_entry,time_at_exit,old.death)~pred_ADA+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="High")) -> cox15
coxph( #IEC High WC
  Surv(time_at_entry,time_at_exit,old.death)~pred_IEC+
    EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+
    strata(age.cat)+strata(MALE),
  data=mcps_fin2 %>% filter(waist_cat=="High")) -> cox16


(mcps_fin2 %>% filter(pred_ADA==1, bmi_cat=="Normal") %>% select(D000_p) %>% table)[2] -> n1
(mcps_fin2 %>% filter(pred_IEC==1, bmi_cat=="Normal") %>% select(D000_p) %>% table)[2] -> n2
(mcps_fin2 %>% filter(pred_ADA==1, bmi_cat=="High") %>% select(D000_p) %>% table)[2] -> n3
(mcps_fin2 %>% filter(pred_IEC==1, bmi_cat=="High") %>% select(D000_p) %>% table)[2] -> n4
(mcps_fin2 %>% filter(pred_ADA==1, waist_cat=="Normal") %>% select(D000_p) %>% table)[2] -> n5
(mcps_fin2 %>% filter(pred_IEC==1, waist_cat=="Normal") %>% select(D000_p) %>% table)[2] -> n6
(mcps_fin2 %>% filter(pred_ADA==1, waist_cat=="High") %>% select(D000_p) %>% table)[2] -> n7
(mcps_fin2 %>% filter(pred_IEC==1, waist_cat=="High") %>% select(D000_p) %>% table)[2] -> n8
(mcps_fin2 %>% filter(pred_ADA==1, bmi_cat=="Normal") %>% select(D000_o) %>% table)[2] -> n9
(mcps_fin2 %>% filter(pred_IEC==1, bmi_cat=="Normal") %>% select(D000_o) %>% table)[2] -> n10
(mcps_fin2 %>% filter(pred_ADA==1, bmi_cat=="High") %>% select(D000_o) %>% table)[2] -> n11
(mcps_fin2 %>% filter(pred_IEC==1, bmi_cat=="High") %>% select(D000_o) %>% table)[2] -> n12
(mcps_fin2 %>% filter(pred_ADA==1, waist_cat=="Normal") %>% select(D000_o) %>% table)[2] -> n13
(mcps_fin2 %>% filter(pred_IEC==1, waist_cat=="Normal") %>% select(D000_o) %>% table)[2] -> n14
(mcps_fin2 %>% filter(pred_ADA==1, waist_cat=="High") %>% select(D000_o) %>% table)[2] -> n15
(mcps_fin2 %>% filter(pred_IEC==1, waist_cat=="High") %>% select(D000_o) %>% table)[2] -> n16

rbind(summary(cox1)$conf.int[1,c(1,3:4)], summary(cox2)$conf.int[1,c(1,3:4)],
      summary(cox3)$conf.int[1,c(1,3:4)], summary(cox4)$conf.int[1,c(1,3:4)],
      summary(cox5)$conf.int[1,c(1,3:4)], summary(cox6)$conf.int[1,c(1,3:4)],
      summary(cox7)$conf.int[1,c(1,3:4)], summary(cox8)$conf.int[1,c(1,3:4)],
      summary(cox9)$conf.int[1,c(1,3:4)],  summary(cox10)$conf.int[1,c(1,3:4)],
      summary(cox11)$conf.int[1,c(1,3:4)], summary(cox12)$conf.int[1,c(1,3:4)],
      summary(cox13)$conf.int[1,c(1,3:4)], summary(cox14)$conf.int[1,c(1,3:4)],
      summary(cox15)$conf.int[1,c(1,3:4)], summary(cox16)$conf.int[1,c(1,3:4)]
      ) %>% as.data.frame() %>% `colnames<-`(c("mean", "low", "up")) %>%
  cbind(
    "N"=paste0("n",1:16) %>% sapply(get),
    "age"=c(rep("Premature",8),rep("Older",8)), "def"=rep(c("ADA","IEC"),8),
    "fat"=rep(c(rep("Normal BMI",2), rep("High BMI",2),
            rep("Normal WC",2), rep("High WC",2)),2)) %>%
  mutate("def"=ordered(def, levels=c("ADA", "IEC")),
         "fat"=ordered(fat, levels=c("Normal BMI", "High BMI",
                                     "Normal WC", "High WC"))) -> base_sf4 

gp1 <- "#5697AF"; gp2 <- "#E7CD4F"
base_sf4 %>% mutate("psize"=1/(mean+1)) %>% filter(age=="Premature") %>%
  ggplot(aes(x=def, y=mean, ymin=low, ymax=up, color=def, size=psize)) +
  geom_hline(yintercept = 1, linetype=2) + geom_pointrange(shape=15) +
  scale_color_manual(values=c(gp1, gp2)) + theme_bw() +
  scale_size_continuous(guide="none", range = c(0.7,1.3)) +
  facet_wrap(~fat, ncol=4) + 
  labs(color="Prediabetes definition",
       y="Mortality Rate Ratio (95%CI)", x="", 
       title="Deaths occurring at ages 35-74 years") +
  scale_y_continuous(trans = scales::log_trans()) +
  theme(legend.position="top", text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text.x = element_text(face = "italic", size=11),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust=0.5, size=13.5))+
  geom_text(aes(x=def,y=mean,label=sprintf("%#.2f", mean)), color="black", nudge_x=0.2, nudge_y=0.05, size=3, fontface="bold.italic") +
  geom_text(aes(x=def,y=low,label=N), color="black", nudge_x=0, nudge_y=-0.025, size=2.75) -> supp4A

base_sf4 %>% mutate("psize"=1/(mean+1)) %>% filter(age=="Older") %>%
  ggplot(aes(x=def, y=mean, ymin=low, ymax=up, color=def, size=psize)) +
  geom_hline(yintercept = 1, linetype=2) + geom_pointrange(shape=15) +
  scale_color_manual(values=c(gp1, gp2)) + theme_bw() +
  scale_size_continuous(guide="none", range = c(0.7,1.3)) +
  facet_wrap(~fat, ncol=4) + 
  labs(color="Prediabetes definition",
       y="Mortality Rate Ratio (95%CI)", x="", 
       title="Deaths occurring at ages 75-84 years") +
  scale_y_continuous(trans = scales::log_trans()) +
  theme(legend.position="top", text=element_text(size=15), axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text.x = element_text(face = "italic", size=11),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust=0.5, size=13.5))+
  geom_text(aes(x=def,y=mean,label=sprintf("%#.2f", mean)), color="black", nudge_x=0.2, nudge_y=0.05, size=3, fontface="bold.italic") +
  geom_text(aes(x=def,y=low,label=N), color="black", nudge_x=0, nudge_y=-0.025, size=2.75) -> supp4B

figs4 <- ggarrange(supp4A,supp4B, nrow=2, ncol=1, common.legend = T,
                    labels=c("A","B"), font.label = list(size = 17))
ggsave(figs4, file="Proyectos/Prediabetes/FSupp4.png", bg="white",
       width=15*1.3, height=15*1.5, units=c("cm"), dpi=600, limitsize = FALSE)

  
#### Cause  2: Cardiac deaths---------------------- ####
remove(mcps_lexis, mcps_fin2)
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D003, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

c1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
c8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause  3: Cerebrovascular deaths-------------- ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D008, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

st1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
st8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause  4: Other vascular deaths--------------- ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D014, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

ov1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ov8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause4.5: Cardiovascular deaths--------------- ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D003+D008+D014, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

cvd1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
cvd2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
cvd3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
cvd4<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause  5: Renal deaths------------------------ ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D019, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

k1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
k8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

summary(k6)

#### Cause  6: Acute diabetic deaths--------------- ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D023, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

d1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
d8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause  7: Neoplastic deaths------------------- ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D040, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

ca1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ca8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause  8: Respiratory deaths------------------ ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D044, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

re1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
re8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause  9: External/ill/other------------------ ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = ext_ill, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

ei1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei7<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
ei8<-coxph(Surv(time_at_entry, time_at_exit, old.wo5)~hba1c_cat +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+factor(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#### Cause 10: Additional classifications---------- ####
#Cerebrovascular and other vascular deaths
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D008+D014, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

so3<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
so4<-coxph(Surv(time_at_entry, time_at_exit, pre.wo5)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
so5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
so6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)

#Hepatobiliary and other GI deaths
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D022+D045, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2$pre.wo5 <- mcps_fin2$premature; mcps_fin2$pre.wo5[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #35-74 w/o first 5y
mcps_fin2$old.wo5 <- mcps_fin2$old.death; mcps_fin2$old.wo5[mcps_fin2$old.death==1 & mcps_fin2$PERSON_YEARS<=5]<-0 #≥75 w/o first 5y
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, old.wo5, premature, pre.wo5)) %>% sapply(table)

eg3<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
eg4<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
eg5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm     +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
eg6<-coxph(Surv(time_at_entry, time_at_exit, old.death)~predm_ice +EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2)
remove(mcps_lexis, mcps_fin2)


#### Figure 3: Cause-specific deaths (40-74y)----------- ####
summ1c <- summary(c1)  ; summ1d <- summary(c2)
summ2c <- summary(st1) ; summ2d <- summary(st2)
summ3c <- summary(ov1) ; summ3d <- summary(ov2)
summ4c <- summary(k1)  ; summ4d <- summary(k2)
summ5c <- summary(d1)  ; summ5d <- summary(d2)
summ7c <- summary(eg3) ; summ7d <- summary(eg4)
summ8c <- summary(ca1) ; summ8d <- summary(ca2)
summ9c <- summary(re1) ; summ9d <- summary(re2)
summ10c <- summary(ei1); summ10d <- summary(ei2)

## Deaths by prediabetes definition ##
deaths_ADA.2 <- mcps_fin %>% group_by(predm) %>% summarise(
  "Cardiac"=sum(D003_p),"Cerebrovascular"=sum(D008_p),"Other_vascular"=sum(D014_p),"Renal"=sum(D019_p),
  "Acute_diabetic"=sum(D023_p), "Hepatobiliary_other_GI"=sum(D045_p)+sum(D022_p),
  "Neoplastic"=sum(D040_p),"Respiratory"=sum(D044_p), "External_ill"=sum(ext_ill_p)) %>%
  pivot_longer(!predm, names_to = "Cause", values_to = "Deaths")
deaths_IEC.2 <- mcps_fin %>% group_by(predm_ice) %>% summarise(
  "Cardiac"=sum(D003_p),"Cerebrovascular"=sum(D008_p),"Other_vascular"=sum(D014_p),"Renal"=sum(D019_p),
  "Acute_diabetic"=sum(D023_p), "Hepatobiliary_other_GI"=sum(D022_p)+sum(D045_p),
  "Neoplastic"=sum(D040_p),"Respiratory"=sum(D044_p), "External_ill"=sum(ext_ill_p)) %>%
  pivot_longer(!predm_ice, names_to = "Cause", values_to = "Deaths")
d.ADA.2 <- with(deaths_ADA.2, paste0(Deaths[predm==1],"/",Deaths[predm==0]))
d.IEC.2 <- with(deaths_IEC.2, paste0(Deaths[predm_ice==1],"/",Deaths[predm_ice==0]))
d.labs <- c("Cardiac",NA,NA,"Cerebrovascular",NA,NA,"Other vascular",NA,NA,"Renal",NA,NA,"Acute diabetic",NA,NA,
            "Gastrointestinal",NA,NA,"Neoplastic",NA,NA,"Respiratory",NA,NA,"External/Ill-defined/Other",NA,NA)

## Build forestplot ##
gp1 <- gpar(fill="#5697AF"); gp2 <- gpar(fill="#E7CD4F"); F.TAB <- rbind(
  c(summ1c$conf.int[1,c(1,3:4)]),c(summ1d$conf.int[1,c(1,3:4)]), c(summ2c$conf.int[1,c(1,3:4)]),c(summ2d$conf.int[1,c(1,3:4)]),
  c(summ3c$conf.int[1,c(1,3:4)]),c(summ3d$conf.int[1,c(1,3:4)]), c(summ4c$conf.int[1,c(1,3:4)]),c(summ4d$conf.int[1,c(1,3:4)]),
  c(summ5c$conf.int[1,c(1,3:4)]),c(summ5d$conf.int[1,c(1,3:4)]), c(summ7c$conf.int[1,c(1,3:4)]),c(summ7d$conf.int[1,c(1,3:4)]),
  c(summ8c$conf.int[1,c(1,3:4)]),c(summ8d$conf.int[1,c(1,3:4)]),
  c(summ9c$conf.int[1,c(1,3:4)]),c(summ9d$conf.int[1,c(1,3:4)]), c(summ10c$conf.int[1,c(1,3:4)]),c(summ10d$conf.int[1,c(1,3:4)])) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper")) %>% mutate("Cause"=c(1,1,2,2,3,3,4,4,5,5,6,6,8,8,9,9,10,10), "Cutoff"=rep(
    c("ADA (≥5.7)","IEC (≥6.0)"),9)) %>% arrange(Cutoff) %>% mutate("Deaths"=c(d.ADA.2, d.IEC.2)) %>% arrange(Cause) %>%
  mutate("HR"=paste0(sprintf("%#.2f", mean), " (",sprintf("%#.2f", lower), "-", sprintf("%#.2f", upper), ")")); rbind(
    NA,F.TAB[1:2,],NA,F.TAB[3:4,],NA,F.TAB[5:6,],NA,F.TAB[7:8,],NA,F.TAB[9:10,],
    NA,F.TAB[11:12,],NA,F.TAB[13:14,],NA,F.TAB[15:16,],NA,F.TAB[17:18,]) %>%
  mutate(Cause=d.labs) %>% forestplot(
    labeltext = c(Cause, Cutoff, Deaths, HR), graph.pos = 4, xlog = T, is.summary = rep(c(T,F,F), 9),
    xlab="RR for death at 35-74 years\n(with vs. without prediabetes)",
    col = fpColors(lines = "black", zero = "black"), shapes_gp = fpShapesGp(box=append(list(gp1), rep(list(gp1,gp1,gp2),9)))) %>%
  fp_add_header(Cause="Cause of death", Cutoff=c("Prediabetes\ndefinition"), Deaths="No. deaths with/\nwithout prediabetes",HR =c("RR (95% CI)")) %>%
  fp_set_style(align="c", txt_gp = fpTxtGp(summary = gpar(fontfamily = "Arial", fontface="bold.italic"), label = gpar(
    fontfamily = "Arial"), ticks = gpar(cex=1.035))) %>% fp_set_zebra_style("grey90","white") %>%
  fp_add_lines(h_2=gpar(lty = 1), h_5=gpar(lty = 1), h_8=gpar(lty = 1), h_11=gpar(lty = 1), h_14=gpar(lty = 1),
               h_17=gpar(lty = 1), h_20=gpar(lty = 1), h_23=gpar(lty = 1), h_26=gpar(lty = 1), h_29=gpar(lty = 1)) %>%
  fp_set_style(txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab  = gpar(cex=1.15, llineheight=0.7, fontface="bold.italic"))) -> forest1; forest1

png("Proyectos/Prediabetes/Figure3.tiff", width=31.75, height=35, units = "cm", res = 300); forest1
dev.off(); print(plot(1))


#### SuppFig5: Cause-specific deaths at ≥75y------------ ####
summ1e <- summary(c5)  ; summ1f <- summary(c6)
summ2e <- summary(st5) ; summ2f <- summary(st6)
summ3e <- summary(ov5) ; summ3f <- summary(ov6)
summ4e <- summary(k5)  ; summ4f <- summary(k6)
summ5e <- summary(d5)  ; summ5f <- summary(d6)
summ7c <- summary(eg5) ; summ7d <- summary(eg6)
summ8e <- summary(ca5) ; summ8f <- summary(ca6)
summ9e <- summary(re5) ; summ9f <- summary(re6)
summ10e <- summary(ei5); summ10f <- summary(ei6)

## Deaths by prediabetes definition ##
deaths_ADA.3 <- mcps_fin %>% group_by(predm) %>% filter(D000_o==1) %>% summarise(
  "Cardiac"=sum(D003_o),"Cerebrovascular"=sum(D008_o),"Other_vascular"=sum(D014_o),"Renal"=sum(D019_o),
  "Acute_diabetic"=sum(D023_o), "Hepatobiliary"=sum(D022_o)+sum(D045_o),
  "Neoplastic"=sum(D040_o),"Respiratory"=sum(D044_o), "External_ill"=sum(ext_ill_o)) %>%
  pivot_longer(!predm, names_to = "Cause", values_to = "Deaths")
deaths_IEC.3 <- mcps_fin %>% group_by(predm_ice) %>% filter(D000_o==1) %>% summarise(
  "Cardiac"=sum(D003_o),"Cerebrovascular"=sum(D008_o),"Other_vascular"=sum(D014_o),"Renal"=sum(D019_o),
  "Acute_diabetic"=sum(D023_o), "Hepatobiliary"=sum(D022_o)+sum(D045_o),
  "Neoplastic"=sum(D040_o),"Respiratory"=sum(D044_o), "External_ill"=sum(ext_ill_o)) %>%
  pivot_longer(!predm_ice, names_to = "Cause", values_to = "Deaths")
d.ADA.3 <- with(deaths_ADA.3, paste0(Deaths[predm==1],"/",Deaths[predm==0]))
d.IEC.3 <- with(deaths_IEC.3, paste0(Deaths[predm_ice==1],"/",Deaths[predm_ice==0]))
d.labs <- c("Cardiac",NA,NA,"Cerebrovascular",NA,NA,"Other vascular",NA,NA,"Renal",NA,NA,"Acute diabetic",NA,NA,
            "Gastrointestinal",NA,NA,"Neoplastic",NA,NA,"Respiratory",NA,NA,"External/Ill-defined/Other",NA,NA)

## Build forestplot ##
gp1 <- gpar(fill="#5697AF"); gp2 <- gpar(fill="#E7CD4F"); F.TAB <- rbind(
  c(summ1e$conf.int[1,c(1,3:4)]),c(summ1f$conf.int[1,c(1,3:4)]), c(summ2e$conf.int[1,c(1,3:4)]),c(summ2f$conf.int[1,c(1,3:4)]),
  c(summ3e$conf.int[1,c(1,3:4)]),c(summ3f$conf.int[1,c(1,3:4)]), c(summ4e$conf.int[1,c(1,3:4)]),c(summ4f$conf.int[1,c(1,3:4)]),
  c(summ5e$conf.int[1,c(1,3:4)]),c(summ5f$conf.int[1,c(1,3:4)]),
  c(summ7c$conf.int[1,c(1,3:4)]),c(summ7d$conf.int[1,c(1,3:4)]), c(summ8e$conf.int[1,c(1,3:4)]),c(summ8f$conf.int[1,c(1,3:4)]),
  c(summ9e$conf.int[1,c(1,3:4)]),c(summ9f$conf.int[1,c(1,3:4)]), c(summ10e$conf.int[1,c(1,3:4)]),c(summ10f$conf.int[1,c(1,3:4)])) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper")) %>% mutate("Cause"=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9), "Cutoff"=rep(
    c("ADA (≥5.7)","IEC (≥6.0)"),9)) %>% arrange(Cutoff) %>% mutate("Deaths"=c(d.ADA.3, d.IEC.3)) %>% arrange(Cause) %>%
  mutate("HR"=paste0(sprintf("%#.2f", mean), " (",sprintf("%#.2f", lower), "-", sprintf("%#.2f", upper), ")")); rbind(
    NA,F.TAB[1:2,],NA,F.TAB[3:4,],NA,F.TAB[5:6,],NA,F.TAB[7:8,],NA,F.TAB[9:10,],
    NA,F.TAB[11:12,],NA,F.TAB[13:14,],NA,F.TAB[15:16,],NA,F.TAB[17:18,]) %>%
  mutate(Cause=d.labs) %>% forestplot(
    labeltext = c(Cause, Cutoff, Deaths, HR), graph.pos = 4, xlog = T, is.summary = rep(c(T,F,F), 9),
    xlab="RR for death at 75-84 years\n(with vs. without prediabetes)",
    col = fpColors(lines = "black", zero = "black"), shapes_gp = fpShapesGp(box=append(list(gp1), rep(list(gp1,gp1,gp2),9)))) %>%
  fp_add_header(Cause="Cause of death", Cutoff=c("Prediabetes\ndefinition"), Deaths="No. deaths with/\nwithout prediabetes",HR =c("RR (95% CI)")) %>%
  fp_set_style(align="c", txt_gp = fpTxtGp(summary = gpar(fontfamily = "Arial", fontface="bold.italic"), label = gpar(
    fontfamily = "Arial"), ticks = gpar(cex=1.035))) %>% fp_set_zebra_style("grey90","white") %>%
  fp_add_lines(h_2=gpar(lty = 1), h_5=gpar(lty = 1), h_8=gpar(lty = 1), h_11=gpar(lty = 1), h_14=gpar(lty = 1),
               h_17=gpar(lty = 1), h_20=gpar(lty = 1), h_23=gpar(lty = 1), h_26=gpar(lty = 1), h_29=gpar(lty = 1)) %>%
  fp_set_style(txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab  = gpar(cex=1.15, llineheight=0.7, fontface="bold.italic"))) -> forest2; forest2

png("Proyectos/Prediabetes/FSupp5.png", width=27, height=40.5, units = "cm", res = 300); forest2
dev.off(); print(plot(1))

#### SuppTab5: Continous HbA1c and Deaths at 35-74y and ≥75y -------------- ####
### All-cause deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D000, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Cardiac deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D003, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5c<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5c<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1c<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1c<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2c<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2c<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3c<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3c<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Cerebrovascular deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D008, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5s<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5s<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1s<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1s<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2s<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2s<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3s<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3s<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Other vascular deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D014, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5ov<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5ov<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1ov<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1ov<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2ov<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2ov<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3ov<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3ov<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Cardiovascular deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D003+D008+D014, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5cvd<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5cvd<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1cvd<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1cvd<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2cvd<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2cvd<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3cvd<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3cvd<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Renal deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D019, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5k<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5k<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1k<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1k<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2k<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2k<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3k<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3k<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Acute diabetic deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D023, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5d<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5d<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1d<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1d<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2d<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2d<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3d<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3d<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Neoplastic deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D040, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5n<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5n<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1n<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1n<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2n<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2n<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3n<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3n<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))

### Respiratory deaths ###
mcps_fin <- mcps_fin %>% mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), exit.status = D044, data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
with(mcps_fin2, list(exit_status, old.death, premature)) %>% sapply(table)

### Cox models ###
#Recode
mcps_fin2$pred_ADA <-mcps_fin2$predm; mcps_fin2$pred_IEC <-mcps_fin2$predm_ice
mcps_fin2$LDL <- mcps_fin2$Clinical_LDL_C; mcps_fin2$HDL <- mcps_fin2$HDL_C
mcps_fin2$SBP <- mcps_fin2$SBP1; mcps_fin2$DBP <- mcps_fin2$DBP1
mcps_fin2$TG <- mcps_fin2$Total_TG
#Sex, municipality
m0.5r<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5r<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle
m0.5_adj1r<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj1r<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids
m0.5_adj2r<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj2r<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5 + SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+strata(age.cat)+cluster(id), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
#Sex, municipality, educational level, lifestyle, BP, lipids, adiposity
m0.5_adj3r<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))
om0.5_adj3r<-coxph(Surv(time_at_entry, time_at_exit, old.death)~BASE_HBA1C.5  +BMI+scale(ICE)+SBP+DBP+LDL+HDL+TG+EDU_LEVEL+factor(PHYSGP)+factor(SMOKEGP)+factor(ALCGP)+COYOACAN+strata(MALE)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% mutate(BASE_HBA1C.5=BASE_HBA1C/0.5))



### Table quick ###
quicktab <- function(a,b,c,d,e,f,g,h,i,x){
  rbind(
    (a %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (b %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (c %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (d %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (e %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (f %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (g %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (h %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f"),
    (i %>% summary)$conf.int[1,c(1,3:4)] %>% sprintf(fmt="%#.2f")) %>%
    as.data.frame %>% `names<-`(LETTERS[1:3]) %>%
    transmute(A=paste0(A, " (", B, ", ", C, ")")) %>% `names<-`(x)}

qt1 <- data.frame("Outcome"=c("All-cause deaths", "Cardiac deaths", "Cerebrovascular deaths", "Other vascular","Cardiovascular deaths",
                 "Renal deaths", "Acute diabetic deaths", "Neoplastic deahts", "Respiratory deaths"),
                 quicktab(m0.5, m0.5c, m0.5s, m0.5ov,m0.5cvd,m0.5k,m0.5d,m0.5n,m0.5r, "Baseline"),
                 quicktab(m0.5_adj1, m0.5_adj1c, m0.5_adj1s, m0.5_adj1ov,m0.5_adj1cvd,m0.5_adj1k,m0.5_adj1d,m0.5_adj1n,m0.5_adj1r, "Lifestyle"),
                 quicktab(m0.5_adj2, m0.5_adj2c, m0.5_adj2s, m0.5_adj2ov,m0.5_adj2cvd,m0.5_adj2k,m0.5_adj2d,m0.5_adj2n,m0.5_adj2r, "BP+Lipids"),
                 quicktab(m0.5_adj3, m0.5_adj3c, m0.5_adj3s, m0.5_adj3ov,m0.5_adj3cvd,m0.5_adj3k,m0.5_adj3d,m0.5_adj3n,m0.5_adj3r, "Adiposity"),
                 quicktab(om0.5, om0.5c, om0.5s, om0.5ov,om0.5cvd,om0.5k,om0.5d,om0.5n,om0.5r, "Baseline"),
                 quicktab(om0.5_adj1, om0.5_adj1c, om0.5_adj1s, om0.5_adj1ov,om0.5_adj1cvd,om0.5_adj1k,om0.5_adj1d,om0.5_adj1n,om0.5_adj1r, "Lifestyle"),
                 quicktab(om0.5_adj2, om0.5_adj2c, om0.5_adj2s, om0.5_adj2ov,om0.5_adj2cvd,om0.5_adj2k,om0.5_adj2d,om0.5_adj2n,om0.5_adj2r, "BP+Lipids"),
                 quicktab(om0.5_adj3, om0.5_adj3c, om0.5_adj3s, om0.5_adj3ov,om0.5_adj3cvd,om0.5_adj3k,om0.5_adj3d,om0.5_adj3n,om0.5_adj3r, "Adiposity")) %>% flextable() %>% align(align = "center",part = "all") %>% autofit()

doc <- read_docx() %>% body_add_flextable(value = qt1, split = TRUE) %>%  body_end_section_landscape() %>%
 print(target = "Proyectos/Prediabetes/tablaS5.1.docx"); remove(mcps_fin2)

#### Table  3: Excess risk (deaths at 40-74y)----------- ####
#1. RR Prediabetes
RR.1A <- (m2_adj1$coefficients[1])  %>% exp; RR.1B <- (m3_adj1$coefficients[1])  %>% exp #All-cause
RR.2A <- (c1$coefficients[1])  %>% exp; RR.2B <- (c2$coefficients[1])  %>% exp #Cardiac
RR.3A <- (st1$coefficients[1]) %>% exp; RR.3B <- (st2$coefficients[1]) %>% exp #Stroke
RR.4A <- (ov1$coefficients[1]) %>% exp; RR.4B <- (ov2$coefficients[1]) %>% exp #Other vascular
RR.5A <- (cvd1$coefficients[1]) %>% exp; RR.5B <- (cvd2$coefficients[1]) %>% exp #Other vascular
RR.6A <- (k1$coefficients[1])  %>% exp; RR.6B <- (k2$coefficients[1])  %>% exp #Renal
RR.7A <- (d1$coefficients[1])  %>% exp; RR.7B <- (d2$coefficients[1])  %>% exp #Acute diabetic
RR.8A <- (eg3$coefficients[1]) %>% exp; RR.8B <- (eg4$coefficients[1]) %>% exp #Gastrointestinal
RR.9A <- (ca1$coefficients[1]) %>% exp; RR.9B <- (ca2$coefficients[1]) %>% exp #Neoplastic
RR.10A <- (re1$coefficients[1]) %>% exp; RR.10B <- (re2$coefficients[1]) %>% exp #Respiratory
RR.11A <- (ei1$coefficients[1]) %>% exp; RR.11B <- (ei2$coefficients[1]) %>% exp #Ext/Ill/Other

#Confidence intervals Prediabetes
CI.1A  <- confint(m2_adj1)[1,] %>% exp; CI.1B  <- confint(m3_adj1)[1,] %>% exp #All cause
CI.2A  <- confint(c1)[1,] %>% exp; CI.2B  <- confint(c2)[1,] %>% exp #Cardiac
CI.3A  <- confint(st1)[1,] %>% exp; CI.3B  <- confint(st2)[1,] %>% exp #Stroke
CI.4A  <- confint(ov1)[1,] %>% exp; CI.4B  <- confint(ov2)[1,] %>% exp #Other vascular
CI.5A  <- confint(cvd1)[1,] %>% exp; CI.5B  <- confint(cvd2)[1,] %>% exp #Other vascular
CI.6A  <- confint(k1)[1,] %>% exp; CI.6B  <- confint(k2)[1,] %>% exp #Renal
CI.7A  <- confint(d1)[1,] %>% exp; CI.7B  <- confint(d2)[1,] %>% exp #Acute diabetic
CI.8A  <- confint(eg3)[1,] %>% exp; CI.8B  <- confint(eg4)[1,] %>% exp #Gastrointestinal
CI.9A  <- confint(ca1)[1,] %>% exp; CI.9B  <- confint(ca2)[1,] %>% exp #Neoplastic
CI.10A  <- confint(re1)[1,] %>% exp; CI.10B  <- confint(re2)[1,] %>% exp #Respiratory
CI.11A  <- confint(ei1)[1,] %>% exp; CI.11B  <- confint(ei2)[1,] %>% exp #Ext/Ill/Other

deaths_ADA.table <- mcps_fin %>% group_by(predm) %>% summarise("All cause"=sum(D000_p),
  "Cardiac"=sum(D003_p),"Cerebrovascular"=sum(D008_p),"Other_vascular"=sum(D014_p),"Cardiovascular"=sum(D003_p+D008_p+D014_p),"Renal"=sum(D019_p),
  "Acute_diabetic"=sum(D023_p), "Hepatobiliary_other_GI"=sum(D045_p)+sum(D022_p),
  "Neoplastic"=sum(D040_p),"Respiratory"=sum(D044_p), "External_ill"=sum(ext_ill_p)) %>%
  pivot_longer(!predm, names_to = "Cause", values_to = "Deaths")
deaths_IEC.table <- mcps_fin %>% group_by(predm_ice) %>% summarise("All cause"=sum(D000_p),
  "Cardiac"=sum(D003_p),"Cerebrovascular"=sum(D008_p),"Other_vascular"=sum(D014_p),"Cardiovascular"=sum(D003_p+D008_p+D014_p),"Renal"=sum(D019_p),
  "Acute_diabetic"=sum(D023_p), "Hepatobiliary_other_GI"=sum(D045_p)+sum(D022_p),
  "Neoplastic"=sum(D040_p),"Respiratory"=sum(D044_p), "External_ill"=sum(ext_ill_p)) %>%
  pivot_longer(!predm_ice, names_to = "Cause", values_to = "Deaths")

#2. TABLE
ER_data <- data.frame(
  (deaths_ADA.table %>% filter(predm==0) %>% select(Deaths)),
  (deaths_ADA.table %>% filter(predm==1) %>% select(Deaths)),
  (paste0("RR.",1:11,"A") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t),
  (paste0("CI.",1:11,"A") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t),
  (deaths_IEC.table %>% filter(predm_ice==0) %>% select(Deaths)),
  (deaths_IEC.table %>% filter(predm_ice==1) %>% select(Deaths)),
  (paste0("RR.",1:11,"B") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t),
  (paste0("CI.",1:11,"B") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t)) %>%
  `colnames<-`(c("Deaths.NoADA","Deaths.ADA","RR.ADA","Lo.ADA","Up.ADA","Deaths.NoIEC","Deaths.IEC","RR.IEC","Lo.IEC","Up.IEC")) %>%
  `rownames<-`(NULL) %>% 
  mutate("RR.IC.ADA"=paste0(round(RR.ADA,2)," (",round(Lo.ADA,2),"-",round(Up.ADA,2),")"),"RR.IC.IEC"=paste0(round(RR.IEC,2)," (",round(Lo.IEC,2),"-",round(Up.IEC,2),")"))%>%
  mutate("AM.ADA"=round(Deaths.ADA/(Deaths.ADA+Deaths.NoADA)*((RR.ADA-1)/(RR.ADA))*100,2)) %>%
  mutate("AM.IEC"=round(Deaths.IEC/(Deaths.IEC+Deaths.NoIEC)*((RR.IEC-1)/(RR.IEC))*100,2)) %>%
  select("Deaths.NoADA","Deaths.ADA","RR.IC.ADA","AM.ADA","Deaths.NoIEC","Deaths.IEC","RR.IC.IEC","AM.IEC")

t1<-ER_data %>% 
  flextable()%>% 
  align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = t1, split = TRUE) %>%
  body_end_section_landscape() %>% print(target = "Proyectos/Prediabetes/tabla2.docx")

#### SuppTabX: Excess risk (deaths at 75-84y)----------- ####
#1. RR Prediabetes
RR.1A <- (m2_adj3$coefficients[1])  %>% exp; RR.1B <- (m3_adj3$coefficients[1])  %>% exp #All-cause
RR.2A <- (c3$coefficients[1])  %>% exp; RR.2B <- (c4$coefficients[1])  %>% exp #Cardiac
RR.3A <- (st3$coefficients[1]) %>% exp; RR.3B <- (st4$coefficients[1]) %>% exp #Stroke
RR.4A <- (ov3$coefficients[1]) %>% exp; RR.4B <- (ov4$coefficients[1]) %>% exp #Other vascular
RR.5A <- (cvd3$coefficients[1]) %>% exp; RR.5B <- (cvd4$coefficients[1]) %>% exp #Other vascular
RR.6A <- (k3$coefficients[1])  %>% exp; RR.6B <- (k4$coefficients[1])  %>% exp #Renal
RR.7A <- (d3$coefficients[1])  %>% exp; RR.7B <- (d4$coefficients[1])  %>% exp #Acute diabetic
RR.8A <- (eg5$coefficients[1]) %>% exp; RR.8B <- (eg6$coefficients[1]) %>% exp #Gastrointestinal
RR.9A <- (ca3$coefficients[1]) %>% exp; RR.9B <- (ca4$coefficients[1]) %>% exp #Neoplastic
RR.10A <- (re3$coefficients[1]) %>% exp; RR.10B <- (re4$coefficients[1]) %>% exp #Respiratory
RR.11A <- (ei3$coefficients[1]) %>% exp; RR.11B <- (ei4$coefficients[1]) %>% exp #Ext/Ill/Other

#Confidence intervals Prediabetes
CI.1A  <- confint(m2_adj3)[1,] %>% exp; CI.1B  <- confint(m3_adj3)[1,] %>% exp #All cause
CI.2A  <- confint(c3)[1,] %>% exp; CI.2B  <- confint(c4)[1,] %>% exp #Cardiac
CI.3A  <- confint(st3)[1,] %>% exp; CI.3B  <- confint(st4)[1,] %>% exp #Stroke
CI.4A  <- confint(ov3)[1,] %>% exp; CI.4B  <- confint(ov4)[1,] %>% exp #Other vascular
CI.5A  <- confint(ov3)[1,] %>% exp; CI.5B  <- confint(ov4)[1,] %>% exp #Other vascular
CI.6A  <- confint(k3)[1,] %>% exp; CI.6B  <- confint(k4)[1,] %>% exp #Renal
CI.7A  <- confint(d3)[1,] %>% exp; CI.7B  <- confint(d4)[1,] %>% exp #Acute diabetic
CI.8A  <- confint(eg5)[1,] %>% exp; CI.8B  <- confint(eg6)[1,] %>% exp #Gastrointestinal
CI.9A  <- confint(ca3)[1,] %>% exp; CI.9B  <- confint(ca4)[1,] %>% exp #Neoplastic
CI.10A  <- confint(re3)[1,] %>% exp; CI.10B  <- confint(re4)[1,] %>% exp #Respiratory
CI.11A  <- confint(ei3)[1,] %>% exp; CI.11B  <- confint(ei4)[1,] %>% exp #Ext/Ill/Other

deaths_ADA.table <- mcps_fin %>% group_by(predm) %>% filter(D000_o==1) %>%
  summarise("All cause"=sum(D000_o),
            "Cardiac"=sum(D003_o),"Cerebrovascular"=sum(D008_o),"Other_vascular"=sum(D014_o),"Cardiovascular"=sum(D003_o+D008_o+D014_o),"Renal"=sum(D019_o),
            "Acute_diabetic"=sum(D023_o), "Hepatobiliary_other_GI"=sum(D045_o)+sum(D022_o),
            "Neoplastic"=sum(D040_o),"Respiratory"=sum(D044_o), "External_ill"=sum(ext_ill_o)) %>%
  pivot_longer(!predm, names_to = "Cause", values_to = "Deaths")
deaths_IEC.table <- mcps_fin %>% group_by(predm_ice) %>% filter(D000_o==1) %>%
  summarise("All cause"=sum(D000_o),
            "Cardiac"=sum(D003_o),"Cerebrovascular"=sum(D008_o),"Other_vascular"=sum(D014_o),"Cardiovascular"=sum(D003_o+D008_o+D014_o),"Renal"=sum(D019_o),
            "Acute_diabetic"=sum(D023_o), "Hepatobiliary_other_GI"=sum(D045_o)+sum(D022_o),
            "Neoplastic"=sum(D040_o),"Respiratory"=sum(D044_o), "External_ill"=sum(ext_ill_o)) %>%
  pivot_longer(!predm_ice, names_to = "Cause", values_to = "Deaths")

#2. TABLE
ER_data <- data.frame(
  (deaths_ADA.table %>% filter(predm==0) %>% select(Deaths)),
  (deaths_ADA.table %>% filter(predm==1) %>% select(Deaths)),
  (paste0("RR.",1:11,"A") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t),
  (paste0("CI.",1:11,"A") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t),
  (deaths_IEC.table %>% filter(predm_ice==0) %>% select(Deaths)),
  (deaths_IEC.table %>% filter(predm_ice==1) %>% select(Deaths)),
  (paste0("RR.",1:11,"B") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t),
  (paste0("CI.",1:11,"B") %>% as.list %>% lapply(get) %>% as.data.frame() %>% t)) %>%
  `colnames<-`(c("Deaths.NoADA","Deaths.ADA","RR.ADA","Lo.ADA","Up.ADA","Deaths.NoIEC","Deaths.IEC","RR.IEC","Lo.IEC","Up.IEC")) %>%
  `rownames<-`(NULL) %>% 
  mutate("RR.IC.ADA"=paste0(round(RR.ADA,2)," (",round(Lo.ADA,2),"-",round(Up.ADA,2),")"),"RR.IC.IEC"=paste0(round(RR.IEC,2)," (",round(Lo.IEC,2),"-",round(Up.IEC,2),")"))%>%
  mutate("AM.ADA"=round(Deaths.ADA/(Deaths.ADA+Deaths.NoADA)*((RR.ADA-1)/(RR.ADA))*100,2)) %>%
  mutate("AM.IEC"=round(Deaths.IEC/(Deaths.IEC+Deaths.NoIEC)*((RR.IEC-1)/(RR.IEC))*100,2)) %>%
  select("Deaths.NoADA","Deaths.ADA","RR.IC.ADA","AM.ADA","Deaths.NoIEC","Deaths.IEC","RR.IC.IEC","AM.IEC")

t1<-ER_data %>% 
  flextable()%>% 
  align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = t1, split = TRUE) %>%
  body_end_section_landscape() %>% print(target = "Proyectos/Prediabetes/SuppTable5.docx")

#### Resurvey ####
mcps_r<-mcps_sens %>% filter(!is.na(R_COYOACAN)) %>%
  filter(ifelse(R_HXDIAB_YR<=YEAR_RECRUITED, 1, 0)==0 | is.na(ifelse(R_HXDIAB_YR<=YEAR_RECRUITED, 1, 0))) %>%
  filter(BASE_CIRR==0) %>% filter(BASE_CKD==0)
nrow(mcps_r)
mcps_r$R_HXDIAB[is.na(mcps_r$R_HXDIAB)]<-0
mcps_r$diab_inc_fin<-ifelse((mcps_r$R_HXDIAB+(mcps_r$R_HBA1C>=6.5))!=0,1,0)
mcps_r$predm<- ifelse(mcps_r$BASE_HBA1C>=5.7,1,0)
mcps_r$predm_ice<- ifelse(mcps_r$BASE_HBA1C>=6,1,0)
mcps_r$r_predm<-ifelse(mcps_r$R_HBA1C>=5.7 & mcps_r$diab_inc_fin==0,1,0)
mcps_r$years<-mcps_r$R_HXDIAB_YR-mcps_r$YEAR_RECRUITED
mcps_r$years[is.na(mcps_r$years)]<-mcps_r$YEAR_RESURVEY[is.na(mcps_r$years)]-mcps_r$YEAR_RECRUITED[is.na(mcps_r$years)]
mcps_r$age_onset<-mcps_r$years+mcps_r$AGE
mcps_r$status_predm<-mcps_r$predm+2*mcps_r$r_predm
mcps_r$status_predm[mcps_r$diab_inc_fin==1]<-5
mcps_r$status_predm<-factor(mcps_r$status_predm, labels = c("No-PreDM", "Regressed to No-PreDM", "Progressed to Pre-DM", "Stayed as Pre-DM", "Progressed to Diabetes"))
mcps_r$status<-mcps_r$r_predm
mcps_r$status[mcps_r$diab_inc_fin==1]<-2
mcps_r$status<-factor(mcps_r$status, labels = c("No-PreDM","Pre-DM","Diabetes"))
mcps_r$predm<-factor(mcps_r$predm, labels = c("No-PreDM","Pre-DM"))
table(mcps_r$R_HXDIAB) %>% prop.table()
table(mcps_r$status_predm, mcps_r$predm)

#### Diabetes progression rates ###
## All-cause mortality, Overall ##
no_predm_overall<- mcps_r %>%
  group_by(diab_inc_fin) %>%
  summarise(cases=n(), time=sum(years)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_overall$time)
cases<- no_predm_overall$cases[no_predm_overall$diab_inc_fin==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(qnorm(.975)*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, Overall ##
predm<- mcps_r %>% filter(predm=="Pre-DM") %>%
  group_by(diab_inc_fin) %>%
  summarise(cases=n(), time=sum(years)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm$time)
cases<- predm$cases[predm$diab_inc_fin==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(qnorm(.975)*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

## All-cause mortality, Overall ##
no_predm<- mcps_r %>% filter(predm=="No-PreDM") %>%
  group_by(diab_inc_fin) %>%
  summarise(cases=n(), time=sum(years)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm$time)
cases<- no_predm$cases[no_predm$diab_inc_fin==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(qnorm(.975)*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

#### Resurvey population table ###
supptab2<-mcps_r %>% 
  dplyr::select(hba1c_categories,AGE,female,COYOACAN,smoke,alcohol,PHYSGP,BASE_HYPERTENSION,BASE_STROKE,
                BASE_CKD,BASE_HEARTATTACK,BMI,ICE,SBP1,DBP1)%>%
  tbl_summary(by = hba1c_categories,missing = "ifany")%>%
  bold_labels()%>%
  add_p()%>%
  add_overall()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>% 
  align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = supptab2, split = TRUE) %>%
  body_end_section_landscape() %>% print(target = "Proyectos/Prediabetes/supptabla2.docx")

#### Transitions of prediabetes ###
mcps_r$predm2<-ifelse(mcps_r$BASE_HBA1C>=5.7, 1, 0)
m2<-coxph(Surv(years, diab_inc_fin)~predm+scale(AGE)+MALE+scale(BMI)+PHYSGP+alcohol+smoke+BASE_HYPERTENSION+scale(ICE)+frailty(COYOACAN), data=mcps_r)
summary(m2)

t1<-timeROC(T = mcps_r$years, 
        delta = mcps_r$diab_inc_fin,
        cause = 1,weighting = "cox",
        marker = mcps_r$predm_ice,
        other_markers = as.matrix(mcps_r$AGE, mcps_r$Sex, mcps_r$BMI, mcps_r$ICE),
        times = seq(2,15, 1))

time1<-data.frame(times=t1$times, auc=t1$AUC, criteria="IEC")

t2<-timeROC(T = mcps_r$years, 
        delta = mcps_r$diab_inc_fin,
        cause = 1,weighting = "cox",
        marker = mcps_r$predm2,
        other_markers = as.matrix(mcps_r$AGE, mcps_r$Sex, mcps_r$BMI, mcps_r$ICE),
        times = seq(2,15,1))

time2<-data.frame(times=t2$times, auc=t2$AUC, criteria="ADA")
time<-rbind(time1, time2)

f4c<-time %>%
  ggplot(aes(x=times, y=auc, color=criteria))+
  geom_point()+geom_line()+theme_bw()+
  ylab("Time-dependent AUROC")+xlab("Years of follow-up")+
  ylim(0.5, 0.7)+scale_x_continuous(breaks = c(2,4,6,8,10,12,14))+
  labs(col="Prediabetes criteria")+theme(legend.position = "top")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A"))

m2<-coxph(Surv(years, diab_inc_fin)~predm_ice+scale(AGE)+MALE+scale(BMI)+PHYSGP+alcohol+smoke+BASE_HYPERTENSION+scale(ICE)+frailty(COYOACAN), data=mcps_r)
summary(m2)

mcps_al <- mcps_r %>% filter(!is.na(status_predm)) %>%
  make_long(predm, status) 
mcps_al$x<-factor(mcps_al$x, labels = c("Baseline", "Resurvey")) 
table(mcps_r$predm, mcps_r$status_predm, useNA = "ifany")

f4a<-ggplot(mcps_al, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_label(size = 3, color = "white", fill = "gray40") +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_manual(values = wes_palette("Zissou1", 3))+
  theme(text = element_text(size = 12))
  
#### Regression to normoglycemia ###
mcps_predm<-mcps_r %>% filter(predm=="Pre-DM") %>% filter(!is.na(r_predm))
mcps_predm$status<-mcps_predm$r_predm+mcps_predm$diab_inc_fin*2
mcps_predm$status[mcps_predm$status==1]<-3
mcps_predm$status[mcps_predm$status==0]<-1
mcps_predm$status[mcps_predm$status==3]<-0
mcps_predm$diab <- NULL
mcps_predm$diab<-factor(mcps_predm$status, labels = c("Censored", "Regression", "Diabetes"))
mcps_predm$diab2<-NULL
mcps_predm$diab2<-mcps_predm$status;mcps_predm$diab2[mcps_predm$diab2==1]<-3;mcps_predm$diab2[mcps_predm$diab2==2]<-1;mcps_predm$diab2[mcps_predm$diab2==3]<-2
mcps_predm$diab2<-factor(mcps_predm$diab2, labels = c("Censored", "Diabetes","Regression"))
table(mcps_predm$diab2)

fgdata_cv <- finegray(Surv(years, diab) ~ ., data=mcps_predm, na.action=na.pass)
m1<-coxph(Surv(fgstart, fgstop, fgstatus)~scale(AGE)+MALE+scale(BMI)+PHYSGP+alcohol+smoke+BASE_HYPERTENSION+scale(ICE)+frailty(COYOACAN), weight=fgwt, data=fgdata_cv)
fgdata_cv2 <- finegray(Surv(years, diab2) ~ ., data=mcps_predm, na.action=na.pass)
m2<-coxph(Surv(fgstart, fgstop, fgstatus)~scale(AGE)+MALE+scale(BMI)+PHYSGP+alcohol+smoke+BASE_HYPERTENSION+scale(ICE)+frailty(COYOACAN), weight=fgwt, data=fgdata_cv2)
summary(m2)

f4b<-jtools::plot_summs(m2, m1,exp=T, model.names = c("Progression", "Regression"),
                        coefs = c("Age" = "scale(AGE)", "BMI" = "scale(BMI)", "WHtR"="scale(ICE)",
                                  "Male sex"="MALE", "Hypertension"="BASE_HYPERTENSION"), 
                        colors = c("#3B9AB2", "#EBCC2A"),legend.title = c("Transition"))+
  xlab("sHR for regression or progression from prediabetes")+theme_bw()+ylab("")+
  theme(legend.position = "top",text = element_text(size = 12))+
  theme(axis.ticks.y = element_blank())

fig3_1<-ggarrange(f4a, f4c, labels = c("A", "B"), widths = c(0.8,1))
fig3<-ggarrange(fig3_1, f4b, labels = c("", "C"), ncol=1)


ggsave(fig3, file="Proyectos/Prediabetes/Figure3.png", bg="transparent",
       width=20, height=20, units=c("cm"), dpi=600, limitsize = FALSE)



g1<-gg_miss_var(mcps_fin[c(1:37, 74:76, 105:133, 138, 682)], show_pct = TRUE, facet = hba1c_categories)
ggsave(g1, file="Proyectos/Prediabetes/Missing_Figure.jpg", bg="transparent",
       width=30, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

2250*2
1850*2+950
