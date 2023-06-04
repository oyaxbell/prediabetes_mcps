# Prediabetes as a risk factor for metabolic disease and cause-specific mortality in 125,442 diabetes-free adults in Mexico City
# Data Analysis: Omar Yaxmehen Bello-Chavolla
# Latest version of Analysis March 2023
# For any question regarding analysis contact Omar Yaxmehen Bello-Chavolla at oyaxbell@yahoo.com.mx

#### Library ####
pacman::p_load(haven, tidyverse, ggpubr, lmtest, nortest, gtools, data.table, caret, glmnet, survival, flextable, blandr, BlandAltmanLeh, corrplot,ggalluvial,coxme,flextable, TimeVTree,ggbreak,simPH,jtools,
               rms, bestNormalize, flexsurv, pROC, timeROC, factoextra, gridExtra,  nhanesA, wesanderson,forestmodel, ggedit,dummy,survivalROC,circlize,forestplot,gtsummary,shadowtext,PropCIs,Epi, popEpi,khroma,
               FactoMineR, fpc, NbClust, ggimage, glmnet, ggsci, survminer, cluster, ggplotify, UpSetR, nortest, viridis, officer, magrittr, coxphw, chorddiag, ggsankey, patchwork, survivalROC, lubridate, naniar)

#### Dataset magament ####
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/MCPS Projects")
mcps<-read.csv("Data/2021-004 MCPS BASELINE.csv") %>% left_join(read.csv("Data/2022-012 MCPS MORTALITY.csv"), by="PATID") %>%
  left_join(read.csv("Data/2022-012 MCPS RESURVEY.csv"), by="PATID") %>% 
  left_join(read.csv("Data/2022-012 MCPS BASELINE_NMR.csv"), by="PATID") 
mcps_fin<-mcps %>% filter(BASE_HBA1C<6.5) %>% filter(BASE_DIABETES==0) %>% 
  filter(DRUG_D1==0) %>% filter(DRUG_D2==0) %>% filter(DRUG_D3==0) %>% filter(DRUG_D4==0)
mcps_sens<-mcps_fin #%>%filter(!is.na(D000))
mcps_fin<-mcps_sens %>% filter(BASE_EMPHYSEMA==0 &
                                 BASE_HEARTATTACK==0 &
                                 BASE_ANGINA==0 &
                                 BASE_ASTHMA==0 &
                                 BASE_BREASTCANCER==0 &
                                 BASE_STROKE==0 &
                                 BASE_CKD==0 &
                                 BASE_CIRR==0 &
                                 BASE_HYPERTENSION==0 &
                                 BASE_LUNGCANCER==0 &
                                 BASE_PROSTATECANCER==0 &
                                 BASE_OTHCANCER==0 &
                                 BASE_STOMCANCER==0 &
                                 BASE_ORALCANCER==0 &
                                 BASE_PAD==0) %>% filter(!is.na(D000))
### Recodify ###
mcps$predm<- ifelse(mcps$BASE_HBA1C>=5.7,1,0)
mcps$predm_ice<- ifelse(mcps$BASE_HBA1C>=6,1,0)
mcps_fin$Sex<-ifelse(mcps_fin$MALE==1, 0, 1)
mcps_fin$predm<- ifelse(mcps_fin$BASE_HBA1C>=5.7,1,0)
mcps_fin$predm_ice<- ifelse(mcps_fin$BASE_HBA1C>=6,1,0)

mcps_fin$BMI<-mcps_fin$WEIGHT/((mcps_fin$HEIGHT/100)^2)
mcps_fin$ICE<-mcps_fin$WAISTC/(mcps_fin$HEIGHT)
mcps_fin$R_HXDIAB[is.na(mcps_fin$R_HXDIAB)]<-0
mcps_fin$hba1c_categories<-ifelse(mcps_fin$BASE_HBA1C<5.7, 0, 1)
mcps_fin$hba1c_categories[mcps_fin$BASE_HBA1C>=6]<-2
mcps_fin$hba1c_categories<-factor(mcps_fin$hba1c_categories, labels = c("<5.7%", "5.7-5.9%", "≥6.0%"))
mcps_fin$hba1c_cat<-mcps_fin$hba1c_categories
mcps_fin$smoke<-mcps_fin$SMOKEGP
mcps_fin$smoke[mcps_fin$SMOKEGP %in% c(3:5)]<-3
mcps_fin$alcohol<-mcps_fin$ALCGP
mcps_fin$alcohol[mcps_fin$ALCGP %in% c(3:5)]<-3
mcps_fin$ASCVD<-ifelse(mcps_fin$D003==1 | mcps_fin$D008==1, 1, 0)
mcps_fin$edad60<-ifelse(mcps_fin$AGE>=60, 1, 0)
mcps_fin$edad60<-factor(mcps_fin$edad60, labels = c("<60 years", "≥60 years"))
s_anthropoage<-function(x){
  # Data frame including Age (years), Sex (coded as "Men" and "Women"), Height (cm), Weight (kg),
  # Waist (cm) and Ethnicity (coded as "White", "Black", "Mexican-American" and "Others")
  x<-as.data.frame(x)
  imc<-x$Weight/(((x$Height)/100)^2); x$tr_imc <- log(imc)
  ice<-x$Waist/x$Height; x$tr_ice <- ice**(1/3)
  ## Load models ##
  model7<-load("Proyectos/Models/7_SAnthropo_M.rda")
  model8<-load("Proyectos/Models/8_SAnthropo_F.rda")
  model9<-load("Proyectos/Models/9_SAge_M.rda")
  model10<-load("Proyectos/Models/10_SAge_F.rda")
  ##Estimate model coefficients based on NHANES models##
  sM1<-1/((exp(coef(gomp1bM1)[1]*120)-1)/((coef(gomp1bM1)[1])))
  b0M1<-coef(gomp1bM1)[2]
  b1M1<-coef(gomp1bM1)[3]
  sW1<-1/((exp(coef(gomp1bF1)[1]*120)-1)/((coef(gomp1bF1)[1])))
  b0W1<-coef(gomp1bF1)[2]
  b1W1<-coef(gomp1bF1)[3]
  ## Loop for S-AnthropoAge estimation ##
  for(i in 1:nrow(x)){
    if(x$Sex[i]=="Women"){
      p1F<-predict(gomp1aF1, newdata = x,type="survival", ci=F, times = c(120))
      pred<-as.numeric(1-p1F$.pred)
      output<-(log(-sW1*log(1-pred))-b0W1)/b1W1}
    else if(x$Sex[i]=="Men") {
      p1M<-predict(gomp1aM1, newdata = x,type="survival", ci=F, times = c(120))
      pred<-as.numeric(1-p1M$.pred)
      output<-(log(-sM1*log(1-pred))-b0M1)/b1M1}
    ##Output
    return(output)}}

#### Recode dataset ####
#Age at risk
mcps_fin$age_risk<-mcps_fin$AGE+mcps_fin$PERSON_YEARS
  
#Premature mortality
mcps_fin$D000_p<-ifelse(mcps_fin$D000==1 & mcps_fin$age_risk<75,1,0)

#Alcohol Intake
mcps_fin$alcohol<-factor(mcps_fin$alcohol,labels = c("Never","Former","Current"))

#Smoking
mcps_fin$smoke<-factor(mcps_fin$smoke,labels = c("Never", "Former","Current"))

#Sex
mcps_fin$female<-ifelse(mcps_fin$MALE==0, 1, 0)

#Physical Activity
mcps_fin$PHYSGP<-factor(mcps_fin$PHYSGP,labels = c("None", "≤2 times a week","≥3 times a week"))

#Labels
setattr(mcps_fin$AGE, "label", "Age, (Years)")
setattr(mcps_fin$female, "label", "Female sex, (%)")
setattr(mcps_fin$COYOACAN, "label", "Residence in Coyoac?n, (%)")
setattr(mcps_fin$EDU_LEVEL, "label", "Educationnal Attainments, (%)")
setattr(mcps_fin$INCOME, "label", "Income, (pesos/month)")
setattr(mcps_fin$smoke, "label", "Smoking, (%)")
setattr(mcps_fin$alcohol, "label", "Alcohol intake, (%)")
setattr(mcps_fin$PHYSGP, "label", "Leisure-time physical activity, (%)")
setattr(mcps_fin$BASE_DIABETES, "label", "Diabetes, (%)")
setattr(mcps_fin$BASE_HYPERTENSION, "label", "Hypertension, (%)")
setattr(mcps_fin$BASE_STROKE, "label", "Stroke, (%)")
setattr(mcps_fin$BASE_CKD, "label", "CKD, (%)")
setattr(mcps_fin$BASE_HEARTATTACK, "label", "Myocardial infarction, (%)")
setattr(mcps_fin$BMI, "label", "Body Mass Index, (kg/m2)")
setattr(mcps_fin$ICE, "label", "Waist-to-Height Ratio, (%)")
setattr(mcps_fin$HIPC, "label", "Hip Circunference, (cm)")
setattr(mcps_fin$SBP1, "label", "Systolic Blood Pressure, (mmHg)")
setattr(mcps_fin$DBP1, "label", "Diastolic Blood Pressure, (mmHg)")
setattr(mcps_fin$BASE_HBA1C, "label", "HbA1c, (%)")
setattr(mcps_fin$predm, "label", "HbA1c ≥5.7%")
setattr(mcps_fin$predm_ice, "label", "HbA1c ≥6.0%")
setattr(mcps_fin$edad60, "label", "Age ≥60")

#### Missing data ####

g1<-gg_miss_var(mcps_fin[,-c(37:73, 77:104,134:137, 139:684)],show_pct = TRUE, facet = hba1c_categories)

ggsave(g1, file="Proyectos/Prediabetes/Missing_Figure.jpg", bg="transparent",
       width=30, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

#### Lexis expansion ####
mcps_fin<- mcps_fin %>% 
  mutate(date_recruited = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))

mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                   "age" = AGE,
                   "time_at_entry" = 0),
      exit = list("period" = date_recruited + PERSON_YEARS),
      exit.status = D000,
      data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T, 1, 0)
#mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0
#mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$PERSON_YEARS<5]<-0

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

table(mcps_fin2$premature)

#### Supplementary Table 1 ####
tab1<-mcps_fin %>% 
  dplyr::select(hba1c_categories,AGE,female,COYOACAN,smoke,alcohol,PHYSGP,BASE_HYPERTENSION,BASE_STROKE,
                BASE_CKD,BASE_HEARTATTACK,BMI,ICE,SBP1,DBP1)%>%
  tbl_summary(by = hba1c_categories,missing = "ifany")%>%
  bold_labels()%>%
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

doc <- read_docx() %>% body_add_flextable(value = tab1, split = TRUE) %>%
  body_end_section_landscape() %>% print(target = "Proyectos/Prediabetes/tabla1.docx")


#### Figure 1 Prevalence of prediabetes ####
prev_ADA<- mcps %>% 
  mutate(category=cut(AGE, breaks=c(-Inf, 40, 45, 50, 55, 60, 65, 70, 75,Inf), 
                      labels=c("<40","40-44","45-49","50-54", "55-59", "60-64", "65-69", "70-74", "≥75"))) %>%
  group_by(category, MALE) %>% summarise(n=n(), predm=sum(predm, na.rm = T)) %>%
  mutate(prev=(predm/n)*100) %>% 
  mutate(Sex=factor(MALE, labels=c("Female", "Male")))

temp<-NULL;L<-NULL;S<-NULL
for(i in 1:nrow(prev_ADA)) {
  temp[[i]]<-exactci(prev_ADA$predm[[i]],prev_ADA$n[[i]],conf.level=0.95)
  L[i]<-temp[[i]]$conf.int[1]
  S[i]<-temp[[i]]$conf.int[2]
}
prev_ADA$lower<-L*100
prev_ADA$upper<-S*100

prevalence_ADA<- prev_ADA %>%
  ggplot(aes(y=prev, x=category, col=Sex, group=Sex))+
  geom_point(size=2)+geom_line(size=1)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+
  xlab("Age groups")+
  ylab("ADA-defined prediabetes prevalence (%)")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A"))+
  theme_bw()+ylim(5, 65)

prev_iec<- mcps %>% 
  mutate(category=cut(AGE, breaks=c(-Inf, 40, 45, 50, 55, 60, 65, 70, 75,Inf), 
                      labels=c("<40","40-44","45-49","50-54", "55-59", "60-64", "65-69", "70-74", "≥75"))) %>%
  group_by(category, MALE) %>% summarise(n=n(), predm=sum(predm_ice, na.rm = T)) %>%
  mutate(prev=(predm/n)*100) %>% mutate(Sex=factor(MALE, labels=c("Female", "Male"))) 

temp<-NULL;L<-NULL;S<-NULL
for(i in 1:nrow(prev_iec)) {
  temp[[i]]<-exactci(prev_iec$predm[[i]],prev_iec$n[[i]],conf.level=0.95)
  L[i]<-temp[[i]]$conf.int[1]
  S[i]<-temp[[i]]$conf.int[2]
}
prev_iec$lower<-L*100
prev_iec$upper<-S*100

prevalence_iec<-prev_iec%>%
  ggplot(aes(y=prev, x=category, col=Sex, group=Sex))+
  geom_point(size=2)+geom_line(size=1)+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2)+
  xlab("Age groups (years)")+ylab("IEC-defined prediabetes prevalence (%)")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A"))+
  theme_bw()+ylim(5, 65)

fig1<-ggarrange(prevalence_ADA, prevalence_iec, labels = c("A", "B"), common.legend = T)

ggsave(fig1, file="Proyectos/Prediabetes/Figure1.jpg", bg="transparent",
         width=25, height=15, units=c("cm"), dpi=600, limitsize = FALSE)

#### Figure 2A - Comparison of death rates by cause ####
## All-cause mortality, Overall ##
no_predm_overall<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D000) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_overall$time)
cases<- no_predm_overall$cases[no_predm_overall$D000==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Overall"
no_predm_overall<-data.frame(lambda, LI, LS, pred, cause)

predm1_overall<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D000) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_overall$time)
cases<- predm1_overall$cases[predm1_overall$D000==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Overall"
predm1_overall<-data.frame(lambda, LI, LS, pred, cause)

predm2_overall<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D000) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_overall$time)
cases<- predm2_overall$cases[predm2_overall$D000==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Overall"
predm2_overall<-data.frame(lambda, LI, LS, pred, cause)

## All-cause mortality, Premature ##
no_predm_premature<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D000_p) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_premature$time)
cases<- no_predm_premature$cases[no_predm_premature$D000_p==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Premature"
no_predm_premature<-data.frame(lambda, LI, LS, pred, cause)

predm1_premature<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D000_p) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_premature$time)
cases<- predm1_premature$cases[predm1_premature$D000_p==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Premature"
predm1_premature<-data.frame(lambda, LI, LS, pred, cause)

predm2_premature<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D000_p) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_premature$time)
cases<- predm2_premature$cases[predm2_premature$D000_p==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Premature"
predm2_premature<-data.frame(lambda, LI, LS, pred, cause)

## Cardiovascular mortality mortality, Overall ##
no_predm_cv<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D003) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_cv$time)
cases<- no_predm_cv$cases[no_predm_cv$D003==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Cardiovascular"
no_predm_cv<-data.frame(lambda, LI, LS, pred, cause)

predm1_cv<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D003) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_cv$time)
cases<- predm1_cv$cases[predm1_cv$D003==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Cardiovascular"
predm1_cv<-data.frame(lambda, LI, LS, pred, cause)

predm2_cv<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D003) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_cv$time)
cases<- predm2_cv$cases[predm2_cv$D003==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Cardiovascular"
predm2_cv<-data.frame(lambda, LI, LS, pred, cause)

## Diabetes-related mortality cause mortality, overall ##
no_predm_db<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D023) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_db$time)
cases<- no_predm_db$cases[no_predm_db$D023==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Diabetes-related"
no_predm_db<-data.frame(lambda, LI, LS, pred, cause)

predm1_db<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D023) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_db$time)
cases<- predm1_db$cases[predm1_db$D023==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Diabetes-related"
predm1_db<-data.frame(lambda, LI, LS, pred, cause)

predm2_db<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D023) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_db$time)
cases<- predm2_db$cases[predm2_db$D023==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Diabetes-related"
predm2_db<-data.frame(lambda, LI, LS, pred, cause)

## Kidney-related mortality, kidney ##
no_predm_kidney<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D019) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_kidney$time)
cases<- no_predm_kidney$cases[no_predm_kidney$D019==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Kidney-related"
no_predm_kidney<-data.frame(lambda, LI, LS, pred, cause)

predm1_kidney<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D019) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_kidney$time)
cases<- predm1_kidney$cases[predm1_kidney$D019==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Kidney-related"
predm1_kidney<-data.frame(lambda, LI, LS, pred, cause)

predm2_kidney<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D019) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_kidney$time)
cases<- predm2_kidney$cases[predm2_kidney$D019==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Kidney-related"
predm2_kidney<-data.frame(lambda, LI, LS, pred, cause)

## Cancer-related cause mortality, Overall ##
no_predm_ca<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D040) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_ca$time)
cases<- no_predm_ca$cases[no_predm_ca$D040==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Cancer-related"
no_predm_ca<-data.frame(lambda, LI, LS, pred, cause)

predm1_ca<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D040) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_ca$time)
cases<- predm1_ca$cases[predm1_ca$D040==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Cancer-related"
predm1_ca<-data.frame(lambda, LI, LS, pred, cause)

predm2_ca<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D040) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_ca$time)
cases<- predm2_ca$cases[predm2_ca$D040==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Cancer-related"
predm2_ca<-data.frame(lambda, LI, LS, pred, cause)

## Stroke mortality, st ##
no_predm_st<- mcps_fin %>%
  filter(hba1c_cat=="<5.7%") %>%
  group_by(D008) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_st$time)
cases<- no_predm_st$cases[no_predm_st$D008==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"<5.7%"
cause<-"Stroke"
no_predm_st<-data.frame(lambda, LI, LS, pred, cause)

predm1_st<- mcps_fin %>%
  filter(hba1c_cat=="5.7-5.9%") %>%
  group_by(D008) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm1_st$time)
cases<- predm1_st$cases[predm1_st$D008==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"5.7-5.9%"
cause<-"Stroke"
predm1_st<-data.frame(lambda, LI, LS, pred, cause)

predm2_st<- mcps_fin %>%
  filter(hba1c_cat=="≥6.0%") %>%
  group_by(D008) %>%
  summarise(cases=n(), time=sum(PERSON_YEARS)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(predm2_st$time)
cases<- predm2_st$cases[predm2_st$D008==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
pred<-"≥6.0%"
cause<-"Stroke"
predm2_st<-data.frame(lambda, LI, LS, pred, cause)

death_rates<-rbind(no_predm_overall, no_predm_premature, no_predm_cv, no_predm_db, no_predm_kidney, no_predm_ca, no_predm_st,
                   predm1_overall, predm1_premature,predm1_cv, predm1_db, predm1_kidney, predm1_ca, predm1_st,
                   predm2_overall, predm2_premature, predm2_cv, predm2_db, predm2_kidney, predm2_ca, predm2_st)

death_rates<-death_rates %>%
  mutate(across(cause, factor, levels=c("Diabetes-related","Stroke","Kidney-related","Cancer-related","Cardiovascular","Premature","Overall"))) %>%
  mutate(across(pred, factor, levels=c("<5.7%","5.7-5.9%", "≥6.0%"))) %>%
  arrange(desc(lambda)) 

n3<-death_rates %>%
  ggplot(aes(y=lambda, fill=pred, x=cause)) + 
  geom_bar(stat="identity", position='dodge', color="black", width=0.8)+
  theme_bw()+ylab("Frequency (%)")+xlab("")+
  theme(axis.ticks.x = element_blank())+
  labs(fill="HbA1c")+
  scale_fill_manual(values = wes_palette("Zissou1", 3))+
  ylab("Deaths per 1,000 person-years")+ 
  scale_y_break(c(5.2,7.8))+ 
  geom_errorbar(aes(ymin=LI, ymax=LS), width=.2, position=position_dodge(.8))+
  coord_flip()+ 
  theme(text = element_text(size = 15), legend.position = "top")

#### Analysis 1 - Prediabetes cut-offs and mortality ####

m0<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C+(COYOACAN)+MALE+strata(age.cat)+cluster(id),data=mcps_fin2)
m1<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m3<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

m0_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m1_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m2_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m3_adj1<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

m0_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~BASE_HBA1C+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m1_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~hba1c_cat+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m2_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
m3_adj2<-coxph(Surv(time_at_entry, time_at_exit, premature)~predm_ice+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

#### Table 2 - All-cause mortality regarding HbA1c definitions ####
#Flex Table Merge
reg.0<-m0 %>%
  tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"),add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.1<-m1 %>%
  tbl_regression(exponentiate = TRUE, include = c("hba1c_cat"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.2<-m2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.3<-m3 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.0_adj1<-m0_adj1 %>%
  tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels

reg.1_adj1<-m1_adj1 %>%
  tbl_regression(exponentiate = TRUE, include = c("hba1c_cat"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.2_adj1<-m2_adj1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.3_adj1<-m3_adj1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.0_adj2<-m0_adj2 %>%
  tbl_regression(exponentiate = TRUE, include = c("BASE_HBA1C"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels

reg.1_adj2<-m1_adj2 %>%
  tbl_regression(exponentiate = TRUE, include = c("hba1c_cat"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.2_adj2<-m2_adj2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

reg.3_adj2<-m3_adj2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

row1 <- tbl_merge(list(reg.0, reg.0_adj1,reg.0_adj2), tab_spanner = c("Model 1", "Model 2","Model 3"))
row2 <- tbl_merge(list(reg.1, reg.1_adj1,reg.1_adj2), tab_spanner = c("Model 1", "Model 2","Model 3"))
row3 <- tbl_merge(list(reg.2, reg.2_adj1,reg.2_adj2), tab_spanner = c("Model 1", "Model 2","Model 3"))
row4 <- tbl_merge(list(reg.3, reg.3_adj1,reg.3_adj2), tab_spanner = c("Model 1", "Model 2","Model 3"))

tab2 <-tbl_stack(list(row1, row2, row3, row4), group_header = c("Continuous", "Categorical", "ADA definition", "IEC definition"))%>%
  as_flex_table() %>%
  align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>% 
  body_end_section_landscape() %>%
  print(target = "Proyectos/Prediabetes/tabla2.docx")

#### Sensitivity analysis 1: Individuals with comorbidities ####

mcps_sens<- mcps_sens %>% filter(!is.na(D000))

### Recodify ###
mcps_sens$Sex<-ifelse(mcps_sens$MALE==1, 0, 1)
mcps_sens$predm<- ifelse(mcps_sens$BASE_HBA1C>=5.7,1,0)
mcps_sens$predm_ice<- ifelse(mcps_sens$BASE_HBA1C>=6,1,0)

mcps_sens$BMI<-mcps_sens$WEIGHT/((mcps_sens$HEIGHT/100)^2)
mcps_sens$ICE<-mcps_sens$WAISTC/(mcps_sens$HEIGHT)
mcps_sens$R_HXDIAB[is.na(mcps_sens$R_HXDIAB)]<-0
mcps_sens$hba1c_categories<-ifelse(mcps_sens$BASE_HBA1C<5.7, 0, 1)
mcps_sens$hba1c_categories[mcps_sens$BASE_HBA1C>=6]<-2
mcps_sens$hba1c_categories<-factor(mcps_sens$hba1c_categories, labels = c("<5.7%", "5.7-5.9%", "≥6.0%"))
mcps_sens$hba1c_cat<-mcps_sens$hba1c_categories
mcps_sens$smoke<-mcps_sens$SMOKEGP
mcps_sens$smoke[mcps_sens$SMOKEGP %in% c(3:5)]<-3
mcps_sens$alcohol<-mcps_sens$ALCGP
mcps_sens$alcohol[mcps_sens$ALCGP %in% c(3:5)]<-3
mcps_sens$YEAR_DOB2<-ym(paste0(mcps_sens$YEAR_RECRUITED-mcps_sens$AGE, mcps_sens$MONTH_RECRUITED))
mcps_sens$YEAR_DOB2<-decimal_date(mcps_sens$YEAR_DOB2)
mcps_sens$n_comorb<-mcps_sens$BASE_EMPHYSEMA+mcps_sens$BASE_HEARTATTACK+mcps_sens$BASE_ANGINA+
  mcps_sens$BASE_ASTHMA+mcps_sens$BASE_STROKE+mcps_sens$BASE_CKD+mcps_sens$BASE_PEP+
  mcps_sens$BASE_CIRR+mcps_sens$BASE_HYPERTENSION+mcps_sens$BASE_DIABETES+mcps_sens$BASE_BREASTCANCER+
  mcps_sens$BASE_CERVCANCER+mcps_sens$BASE_LUNGCANCER+ mcps_sens$BASE_ORALCANCER+mcps_sens$BASE_OTHCANCER+
  mcps_sens$BASE_PAD+mcps_sens$BASE_PROSTATECANCER+mcps_sens$BASE_STOMCANCER

mcps_lexis <- Lexis(entry = list("period"=AGE+YEAR_DOB2, "age"=AGE), 
                    exit = list("age" = AGE+PERSON_YEARS), 
                    exit.status = D000, 
                    data = mcps_sens)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T, 1, 0)
mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0

m0<-coxph(Surv(lex.dur, premature)~BASE_HBA1C+(COYOACAN)+MALE+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m1<-coxph(Surv(lex.dur, premature)~hba1c_cat+MALE+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m2<-coxph(Surv(lex.dur, premature)~predm+MALE+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m3<-coxph(Surv(lex.dur, premature)~predm_ice+MALE+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)

m0_adj1<-coxph(Surv(lex.dur, premature)~BASE_HBA1C+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+MALE+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m1_adj1<-coxph(Surv(lex.dur, premature)~hba1c_cat+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m2_adj1<-coxph(Surv(lex.dur, premature)~predm+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m3_adj1<-coxph(Surv(lex.dur, premature)~predm_ice+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)

m0_adj2<-coxph(Surv(lex.dur, premature)~BASE_HBA1C+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+MALE+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m1_adj2<-coxph(Surv(lex.dur, premature)~hba1c_cat+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m2_adj2<-coxph(Surv(lex.dur, premature)~predm+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m3_adj2<-coxph(Surv(lex.dur, premature)~predm_ice+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+n_comorb+cluster(lex.id)+strata(age.cat), data=mcps_fin2)

#### Sensitivity analysis 2: Premature mortality ####
mcps_fin$YEAR_DOB2<-ym(paste0(mcps_fin$YEAR_RECRUITED-mcps_fin$AGE, mcps_fin$MONTH_RECRUITED))
mcps_fin$YEAR_DOB2<-decimal_date(mcps_fin$YEAR_DOB2)

mcps_lexis <- Lexis(entry = list("period"=AGE+YEAR_DOB2, "age"=AGE), 
                    exit = list("age" = AGE+PERSON_YEARS), 
                    exit.status = D000, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$ageout<-mcps_fin2$age+mcps_fin2$lex.dur
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T, 1, 0)
mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0

m0<-coxph(Surv(lex.dur, premature)~BASE_HBA1C+(COYOACAN)+MALE+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m1<-coxph(Surv(lex.dur, premature)~hba1c_cat+MALE+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m2<-coxph(Surv(lex.dur, premature)~predm+MALE+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m3<-coxph(Surv(lex.dur, premature)~predm_ice+MALE+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)

m0_adj1<-coxph(Surv(lex.dur, premature)~BASE_HBA1C+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+MALE+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m1_adj1<-coxph(Surv(lex.dur, premature)~hba1c_cat+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m2_adj1<-coxph(Surv(lex.dur, premature)~predm+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m3_adj1<-coxph(Surv(lex.dur, premature)~predm_ice+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)


m0_adj2<-coxph(Surv(lex.dur, premature)~BASE_HBA1C+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+MALE+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m1_adj2<-coxph(Surv(lex.dur, premature)~hba1c_cat+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m2_adj2<-coxph(Surv(lex.dur, premature)~predm+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)
m3_adj2<-coxph(Surv(lex.dur, premature)~predm_ice+BMI+scale(ICE)+EDU_LEVEL+factor(alcohol)+factor(smoke)+(MALE)+(COYOACAN)+cluster(lex.id)+strata(age.cat), data=mcps_fin2)

#### Analysis 2 - Cardiovascular mortality ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS),
                    exit.status = D003, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm_ice+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

#### Analysis 3 - Diabetes-related mortality ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS),
                    exit.status = D023, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

d1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
d2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm_ice+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

#### Analysis 4 - Kidney-related mortality ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), 
                    exit.status = D019, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

k1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
k2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm_ice+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

#### Analysis 5 - Cancer-related mortality ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS),
                    exit.status = D040, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

ca1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
ca2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm_ice+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

#### Analysis 6 - Stroke-related mortality ####
mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS),
                    exit.status = D008, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

st1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)
st2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~predm_ice+MALE+BMI+scale(ICE)+EDU_LEVEL+factor(PHYSGP)+factor(alcohol)+factor(smoke)+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2)

#### Table 2 - Cause-specific mortality regarding HbA1c definitions ####

#Flex Table Merge

c.1<-c1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

c.2<-c2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

d.1<-d1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

d.2<-d2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

k.1<-k1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

k.2<-k2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

ca.1<-ca1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

ca.2<-ca2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

st.1<-st1 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

st.2<-st2 %>%
  tbl_regression(exponentiate = TRUE, include = c("predm_ice"), add_estimate_to_reference_rows = T, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
  bold_labels() %>%
  italicize_levels()

row1 <- tbl_merge(list(c.1, c.2), tab_spanner = c("ADA definition", "IEC definition"))
row2 <- tbl_merge(list(d.1, d.2), tab_spanner = c("ADA definition", "IEC definition"))
row3 <- tbl_merge(list(k.1, k.2), tab_spanner = c("ADA definition", "IEC definition"))
row4 <- tbl_merge(list(ca.1, ca.2), tab_spanner = c("ADA definition", "IEC definition"))
row5 <- tbl_merge(list(st.1, st.2), tab_spanner = c("ADA definition", "IEC definition"))

tab3 <-tbl_stack(list(row1, row2, row3, row4, row5), group_header = c("Cardiovascular", "Diabetes-related", "Kidney-related", "Cancer-related", "Stroke"))%>%
  as_flex_table() %>%
  align(align = "center",part = "all") %>% 
  autofit()

doc <- read_docx() %>% body_add_flextable(value = tab3, split = TRUE) %>% 
  print(target = "Proyectos/Prediabetes/tabla3.docx")

#### Figure 2B - Interaction between Age*HbA1c for all-cause mortality ####

mcps_lexis <- Lexis(entry = list("period" = date_recruited,
                                 "age" = AGE,
                                 "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), 
                    exit.status = D000, 
                    data = mcps_fin)
mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")

mcps_fin2<- mcps_fin2 %>%
  rename(id = lex.id,
         fu_period = lex.dur,
         entry_status = lex.Cst,
         exit_status = lex.Xst) %>% 
  mutate(time_at_exit = time_at_entry + fu_period,
         age_at_exit = age + fu_period)

# Run basic model
m0.1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~hba1c_cat+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% filter(AGE>=65))
m0.2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~hba1c_cat+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% filter(between(AGE, 50, 65)))
m0.3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~hba1c_cat+MALE+(COYOACAN)+cluster(id)+strata(age.cat), data=mcps_fin2 %>% filter(AGE<50))

f2b<-jtools::plot_summs(m0.3, m0.2, m0.1, exp=T, model.names = c("<50", "50-64", "≥65"),
                   coefs = c("HbA1c 5.7-5.9%" = "hba1c_cat5.7-5.9%", "HbA1c ≥6.0%" = "hba1c_cat≥6.0%"),
                   legend.title = c("Age at baseline (yrs)"), colors = wes_palette("Zissou1", 4))+
  xlab("HR for all-cause mortality")+ylab("")+theme_bw()+
  theme(legend.position = "top",text = element_text(size = 15))+
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5), axis.ticks.y = element_blank())

#### Build Figure 2 ####

fig2<-n3+f2b+plot_layout(widths = c(1,0.8))+ 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))

ggsave(fig2, file="Proyectos/Prediabetes/Figure2.jpg", bg="transparent",
       width=35, height=15, units=c("cm"), dpi=600, limitsize = FALSE)

#### Resampling ####
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

#### Diabetes progression rates ####
## All-cause mortality, Overall ##
no_predm_overall<- mcps_r %>%
  group_by(diab_inc_fin) %>%
  summarise(cases=n(), time=sum(years)) %>% drop_na()
#Luego extraes el tiempo de cada caso
time<- sum(no_predm_overall$time)
cases<- no_predm_overall$cases[no_predm_overall$diab_inc_fin==1]
lambda<-(cases/time)*1000
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
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
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
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
lambda.se<-(sqrt(cases/time^2))*1000; LI<-lambda-(1.96*lambda.se); 
LS<-lambda+(qnorm(1-0.05/2)*lambda.se)
c(lambda, LI, LS)

#### Supplementary Table 2 ####
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

#### Transitions of prediabetes ####
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
  
#### Regression to normoglycemia ####
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


ggsave(fig3, file="Proyectos/Prediabetes/Figure3.jpg", bg="transparent",
       width=20, height=20, units=c("cm"), dpi=600, limitsize = FALSE)

  