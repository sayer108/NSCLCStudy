setwd("~/Desktop/NSCLCGithub")
####Packages needed for analysis
library(glmulti)
library(survcomp)
library(finalfit)
library(survminer)
library(caret)
library(pROC)
library(ROCR)
library(MLmetrics)

##Data sets used
#All patients combined into 1 file
all_data<- read.csv("combined.data.csv")
#Training data subset
train_data_whole<- read.csv("training.data.csv")
#testing data subset
test_data_whole<- read.csv("testing.data.csv")
#Machine learning results
OS_Figure_2E<- read.csv("Figure_3_OS.csv")
PFS_Figure_2F<- read.csv("Figure_3_PFS.csv")

#########Figure 1, Kaplan Meier Survival analysis plots
##Kaplan Meier survival curves evaluating the role of immune-related adverse events (irAE) 
##and history of tyrosine kinase inhibitor (TKI) therapy prior to immune checkpoint inhibitor (ICI) 
##therapy on overall survival (OS) and real-world progression free survival (rwPFS). 

#(A) and (B) compare patients experiencing an irAE to those that did not,
#evaluating OS and rwPFS, respectively.
#(C) and (D) compare patients experiencing an early474
#irAE to those experiencing a late-irAE event, evaluating OS and rwPFS, respectively
#(E) and 475 (F) compares patients who received a TKI prior to 
#ICI therapy to those that did not, evaluating OS and rwPFS, respectively.
#Figure 1, Panel A
OS_irAE_yn<-  ggsurvplot(
  survfit(Surv(I(
    os_days_update
  )/30, survival_flag) ~ irAE.1_pos_neg, all_data),
  pval = T,
  palette = "jama",
  pval.coord = c(0, 0.03),
  ylab = "OS Probability",
  xlab= "Time in Months",
  break.time.by=3,      
  risk.table = T,
  xlim= c(0,24),
  legend= "right",
  legend.title= "",
  legend.labs= c("No irAE","irAE"),
  fontsize=3,
  font.tickslab= c(10, "bold"),
  ggtheme= theme_classic2(base_size = 10),
  tables.y.text= F,
  tables.theme = theme_cleantable()
)
###Confidence Interval and P-value statistics
#p-value
coxph(Surv(os_days_update/30, survival_flag)~ irAE.1_pos_neg,data= all_data)
#confidence interval
survfit(Surv(I(os_days_update),  survival_flag) ~ irAE.1_pos_neg, all_data)

#Figure 1, Panel B
PFS_irAE_yn<-  
  ggsurvplot(
    survfit(Surv(I(
      pfs_days
    )/30,  pfs_flag) ~ irAE.1_pos_neg, all_data),
    pval = T,
    palette = "jama",
    pval.coord = c(0, 0.03),
    ylab = "rwPFS Probability",
    xlab= "Time in Months",
    break.time.by=3,      
    risk.table = T,
    xlim= c(0,24),
    legend= "none",
    fontsize=3,
    font.tickslab= c(10, "bold"),
    ggtheme= theme_classic2(base_size = 10),
    tables.y.text= F,
    tables.theme = theme_cleantable()
  )
##confidence interval
survfit(Surv(I(pfs_days/30),  pfs_flag) ~ irAE.1_pos_neg, data= all_data)
#p-value
coxph(Surv(pfs_days/30, pfs_flag)~ irAE.1_pos_neg,data= all_data)

#Figure 1, Panel C
OS_irAE_earlylate<- 
  ggsurvplot(
    survfit(Surv(I(
      os_days_update
    )/30, survival_flag) ~ irAE.1_EarlyLate, all_data),
    pval = T,
    palette = "jama",
    pval.coord = c(0, 0.03),
    ylab = "OS Probability",
    xlab= "Time in Months",
    break.time.by=3,      
    risk.table = T,
    xlim= c(0,24),
    legend= "right",
    legend.title= "",
    legend.labs= c("Early irAE","Late irAE"),
    fontsize=3,
    font.tickslab= c(10, "bold"),
    ggtheme= theme_classic2(base_size = 10),
    tables.y.text= F,
    tables.theme = theme_cleantable()
  )
#Confidence Interval
survfit(Surv(I(os_days_update/30),  survival_flag) ~ irAE.1_EarlyLate, data= all_data)
#P-value
coxph(Surv(os_days_update/30, survival_flag)~ irAE.1_EarlyLate, data= all_data)

#Figure 1, Panel D
PFS_irAE_earlylate<- 
  ggsurvplot(
    survfit(Surv(I(
      pfs_days
    )/30,  pfs_flag) ~ irAE.1_EarlyLate, all_data),
    pval = T,
    palette = "jama",
    pval.coord = c(0, 0.03),
    ylab = "rwPFS Probability",
    xlab= "Time in Months",
    break.time.by=3,      
    risk.table = T,
    xlim= c(0,24),
    legend= "none",
    fontsize=3,
    font.tickslab= c(10, "bold"),
    ggtheme= theme_classic2(base_size = 10),
    tables.y.text= F,
    tables.theme = theme_cleantable()
  )
#Confidence Interval
survfit(Surv(I(pfs_days/30),  pfs_flag) ~ irAE.1_EarlyLate, all_data)
#P-value
coxph(Surv(pfs_days/30, pfs_flag)~ irAE.1_EarlyLate,data= all_data)

#Figure 1, Panel E
OS_Prior_TKI<- 
  ggsurvplot(
    survfit(Surv(I(
      os_days_update
    )/30, survival_flag) ~ prior_tki, all_data),
    pval = T,
    palette = "jama",
    pval.coord = c(0, 0.03),
    ylab = "OS Probability",
    xlab= "Time in Months",
    break.time.by=3,      
    risk.table = T,
    xlim= c(0,24),
    legend= "right",
    legend.title= "",
    legend.labs= c("No Prior TKI","Prior TKI"),
    fontsize=3,
    font.tickslab= c(10, "bold"),
    ggtheme= theme_classic2(base_size = 10),
    tables.y.text= F,
    tables.theme = theme_cleantable()
  )
#Confidence Interval
survfit(Surv(I(os_days_update/30), survival_flag) ~ priortki, all_data)
#P-Value
coxph(Surv(os_days_update/30, survival_flag)~ priortki, data= all_data)

#Figure 1, Panel F
PFS_Prior_TKI<- 
  ggsurvplot(
    survfit(Surv(I(
      pfs_days
    )/30,  pfs_flag) ~ prior_tki, all_data),
    pval = T,
    palette = "jama",
    pval.coord = c(0, 0.03),
    ylab = "rwPFS Probability",
    xlab= "Time in Months",
    break.time.by=3,      
    risk.table = T,
    xlim= c(0,24),
    legend= "none",
    legend.title= "",
    legend.labs= c("No Prior TKI","Prior TKI"),
    fontsize=3,
    font.tickslab= c(10, "bold"),
    ggtheme= theme_classic2(base_size = 10),
    tables.y.text= F,
    tables.theme = theme_cleantable()
  )
#Confidence Interval
survfit(Surv(I(pfs_days/30),  pfs_flag) ~ priortki, data= all_data)
#P-Value
coxph(Surv(pfs_days/30, pfs_flag)~ priortki,data= all_data)

#####This creates the multi-panel figure utilized in the paper, attaching labeling and formatting
##A, OS irAE YN
A_Figure1<- OS_irAE_yn
A_Figure1$plot<- A_Figure1$plot + labs(tag= "A") + theme(plot.tag.position = "topleft")
##B, rwPFS irAE YN
B_Figure1<- PFS_irAE_yn
B_Figure1$plot<- B_Figure1$plot + labs(tag= "B") + theme(plot.tag.position = "topleft")
##C, OS earlylate irAE
C_Figure1<- OS_irAE_earlylate
C_Figure1$plot<- C_Figure1$plot + labs(tag= "C") + theme(plot.tag.position = "topleft")
##D, rwPFS earlylate irAE
D_Figure1<- PFS_irAE_earlylate
D_Figure1$plot<- D_Figure1$plot + labs(tag= "D") + theme(plot.tag.position = "topleft")
##E, OS Prior TKI
E_Figure1<- OS_Prior_TKI
E_Figure1$plot<- E_Figure1$plot + labs(tag= "E") + theme(plot.tag.position = "topleft")
##F, PFS Prior TKI
F_Figure1<- PFS_Prior_TKI
F_Figure1$plot<- F_Figure1$plot + labs(tag= "F") + theme(plot.tag.position = "topleft")
###This creates the final multi-panel figure utilized in the paper with the 6 plots
Combined_Figure1<- list(A_Figure1, C_Figure1, E_Figure1, B_Figure1, D_Figure1, F_Figure1)
Figure1_Plot<- arrange_ggsurvplots(Combined_Figure1, ncol = 2, nrow = 3)
ggsave("Figure1.png", plot= Figure1_Plot, units= "in", height= 10, width = 8 , dpi= 300)



####This code implements various logistic regression analyses utilized in the paper, with the results referenced
####in Table 2, Figure 2, and supplementary figure 3
####Table 2 Description
####Evaluation of model performance Metrics used to evaluate performance of various
####models implemented are accuracy, Area Under the Receiver Operator Curve (AUROC), Area
####Under Precision-Recall Curve (AUPRC) and F-1 score. All measures range between 0-1.

####Figure 2 Description
####Kaplan Meier survival curves comparing patients assigned high survival scores
####compared to low survival scores for baseline and optimized logistic regression models in
####addition to elastic-net logistic regression model for overall survival (OS) and real-world
####progression free survival (rwPFS)
####A and B, Prediction results of baseline logistic regression models (A and B) for OS and rwPFS
####C and D, Prediction results of optimized logistic regression models (C and D)
####E and F, Prediction results of elastic-net logistic regression model (E and F)
####Panels A, C and E are OS, with B, D, and F being rwPFS
####Patients with survival scores higher than the median cut-off had significantly longer OS and rwPFS 
####than those with lower survival scores. P-values less than 0.05 indicates a statistically significant 
####probability from the log rank statistical test that the high and low score groups have different 
####survival probability.

####Supplementary Table 3 Description
##Coefficients of baseline and optimized logistic regression models 
##predicting 1-year overall survival (OS) and 6-month real-world progression free survival (rwPFS). 
##Significant variables with P-values < 0.05 are noted 
##with * symbols in the table; variables excluded in optimized logistic regression models are left blank in the 
##table. Variables having p-values with 2 asterisks (**) have P-values < 0.01, and 3 asterisks (***) 
##have p-values < 0.001. 



####Formatted subset data for training and testing data sets
####OS_1yr_test_data is the training data set used to create the logistic regression model for 1 year OS
####OS_1yr_test_data is the testing data set used to make predictions of 1 year OS with the created model 
train_data_whole$y<- train_data_whole$os_1yr_flag_surv
test_data_whole$y<- test_data_whole$os_1yr_flag_surv
OS_1yr_train_data<- subset(train_data_whole, select= c(male, Age.at.diagnosis, Adenocarcinoma, Squamous, White, Asian, Black, Smoking, StageIorII, StageIII, StageIV,
                                                antipd1, antipdl1, priortki, irae_pos_1yr,irae_cat_1yr, early_irae, pdl1_neg, pdl1_149, pdl1_50100, y))
OS_1yr_test_data<- subset(test_data_whole, select= c(male, Age.at.diagnosis, Adenocarcinoma, Squamous, White, Asian, Black, Smoking, StageIorII, StageIII, StageIV,
                                              antipd1, antipdl1, priortki, irae_pos_1yr,irae_cat_1yr, early_irae, pdl1_neg, pdl1_149, pdl1_50100, y))
##Creation of baseline 1 year overall survival logistic regression model
##Baseline incorporates all variables in the data subset
Baseline_Model_1yr_Log<- glm(y~., data= OS_1yr_train_data)
##Baseline logistic regression model for 1 year overall survival metrics with coefficients and P-values
##Supplementary Figure 3 coefficients and P-values
summary(Baseline_Model_1yr_Log)
##Baseline 1 year overall survival logistic regression model predictions
OS_baseline_survival_prediction<- predict(Baseline_Model_1yr_Log, OS_1yr_test_data, type= "response")
##Baseline 1 year overall survival logistic regression model peformance metrics (Table 2)
confusionMatrix(data= as.factor(as.numeric(OS_baseline_survival_prediction>.5)), reference= as.factor(OS_1yr_test_data$y), mode= "everything", positive = "1")
confusionMatrix(data= as.factor(as.numeric(OS_baseline_survival_prediction>.5)), reference= as.factor(OS_1yr_test_data$y), mode= "prec_recall", positive = "1")
baseline_1yr_roc_model<- roc(OS_1yr_test_data$y, OS_baseline_survival_prediction)
PRAUC(OS_baseline_survival_prediction, OS_1yr_test_data$y)
auc(baseline_1yr_roc_model)

##Creation of optimized 1 year overall survival logistic regression model
##optimized selects fewer variables based on an AIC minimization strategy
Optimized_Model_1yr_Log<- step(Baseline_Model_1yr_Log, direction = c("both"))
##optimized logistic regression model for 1 year overall survival metrics with coefficients and P-values
##Supplentary Figure 3 coefficients and P-values
summary(Optimized_Model_1yr_Log)
##Optimized 1 year overall survival logistic regression model predictions
OS_optimize_survival_prediction<- predict(Optimized_Model_1yr_Log, OS_1yr_test_data, type= "response")
##Optimized 1 year overall survival logistic regression model peformance metrics (Table 2)
confusionMatrix(data= as.factor(as.numeric(OS_optimize_survival_prediction>.5)), reference= as.factor(OS_1yr_test_data$y), mode= "everything", positive = "1")
confusionMatrix(data= as.factor(as.numeric(OS_optimize_survival_prediction>.5)), reference= as.factor(OS_1yr_test_data$y), mode= "prec_recall", positive = "1")
optimized_1yr_roc_model<- roc(OS_1yr_test_data$y, OS_optimize_survival_prediction)
PRAUC(OS_optimize_survival_prediction, OS_1yr_test_data$y)
auc(optimized_1yr_roc_model)

#########Creation of 6month predictive models###############
####PFS_6mo_train_data is the training data set used to create the logistic regression model for 6 month rwPFS
####PFS_6mo_test_data is the testing data set used to make predictions of 6 month rwPFS with the created model 
train_data_whole$z<- train_data_whole$pfs_6months_flag_surv
test_data_whole$z<- test_data_whole$pfs_6months_flag_surv
PFS_6mo_train_data<- subset(train_data_whole, select= c(male, Age.at.diagnosis, Adenocarcinoma, Squamous, White, Asian, Black, Smoking, StageIorII, StageIII, StageIV,
                                                 antipd1, antipdl1, priortki, irae_pos_6mo,irae_cat_6mo, early_irae, pdl1_neg, pdl1_149, pdl1_50100, z))
PFS_6mo_test_data<- subset(test_data_whole, select= c(male, Age.at.diagnosis, Adenocarcinoma, Squamous, White, Asian, Black, Smoking, StageIorII, StageIII, StageIV,
                                               antipd1, antipdl1, priortki, irae_pos_6mo,irae_cat_6mo, early_irae, pdl1_neg, pdl1_149, pdl1_50100, z))
##creation of Baseline 6 month progression free survival logistic regression model
##Baseline incorporates all variables in the data subset
Baseline_Model_6mo_Log<- glm(z~., data= PFS_6mo_train_data)
##Baseline logistic regression model for 1 year overall survival metrics with coefficients and P-values
##Supplentary Figure 3 coefficients and P-values
summary(Baseline_Model_1yr_Log)
##Baseline 6 month progression free survival logistic regression model predictions
PFS_baseline_survival_prediction<- predict(Baseline_Model_6mo_Log, PFS_6mo_test_data, type= "response")
##Baseline 6 month progression free survival logistic regression model peformance metrics (Table 2)
confusionMatrix(data= as.factor(as.numeric(PFS_baseline_survival_prediction>.5)), reference= as.factor(PFS_6mo_test_data$z), mode= "everything", positive = "1")
confusionMatrix(data= as.factor(as.numeric(PFS_baseline_survival_prediction>.5)), reference= as.factor(PFS_6mo_test_data$z), mode= "prec_recall", positive = "1")
baseline_6mo_roc_model<- roc(PFS_6mo_test_data$z, PFS_baseline_survival_prediction)
PRAUC(PFS_baseline_survival_prediction, PFS_6mo_test_data$z)
auc(baseline_6mo_roc_model)

##Optimized 6 month progression free survival logistic regression model
##optimized selects fewer variables based on an AIC minimization strategy
Optimized_Model_6mo_Log<- step(Baseline_Model_6mo_Log, direction = c("both"))
##optimized logistic regression model for 6 month progression free survival metrics with coefficients and P-values
##Supplentary Figure 3 coefficients and P-values
summary(Optimized_Model_6mo_Log)
##Optimized 6 month progression free survival logistic regression model predictions
PFS_optimize_survival_prediction<- predict(Optimized_Model_6mo_Log, PFS_6mo_test_data, type= "response")
##Optimized 6 month progression free survival logistic regression model peformance metrics (Table 2)
confusionMatrix(data= as.factor(as.numeric(PFS_optimize_survival_prediction>.5)), reference= as.factor(PFS_6mo_test_data$z), mode= "everything", positive = "1")
confusionMatrix(data= as.factor(as.numeric(PFS_optimize_survival_prediction>.5)), reference= as.factor(PFS_6mo_test_data$z), mode= "prec_recall", positive = "1")
optimized_6mo_roc_model<- roc(PFS_6mo_test_data$z, PFS_optimize_survival_prediction)
PRAUC(PFS_optimize_survival_prediction, PFS_6mo_test_data$z)
auc(optimized_6mo_roc_model)



######Creation of Figure 2
##Panel A, Figure 2
##l year OS baseline logistic regression model high and low score assignment and kaplan meier curve
##Create and format data frame for analysis titled Baseline_1yr_risk, assigned survival scores from previous logistic regression analysis
Baseline_1yr_risk<- as.data.frame(cbind(OS_baseline_survival_prediction, test_data_whole$os_days_update, test_data_whole$VitalStatus))
colnames(Baseline_1yr_risk)<- c("Risk", "os_days_update", "VitalStatus")
Baseline_1yr_risk$Risk<- as.numeric(Baseline_1yr_risk$Risk)
Baseline_1yr_risk$os_days_update<- as.numeric(Baseline_1yr_risk$os_days_update)
##Assigning survival scores as "High or Low" based on median split
Baseline_1yr_risk$cat<- ifelse(Baseline_1yr_risk$Risk> median(Baseline_1yr_risk$Risk), "High Risk", "Low Risk")
Baseline_1yr_risk$survival_flag<- ifelse(Baseline_1yr_risk$VitalStatus== "Deceased", 1,0)
Baseline_1yr_Risk_KMC<- ggsurvplot(
  survfit(Surv(I(
    os_days_update
  )/30, survival_flag) ~ cat, Baseline_1yr_risk),
  pval = T,
  palette = "jama",
  pval.coord = c(0, 0.03),
  ylab = "OS Probability",
  xlab= "Time in Months",
  break.time.by=3,      
  risk.table = T,
  xlim= c(0,24),
  legend= "right",
  legend.title= "",
  legend.labs= c("High Score","Low Score"),
  fontsize=3,
  font.tickslab= c(10, "bold"),
  ggtheme= theme_classic2(base_size = 10),
  tables.y.text= F,
  tables.theme = theme_cleantable()
)
#Confidence Interval
survfit(Surv(I(os_days_update/30), survival_flag) ~ cat, data= Baseline_1yr_risk)
#P-Value
coxph(Surv(os_days_update/30, survival_flag)~ cat, data= Baseline_1yr_risk)



##Panel B, Figure 2
##baseline progression free survival logistic regression model scores
##Create and format data frame for analysis titled Baseline_1yr_risk, assigned survival scores from previous logistic regression analysis
Baseline_6mo_risk<- as.data.frame(cbind(PFS_baseline_survival_prediction, test_data_whole$pfs_days, test_data_whole$pfs_flag))
colnames(Baseline_6mo_risk)<- c("Risk", "pfs_days", "pfs_flag")
Baseline_6mo_risk$Risk<- as.numeric(Baseline_6mo_risk$Risk)
Baseline_6mo_risk$os_days_update<- as.numeric(Baseline_6mo_risk$pfs_days)
##Assigning survival scores as "High or Low" based on median split
Baseline_6mo_risk$cat<- ifelse(Baseline_6mo_risk$Risk> median(Baseline_6mo_risk$Risk), "High Risk", "Low Risk")
Baseline_6mo_Risk_KMC<- ggsurvplot(
  survfit(Surv(I(
    pfs_days
  )/30, pfs_flag) ~ cat, Baseline_6mo_risk),
  pval = T,
  palette = "jama",
  xlab="",
  pval.coord = c(0, 0.03),
  linetype = 2,
  ylab = "",
  break.time.by=2,      
  risk.table = F,
  xlim= c(0,24),
  legend= "none",
  ggtheme= theme_classic2(base_size = 10)
)
#Confidence Interval
survfit(Surv(I(pfs_days/30), pfs_flag) ~ cat, data= Baseline_6mo_risk)
#P-Value
coxph(Surv(pfs_days/30, pfs_flag)~ cat, data= Baseline_6mo_risk)

##Panel C, Figure 2
##l year OS optimized logistic regression model high and low score assignment and kaplan meier curve
##Create and format data frame for analysis titled Baseline_1yr_risk, assigned survival scores from previous logistic regression analysis
Optimized_1yr_risk<- as.data.frame(cbind(OS_optimize_survival_prediction, test_data_whole$os_days_update, test_data_whole$VitalStatus))
colnames(Optimized_1yr_risk)<- c("Risk", "os_days_update", "VitalStatus")
Optimized_1yr_risk$Risk<- as.numeric(Optimized_1yr_risk$Risk)
Optimized_1yr_risk$os_days_update<- as.numeric(Optimized_1yr_risk$os_days_update)
##Assigning survival scores as "High or Low" based on median split
Optimized_1yr_risk$cat<- ifelse(Optimized_1yr_risk$Risk> median(Optimized_1yr_risk$Risk), "High Risk", "Low Risk")
Optimized_1yr_risk$survival_flag<- ifelse(Optimized_1yr_risk$VitalStatus== "Deceased", 1,0)
Optimized_1year_Risk_KMC<- ggsurvplot(
  survfit(Surv(I(
    os_days_update
  )/30, survival_flag) ~ cat, Optimized_1yr_risk),
  pval = T,
  palette = "jama",
  pval.coord = c(0, 0.03),
  ylab = "OS Probability",
  xlab= "Time in Months",
  break.time.by=3,      
  risk.table = T,
  xlim= c(0,24),
  legend= "right",
  legend.title= "",
  legend.labs= c("High Score","Low Score"),
  fontsize=3,
  font.tickslab= c(10, "bold"),
  ggtheme= theme_classic2(base_size = 10),
  tables.y.text= F,
  tables.theme = theme_cleantable()
)
#Confidence Interval
survfit(Surv(I(os_days_update/30), survival_flag) ~ cat, data= Optimized_1yr_risk)
#P-Value
coxph(Surv(os_days_update/30, survival_flag)~ cat, data= Optimized_1yr_risk)

##Panel D
##optimized progression free survival logistic regression model scores
##Create and format data frame for analysis titled Baseline_1yr_risk, assigned survival scores from previous logistic regression analysis
Optimized_6mo_risk<- as.data.frame(cbind(PFS_optimize_survival_prediction, test_data_whole$pfs_days, test_data_whole$pfs_flag))
colnames(Optimized_6mo_risk)<- c("Risk", "pfs_days", "pfs_flag")
Optimized_6mo_risk$Risk<- as.numeric(Optimized_6mo_risk$Risk)
Optimized_6mo_risk$os_days_update<- as.numeric(Optimized_6mo_risk$pfs_days)
##Assigning survival scores as "High or Low" based on median split
Optimized_6mo_risk$cat<- ifelse(Optimized_6mo_risk$Risk> median(Optimized_6mo_risk$Risk), "High Risk", "Low Risk")
Optimized_6mo_Risk_KMC<- ggsurvplot(
  survfit(Surv(I(
    pfs_days
  )/30, pfs_flag) ~ cat, Optimized_6mo_risk),
  pval = T,
  palette = "jama",
  xlab="",
  pval.coord = c(0, 0.03),
  linetype = 2,
  ylab = "",
  break.time.by=2,      
  risk.table = F,
  xlim= c(0,24),
  legend= "none",
  ggtheme= theme_classic2(base_size = 10)
)
#Confidence Interval
survfit(Surv(I(pfs_days/30), pfs_flag) ~ cat, data= Optimized_6mo_risk)
#P-Value
coxph(Surv(pfs_days/30, pfs_flag)~ cat, data= Optimized_6mo_risk)

##Panel E, Figure 2 
##Machine learning predictions with high and low score groups with overall survival
OS_Figure_2E$HighLow<- ifelse(OS_Figure_2E$survival_score>median(OS_Figure_2E$survival_score), "High Survival Score", "Low Survival Score")
OS_Figure_2E_Plot<-  ggsurvplot(
  survfit(Surv(I(
    surv_months
  ), deceased) ~ HighLow, OS_Figure_2E),
  pval = T,
  palette = "jama",
  pval.coord = c(0, 0.03),
  ylab = "OS Probability",
  xlab= "Time in Months",
  break.time.by=3,      
  risk.table = T,
  xlim= c(0,24),
  legend= "right",
  legend.title= "",
  legend.labs= c("High Score","Low Score"),
  fontsize=3,
  font.tickslab= c(10, "bold"),
  ggtheme= theme_classic2(base_size = 10),
  tables.y.text= F,
  tables.theme = theme_cleantable()
)
#Confidence Interval
survfit(Surv(I(surv_months),  deceased) ~ HighLow, data= OS_Figure_2E)
#P-Value
coxph(Surv(surv_months, deceased)~ HighLow,data= OS_Figure_2E)

##Panel F, Figure 2
##Machine learning predictions with high and low score groups with overall survival
PFS_Figure_2F$HighLow<- ifelse(PFS_Figure_2F$survival_score>median(PFS_Figure_2F$survival_score), "High Survival Score", "Low Survival Score")
PFS_Figure_2F_Plot<-  ggsurvplot(
  survfit(Surv(I(
    surv_months
  ), progressed) ~ HighLow, PFS_Figure_2F),
  pval = T,
  palette = "jama",
  pval.coord = c(0, 0.03),
  ylab = "rwPFS Probability",
  xlab= "Time in Months",
  break.time.by=3,      
  risk.table = T,
  xlim= c(0,24),
  legend= "none",
  fontsize=3,
  font.tickslab= c(10, "bold"),
  ggtheme= theme_classic2(base_size = 10),
  tables.y.text= F,
  tables.theme = theme_cleantable()
)
#Confidence Interval
survfit(Surv(I(surv_months),  progressed) ~ HighLow, data= PFS_Figure_2F)
#P-Value
coxph(Surv(surv_months, progressed)~ HighLow,data= PFS_Figure_2F)

##Combining and formating survival curves for figure 2
A_Figure2<- Baseline_1yr_Risk_KMC
A_Figure2$plot<- A_Figure2$plot + labs(tag= "A") + theme(plot.tag.position = "topleft")

B_Figure2<- Baseline_6mo_Risk_KMC
B_Figure2$plot<- B_Figure2$plot + labs(tag= "B") + theme(plot.tag.position = "topleft")

C_Figure2<- Optimized_1year_Risk_KMC
C_Figure2$plot<- C_Figure2$plot + labs(tag = "C") + theme(plot.tag.position = "topleft")

D_Figure2<- Optimized_6mo_Risk_KMC
D_Figure2$plot<- D_Figure2$plot + labs(tag= "D") + theme(plot.tag.position = "topleft")

E_Figure2<- OS_Figure_2E_Plot
E_Figure2$plot<- E_Figure2$plot + labs(tag= "E") + theme(plot.tag.position = "topleft")

F_Figure2<- PFS_Figure_2F_Plot
F_Figure2$plot<- F_Figure2$plot + labs(tag= "F") + theme(plot.tag.position = "topleft")

Figure2_Combined<- list(A_Figure2, C_Figure2, E_Figure2, B_Figure2, D_Figure2, F_Figure2)
Figure2_Plot<- arrange_ggsurvplots(Figure2_Combined, ncol = 2, nrow = 3)
ggsave("Figure2.png", plot= Figure2_Plot, units= "in", height= 10, width = 8 , dpi= 300)




##Supplementary Figure 5: Forest plot of hazard ratios for multivariate 
##Cox-proportional hazards model for overall survival.
##From left to right, the variable is listed followed by its hazard ratio and 95% confidence interval 
##in parenthesis below it. The graphic representation has square points representing the hazard ratio value 
##and bracketed lines representing its 95% confidence interval. P-values are shown on the right with significant 
##variables. * P-value < 0.05; ** P-values < 0.01; ***P-values < 0.001. 
##Creation of data frame for analysis
OS_COXPH_Dataframe<- subset(all_data, select= c(male, Age.at.diagnosis, Adenocarcinoma, Squamous, White, Asian, Black, Smoking, StageIorII, StageIII, StageIV,
                                                 antipd1, antipdl1, priortki, irae_pos,irae_cat, early_irae, pdl1_neg, pdl1_149, pdl1_50100, os_days_update, survival_flag))
##Format variables
OS_COXPH_Dataframe$male<- as.numeric(OS_COXPH_Dataframe$male)
OS_COXPH_Dataframe$Adenocarcinoma<- as.numeric(OS_COXPH_Dataframe$Adenocarcinoma)
OS_COXPH_Dataframe$Squamous<- as.numeric(OS_COXPH_Dataframe$Squamous)
OS_COXPH_Dataframe$White<- as.numeric(OS_COXPH_Dataframe$White)
OS_COXPH_Dataframe$Asian<- as.numeric(OS_COXPH_Dataframe$Asian)
OS_COXPH_Dataframe$Black<- as.numeric(OS_COXPH_Dataframe$Black)
OS_COXPH_Dataframe$Smoking<- as.numeric(OS_COXPH_Dataframe$Smoking)
OS_COXPH_Dataframe$StageIorII<- as.numeric(OS_COXPH_Dataframe$StageIorII)
OS_COXPH_Dataframe$StageIII<- as.numeric(OS_COXPH_Dataframe$StageIII)
OS_COXPH_Dataframe$StageIV<- as.numeric(OS_COXPH_Dataframe$StageIV)
OS_COXPH_Dataframe$antipd1<- as.numeric(OS_COXPH_Dataframe$antipd1)
OS_COXPH_Dataframe$antipdl1<- as.numeric(OS_COXPH_Dataframe$antipdl1)
OS_COXPH_Dataframe$priortki<- as.numeric(OS_COXPH_Dataframe$priortki)
OS_COXPH_Dataframe$irae_pos<- as.numeric(OS_COXPH_Dataframe$irae_pos)
OS_COXPH_Dataframe$irae_cat<- as.numeric(OS_COXPH_Dataframe$irae_cat)
OS_COXPH_Dataframe$early_irae<- as.numeric(OS_COXPH_Dataframe$early_irae)
OS_COXPH_Dataframe$pdl1_neg<- as.numeric(OS_COXPH_Dataframe$pdl1_neg)
OS_COXPH_Dataframe$pdl1_149<- as.numeric(OS_COXPH_Dataframe$pdl1_149)
OS_COXPH_Dataframe$pdl1_50100<- as.numeric(OS_COXPH_Dataframe$pdl1_50100)
OS_coxmodel<- coxph(Surv(os_days_update, survival_flag)~ .,data= OS_COXPH_Dataframe)
##Create COXPH model
OS_Forest<- ggforest(OS_coxmodel, cpositions = c(0.02, -0.20, .35), fontsize = 1) 
##Print results
OS_Forest



#####Supplemental Figure 6, Forest plot of hazard ratios for multivariate 
##Cox-proportional hazards model for real-world progression-free survival (rwPFS)
##From left to right, the variable is listed followed by its hazard ratio and 95% confidence interval 
##in parenthesis below it. The graphic representation has square points representing the hazard ratio value
## and bracketed lines representing its 95% confidence interval. P-values are shown on the right with significant 
##variables.  * P-value < 0.05; ** P-values < 0.01; ***P-values < 0.001. 
##Creation of data frame for analysis
PFS_COXPH_Dataframe<- subset(all_data, select= c(male, Age.at.diagnosis, Adenocarcinoma, Squamous, White, Asian, Black, Smoking, StageIorII, StageIII, StageIV,
                                                  antipd1, antipdl1, priortki, irae_pos,irae_cat, early_irae, pdl1_neg, pdl1_149, pdl1_50100, pfs_days, pfs_flag))
##Format variables
PFS_COXPH_Dataframe$male<- as.numeric(PFS_COXPH_Dataframe$male)
PFS_COXPH_Dataframe$Adenocarcinoma<- as.numeric(PFS_COXPH_Dataframe$Adenocarcinoma)
PFS_COXPH_Dataframe$Squamous<- as.numeric(PFS_COXPH_Dataframe$Squamous)
PFS_COXPH_Dataframe$White<- as.numeric(PFS_COXPH_Dataframe$White)
PFS_COXPH_Dataframe$Asian<- as.numeric(PFS_COXPH_Dataframe$Asian)
PFS_COXPH_Dataframe$Black<- as.numeric(PFS_COXPH_Dataframe$Black)
PFS_COXPH_Dataframe$Smoking<- as.numeric(PFS_COXPH_Dataframe$Smoking)
PFS_COXPH_Dataframe$StageIorII<- as.numeric(PFS_COXPH_Dataframe$StageIorII)
PFS_COXPH_Dataframe$StageIII<- as.numeric(PFS_COXPH_Dataframe$StageIII)
PFS_COXPH_Dataframe$StageIV<- as.numeric(PFS_COXPH_Dataframe$StageIV)
PFS_COXPH_Dataframe$antipd1<- as.numeric(PFS_COXPH_Dataframe$antipd1)
PFS_COXPH_Dataframe$antipdl1<- as.numeric(PFS_COXPH_Dataframe$antipdl1)
PFS_COXPH_Dataframe$priortki<- as.numeric(PFS_COXPH_Dataframe$priortki)
PFS_COXPH_Dataframe$irae_pos<- as.numeric(PFS_COXPH_Dataframe$irae_pos)
PFS_COXPH_Dataframe$irae_cat<- as.numeric(PFS_COXPH_Dataframe$irae_cat)
PFS_COXPH_Dataframe$early_irae<- as.numeric(PFS_COXPH_Dataframe$early_irae)
PFS_COXPH_Dataframe$pdl1_neg<- as.numeric(PFS_COXPH_Dataframe$pdl1_neg)
PFS_COXPH_Dataframe$pdl1_149<- as.numeric(PFS_COXPH_Dataframe$pdl1_149)
PFS_COXPH_Dataframe$pdl1_50100<- as.numeric(PFS_COXPH_Dataframe$pdl1_50100)
PFS_coxmodel<- coxph(Surv(pfs_days, pfs_flag)~ .,data= PFS_COXPH_Dataframe)
##Create COXPH model
PFS_Forest<- ggforest(PFS_coxmodel, cpositions = c(0.02, -0.20, .35), fontsize = 1) 
##Print results
PFS_Forest




