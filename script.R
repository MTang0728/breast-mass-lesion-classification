.libPaths(c('/usr/local/lib/R/library', .libPaths()))
library(ggplot2)
library(mice)
library(VIM)
library(rms)
library(arm)
library(pROC)
library(e1071)
library(caret)
library(lattice)
library(gridExtra)

# import data
data <- read.csv('./mammographic_masses.data', 
                 col.names = c('birads', 'age', 'shape', 'margin', 'density', 'severity'),
                 na.strings = '?')
# inspect data first
str(data)
summary(data)
## for birads, regular range is from 0 to 6, where 0 indicates an incomplete test. Therefore, I will
## replace 0 with NA to indicate missing values. There is one observation where its birads is 55, I will
## also replace it with NA.
data$birads[data$birads == 0 | data$birads == 55] <- NA
# convert data type
data$birads <- factor(data$birads)
data$shape <- factor(data$shape, labels = c('round', 'oval', 'lobular', 'irregular'))
data$margin <- factor(data$margin, labels = c('circumscribed',
                                              'microlobulated',
                                              'obscured',
                                              'ill-defined',
                                              'spiculated'))
data$density <- factor(data$density, labels = c('high', 'iso', 'low', 'fat-containing'))
# create a factorized version of the response variable
data$severity_fac <- factor(data$severity)



# set random seed for reproducibility
set.seed(1993)

## data imputation ##
# visualize missing data 1
md.pattern(data, rotate.names = TRUE)
# visualize missing data 2
aggr(data,col=c("lightblue3","darkred"),numbers=TRUE,sortVars=TRUE,
     labels=names(data),cex.axis=.7,gap=3, ylab=c("Proportion missing","Missingness pattern"))
# apply multiple imputation and make 10 copies
data_imp <- mice(data, m = 10,
                 defaultMethod = c('norm', 'logreg', 'polyreg', 'polr'),
                 print = F)
## used norm for continuous variable instead of pmm to impute with values potentially outside of the current range
# plot imputed and observed values for age
stripplot(x = data_imp, age ~ .imp, col=c("grey","darkred"), pch = c(1, 20), main = 'Age')
## imputed values look pretty evenly distributed across all 10 copies of the original dataset. All imputed values are within the range of the original non-missing values
# plot density of different predictor
densityplot(data_imp)
## the density plot of the observed values is unimodal. The density plot of imputed values from each copy generally only has 1 peak. Note that even though some density plots appear to be multi-modad, the valleys in these curves could come from the fact that there are only 5 missing values to be imputed. This sparsity makes it hard to construct a smooth unimodal curve when when the imputed values span a wide range. 



# take out 1 copy for further analysis
data_imp_7 <- complete(data_imp, 7)
# fit a logistic model on the data_imp_7
full_formula <-formula(severity ~ birads + age + shape + margin + density)
data_reg_7 <- glm(full_formula, data = data_imp_7, family = binomial)
summary(data_reg_7)
# assess the independence and whole fit
rawresid_7 <- residuals(data_reg_7,"resp")
# binned residual plots
#binned residual plots
binnedplot(x=fitted(data_reg_7),y=rawresid_7,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
## There are 2 points outside the red band
# check if inverse logit function fits age
binnedplot(x=data_imp_7$age,y=rawresid_7,xlab="Age",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy") 
## no patterns, there are 2 points outside the red band

# let's do the confusion matrix with .5 threshold
ConfMat_7 <- confusionMatrix(as.factor(ifelse(fitted(data_reg_7) >= 0.5, "1","0")),
                             data_imp_7$severity_fac,positive = "1")
# look at couple evaluation metrics
ConfMat_7$table
ConfMat_7$overall["Accuracy"]; ## accuracy is 0.840625
ConfMat_7$byClass[c("Sensitivity","Specificity")] ## 0.8153153   0.8624031

# let's repeat with the marginal percentage in the data
ConfMat_7_mean <- confusionMatrix(as.factor(ifelse(fitted(data_reg_7) >= 
                                                     mean(data_imp_7$severity), "1","0")),
                                  data_imp_7$severity_fac,positive = "1")
# look at couple evaluation metrics
ConfMat_7_mean$table
ConfMat_7_mean$overall["Accuracy"]; ## accuracy is 0.8427083
ConfMat_7_mean$byClass[c("Sensitivity","Specificity")] ## 0.8400901   0.8449612 
#look at ROC curve
roc(data_imp_7$severity,fitted(data_reg_7),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3") 
## AUC: 0.9127     threshold: 0.399, sensitivity: 0.828, specificity: 0.869



## perform EDA on this copy
# check continuous predictors
ggplot(data_imp_7,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Age vs Severity",
       x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") ## different median and distribution, could be significant

# check categorical predictors
# bi-rads
table(data_imp_7$birads) ## not super unbalanced
apply(table(data_imp_7[,c("severity_fac","birads")])/sum(table(data_imp_7[,c("severity_fac","birads")])),
      2,function(x) x/sum(x)) ## ratio of severity varies across birads, could be significant
# shape
table(data_imp_7$shape) ## balanced
apply(table(data_imp_7[,c("severity_fac","shape")])/sum(table(data_imp_7[,c("severity_fac","shape")])),
      2,function(x) x/sum(x)) ## ratio of severity varies across shape, could be significant
# margin
table(data_imp_7$margin) ## not super unbalanced
apply(table(data_imp_7[,c("severity_fac","margin")])/sum(table(data_imp_7[,c("severity_fac","margin")])),
      2,function(x) x/sum(x)) ## ratio of severity varies across margin, could be significant
# density
table(data_imp_7$density) ## not super unbalanced
apply(table(data_imp_7[,c("severity_fac","density")])/sum(table(data_imp_7[,c("severity_fac","density")])),
      2,function(x) x/sum(x)) # ratio of severity consistent across density, might NOT be significant


## check interaction
# severity and age by birads
g1 <- ggplot(data_imp_7,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Birads",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ birads,ncol=3)  ## trend consistent across all levels, no interaction
# severity and age by shape
g2 <- ggplot(data_imp_7,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Shape",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ shape,ncol=2)  ## trend consistent across all levels, no interaction
# severity and age by margin
g3 <- ggplot(data_imp_7,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Margin",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ margin,ncol=3) ## trend consistent across all levels, no interaction
# severity and age by density
g4 <- ggplot(data_imp_7,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Density",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ density,ncol=2) ## trend consistent across all levels, no interaction
grid.arrange(g1, g2, g3, g4, ncol = 4)

# severity and density by birads
table(data_imp_7$density, data_imp_7$birads) ## 0 samples in certain combinations
# severity and density by shape
table(data_imp_7$density, data_imp_7$shape) ## 0 samples in certain combinations
# severity and density by margin
table(data_imp_7$density, data_imp_7$margin) ## 0 samples in certain combinations

# severity and margin by birads
table(data_imp_7$margin, data_imp_7$birads) ## 0 samples in certain combinations
# severity and margin by shape
table(data_imp_7$margin, data_imp_7$shape) ## sample size too small in certain combinations

# severity and shape by birads
table(data_imp_7$shape, data_imp_7$birads) ## 0 samples in certain combinations

## in conclusion, through visual inspection, all main effects except for density should be controlled for. No interactions need to be added.



# build a null model that contains only the effect i care about: shape and density
reg_7_null <- glm(severity ~ shape + density, data = data_imp_7, family = binomial)
# use data_reg_7 as the full model
reg_7_full <- data_reg_7
# perform model selection
# do forward selection first
model_forward <- step(reg_7_null, scope = formula(reg_7_full), direction = 'forward', trace = 0)
model_forward$call ## shape + density + birads + age + margin
# repeat the process using stepwise selection
model_stepwise <- step(reg_7_null, scope = formula(reg_7_full),direction="both", trace=0)
model_stepwise$call ## has shape + birads + age + margin
# repeat the process using backward selection
model_backward <- step(reg_7_full, direction = 'backward', trace = 0)
model_backward$call ## has birads + age + shape + margin
# since forward selection methods produced different results, check which one to keep
model_select_7 <- glm(formula = severity ~ shape + birads + age + margin, family = binomial, 
                      data = data_imp_7)
anova(model_select_7, reg_7_full, test = 'Chisq')
## since the p-value of F-test is 0.2133, which is greater than 0.1. Thus we fail to reject the null hypothesis, and will keep the result of stepwise and backward selections.

summary(model_select_7)
# assess the independence and whole fit
rawresid_7_final <- residuals(model_select_7,"resp")
# binned residual plots
#binned residual plots
binnedplot(x=fitted(model_select_7),y=rawresid_7_final,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
## There are 3 points outside the red band
# check if inverse logit function fits age
binnedplot(x=data_imp_7$age,y=rawresid_7_final,xlab="Age",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy") 
## no patterns, there are 2 points outside the red band

# let's do the confusion matrix with .5 threshold
ConfMat_7_final <- confusionMatrix(as.factor(ifelse(fitted(model_select_7) >= 0.5, "1","0")),
                                   data_imp_7$severity_fac,positive = "1")
# look at couple evaluation metrics
ConfMat_7_final$table
ConfMat_7_final$overall["Accuracy"]; ## accuracy is 0.840625
ConfMat_7_final$byClass[c("Sensitivity","Specificity")] ## 0.8175676   0.8604651 

# let's repeat with the marginal percentage in the data
ConfMat_7_mean_final <- confusionMatrix(as.factor(ifelse(fitted(model_select_7) 
                                                         >= mean(data_imp_7$severity), "1","0")),
                                        data_imp_7$severity_fac,positive = "1")
# look at couple evaluation metrics
ConfMat_7_mean_final$table
ConfMat_7_mean_final$overall["Accuracy"]; ## accuracy is 0.8427083
ConfMat_7_mean_final$byClass[c("Sensitivity","Specificity")] ## 0.8333333   0.8507752
#look at ROC curve
roc(data_imp_7$severity,fitted(model_select_7),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3") 
## AUC: 0.911     threshold: 0.409, sensitivity: 0.831, specificity: 0.865



# take another copy
data_imp_2 <- complete(data_imp, 2)
# fit a logistic model on data_imp_2 using all main effects
data_reg_2 <- glm(full_formula, data = data_imp_2, family = binomial)
summary(data_reg_2) ## result is similar to what the 7th copy produced
# fit a logistic model on data_imp_2 using the output of model selection
model_select_2 <- glm(formula = severity ~ shape + birads + age + margin, family = binomial, 
                      data = data_imp_2)
anova(model_select_2, data_reg_2, test = 'Chisq')
## since the p-value of F-test is 0.632, which is greater than 0.1. Thus we fail to reject the null hypothesis, and will keep the result of stepwise and backward selections. This finding is consistent with what the 7th copy produced. 
# check the model summary
summary(model_select_2)
## the coefficients are similar to what we have seem previously. The predictors that we have previously determined to be statistically significant are still significant in this copy. Except that now marginmicrolobulated is also statistically significant at level 0.1
# assess the independence and whole fit
rawresid_2 <- residuals(model_select_2,"resp")
# binned residual plots
#binned residual plots
binnedplot(x=fitted(model_select_2),y=rawresid_2,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
## There are 3 points outside the red band
# check if inverse logit function fits age
binnedplot(x=data_imp_2$age,y=rawresid_2,xlab="Age",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy") 
## no patterns, there are 2 points outside the red band

# let's do the confusion matrix with .5 threshold
ConfMat_2 <- confusionMatrix(as.factor(ifelse(fitted(model_select_2) >= 0.5, "1","0")),
                            data_imp_2$severity_fac,positive = "1")
# look at couple evaluation metrics
ConfMat_2$table
ConfMat_2$overall["Accuracy"]; ## accuracy is 0.8322917
ConfMat_2$byClass[c("Sensitivity","Specificity")] ## 0.8085586   0.8527132
# let's repeat with the marginal percentage in the data
ConfMat_2_mean <- confusionMatrix(as.factor(ifelse(fitted(model_select_2) >=
                                                     mean(data_imp_2$severity), "1","0")),
                                  data_imp_2$severity_fac,positive = "1")
# look at couple evaluation metrics
ConfMat_2_mean$table
ConfMat_2_mean$overall["Accuracy"]; ## accuracy is 0.8395833
ConfMat_2_mean$byClass[c("Sensitivity","Specificity")] ## 0.8400901   0.8391473
#look at ROC curve
roc(data_imp_2$severity,fitted(model_select_2),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3") ## 
## AUC: 0.9104     threshold: 0.333, sensitivity: 0.800, specificity: 0.885

### in general, all findings are pretty similar and consistent across 2 different copies of the original dataset


## apply the same model on all 10 datasets
data_reg_imp <- with(data = data_imp, 
                     glm(severity ~ shape + birads + age + margin, family = binomial))
# results for the 2nd and 7th dataset
data_reg_imp[[4]][[2]]
data_reg_imp[[4]][[7]]
#get the MI inferences based on the Rubin combining rules
data_reg <- pool(data_reg_imp)
summary(data_reg)

