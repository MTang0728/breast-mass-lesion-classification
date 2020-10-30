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
set.seed(40)

## data imputation ##
# visualize missing data 1
md.pattern(data, rotate.names = TRUE)
# visualize missing data 2
aggr(data,col=c("lightblue3","darkred"),numbers=TRUE,sortVars=TRUE,
     labels=names(data),cex.axis=.7,gap=3, ylab=c("Proportion missing","Missingness pattern"))
# apply multiple imputation and make 50 copies
data_imp <- mice(data, m = 50,
                 defaultMethod = c('norm', 'logreg', 'polyreg', 'polr'),
                 print = F)
# plot imputed and observed values for age
stripplot(x = data_imp, age ~ .imp, col=c("grey","darkred"), pch = c(1, 20), main = 'Age')
# plot density of different predictor
densityplot(data_imp)
## take out 2 copies for further analysis
data_imp_7 <- complete(data_imp, 7)
data_imp_28 <- complete(data_imp, 28)


# fit a logistic model on the data_imp_7
data_reg_7 <- glm(severity ~ birads + age + shape + margin + density, data = data_imp_7)
summary(data_reg_7)
# assess the independence and whole fit
rawresid7 <- residuals(data_reg_7,"resp")
# binned residual plots
#binned residual plots
binnedplot(x=fitted(data_reg_7),y=rawresid7,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
## There are 4 points outside the red band
# check if inverse logit function fits age
binnedplot(x=data_imp_7$age,y=rawresid7,xlab="Age",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy") 
## no patterns, pretty much all within the red line

# let's do the confusion matrix with .5 threshold
Conf_mat_7 <- confusionMatrix(as.factor(ifelse(fitted(data_reg_7) >= 0.5, "1","0")),
                              data_imp_28$severity_fac,positive = "1")
# look at couple evaluation metrics
Conf_mat_7$table
Conf_mat_7$overall["Accuracy"];
Conf_mat_7$byClass[c("Sensitivity","Specificity")]
# let's repeat with the marginal percentage in the data
Conf_mat_7 <- confusionMatrix(as.factor(ifelse(fitted(data_reg_7) >= mean(data_imp_7$severity), "1","0")),
                              data_imp_28$severity_fac,positive = "1")
# look at couple evaluation metrics
Conf_mat_7$table
Conf_mat_7$overall["Accuracy"];
Conf_mat_7$byClass[c("Sensitivity","Specificity")]
#look at ROC curve
roc(data_imp_7$severity,fitted(data_reg_7),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3") ## 

# fit a logistic model on the data_imp_28
data_reg_28 <- glm(severity ~ birads + age + shape + margin + density, data = data_imp_28)
summary(data_reg_28)
# assess the independence and whole fit
rawresid28 <- residuals(data_reg_28,"resp")
# binned residual plots
#binned residual plots
binnedplot(x=fitted(data_reg_28),y=rawresid28,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
## There are 3 points outside the red band
# check if inverse logit function fits age
binnedplot(x=data_imp_28$age,y=rawresid28,xlab="Age",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy") 
## no patterns, pretty much all within the red line

# let's do the confusion matrix with .5 threshold
Conf_mat_28 <- confusionMatrix(as.factor(ifelse(fitted(data_reg_28) >= 0.5, "1","0")),
                               data_imp_28$severity_fac,positive = "1")
# look at couple evaluation metrics
Conf_mat_28$table
Conf_mat_28$overall["Accuracy"];
Conf_mat_28$byClass[c("Sensitivity","Specificity")]
# let's repeat with the marginal percentage in the data
Conf_mat_28 <- confusionMatrix(as.factor(ifelse(fitted(data_reg_28) >= mean(data_imp_28$severity), "1","0")),
                               data_imp_28$severity_fac,positive = "1")
# look at couple evaluation metrics
Conf_mat_28$table
Conf_mat_28$overall["Accuracy"];
Conf_mat_28$byClass[c("Sensitivity","Specificity")]
#look at ROC curve
roc(data_imp_28$severity,fitted(data_reg_28),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3") ## 

### in general, inputation is ok, will proceed with inputation.

## perform EDA on one copy
data_imp_13 <- complete(data_imp, 13)
# check continuous predictor: age
ggplot(data_imp_13,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Age vs Severity",
       x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") ## different median and distribution, could be significant

# check categorical predictors
# bi-rads
table(data_imp_13$birads)
apply(table(data_imp_13[,c("severity_fac","birads")])/sum(table(data_imp_13[,c("severity_fac","birads")])),
      2,function(x) x/sum(x)) ## ratio of severity varies across birads, could be significant
# shape
table(data_imp_13$shape)
apply(table(data_imp_13[,c("severity_fac","shape")])/sum(table(data_imp_13[,c("severity_fac","shape")])),
      2,function(x) x/sum(x)) ## ratio of severity varies across shape, could be significant
# margin
table(data_imp_13$margin)
apply(table(data_imp_13[,c("severity_fac","margin")])/sum(table(data_imp_13[,c("severity_fac","margin")])),
      2,function(x) x/sum(x)) ## ratio of severity varies across margin, could be significant
# density
table(data_imp_13$density)
apply(table(data_imp_13[,c("severity_fac","density")])/sum(table(data_imp_13[,c("severity_fac","density")])),
      2,function(x) x/sum(x)) # ratio of severity consistent across density, might NOT be significant

## check interaction
# severity and age by margin
ggplot(data_imp_13,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Margin",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ margin,ncol=2) ## no interaction
# severity and shape by margin
table(data_imp_13$shape, data_imp_13$margin)
apply(table(data_imp_13[data_imp_13$margin == "circumscribed",c("severity_fac","shape")])/sum(table(data_imp_13[data_imp_13$margin == "circumscribed",c("severity_fac","shape")])),2,function(x) x/sum(x))
apply(table(data_imp_13[data_imp_13$margin == "microlobulated",c("severity_fac","shape")])/sum(table(data_imp_13[data_imp_13$margin == "microlobulated",c("severity_fac","shape")])),2,function(x) x/sum(x))
apply(table(data_imp_13[data_imp_13$margin == "obscured",c("severity_fac","shape")])/sum(table(data_imp_13[data_imp_13$margin == "obscured",c("severity_fac","shape")])),2,function(x) x/sum(x))
apply(table(data_imp_13[data_imp_13$margin == "ill-defined",c("severity_fac","shape")])/sum(table(data_imp_13[data_imp_13$margin == "ill-defined",c("severity_fac","shape")])),2,function(x) x/sum(x))
apply(table(data_imp_13[data_imp_13$margin == "spiculated",c("severity_fac","shape")])/sum(table(data_imp_13[data_imp_13$margin == "spiculated",c("severity_fac","shape")])),2,function(x) x/sum(x))
## patterns vary across margin, but not enough samples in certain combinations

# severity and birads by margin
table(data_imp_13$birads, data_imp_13$margin) ## 0 samples in certain combinations

# severity and density by margin
table(data_imp_13$density, data_imp_13$margin) ## 0 samples in certain combinations

# severity and age by birads
ggplot(data_imp_13,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Margin",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ birads,ncol=2) ## no interaction
# severity and shape by birads
table(data_imp_13$shape, data_imp_13$birads) ## 0 samples in certain combinations
# severity and density by birads
table(data_imp_13$density, data_imp_13$birads) ## 0 samples in certain combinations

# severity and age by shape
ggplot(data_imp_13,aes(x=severity_fac, y=age, fill=severity_fac)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Blues") +
  labs(title="Severity vs Age by Margin",x="Severity",y="Age") + 
  theme_classic() + theme(legend.position="none") +
  facet_wrap( ~ shape,ncol=2) ## no interaction
# severity and density by shape
table(data_imp_13$density, data_imp_13$shape) ## 0 samples in certain combinations

## in conclusion, thorugh visual inspection, all main effects except for density should be controlled for. No interactions need to be added.





















