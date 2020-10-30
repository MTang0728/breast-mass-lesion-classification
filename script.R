.libPaths(c('/usr/local/lib/R/library', .libPaths()))
library(ggplot2)
library(mice)
library(VIM)
library(lattice)
library(gridExtra)

# import data
data <- read.csv('./mammographic_masses.data', 
                 col.names = c('birads', 'age', 'shape', 'margin', 'density', 'severity'),
                 na.strings = '?')
str(data)
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
severity_fac <- factor(data$severity)
# set random seed
set.seed(40)

## data imputation ##
# visualize missing data 1
md.pattern(data, rotate.names = TRUE)
# visualize missing data 2
aggr(data,col=c("lightblue3","darkred"),numbers=TRUE,sortVars=TRUE,
     labels=names(data),cex.axis=.7,gap=3, ylab=c("Proportion missing","Missingness pattern"))
# visualize patterns of missing data 3
marginplot([,c("diameter","age")],col=c("lightblue3","darkred"),cex.numbers=1.2,pch=19)