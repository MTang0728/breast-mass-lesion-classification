# breast-mass-lesion-classification

## Intro:
Breast cancer is the most commonly diagnosed cancer in women in the United States. Despite the low mortality rate of breast cancer, it is desired to detect and introduce intervention in the early stage to prevent further complications. The most effective method for screening breast cancer is mammography, which involves taking an x-ray image of the breast. Breast biopsy is also used in addition to mammography on suspicious lesion but is often an unpleasant process. Therefore, there has been a high demand in utilizing machine learning on mammography data to assist diagnosis. This analysis aims to determine the effects that several factors, including a breast lesion’s physical attributes have on a lesion’s severity. My final findings suggest that there are significant relationships between the odds of developing a malignant breast lesion and BI-RADS score, its mass shape and mass margin, as well as patient age, which has the strongest association with malignancy.

For details on the analysis and the findings, please refer the [report](./Resources/Report.pdf).

## Data
The data used used in this analysis is the [Mammography Mass Data Set](http://archive.ics.uci.edu/ml/datasets/Mammographic+Mass) sourced from the UCI Machine Learning Repository. This dataset contains 961 instances and 6 attributes including the response. There are 445 positive cases and 516 negative cases. The attribute detail is summarized below:

Name | Description | Type | Missing
-----|-------------|------|--------
BI-RADS | BI-RADS assessment ranging from 1 (definitely benign) to 6 (biopsy proven malignancy) | Ordinal | 5
Age | Patient age | integer | 5
Shape | Lesion mass shape, includes round, oval, lobular and irregular | nominal | 31
Margin | Lesion mass margin, includes circumscribed, microlobulated, obscured, ill-defined and spiculated | Nominal | 48
Density | Lesion mass density, includes high, iso, low and fat-containing | ordinal | 76
**Severity** | the response predictor that represents lesion severity, can be either benign or malignant | Binary | 0

## Statistical Methods Used:
- **Multiple Imputation** for computing missing values
- **Exploratory Data Analysis (EDA)** for exploring predictors and interactions to be controlled for
- **Logistic Regression** for fitting data
- **ANOVA & Chi-Square Test** for model and predictor selection
- **Accuracy, Sensitivity, Specificity, ROC** for evaluating model performance
- **P-value (Statistical Significance)** for model interpretation

## Tools Used:
- RStudio
