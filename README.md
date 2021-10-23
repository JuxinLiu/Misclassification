# JointMisclassification
# Oct 22, 2021

For situations where the response variable and the covariate are simultaneously subject to misclassification errors, often an assumption of independent misclassification errors is used without justification. The aim of this project is to show the harmful consequences of inappropriate adjustment for such joint misclassification errors. In particular, we focus on the wrong adjustment by ignoring the dependence between the misclassification process of the response variable and that of the covariate.

When both variables are binary, we consider four different types of joint misclassification pattern, i.e., both variables subject to differential errors, both variables subject to nondifferential errors, one subject to differential errors while the other subject to nondifferential errors. For each type of joint misclassification pattern, the corresponding subfolder in folder entitled "binary" contains the R code for generating synthetic data and estimating model parameters for each of naive model (without any adjustment), independent model (wrong assumption), and dependent model (data generative model).

When both variables are trinary and ordinal, we only consider both subject to nondifferential misclassification errors. In the folder entitled "categorical", the R code for generating synthetic data and estimating model parameters for each of naive model (without any adjustment), independent model (wrong assumption), and dependent model (data generative model).
