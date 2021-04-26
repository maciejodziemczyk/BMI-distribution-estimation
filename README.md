# BMI distribution estimation
Project created for Analytical Tools Programing (org. Programowanie Narzędzi Analitycznych) at WNE UW

Language: Polish - report and code comments 
Semester: I (MA studies)

## About
The main objective of this project was to train distributions estimations methods, statistical hypotesis veryfing and R programming language basics learned during classes. The project is about Body Mass Index (BMI) distribution estimation on dataset found on kaggle.com for the whole sample and for mans and women separately. The idea was to assume some theoretical distributions and find proper parameters via Maximum Likelihood Estimation (MLE) and General Method of Moments (GMM). Assumed theoretical distributions were Normal (Gaussian) and 2-parametric Weibull. Godness of fit and estimation methods comparison was verified with quantile-quantile plots (empirical distribution vs theoritical with estimated parameters). Research hypotesis was:

  1. 37.5 distribution mean with 2.5 standard deviation which means secondary obesity (Normal distribution)
  2. scale parameter k=3 (Weibull distribution)
  
Joint hypotesis was tested via Likelihood Ratio test (MLE) and Wald test (GMM). Simple hypotesis was tested via z test.

Findings:
  
  *Both Normal and Weibull distributions can be considered while BMI modelling, what is impressive because that proves Weibull distribution flexibility (it is considered for many tasks - eqrthquakes magnitudes, failure rates or winds modelling),
  *no basics for accepting 1. hypotesis what is good because proves good society physical condition, 2. hypotesis was rejected too,
  *Mans has worse physical condition than women in terms of BMI
  
