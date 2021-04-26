# BMI distribution estimation
Project created for *Analytical Tools Programing* (org. *Programowanie Narzędzi Analitycznych*) classes at WNE UW.

Language: 
 * Polish - classes, report and code comments.

Semester: I (MA studies).

## About
The main objective of this project was to apply distributions estimations methods, statistical hypothesis verification and practice programming in R, learned during classes. The project is about Body Mass Index (BMI) distribution estimation on dataset found on kaggle.com for the whole sample and for men and women separately. The idea was to assume some theoretical distributions and find proper parameters via Maximum Likelihood Estimation (MLE) and General Method of Moments (GMM). Assumed theoretical distributions were Normal (Gaussian) and 2-parametric Weibull. Godness of fit verification and estimation methods comparison were carried out using quantile-quantile plots (empirical distribution vs theoritical with estimated parameters). Research hypothesis was:

1. 37.5 distribution mean with 2.5 standard deviation which means obesity of the second degree (Normal distribution).
2. Scale parameter k=3 (Weibull distribution).
  
Joint hypothesis was tested via Likelihood Ratio test (MLE) and Wald test (GMM). Simple hypothesis was tested via z test.

Findings:  
* both Normal and Weibull distributions can be considered while BMI modelling, what is impressive because that proves Weibull distribution flexibility (it is considered for many tasks - eqrthquakes magnitudes, failure rates or winds modelling),
* no basics for accepting 1st hypothesis what is good because proves good society physical condition, 2nd hypothesis was rejected too,
* men has worse physical condition than women in terms of BMI and WHO guidelines.


## Repository description
* BMI.csv - data in .csv format,
* Projekt PNA Maciej Odziemczyk [388581].pdf - project report (in Polish),
* Projekt PNA Maciej Odziemczyk [388581].R - project code.

## Technologies
 * R (simple visualizations, estimations performed with written functions (MLE, GMM), tests written from scratch),
 * Word (report writing).

## Author
Maciej Odziemczyk

## Notes
* That was one of my first projects at WNE UW,
* to run the code pwd must be specified.
