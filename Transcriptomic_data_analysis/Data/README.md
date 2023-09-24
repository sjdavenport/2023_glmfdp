## General info

data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42057
paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3824029/
supplement: https://www.atsjournals.org/doi/suppl/10.1165/rcmb.2012-0230OC/suppl_file/bahr_data-supplement.pdf

## Data description

> 42 control subjects and 94 subjects with varying severity of COPD had PBMC gene expression profiles generated. All subjects are non-hispanic white, current or former smokers.

## Data files

- Both files are CSV files.
- The first row contains row names.
- The first column contains column names.

### GSE42057_X.csv

136 x 12 design matrix: each row is a subject and each column is a variable of interest.

> head  GSE42057_X.csv
group,age_enroll,ats_packyears,bmi,distwalked,fev1_fvc_utah,fev1pp_utah,finalgold,gender,parentalcopd,smokcignow
GSM1031549,COPD,64.7,38,27.71,1600,0.46,37,3,1,1,0
GSM1031550,COPD,66.2,18,32.64,1805,0.65,54,2,1,2,0
GSM1031551,COPD,65.5,63,26.58,1620,0.49,53,2,0,0,0
GSM1031552,COPD,75.8,41,23.62,1010,0.43,32,3,0,0,0
GSM1031553,Control,61.2,39,30.82,1641,0.78,94,0,1,1,0
GSM1031554,Control,50.6,34.6,26.87,1460,0.76,96,0,0,2,1
GSM1031555,Control,60.6,50,24.4,2100,0.71,101,0,0,0,0
GSM1031556,COPD,55.8,46.8,17.92,538,0.69,60,2,0,2,0
GSM1031557,COPD,46.3,20,26.01,2109,0.66,84,1,1,0,0

### GSE42057_Y.csv

12531 x 136 gene expression data matrix: each row is a gene and each column is a subject.

## Comments on the model of Bahr et al vs ours

1) I've looked at these fev1* variables and they seem to be used for *defining* the COPD condition:

> COPD is defined as post–bronchodilator airflow obstruction with a ratio of forced expiratory volume in 1 second (FEV1) to forced expiratory volume (FVC) of less than 0.70, and is further subdivided into Stages 1–4, based on Global Initiative for Chronic Obstructive Lung Disease (GOLD) guidelines (24). Subjects with FEV1/FVC greater than or equal to 70 and FEV1 percentage of less than 80% post–bronchodilator were considered GOLD unclassified (25).

Therefore, it seems important to exclude the fev1 variables from the design, as they are probably quite correlated to the group defining variable.

2) Another point I have noted is that they do not run the same type of linear regression as we do. Indeed, section "Microarray Studies of PBMCs" states:

> separate linear models for the remaining 12,531 transcripts were fit for post-bronchodilator FEV1 percentage and FEV1/FVC across the 136 subjects. The total number of transcripts whose relative abundance was associated with FEV1 percentage, while controlling for demographic and smoking covariates, was 1,090. The total number of transcripts whose relative abundance was associated with FEV1/FVC while controlling for covariates was 1,745.

So what I understand is that they use each of the fev1 variables (so that's two sets of linear regressions) *as the response variable Y*, and for each gene v they run a linear regression of Y on X = covariates + the expression profile of v. That is, rather than looking at the model:

Y(v) = X beta(v) + epsilon(v)  (*)

where Y(v) is the expression profile of gene v and X the fixed matrix of covariates, they are instead looking at

Y^1 = X(v) beta^1(v) + epsilon(v)
Y^2 = X(v) beta^2(v) + epsilon(v)

where Y^1 and Y^2 are the two fev1 variables, and X(v) is the design matrix combining all covariates and gene v.

So they are looking at a different model and testing different assumptions. However, for me it is fine if we do (*) as proposed in the draft.

### Update (2021-11-22)

I've been able to reproduce the p-values reported in the paper using the linear model described in the abstract: "Multiple linear regression with adjustment for covariates (gender, age, body mass index, family history, smoking status, and pack-years)", that is:

Y = fev1_fvc_utah

X = expression + age_enroll + ats_packyears + bmi + gender + parentalcopd + smokcignow


This gives 1,745 genes selected by the BH procedure at the 0.05 target FDR level. 

NB: In the above model, an intercept is included.

Of course we get the same p-values if we swap "fev1_fvc_utah" and "expression", and we can test this model using your implementation.

Caveat: the variable dist_walked has been (silently?) excluded.  It turns out that this covariate is highly correlated with the fev* variables, so if one includes it in the model, the association between gene expression and fev* appears much weaker.

Conclusion: it seems that we can look at this model, with two contrasts that test the effect of (i) the fev1_fvc_utah variable and (ii) some other covariate (e.g. BMI? I guess that any other covariate will have 0 or very weak signal). 
