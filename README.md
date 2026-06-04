# Adjusting for Outcome Reporting Bias in Meta-analysis: A Multiple Imputation Approach

## Introduction
Selective reporting of outcomes in clinical trial reports can be described as the reporting of a subset of the originally recorded outcome variables in the final report. Selective outcome reporting can create outcome reporting bias (ORB) when the decision to report results within published studies is influenced by the significance or direction of the results [^1]. Previous reviews showed that statistically significant outcomes were more likely to be fully reported than non-significant outcomes [^2]. ORB potentially undermines the credibility and validity of meta-analyses [^3] and contributes to research waste by distorting overall treatment effects [^4]. ORB can be viewed as a missing data problem where unreported outcomes introduce bias. Despite the serious implications ORB poses, it remains an underrecognized issue, with only a few adjustment methods available.

## Methods
We propose an approach that addresses unreported outcomes in meta-analyses through multiple imputation. The imputed data are reweighted using importance sampling to pro-
vide an adjusted estimate of the treatment effect, building on existing methods for selection bias from the literature [^5]. To assess the impact of ORB in meta-analyses of clinical trials, we apply our proposed methodology to real clinical data affected by ORB. We compare the imputation of unreported outcomes in the univariate meta-analysis to the multivariate meta-analysis. Additionally, we conduct a simulation study to evaluate the method’s performance, focusing on treatment effect estimation across varying degrees of selective non-reporting.

## Repository Organization
Below is a description of the different folders and their respective contents. This ORBproject repository contains the code and source files of the paper "Adjusting for Outcome Reporting Bias in Meta-analysis: a Multiple Imputation Approach". 

### Paper Folder
The folder contains the file ImpSimORB.Rnw file, which can be compiled in R, so as to obtain the ImpSimORB.pdf file, also included in this folder. The pdf file is the main file of the manuscript for this project, dynamically producing text, tables, and figures obtained from the analyses. The folder additionally includes the bibliography file references.bib, as well as the figures in the respective folder.

[^1]: Cynthia MC Lemmens, Suzan van Amerongen, Eva M Strijbis, and Joep Killestein. Outcome reporting bias in clinical trials researching disease-modifying therapy in patients with multiple sclerosis. Neurology, 102(6):e208032, 2024. doi: https://doi.org/10.1212/WNL.0000000000208032.
[^2]: Kerry Dwan, Carrol Gamble, Paula R Williamson, Jamie J Kirkham, and Reporting Bias Group. Systematic review of the empirical evidence of study publication bias and outcome
reporting bias—an updated review. PloS one, 8(7):e66844, 2013. doi: https://doi.org/10.1371/journal.pone.0066844.
[^3]: John PA Ioannidis. Clinical trials: what a waste. BMJ, 349, 2014. doi: https://doi.org/10.1136/bmj.g7089.
[^4]: Elizabeth T Thomas and Carl Heneghan. Catalogue of bias: selective outcome reporting bias. BMJ Evidence-Based Medicine, 27(6):370–372, 2022. doi: https://doi.org/10.1136/
bmjebm-2021-111845.
[^5]: James Carpenter, Gerta Rücker, and Guido Schwarzer. Assessing the sensitivity of meta-analysis to selection bias: a multiple imputation approach. Biometrics, 67(3):1066–1072, 2011. doi: https://doi.org/10.1111/j.1541-0420.2010.01498.x
