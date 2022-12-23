This repository includes the data and codes for the manuscript "Predicting survival of NSCLC patients treated with immune checkpoint inhibitors: Impact and timing of immune-related adverse events and prior tyrosine kinase inhibitor therapy" by Sayer et al.

##Description of the data files

	1. `combined.data.csv` dataset contains all data from 354 patients used
	2. `training.data.csv` dataset contains data from 283 patients used as the training datasets for the predictive modeling
	3. `testing.data.csv` dataset contains data from 71 patients used as the testing datasets for the predictive modeling
	4. `NSCLC Data Dictionary` contains the description of the variables used in the datasets


##Description of the code files
	1. ```NSCLC_code.R``` contains the R codes for generating the figures in the manuscript using the datasets. Install the following R packages, which can be obtained using either the ```install.packages()``` function in R or via the [Bioconductor framework](http://www.bioconductor.org).
	
		* glmulti
		* survcomp
		* finalfit
		* caret
		* pROC
		* ROCR
		* MLmetrics


## Contact information

* Moom Roosan, PharmD, PhD [roosan@chapman.edu]
* Ravi Salgia, MD, PhD [rsalgia@coh.org]
