# Time-dependent impact of immune-related adverse events on survival in NSCLC patients treated with immunotherapy

In this study, we utilized survival analysis and machine learning modeling to investigate how the presence of immune-related adverse events (irAE) and their timings may have an impact on clinical outcome on a cohort of non-small cell lung cancer patients treated with immune checkpoint inhibitors. 

## Dependencies

The following python package versions were used to conduct the machine learning analysis:
Python (v 3.9.7), scikit-learn (v 1.0.2), lifelines (v 0.26.4), and shap (v 0.39.0)

## How to Run the code

Place training and testing datasets within project folder. Once completed, run steps 01, 02, 03, and 04 sequentially to reproduce the analysis.

## File Descriptions

1. `step01_data_cleaning.ipynb`: Data wrangling
2. `step02_exploratory_analysis.ipynb`: Generation of heatmap to characterize the association between patient demographic features and survival outcome in the training data set
3. `step03_predictive_modeling_os.ipynb`: Machine learning modeling of OS after 1 year with SHAP model explanations
4. `step04_predictive_modeling_pfs.ipynb`: Machine learning modeling of PFS after 6 months with SHAP model explanations

