# Ultrasonic Texture Analysis for Acute Myocardial Infarction Risk Stratification: A Pilot Study


### Background:
Current risk stratification tools for acute myocardial infarction (AMI) have limitations, particularly in predicting mortality. This study utilizes cardiac ultrasound radiomics (i.e., ultrasomics) to risk stratify AMI patients when predicting all-cause mortality.

### Methods:
The study included 197 patients: a) retrospective internal cohort (n=155) of non-ST-elevation myocardial infarction (n=63) and ST-elevation myocardial infarction (n=92) patients, and b) external cohort from the multicenter Door-To-Unload in ST-segment–elevation myocardial infarction [DTU-STEMI] Pilot Trial (n=42). Echocardiography images of apical 2, 3, and 4-chamber were processed through an automated deep-learning pipeline to extract ultrasomic features. Unsupervised machine learning (topological data analysis) generated AMI clusters followed by a supervised classifier to generate individual predicted probabilities. Validation included assessing the incremental value of predicted probabilities over the Global Registry of Acute Coronary Events (GRACE) risk score 2.0 to predict 1-year all-cause mortality in the internal cohort and infarct size in the external cohort.

### Results:
Three phenogroups were identified: Cluster A (high-risk), Cluster B (intermediate-risk), and Cluster C (low-risk). Cluster A patients had decreased LV ejection fraction (P=0.004) and global longitudinal strain (P=0.027) and increased mortality at 1-year (log rank P=0.049). Ultrasomics features alone (C-Index: 0.74 vs. 0.70, P=0.039) and combined with global longitudinal strain (C-Index: 0.81 vs. 0.70, P<0.001) increased prediction of mortality beyond the GRACE 2.0 score. In the DTU-STEMI clinical trial, Cluster A was associated with larger infarcts size (>10% LV mass, P=0.003), compared to remaining clusters.

### Conclusions:
Ultrasomics-based phenogroup clustering, augmented by TDA and supervised machine learning, provides a novel approach for AMI risk stratification.


### Figure
![alt text](https://github.com/qahathaway/AMI_Phenogroups/blob/main/Figure1.jpg)

#### Study Design and Overview. (A) The internal validation patient cohort included patients with presenting with non-ST-elevation myocardial infarction (NSTEMI, n=63) and ST-elevation myocardial infarction (STEMI, n=92) who underwent echocardiography with views of the Apical 2-Chamber (A2C), Apical 3-Chamber (A3C), and Apical 4-Chamber (A4C). (B) Ultrasomics features were extracted using echocv and pyradiomics (v3.0.1). TDAView was used to cluster patients into three phenogroups: Cluster A, Cluster B, and Cluster C. The identified phenogroups were used to develop individual patient predicted probability of cluster assignment using a supervised machine learning classifier. (C) The generated probabilities from the supervised classifier were used to predict mortality and illustrate the incremental value of ultrasomics features over GRACE 2.0. Ultrasomics features were also extracted from the external validation group and applied to the supervised machine learning classifier to produce class labels (i.e., Cluster A, B, and C). The external validation phenogroups were used to predict findings on cardiac magnetic resonance, including acute infarct size.

### Figure
![alt text](https://github.com/qahathaway/AMI_Phenogroups/blob/main/FigureS1.jpg)

#### Illustration of Left Ventricular Ultrasomics in Apical 2-Chamber (A2C), Apical 3-Chamber (A3C), and Apical 4-Chamber (A4C). (A) Representative images for region of interest (ROI) placement using semantic segmentation through echocv. (B) Example of ultrasomics features extracted from the Python package pyradiomics (v3.0.1), including 1st order (n=18), shape-based (n=9), and texture-based (n=73) features. GLCM = gray-level cooccurrence matrix, GLDM = gray-level difference matrix, NGTDM = Neighborhood gray-tone difference matrix, GLRLM = gray-level run-length, GLSZM = gray-level size zone matrix.
