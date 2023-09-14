# ARM

We developed an advanced ARM framework and systematically investigated the interactions between microbial species using a public large-scale uniformly processed human microbiome data from the curatedMetagenomic Database together with ARM. 



In the rule mining part of this project, we obtain the relative abundance table from curatedMetagenomicData as our input. Then,

1. In Transaction_table.R, we retrieve the data from the curatedMetagenomicData package and filter the samples and taxa. And we partition the samples into groups with respect to their disease status (Healthy, IBD, CRC, IGT, T2D). Then, we transform the abundance table into transition table based on whether the taxa is present or absent in the sample. Here, we have the output in /ARM/Interim/transaction*.csv

2. In FP_growth_rule_mining.ipynb, we performed rule mining framework with FP-growth. The input is the transaction tables we generated above, and the output file was saved as /ARM/Interim/python_rules*.csv

3. In Clean_rules.R, we find the universal rules and do a visualization of top rules. The inputs are the rules mined above, and the output is in /ARM/Interim/univ_rule*.csv. The visualization step requires rule_viz.R which containing functions to plot our rules. The visualization result is shown in Figure 2.

4. In Freqitemset_Nestedness_analysis.R, we analyzed nestedness of the frequent item sets using based on the classical NODF measure and plot them for comparison. The inputs are the universal rule sets from step 3, and the outputs are NODF values and the plots shown in Figure 3.

5. In Rule_Comparison.R, we compare the universal rules between healthy and disease samples. The inputs are also the universal rule sets from step 3. We also generated plots shown in Figure 4.


In the classification part, there are four each classification model : 

1. Logit - logistic regression with LASSO
2. RF - Random Forest
3. SVM - Support Vector Machine
4. XGB - eXtreme Gradient Boosting decision trees 

And for each method, we tried two kinds of input to perform the classification:

1. RA - input as relative abundance table
2. CPAR - input as selected features from relative abundance table using predictive rules

Each R script is corresponding to the methods and input types.

Also, in functions_CPAR_classification.R includes the necessary function for *CPAR.R files to perform feature selection.

Finally, the performance was plotted using Performance_Plot.R. The input has the format as in /ARM?Interim/performace_list_final.RData, and the output is shown as in Figure 5b.


![alt text](https://github.com/Yiyannnnn/ARM/blob/0a03f81de4dd137acf607f5cc66e726810f88e05/Fig1_ARMworkflow.pdf)
