
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