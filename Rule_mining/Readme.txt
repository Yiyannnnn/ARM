{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset134 PingFangSC-Regular;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 In the rule mining part of this project, we obtain the relative abundance table from curatedMetagenomicData as our input. Then,\
\
1. In Transaction_table
\f1 .R, we retrieve the data from the 
\f0 curatedMetagenomicData package and filter the samples and taxa. And we partition the samples into groups with respect to their disease status (Healthy, IBD, CRC, IGT, T2D). Then, we 
\f1 transform
\f0  
\f1 the abundance table into transition table based on whether the taxa is present or absent in the sample. Here, we have the output in /ARM/Interim/transaction*.csv\
\
2. In FP_growth_rule_mining.ipynb, we performed rule mining framework with FP-growth. The input is the transaction tables we generated above, and the output file was saved as /ARM/Interim/python_rules*.csv\
\
3. In Clean_rules.R, we find the universal rules and do a visualization of top rules. The inputs are the rules mined above, and the output is in /ARM/Interim/univ_rule*.csv. The visualization step requires rule_viz.R which containing functions to plot our rules. The visualization result is shown in Figure 2.\
\
4. In Freqitemset_Nestedness_analysis.R, we analyzed nestedness of the frequent item sets using based on the classical NODF measure and plot them for comparison. The inputs are the universal rule sets from step 3, and the outputs are NODF values and the plots shown in Figure 3.\
\
5. In Rule_Comparison.R, we compare the universal rules between healthy and disease samples. The inputs are also the universal rule sets from step 3. We also generated plots shown in Figure 4.\
\
}