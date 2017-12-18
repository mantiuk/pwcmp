# pwcmp 

This is a set of matlab functions for scaling of pairwise comparison experiment results based on Thurstone's model V assumptions.

The main features:

* The scaling can work with imbalanced and incomplete data, in which not all pairs are compared and some pairs are compared more often than the others.

* Additional priors reduce bias due to the non-linear nature of the problem.

* Outlier rejection to screen observers that perform differently than the rest.

* The code can compute confidence intervals using bootstrapping. 

## Usage


The code for this example can be found in the "examples" folder, under the name of video_TMO_analysis_example, in which we analyse the data from a video tone mapping evaluation project.

We recommend to keep the data in a tabulated format, such as comma-separated-files (CSV), in which each condition is described by meaningful labels. Such files are easy to read with any software and can be easily interpreted even long after the data have been collected. The following table shows a few rows from the analysed dataset.


| Observer | Session | Scene | Condition_1 | Condition_2 | Selection |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 1 | 1 | Window	| TMO_Camera	| Ferwerda96 | 1 |
| 1 | 1 | Exhibition | Ronan12 | Irawan05	| 2 |
| 1 | 1 | Corridor | Irawan05 | Ferwerda96 |	1 |
| 2 | 2 | Corridor |	Ronan12 | TMO_Camera | 2|


The first step is to convert the answers from the table into a set of comparison matrices M, one matrix per each observer. In such a matrix, columns and rows correspond to compared conditions and matrix value c_{ij}=n means that condition O_i was n times selected as better than condition O_j. If there is a reference condition, such as a non-distorted image, it should be put in the matrix as the first condition in the first row and column. The first condition will be assigned a fixed quality value of 0. 

The second step is to perform outlier analysis to detect potential observers who performed very differently from the rest. The function to perform this analysis is

```

[L,L_dist]=pw_outlier_analysis(M)

```

which receives a matrix M with the responses per observer and returns the likelihood L of observing the data of each observer and a inter-quartile-normalised score L_dist, which indicates the observers that should be further investigated. Since there is no objective threshold that could distinguish outliers with high confidence, we advise to investigate all observers whose L_dist score is close or above the customary threshold of 1.5.

The results for the 18 observers in the analysed dataset indicate that there is one observer with a score of 2.72, which requires further attention.

To compare the answers of the indicated observer (observer number n_obs) to the rest of observers, we use the function 


```

compare_probs_observer(M,n_obs)

```

which plots the probabilities of selecting one condition over all others. Note that this presentation of the data does not involve scaling, which could obscure the patterns that are specific to an outlier. The black circles in the plot represent the answers of the potential outlier. We recommend performing a detailed per-observer analysis, rather than using an arbitrary measure to exclude observers. 

Once we are confident there are no outliers in the dataset, we can scale the results and compute confidence intervals using
```

[jod, stats]=pw_scale_bootstrp(M)

```

The function expects the same matrix of comparison per observer M as the outlier analysis and returns the scaling solution and a set of statistics. Confidence intervals represent the range in which the estimated quality values lie with 95\% confidence. The confidence intervals, however, should not be used to infer statistical significance of the difference. 

Statistical tests can be performed using the function 

```

pw_plot_ranking_triangles(jod,stats)

```
which produces a plot that can be used to interpret data in terms of statistical significance. The continuous lines in this plot indicate statistically significant difference between the pair of conditions and the dashed lines indicate the lack of evidence for statistically significant difference.

If we do not have access to the comparison matrix per observer or we do not want to analyse confidence intervals the following scaling function can be used instead


```
% Simple example showing how to execute scaling method

% Comparison matrix. This is an incomplete design in which 1-st conditon was
% compared with the 2nd and 2nd with the 3rd. In both cases 75 observers
% seleted one condition and 25 the other.
D = [ 0   25  0;
      75  0  25;
      0   75 0  ];
   
Q = pw_scale( D );

display( Q )
```


See more examples in "examples" directory.

## Revision history

v0.1 (14.12.2016) Internal release

## Literature

There is a number of papers describing the technique. When using the code, please cite [1]: 

[1] M. Perez-Ortiz and R. K. Mantiuk, “A practical guide and software for analysing pairwise comparison experiments”, arXiv Stat.AP, 2017, accessible at https://arxiv.org/abs/1712.03686

## License

MIT
