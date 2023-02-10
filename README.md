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

compare_probs_observer(M, n_obs)

```

which plots the probabilities of selecting one condition over all others. Note that this presentation of the data does not involve scaling, which could obscure the patterns that are specific to an outlier. The black circles in the plot represent the answers of the potential outlier. We recommend performing a detailed per-observer analysis, rather than using an arbitrary measure to exclude observers. 

Once we are confident there are no outliers in the dataset, we can scale the results and compute confidence intervals using
```

[jod, stats]=pw_scale_bootstrp(M, boostrap_samples, options)

```

The function expects the same matrix of comparison per observer `M` as the outlier analysis and returns the scaling solution and a set of statistics. Confidence intervals represent the range in which the estimated quality values lie with 95\% confidence. The confidence intervals, however, should not be used to infer statistical significance of the difference. Use at least 500-1000 `bootstrap_samples` to obtain accurate confidence intervals. As bootstrapping takes time, you can pass 0 for `bootstrap_samples` to disable computation of confidence intervals. 

The cell-array `options` allows to select the type of prior, regularization, set the confidence level, and display options. Refer to the documentation of the function for details. The two most relevant options are:

* `prior` - The prior that ensures that the distance between pairs of conditions is finite. It helps to reduce the estimation error when the number of measurements is low. The default option is to use the `gaussian` prior, which is slower but more accurate than the `bounded` prior. Use `bounded` when scaling (or bootstrapping) takes too much time. Set to `none` to disable prior. 

* `regularization` - Since the quality scores in pairwise comparisons are relative and the absolute value cannot be determined, it is necessary to make an assumption how the absolute values are fixed in the optimization. The two options are  `mean0` (default) that makes mean JOD value equal to 0 and `fix0` which fixes the score of the first condition to 0. `mean0` results in a reduced overall estimation error as compared to `fix0`. However, `fix0` can be used to ensure that the estimation error is the smallest near the first condition. 

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

## Related projects

* [ASAP: Active Sampling for Pairwise Comparisons](https://github.com/gfxdisp/asap) - an effective method for reducing the number of required pairwise comparisons. 

## Revision history

* 2 Feb 2017 Simplified cost function makes computation faster without alterning results 

* 11 Dec 2016 Initial publicly available version

* 29 Apr 2020 Added `regularization` options: `mean0` and `fix0`. Cleaned up `ex4_lf.m`. 

* 16 Sep 2022 Added `pw_scale_table()` fuction. Fixed issues with Kolmogorov-Smirnov test.

* 10 Jan 2023 The code is now vectorized and provides gradients for optimization. The scaling is much faster now. [Thanks to Zhongshi Jiang]

## Literature

There is a number of papers describing the technique. When using the code, please cite [1]: 

[1] M. Perez-Ortiz and R. K. Mantiuk, “A practical guide and software for analysing pairwise comparison experiments”, arXiv Stat.AP, 2017, accessible at https://arxiv.org/abs/1712.03686

The scaling method and code has been used in the following papers (and many more): 

[2] K. Karaduzovic-Hadziabdic, J. H. Telalovic and R. K. Mantiuk, “Subjective and Objective Evaluation of Multi-exposure High Dynamic Range Image Deghosting Methods”, In Eurographics - Short Papers, 2016, pp. 29-32.

[3] G. Eilertsen, R. K. Mantiuk and J. Unger, “Real-time noise-aware tone mapping”. ACM Transactions on Graphics, 2015, 34(6), pp. 1-15.

[4] P. Vangorp, R. K. Mantiuk, B. Bazyluk, K. Myszkowski, R. Mantiuk, S. J. Watt and H.-P. Seidel, “Depth from HDR: depth induction or increased realism?” In ACM Symposium on Applied Perception - SAP, 2014, pp. 71-78. ACM Press.

[5] R. Wanat and R. K. Mantiuk, “Simulating and compensating changes in appearance between day and night vision”, ACM Transactions on Graphics (Proc. of SIGGRAPH), 2014, 33(4):147.

[6] M.H. Kedjar, G. Ward, H. Yoo, A. Soudi, T. Akhavan and C. Vazquez, “A Unified Color and Contrast Age-Dependent Visual Content Adaptation”, International Conference on Image Analysis and Processing ‑ ICIAP, 2017, pp. 765-778.

[7] G. Eilertsen, J. Kronander, G. Denes, R.K. Mantiuk and J. Unger, “HDR image reconstruction from a single exposure using deep CNNs”, ACM Transaction on Graphics (Proc. SIGGRAPH Asia), 2017, 36(6):178.

[8] E. Zerman, V. Hulusic, G. Valenzise, R.K. Mantiuk and F. Dufaux, “Effect of color space on high dynamic range video compression performance”, 9th International Conference on Quality of Multimedia Experience, QoMEX, 2017.

[9] V.K. Adhikarla, M. Vinkler, D. Sumin, et al., “Towards a quality metric for dense light fields”, In Proceedings of Computer Vision and Pattern Recognition (CVPR), 2017, pp 58-67.

[10] G. Eilertsen, J. Unger and R.K. Mantiuk, “Evaluation of Tone Mapping Operators for HDR Video”, High Dynamic Range Video (book), 2016, pp. 185-207.

[11] Chubarau, A., Akhavan, T., Yoo, H., Mantiuk, R. K., & Clark, J. (2020). Perceptual Image Quality Assessment for Various Viewing Conditions and Display Systems. Human Vision and Electronic Imaging.

[12] Zhong, F., Koulieris, G. A., Drettakis, G., Banks, M. S., Chambe, M., Durand, F., & Mantiuk, R. K. (2019). DiCE: Dichoptic Contrast Enhancement for VR and Stereo Displays. ACM Transactions on Graphics (TOG), 38(6), 1–13. https://doi.org/10.1145/3355089.3356552

[13] Denes, G., Maruszczyk, K., Ash, G., & Mantiuk, R. K. (2019). Temporal Resolution Multiplexing: Exploiting the limitations of spatio-temporal vision for more efficient VR rendering. IEEE Transactions on Visualization and Computer Graphics, 25(5), 2072–2082. https://doi.org/10.1109/TVCG.2019.2898741

[14] Perez-Ortiz, M., Mikhailiuk, A., Zerman, E., Hulusic, V., Valenzise, G., & Mantiuk, R. K. (2019). From pairwise comparisons and rating to a unified quality scale. IEEE Transactions on Image Processing, 1–1. https://doi.org/10.1109/tip.2019.2936103

[15] Mikhailiuk, A., Pérez-Ortiz, M., & Mantiuk, R. K. (2018). Psychometric scaling of TID2013 dataset. International Conference on Quiality of Multimedia Experience (QoMEX). https://doi.org/10.1109/QoMEX.2018.8463376

[16] Adams, W. J., Kucukoglu, G., Landy, M. S., & Mantiuk, R. K. (2018). Naturally glossy: Gloss perception, illumination statistics, and tone mapping. Journal of Vision, 18(13), 4. https://doi.org/10.1167/18.13.4

[17] Mantiuk, Rafał K., Gyorgy Denes, Alexandre Chapiro, Anton Kaplanyan, Gizem Rufo, Romain Bachy, Trisha Lian, and Anjul Patney. “FovVideoVDP : A Visible Difference Predictor for Wide Field-of-View Video.” ACM Transaction on Graphics 40, no. 4 (2021): 49. https://doi.org/10.1145/3450626.3459831.

[18] Hanji, Param, Rafał K. Mantiuk, Gabriel Eilertsen, Saghi Hajisharif, and Jonas Unger. “Comparison of Single Image HDR Reconstruction Methods — the Caveats of Quality Assessment.” In SIGGRAPH ’22 Conference Proceedings. Association for Computing Machinery, 2022. https://doi.org/10.1145/3528233.3530729.

A different version of the software was used in each of the above mentioned publications, which could result in small differences in the scaled results.

## License

MIT
