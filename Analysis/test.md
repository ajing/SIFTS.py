---
title: "sup_fig"
author: "ajing"
date: "Tuesday, March 24, 2015"
output: html_document
---


#```{r setup, include=FALSE}
#opts_chunk$set(dev = 'pdf')
#```



Create aa_change_preference function


# Amino acid change after mutation
## All locations
![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

##On surface (with binding sites)
![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

##On binding sites

```
##       aa_change           vartypes      estimate       ci_low       
##  ALA to ASP:  3   Disease     :158   Min.   :  0   Min.   :0.00000  
##  ALA to GLU:  3   Polymorphism:158   1st Qu.:  0   1st Qu.:0.00000  
##  ALA to GLY:  3   Unclassified:158   Median :  1   Median :0.03253  
##  ALA to PRO:  3                      Mean   :Inf   Mean   :0.16531  
##  ALA to SER:  3                      3rd Qu.:  2   3rd Qu.:0.17140  
##  ALA to THR:  3                      Max.   :Inf   Max.   :2.20559  
##  (Other)   :456                                                     
##      ci_up         pvalue        
##  Min.   :  1   Min.   :0.001318  
##  1st Qu.:  4   1st Qu.:0.252103  
##  Median : 10   Median :0.573674  
##  Mean   :Inf   Mean   :0.578123  
##  3rd Qu.:105   3rd Qu.:1.000000  
##  Max.   :Inf   Max.   :1.000000  
## 
```

```
## Error in matrix(value, n, p): 'data' must be of a vector type, was 'NULL'
```

In protein core

```
## Error in matrix(value, n, p): 'data' must be of a vector type, was 'NULL'
```

# Size change after mutation
![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

## For protein core
![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

## For surface (with binding site)
![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

## For binding site

```
##       volumn_change     vartypes  estimate    ci_low     ci_up
## 1   median to large      Disease 0.6078145 0.4080431 0.9040288
## 2   median to small      Disease 0.8180199 0.4912494 1.3690729
## 3    large to large      Disease 1.4712091 0.9734697 2.2511115
## 4   small to median      Disease 1.2267921 0.7904466 1.9256721
## 5    small to large      Disease 1.1203329 0.6855239 1.8546639
## 6  median to median      Disease 0.8356053 0.5895596 1.1868511
## 7    large to small      Disease 1.6283609 1.0522930 2.5608420
## 8   large to median      Disease 0.7309206 0.4927867 1.0852598
## 9    small to small      Disease 1.2308042 0.7479822 2.0575977
## 10  median to large Polymorphism 2.3168669 1.5008108 3.5434272
## 11  median to small Polymorphism 0.5933546 0.2661739 1.1937078
## 12   large to large Polymorphism 0.7097826 0.4050013 1.1914879
## 13  small to median Polymorphism 0.6201629 0.3253047 1.1106095
## 14   small to large Polymorphism 0.7920292 0.4017593 1.4630894
## 15 median to median Polymorphism 1.6326926 1.0910368 2.4177281
## 16   large to small Polymorphism 0.4653678 0.2341434 0.8561169
## 17  large to median Polymorphism 1.4321322 0.8997655 2.2377014
## 18   small to small Polymorphism 0.5279532 0.2379210 1.0556225
## 19  median to large Unclassified 0.7791219 0.4494787 1.2977119
## 20  median to small Unclassified 1.9341298 1.0939965 3.3354868
## 21   large to large Unclassified 0.7791219 0.4494787 1.2977119
## 22  small to median Unclassified 1.1352886 0.6616243 1.8860450
## 23   small to large Unclassified 1.0558473 0.5650870 1.8801796
## 24 median to median Unclassified 0.7500435 0.4648595 1.1770992
## 25   large to small Unclassified 0.9464313 0.5491258 1.5734760
## 26  large to median Unclassified 1.0845421 0.6590872 1.7374001
## 27   small to small Unclassified 1.2576567 0.6906734 2.2012299
##          pvalue
## 1  0.0108365318
## 2  0.4561335630
## 3  0.0627963299
## 4  0.3975877936
## 5  0.7224976396
## 6  0.3015997963
## 7  0.0247201901
## 8  0.1195245395
## 9  0.4071567750
## 10 0.0001186921
## 11 0.1736600467
## 12 0.1937987299
## 13 0.1229903773
## 14 0.5653297344
## 15 0.0157998407
## 16 0.0092928806
## 17 0.1246633182
## 18 0.0816156675
## 19 0.4023407200
## 20 0.0151240517
## 21 0.4023407200
## 22 0.6064850456
## 23 0.8849465483
## 24 0.2081355715
## 25 0.9012746537
## 26 0.7224039390
## 27 0.3859480584
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 




# hydrophibicity change after mutation


## For all locations
![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 

## For surface
![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png) 

## For binding site only
![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png) 


# Linear model

```
## 
## Call:
## glm(formula = is_disease ~ abs(hydro_change) + abs(size_change), 
##     family = "binomial")
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.5088  -1.1512   0.9462   1.1474   1.3298  
## 
## Coefficients:
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       -0.3862380  0.0510458  -7.566 3.83e-14 ***
## abs(hydro_change)  0.0035728  0.0005937   6.017 1.77e-09 ***
## abs(size_change)   0.0055358  0.0006429   8.611  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 9958.9  on 7188  degrees of freedom
## Residual deviance: 9856.7  on 7186  degrees of freedom
## AIC: 9862.7
## 
## Number of Fisher Scoring iterations: 4
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -0.35140 -0.12440  0.04665  0.06454  0.27230  0.75190
```

```
## [1] 0.5643344
```

```
## [1] 0.4356656
```

## With only disease and polymorphism

```
## 
## Call:
## glm(formula = is_disease ~ abs(hydro_change) + abs(size_change), 
##     family = "binomial")
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.6899  -1.4558   0.8478   0.9061   0.9596  
## 
## Coefficients:
##                    Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       0.5162465  0.0634145   8.141 3.93e-16 ***
## abs(hydro_change) 0.0006213  0.0007099   0.875    0.381    
## abs(size_change)  0.0038335  0.0007852   4.882 1.05e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 6882.1  on 5472  degrees of freedom
## Residual deviance: 6857.6  on 5470  degrees of freedom
## AIC: 6863.6
## 
## Number of Fisher Scoring iterations: 4
```

```
## [1] 0.6775078
```

```
## [1] 7189   41
```

```
## 
## Attaching package: 'dplyr'
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## Source: local data frame [3 x 5]
## 
##        VarType size.mean size.abs.mean hydro.mean hydro.abs.mean
## 1      Disease -3.758900      53.47627   2.179612       51.07120
## 2 Polymorphism  1.601399      44.42308  -7.571678       41.42016
## 3 Unclassified -3.093484      48.10198  -4.593201       50.81983
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png) ![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-2.png) 

```
## [1] "with abs"
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-3.png) ![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-4.png) 

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.5367  0.6446  0.7104  0.7463  0.8310  1.1540
```

```
## 
## FALSE  TRUE 
##  3481  3708
```

```
## [1] 3708
```

```
## [1] 5473
```

```
## [1] 3708
```

```
## [1] 5473
```

