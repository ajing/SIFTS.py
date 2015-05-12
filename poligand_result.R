#Some results

> fish_bs(protein_annotate_withsnp_site, "Disease", "center")

Fisher's Exact Test for Count Data

data:
p-value = 0.5528
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.7039081 1.2002189
sample estimates:
odds ratio
0.9190559

> fish_bs(protein_annotate_withsnp_site, "Polymorphism", "center")

Fisher's Exact Test for Count Data

data:
  p-value = 1.285e-13
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  2.392356 4.764285
sample estimates:
  odds ratio
3.363529

> fish_bs(protein_annotate_withsnp_site, "Unclassified", "center")

Fisher's Exact Test for Count Data

data:
p-value = 1.097e-11
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.2044393 0.4383788
sample estimates:
odds ratio
0.301712