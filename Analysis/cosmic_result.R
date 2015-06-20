> table(protein_annotate_cancer$VarType, protein_annotate_cancer$location)

Binding Site   Core Surface
Cancer       314369  39338  324779
> table(protein_annotate_cancer$location)

Binding Site         Core      Surface 
451184       540490      1701625 


> get_stat_eachtype(protein_annotate_cancer, "Cancer")
[1] "Observed"

Binding Site         Core      Surface 
314369        39338       324779 
[1] 678486
[1] "Expected:"

Binding Site         Core      Surface 
113660.6     136158.3     428667.1 
[1] "O/E"

Binding Site         Core      Surface 
2.7658569    0.2889138    0.7576485


> odds_ratio_stat(protein_annotate_cancer, "Cancer")
$data

Not Core   Core   Total
Not Cancer  1513661 501152 2014813
Cancer       639148  39338  678486
Total       2152809 540490 2693299

$measure
odds ratio with 95% C.I.
estimate     lower     upper
Not Cancer 1.000000        NA        NA
Cancer     0.185896 0.1839229 0.1878903

$p.value
two-sided
midp.exact fisher.exact chi.square
Not Cancer         NA           NA         NA
Cancer              0            0          0

#####################
Binding Site Surface   Total
Not Cancer       136815 1376846 1513661
Cancer           314369  324779  639148
Total            451184 1701625 2152809

$measure
odds ratio with 95% C.I.
estimate    lower     upper
Not Cancer 1.0000000       NA        NA
Cancer     0.1026589 0.101901 0.1034225

$p.value
two-sided
midp.exact fisher.exact chi.square
Not Cancer         NA           NA         NA
Cancer              0            0          0

############
Binding Site   Core  Total
Not Cancer       136815 501152 637967
Cancer           314369  39338 353707
Total            451184 540490 991674

$measure
odds ratio with 95% C.I.
estimate      lower      upper
Not Cancer 1.00000000         NA         NA
Cancer     0.03416149 0.03375174 0.03457622

$p.value
two-sided
midp.exact fisher.exact chi.square
Not Cancer         NA           NA         NA
Cancer              0            0          0

