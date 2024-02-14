Feb 14 by Hongru

- remove blood cells from DasGupta study
- remove the duplicated sequencing runs (5 runs) from Andrews 2024 study (confirmed the sequencing run names with Rachel)
- move the C41/P3TLH from Andrews2022/Toronto_2 back to MacParland2018/Toronto_1 study
- covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

```
#                         sample sequencing_run donor library   cell
# Mullen                      23             23    23      23 137360
# DasGupta                    23             23    23      23  49921
# Grün                        4              4     4       4  29242
# Henderson                   29             29    17      29 105525
# Scott                       37             37    19      37  59201
# MacParland2018/Toronto1      4              4     4       4   6603
# Andrews2022/Toronto2        10             10     4      10  80130
# Andrews2024/Toronto3        26             26    21      26 176698

# sample sequencing_run          donor        library           cell
#    156            156            115            156         644680
```
