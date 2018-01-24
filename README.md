
Data and code for analyses described in the paper "Among-strain consistency in the pace and shape of senescence in duckweed" by Barks et al. (_Journal of Ecology_, in press).

&nbsp;

### Description of data files

Column | Description
-----------------------------|--------------------------------------------------
_data_demographic.csv_ |
`id`                         | frond identifier
`site`                       | site identifier
`strain`                     | strain identifier
`block`                      | temporal block identifier (a or b)
`tray`                       | tray identifier
`date_birth`                 | date of frond birth (yyyy-mm-dd)
`discard`                    | was frond discarded from experiment? (T/F)
`uncertain_repro`            | was reproduction uncertain at some point in frond's life? (T/F)
`area`                       | frond surface area (cm^2)
`Jun_01_2014`:`Aug_31_2014`  | number of offspring detached from frond on given date
&nbsp;                       |
_data_areas_supplementary.csv_ |
`site`                       | site identifier
`strain`                     | strain identifier
`block`                      | spatial block (a or b)
`area`                       | frond surface area (cm^2)
&nbsp;                       |
_data_sites.csv_ |
`site`                       | site identifier
`dropped`                    | was site dropped from main common garden experiment? (T/F)
`replicated`                 | was site replicated in main common garden experiment? (T/F)
`lat`                        | latitude (decimal degrees)
`lon`                        | longitude (decimal degrees)
`conductivity`               | surface water conductivity (μS/cm at 25 C)
`tdn`                        | surface water total dissolved nitrogen (μg/L)
`tdp`                        | surface water total dissolved phosphorus (μg/L)
`degree_days`                | degree-days above 10 C
`climate_station`            | nearest Environmental Canada climate station
`station_dist`               | distance to nearest climate station (km)
`nutrient_pc1`               | first principle component of ln-tdn and ln-tdp

&nbsp;

### Session Info

```
R version 3.4.2 (2017-09-28)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Sierra 10.12.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] maps_3.2.0           bindrcpp_0.2         rstan_2.16.2         StanHeaders_2.16.0-1
 [5] loo_1.1.0            mgcv_1.8-22          nlme_3.1-131         raster_2.6-7        
 [9] sp_1.2-5             geosphere_1.5-7      ape_5.0              RColorBrewer_1.1-2  
[13] gridExtra_2.3        ggplot2_2.2.1        data.table_1.10.4-3  tidyr_0.7.2         
[17] broom_0.4.3          tibble_1.3.4         dplyr_0.7.4          readr_1.1.1         

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.14       compiler_3.4.2     plyr_1.8.4         bindr_0.1         
 [5] tools_3.4.2        digest_0.6.12      gtable_0.2.0       lattice_0.20-35   
 [9] pkgconfig_2.0.1    rlang_0.1.4        Matrix_1.2-11      psych_1.7.8       
[13] mapproj_1.2-5      yaml_2.1.15        parallel_3.4.2     stringr_1.2.0     
[17] hms_0.4.0          tidyselect_0.2.3   stats4_3.4.2       inline_0.3.14     
[21] glue_1.2.0         R6_2.2.2           foreign_0.8-69     reshape2_1.4.2    
[25] purrr_0.2.4        magrittr_1.5       codetools_0.2-15   matrixStats_0.52.2
[29] scales_0.5.0       assertthat_0.2.0   mnormt_1.5-5       colorspace_1.3-2  
[33] labeling_0.3       stringi_1.1.6      lazyeval_0.2.1     munsell_0.4.3 
```

