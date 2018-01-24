
Data and code for analyses described in the paper "Among-strain consistency in the pace and shape of senescence in duckweed", by Patrick Barks et al. (Journal of Ecology, in press).

### Description of data files

##### 1. data_demographic.csv

Column | Description
-----------------------------|--------------------------------------------------
data_demographic.csv         | 
`id`                         | frond identifier
`site`                       | site identifier
`strain`                     | strain identifier
`block`                      | temporal block identifier (`a` or `b`)
`tray`                       | tray identifier
`date_birth`                 | date of frond 'birth'
`discard`                    | was frond discarded from experiment? (T/F)
`uncertain_repro`            | was there uncertainty in tracking reproduction at some point in frond's life? (T/F)
`area`                       | frond surface area (cm^2)
`Jun.01.2014`:`Aug.31.2014`  | number of offspring detached from frond on given date
                             | 
data_areas_supplementary.csv | 
`site`                       | site abbreviation
`strain`                     | strain abbreviation
`block`                      | spatial block (`a` or `b`)
`area`                       | frond surface area (cm^2)
                             | 
data_sites.csv               | 
`site`                       | site abbreviation
`dropped`                    | was site dropped from main common garden experiment? (T/F)
`replicated`                 | was site replicated in main common garden experiment? (T/F)
`lat`                        | latitude (decimal degrees)
`lon`                        | longitude (decimal degrees)
`conductivity`               | surface water conductivity (μS / cm at 25 C)
`tdn`                        | surface water total dissolved nitrogen (μg / L)
`tdp`                        | surface water total dissolved phosphorus (μg / L)
`degree_days`                | degree-days above 10 C
`climate_station`            | nearest Environmental Canada climate station
`station_dist`               | distance to nearest climate station (km)
`nutrient_pc1`               | first principle component of ln-tdn and ln-tdp
