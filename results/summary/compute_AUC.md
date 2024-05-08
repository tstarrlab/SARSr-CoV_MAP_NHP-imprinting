Compute per-barcode binding to vaccinated mouse sera samples
================
Tyler Starr
5/12/2022

This notebook reads in per-barcode counts from `count_variants.ipynb`
for sera-binding titration experiments, computes functional scores for
RBD binding values via delta-AUC metrics, and does some basic QC on
variant binding functional scores.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$sera_delta_AUC_dir)){
  dir.create(file.path(config$sera_delta_AUC_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.8 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/linux-rocky8-nehalem/gcc-8.5.0/intel-oneapi-mkl-2021.4.0-h43nkmwzvaltaa6ii5l7n6e7ruvjbmnv/mkl/2021.4.0/lib/intel64/libmkl_rt.so.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [5] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ##  [9] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## [13] knitr_1.37       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.30        haven_2.4.3      colorspace_2.0-3
    ##  [5] vctrs_0.5.2      generics_0.1.2   htmltools_0.5.2  utf8_1.2.2      
    ##  [9] rlang_1.0.6      pillar_1.7.0     glue_1.6.2       withr_2.5.0     
    ## [13] DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1    
    ## [17] lifecycle_1.0.3  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_1.0.2      evaluate_0.15    tzdb_0.2.0       fastmap_1.1.0   
    ## [25] fansi_1.0.2      broom_0.7.12     Rcpp_1.0.11      backports_1.4.1 
    ## [29] scales_1.2.1     jsonlite_1.8.7   fs_1.5.2         hms_1.1.1       
    ## [33] digest_0.6.29    stringi_1.7.6    grid_4.1.3       cli_3.6.0       
    ## [37] tools_4.1.3      magrittr_2.0.2   crayon_1.5.0     pkgconfig_2.0.3 
    ## [41] ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0 
    ## [45] rstudioapi_0.13  assertthat_0.2.1 rmarkdown_2.13   httr_1.4.7      
    ## [49] R6_2.5.1         compiler_4.1.3

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#read file giving count of each barcode in each sort partition
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F))

#read in barcode-variant lookup tables
dt <- data.table(read.csv(file=config$codon_variant_table,stringsAsFactors=F))

#rename some targets
dt[target=="XBB15",target:="SARS-CoV-2_XBB15"]
dt[target=="SARS-CoV-1_2693",target:="SARS-CoV-1_Urbani"]
dt[target=="SARS-CoV-1_Urbani_HP03L",target:="SARS-CoV-1_Urbani"]
dt[target=="Wuhan_Hu_1",target:="SARS-CoV-2_WH1"]

setkey(dt,barcode,library)

#eliminate duplicated barcodes
duplicates <- dt[duplicated(dt,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplciates, and then remove
dt[,duplicate:=FALSE]
for(i in 1:nrow(duplicates)){
  dt[library==duplicates[i,library] & barcode==duplicates[i,barcode],duplicate:=TRUE]
}
dt <- dt[duplicate==FALSE,]; dt[,duplicate:=NULL]

dt <- merge(counts, dt, by=c("library","barcode"));rm(counts); rm(duplicates)

samples_D5338_3 <- data.frame(sample=sort(unique(paste(rep("D5338-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D5338-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6391_3 <- data.frame(sample=sort(unique(paste(rep("D6391-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D6391-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6404_3 <- data.frame(sample=sort(unique(paste(rep("D6404-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D6404-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5733_3 <- data.frame(sample=sort(unique(paste(rep("D5733-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D5733-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6343_3 <- data.frame(sample=sort(unique(paste(rep("D6343-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D6343-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6271_3 <- data.frame(sample=sort(unique(paste(rep("D6271-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D6271-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5220_3 <- data.frame(sample=sort(unique(paste(rep("D5220-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D5220-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5417_3 <- data.frame(sample=sort(unique(paste(rep("D5417-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D5417-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5379_3 <- data.frame(sample=sort(unique(paste(rep("D5379-3",3),formatC(barcode_runs[barcode_runs$sample_type=="D5379-3","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5338_4 <- data.frame(sample=sort(unique(paste(rep("D5338-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D5338-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6391_4 <- data.frame(sample=sort(unique(paste(rep("D6391-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D6391-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6404_4 <- data.frame(sample=sort(unique(paste(rep("D6404-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D6404-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5733_4 <- data.frame(sample=sort(unique(paste(rep("D5733-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D5733-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6343_4 <- data.frame(sample=sort(unique(paste(rep("D6343-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D6343-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D6271_4 <- data.frame(sample=sort(unique(paste(rep("D6271-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D6271-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5220_4 <- data.frame(sample=sort(unique(paste(rep("D5220-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D5220-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5417_4 <- data.frame(sample=sort(unique(paste(rep("D5417-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D5417-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))

samples_D5379_4 <- data.frame(sample=sort(unique(paste(rep("D5379-4",3),formatC(barcode_runs[barcode_runs$sample_type=="D5379-4","concentration"], width=2,flag="0"),sep="_"))),conc=c(1/100,1/10000,0))
```

Convert from Illumina read counts to estimates of the number of cells
that were sorted into a bin, and add some other useful information to
our data frame.

``` r
#for each bin, normalize the read counts to the observed ratio of cell recovery among bins
for(i in 1:nrow(barcode_runs)){
  lib <- as.character(barcode_runs$library[i])
  bin <- as.character(barcode_runs$sample[i])
  ratio <- sum(dt[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
  if(ratio<1){ #if there are fewer reads from a FACS bin than cells sorted
    dt[library==lib & sample==bin, count.norm := as.numeric(count)] #don't normalize cell counts, make count.norm the same as count
    print(paste("reads < cells for",lib,bin,", un-normalized (ratio",ratio,")")) #print to console to inform of undersampled bins
  }else{
    dt[library==lib & sample==bin, count.norm := as.numeric(count/ratio)] #normalize read counts by the average read:cell ratio, report in new "count.norm" column
    print(paste("read:cell ratio for",lib,bin,"is",ratio))
  }
}
```

    ## [1] "read:cell ratio for pool1 D5338-3_01_bin1 is 5.14407662601061"
    ## [1] "read:cell ratio for pool1 D5338-3_01_bin2 is 2.4260395586215"
    ## [1] "read:cell ratio for pool1 D5338-3_01_bin3 is 2.16363100324304"
    ## [1] "read:cell ratio for pool1 D5338-3_01_bin4 is 2.91819234725809"
    ## [1] "read:cell ratio for pool1 D5338-3_02_bin1 is 2.22853690702828"
    ## [1] "read:cell ratio for pool1 D5338-3_02_bin2 is 2.23004509568558"
    ## [1] "read:cell ratio for pool1 D5338-3_02_bin3 is 1.96382885585045"
    ## [1] "read:cell ratio for pool1 D5338-3_02_bin4 is 1.8798411122145"
    ## [1] "read:cell ratio for pool1 D5338-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5338-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5338-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5338-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6391-3_01_bin1 is 4.0744907525782"
    ## [1] "read:cell ratio for pool1 D6391-3_01_bin2 is 4.16023186316469"
    ## [1] "read:cell ratio for pool1 D6391-3_01_bin3 is 2.66023918180342"
    ## [1] "read:cell ratio for pool1 D6391-3_01_bin4 is 2.75201632520459"
    ## [1] "read:cell ratio for pool1 D6391-3_02_bin1 is 2.70001351835294"
    ## [1] "read:cell ratio for pool1 D6391-3_02_bin2 is 2.70020948317322"
    ## [1] "read:cell ratio for pool1 D6391-3_02_bin3 is 2.43978842614886"
    ## [1] "read:cell ratio for pool1 D6391-3_02_bin4 is 2.25912799397523"
    ## [1] "read:cell ratio for pool1 D6391-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6391-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6391-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6391-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6404-3_01_bin1 is 4.43490516053831"
    ## [1] "read:cell ratio for pool1 D6404-3_01_bin2 is 2.57414461470387"
    ## [1] "read:cell ratio for pool1 D6404-3_01_bin3 is 2.38335775936145"
    ## [1] "read:cell ratio for pool1 D6404-3_01_bin4 is 1.86778720883187"
    ## [1] "read:cell ratio for pool1 D6404-3_02_bin1 is 2.58149184164969"
    ## [1] "read:cell ratio for pool1 D6404-3_02_bin2 is 2.55027929994478"
    ## [1] "read:cell ratio for pool1 D6404-3_02_bin3 is 2.46006781384776"
    ## [1] "read:cell ratio for pool1 D6404-3_02_bin4 is 2.79188123411734"
    ## [1] "read:cell ratio for pool1 D6404-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6404-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6404-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6404-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5733-3_01_bin1 is 3.06015781017245"
    ## [1] "read:cell ratio for pool1 D5733-3_01_bin2 is 2.89376386373103"
    ## [1] "read:cell ratio for pool1 D5733-3_01_bin3 is 2.68392463866896"
    ## [1] "read:cell ratio for pool1 D5733-3_01_bin4 is 2.68451269635343"
    ## [1] "read:cell ratio for pool1 D5733-3_02_bin1 is 2.20726288098519"
    ## [1] "read:cell ratio for pool1 D5733-3_02_bin2 is 2.46463539941033"
    ## [1] "read:cell ratio for pool1 D5733-3_02_bin3 is 3.3465282590889"
    ## [1] "reads < cells for pool1 D5733-3_02_bin4 , un-normalized (ratio 0.716216216216216 )"
    ## [1] "read:cell ratio for pool1 D5733-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5733-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5733-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5733-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6343-3_01_bin1 is 3.97342229902714"
    ## [1] "read:cell ratio for pool1 D6343-3_01_bin2 is 3.63844562012225"
    ## [1] "read:cell ratio for pool1 D6343-3_01_bin3 is 2.61807097359232"
    ## [1] "read:cell ratio for pool1 D6343-3_01_bin4 is 2.4358712277782"
    ## [1] "read:cell ratio for pool1 D6343-3_02_bin1 is 2.79193452714943"
    ## [1] "read:cell ratio for pool1 D6343-3_02_bin2 is 2.64437857948069"
    ## [1] "read:cell ratio for pool1 D6343-3_02_bin3 is 2.73817029244014"
    ## [1] "read:cell ratio for pool1 D6343-3_02_bin4 is 1.90469208211144"
    ## [1] "read:cell ratio for pool1 D6343-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6343-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6343-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6343-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6271-3_01_bin1 is 4.38861839460401"
    ## [1] "read:cell ratio for pool1 D6271-3_01_bin2 is 3.23181363576321"
    ## [1] "read:cell ratio for pool1 D6271-3_01_bin3 is 2.55321299034832"
    ## [1] "read:cell ratio for pool1 D6271-3_01_bin4 is 2.52402655835786"
    ## [1] "read:cell ratio for pool1 D6271-3_02_bin1 is 2.52945930231849"
    ## [1] "read:cell ratio for pool1 D6271-3_02_bin2 is 2.78267261916355"
    ## [1] "read:cell ratio for pool1 D6271-3_02_bin3 is 2.47740031692377"
    ## [1] "read:cell ratio for pool1 D6271-3_02_bin4 is 5.89067443198159"
    ## [1] "read:cell ratio for pool1 D6271-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6271-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6271-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6271-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5220-3_01_bin1 is 3.30000645545378"
    ## [1] "read:cell ratio for pool1 D5220-3_01_bin2 is 2.60510630598028"
    ## [1] "read:cell ratio for pool1 D5220-3_01_bin3 is 2.1620365525351"
    ## [1] "read:cell ratio for pool1 D5220-3_01_bin4 is 3.69919557638079"
    ## [1] "read:cell ratio for pool1 D5220-3_02_bin1 is 2.78319790807323"
    ## [1] "read:cell ratio for pool1 D5220-3_02_bin2 is 2.95533497936712"
    ## [1] "read:cell ratio for pool1 D5220-3_02_bin3 is 2.42650625279775"
    ## [1] "read:cell ratio for pool1 D5220-3_02_bin4 is 1.39859693877551"
    ## [1] "read:cell ratio for pool1 D5220-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5220-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5220-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5220-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5417-3_01_bin1 is 3.74009357900024"
    ## [1] "read:cell ratio for pool1 D5417-3_01_bin2 is 2.61620950885202"
    ## [1] "read:cell ratio for pool1 D5417-3_01_bin3 is 3.33938501433312"
    ## [1] "read:cell ratio for pool1 D5417-3_01_bin4 is 2.9377023488995"
    ## [1] "read:cell ratio for pool1 D5417-3_02_bin1 is 2.91395659481083"
    ## [1] "read:cell ratio for pool1 D5417-3_02_bin2 is 2.94607169764196"
    ## [1] "read:cell ratio for pool1 D5417-3_02_bin3 is 2.44133049506645"
    ## [1] "read:cell ratio for pool1 D5417-3_02_bin4 is 11.8817188335214"
    ## [1] "read:cell ratio for pool1 D5417-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5417-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5417-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5417-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5379-3_01_bin1 is 3.87536165176223"
    ## [1] "read:cell ratio for pool1 D5379-3_01_bin2 is 3.67441631440995"
    ## [1] "read:cell ratio for pool1 D5379-3_01_bin3 is 2.89292945846248"
    ## [1] "read:cell ratio for pool1 D5379-3_01_bin4 is 2.34539203066172"
    ## [1] "read:cell ratio for pool1 D5379-3_02_bin1 is 3.26801042970922"
    ## [1] "read:cell ratio for pool1 D5379-3_02_bin2 is 2.9314087885025"
    ## [1] "read:cell ratio for pool1 D5379-3_02_bin3 is 2.83912186085839"
    ## [1] "read:cell ratio for pool1 D5379-3_02_bin4 is 1.79285184269107"
    ## [1] "read:cell ratio for pool1 D5379-3_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5379-3_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5379-3_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5379-3_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5338-4_01_bin1 is 11.6389716971206"
    ## [1] "read:cell ratio for pool1 D5338-4_01_bin2 is 6.44946027221054"
    ## [1] "read:cell ratio for pool1 D5338-4_01_bin3 is 2.58675206942105"
    ## [1] "read:cell ratio for pool1 D5338-4_01_bin4 is 2.54923525671529"
    ## [1] "read:cell ratio for pool1 D5338-4_02_bin1 is 2.27041383902171"
    ## [1] "read:cell ratio for pool1 D5338-4_02_bin2 is 2.20925902874143"
    ## [1] "read:cell ratio for pool1 D5338-4_02_bin3 is 2.03616746642606"
    ## [1] "read:cell ratio for pool1 D5338-4_02_bin4 is 2.78067814009029"
    ## [1] "read:cell ratio for pool1 D5338-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5338-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5338-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5338-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6391-4_01_bin1 is 10.8125590991146"
    ## [1] "read:cell ratio for pool1 D6391-4_01_bin2 is 8.52101973281039"
    ## [1] "read:cell ratio for pool1 D6391-4_01_bin3 is 2.89377156692658"
    ## [1] "read:cell ratio for pool1 D6391-4_01_bin4 is 3.64389395440608"
    ## [1] "read:cell ratio for pool1 D6391-4_02_bin1 is 3.16958118988884"
    ## [1] "read:cell ratio for pool1 D6391-4_02_bin2 is 2.78720681713814"
    ## [1] "read:cell ratio for pool1 D6391-4_02_bin3 is 2.50946533406331"
    ## [1] "read:cell ratio for pool1 D6391-4_02_bin4 is 2.49615590034514"
    ## [1] "read:cell ratio for pool1 D6391-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6391-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6391-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6391-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6404-4_01_bin1 is 9.17018806007825"
    ## [1] "read:cell ratio for pool1 D6404-4_01_bin2 is 10.0261999179543"
    ## [1] "read:cell ratio for pool1 D6404-4_01_bin3 is 7.13123210877863"
    ## [1] "read:cell ratio for pool1 D6404-4_01_bin4 is 2.08708958813543"
    ## [1] "read:cell ratio for pool1 D6404-4_02_bin1 is 2.44919843269067"
    ## [1] "read:cell ratio for pool1 D6404-4_02_bin2 is 2.48930026694698"
    ## [1] "read:cell ratio for pool1 D6404-4_02_bin3 is 2.61381871657572"
    ## [1] "read:cell ratio for pool1 D6404-4_02_bin4 is 2.98782984776533"
    ## [1] "read:cell ratio for pool1 D6404-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6404-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6404-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6404-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5733-4_01_bin1 is 8.19648205011562"
    ## [1] "read:cell ratio for pool1 D5733-4_01_bin2 is 5.46669462500874"
    ## [1] "read:cell ratio for pool1 D5733-4_01_bin3 is 4.36844648909227"
    ## [1] "read:cell ratio for pool1 D5733-4_01_bin4 is 2.49909440428023"
    ## [1] "read:cell ratio for pool1 D5733-4_02_bin1 is 2.197057972339"
    ## [1] "read:cell ratio for pool1 D5733-4_02_bin2 is 2.36509942895405"
    ## [1] "read:cell ratio for pool1 D5733-4_02_bin3 is 2.50816728146797"
    ## [1] "read:cell ratio for pool1 D5733-4_02_bin4 is 2.90472367792039"
    ## [1] "read:cell ratio for pool1 D5733-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5733-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5733-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5733-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6343-4_01_bin1 is 9.02855546588863"
    ## [1] "read:cell ratio for pool1 D6343-4_01_bin2 is 7.33686696141479"
    ## [1] "read:cell ratio for pool1 D6343-4_01_bin3 is 3.57402633670598"
    ## [1] "read:cell ratio for pool1 D6343-4_01_bin4 is 2.82135215970721"
    ## [1] "read:cell ratio for pool1 D6343-4_02_bin1 is 2.36345465145332"
    ## [1] "read:cell ratio for pool1 D6343-4_02_bin2 is 2.50195354870943"
    ## [1] "read:cell ratio for pool1 D6343-4_02_bin3 is 2.32245223564545"
    ## [1] "read:cell ratio for pool1 D6343-4_02_bin4 is 2.63242171241094"
    ## [1] "read:cell ratio for pool1 D6343-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6343-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6343-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6343-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D6271-4_01_bin1 is 1.33214479497713"
    ## [1] "read:cell ratio for pool1 D6271-4_01_bin2 is 3.31438870235039"
    ## [1] "read:cell ratio for pool1 D6271-4_01_bin3 is 4.28998149616706"
    ## [1] "read:cell ratio for pool1 D6271-4_01_bin4 is 2.27856555715857"
    ## [1] "read:cell ratio for pool1 D6271-4_02_bin1 is 3.36837281526194"
    ## [1] "read:cell ratio for pool1 D6271-4_02_bin2 is 2.70850823943079"
    ## [1] "read:cell ratio for pool1 D6271-4_02_bin3 is 2.2252301644449"
    ## [1] "read:cell ratio for pool1 D6271-4_02_bin4 is 3.01140276725693"
    ## [1] "read:cell ratio for pool1 D6271-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D6271-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D6271-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D6271-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5220-4_01_bin1 is 5.17605089273162"
    ## [1] "read:cell ratio for pool1 D5220-4_01_bin2 is 3.02066423720391"
    ## [1] "read:cell ratio for pool1 D5220-4_01_bin3 is 3.46360630748473"
    ## [1] "read:cell ratio for pool1 D5220-4_01_bin4 is 2.54653162729434"
    ## [1] "read:cell ratio for pool1 D5220-4_02_bin1 is 2.76709189189189"
    ## [1] "read:cell ratio for pool1 D5220-4_02_bin2 is 2.67533136609709"
    ## [1] "read:cell ratio for pool1 D5220-4_02_bin3 is 3.48980163727702"
    ## [1] "read:cell ratio for pool1 D5220-4_02_bin4 is 2.69194378896738"
    ## [1] "read:cell ratio for pool1 D5220-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5220-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5220-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5220-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5417-4_01_bin1 is 2.39802880970432"
    ## [1] "read:cell ratio for pool1 D5417-4_01_bin2 is 2.26895010189994"
    ## [1] "read:cell ratio for pool1 D5417-4_01_bin3 is 7.54376874829561"
    ## [1] "read:cell ratio for pool1 D5417-4_01_bin4 is 3.49083503509812"
    ## [1] "read:cell ratio for pool1 D5417-4_02_bin1 is 4.23162303885882"
    ## [1] "read:cell ratio for pool1 D5417-4_02_bin2 is 2.43955551553958"
    ## [1] "read:cell ratio for pool1 D5417-4_02_bin3 is 2.45631105521251"
    ## [1] "read:cell ratio for pool1 D5417-4_02_bin4 is 2.98335677344048"
    ## [1] "read:cell ratio for pool1 D5417-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5417-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5417-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5417-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"
    ## [1] "read:cell ratio for pool1 D5379-4_01_bin1 is 4.01666468748763"
    ## [1] "read:cell ratio for pool1 D5379-4_01_bin2 is 2.74770733652313"
    ## [1] "read:cell ratio for pool1 D5379-4_01_bin3 is 4.85336224889121"
    ## [1] "read:cell ratio for pool1 D5379-4_01_bin4 is 3.21527997496256"
    ## [1] "read:cell ratio for pool1 D5379-4_02_bin1 is 2.72350563773398"
    ## [1] "read:cell ratio for pool1 D5379-4_02_bin2 is 2.9805970629205"
    ## [1] "read:cell ratio for pool1 D5379-4_02_bin3 is 3.49916906227904"
    ## [1] "read:cell ratio for pool1 D5379-4_02_bin4 is 3.13718530127202"
    ## [1] "read:cell ratio for pool1 D5379-4_03_bin1 is 2.60447668505432"
    ## [1] "read:cell ratio for pool1 D5379-4_03_bin2 is 2.74591737283407"
    ## [1] "read:cell ratio for pool1 D5379-4_03_bin3 is 1.05392156862745"
    ## [1] "reads < cells for pool1 D5379-4_03_bin4 , un-normalized (ratio 0.525547445255474 )"

``` r
#annotate each barcode as to whether it's a homolog variant, wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
dt[,variant_class:=as.character(NA)]
dt[n_codon_substitutions==0, variant_class := "wildtype"]
dt[n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

#cast the data frame into wide format
dt <- dcast(dt, library + sublibrary + barcode + target + variant_class + aa_substitutions + n_aa_substitutions ~ sample, value.var="count.norm")
```

## Calculating mean bin for each barcode at each sample concentration

Next, for each barcode at each of the sera concentrations, calculate the
“mean bin” response variable. This is calculated as a simple mean, where
the value of each bin is the integer value of the bin (bin1=unbound,
bin4=highly bound) – because of how bins are defined, the mean
fluorescence of cells in each bin are equally spaced on a log-normal
scale, so mean bin correlates with simple mean fluorescence.

We do not use the fluorescence boundaries of the FACS bins in our
calculations here, but we provide them for posterity’s sake below.

\`\`

``` r
#function that returns mean bin and sum of counts for four bins cell counts. Includes cutoffs for bimodal sample splits to filter out
calc.meanbin <- function(vec, split13filter=0.4, split24filter=0.4, split14filter=0.2){
  total <- sum(vec)
  if(is.na(total) | (vec[1] > split13filter*total & vec[3] > split13filter*total) | (vec[2] > split24filter*total & vec[4] > split24filter*total) | (vec[1] > split14filter*total & vec[4] > split14filter*total)){
    return(list(as.numeric(NA),as.numeric(NA)))
  }else{
    return( list(as.numeric((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4])), as.numeric(total)) )
  }
}
  

#iterate through titration samples, compute mean_bin and total_count for each barcode variant
#D5338_3
for(i in 1:nrow(samples_D5338_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5338_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5338_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5338_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5338_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5338_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5338_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6391_3
for(i in 1:nrow(samples_D6391_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6391_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6391_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6391_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6391_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6391_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6391_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6404_3
for(i in 1:nrow(samples_D6404_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6404_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6404_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6404_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6404_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6404_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6404_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5733_3
for(i in 1:nrow(samples_D5733_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5733_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5733_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5733_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5733_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5733_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5733_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6343_3
for(i in 1:nrow(samples_D6343_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6343_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6343_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6343_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6343_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6343_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6343_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6271_3
for(i in 1:nrow(samples_D6271_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6271_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6271_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6271_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6271_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6271_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6271_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5220_3
for(i in 1:nrow(samples_D5220_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5220_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5220_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5220_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5220_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5220_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5220_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5417_3
for(i in 1:nrow(samples_D5417_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5417_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5417_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5417_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5417_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5417_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5417_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5379_3
for(i in 1:nrow(samples_D5379_3)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5379_3[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5379_3[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5379_3[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5379_3[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5379_3[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5379_3[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}



#D5338_4
for(i in 1:nrow(samples_D5338_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5338_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5338_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5338_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5338_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5338_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5338_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6391_4
for(i in 1:nrow(samples_D6391_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6391_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6391_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6391_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6391_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6391_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6391_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6404_4
for(i in 1:nrow(samples_D6404_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6404_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6404_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6404_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6404_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6404_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6404_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5733_4
for(i in 1:nrow(samples_D5733_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5733_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5733_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5733_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5733_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5733_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5733_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6343_4
for(i in 1:nrow(samples_D6343_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6343_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6343_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6343_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6343_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6343_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6343_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D6271_4
for(i in 1:nrow(samples_D6271_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D6271_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D6271_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D6271_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D6271_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D6271_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D6271_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5220_4
for(i in 1:nrow(samples_D5220_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5220_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5220_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5220_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5220_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5220_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5220_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5417_4
for(i in 1:nrow(samples_D5417_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5417_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5417_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5417_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5417_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5417_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5417_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}

#D5379_4
for(i in 1:nrow(samples_D5379_4)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_D5379_4[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_D5379_4[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_D5379_4[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_D5379_4[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_D5379_4[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_D5379_4[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}
```

## Calculate per-bc AUC metrics

We will calculate a simple AUC metric across each barcode’s titration
series. We will also include a minimum cell count that is required for a
meanbin estimate to be used in the titration fit, and a minimum number
of concentrations with determined meanbin that is required for a
titration to be reported.

``` r
#For QC and filtering, output columns giving the average number of cells that were sampled for a barcode across the 9 sample concentrations, and a value for the number of meanbin estimates that were removed for being below the # of cells cutoff
cutoff <- 2

#D5338-3
dt[,`D5338-3_avgcount` := mean(c(`D5338-3_01_totalcount`,`D5338-3_02_totalcount`,`D5338-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5338-3_min_cell_filtered` := sum(c(c(`D5338-3_01_totalcount`,`D5338-3_02_totalcount`,`D5338-3_03_totalcount`)<cutoff,
                                     is.na(c(`D5338-3_01_totalcount`,`D5338-3_02_totalcount`,`D5338-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6391-3
dt[,`D6391-3_avgcount` := mean(c(`D6391-3_01_totalcount`,`D6391-3_02_totalcount`,`D6391-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6391-3_min_cell_filtered` := sum(c(c(`D6391-3_01_totalcount`,`D6391-3_02_totalcount`,`D6391-3_03_totalcount`)<cutoff,
                                     is.na(c(`D6391-3_01_totalcount`,`D6391-3_02_totalcount`,`D6391-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6404-3
dt[,`D6404-3_avgcount` := mean(c(`D6404-3_01_totalcount`,`D6404-3_02_totalcount`,`D6404-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6404-3_min_cell_filtered` := sum(c(c(`D6404-3_01_totalcount`,`D6404-3_02_totalcount`,`D6404-3_03_totalcount`)<cutoff,
                                     is.na(c(`D6404-3_01_totalcount`,`D6404-3_02_totalcount`,`D6404-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5733-3
dt[,`D5733-3_avgcount` := mean(c(`D5733-3_01_totalcount`,`D5733-3_02_totalcount`,`D5733-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5733-3_min_cell_filtered` := sum(c(c(`D5733-3_01_totalcount`,`D5733-3_02_totalcount`,`D5733-3_03_totalcount`)<cutoff,
                                     is.na(c(`D5733-3_01_totalcount`,`D5733-3_02_totalcount`,`D5733-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6343-3
dt[,`D6343-3_avgcount` := mean(c(`D6343-3_01_totalcount`,`D6343-3_02_totalcount`,`D6343-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6343-3_min_cell_filtered` := sum(c(c(`D6343-3_01_totalcount`,`D6343-3_02_totalcount`,`D6343-3_03_totalcount`)<cutoff,
                                     is.na(c(`D6343-3_01_totalcount`,`D6343-3_02_totalcount`,`D6343-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6271-3
dt[,`D6271-3_avgcount` := mean(c(`D6271-3_01_totalcount`,`D6271-3_02_totalcount`,`D6271-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6271-3_min_cell_filtered` := sum(c(c(`D6271-3_01_totalcount`,`D6271-3_02_totalcount`,`D6271-3_03_totalcount`)<cutoff,
                                     is.na(c(`D6271-3_01_totalcount`,`D6271-3_02_totalcount`,`D6271-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5220-3
dt[,`D5220-3_avgcount` := mean(c(`D5220-3_01_totalcount`,`D5220-3_02_totalcount`,`D5220-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5220-3_min_cell_filtered` := sum(c(c(`D5220-3_01_totalcount`,`D5220-3_02_totalcount`,`D5220-3_03_totalcount`)<cutoff,
                                     is.na(c(`D5220-3_01_totalcount`,`D5220-3_02_totalcount`,`D5220-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5417-3
dt[,`D5417-3_avgcount` := mean(c(`D5417-3_01_totalcount`,`D5417-3_02_totalcount`,`D5417-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5417-3_min_cell_filtered` := sum(c(c(`D5417-3_01_totalcount`,`D5417-3_02_totalcount`,`D5417-3_03_totalcount`)<cutoff,
                                     is.na(c(`D5417-3_01_totalcount`,`D5417-3_02_totalcount`,`D5417-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5379-3
dt[,`D5379-3_avgcount` := mean(c(`D5379-3_01_totalcount`,`D5379-3_02_totalcount`,`D5379-3_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5379-3_min_cell_filtered` := sum(c(c(`D5379-3_01_totalcount`,`D5379-3_02_totalcount`,`D5379-3_03_totalcount`)<cutoff,
                                     is.na(c(`D5379-3_01_totalcount`,`D5379-3_02_totalcount`,`D5379-3_03_totalcount`))),na.rm=T),by=c("library","barcode")]
#D5338-4
dt[,`D5338-4_avgcount` := mean(c(`D5338-4_01_totalcount`,`D5338-4_02_totalcount`,`D5338-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5338-4_min_cell_filtered` := sum(c(c(`D5338-4_01_totalcount`,`D5338-4_02_totalcount`,`D5338-4_03_totalcount`)<cutoff,
                                     is.na(c(`D5338-4_01_totalcount`,`D5338-4_02_totalcount`,`D5338-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6391-4
dt[,`D6391-4_avgcount` := mean(c(`D6391-4_01_totalcount`,`D6391-4_02_totalcount`,`D6391-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6391-4_min_cell_filtered` := sum(c(c(`D6391-4_01_totalcount`,`D6391-4_02_totalcount`,`D6391-4_03_totalcount`)<cutoff,
                                     is.na(c(`D6391-4_01_totalcount`,`D6391-4_02_totalcount`,`D6391-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6404-4
dt[,`D6404-4_avgcount` := mean(c(`D6404-4_01_totalcount`,`D6404-4_02_totalcount`,`D6404-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6404-4_min_cell_filtered` := sum(c(c(`D6404-4_01_totalcount`,`D6404-4_02_totalcount`,`D6404-4_03_totalcount`)<cutoff,
                                     is.na(c(`D6404-4_01_totalcount`,`D6404-4_02_totalcount`,`D6404-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5733-4
dt[,`D5733-4_avgcount` := mean(c(`D5733-4_01_totalcount`,`D5733-4_02_totalcount`,`D5733-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5733-4_min_cell_filtered` := sum(c(c(`D5733-4_01_totalcount`,`D5733-4_02_totalcount`,`D5733-4_03_totalcount`)<cutoff,
                                     is.na(c(`D5733-4_01_totalcount`,`D5733-4_02_totalcount`,`D5733-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6343-4
dt[,`D6343-4_avgcount` := mean(c(`D6343-4_01_totalcount`,`D6343-4_02_totalcount`,`D6343-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6343-4_min_cell_filtered` := sum(c(c(`D6343-4_01_totalcount`,`D6343-4_02_totalcount`,`D6343-4_03_totalcount`)<cutoff,
                                     is.na(c(`D6343-4_01_totalcount`,`D6343-4_02_totalcount`,`D6343-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D6271-4
dt[,`D6271-4_avgcount` := mean(c(`D6271-4_01_totalcount`,`D6271-4_02_totalcount`,`D6271-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D6271-4_min_cell_filtered` := sum(c(c(`D6271-4_01_totalcount`,`D6271-4_02_totalcount`,`D6271-4_03_totalcount`)<cutoff,
                                     is.na(c(`D6271-4_01_totalcount`,`D6271-4_02_totalcount`,`D6271-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5220-4
dt[,`D5220-4_avgcount` := mean(c(`D5220-4_01_totalcount`,`D5220-4_02_totalcount`,`D5220-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5220-4_min_cell_filtered` := sum(c(c(`D5220-4_01_totalcount`,`D5220-4_02_totalcount`,`D5220-4_03_totalcount`)<cutoff,
                                     is.na(c(`D5220-4_01_totalcount`,`D5220-4_02_totalcount`,`D5220-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5417-4
dt[,`D5417-4_avgcount` := mean(c(`D5417-4_01_totalcount`,`D5417-4_02_totalcount`,`D5417-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5417-4_min_cell_filtered` := sum(c(c(`D5417-4_01_totalcount`,`D5417-4_02_totalcount`,`D5417-4_03_totalcount`)<cutoff,
                                     is.na(c(`D5417-4_01_totalcount`,`D5417-4_02_totalcount`,`D5417-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]

#D5379-4
dt[,`D5379-4_avgcount` := mean(c(`D5379-4_01_totalcount`,`D5379-4_02_totalcount`,`D5379-4_03_totalcount`),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,`D5379-4_min_cell_filtered` := sum(c(c(`D5379-4_01_totalcount`,`D5379-4_02_totalcount`,`D5379-4_03_totalcount`)<cutoff,
                                     is.na(c(`D5379-4_01_totalcount`,`D5379-4_02_totalcount`,`D5379-4_03_totalcount`))),na.rm=T),by=c("library","barcode")]


#function that calculates an AUC metric across two log-spaced points and a zero-point for substraction, including an option to filter below certain thresholds for average cells across all samples, and number of samples below a cutoff of cells
fit.auc <- function(x.vals,y.vals,zero.val,count.vals,zero.count.val,min.cfu=cutoff){
  if(sum(!is.na(y.vals))==2 & !is.na(zero.val)){
    if(sum(count.vals > min.cfu) == length(count.vals) & zero.count.val > min.cfu){
      y.bg <- y.vals - zero.val
      y.bg[y.bg<0] <- 0
      auc <- sum(diff(rev(x.vals)) * (head(rev(y.bg),-1)+tail(rev(y.bg),-1)))/2 #reverse order, I supply high to low but I want low to high
      return(auc)
    }else{
      return(as.numeric(NA))
    }
  }else{
    return(as.numeric(NA))
  }
}

#fit auc to D5338_3 sera data for each barcode
dt[,c("D5338-3_AUC") := fit.auc(x.vals=log10(samples_D5338_3$conc[1:2]),
                              y.vals=c(`D5338-3_01_meanbin`,`D5338-3_02_meanbin`),
                      zero.val=`D5338-3_03_meanbin`,
                      count.vals=c(`D5338-3_01_totalcount`,`D5338-3_02_totalcount`),
                      zero.count.val=`D5338-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6391_3 sera data for each barcode
dt[,c("D6391-3_AUC") := fit.auc(x.vals=log10(samples_D6391_3$conc[1:2]),
                              y.vals=c(`D6391-3_01_meanbin`,`D6391-3_02_meanbin`),
                      zero.val=`D6391-3_03_meanbin`,
                      count.vals=c(`D6391-3_01_totalcount`,`D6391-3_02_totalcount`),
                      zero.count.val=`D6391-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6404_3 sera data for each barcode
dt[,c("D6404-3_AUC") := fit.auc(x.vals=log10(samples_D6404_3$conc[1:2]),
                              y.vals=c(`D6404-3_01_meanbin`,`D6404-3_02_meanbin`),
                      zero.val=`D6404-3_03_meanbin`,
                      count.vals=c(`D6404-3_01_totalcount`,`D6404-3_02_totalcount`),
                      zero.count.val=`D6404-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5733_3 sera data for each barcode
dt[,c("D5733-3_AUC") := fit.auc(x.vals=log10(samples_D5733_3$conc[1:2]),
                              y.vals=c(`D5733-3_01_meanbin`,`D5733-3_02_meanbin`),
                      zero.val=`D5733-3_03_meanbin`,
                      count.vals=c(`D5733-3_01_totalcount`,`D5733-3_02_totalcount`),
                      zero.count.val=`D5733-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6343_3 sera data for each barcode
dt[,c("D6343-3_AUC") := fit.auc(x.vals=log10(samples_D6343_3$conc[1:2]),
                              y.vals=c(`D6343-3_01_meanbin`,`D6343-3_02_meanbin`),
                      zero.val=`D6343-3_03_meanbin`,
                      count.vals=c(`D6343-3_01_totalcount`,`D6343-3_02_totalcount`),
                      zero.count.val=`D6343-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6271_3 sera data for each barcode
dt[,c("D6271-3_AUC") := fit.auc(x.vals=log10(samples_D6271_3$conc[1:2]),
                              y.vals=c(`D6271-3_01_meanbin`,`D6271-3_02_meanbin`),
                      zero.val=`D6271-3_03_meanbin`,
                      count.vals=c(`D6271-3_01_totalcount`,`D6271-3_02_totalcount`),
                      zero.count.val=`D6271-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5220_3 sera data for each barcode
dt[,c("D5220-3_AUC") := fit.auc(x.vals=log10(samples_D5220_3$conc[1:2]),
                              y.vals=c(`D5220-3_01_meanbin`,`D5220-3_02_meanbin`),
                      zero.val=`D5220-3_03_meanbin`,
                      count.vals=c(`D5220-3_01_totalcount`,`D5220-3_02_totalcount`),
                      zero.count.val=`D5220-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5417_3 sera data for each barcode
dt[,c("D5417-3_AUC") := fit.auc(x.vals=log10(samples_D5417_3$conc[1:2]),
                              y.vals=c(`D5417-3_01_meanbin`,`D5417-3_02_meanbin`),
                      zero.val=`D5417-3_03_meanbin`,
                      count.vals=c(`D5417-3_01_totalcount`,`D5417-3_02_totalcount`),
                      zero.count.val=`D5417-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5379_3 sera data for each barcode
dt[,c("D5379-3_AUC") := fit.auc(x.vals=log10(samples_D5379_3$conc[1:2]),
                              y.vals=c(`D5379-3_01_meanbin`,`D5379-3_02_meanbin`),
                      zero.val=`D5379-3_03_meanbin`,
                      count.vals=c(`D5379-3_01_totalcount`,`D5379-3_02_totalcount`),
                      zero.count.val=`D5379-3_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5338_3 sera data for each barcode
dt[,c("D5338-4_AUC") := fit.auc(x.vals=log10(samples_D5338_3$conc[1:2]),
                              y.vals=c(`D5338-4_01_meanbin`,`D5338-4_02_meanbin`),
                      zero.val=`D5338-4_03_meanbin`,
                      count.vals=c(`D5338-4_01_totalcount`,`D5338-4_02_totalcount`),
                      zero.count.val=`D5338-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6391_3 sera data for each barcode
dt[,c("D6391-4_AUC") := fit.auc(x.vals=log10(samples_D6391_3$conc[1:2]),
                              y.vals=c(`D6391-4_01_meanbin`,`D6391-4_02_meanbin`),
                      zero.val=`D6391-4_03_meanbin`,
                      count.vals=c(`D6391-4_01_totalcount`,`D6391-4_02_totalcount`),
                      zero.count.val=`D6391-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6404_3 sera data for each barcode
dt[,c("D6404-4_AUC") := fit.auc(x.vals=log10(samples_D6404_3$conc[1:2]),
                              y.vals=c(`D6404-4_01_meanbin`,`D6404-4_02_meanbin`),
                      zero.val=`D6404-4_03_meanbin`,
                      count.vals=c(`D6404-4_01_totalcount`,`D6404-4_02_totalcount`),
                      zero.count.val=`D6404-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5733_3 sera data for each barcode
dt[,c("D5733-4_AUC") := fit.auc(x.vals=log10(samples_D5733_3$conc[1:2]),
                              y.vals=c(`D5733-4_01_meanbin`,`D5733-4_02_meanbin`),
                      zero.val=`D5733-4_03_meanbin`,
                      count.vals=c(`D5733-4_01_totalcount`,`D5733-4_02_totalcount`),
                      zero.count.val=`D5733-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6343_3 sera data for each barcode
dt[,c("D6343-4_AUC") := fit.auc(x.vals=log10(samples_D6343_3$conc[1:2]),
                              y.vals=c(`D6343-4_01_meanbin`,`D6343-4_02_meanbin`),
                      zero.val=`D6343-4_03_meanbin`,
                      count.vals=c(`D6343-4_01_totalcount`,`D6343-4_02_totalcount`),
                      zero.count.val=`D6343-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D6271_3 sera data for each barcode
dt[,c("D6271-4_AUC") := fit.auc(x.vals=log10(samples_D6271_3$conc[1:2]),
                              y.vals=c(`D6271-4_01_meanbin`,`D6271-4_02_meanbin`),
                      zero.val=`D6271-4_03_meanbin`,
                      count.vals=c(`D6271-4_01_totalcount`,`D6271-4_02_totalcount`),
                      zero.count.val=`D6271-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5220_3 sera data for each barcode
dt[,c("D5220-4_AUC") := fit.auc(x.vals=log10(samples_D5220_3$conc[1:2]),
                              y.vals=c(`D5220-4_01_meanbin`,`D5220-4_02_meanbin`),
                      zero.val=`D5220-4_03_meanbin`,
                      count.vals=c(`D5220-4_01_totalcount`,`D5220-4_02_totalcount`),
                      zero.count.val=`D5220-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5417_3 sera data for each barcode
dt[,c("D5417-4_AUC") := fit.auc(x.vals=log10(samples_D5417_3$conc[1:2]),
                              y.vals=c(`D5417-4_01_meanbin`,`D5417-4_02_meanbin`),
                      zero.val=`D5417-4_03_meanbin`,
                      count.vals=c(`D5417-4_01_totalcount`,`D5417-4_02_totalcount`),
                      zero.count.val=`D5417-4_03_totalcount`),
   by=c("library","barcode")]

#fit auc to D5379_3 sera data for each barcode
dt[,c("D5379-4_AUC") := fit.auc(x.vals=log10(samples_D5379_3$conc[1:2]),
                              y.vals=c(`D5379-4_01_meanbin`,`D5379-4_02_meanbin`),
                      zero.val=`D5379-4_03_meanbin`,
                      count.vals=c(`D5379-4_01_totalcount`,`D5379-4_02_totalcount`),
                      zero.count.val=`D5379-4_03_totalcount`),
   by=c("library","barcode")]

#save temp data file for downstream troubleshooting
save(dt,file=paste(config$sera_delta_AUC_dir,"/dt.temp.Rda",sep=""))
```

## QC and sanity checks

We will do some QC to make sure we got good titration curves for most of
our library barcodes.

Let’s visualize the AUC binding measurements as violin plots for the
different wildtype targets, for each serum metric.

``` r
p1 <- ggplot(dt[!is.na(`D5338-3_AUC`),],aes(x=variant_class,y=`D5338-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5338-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5338-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5338-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6391-3_AUC`),],aes(x=variant_class,y=`D6391-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6391-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6391-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6391-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6404-3_AUC`),],aes(x=variant_class,y=`D6404-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6404-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6404-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6404-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5733-3_AUC`),],aes(x=variant_class,y=`D5733-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5733-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5733-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5733-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6343-3_AUC`),],aes(x=variant_class,y=`D6343-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6343-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6343-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6343-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6271-3_AUC`),],aes(x=variant_class,y=`D6271-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6271-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6271-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6271-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5220-3_AUC`),],aes(x=variant_class,y=`D5220-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5220-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5220-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5220-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5417-3_AUC`),],aes(x=variant_class,y=`D5417-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5417-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5417-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5417-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5379-3_AUC`),],aes(x=variant_class,y=`D5379-3_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5379-3 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5379-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5379-3.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5338-4_AUC`),],aes(x=variant_class,y=`D5338-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5338-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5338-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5338-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6391-4_AUC`),],aes(x=variant_class,y=`D6391-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6391-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6391-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6391-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6404-4_AUC`),],aes(x=variant_class,y=`D6404-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6404-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6404-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6404-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5733-4_AUC`),],aes(x=variant_class,y=`D5733-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5733-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5733-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5733-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6343-4_AUC`),],aes(x=variant_class,y=`D6343-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6343-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6343-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6343-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D6271-4_AUC`),],aes(x=variant_class,y=`D6271-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6271-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D6271-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D6271-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5220-4_AUC`),],aes(x=variant_class,y=`D5220-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5220-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5220-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5220-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5417-4_AUC`),],aes(x=variant_class,y=`D5417-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5417-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5417-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5417-4.pdf",sep="")))
```

``` r
p1 <- ggplot(dt[!is.na(`D5379-4_AUC`),],aes(x=variant_class,y=`D5379-4_AUC`))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5379-4 sera AUC")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~library+target,nrow=10)

grid.arrange(p1,ncol=1)
```

<img src="compute_AUC_files/figure-gfm/binding_distribution_vioplot_D5379-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/violin-plot_AUC-by-target_D5379-4.pdf",sep="")))
```

## Save barcode-level metrics

In the next script, we will collapse bcs down to final
mutant/variant-level phenotypes, integrate things like expression
effects of variants, and visualize final phenotypes.

``` r
dt[,.(library,sublibrary,barcode,target,variant_class,aa_substitutions,n_aa_substitutions,
     `D5338-3_avgcount`,`D5338-3_AUC`,
     `D6391-3_avgcount`,`D6391-3_AUC`,
     `D6404-3_avgcount`,`D6404-3_AUC`,
     `D5733-3_avgcount`,`D5733-3_AUC`,
     `D6343-3_avgcount`,`D6343-3_AUC`,
     `D6271-3_avgcount`,`D6271-3_AUC`,
     `D5220-3_avgcount`,`D5220-3_AUC`,
     `D5417-3_avgcount`,`D5417-3_AUC`,
     `D5379-3_avgcount`,`D5379-3_AUC`,
     `D5338-4_avgcount`,`D5338-4_AUC`,
     `D6391-4_avgcount`,`D6391-4_AUC`,
     `D6404-4_avgcount`,`D6404-4_AUC`,
     `D5733-4_avgcount`,`D5733-4_AUC`,
     `D6343-4_avgcount`,`D6343-4_AUC`,
     `D6271-4_avgcount`,`D6271-4_AUC`,
     `D5220-4_avgcount`,`D5220-4_AUC`,
     `D5417-4_avgcount`,`D5417-4_AUC`,
     `D5379-4_avgcount`,`D5379-4_AUC`)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$sera_delta_AUC_file, row.names=F)
```

## Plot representative binding curves

Want to illustrate representative binding curves from which AUCs were
measured. Will do, for each of the six sera, binding curves for wildtype
PRD-0038,, SARS-CoV-2 WH1, XBB.1.5, SARS-CoV-1 Urbani

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-3_03_totalcount` > cutoff,`D5338-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-3_02_totalcount` > cutoff,`D5338-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-3_01_totalcount` > cutoff,`D5338-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5338-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5338-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-3_03_totalcount` > cutoff,`D6391-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-3_02_totalcount` > cutoff,`D6391-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-3_01_totalcount` > cutoff,`D6391-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6391-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6391-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-3_03_totalcount` > cutoff,`D6404-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-3_02_totalcount` > cutoff,`D6404-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-3_01_totalcount` > cutoff,`D6404-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6404-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6404-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-3_03_totalcount` > cutoff,`D5733-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-3_02_totalcount` > cutoff,`D5733-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-3_01_totalcount` > cutoff,`D5733-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5733-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5733-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-3_03_totalcount` > cutoff,`D6343-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-3_02_totalcount` > cutoff,`D6343-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-3_01_totalcount` > cutoff,`D6343-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6343-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6343-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-3_03_totalcount` > cutoff,`D6271-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-3_02_totalcount` > cutoff,`D6271-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-3_01_totalcount` > cutoff,`D6271-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6271-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6271-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-3_03_totalcount` > cutoff,`D5220-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-3_02_totalcount` > cutoff,`D5220-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-3_01_totalcount` > cutoff,`D5220-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5220-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5220-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-3_03_totalcount` > cutoff,`D5417-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-3_02_totalcount` > cutoff,`D5417-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-3_01_totalcount` > cutoff,`D5417-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5417-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5417-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-3, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-3, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-3, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-3, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-3, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-3, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-3_03_totalcount` > cutoff,`D5379-3_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-3_02_totalcount` > cutoff,`D5379-3_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-3_01_totalcount` > cutoff,`D5379-3_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5379-3-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5379-3.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5338-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-4_03_totalcount` > cutoff,`D5338-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-4_02_totalcount` > cutoff,`D5338-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5338-4_01_totalcount` > cutoff,`D5338-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5338-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5338-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6391-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-4_03_totalcount` > cutoff,`D6391-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-4_02_totalcount` > cutoff,`D6391-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6391-4_01_totalcount` > cutoff,`D6391-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6391-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6391-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6404-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-4_03_totalcount` > cutoff,`D6404-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-4_02_totalcount` > cutoff,`D6404-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6404-4_01_totalcount` > cutoff,`D6404-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6404-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6404-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5733-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-4_03_totalcount` > cutoff,`D5733-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-4_02_totalcount` > cutoff,`D5733-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5733-4_01_totalcount` > cutoff,`D5733-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5733-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5733-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6343-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-4_03_totalcount` > cutoff,`D6343-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-4_02_totalcount` > cutoff,`D6343-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6343-4_01_totalcount` > cutoff,`D6343-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6343-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6343-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D6271-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-4_03_totalcount` > cutoff,`D6271-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-4_02_totalcount` > cutoff,`D6271-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D6271-4_01_totalcount` > cutoff,`D6271-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D6271-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D6271-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5220-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-4_03_totalcount` > cutoff,`D5220-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-4_02_totalcount` > cutoff,`D5220-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5220-4_01_totalcount` > cutoff,`D5220-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5220-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5220-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5417-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-4_03_totalcount` > cutoff,`D5417-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-4_02_totalcount` > cutoff,`D5417-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5417-4_01_totalcount` > cutoff,`D5417-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5417-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5417-4.pdf",sep="")))
```

``` r
par(mfrow=c(1,6))

#SARS-CoV-2_WH1
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-4, SARS-CoV-2_WH1",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_WH1" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-2_Omicron-XBB15
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-4, SARS-CoV-2_Omicron-XBB15",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-2_Omicron-XBB15" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#PRD-0038
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-4, PRD-0038",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="PRD-0038" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#SARS-CoV-1_Urbani
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-4, SARS-CoV-1_Urbani",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="SARS-CoV-1_Urbani" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RsYN04
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-4, RsYN04",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RsYN04" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))

#RmYN02
#set empty plot window
plot(NULL,NULL,xlim=c(10^-6,10^-2),ylim=c(1,4),log="x",main="D5379-4, RmYN02",ylab="serum binding (mean FACS bin)",xlab="serum dilution")
#put in faint points per-replicate
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`]
points(rep(10^-6,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`]
points(rep(10^-4,length(y)),y,pch=16,col="#7f7f7f02")
y <- dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`]
points(rep(10^-2,length(y)),y,pch=16,col="#7f7f7f02")
#put in black line for the average
y3 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-4_03_totalcount` > cutoff,`D5379-4_03_meanbin`],na.rm=T)
points(10^-6,y3,pch=16,col="black")
y4 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-4_02_totalcount` > cutoff,`D5379-4_02_meanbin`],na.rm=T)
points(10^-4,y4,pch=16,col="black")
y5 <- mean(dt[target=="RmYN02" & sublibrary=="lib61_SARSr-wts" & `D5379-4_01_totalcount` > cutoff,`D5379-4_01_meanbin`],na.rm=T)
points(10^-2,y5,pch=16,col="black")
#connect average points with lines
lines(c(10^-4, 10^-2),c(y4, y5),lwd=1.5,col="black",lty=2)
abline(h=y3,lty=2,col="gray50")
legend("topleft",bty="n",cex=1,legend=paste("mean AUC: "))
```

<img src="compute_AUC_files/figure-gfm/binding_curves_D5379-4-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$sera_delta_AUC_dir,"/representative-plots_D5379-4.pdf",sep="")))
```
