Collapse barcodes to final per-RBD/mutant phenotype scores for the
wildtype sarbecoviruses pool
================
Tyler Starr
05/03/2024

- <a href="#setup" id="toc-setup">Setup</a>
- <a href="#calculate-per-variant-mean-scores"
  id="toc-calculate-per-variant-mean-scores">Calculate per-variant mean
  scores</a>
- <a href="#heatmaps" id="toc-heatmaps">Heatmaps!</a>

This notebook reads in the per-barcode sera binding values and
previously measured expression for sarbecovirus homologs pool. It
synthesizes these two sets of results and calculates the final ‘mean’
phenotypes for each variant, and generates some coverage and QC
analyses.

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
if(!file.exists(config$final_variant_scores_dir)){
  dir.create(file.path(config$final_variant_scores_dir))
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

Read in tables of per-barcode AUC

``` r
#read in data, keep just the wildtypes (which should be lib47, but also Omicron_BA2 and the other wt barcodes in the mutant DMS pools)
dt <- data.table(read.csv(config$sera_delta_AUC_file),stringsAsFactors=F)[variant_class=="wildtype" & sublibrary=="lib61_SARSr-wts",]

#check all targets are in the targets_ordered config list
unique(dt$target) %in% config$targets_ordered
```

    ##   [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [91] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
#assign target as a factor in my desired overall plotting order
dt[,target := factor(dt$target,levels=config$targets_ordered)]

#remove substitutiosn columns as these are all widltype
dt[,c("aa_substitutions","n_aa_substitutions") := NULL]

#read in previously measured expression measurements for these variants? Note, don't have this metric for the pool7 additions to this v2 pool
dt_expr <- data.table(read.csv(config$SARSr_lib47_mut_bind_expr),stringsAsFactors=F)
dt_expr[target=="SARS-CoV-2",target:="SARS-CoV-2_WH1"]
dt_expr[target=="SARS-CoV-1_Urbani_HP03L",target:="SARS-CoV-1_Urbani"]
```

## Calculate per-variant mean scores

Unfiltered, look at distribution of AUC scores

First, third draw timepoints

``` r
p1 <- ggplot(dt[!is.na(D5338.3_AUC),],aes(x=target,y=D5338.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5338.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p2 <- ggplot(dt[!is.na(D6391.3_AUC),],aes(x=target,y=D6391.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6391.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p3 <- ggplot(dt[!is.na(D6404.3_AUC),],aes(x=target,y=D6404.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6404.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p4 <- ggplot(dt[!is.na(D5733.3_AUC),],aes(x=target,y=D5733.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5733.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p5 <- ggplot(dt[!is.na(D6343.3_AUC),],aes(x=target,y=D6343.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6343.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p6 <- ggplot(dt[!is.na(D6271.3_AUC),],aes(x=target,y=D6271.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6271.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p7 <- ggplot(dt[!is.na(D5220.3_AUC),],aes(x=target,y=D5220.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5220.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p8 <- ggplot(dt[!is.na(D5417.3_AUC),],aes(x=target,y=D5417.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5417.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p9 <- ggplot(dt[!is.na(D5379.3_AUC),],aes(x=target,y=D5379.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5379.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=1)
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/unfiltered_AUCs_3-1.png" style="display: block; margin: auto;" />

Next, fourth draw time point

``` r
p1 <- ggplot(dt[!is.na(D5338.4_AUC),],aes(x=target,y=D5338.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5338.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p2 <- ggplot(dt[!is.na(D6391.4_AUC),],aes(x=target,y=D6391.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6391.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p3 <- ggplot(dt[!is.na(D6404.4_AUC),],aes(x=target,y=D6404.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6404.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p4 <- ggplot(dt[!is.na(D5733.4_AUC),],aes(x=target,y=D5733.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5733.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p5 <- ggplot(dt[!is.na(D6343.4_AUC),],aes(x=target,y=D6343.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6343.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p6 <- ggplot(dt[!is.na(D6271.4_AUC),],aes(x=target,y=D6271.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6271.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p7 <- ggplot(dt[!is.na(D5220.4_AUC),],aes(x=target,y=D5220.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5220.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p8 <- ggplot(dt[!is.na(D5417.4_AUC),],aes(x=target,y=D5417.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5417.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p9 <- ggplot(dt[!is.na(D5379.4_AUC),],aes(x=target,y=D5379.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5379.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=1)
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/unfiltered_AUCs_4-1.png" style="display: block; margin: auto;" />

Let’s add a variable that flags the top and bottom 2.5% of expression
scores for each variant, and see how violin plots look when censoring
these most extreme 5% of expressed barcodes

``` r
#D5338.3
dt[,D5338.3_censor_lower:=quantile(D5338.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5338.3_censor_upper:=quantile(D5338.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6391.3
dt[,D6391.3_censor_lower:=quantile(D6391.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6391.3_censor_upper:=quantile(D6391.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6404.3
dt[,D6404.3_censor_lower:=quantile(D6404.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6404.3_censor_upper:=quantile(D6404.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5733.3
dt[,D5733.3_censor_lower:=quantile(D5733.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5733.3_censor_upper:=quantile(D5733.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6343.3
dt[,D6343.3_censor_lower:=quantile(D6343.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6343.3_censor_upper:=quantile(D6343.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6271.3
dt[,D6271.3_censor_lower:=quantile(D6271.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6271.3_censor_upper:=quantile(D6271.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5220.3
dt[,D5220.3_censor_lower:=quantile(D5220.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5220.3_censor_upper:=quantile(D5220.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5417.3
dt[,D5417.3_censor_lower:=quantile(D5417.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5417.3_censor_upper:=quantile(D5417.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5379.3
dt[,D5379.3_censor_lower:=quantile(D5379.3_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5379.3_censor_upper:=quantile(D5379.3_AUC,0.975,na.rm=T,type=7),by=c("library","target")]



p1 <- ggplot(dt[!is.na(D5338.3_AUC) & D5338.3_AUC >= D5338.3_censor_lower & D5338.3_AUC <= D5338.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5338.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5338.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p2 <- ggplot(dt[!is.na(D6391.3_AUC) & D6391.3_AUC >= D6391.3_censor_lower & D6391.3_AUC <= D6391.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6391.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6391.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p3 <- ggplot(dt[!is.na(D6404.3_AUC) & D6404.3_AUC >= D6404.3_censor_lower & D6404.3_AUC <= D6404.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6404.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6404.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p4 <- ggplot(dt[!is.na(D5733.3_AUC) & D5733.3_AUC >= D5733.3_censor_lower & D5733.3_AUC <= D5733.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5733.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5733.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p5 <- ggplot(dt[!is.na(D6343.3_AUC) & D6343.3_AUC >= D6343.3_censor_lower & D6343.3_AUC <= D6343.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6343.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6343.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p6 <- ggplot(dt[!is.na(D6271.3_AUC) & D6271.3_AUC >= D6271.3_censor_lower & D6271.3_AUC <= D6271.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6271.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6271.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p7 <- ggplot(dt[!is.na(D5220.3_AUC) & D5220.3_AUC >= D5220.3_censor_lower & D5220.3_AUC <= D5220.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5220.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5220.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p8 <- ggplot(dt[!is.na(D5417.3_AUC) & D5417.3_AUC >= D5417.3_censor_lower & D5417.3_AUC <= D5417.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5417.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5417.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p9 <- ggplot(dt[!is.na(D5379.3_AUC) & D5379.3_AUC >= D5379.3_censor_lower & D5379.3_AUC <= D5379.3_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5379.3_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5379.3 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))


grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=1)
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/censor_2.5_AUCs.3-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_vioplots_AUC-cens_time.3.pdf",sep="")))
```

``` r
#D5338.4
dt[,D5338.4_censor_lower:=quantile(D5338.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5338.4_censor_upper:=quantile(D5338.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6391.4
dt[,D6391.4_censor_lower:=quantile(D6391.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6391.4_censor_upper:=quantile(D6391.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6404.4
dt[,D6404.4_censor_lower:=quantile(D6404.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6404.4_censor_upper:=quantile(D6404.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5733.4
dt[,D5733.4_censor_lower:=quantile(D5733.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5733.4_censor_upper:=quantile(D5733.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6343.4
dt[,D6343.4_censor_lower:=quantile(D6343.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6343.4_censor_upper:=quantile(D6343.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D6271.4
dt[,D6271.4_censor_lower:=quantile(D6271.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D6271.4_censor_upper:=quantile(D6271.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5220.4
dt[,D5220.4_censor_lower:=quantile(D5220.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5220.4_censor_upper:=quantile(D5220.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5417.4
dt[,D5417.4_censor_lower:=quantile(D5417.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5417.4_censor_upper:=quantile(D5417.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]

#D5379.4
dt[,D5379.4_censor_lower:=quantile(D5379.4_AUC,0.025,na.rm=T,type=7),by=c("library","target")]
dt[,D5379.4_censor_upper:=quantile(D5379.4_AUC,0.975,na.rm=T,type=7),by=c("library","target")]



p1 <- ggplot(dt[!is.na(D5338.4_AUC) & D5338.4_AUC >= D5338.4_censor_lower & D5338.4_AUC <= D5338.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5338.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5338.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p2 <- ggplot(dt[!is.na(D6391.4_AUC) & D6391.4_AUC >= D6391.4_censor_lower & D6391.4_AUC <= D6391.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6391.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6391.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p3 <- ggplot(dt[!is.na(D6404.4_AUC) & D6404.4_AUC >= D6404.4_censor_lower & D6404.4_AUC <= D6404.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6404.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6404.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p4 <- ggplot(dt[!is.na(D5733.4_AUC) & D5733.4_AUC >= D5733.4_censor_lower & D5733.4_AUC <= D5733.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5733.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5733.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p5 <- ggplot(dt[!is.na(D6343.4_AUC) & D6343.4_AUC >= D6343.4_censor_lower & D6343.4_AUC <= D6343.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6343.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6343.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p6 <- ggplot(dt[!is.na(D6271.4_AUC) & D6271.4_AUC >= D6271.4_censor_lower & D6271.4_AUC <= D6271.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D6271.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D6271.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p7 <- ggplot(dt[!is.na(D5220.4_AUC) & D5220.4_AUC >= D5220.4_censor_lower & D5220.4_AUC <= D5220.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5220.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5220.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p8 <- ggplot(dt[!is.na(D5417.4_AUC) & D5417.4_AUC >= D5417.4_censor_lower & D5417.4_AUC <= D5417.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5417.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5417.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p9 <- ggplot(dt[!is.na(D5379.4_AUC) & D5379.4_AUC >= D5379.4_censor_lower & D5379.4_AUC <= D5379.4_censor_upper | target %in% config$targets_low_bc,],aes(x=target,y=D5379.4_AUC))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("D5379.4 sera binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))


grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=1)
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/censor_2.5_AUCs.4-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_vioplots_AUC-cens_time.4.pdf",sep="")))
```

Calculate the mean per variant, the standard deviation, and the number
of (post-filter) barcodes on which a variant score was determined

``` r
#apply the censors to NA out the phenotypes outside the range
dt[!(target %in% config$targets_low_bc) & (D5338.3_AUC < D5338.3_censor_lower | D5338.3_AUC > D5338.3_censor_upper), D5338.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6391.3_AUC < D6391.3_censor_lower | D6391.3_AUC > D6391.3_censor_upper), D6391.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6404.3_AUC < D6404.3_censor_lower | D6404.3_AUC > D6404.3_censor_upper), D6404.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5733.3_AUC < D5733.3_censor_lower | D5733.3_AUC > D5733.3_censor_upper), D5733.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6343.3_AUC < D6343.3_censor_lower | D6343.3_AUC > D6343.3_censor_upper), D6343.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6271.3_AUC < D6271.3_censor_lower | D6271.3_AUC > D6271.3_censor_upper), D6271.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5220.3_AUC < D5220.3_censor_lower | D5220.3_AUC > D5220.3_censor_upper), D5220.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5417.3_AUC < D5417.3_censor_lower | D5417.3_AUC > D5417.3_censor_upper), D5417.3_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5379.3_AUC < D5379.3_censor_lower | D5379.3_AUC > D5379.3_censor_upper), D5379.3_AUC:=NA]

dt[!(target %in% config$targets_low_bc) & (D5338.4_AUC < D5338.4_censor_lower | D5338.4_AUC > D5338.4_censor_upper), D5338.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6391.4_AUC < D6391.4_censor_lower | D6391.4_AUC > D6391.4_censor_upper), D6391.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6404.4_AUC < D6404.4_censor_lower | D6404.4_AUC > D6404.4_censor_upper), D6404.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5733.4_AUC < D5733.4_censor_lower | D5733.4_AUC > D5733.4_censor_upper), D5733.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6343.4_AUC < D6343.4_censor_lower | D6343.4_AUC > D6343.4_censor_upper), D6343.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D6271.4_AUC < D6271.4_censor_lower | D6271.4_AUC > D6271.4_censor_upper), D6271.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5220.4_AUC < D5220.4_censor_lower | D5220.4_AUC > D5220.4_censor_upper), D5220.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5417.4_AUC < D5417.4_censor_lower | D5417.4_AUC > D5417.4_censor_upper), D5417.4_AUC:=NA]
dt[!(target %in% config$targets_low_bc) & (D5379.4_AUC < D5379.4_censor_lower | D5379.4_AUC > D5379.4_censor_upper), D5379.4_AUC:=NA]


dt[,mean_D5338.3_AUC:=mean(D5338.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5338.3_AUC:=sd(D5338.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5338.3_AUC:=sum(!is.na(D5338.3_AUC)),by=c("library","target")]

dt[,mean_D6391.3_AUC:=mean(D6391.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6391.3_AUC:=sd(D6391.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6391.3_AUC:=sum(!is.na(D6391.3_AUC)),by=c("library","target")]

dt[,mean_D6404.3_AUC:=mean(D6404.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6404.3_AUC:=sd(D6404.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6404.3_AUC:=sum(!is.na(D6404.3_AUC)),by=c("library","target")]

dt[,mean_D5733.3_AUC:=mean(D5733.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5733.3_AUC:=sd(D5733.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5733.3_AUC:=sum(!is.na(D5733.3_AUC)),by=c("library","target")]

dt[,mean_D6343.3_AUC:=mean(D6343.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6343.3_AUC:=sd(D6343.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6343.3_AUC:=sum(!is.na(D6343.3_AUC)),by=c("library","target")]

dt[,mean_D6271.3_AUC:=mean(D6271.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6271.3_AUC:=sd(D6271.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6271.3_AUC:=sum(!is.na(D6271.3_AUC)),by=c("library","target")]

dt[,mean_D5220.3_AUC:=mean(D5220.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5220.3_AUC:=sd(D5220.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5220.3_AUC:=sum(!is.na(D5220.3_AUC)),by=c("library","target")]

dt[,mean_D5417.3_AUC:=mean(D5417.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5417.3_AUC:=sd(D5417.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5417.3_AUC:=sum(!is.na(D5417.3_AUC)),by=c("library","target")]

dt[,mean_D5379.3_AUC:=mean(D5379.3_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5379.3_AUC:=sd(D5379.3_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5379.3_AUC:=sum(!is.na(D5379.3_AUC)),by=c("library","target")]

dt[,mean_D5338.4_AUC:=mean(D5338.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5338.4_AUC:=sd(D5338.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5338.4_AUC:=sum(!is.na(D5338.4_AUC)),by=c("library","target")]

dt[,mean_D6391.4_AUC:=mean(D6391.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6391.4_AUC:=sd(D6391.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6391.4_AUC:=sum(!is.na(D6391.4_AUC)),by=c("library","target")]

dt[,mean_D6404.4_AUC:=mean(D6404.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6404.4_AUC:=sd(D6404.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6404.4_AUC:=sum(!is.na(D6404.4_AUC)),by=c("library","target")]

dt[,mean_D5733.4_AUC:=mean(D5733.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5733.4_AUC:=sd(D5733.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5733.4_AUC:=sum(!is.na(D5733.4_AUC)),by=c("library","target")]

dt[,mean_D6343.4_AUC:=mean(D6343.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6343.4_AUC:=sd(D6343.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6343.4_AUC:=sum(!is.na(D6343.4_AUC)),by=c("library","target")]

dt[,mean_D6271.4_AUC:=mean(D6271.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D6271.4_AUC:=sd(D6271.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D6271.4_AUC:=sum(!is.na(D6271.4_AUC)),by=c("library","target")]

dt[,mean_D5220.4_AUC:=mean(D5220.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5220.4_AUC:=sd(D5220.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5220.4_AUC:=sum(!is.na(D5220.4_AUC)),by=c("library","target")]

dt[,mean_D5417.4_AUC:=mean(D5417.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5417.4_AUC:=sd(D5417.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5417.4_AUC:=sum(!is.na(D5417.4_AUC)),by=c("library","target")]

dt[,mean_D5379.4_AUC:=mean(D5379.4_AUC,na.rm=T),by=c("library","target")]
dt[,sd_D5379.4_AUC:=sd(D5379.4_AUC,na.rm=T),by=c("library","target")]
dt[,n_bc_D5379.4_AUC:=sum(!is.na(D5379.4_AUC)),by=c("library","target")]
```

Collapse down to tables reporting just the summary statistics for each
genotype.

``` r
dt_final <- dt[,.(library,target,variant_class,
                  mean_D5338.3_AUC, sd_D5338.3_AUC, n_bc_D5338.3_AUC,
                  mean_D6391.3_AUC, sd_D6391.3_AUC, n_bc_D6391.3_AUC,
                  mean_D6404.3_AUC, sd_D6404.3_AUC, n_bc_D6404.3_AUC,
                  mean_D5733.3_AUC, sd_D5733.3_AUC, n_bc_D5733.3_AUC,
                  mean_D6343.3_AUC, sd_D6343.3_AUC, n_bc_D6343.3_AUC,
                  mean_D6271.3_AUC, sd_D6271.3_AUC, n_bc_D6271.3_AUC,
                  mean_D5220.3_AUC, sd_D5220.3_AUC, n_bc_D5220.3_AUC,
                  mean_D5417.3_AUC, sd_D5417.3_AUC, n_bc_D5417.3_AUC,
                  mean_D5379.3_AUC, sd_D5379.3_AUC, n_bc_D5379.3_AUC,
                  mean_D5338.4_AUC, sd_D5338.4_AUC, n_bc_D5338.4_AUC,
                  mean_D6391.4_AUC, sd_D6391.4_AUC, n_bc_D6391.4_AUC,
                  mean_D6404.4_AUC, sd_D6404.4_AUC, n_bc_D6404.4_AUC,
                  mean_D5733.4_AUC, sd_D5733.4_AUC, n_bc_D5733.4_AUC,
                  mean_D6343.4_AUC, sd_D6343.4_AUC, n_bc_D6343.4_AUC,
                  mean_D6271.4_AUC, sd_D6271.4_AUC, n_bc_D6271.4_AUC,
                  mean_D5220.4_AUC, sd_D5220.4_AUC, n_bc_D5220.4_AUC,
                  mean_D5417.4_AUC, sd_D5417.4_AUC, n_bc_D5417.4_AUC,
                  mean_D5379.4_AUC, sd_D5379.4_AUC, n_bc_D5379.4_AUC
                  )]

dt_final <- unique(dt_final); setkey(dt_final, target)
```

Let’s look how SEM is distributed. Can see that SEM is generally very,
very low. Also that it doesn’t really have a relationship with the AUC
metric, which is good.

``` r
par(mfrow=c(3,3))
#D5338.3
x <- dt_final[,mean_D5338.3_AUC]; y <- dt_final[,sd_D5338.3_AUC/sqrt(n_bc_D5338.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5338.3 sera")

#D6391.3
x <- dt_final[,mean_D6391.3_AUC]; y <- dt_final[,sd_D6391.3_AUC/sqrt(n_bc_D6391.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6391.3 sera")

#D6404.3
x <- dt_final[,mean_D6404.3_AUC]; y <- dt_final[,sd_D6404.3_AUC/sqrt(n_bc_D6404.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6404.3 sera")

#D5733.3
x <- dt_final[,mean_D5733.3_AUC]; y <- dt_final[,sd_D5733.3_AUC/sqrt(n_bc_D5733.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5733.3 sera")

#D6343.3
x <- dt_final[,mean_D6343.3_AUC]; y <- dt_final[,sd_D6343.3_AUC/sqrt(n_bc_D6343.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6343.3 sera")

#D6271.3
x <- dt_final[,mean_D6271.3_AUC]; y <- dt_final[,sd_D6271.3_AUC/sqrt(n_bc_D6271.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6271.3 sera")

#D5220.3
x <- dt_final[,mean_D5220.3_AUC]; y <- dt_final[,sd_D5220.3_AUC/sqrt(n_bc_D5220.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5220.3 sera")

#D5417.3
x <- dt_final[,mean_D5417.3_AUC]; y <- dt_final[,sd_D5417.3_AUC/sqrt(n_bc_D5417.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5417.3 sera")

#D5379.3
x <- dt_final[,mean_D5379.3_AUC]; y <- dt_final[,sd_D5379.3_AUC/sqrt(n_bc_D5379.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5379.3 sera")
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/plot_SEMs.3-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/SEM-v-AUC.3.pdf",sep=""),useDingbats=F))
```

``` r
par(mfrow=c(3,3))
#D5338.4
x <- dt_final[,mean_D5338.4_AUC]; y <- dt_final[,sd_D5338.4_AUC/sqrt(n_bc_D5338.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5338.4 sera")

#D6391.4
x <- dt_final[,mean_D6391.4_AUC]; y <- dt_final[,sd_D6391.4_AUC/sqrt(n_bc_D6391.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6391.4 sera")

#D6404.4
x <- dt_final[,mean_D6404.4_AUC]; y <- dt_final[,sd_D6404.4_AUC/sqrt(n_bc_D6404.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6404.4 sera")

#D5733.4
x <- dt_final[,mean_D5733.4_AUC]; y <- dt_final[,sd_D5733.4_AUC/sqrt(n_bc_D5733.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5733.4 sera")

#D6343.4
x <- dt_final[,mean_D6343.4_AUC]; y <- dt_final[,sd_D6343.4_AUC/sqrt(n_bc_D6343.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6343.4 sera")

#D6271.4
x <- dt_final[,mean_D6271.4_AUC]; y <- dt_final[,sd_D6271.4_AUC/sqrt(n_bc_D6271.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D6271.4 sera")

#D5220.4
x <- dt_final[,mean_D5220.4_AUC]; y <- dt_final[,sd_D5220.4_AUC/sqrt(n_bc_D5220.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5220.4 sera")

#D5417.4
x <- dt_final[,mean_D5417.4_AUC]; y <- dt_final[,sd_D5417.4_AUC/sqrt(n_bc_D5417.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5417.4 sera")

#D5379.4
x <- dt_final[,mean_D5379.4_AUC]; y <- dt_final[,sd_D5379.4_AUC/sqrt(n_bc_D5379.4_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="variant AUC",ylab="SEM",main="D5379.4 sera")
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/plot_SEMs.4-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/SEM-v-AUC.4.pdf",sep=""),useDingbats=F))
```

Let’s also look at how standard error of a within-replicate mean varies
with the number of barcodes

``` r
par(mfrow=c(3,3))
#D5338.3
x <- dt_final[,n_bc_D5338.3_AUC]; y <- dt_final[,sd_D5338.3_AUC/sqrt(n_bc_D5338.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D5338.3 sera")

#D6391.3
x <- dt_final[,n_bc_D6391.3_AUC]; y <- dt_final[,sd_D6391.3_AUC/sqrt(n_bc_D6391.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D6391.3 sera")

#D6404.3
x <- dt_final[,n_bc_D6404.3_AUC]; y <- dt_final[,sd_D6404.3_AUC/sqrt(n_bc_D6404.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D6404.3 sera")

#D5733.3
x <- dt_final[,n_bc_D5733.3_AUC]; y <- dt_final[,sd_D5733.3_AUC/sqrt(n_bc_D5733.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D5733.3 sera")

#D6343.3
x <- dt_final[,n_bc_D6343.3_AUC]; y <- dt_final[,sd_D6343.3_AUC/sqrt(n_bc_D6343.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D6343.3 sera")

#D6271.3
x <- dt_final[,n_bc_D6271.3_AUC]; y <- dt_final[,sd_D6271.3_AUC/sqrt(n_bc_D6271.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D6271.3 sera")

#D5220.3
x <- dt_final[,n_bc_D5220.3_AUC]; y <- dt_final[,sd_D5220.3_AUC/sqrt(n_bc_D5220.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D5220.3 sera")

#D5417.3
x <- dt_final[,n_bc_D5417.3_AUC]; y <- dt_final[,sd_D5417.3_AUC/sqrt(n_bc_D5417.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D5417.3 sera")

#D5379.3
x <- dt_final[,n_bc_D5379.3_AUC]; y <- dt_final[,sd_D5379.3_AUC/sqrt(n_bc_D5379.3_AUC)]; plot(x,y,pch=16,col="#00000090",xlab="number barcodes",ylab="SEM",main="D5379.3 sera")
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/plot_sterr_v_n.3-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib47_SEM-v-n-bc.pdf",sep=""),useDingbats=F))
```

Add in the previously-measured expression values.

``` r
dt_final[,expr:=as.numeric(NA)]
for(i in 1:nrow(dt_final)){
  bg <- as.character(dt_final[i,target])
  if(bg %in% dt_expr$target){
    dt_final[i,expr := dt_expr[target==bg,expression]]
  }
}

#compute a delta_expr relative to the median
dt_final[,delta_expr := expr - median(dt_final$expr,na.rm=T)]

#use the normalization values applied to the single-mut DMS data in the corresponding lib40 notebook? The normalizations for each serum are given in that notebook text and manually entered here. We might also imagine using the average across the six coefficients
dt_final[,D5338.3_normAUC := mean_D5338.3_AUC - (0.48648*delta_expr)]
dt_final[,D6391.3_normAUC := mean_D6391.3_AUC - (0.45686*delta_expr)]
dt_final[,D6404.3_normAUC := mean_D6404.3_AUC - (0.53314*delta_expr)]
dt_final[,D5733.3_normAUC := mean_D5733.3_AUC - (0.6619*delta_expr)]
dt_final[,D6343.3_normAUC := mean_D6343.3_AUC - (0.63838*delta_expr)]
dt_final[,D6271.3_normAUC := mean_D6271.3_AUC - (0.48909*delta_expr)]
dt_final[,D5220.3_normAUC := mean_D5220.3_AUC - (0.53003*delta_expr)]
dt_final[,D5417.3_normAUC := mean_D5417.3_AUC - (0.6391*delta_expr)]
dt_final[,D5379.3_normAUC := mean_D5379.3_AUC - (0.44811*delta_expr)]

dt_final[,D5338.4_normAUC := mean_D5338.4_AUC - (0.47022*delta_expr)]
dt_final[,D6391.4_normAUC := mean_D6391.4_AUC - (0.49721*delta_expr)]
dt_final[,D6404.4_normAUC := mean_D6404.4_AUC - (0.43672*delta_expr)]
dt_final[,D5733.4_normAUC := mean_D5733.4_AUC - (0.4110*delta_expr)]
dt_final[,D6343.4_normAUC := mean_D6343.4_AUC - (0.42439*delta_expr)]
dt_final[,D6271.4_normAUC := mean_D6271.4_AUC - (0.42173*delta_expr)]
dt_final[,D5220.4_normAUC := mean_D5220.4_AUC - (0.43962*delta_expr)]
dt_final[,D5417.4_normAUC := mean_D5417.4_AUC - (0.46937*delta_expr)]
dt_final[,D5379.4_normAUC := mean_D5379.4_AUC - (0.33859*delta_expr)]

par(mfrow=c(9,4))
plot(dt_final$delta_expr,dt_final$mean_D5338.3_AUC,xlab="expr relative to median",ylab="AUC, D5338.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5338.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5338.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5338.4_AUC,xlab="expr relative to median",ylab="AUC, D5338.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5338.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5338.4 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6391.3_AUC,xlab="expr relative to median",ylab="AUC, D6391.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6391.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6391.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6391.4_AUC,xlab="expr relative to median",ylab="AUC, D6391.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6391.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6391.4 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6404.3_AUC,xlab="expr relative to median",ylab="AUC, D6404.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6404.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6404.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6404.4_AUC,xlab="expr relative to median",ylab="AUC, D6404.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6404.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6404.4 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5733.3_AUC,xlab="expr relative to median",ylab="AUC, D5733.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5733.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5733.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5733.4_AUC,xlab="expr relative to median",ylab="AUC, D5733.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5733.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5733.4 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6343.3_AUC,xlab="expr relative to median",ylab="AUC, D6343.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6343.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6343.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6343.4_AUC,xlab="expr relative to median",ylab="AUC, D6343.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6343.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6343.4 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6271.3_AUC,xlab="expr relative to median",ylab="AUC, D6271.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6271.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6271.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D6271.4_AUC,xlab="expr relative to median",ylab="AUC, D6271.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D6271.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D6271.4 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5220.3_AUC,xlab="expr relative to median",ylab="AUC, D5220.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5220.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5220.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5220.4_AUC,xlab="expr relative to median",ylab="AUC, D5220.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5220.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5220.4 sera",pch=16)


plot(dt_final$delta_expr,dt_final$mean_D5417.3_AUC,xlab="expr relative to median",ylab="AUC, D5417.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5417.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5417.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5417.4_AUC,xlab="expr relative to median",ylab="AUC, D5417.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5417.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5417.4 sera",pch=16)


plot(dt_final$delta_expr,dt_final$mean_D5379.3_AUC,xlab="expr relative to median",ylab="AUC, D5379.3 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5379.3_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5379.3 sera",pch=16)

plot(dt_final$delta_expr,dt_final$mean_D5379.4_AUC,xlab="expr relative to median",ylab="AUC, D5379.4 sera",pch=16)
plot(dt_final$delta_expr,dt_final$D5379.4_normAUC,xlab="expr relative to median",ylab="expr-normalized AUC, D5379.4 sera",pch=16)
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/add_expr-1.png" style="display: block; margin: auto;" />

Probably don’t need to do expression-normalization

Filter out the two backgrounds that were completely non-expressing. Most
barcodes were purged before the affinity measurements for these
backgrounds, so the affinities are determined from few barcodes and are
just generally unreliable because these are poorly folded/expressing
variants. (E.g. could see very high standard deviations)

``` r
dt_final <- dt_final[!(target %in% c("HKU3-8","AncSARS1a_alt")),]
```

Order factor variables for plotting

``` r
#order target by order given in config
dt_final$target <- factor(dt_final$target,levels=config$targets_ordered)

#rename some columns for convenience
setnames(dt_final,"mean_D5338.3_AUC","D5338.3_AUC")
setnames(dt_final,"mean_D6391.3_AUC","D6391.3_AUC")
setnames(dt_final,"mean_D6404.3_AUC","D6404.3_AUC")
setnames(dt_final,"mean_D5733.3_AUC","D5733.3_AUC")
setnames(dt_final,"mean_D6343.3_AUC","D6343.3_AUC")
setnames(dt_final,"mean_D6271.3_AUC","D6271.3_AUC")
setnames(dt_final,"mean_D5220.3_AUC","D5220.3_AUC")
setnames(dt_final,"mean_D5417.3_AUC","D5417.3_AUC")
setnames(dt_final,"mean_D5379.3_AUC","D5379.3_AUC")

setnames(dt_final,"mean_D5338.4_AUC","D5338.4_AUC")
setnames(dt_final,"mean_D6391.4_AUC","D6391.4_AUC")
setnames(dt_final,"mean_D6404.4_AUC","D6404.4_AUC")
setnames(dt_final,"mean_D5733.4_AUC","D5733.4_AUC")
setnames(dt_final,"mean_D6343.4_AUC","D6343.4_AUC")
setnames(dt_final,"mean_D6271.4_AUC","D6271.4_AUC")
setnames(dt_final,"mean_D5220.4_AUC","D5220.4_AUC")
setnames(dt_final,"mean_D5417.4_AUC","D5417.4_AUC")
setnames(dt_final,"mean_D5379.4_AUC","D5379.4_AUC")
```

## Heatmaps!

Output heatmaps illustrating all wildtype variants with separate columns
for each serum. Do with both the raw AUC and the expression-normalized
metric.

``` r
#make temp long-form data frame
temp1 <- data.table::melt(dt_final[,.(target,
                                      D5338.3_AUC,D5338.4_AUC,
                                      D6391.3_AUC,D6391.4_AUC,
                                      D6404.3_AUC,D6404.4_AUC,
                                      D5733.3_AUC,D5733.4_AUC,
                                      D6343.3_AUC,D6343.4_AUC,
                                      D6271.3_AUC,D6271.4_AUC,
                                      D5220.3_AUC,D5220.4_AUC,
                                      D5417.3_AUC,D5417.4_AUC,
                                      D5379.3_AUC,D5379.4_AUC)],
                          id.vars=c("target"),
                          measure.vars=c("D5338.3_AUC","D5338.4_AUC",
                                      "D6391.3_AUC","D6391.4_AUC",
                                      "D6404.3_AUC","D6404.4_AUC",
                                      "D5733.3_AUC","D5733.4_AUC",
                                      "D6343.3_AUC","D6343.4_AUC",
                                      "D6271.3_AUC","D6271.4_AUC",
                                      "D5220.3_AUC","D5220.4_AUC",
                                      "D5417.3_AUC","D5417.4_AUC",
                                      "D5379.3_AUC","D5379.4_AUC"),
                          variable.name="sera",value.name="AUC")

#reorder mouse sera for display
temp1$sera <- factor(temp1$sera,levels=c("D5379.4_AUC","D5379.3_AUC",
                                      "D5417.4_AUC","D5417.3_AUC",
                                      "D5220.4_AUC","D5220.3_AUC",
                                      "D6271.4_AUC","D6271.3_AUC",
                                      "D6343.4_AUC","D6343.3_AUC",
                                      "D5733.4_AUC","D5733.3_AUC",
                                      "D6404.4_AUC","D6404.3_AUC",
                                      "D6391.4_AUC","D6391.3_AUC",
                                      "D5338.4_AUC","D5338.3_AUC"))

p1 <- ggplot(temp1,aes(target,sera))+geom_tile(aes(fill=AUC),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(0,6),values=c(0,2/6,6/6),na.value="gray40")+ #effective range 2 to 6
  #scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
  #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
  labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

grid.arrange(p1,nrow=1)
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/heatmap_wildtypes_all_AUC-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_heatmap_AUC_all_wildtypes.pdf",sep="")))
```

Showing just extant sarbs.

``` r
#make temp long-form data frame
extant <- c(config$EurAf_extant,config$RsYN04_extant,config$SARS2_extant,config$SARS1_extant,config$Clade2_extant)

temp2 <- temp1[target %in% extant,];temp2$target <- factor(temp2$target,levels=extant)

p1 <- ggplot(temp2,aes(target,sera))+geom_tile(aes(fill=AUC),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(0,6),values=c(0,1/6,6/6),na.value="gray40")+ #effective range 1 to 6
  #scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
  #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
  labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

p1
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/heatmap_wildtypes_extants_AUC-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_heatmap_AUC_extant_wildtypes.pdf",sep="")))
```

Want to average the values across monkeys from the same condition

``` r
#compute averages across conditions

dt_final[,mean_WH1mono_dose3_AUC:=mean(c(D5338.3_AUC,D6391.3_AUC,D6404.3_AUC),na.rm=T),by=c("library","target")]
dt_final[,mean_WH1mono_dose4_AUC:=mean(c(D5338.4_AUC,D6391.4_AUC,D6404.4_AUC),na.rm=T),by=c("library","target")]

dt_final[,mean_tri_dose3_AUC:=mean(c(D5733.3_AUC,D6343.3_AUC,D6271.3_AUC),na.rm=T),by=c("library","target")]
dt_final[,mean_tri_dose4_AUC:=mean(c(D5733.4_AUC,D6343.4_AUC,D6271.4_AUC),na.rm=T),by=c("library","target")]

dt_final[,mean_tetra_dose3_AUC:=mean(c(D5220.3_AUC,D5417.3_AUC,D5379.3_AUC),na.rm=T),by=c("library","target")]
dt_final[,mean_tetra_dose4_AUC:=mean(c(D5220.4_AUC,D5417.4_AUC,D5379.4_AUC),na.rm=T),by=c("library","target")]


#make temp long-form data frame
temp3 <- data.table::melt(dt_final[,.(target,
                                      mean_WH1mono_dose3_AUC, mean_WH1mono_dose4_AUC,
                                      mean_tri_dose3_AUC, mean_tri_dose4_AUC,
                                      mean_tetra_dose3_AUC, mean_tetra_dose4_AUC)],
                          id.vars=c("target"),
                          measure.vars=c("mean_WH1mono_dose3_AUC","mean_WH1mono_dose4_AUC",
                                      "mean_tri_dose3_AUC","mean_tri_dose4_AUC",
                                      "mean_tetra_dose3_AUC","mean_tetra_dose4_AUC"),
                          variable.name="sera",value.name="AUC")

#reorder mouse sera for display
temp3$sera <- factor(temp3$sera,levels=c("mean_tetra_dose4_AUC","mean_tetra_dose3_AUC",
                                         "mean_tri_dose4_AUC","mean_tri_dose3_AUC",
                                         "mean_WH1mono_dose4_AUC","mean_WH1mono_dose3_AUC"))

#make temp long-form data frame
extant <- c(config$EurAf_extant,config$RsYN04_extant,config$SARS2_extant,config$SARS1_extant,config$Clade2_extant)

temp3 <- temp3[target %in% extant,];temp3$target <- factor(temp3$target,levels=extant)

p1 <- ggplot(temp3,aes(target,sera))+geom_tile(aes(fill=AUC),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(1,6),values=c(0,1/5,5/5),na.value="gray40")+ #effective range 2 to 6
  #scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
  #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
  labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

p1
```

<img src="collapse_barcodes_SARSr-wts_files/figure-gfm/heatmap_wildtypes_extants_condition-averaged_AUC-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_heatmap_AUC_extant_wildtypes_condition-averaged.pdf",sep="")))
```

Save output file.

``` r
dt_final %>%
  mutate_if(is.numeric, round, digits=5) %>%
  write.csv(file=config$final_variant_scores_wts_file, row.names=F,quote=F)
```
