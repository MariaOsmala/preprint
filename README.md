# PREPRINT

This package is a PRobabilistic Enhancer PRedictIoN Tool PREPRINT.

## Installation

Clone the project to the desired directory, this will be your path_to_dir

```
path_to_dir=../preprint/
```
The codes are in folder

```
.../preprint/code/
```

## Create directories

In file create_dirs.sh, specify `path_to_dir=.../preprint`. Then run:
```
bash create_dirs.sh
```
## Download data
see Additional File 2.xlsx for links to download files. 
Downloads the data into `preprint/`. File `create_dirs.sh` has instructions where to put the downloaded data.

```
bash download_data.sh
```

## Required softwares and tools

### Linux environment

The codes have been run and tested using the following linux environment:
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

```
module list 
Currently Loaded Modules:
  1) SAMtools/1.3.1-goolf-triton-2016a  23) sqlite/3.23.1         45) libxext/1.3.3
  2) PCRE/8.38-goolf-triton-2016b       24) tcl/8.6.8             46) pixman/0.34.0
  3) GLib/2.48.0-goolf-triton-2016b     25) kbproto/1.0.7         47) cairo/1.14.12-python3
  4) xz/5.2.4                           26) xextproto/7.3.0       48) libjpeg-turbo/1.5.90
  5) openssl/1.0.2k                     27) libxdmcp/1.1.2        49) libtiff/4.0.9
  6) libssh2/1.8.0                      28) libpthread-stubs/0.4  50) icu4c/60.1
  7) curl/7.60.0                        29) xproto/7.0.31         51) libice/1.0.9
  8) pcre/8.42                          30) libxau/1.0.8          52) libsm/1.2.2
  9) libbsd/0.8.6                       31) libxcb/1.13           53) libxt/1.1.5
 10) expat/2.2.5                        32) libx11/1.6.5          54) libxft/2.3.2
 11) gettext/0.19.8.1                   33) tk/8.6.8              55) harfbuzz/1.4.6-python3
 12) ncurses/6.1                        34) python/3.6.3          56) gobject-introspection/1.49.2-python3
 13) readline/7.0                       35) libffi/3.2.1          57) pango/1.41.0-python3
 14) gdbm/1.14.1                        36) glib/2.56.1-python3   58) openblas/0.3.2
 15) perl/5.26.2                        37) jdk/8u181-b13         59) r/3.4.3-python3
 16) libiconv/1.15                      38) libpng/1.6.34         60) matlab/r2012a
 17) git/2.18.0                         39) freetype/2.7.1        61) gmp/6.1.2
 18) GCC/8.2.0-2.31.1                   40) libxml2/2.9.8         62) mpfr/4.0.1
 19) anaconda3/latest                   41) font-util/1.3.1       63) mpc/1.1.0
 20) Boost/1.61.0-iomkl-triton-2017a    42) fontconfig/2.12.3     64) isl/0.19
```

### Required R packages

* Install spp and accelerometry_2.2.5

```
cd $path_to_dir/softwares/
git clone https://github.com/kundajelab/phantompeakqualtools
cd phantompeakqualtools
R
install.packages("snow", repos="http://cran.us.r-project.org")
install.packages("snowfall", repos="http://cran.us.r-project.org")
install.packages("bitops", repos="http://cran.us.r-project.org")
install.packages("caTools", repos="http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools",suppressUpdates=TRUE)
install.packages("./spp_1.14.tar.gz")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/accelerometry/accelerometry_2.2.5.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

* Other required Cran and Bioconductor packages

```
sessionInfo()
R version 3.4.3 (2017-11-30)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /share/apps/spack/software/openblas/0.3.2/7rbfsnk/lib/libopenblasp-r0.3.2.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] gdata_2.18.0                      accelerometry_2.2.5              
 [3] stringr_1.4.0                     reshape2_1.4.3                   
 [5] gridExtra_2.3                     MASS_7.3-51.4                    
 [7] pryr_0.1.4                        ggplot2_3.2.0                    
 [9] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.46.0                  
[11] RColorBrewer_1.1-2                circlize_0.4.6                   
[13] GetoptLong_0.1.7                  doParallel_1.0.14                
[15] iterators_1.0.10                  foreach_1.4.4                    
[17] bitops_1.0-6                      rtracklayer_1.38.3               
[19] ShortRead_1.36.1                  GenomicAlignments_1.14.2         
[21] SummarizedExperiment_1.8.1        DelayedArray_0.4.1               
[23] matrixStats_0.54.0                Biobase_2.38.0                   
[25] Rsamtools_1.30.0                  GenomicRanges_1.30.3             
[27] GenomeInfoDb_1.14.0               Biostrings_2.46.0                
[29] XVector_0.18.0                    IRanges_2.12.0                   
[31] S4Vectors_0.16.0                  BiocParallel_1.12.0              
[33] BiocGenerics_0.24.0               optparse_1.6.2                   

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1             lattice_0.20-35        gtools_3.8.1          
 [4] assertthat_0.2.1       R6_2.4.0               plyr_1.8.4            
 [7] pillar_1.4.1           GlobalOptions_0.1.0    zlibbioc_1.24.0       
[10] rlang_0.4.0            lazyeval_0.2.2         Matrix_1.2-12         
[13] RCurl_1.95-4.12        munsell_0.5.0          compiler_3.4.3        
[16] pkgconfig_2.0.2        shape_1.4.4            tidyselect_0.2.5      
[19] tibble_2.1.3           GenomeInfoDbData_1.0.0 codetools_0.2-15      
[22] XML_3.98-1.20          crayon_1.3.4           dplyr_0.8.1           
[25] withr_2.1.2            grid_3.4.3             gtable_0.3.0          
[28] magrittr_1.5           scales_1.0.0           stringi_1.4.3         
[31] hwriter_1.3.2          getopt_1.20.3          latticeExtra_0.6-28   
[34] rjson_0.2.20           tools_3.4.3            glue_1.3.1            
[37] purrr_0.3.2            colorspace_1.4-1      


```

### bedtools2
* bedtools2 v2.28.0 downloaded from `https://github.com/arq5x/bedtools2`. 
* Extract to `path_to_dir/softwares/`, `cd` to `$path_to_dir/softwares/bedtools2` 
* Run `make`
* Add bedtools to path: `$PATH=$path_to_dir/softwares/bedtools2/bin:$PATH`

### Bowtie2

* Bowtie2 is downloaded as precombiled binaries from https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/
* bowtie2-2.3.3.1 extracted to ```$path_to_dir/softwares/```
* The bowtie2 index: ```wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz```
* Path to index ```$path_to_dir/softwares/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome```

### libsvm
* Download from ```https://github.com/cjlin1/libsvm/archive/v322.zip``` and put here ```$path_to_dir/softwares/``

* ```cd $path_to_dir/softwares/libsvm-3.22/```
* ```make```
* Install gnuplot into ```$path_to_dir/softwares/bin```
* Change the line 19 in file ```$path_to_dir/softwares/libsvm-3.22/tools/easy.py```

```
19    gnuplot_exe = "$path_to_dir/softwares/bin/gnuplot"
```
* ```export PATH=$PATH:$path_to_dir/softwares/libsvm-3.22/```
* ``` export PATH=$path_to_dir/softwares/bin:$PATH``` 


### RFECS
Codes available here:
```
https://github.com/MariaOsmala/RFECS_BUG_FIXES.git

```

### Utils
Download utils from ```http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/```
particularly, ```bedGraphToBigWig```, ```bedToBigBed```, ```fetchChromSizes```, ```hg19.chrom.sizes```
```
mkdir $path_to_dir/softwares/utils
cd $path_to_dir/softwares/utils
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
```
To get the chrom sizes
```
bash fetchChromSizes hg19 > hg19.chrom.sizes
```

## [Preprocessing steps](PreprocessingSteps.md)

## [Model training and prediction, lists of predicted enhancers and instructions to visualize results in genome browser](TrainingPrediction.md)


## Citation:

@Article{Osmala2020,
    AUTHOR = "Maria Osmala and Harri Lähdesmäki",
    TITLE = "Enhancer prediction in the human genome by probabilistic modelling of the chromatin feature patterns",
    VOLUME = {},
    PAGES = "",
    ADDRESS= "",
    ORGANIZATION = "",
    PUBLISHER="",
    MONTH = "",
    YEAR = }



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


## Contact

Maria Osmala, MSc
PhD student
Aalto University School of Science
Department of Computer Science
Email: firstname.surname@aalto.fi
Home Page: https://people.aalto.fi/maria.osmala

Harri Lähdesmäki, D. Sc. (Tech)
Associate Professor
Aalto University School of Science
Department of Computer Science
Email: firstname.surname@aalto.fi
Home Page: http://users.ics.aalto.fi/harrila



===========================================================


## Acknowledgments

* Thank you for Nisha Rajagopal and Göcken Eraslan providing support for using RFECS.

