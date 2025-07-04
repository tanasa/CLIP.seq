#!/bin/bash

uropa \
--bed 16h_peaks_intersection_all_4_sets.bed \
--gtf gencode.v48.basic.annotation.gtf \
-o output \
--summary

# other resources : 

# https://bioconductor.org/packages/release/bioc/vignettes/RCAS/inst/doc/RCAS.metaAnalysis.vignette.html
# annotatr : https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html
# regioneR: https://bioconductor.org/packages/devel/bioc/vignettes/regioneR/inst/doc/regioneR.html
# ChIPseeker : https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html
# ChIPpeakAnno : https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html
