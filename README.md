# PMC_detect_circRNA

<h2>Overview</h2>

<p>All pipelines and scripts used in the Master's project <b>Biomarker potential of neuroblastoma circRNAs</b> are deposited in this repository.

<p>This repository has two part:</p>
<ul>
    <li>Shell</li>
    <li>R</li>
</ul>

<p>Shell pipelines and scripts are deposited in the bash-folder. R-pipelines and scripts are deposited in the R-folder.<p>

<h2>Software and algorithms</h2>

|**Name**              |**Version**|**Source**                                            |**Identifier**                                                                                                          |
|----------------------|-----------|------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|
|AnnotationHub         |3.0.1      |Morgan & Shepherd, 2021                               |https://bioconductor.org/<br>packages/release/bioc/<br>vignettes/AnnotationHub/<br>inst/doc/AnnotationHub.html          |
|apeglm                |1.14.0     |The University of North Carolina at Chapel Hill       |DOI: 10.18129/B9.bioc.apeglm                                                                                            |
|bedtools              |2.25.0     |Quinlan and Hall, 2010                                |RRID:SCR\_006646                                                                                                        |
|biomaRt               |2.48.3     |Durinck et al., 2009                                  |RRID:SCR\_019214                                                                                                        |
|Burrow-Wheeler Aligner|0.7.17     |Li, 2013                                              |RRID:SCR\_010910                                                                                                        |
|Cairo                 |1.5-12.2   |R Project                                             |http://www.rforge.net/Cairo/                                                                                            |
|circlize              |0.4.13     |Gu, 2014                                              |RRID:SCR\_002141                                                                                                        |
|CIRI2                 |2.0.6      |Gao et al., 2018                                      |https://ciri-cookbook.readthedocs.io/en/latest/                                                                         |
|CIRIquant             |1.1.2      |Zhang et al., 2020                                    |https://ciri-cookbook.readthedocs.io/en/latest/                                                                         |
|ComplexHeatmap        |2.8.0      |Gu, 2016                                              |RRID:SCR\_017270                                                                                                        |
|Cutadapt              |3.4        |National Bioinformatics Infrastructure Sweden         |RRID:SCR\_011841                                                                                                        |
|DESeq2                |1.32.0     |Love et al., 2014                                     |RRID:SCR\_015687                                                                                                        |
|EnhancedVolcano       |1.10.0     |Blighe et al., 2021                                   |RRID:SCR\_018931                                                                                                        |
|ensembldb             |2.16.4     |Rainer et al., 2019                                   |RRID:SCR\_019103                                                                                                        |
|factoextra            |1.0.7      |Alboukadel Kassambara & Fabian Mundt                  |RRID:SCR\_016692                                                                                                        |
|FastQC                |0.11.9     |Brabaham Institute                                    |RRID:SCR\_014583                                                                                                        |
|GenomicFeatures       |1.44.2     |Lawrence et al., 2013                                 |RRID:SCR\_016960                                                                                                        |
|GenomicRanges         |1.44.0     |Lawrence et al., 2013                                 |RRID:SCR\_000025                                                                                                        |
|ggbeeswarm            |0.6.0      |Erik Clarke & Scott Sherril-Mix                       |https://cran.r-project.org/web/packages/ggbeeswarm/<br>vignettes/usageExamples.pdf                                      |
|ggpubr                |0.4.0      |Alboukadel Kassambara                                 |RRID:SCR\_021139                                                                                                        |
|ggrepel               |0.9.1      |Kamil Slowikowski                                     |RRID:SCR\_017393                                                                                                        |
|gridExtra             |2.3        |Baptiste Auguie                                       |https://cran.r-project.org/web/packages/gridExtra/<br>vignettes/arrangeGrob.html                                        |
|HISAT2                |2.2.1      |Johns Hopkins University                              |RRID:SCR\_015530                                                                                                        |
|magick                |2.7.3      |Jeroen Ooms                                           |https://cran.r-project.org/web/packages/magick/vignettes/<br>intro.html                                                 |
|matrixStats           |0.61.0     |https://github.com/<br>HenrikBengtsson/<br>matrixStats|https://cran.rstudio.com/<br>web/packages/matrixStats/<br>vignettes/matrixStats-methods.html                            |
|NbClust               |3.0        |Charrad et al., 2014                                  |https://sites.google.com/<br>site/malikacharrad/<br>research/<br>nbclust-package/                                       |
|R                     |4.1.1      |R Project                                             |RRID:SCR\_001905                                                                                                        |
|RColorBrewer          |1.1-2      |Erich Neuwirth                                        |https://www.r-graph-gallery.com/<br>38-rcolorbrewers-palettes.html                                                      |
|rtracklayer           |1.52.1     |Lawrence et al., 2009                                 |RRID:SCR\_021325                                                                                                        |
|RVenn                 |1.1.0      |Turgut Yigit Akyol                                    |https://cran.r-project.org/web/packages/RVenn/<br>vignettes/vignette.html                                               |
|salmon                |1.4.0      |Patro et al., 2017                                    |https://combine-lab.github.io/salmon/                                                                                   |
|samtools              |1.12       |Wellcome Sanger Institute                             |RRID:SCR\_002105                                                                                                        |
|SeqKit                |0.16.1     |Shen et al., 2016                                     |http://bioinf.shenwei.me/seqkit/                                                                                        |
|StringTie             |2.1.5      |Johns Hopkins University                              |RRID:SCR\_016323                                                                                                        |
|Subread package       |2.0.2      |University of Melbourne                               |RRID:SCR\_009803                                                                                                        |
|SummarizedExperiment  |1.22.0     |Morgan et al., 2021                                   |https://bioconductor.org/packages/release/<br>bioc/vignettes/SummarizedExperiment/inst/<br>doc/SummarizedExperiment.html|
|tidyHeatmap           |1.3.1      |Mangiola et al., 2020                                 |https://cran.r-project.org/<br>web/packages/<br>tidyHeatmap/<br>vignettes/<br>introduction.html                         |
|tidyverse             |1.3.1      |R Studio                                              |RRID:SCR\_019186                                                                                                        |
|Trim Galore           |0.66       |Brabaham Institute                                    |RRID:SCR\_011847                                                                                                        |
|tximport              |1.20.0     |Soneson et al., 2015                                  |RRID:SCR\_016752                                                                                                        |
|tximportData          |1.20.0     |Love, 2021                                            |DOI: 10.18129/B9.bioc.tximportData                                                                                      |
|UpSetR                |1.4.0      |Conway et al., 2017                                   |https://cran.r-project.org/web/packages/UpSetR/<br>vignettes/basic.usage.html                                           |
|viridis               |0.6.2      |https://github.com/<br>sjmgarnier/viridis/            |RRID:SCR\_016696                                                                                                        |
<h2>Contact</h2>
<p>For questions don't hesitate to send me a message at:
<br>
    <a href="mailto:q.c.lin@students.uu.nl">q.c.lin@students.uu.nl</a>
<br>or<br>
    <a href="mailto:q.c.lin-2@prinsesmaximacentrum.nl">q.c.lin-2@prinsesmaximacentrum.nl</a> (Non-active from 21/01/2022)
</p>