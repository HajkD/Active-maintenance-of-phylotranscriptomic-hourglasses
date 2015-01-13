# Reproducible Scripts for the Publication: 

Drost HG, Gabel A, Grosse I, Quint M (2015). __Evidence for active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis__ Mol. Biol. Evol. (In Review)

## Performing Phylostratigraphy

__Phylostratigraphy__ was introduced by <a href="http://www.sciencedirect.com/science/article/pii/S0168952507002995">Domazet-Lo&scaron;o et al. in 2007</a> to trace the evolutionary origin of protein coding genes. It was performed by using the Perl script `createPsMap.pl`. The resulting phylostratigraphic map stores the phylostratum in the first column and the corresponding gene id in the second column.

For creating the phylostratigraphic map the following steps have to be done:
    
1) Make sure that BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.21/) is installed on your machine.

2) Download the sequence database <a href="http://msbi.ipb-halle.de/download/phyloBlastDB_Drost_Gabel_Grosse_Quint.fa.tbz">phyloBlastDB_Drost_Gabel_Grosse_Quint.fa</a> used for BLAST searches and unpack it (`tar xfvj phyloBlastDB_Drost_Gabel_Grosse_Quint.fa.tbz`).

3) Make sure that the header of your FASTA-files (e.g. Athaliana_167_protein.fa) fullfills the following specification:<br />
  <code>>GeneID | [organism_name] | [taxonomy]</code><br />
  Notice, the taxonomy begins after the node "Cellular organisms" e.g.
```{terminal}
>NP_146894.1 | [Aeropyrum pernix] | [Archaea; Crenarchaeota; Thermoprotei; Desulfurococcales; Desulfurococcaceae; Aeropyrum]
or
>YP_001514406.1 | [Acaryochloris marina MBIC11017] | [Bacteria; Cyanobacteria; Oscillatoriophycideae; Chroococcales; Acaryochloris; Acaryochloris marina]
or
>ATCG00500.1|PACid:19637947 | [Arabidopsis thaliana] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis]
```

4) Use the following command to start the Perl script
```terminal
perl createPsMap.pl -i Athaliana_167_protein_with_new_Header.fa -d phyloBlastDB_Drost_Gabel_Grosse_Quint.fa -p BLAST_Athaliana 
                    -r athaliana_blast_results -t 30 -a 64             
Arguments:
-i,--input          input file name in FASTA format
-d,--database       BLAST sequence database name
-p,--prefix         Prefix for generated mysql-files containing BLAST results
-r,--resultTable    mysql table name
-t,--threshold      threshold for sequence length (Default 30 amino acids)
-a                  threads for BLAST searches
-e,--evalue         e-value threshold for BLAST 
```

## Performing Divergence Stratigraphy

__Divergence Stratigraphy__ is the process of quantifying the selection pressure (in terms of amino acid sequence divergence) acting on orthologous genes between closely related species. The resulting sequence divergence map (short divergence map), stores the divergence stratum in the first column and the query_id of inferred orthologous genes in the second column ([Quint et al., 2012](http://www.nature.com/nature/journal/v490/n7418/full/nature11394.html); Drost et al., 2015).

Following steps are performed to obtain a standard divergence map based on the `divergence_stratigraphy()` function:

1) Orthology Inference using BLAST reciprocal best hit ("RBH") based on blastp

2) Pairwise global amino acid alignments of orthologous genes using the [Needleman-Wunsch algorithm](http://www.sciencedirect.com/science/article/pii/0022283670900574)

3) Codon alignments of orthologous genes using [PAL2NAL](http://www.bork.embl.de/pal2nal/)

4) dNdS estimation using [Comeron's method (1995)](http://link.springer.com/article/10.1007/BF00173196)

5) Assigning estimated dNdS values to divergence strata (deciles of all dNdS values)

When using the `divergence_stratigraphy()` function implemented in `orthologr` it is assumed that you have BLAST installed on your machine.


```r

# install package 'orthologr' from: https://github.com/HajkD/orthologr
install.packages("devtools") # note for wondows installation see https://github.com/HajkD/orthologr for details
devtools::install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)
library(orthologr)

install.packages("myTAI")
library(myTAI)

install.packages("gdata")
library(gdata)

# install the Biostrings package from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")

```


### Retrieving Divergence Maps for D. rerio, D. melanogaster, and A. thaliana

To perform __Divergence Stratigraphy__ using `orthologr` you need the following prerequisites

* a CDS file covering all protein coding genes of the query organism of interest
* a CDS file covering all protein coding genes of the subject organism of interest

First, the CDS sequences for all protein coding genes of all query and subject organisms need to
be obtained.

The CDS retrieval can be done using a `Terminal` or by manual downloading the files

### For D. rerio
* ftp://ftp.ensembl.org/pub/release-77/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/takifugu_rubripes/cds/Takifugu_rubripes.FUGU4.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/astyanax_mexicanus/cds/Astyanax_mexicanus.AstMex102.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/xiphophorus_maculatus/cds/Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/gadus_morhua/cds/Gadus_morhua.gadMor1.cds.all.fa.gz


### For D. melanogaster
* ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.53_FB2013_05/fasta/dmel-all-CDS-r5.53.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r1.4_FB2012_03/fasta/dsim-all-CDS-r1.4.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_yakuba/dyak_r1.3_FB2011_08/fasta/dyak-all-CDS-r1.3.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_persimilis/dper_r1.3_FB2014_03/fasta/dper-all-CDS-r1.3.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_r1.2_FB2012_01/fasta/dvir-all-CDS-r1.2.fasta.gz


### For A. thaliana 
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Athaliana/annotation/Athaliana_167_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Alyrata/annotation/Alyrata_107_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Cpapaya/annotation/Cpapaya_113_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Crubella/annotation/Crubella_183_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Brapa/annotation/Brapa_197_cds.fa.gz

This is an example shell script how to download the CDS files listed above:

```shell

# download CDS file of A. thaliana
curl ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Athaliana/annotation/Athaliana_167_cds.fa.gz -o Athaliana_167_cds.fa.gz

# unzip the fasta file
gunzip -d Athaliana_167_cds.fa.gz

# download CDS file of A. lyrata

curl ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Alyrata/annotation/Alyrata_107_cds.fa.gz
-o Alyrata_107_cds.fa.gz

# unzip the fasta file
gunzip -d Alyrata_107_cds.fa.gz

```
When the download is finished you need to unzip the files and
then start R to perform the following analyses. Make sure you
have all fasta files stored in the current working directory of
your R session.

When using the `divergence_stratigraphy()` function, you can specify the `comp_cores`
argument to perform parallel computations on a multicore machine.

```r

library(orthologr)

### For Danio rerio

# compute the divergence map of D. rerio vs. T. rubripes
Drerio_vs_Trubripes_DM <- divergence_stratigraphy(
                         query_file      = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file    = "Takifugu_rubripes.FUGU4.cds.all.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         



# compute the divergence map of D. rerio vs. A. mexicanus
Drerio_vs_Amexicanus_DM <- divergence_stratigraphy(
                         query_file      = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file    = "Astyanax_mexicanus.AstMex102.cds.all.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         
                         
                         

# compute the divergence map of D. rerio vs. X. maculatus
Drerio_vs_Xmaculatus_DM <- divergence_stratigraphy(
                         query_file      = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file    = "Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )




# compute the divergence map of D. rerio vs. G. morhua
Drerio_vs_Gmorhua_DM <- divergence_stratigraphy(
                         query_file      = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file    = "Gadus_morhua.gadMor1.cds.all.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )




### For Drosophila melanogaster

# compute the divergence map of D. melanogaster vs. D. simulans
Dmel_vs_Dsim_DM <- divergence_stratigraphy(
                         query_file      = "dmel-all-CDS-r5.53.fasta",
                         subject_file    = "dsim-all-CDS-r1.4.fasta",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         


# compute the divergence map of D. melanogaster vs. D. persimilis
Dmel_vs_Dper_DM <- divergence_stratigraphy(
                         query_file      = "dmel-all-CDS-r5.53.fasta",
                         subject_file    = "dper-all-CDS-r1.3.fasta",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         
                         

# compute the divergence map of D. melanogaster vs. D. yakuba
Dmel_vs_Dyak_DM <- divergence_stratigraphy(
                         query_file      = "dmel-all-CDS-r5.53.fasta",
                         subject_file    = "dyak-all-CDS-r1.3.fasta",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )



# compute the divergence map of D. melanogaster vs. D. virilis
Dmel_vs_Dvir_DM <- divergence_stratigraphy(
                         query_file      = "dmel-all-CDS-r5.53.fasta",
                         subject_file    = "dvir-all-CDS-r1.2.fasta",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         


### For Arabidopsis thaliana

# compute the divergence map of A. thaliana vs. A. lyrata
Ath_vs_Aly_DM <- divergence_stratigraphy(
                         query_file      = "Athaliana_167_cds.fa",
                         subject_file    = "Alyrata_107_cds.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         
  
  
# compute the divergence map of A. thaliana vs. C. rubella
Ath_vs_Crubella_DM <- divergence_stratigraphy(
                         query_file      = "Athaliana_167_cds.fa",
                         subject_file    = "Crubella_183_cds.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         


# compute the divergence map of A. thaliana vs. B. rapa
Ath_vs_Brapa_DM <- divergence_stratigraphy(
                         query_file      = "Athaliana_167_cds.fa",
                         subject_file    = "Brapa_197_cds.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )
                         
                         

# compute the divergence map of A. thaliana vs. C. papaya
Ath_vs_Cpapaya_DM <- divergence_stratigraphy(
                         query_file      = "Athaliana_167_cds.fa",
                         subject_file    = "Cpapaya_113_cds.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )

```

## Mapping Gene IDs

It is now assumed that the Divergence Map of interest and the corresponding gene expression data set
are joined. For this purpose the `MatchMap()` function implemented in the `myTAI` package can be used.
See `?myTAI::MatchMap` for details. However, the `MatchMap()` function can only deal with identical gene ids
present in the Phylo/Divergence-Maps and the corresponding gene expression set. 


## Reading PhyloExpressionSets and DivergenceExpressionSets

After performing __Phylostratigraphy__ and __Divergence Stratigraphy__ PhyloExpressionSets and
DivergenceExpressionSets can be obtained by matching the corresponding _phylostratigraphic maps_ and _divergence maps_ of _Danio rerio_, _Drosophila melanogaster_, and _Arabidopsis thaliana_ with the corresponding transcriptome data sets covering
the embryogenesis of _Danio rerio_, _Drosophila melanogaster_, and _Arabidopsis thaliana_ (see Methods in Drost et al., 2015 for details)


The following data sets (Supplementary table S3 and S5) can be downloaded here : [Supplementary table S3.xls](http://figshare.com/articles/Supplementary_table_S3/1244948) and [Supplementary table S5.xls](http://figshare.com/articles/Supplementary_table_S5/1244950).

```r

library(gdata)

### read data sets

## read PhyloExpressionSets
Drerio_PhyloExpressionSet <- read.xls("Supplementary table S3.xls",sheet = 1)
Dmelanogaster_PhyloExpressionSet <- read.xls("Supplementary table S3.xls",sheet = 2)
Athaliana_PhyloExpressionSet <- read.xls("Supplementary table S3.xls",sheet = 3)


## read DivergenceExpressionSets

# Danio rerio

# D. rerio vs. A. mexicanus
Drerio_vs_Amex_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 9)
# D. rerio vs. F. rubripes
Drerio_vs_Frubripes_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 10)
# D. rerio vs. X. maculatus
Drerio_vs_Xmac_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 11)
# D. rerio vs. G. morhua
Drerio_vs_Gmor_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 12)



# Drosophila melanogaster

# D. melanogaster vs. D. simulans
Dmel_Dsim_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 1)
# D. melanogaster vs. D. yakuba
Dmel_Dyak_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 2)
# D. melanogaster vs. D. persimilis
Dmel_Dper_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 3)
# D. melanogaster vs. D. virilis
Dmel_Dvir_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 4)



# Arabidopsis thaliana
# A. thaliana vs A. lyrata
Ath_Aly_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 5)
# A thaliana vs. T. halophila
Ath_Brapa_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 8)
# A thaliana vs. C. rubella
Ath_Crub_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 7)
# A thaliana vs. C. papaya 
Ath_Cpapaya_DivergenceExpressionSet <- read.xls("Supplementary table S5.xls",sheet = 6)

```

## Generating Figures

First load the following packages into the work space:

```r

library(myTAI)

```


### Figure 2

```r
svg("Fig2.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(8,1,0))

PlotPattern(Drerio_PhyloExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:18, mid = 19:36, late = 37:40),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TAI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Dmelanogaster_PhyloExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TAI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Athaliana_PhyloExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TAI", mgp = c(3,0.5,0), cex.lab = 1.5)

dev.off()

```

### Figure 3

```r
svg("Fig3.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(8,1,0))

PlotCorrelation(Drerio_PhyloExpressionSet, Drerio_vs_Amex_DivergenceExpressionSet , 
                method = "kendall", main.text = "")

PlotCorrelation(Dmelanogaster_PhyloExpressionSet, Dmel_Dsim_DivergenceExpressionSet, 
                method = "kendall", main.text = "")

PlotCorrelation(Athaliana_PhyloExpressionSet, Ath_Aly_DivergenceExpressionSet, 
                method = "kendall", main.text = "")

dev.off()
```


### Figure 4

```r
svg("Fig4.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(8,1,0))

PlotPattern(Drerio_vs_Amex_DivergenceExpressionSet[ , 1:42] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:18, mid = 19:36, late = 37:40),
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Drerio_vs_Frubripes_DivergenceExpressionSet : ",nrow(Drerio_vs_Frubripes_DivergenceExpressionSet), " genes.")
cat("\n")

PlotPattern(Dmel_Dsim_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Dmel_Dsim_DivergenceExpressionSet : ",nrow(Dmel_Dsim_DivergenceExpressionSet), " genes.")
cat("\n")

PlotPattern(Ath_Aly_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Ath_Aly_DivergenceExpressionSet : ",nrow(Ath_Aly_DivergenceExpressionSet), " genes.")
cat("\n")

dev.off()

```

## Supplementary Figures


### Suppl_Figure 2

```r
svg("S2.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.5))
par(mai = c(1.4,0.6,0.5,0.5))
par(mgp = c(5,1,0))

PlotPattern(Drerio_vs_Frubripes_DivergenceExpressionSet[ , 1:42] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:18, mid = 19:36, late = 37:40),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", main = "D. rerio vs T. rubripes", las = 3, cex.lab = 1.5, cex.axis = 1)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(2,0.5,0), cex.lab = 1.5)

cat(paste0("Drerio_vs_Amex_DivergenceExpressionSet : ",nrow(Drerio_vs_Amex_DivergenceExpressionSet), " genes."))
cat("\n")


PlotPattern(Drerio_vs_Xmac_DivergenceExpressionSet[ , 1:42] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:18, mid = 19:36, late = 37:40),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", main = "D. rerio vs X. maculatus", las = 3, cex.lab = 1.5, cex.axis = 1)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(2,0.5,0), cex.lab = 1.5)

cat(paste0("Drerio_vs_Xmac_DivergenceExpressionSet : ",nrow(Drerio_vs_Xmac_DivergenceExpressionSet), " genes."))
cat("\n")


PlotPattern(Drerio_vs_Gmor_DivergenceExpressionSet[ , 1:42] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:18, mid = 19:36, late = 37:40),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", main = "D. rerio vs G. morhua", las = 3, cex.lab = 1.4, cex.axis = 1)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("D")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(2,0.5,0), cex.lab = 1.5)

cat(paste0("Drerio_vs_Gmor_DivergenceExpressionSet : ",nrow(Drerio_vs_Gmor_DivergenceExpressionSet), " genes."))

dev.off()



```


### Suppl_Figure 3

```r
svg("S3.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(8,1,0))

PlotPattern(Dmel_Dyakuba_DivergenceExpressionSet[ , 1:14] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", main = "D. melanogaster vs D. yakuba",las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Dmel_Dyakuba_DivergenceExpressionSet : ",nrow(Dmel_Dyakuba_DivergenceExpressionSet), " genes."))
cat("\n")


PlotPattern(Dmel_Dper_DivergenceExpressionSet[ , 1:14] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", main = "D. melanogaster vs D. persimilis", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Dmel_Dper_DivergenceExpressionSet : ",nrow(Dmel_Dper_DivergenceExpressionSet), " genes."))
cat("\n")

PlotPattern(Dmel_Dvir_DivergenceExpressionSet[ , 1:14] , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", main = "D. melanogaster vs D. virilis", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Dmel_Dvir_DivergenceExpressionSet : ",nrow(Dmel_Dvir_DivergenceExpressionSet), " genes."))
cat("\n")

dev.off()
```

### Suppl_Figure 4

```r
svg("S4.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(8,1,0))

PlotPattern(Ath_Crub_DivergenceExpressionSet , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", main = "A. thaliana vs C. rubella", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Ath_Crub_DivergenceExpressionSet : ",nrow(Ath_Crub_DivergenceExpressionSet), " genes."))
cat("\n")



PlotPattern(Ath_Brapa_DivergenceExpressionSet , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", main = "A. thaliana vs B. rapa", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Ath_Brapa_DivergenceExpressionSet : ",nrow(Ath_Brapa_DivergenceExpressionSet), " genes."))
cat("\n")


PlotPattern(Ath_Cpapaya_DivergenceExpressionSet , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", main = "A. thaliana vs C. papaya", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

cat(paste0("Ath_Cpapaya_DivergenceExpressionSet : ",nrow(Ath_Cpapaya_DivergenceExpressionSet), " genes."))
cat("\n")


dev.off()

```


### Suppl_Figure 5

```r

svg("S5.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(3,1,0))

PlotCorrelation(Drerio_PhyloExpressionSet, Drerio_vs_Frubripes_DivergenceExpressionSet , 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

PlotCorrelation(Drerio_PhyloExpressionSet, Drerio_vs_Xmac_DivergenceExpressionSet, 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

PlotCorrelation(Drerio_PhyloExpressionSet, Drerio_vs_Gmor_DivergenceExpressionSet, 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

dev.off()

```

### Suppl_Figure 6

```r
svg("S6.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(3,1,0))

PlotCorrelation(Dmelanogaster_PhyloExpressionSet, Dmel_Dyakuba_DivergenceExpressionSet , 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

PlotCorrelation(Dmelanogaster_PhyloExpressionSet, Dmel_Dper_DivergenceExpressionSet, 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

PlotCorrelation(Dmelanogaster_PhyloExpressionSet, Dmel_Dvir_DivergenceExpressionSet, 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

dev.off()

```


### Suppl_Figure 7

```r
svg("S7.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(3,1,0))

PlotCorrelation(Athaliana_PhyloExpressionSet, Ath_Crub_DivergenceExpressionSet , 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

PlotCorrelation(Athaliana_PhyloExpressionSet, Ath_Brapa_DivergenceExpressionSet, 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

PlotCorrelation(Athaliana_PhyloExpressionSet, Ath_Cpapaya_DivergenceExpressionSet, 
                method = "kendall", main.text = "",cex.axis = 1.7,cex.lab = 1.7,xlab = "Phylostratum",
                ylab = "Divergencestratum",cex.main = 1.7)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.7,inset = c(-0.08,-0.15))
box()

dev.off()

```

### Suppl_Figure 8

```r
# read embryo defective genes
Ath <- read.csv("EmbryoDefective_Expression.csv", sep = ";", header = TRUE)

Ath_MeanProfile <- colMeans(Ath[ , 2:8])

std_error <- function(x){sd(x)/sqrt(length(x))}

svg("S8A.svg",width = 8,height = 6)

plot(Ath_MeanProfile,
    type = "l", lwd = 6,cex.axis = 1.3, cex.lab = 1.3, ylim = c(2600,4600),
    xaxt = "n",xlab = "Ontogeny", ylab = "Expression Level")

axis(1, 1:7, names(Ath[ , 2:8]),cex.axis = 1.3, cex.lab = 1.3)

lines(Ath_MeanProfile + apply(Ath[ , 2:8],2,std_error), col = "lightgrey", lwd = 6)
lines(Ath_MeanProfile - apply(Ath[ , 2:8],2,std_error), col = "lightgrey", lwd = 6)
legend("bottom",legend = c("Mean EDG Expression","std. error"), fill = c("black","lightgrey"), bty = "n", cex = 2)
par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 2,inset = c(-0.08,-0.15))

dev.off()


```

Performing Statistical Tests

```r
# install.packages(dunn.test)
library("dunn.test")

# perform a Kruskal-Wallis Rank Sum Test
kruskal.test(Ath[ , 2:8])

# perform Dunn's test of multiple comparisons using rank sums
dunn.test(Ath[ , 2:8], method = "bh")

```


### Suppl_Figure 9

```r
svg("S8.svg",width = 8,height = 5)
PlotPattern(Drerio_vs_Amex_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest",
            modules = list(early = 1:18, mid = 19:36, late = 37:61), shaded.area = TRUE, 
            type = "l", lwd = 9,xlab = "Ontogeny", ylab = "TDI", main = "D. rerio vs A. mexicanus")

dev.off()
```


### Suppl_Figure 10

This figure shall illustrate how the test statistic for the `Reductive Hourglass Test` is build.
Here we only show an example for `Athaliana_PhyloExpressionSet`, but the same procedure can analogously be
applied to all other `PhyloExpressionSets` and `DivergenceExpressionSets`.

```r
library(myTAI)

rhScore_real <- rhScore(TAI(Athaliana_PhyloExpressionSet),
                       method = "min",scoringMethod = "mean-mean", 
                       early = 1:2, mid = 3:5, late = 6:7)
                       
rhScores <- apply(bootMatrix(Athaliana_PhyloExpressionSet, 10000), 1 , rhScore,
                   method = "min",scoringMethod = "mean-mean",
                   early = 1:2, mid = 3:5, late = 6:7)
                   
norm_MME <- fitdistrplus::fitdist(rhScores,distr = "norm",method = "mme")

# estimate mean:
mean_est <- norm_MME$estimate[1]
# estimate sd:
sd_est <- norm_MME$estimate[2]


svg("S9.svg",width = 9,height = 5)

norm_distr <- function(x){ return(dnorm(x = x,mean = mean_est,sd = sd_est)) }

# plot the density function and the histogram of rhScores
curve(norm_distr,
      xlim = c(min(rhScores),max(c(rhScores,rhScore_real))),
      col = "steelblue",lwd = 5,xlab = "scores", ylab="Frequency")

# plot the histogram of rhScores
hist(rhScores,prob = TRUE,add = TRUE, breaks = 100,main = "10,000 permutations")
rug(rhScores)

# plot a red line at the position where we can find the real variance
abline(v = rhScore_real, lwd = 5, col = "red")


dev.off()
```



