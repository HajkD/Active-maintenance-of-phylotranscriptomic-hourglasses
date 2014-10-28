# Reproducible Scripts for the publication: 

Drost HG, Gabel A, Grosse I, Quint M (2014). __Active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis__ Mol. Biol. Evol. (In Review)


## Performing Divergence Stratigraphy

__Divergence Stratigraphy__ is the process of quantifying the selection pressure (in terms of amino acid sequence divergence) acting on orthologous genes between closely related species. The resulting sequence divergence map (short divergence map), stores the divergence stratum in the first column and the query_id of inferred orthologous genes in the second column ([Quint et al., 2012](http://www.nature.com/nature/journal/v490/n7418/full/nature11394.html); Drost et al., 2014).

Following steps are performed to obtain a standard divergence map based on the `divergence_stratigraphy()` function:

1) Orthology Inference using BLAST reciprocal best hit ("RBH") based on blastp

2) Pairwise global amino acid alignments of orthologous genes using the [Needleman-Wunsch algorithm](http://www.sciencedirect.com/science/article/pii/0022283670900574)

3) Codon alignments of orthologous genes using [PAL2NAL](http://www.bork.embl.de/pal2nal/)

4) dNdS estimation using [Comeron's method (1995)](http://link.springer.com/article/10.1007/BF00173196)

5) Assigning estimated dNdS values to divergence strata (deciles of all dNdS values)


```r

# install package 'orthologr' from: https://github.com/HajkD/orthologr
# library(devtools)
# install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)
library(orthologr)

# install package 'myTAI': https://github.com/HajkD/myTAI
# library(devtools)
# install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)
library(myTAI)

# install.packages("gdata")
library(gdata)

```


### Retrieving divergence maps for D. rerio, D. melanogaster, and A. thaliana

To perform __Divergence Stratigraphy__ using `orthologr` you need the following prerequisites

* a CDS file covering all protein coding genes of the query organism of interest
* a CDS file covering all protein coding genes of the subject organism of interest

First, the CDS sequences for all protein coding genes of all query and subject organisms need to
be obtained.

The CDS retrieval can be done using a `Terminal` or by manual downloading the files

### For D. rerio
* Danio_rerio.Zv9.75.cds.all.fa
* Takifugu_rubripes.FUGU4.cds.all.fa
* Astyanax_mexicanus.AstMex102.cds.all.fa
* Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa
* Oryzias_latipes.MEDAKA1.cds.all.fa
* Gadus_morhua.gadMor1.cds.all.fa

### For D. melanogaster
* dmel-all-CDS-r5.53.fasta
* dsim-all-CDS-r1.4.fasta
* dpse-all-CDS-r3.1.fasta
* dper-all-CDS-r1.3.fasta
* dvir-all-CDS-r1.2.fasta

### For A. thaliana 
* Athaliana_167_cds.fa
* Alyrata_107_cds.fa
* Brapa_197_cds.fa
* Crubella_183_cds.fa
* Thalophila_173_cds.fa


```shell

# download CDS file of A. thaliana
curl ftp://ftp.ensemblgenomes.org/pub/
plants/release-23/fasta/arabidopsis_thaliana/
cds/Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz 
-o Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz

# unzip the fasta file
gunzip -d Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz

# download CDS file of A. lyrata

curl ftp://ftp.ensemblgenomes.org/pub/plants/
release-23/fasta/arabidopsis_lyrata/cds/
Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz 
-o Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz

# unzip the fasta file
gunzip -d Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz

```
When the download is finished you need to unzip the files and
then start R to perform the following analyses. Make sure you
have all fasta files stored in the current working directory of
your R session.

```r
library(orthologr)

# compute the divergence map of A. thaliana vs. A. lyrata
Ath_vs_Aly_DM <- divergence_stratigraphy(
                         query_file = "Arabidopsis_thaliana.TAIR10.23.cds.all.fa",
                         subject_file = "Arabidopsis_lyrata.v.1.0.23.cds.all.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )

```

## Reading PhyloExpressionSets and DivergenceExpressionSets

After performing __Phylostratigraphy__ and __Divergence Stratigraphy__ PhyloExpressionSets and
DivergenceExpressionSets can be obtained by matching the corresponding _phylostratigraphic maps_ and _divergence maps_ of _Danio rerio_, _Drosophila melanogaster_, and _Arabidopsis thaliana_ with the corresponding transcriptome data sets covering
the embryogenesis of _Danio rerio_, _Drosophila melanogaster_, and _Arabidopsis thaliana_ (see Methods in Drost et al., 2014 for details)

```r

### read data sets

## read PhyloExpressionSets
Drerio_PhyloExpressionSet <- read.xls()
Dmelanogaster_PhyloExpressionSet <- read.xls()
Athaliana_PhyloExpressionSet <- read.xls()


## read DivergenceExpressionSets

# Danio rerio

# D. rerio vs. F. rubripes
Drerio_vs_Frubripes_DivergenceExpressionSet <- read.xls()
# D. rerio vs. A. mexicanus
Drerio_vs_Amex_DivergenceExpressionSet <- read.xls()
# D. rerio vs. X. maculatus
Drerio_vs_Xmac_DivergenceExpressionSet <- read.xls()
# D. rerio vs. O. latipes
Drerio_vs_Olat_DivergenceExpressionSet <- read.xls()
# D. rerio vs. G. morhua
Drerio_vs_Gmor_DivergenceExpressionSet <- read.xls()



# Drosophila melanogaster

# D. melanogaster vs. D. simulans
Dmel_Dsim_DivergenceExpressionSet <- read.xls()
# D. melanogaster vs. D. ananassae
Dmel_Dana_DivergenceExpressionSet <- read.xls()
# D. melanogaster vs. D. pseudoobscura
Dmel_Dpse_DivergenceExpressionSet <- read.xls()
# D. melanogaster vs. D. persimilis
Dmel_Dper_DivergenceExpressionSet <- read.xls()
# D. melanogaster vs. D. virilis
Dmel_Dvir_DivergenceExpressionSet <- read.xls()



# Arabidopsis thaliana
# A. thaliana vs A. lyrata
Ath_Aly_DivergenceExpressionSet <- read.xls()
# A thaliana vs. B. rapa 
Ath_Bra_DivergenceExpressionSet <- read.xls()
# A thaliana vs. T. halophila
Ath_Tha_DivergenceExpressionSet <- read.xls()
# A thaliana vs. C. rubella
Ath_Crub_DivergenceExpressionSet <- read.xls()


```

## Testing the statistical significance of observed phylotranscriptomics patterns




## Generating Figures

First load the following packages into the workspace:

```r

library(myTAI)

library(gdata)

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
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TAI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Dmelanogaster_PhyloExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TAI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Athaliana_PhyloExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
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

PlotCorrelation(Drerio_PhyloExpressionSet, Drerio_vs_Frubripes_DivergenceExpressionSet , 
                method = "kendall", main.text = "Kendall's ")

PlotCorrelation(Dmelanogaster_PhyloExpressionSet, Dmel_Dsim_DivergenceExpressionSet, 
                method = "kendall", main.text = "Kendall's ")

PlotCorrelation(Athaliana_PhyloExpressionSet, Ath_Aly_DivergenceExpressionSet, 
                method = "kendall", main.text = "Kendall's ")

dev.off()
```


### Figure 4

```r
svg("Fig4.svg",width = 16.9,height = 5)
par(mfrow = c(1,3))
par(mar = c(1.5, 0.5, 0.5, 0.1))
par(mai = c(1.4,0.6,0.5,0.1))
par(mgp = c(8,1,0))

PlotPattern(Drerio_vs_Frubripes_DivergenceExpressionSet , TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:18, mid = 19:36, late = 37:40),
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Dmel_Dsim_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Ath_Aly_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = FALSE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

dev.off()

```

