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

# install package 'myTAI': install.packages("myTAI")
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
* ftp://ftp.ensembl.org/pub/release-77/fasta/danio_rerio/cds/Danio_rerio.Zv9.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/takifugu_rubripes/cds/Takifugu_rubripes.FUGU4.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/astyanax_mexicanus/cds/Astyanax_mexicanus.AstMex102.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/xiphophorus_maculatus/cds/Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/oryzias_latipes/cds/Oryzias_latipes.MEDAKA1.cds.all.fa.gz
* ftp://ftp.ensembl.org/pub/release-77/fasta/gadus_morhua/cds/Gadus_morhua.gadMor1.cds.all.fa.gz


### For D. melanogaster
* ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.53_FB2013_05/fasta/dmel-all-CDS-r5.53.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r1.4_FB2012_03/fasta/dsim-all-CDS-r1.4.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_pseudoobscura/dpse_r3.1_FB2013_02/fasta/dpse-all-CDS-r3.1.fasta.gz
* ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_r1.2_FB2012_01/fasta/dvir-all-CDS-r1.2.fasta.gz


### For A. thaliana 
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Athaliana/annotation/Athaliana_167_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Alyrata/annotation/Alyrata_107_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Brapa/annotation/Brapa_197_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Crubella/annotation/Crubella_183_cds.fa.gz
* ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Thalophila/annotation/Thalophila_173_cds.fa.gz

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

```r

library(orthologr)

### For Danio rerio

# compute the divergence map of D. rerio vs. T. rubripes
Drerio_vs_Trubripes_DM <- divergence_stratigraphy(
                         query_file = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file = "Takifugu_rubripes.FUGU4.cds.all.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         



# compute the divergence map of D. rerio vs. A. mexicanus
Drerio_vs_Amexicanus_DM <- divergence_stratigraphy(
                         query_file = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file = "Astyanax_mexicanus.AstMex102.cds.all.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         
                         
                         

# compute the divergence map of D. rerio vs. X. maculatus
Drerio_vs_Xmaculatus_DM <- divergence_stratigraphy(
                         query_file = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file = "Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )



# compute the divergence map of D. rerio vs. O. latipes
Drerio_vs_Olatipes_DM <- divergence_stratigraphy(
                         query_file = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file = "Oryzias_latipes.MEDAKA1.cds.all.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )



# compute the divergence map of D. rerio vs. G. morhua
Drerio_vs_Gmorhua_DM <- divergence_stratigraphy(
                         query_file = "Danio_rerio.Zv9.cds.all.fa",
                         subject_file = "Gadus_morhua.gadMor1.cds.all.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )




### For Drosophila melanogaster

# compute the divergence map of D. melanogaster vs. D. simulans
Dmel_vs_Dsim_DM <- divergence_stratigraphy(
                         query_file = "dmel-all-CDS-r5.53.fasta",
                         subject_file = "dsim-all-CDS-r1.4.fasta",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         


# compute the divergence map of D. melanogaster vs. D. pseudoobscura
Dmel_vs_Dpse_DM <- divergence_stratigraphy(
                         query_file = "dmel-all-CDS-r5.53.fasta",
                         subject_file = "dpse-all-CDS-r3.1.fasta",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         
                         

# compute the divergence map of D. melanogaster vs. D. virilis
Dmel_vs_Dvir_DM <- divergence_stratigraphy(
                         query_file = "dmel-all-CDS-r5.53.fasta",
                         subject_file = "dvir-all-CDS-r1.2.fasta",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         


### For Arabidopsis thaliana

# compute the divergence map of A. thaliana vs. A. lyrata
Ath_vs_Aly_DM <- divergence_stratigraphy(
                         query_file = "Athaliana_167_cds.fa",
                         subject_file = "Alyrata_107_cds.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         


# compute the divergence map of A. thaliana vs. B. rapa
Ath_vs_Brapa_DM <- divergence_stratigraphy(
                         query_file = "Athaliana_167_cds.fa",
                         subject_file = "Brapa_197_cds.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
  
  
  
# compute the divergence map of A. thaliana vs. C. rubella
Ath_vs_Crubella_DM <- divergence_stratigraphy(
                         query_file = "Athaliana_167_cds.fa",
                         subject_file = "Crubella_183_cds.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )
                         

# compute the divergence map of A. thaliana vs. T. halophila
Ath_vs_Thalophila_DM <- divergence_stratigraphy(
                         query_file = "Athaliana_167_cds.fa",
                         subject_file = "Thalophila_173_cds.fa",
                         eval = "1E-5", ortho_detection = "RBH",
                         comp_cores = 1, quiet = TRUE, clean_folders = TRUE )

```

## Reading PhyloExpressionSets and DivergenceExpressionSets

After performing __Phylostratigraphy__ and __Divergence Stratigraphy__ PhyloExpressionSets and
DivergenceExpressionSets can be obtained by matching the corresponding _phylostratigraphic maps_ and _divergence maps_ of _Danio rerio_, _Drosophila melanogaster_, and _Arabidopsis thaliana_ with the corresponding transcriptome data sets covering
the embryogenesis of _Danio rerio_, _Drosophila melanogaster_, and _Arabidopsis thaliana_ (see Methods in Drost et al., 2014 for details)

```r

library(gdata)

### read data sets

## read PhyloExpressionSets
Drerio_PhyloExpressionSet <- read.xls("TAI_computation_data.xls",sheet = 1)
Dmelanogaster_PhyloExpressionSet <- read.xls("TAI_computation_data.xls",sheet = 2)
Athaliana_PhyloExpressionSet <- read.xls("TAI_computation_data.xls",sheet = 3)


## read DivergenceExpressionSets

# Danio rerio

# D. rerio vs. F. rubripes
Drerio_vs_Frubripes_DivergenceExpressionSet <- read.csv("Danio_Takifugu_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. rerio vs. A. mexicanus
Drerio_vs_Amex_DivergenceExpressionSet <- read.csv("Danio_Astyanax_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. rerio vs. X. maculatus
Drerio_vs_Xmac_DivergenceExpressionSet <- read.csv("Danio_Xiphophorus_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. rerio vs. O. latipes
Drerio_vs_Olat_DivergenceExpressionSet <- read.csv("Danio_Oryzias_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. rerio vs. G. morhua
Drerio_vs_Gmor_DivergenceExpressionSet <- read.csv("Danio_Gadus_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)



# Drosophila melanogaster

# D. melanogaster vs. D. simulans
Dmel_Dsim_DivergenceExpressionSet <- read.csv("dmel_dsim_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. melanogaster vs. D. pseudoobscura
Dmel_Dpse_DivergenceExpressionSet <- read.csv("dmel_dpse_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. melanogaster vs. D. persimilis
Dmel_Dper_DivergenceExpressionSet <- read.csv("dmel_dper_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# D. melanogaster vs. D. virilis
Dmel_Dvir_DivergenceExpressionSet <- read.csv("dmel_dvir_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)



# Arabidopsis thaliana
# A. thaliana vs A. lyrata
Ath_Aly_DivergenceExpressionSet <- read.csv("Athaliana_Alyrata_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# A thaliana vs. B. rapa 
Ath_Bra_DivergenceExpressionSet <- read.csv("Athaliana_Brapa_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# A thaliana vs. T. halophila
Ath_Tha_DivergenceExpressionSet <- read.csv("Athaliana_Thalophila_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)
# A thaliana vs. C. rubella
Ath_Crub_DivergenceExpressionSet <- read.csv("Athaliana_Crubella_RBH_decil_kaks_expFile.csv", sep = ";", header = TRUE)


```

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
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkblue", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)


par(xpd = TRUE)
legend("topleft",legend = expression(bold("A")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Dmel_Dsim_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:3, mid = 4:5, late = 6:12),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "magenta", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("B")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

PlotPattern(Ath_Aly_DivergenceExpressionSet, TestStatistic = "ReductiveHourglassTest", 
            permutations = 10000, modules = list(early = 1:2, mid = 3:5, late = 6:7),
            shaded.area = TRUE, p.value = TRUE, y.ticks = 5, type = "l", lwd = 6, col = "darkgreen", 
            ylab = "", xlab = "Ontogeny", las = 3, cex.lab = 1.5, cex.axis = 1.5)

par(xpd = TRUE)
legend("topleft",legend = expression(bold("C")),bty = "n",cex = 1.5,inset = c(-0.08,-0.15))
box()
title(ylab = "TDI", mgp = c(3,0.5,0), cex.lab = 1.5)

dev.off()

```

