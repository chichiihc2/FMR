# Data

The real-data analysis (Section 5) uses data of the Alzheimer's Disease Neuroimaging
Initiative (ADNI), which cannot be redistributed. Apply for access at
<https://adni.loni.usc.edu> and download the files below from the LONI Image and Data Archive
(IDA). Nothing under `data/raw/` or `data/processed/` is ever committed (see `.gitignore`).

## 1. What to put in `data/raw/`

```
data/raw/
├── filtered_ADNI_808_plink2.bed       ADNI whole-genome sequencing, 808 subjects, after QC
├── filtered_ADNI_808_plink2.bim       (PLINK 1 binary, GRCh37)
├── filtered_ADNI_808_plink2.fam       IID = ADNI PTID, e.g. 003_S_1057
├── ADNI_Gene_Expression_Profile.csv   ADNI blood gene expression (Affymetrix U219), as downloaded
├── data_cleaned_demo.csv              phenotypes of Section 5.1 (see below)
└── ADNIMERGE2/                        the ADNIMERGE2 R data package (see below)
```

**Genotypes.** The PLINK files are the ADNI WGS VCF
(`ADNI.808_indiv.minGQ_21.pass_combined.vcf.gz`, "ADNI WGS" in IDA's genetic data) after
the QC of the supplement (Section S.4), done by `realdata/prep/wgs_qc/01_filter_vcf.sh` and
`02_qc_plink.sh`. If you put the VCF in `data/raw/` instead of the PLINK files,
`run_realdata.R` runs the QC first (needs `vcftools` and `plink2`).

**Gene expression.** `ADNI_Gene_Expression_Profile.csv` exactly as downloaded (probe level;
8 header rows of sample information, then one row per probe set).

**Phenotypes of Section 5.1** (`data_cleaned_demo.csv`): one row per visit; the first row of
each `PTID` is used. Columns used: `PTID`, `Age`, `Gender`, and the outcomes `ADAS13`, `MMSE`,
`PTAU`, `FDG`, `ABETA`. It is made by `realdata/prep/R/01_clean_phenotypes.R` from three IDA
downloads; if you put these in `data/raw/` instead, `run_realdata.R` makes the table first:
`ADNI_adnimerge_20170629_QT-freeze.csv` (ADNIMERGE, 2017-06-29 freeze),
`ADNINIGHTINGALE_20210219_04Jan2024.csv`, and `Demo.csv`.

**Phenotypes of Section 5.2** (`ADNIMERGE2/`): the ADNIMERGE2 R data package, version 0.1.2
(IDA download `adnimerge_download.zip`, which contains `ADNIMERGE2.tar.gz`). Put the zip, the
tar.gz or the unpacked folder in `data/raw/`; the tables are read with `load()`, and the
package does not need to be installed. Tables used (`data/raw/ADNIMERGE2/data/`): ADAS,
ADSP_PHC_PET_Amyloid_Simple, FAQ, PTDEMOG, VITALS.

`PTID.csv` (subject IDs in `.fam` order) is made from the `.fam` file automatically.

## 2. Checking the data release

`data/checksums.md5` holds the MD5 sums of the files used for the paper. `run_realdata.R`
compares them and warns on a mismatch: another ADNI release (in particular of ADNIMERGE2)
changes the sample and the estimates. By hand: `md5sum -c data/checksums.md5` (Linux) or
`md5 -r` (macOS).

## 3. Public data

* `data/external/gene_locations37.csv`: gene coordinates (GRCh37), used to assign variants to
  genes. Included.
* `data/external/eqtlgen/`: eQTLGen cis-eQTL summary statistics (Vosa et al., 2021; about
  4 GB), downloaded by `realdata/prep/eqtlgen/00_download_eqtlgen.sh` from
  <https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/>.

## 4. What the pipeline writes to `data/processed/` (about 60 GB)

| File | Script |
|---|---|
| `filtered_ADNI_808_plink2/pca50_ADNI_808_plink2.{eigenvec,csv}`: genotype PCs | `realdata/prep/wgs_qc/03_pca_plink.sh`, `realdata/prep/R/04_convert_pca_eigenvec.R` |
| `filtered_ADNI_808_plink2/{bim.csv, gene_snps.csv}`: variant to gene map | `realdata/prep/R/03_assign_snps_to_genes.R` |
| `Genes_pca6/<gene>.csv`: instrument of each gene, adjusted for PC1-PC6 (56 GB) | `realdata/prep/R/05_adjust_FIV_with_PCA.R` |
| `ADNI_Gene_Expression.csv`, `Gene_expression_per_gene/<gene>.csv` | `realdata/prep/R/02_clean_gene_expression.R`, `06_split_expression_per_gene.R` |
| `eqtlgen/eqtlgen_adni_cis.csv.gz` | `realdata/prep/eqtlgen/01_build_eqtlgen_gene_tables.R` |
| `apoe/`: APOE-region genotypes, PCs, relatedness, phenotypes | `realdata/apoe/genotypes.sh`, `realdata/apoe/phenotypes.R` |

Scenario S1 of the simulation study uses `Genes_pca6/DLGAP2.csv`,
`Gene_expression_per_gene/DLGAP2.csv`, `data_cleaned_demo.csv` and `PTID.csv`; its Table S.1
cells also use `filtered_ADNI_808_plink2/pca50_ADNI_808_plink2.csv`.
