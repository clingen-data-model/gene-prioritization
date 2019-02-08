# GenePrioritization
Python script to generate Gene Prioritization Reports for ClinGen curation.


## About this project
This script generates an Excel file with a row for each gene that has one of the following:
#a) at least one ClinVar 1-star variant
#b) a + or # OMIMID
#c) a MANE RefSeq transcript
#d) a VCEP
#e) a GCEP
#f) a Clinical Validity curation
#g) a Gene Dosage curation
#h) an Actionability curation

This script uses the following files:
#'Gene2GeneIDSyns.txt' #From HGNC genenames.org downloads (Gene Symbol, Synonyms, NCBI GeneID)
#'gene_info_human.txt' #From: https://ftp.ncbi.nih.gov/gene/DATA/gene_info extracted '9606'
#'test_version_20181206.txt' #From https://ftp.ncbi.nih.gov/pub/GTR/data/_README.html
#'variant_summary.txt' #From ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
#'mimTitles.txt' #From https://www.omim.org/downloads/
#'mim2gene_medgen.txt' #From https://www.omim.org/downloads/
#'MANE.GRCh38.v0.5.summary.txt' #From ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
#'VCEPGeneList.txt' #Internal file on DropBox Shared U41/ClinVar/ClinVarReports
#'ClinGen-Gene-Disease-Summary-2019-02-04.csv' #From https://search.clinicalgenome.org/kb/ - Manually removed header rows
#'ClinGen-Dosage-Sensitivity-2019-02-04.csv' #From https://search.clinicalgenome.org/kb/ - Manually removed header rows
#'ClinGen-Clinical-Actionability-2019-02-04.csv' #From https://search.clinicalgenome.org/kb/ - Manually removed header rows
#'CGD.txt' #From https://research.nhgri.nih.gov/CGD/download/
#'curations_export_at_2019-02-01_11_52_06.csv' #From Gene Tracker download (needs log in)

**Gene_Prioritization_Report.py** -
  * \Tab#GTR_REPORT: main file
  * \Tab#GTR_STATS_AllCGDs: Stats of genes per test for all disease areas
  * \Tab#GTR_STATS_SingleCGD: Stats of genes per test for ONLY A SINGLE disease area
