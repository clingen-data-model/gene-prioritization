# GenePrioritization
Python script to generate Gene Prioritization Reports for ClinGen curation.

## About this project
This script generates an Excel file with a row for each gene that has one of the following:
* a) at least one test in the GTR
* b) at least one ClinVar variant
* c) a + or # OMIMID
* d) a MANE RefSeq transcript
* e) a VCEP
* f) a GCEP
* g) a Clinical Validity curation
* h) a Gene Dosage curation
* i) an Actionability curation

This script uses the following files:
  * 'Gene2GeneIDSyns.txt' - From https://www.genenames.org/download/custom/ (Gene Symbol, Previous Gene Symbol, Synonyms, NCBI GeneID, Approved)
  * 'Homo_sapiens.gene_info.gz' - From https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/
  * 'test_version.gz' - From https://ftp.ncbi.nih.gov/pub/GTR/data/
  * 'variant_summary.txt' - From ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
  * 'mimTitles.txt' - From https://www.omim.org/downloads/
  * 'mim2gene_medgen.txt' - From ftp://ftp.ncbi.nih.gov/gene/DATA/
  * 'MANE.GRCh38.v0.5.summary.txt' - From ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
  * 'VCEPGeneList.txt' - Internal file on DropBox Shared U41/ClinVar/ClinVarReports
  * 'ClinGen-Gene-Disease-Summary-2019-03-14.csv' - From https://search.clinicalgenome.org/kb/gene-validity.csv - Manually removed header rows
  * 'ClinGen-Dosage-Sensitivity-2019-03-14.csv' - From https://search.clinicalgenome.org/kb/gene-dosage.csv -  Manually removed header rows
  * 'ClinGen-Clinical-Actionability-2019-03-14.csv' - From https://search.clinicalgenome.org/kb/actionability.csv - Manually removed header rows
  * 'CGD.txt.gz' - From https://research.nhgri.nih.gov/CGD/download/txt/
  * 'curations_export_at_2019-03-14_09_25_28.csv' - From https://clingen.sirs.unc.edu/#/curations/export (needs log in)
  * 'Concert_Genetics_1000_Most_Popular_Genes.txt' - Obtained internally

**Gene_Prioritization_Report.py**
  * Tab#GTR_REPORT: main file
  * Tab#GTR_STATS_AllCGDs: Stats of genes per test for all disease areas
  * Tab#GTR_STATS_SingleCGD: Stats of genes per test for ONLY A SINGLE disease area

## How to run these scripts
All scripts are run as 'python3 *filename.py*'.
