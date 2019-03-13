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
  * 'Gene2GeneIDSyns.txt' - From HGNC genenames.org downloads (Gene Symbol, Synonyms, NCBI GeneID)
  * 'gene_info_human.txt' - From: https://ftp.ncbi.nih.gov/gene/DATA/gene_info extracted '9606'
  * 'test_version_20181206.txt' - From https://ftp.ncbi.nih.gov/pub/GTR/data/_README.html
  * 'variant_summary.txt' - From ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
  * 'mimTitles.txt' - From https://www.omim.org/downloads/
  * 'mim2gene_medgen.txt' - From https://www.omim.org/downloads/
  * 'MANE.GRCh38.v0.5.summary.txt' - From ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
  * 'VCEPGeneList.txt' - Internal file on DropBox Shared U41/ClinVar/ClinVarReports
  * 'ClinGen-Gene-Disease-Summary-2019-02-04.csv' - From https://search.clinicalgenome.org/kb/ - Manually removed header rows
  * 'ClinGen-Dosage-Sensitivity-2019-02-04.csv' - From https://search.clinicalgenome.org/kb/ - Manually removed header rows
  * 'ClinGen-Clinical-Actionability-2019-02-04.csv' - From https://search.clinicalgenome.org/kb/ - Manually removed header rows
  * 'CGD.txt' - From https://research.nhgri.nih.gov/CGD/download/
  * 'curations_export_at_2019-02-01_11_52_06.csv' - From Gene Tracker download (needs log in)
  * 'Concert_Genetics_1000_Most_Popular_Genes.txt' - Obtained internally

**Gene_Prioritization_Report.py**
  * Tab#GTR_REPORT: main file
  * Tab#GTR_STATS_AllCGDs: Stats of genes per test for all disease areas
  * Tab#GTR_STATS_SingleCGD: Stats of genes per test for ONLY A SINGLE disease area

## How to run these scripts
All scripts are run as 'python3 *filename.py*'.
