#Must run as >python3 Gene_Prioritization_Report.py
#This script outputs an Excel file with a row for each gene that has one of the following:
#a) at least one ClinVar 1-star variant
#b) a + or # OMIMID
#c) a MANE RefSeq transcript
#d) a VCEP
#e) a GCEP
#f) a Clinical Validity curation
#g) a Gene Dosage curation
#h) an Actionability curation
#This script also outputs two stats tabs based on number of tests per gene and disease area(s)

import re
import pprint
import csv
import xlsxwriter

gtrHash = {}
omimHash = {}
osystems = []
diseases = {}
diseasesOne = {}

def create_gene2idlookup1(infile):
    '''This function makes a lookup of GeneSym to GeneID from HGNC'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()

            if not line.startswith('Approved'): #ignore first line
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    geneSym = col[0].upper()
                    synonyms = col[1]
                    geneID = col[2].rstrip('\n\r')
                    if geneID and geneID != '':
                        geneID = int(geneID)
                        gtrHash[geneSym] = {'GeneID':geneID}
                        if synonyms != '':
                            gtrHash[geneSym].update({'Synonyms':col[1].split(', ')})

    input.close()
    return(gtrHash)


def create_gene2idlookup2(infile):
    '''This function adds to the lookup, GeneSym to GeneID for additional NCBI Genes'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()

            if not line.startswith('#'): #ignore first line
                col = re.split(r'\t', line) #split on tabs
                if col[0] != '':
                    geneID = int(col[1])
                    geneSym = col[2].upper()
                    if geneID and geneID != '' and geneSym not in gtrHash.keys():
                        gtrHash[geneSym] = {'GeneID':geneID}

    input.close()
    return(gtrHash)


def gene_convert(gene):
    '''This function converts geneSyms in the GTR that are not in the gtrHash'''

    if gene == 'GS1-259H13.2':
        return('TMEM225B')
    elif gene == 'NARR':
        return('RAB34')
    elif gene == 'CH17-360D5.1':
        return('NPY4R2')
    elif gene == 'CRHR1-IT1-CRHR1':
        return('LINC02210-CRHR1')
    elif gene == 'HMP19':
        return('NSG2')
    elif gene == 'WTH3DI':
        return('RAB6D')
    elif gene == 'LOC200726':
        return('FAM237A')
    elif gene == 'LOC100132146':
        return('FAM240A')
    elif gene == 'LOC100996693':
        return('SPEGNB')
    elif gene == 'LOC100130705':
        return('ATP6V1FNB')
    elif gene == 'LOC81691':
        return('REXO5')
    elif gene == 'LOC100505478':
        return('TEX48')
    elif gene == 'LOC100129216':
        return('DEFB131B')
    elif gene == 'LOC653602':
        return('FAM84A')
    elif gene == 'CCNA2':
        return('CCNA2')
    else:
        return(gene)


def gene_check(gene):
    '''This function checks if the geneSym in the GTR is currently approved or assigns based on unique synonym'''

    if gene.upper() in gtrHash.keys():
        return(gene)
    else:
        count = 0
        for record in gtrHash:
            if 'Synonyms' in gtrHash[record].keys():
                for syn in gtrHash[record]['Synonyms']:
                    if gene == syn.upper():
                        count +=  1
                        gene = record
        if count == 1:
            return(gene)
        else:
            return('Not found: ' + gene)


def create_gtrHash(infile):
    '''This function makes a hash of metadata for each GeneSym in the GTR registry'''

    unique_tests = []

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()

            if not line.startswith('test'): #ignore first line
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    gtr = col[0]
                    lab_name = col[1]
                    lab_name = lab_name.strip('\'"')
                    test_name = col[11]
                    test_name = test_name.strip('\'"')
                    current_lab_test = lab_name + ' | ' + test_name
                    conditions = col[15].split('|') #pipe separated
                    indications = col[16].split('|') #pipe separated
                    geneSymList = col[21].split('|') #pipe separated
                    for geneSym in geneSymList:
                        if geneSym == 'BCR-ABL':
                            geneSymList.remove(geneSym)
                            geneSymList.append('BCR')
                            geneSymList.append('ABL1')
                    current_version = int(col[23])
                    current_test = int(col[24])
                    current_lab = int(col[26])

                    if current_version == 1 and current_test == 1 and current_lab == 1 and len(geneSymList) != 0 \
                       and current_lab_test not in unique_tests: #ensure test and lab are current and genelist not empty
                        unique_tests.append(current_lab_test)
                        for geneSym in geneSymList:
                            if geneSym != '':
                                geneSym = gene_convert(geneSym)
                                geneSym = gene_check(geneSym.upper())

                                if 'Not found' not in geneSym:
                                    if 'geneCount' not in gtrHash[geneSym].keys():
                                        gtrHash[geneSym].update({'geneCount': 0})

                                    gtrHash[geneSym]['geneCount'] += 1

                                    if 'Conditions' not in gtrHash[geneSym].keys():
                                        gtrHash[geneSym].update({'Conditions': []})

                                    for condition in conditions:
                                        condition = condition.strip('\'"')
                                        if condition not in gtrHash[geneSym]['Conditions']:
                                            gtrHash[geneSym]['Conditions'].append(condition)

                                    if 'Indications' not in gtrHash[geneSym].keys():
                                        gtrHash[geneSym].update({'Indications': []})

                                    for indication in indications:
                                        indication = indication.strip('\'"')
                                        if indication not in gtrHash[geneSym]['Indications']:
                                            gtrHash[geneSym]['Indications'].append(indication)

                                    if 'LabNames' not in gtrHash[geneSym].keys():
                                        gtrHash[geneSym].update({'LabNames': []})

                                    gtrHash[geneSym]['LabNames'].append(lab_name)

                                    if 'TestNames' not in gtrHash[geneSym].keys():
                                        gtrHash[geneSym].update({'TestNames': []})

                                    gtrHash[geneSym]['TestNames'].append(test_name)
                                else:
                                    print(geneSym)

    input.close()
    return(gtrHash)


def add_VarCounts(infile):
    '''This function adds ClnVar varID counts to the gtr Hash'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()

            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '' and 'GRCh38' in col[16]: #ignore empty lines
                    geneID = int(col[3])
                    geneSym = col[4].upper()
                    clinSig = col[6].split(', ')
                    revStat = col[24]
                    varID = int(col[30].rstrip('\n\r'))

                    if geneSym in gtrHash.keys():
                        if 'VarCount' not in gtrHash[geneSym].keys():
                            gtrHash[geneSym].update({'VarCount': 0})
                            gtrHash[geneSym].update({'VarCountP': 0})
                            gtrHash[geneSym].update({'VarCountC': 0})

                        gtrHash[geneSym].update({'GeneID': geneID})

                        if 'no assertion criteria provided' not in revStat and \
                           'no interpretation for the single variant' not in revStat and \
                           'no assertion provided' not in revStat:
                           gtrHash[geneSym]['VarCount'] += 1

                        if 'no assertion criteria provided' not in revStat and \
                           'no interpretation for the single variant' not in revStat and \
                           'no assertion provided' not in revStat and ('Pathogenic' in clinSig or \
                           'Likely pathogenic' in clinSig or 'Pathogenic/Likely pathogenic' in clinSig or \
                           'drug response' in clinSig or 'risk factor' in clinSig):
                            gtrHash[geneSym]['VarCountP'] += 1

                        if 'no assertion criteria provided' not in revStat and \
                           'no interpretation for the single variant' not in revStat and \
                           'no assertion provided' not in revStat and 'Conflicting interpretations of pathogenicity' in clinSig:
                            gtrHash[geneSym]['VarCountC'] += 1

    input.close()
    return(gtrHash)


def create_omimHash(infile):
    '''This function makes a list # and + OMIM IDs'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()
            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if col[0] == 'Number Sign' or col[0] == 'Plus':
                    omim = col[1]
                    disease = col[2].split(';', 1)[0]
                    omimHash[omim] = {'Disease':disease}

    input.close()
    return(omimHash)


def add_OMIMID(infile):
    '''This function adds the OMIM ID to the gtr Hash'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()
            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '' and col[0] in omimHash and '-' not in col[1]: #ignore empty lines

                    MIM = col[0]
                    geneID = int(col[1])

                    for geneSym in gtrHash:
                        if 'GeneID' in gtrHash[geneSym].keys() and gtrHash[geneSym]['GeneID'] == geneID:
                            if 'MIM#' not in gtrHash[geneSym].keys():
                                gtrHash[geneSym].update({'MIM#': []})

                            gtrHash[geneSym]['MIM#'].append(MIM)

                            if 'Diseases' not in gtrHash[geneSym].keys():
                                gtrHash[geneSym].update({'Diseases': []})

                            gtrHash[geneSym]['Diseases'].append(omimHash[MIM]['Disease'])

    input.close()
    return(gtrHash)


def add_MANEflag(infile):
    '''This function adds the MANE flag to the gtr Hash'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()
            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines

                    geneID = int(col[0].split(':', 1)[1])
                    geneSym = col[3].rstrip('\n\r')
                    geneSym = gene_check(geneSym.upper())
                    RefSeq = col[5]

                    if 'GeneID' in gtrHash[geneSym].keys() and gtrHash[geneSym]['GeneID'] == geneID:
                        gtrHash[geneSym].update({'MANEstatus': RefSeq})

    input.close()
    return(gtrHash)


def add_VCEPs(infile):
    '''This function adds the VCEP name to the GTR hash'''

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()
            if not line.startswith('Gene'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    geneSym = col[0].upper()
                    VCEP = col[1].rstrip('\n\r')

                    gtrHash[geneSym].update({'VCEP': VCEP})

    input.close()
    return(gtrHash)


def add_GCEPs(infile):
    '''This function adds the GCEPs and status to the GTR hash'''

    with open(infile, newline='') as input:
        reader = csv.DictReader(input,dialect='excel')
        rows = []

        for row in reader:
            rows.append(row)
            geneSym = row['Gene Symbol'].rstrip()
            geneSym = geneSym.split(' ', 1)[0]
            if geneSym == 'NKX2.2':
                geneSym = 'NKX2-2'
            if geneSym == 'DHPR':
                geneSym = 'QDPR'
            if geneSym == 'PLKHM2':
                geneSym = 'PLEKHM2'
            if geneSym == '3-PGDH':
                geneSym = 'PHGDH'
            if geneSym == 'SYN4':
                geneSym = 'SNTG1'
            if geneSym != 'GCS' and '(' not in geneSym and geneSym != 'AGL1' and \
               geneSym != 'PCC' and geneSym != 'MAT' and geneSym != '3PSPH' and \
               geneSym != 'ALD18A1':
                geneSym = gene_check(geneSym.upper())

                if row['Expert Panel'] != '':
                    if 'GCEP' not in gtrHash[geneSym].keys():
                        gtrHash[geneSym].update({'GCEP': []})

                    gtrHash[geneSym]['GCEP'].append(row['Expert Panel'])

                    if row['Curation Approved date'] != '':
                        status = 'Approved'
                    elif row['Curation Provisional date'] != '':
                        status = 'Provisional'
                    elif row['Curation In Progress date'] != '':
                        status = 'In progress'
                    else:
                        status = 'Not started'

                    if 'GCEP_status' not in gtrHash[geneSym].keys():
                        gtrHash[geneSym].update({'GCEP_status': []})

                    gtrHash[geneSym]['GCEP_status'].append(status)

    input.close()
    return(gtrHash)


def add_GeneValidity(infile):
    '''This function adds the GeneValidity status to the GTR hash'''

    with open(infile, newline='') as input:
        reader = csv.DictReader(input,dialect='excel')
        rows = []

        for row in reader:
            rows.append(row)
            geneSym = row['GENE SYMBOL']
            geneSym = gene_check(geneSym.upper())
            string = row['CLASSIFICATION'] + ' (' + row['DISEASE LABEL'] + ')'
            if 'Validity' not in gtrHash[geneSym].keys():
                gtrHash[geneSym].update({'Validity': []})

            gtrHash[geneSym]['Validity'].append(string)

    input.close()
    return(gtrHash)


def MyFn(clinSig):
    '''This function orders the Clinical Validity classification based on ranking'''

    if 'Definitive' in clinSig:
        return 1
    if 'Strong' in clinSig:
        return 2
    if 'Moderate' in clinSig:
        return 3
    if 'Limited' in clinSig:
        return 4
    if 'No Reported Evidence' in clinSig:
        return 5
    if 'Disputed' in clinSig:
        return 6
    if 'Refuted' in clinSig:
        return 7
    else:
        return 8


def add_GeneDosage(infile):
    '''This function adds the GeneDosage status to the GTR hash'''

    with open(infile, newline='') as input:
        reader = csv.DictReader(input,dialect='excel')
        rows = []

        for row in reader:
            rows.append(row)
            geneSym = row['GENE SYMBOL']
            geneSym = gene_check(geneSym.upper())

            if row['HAPLOINSUFFICIENCY'] != '':
                haplo = row['HAPLOINSUFFICIENCY']
            else:
                haplo = 'N/A'

            if row['TRIPLOSENSITIVITY'] != '':
                triplo = row['TRIPLOSENSITIVITY']
            else:
                triplo = 'N/A'

            string = 'Haplo: ' + haplo + ' | Triplo: ' + triplo

            gtrHash[geneSym].update({'Dosage': string})

    input.close()
    return(gtrHash)


def add_Actionability(infile):
    '''This function adds the Actionability status to the GTR hash'''

    with open(infile, newline='') as input:
        reader = csv.DictReader(input,dialect='excel')
        rows = []

        for row in reader:
            rows.append(row)
            geneSymList = re.split(r',', row['GENE SYMBOL']) #split on commas
            for geneSym in geneSymList:
                geneSym = geneSym.strip('\n\r')
                if geneSym != '':
                    geneSym = gene_check(geneSym.upper())
                    gtrHash[geneSym].update({'Actionability': row['CONDITION']})

    input.close()
    return(gtrHash)


def add_CGD(infile):
    '''This function adds the manifestation categories of CGD status to the GTR hash'''

    global osystems

    with open(infile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()
            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    geneSym = col[0].rstrip('\n\r')
                    geneSym = gene_check(geneSym.upper())
                    geneID = col[2]
                    manifestations = col[7]
                    system = [x.strip() for x in col[7].split(';') if x != '']
                    osystems.extend(system)
                    gtrHash[geneSym].update({'Manifestations': manifestations})

    osystems = sorted(set(osystems))
    input.close()
    return(gtrHash)


def get_stats():
    '''This function adds the Stat counts to the diseases and diseasesOne hashes'''

    tests = ['TestAll', 'Test10', 'Test20', 'Test50', 'Test75', 'Test100']
    i = 1
    for test in tests:
        if test not in diseases.keys():
            diseases[test] = {}
            diseasesOne[test] = {}
        for organ in osystems:
            diseases[test][organ] = {'Count':0}
            diseasesOne[test][organ] = {'Count':0}
            i += 1

    for geneSym in gtrHash:
        if 'Manifestations' in gtrHash[geneSym].keys():
            manifestations = [x.strip() for x in gtrHash[geneSym]['Manifestations'].split(';') if x != '']
            for organ in osystems:
                for mani in manifestations:
                    if mani == organ:
                        if 'geneCount' in gtrHash[geneSym].keys():
                            diseases['TestAll'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 10:
                            diseases['Test10'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 20:
                            diseases['Test20'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 50:
                            diseases['Test50'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 75:
                            diseases['Test75'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 100:
                            diseases['Test100'][organ]['Count'] += 1

            if len(manifestations) == 1:
                for organ in osystems:
                    if manifestations[0] == organ:
                        if 'geneCount' in gtrHash[geneSym].keys():
                            diseasesOne['TestAll'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 10:
                            diseasesOne['Test10'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 20:
                            diseasesOne['Test20'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 50:
                            diseasesOne['Test50'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 75:
                            diseasesOne['Test75'][organ]['Count'] += 1
                        if 'geneCount' in gtrHash[geneSym].keys() and gtrHash[geneSym]['geneCount'] >= 100:
                            diseasesOne['Test100'][organ]['Count'] += 1

    return(diseases, diseasesOne)


def print_allgenes(outfile):
    '''This function creates the Excel output file'''

    workbook = xlsxwriter.Workbook(outfile)
    worksheet = workbook.add_worksheet('GTR_REPORT')
    worksheet1 = workbook.add_worksheet('GTR_STATS_AllCGDs')
    worksheet2 = workbook.add_worksheet('GTR_STATS_SingleCGD')

    worksheet.write(0, 0, 'GeneSym')
    worksheet.write(0, 1, 'NCBIGeneID')
    worksheet.write(0, 2, '#Tests with Gene')
    worksheet.write(0, 3, 'Indication(s)')
    worksheet.write(0, 4, 'Condition(s)')
    worksheet.write(0, 5, 'LabTest(s)')
    worksheet.write(0, 6, 'Lab(s)')
    worksheet.write(0, 7, '#ClinVar Variants (1-star or above)')
    worksheet.write(0, 8, '#P/LP, Drug response or Risk Factor ClinVar Variants (1-star or above)')
    worksheet.write(0, 9, '#ClinVar Variants (1-star or above with Conflicting interpretations of pathogenicity)')
    worksheet.write(0, 10, 'MANE SELECT RefSeq')
    worksheet.write(0, 11, '# or + OMIM IDs')
    worksheet.write(0, 12, '# or + OMIM Diseases')
    worksheet.write(0, 13, 'VCEP')
    worksheet.write(0, 14, 'GCEP(s)')
    worksheet.write(0, 15, 'Highest GCEP curation status')
    worksheet.write(0, 16, 'Gene Validity Curation(s)')
    worksheet.write(0, 17, 'Gene Dosage Curation')
    worksheet.write(0, 18, 'Actionability Curation')
    worksheet.write(0, 19, 'CGD manifestations')
    i = 20
    for organ in osystems:
        worksheet.write(0, i, organ + ' (CGD)')
        i += 1

    worksheet.write(0, i, 'Total disease areas (CGD)')

    row = 1

    for geneSym in gtrHash:
        if ('geneCount' in gtrHash[geneSym] and gtrHash[geneSym]['geneCount'] != 0) or 'MIM#' in gtrHash[geneSym] or 'MANEstatus' in gtrHash[geneSym] or \
           ('VarCount' in gtrHash[geneSym] and gtrHash[geneSym]['VarCount'] != 0) or 'Validity' in gtrHash[geneSym] or 'Dosage' in gtrHash[geneSym] or \
           'Actionability' in gtrHash[geneSym] or 'GCEP' in gtrHash[geneSym]:
            worksheet.write(row, 0, geneSym)
            worksheet.write(row, 1, gtrHash[geneSym]['GeneID'])
            if 'geneCount' in gtrHash[geneSym].keys():
                worksheet.write(row, 2, gtrHash[geneSym]['geneCount'])
            else:
                worksheet.write(row, 2, 0)

            if 'Indications' in gtrHash[geneSym].keys():
                worksheet.write(row, 3, '| '.join(sorted(set(gtrHash[geneSym]['Indications']))))
            else:
                worksheet.write(row, 3, '-')

            if 'Conditions' in gtrHash[geneSym].keys():
                worksheet.write(row, 4, ' | '.join(sorted(set(gtrHash[geneSym]['Conditions']))))
            else:
                worksheet.write(row, 4, '-')

            if 'TestNames' in gtrHash[geneSym].keys():
                worksheet.write(row, 5, ' | '.join(sorted(set(gtrHash[geneSym]['TestNames']))))
            else:
                worksheet.write(row, 5, '-')

            if 'LabNames' in gtrHash[geneSym].keys():
                worksheet.write(row, 6, ' | '.join(sorted(set(gtrHash[geneSym]['LabNames']))))
            else:
                worksheet.write(row, 6, '-')

            if 'VarCount' in gtrHash[geneSym].keys():
                worksheet.write(row, 7, gtrHash[geneSym]['VarCount'])
            else:
                worksheet.write(row, 7, 0)

            if 'VarCountP' in gtrHash[geneSym].keys():
                worksheet.write(row, 8, gtrHash[geneSym]['VarCountP'])
            else:
                worksheet.write(row, 8, 0)

            if 'VarCountC' in gtrHash[geneSym].keys():
                worksheet.write(row, 9, gtrHash[geneSym]['VarCountC'])
            else:
                worksheet.write(row, 9, 0)

            if 'MANEstatus' in gtrHash[geneSym].keys():
                worksheet.write(row, 10, gtrHash[geneSym]['MANEstatus'])
            else:
                worksheet.write(row, 10, '-')

            if 'MIM#' in gtrHash[geneSym].keys():
                worksheet.write(row, 11, ' | '.join(sorted(set(gtrHash[geneSym]['MIM#']))))
            else:
                worksheet.write(row, 11, '-')

            if 'Diseases' in gtrHash[geneSym].keys():
                worksheet.write(row, 12, ' | '.join(sorted(set(gtrHash[geneSym]['Diseases']))))
            else:
                worksheet.write(row, 12, '-')

            if 'VCEP' in gtrHash[geneSym].keys():
                worksheet.write(row, 13, gtrHash[geneSym]['VCEP'])
            else:
                worksheet.write(row, 13, '-')

            if 'GCEP' in gtrHash[geneSym].keys():
                worksheet.write(row, 14, ' | '.join(sorted(set(gtrHash[geneSym]['GCEP']))))
            else:
                worksheet.write(row, 14, '-')

            if 'GCEP_status' in gtrHash[geneSym].keys():
                if 'Approved' in gtrHash[geneSym]['GCEP_status']:
                    worksheet.write(row, 15, 'Approved')
                elif 'Provisional' in gtrHash[geneSym]['GCEP_status']:
                    worksheet.write(row, 15, 'Provisional')
                elif 'In progress' in gtrHash[geneSym]['GCEP_status']:
                    worksheet.write(row, 15, 'In progress')
                else:
                    worksheet.write(row, 15, 'Not started')
            else:
                worksheet.write(row, 15, '-')

            if 'Validity' in gtrHash[geneSym].keys():
                worksheet.write(row, 16, ' | '.join(sorted(gtrHash[geneSym]['Validity'], key=MyFn)))
            else:
                worksheet.write(row, 16, '-')

            if 'Dosage' in gtrHash[geneSym].keys():
                worksheet.write(row, 17, gtrHash[geneSym]['Dosage'])
            else:
                worksheet.write(row, 17, '-')

            if 'Actionability' in gtrHash[geneSym].keys():
                worksheet.write(row, 18, gtrHash[geneSym]['Actionability'])
            else:
                worksheet.write(row, 18, '-')

            if 'Manifestations' in gtrHash[geneSym].keys():
                worksheet.write(row, 19, gtrHash[geneSym]['Manifestations'])
            else:
                worksheet.write(row, 19, '-')

            i = 20
            count = 0
            if 'Manifestations' in gtrHash[geneSym].keys():
                manifestations = [x.strip() for x in gtrHash[geneSym]['Manifestations'].split(';') if x != '']
                for organ in osystems:
                    p2file = 0
                    for mani in manifestations:
                        if mani == organ:
                            p2file = 1

                    if p2file == 1:
                        worksheet.write(row, i, 1)
                        count += 1
                        i += 1
                    else:
                        worksheet.write(row, i, 0)
                        i += 1
            else:
                for organ in osystems:
                    worksheet.write(row, i, '-')
                    i += 1

            worksheet.write(row, i, count)

            row += 1

    print_stats(workbook, worksheet1, diseases)
    print_stats(workbook, worksheet2, diseasesOne)

    workbook.close()


def print_stats(workbook, worksheet, diseaseHash):
    '''This function writes the stats for gene tests per disease areas'''

    worksheet.write(0, 0, '#Tests with >= number of genes below:')

    j = 1
    for organ in osystems:
        worksheet.write(0, j, organ + ' (CGD)')
        j += 1

    worksheet.write(0, j, 'Total disease areas (CGD)')

    row = 1
    for test in diseaseHash:
        worksheet.write(row, 0, test)
        total = 0
        k = 1
        for organ in osystems:
            worksheet.write(row, k, diseaseHash[test][organ]['Count'])
            total += diseaseHash[test][organ]['Count']
            k += 1

        worksheet.write(row, k, total)

        row += 1


def main():

    sourceFile0 = 'Gene2GeneIDSyns.txt' #From HGNC genenames.org downloads (Gene Symbol, Synonyms, NCBI GeneID)
    sourceFile1 = 'gene_info_human.txt' #From: https://ftp.ncbi.nih.gov/gene/DATA/gene_info extracted '9606'
    inputFile1 = 'test_version_20181206.txt' #From https://ftp.ncbi.nih.gov/pub/GTR/data/_README.html
    inputFile2 = 'variant_summary.txt' #From ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ 1/29/19
    inputFile3 = 'mimTitles.txt' #From https://www.omim.org/downloads/
    inputFile4 = 'mim2gene_medgen.txt' #From https://www.omim.org/downloads/
    inputFile5 = 'MANE.GRCh38.v0.5.summary.txt' #From ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
    inputFile6 = 'newVCEPGeneList.txt' #Internal file on DropBox Shared U41/ClinVar/ClinVarReports
    inputFile7 = 'ClinGen-Gene-Disease-Summary-2019-02-04.csv' #From https://search.clinicalgenome.org/kb/ - Manually removed header rows
    inputFile8 = 'ClinGen-Dosage-Sensitivity-2019-02-04.csv' #From https://search.clinicalgenome.org/kb/ - Manually removed header rows
    inputFile9 = 'ClinGen-Clinical-Actionability-2019-02-04.csv' #From https://search.clinicalgenome.org/kb/ - Manually removed header rows
    inputFile10 ='CGD.txt' #From https://research.nhgri.nih.gov/CGD/download/
    inputFile11 = 'curations_export_at_2019-02-01_11_52_06.csv' #From Gene Tracker download (needs log in)

    outputFile = 'Gene_Prioritization_Report_Feb_2019.xlsx'

    create_gene2idlookup1(sourceFile0)
    create_gene2idlookup2(sourceFile1)
    create_gtrHash(inputFile1)
    add_VarCounts(inputFile2)
    create_omimHash(inputFile3)
    add_OMIMID(inputFile4)
    add_MANEflag(inputFile5)
    add_VCEPs(inputFile6)
    add_GeneValidity(inputFile7)
    add_GeneDosage(inputFile8)
    add_Actionability(inputFile9)
    add_CGD(inputFile10)
    add_GCEPs(inputFile11)
    get_stats()

    print_allgenes(outputFile)

main()
