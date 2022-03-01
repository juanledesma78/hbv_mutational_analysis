__version__ = '1.0'
__date__ = '24/02/2022'
__author__ = 'Juan Ledesma'

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from glob import glob
import os
import shutil
import pandas as pd
import numpy as np
import mutation_utilities as mu
import argparse
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

program_specs = argparse.ArgumentParser(
    prog='hbv_mutation_genotyping.py',
    usage='%(prog)s -dir <fasta file location>',
    description='The script %(prog)s takes as an argument the RUN location where\
                a folder "Consensus" contains the fasta sequences to be analysed\
                and a csv file with their genotyping results are located.\
                It returns an excel report with the mutations found for each\
                sequence, a log with details about the sequences analysed and\
                alignmnets for the user to inspect if needed.')
program_specs.add_argument('-dir', 
                            help='path to the directory where the fasta files are\
                                  This directory MUST contain a csv file with the genotyping results',
                            required=True)
args = program_specs.parse_args()

path_to_dir = os.path.join(os.getcwd(), args.dir) #E1647/Consensus
path_to_input = os.path.join(path_to_dir, 'Consensus')

if os.path.isdir(path_to_input):
    path_to_alns = os.path.join(path_to_input, 'alignments')
    path_to_trimmed = os.path.join(path_to_alns, 'raw_trimmed')
    path_to_tmp = os.path.join(path_to_input, 'tmp')
    run_name = args.dir.split('/')[0]
    if os.path.exists(path_to_tmp): #it will remove tmp in uncomplete analysis
        shutil.rmtree(path_to_tmp) 
    if os.path.exists(path_to_alns):# and os.path.exists(path_to_tmp):
        shutil.rmtree(path_to_alns)
        os.mkdir(path_to_alns)
        os.mkdir(path_to_trimmed) # for checking not inframe sequences
        os.mkdir(path_to_tmp) #used for temp files
    else:
        os.mkdir(path_to_alns)
        os.mkdir(path_to_trimmed) 
        os.mkdir(path_to_tmp) 

else: 
    print(f'''The input {args.dir} is NOT a DIRECTORY.
Please enter the proper path to the RUN directory.
    ''')
    raise NotADirectoryError(path_to_input)
    

""" import the reference sequences align them to each consensus sequence 
    in individual files"""

xprecore_refseq = SeqIO.read('refseqs/LC360507.1_WG_xprecore.fas','fasta')
spol_refseqs = {    'A':'AY128092.1_Genotype_A_polymerase_surface', 
                    'B':'AY167089.1_Genotype_B_polymerase_surface',
                    'C':'KM999990.1_Genotype_C_polymerase_surface',
                    'D':'X65257.1_Genotype_D_polymerase_surface',
                    'E':'X75657.1_Genotype_E_polymerase_surface',
                    'F':'KX264497.1_Genotype_F_polymerase_surface',
                    'G':'AF160501.1_Genotype_G_polymerase_surface',
                    'H':'AY090457.1_Genotype_H_polymerase_surface',
                    'I':'AF241411.1_Genotype_I_polymerase_surface',
                    'J':'AB486012.1_Genotype_J_polymerase_surface'
                        }

try:
    csv_file = os.path.join(path_to_input, 'Sample_list.csv')
    genotype_results = pd.read_csv(csv_file)
    if genotype_results.iloc[0]['Seq ID'] and genotype_results.iloc[0]['Genotype']:

        
        """It takes a fasta file and the corresponding reference sequence above
        and creates a SeqRecord that will be used to do a initial pair alignment 
        using CLustalw, only for trimming the amplicons at the proper length"""

        key = 0
        id_map ={} # to prevent issues with long IDs for the sequences
        report ={}
        xc_fasta_input = [] # log
        sp_fasta_input = [] # log
        csv_seq_list = []
        
        for fasta_file in glob(os.path.join(path_to_input,'*.fas')):
            sequence= SeqIO.read(fasta_file, "fasta")
            name = (sequence.id)[:-2]
            key +=  1
            id_map[str(key)] = name
            report[name] = {'surface':{}, 'polymeraseRT':{}, 'xregion':{}, 'precore':{}}    

            ''' XCORE AMPLICONS '''

            if 'xc' in sequence.id:
                xc_fasta_input.append(sequence.id)
                xc_SeqRecord = []
                xc_SeqRecord.append(xprecore_refseq)
                sequence.id = str(key)
                xc_SeqRecord.append(sequence)
                output_name = os.path.join(path_to_tmp, f'{name}_xc_refseq.fas')
                SeqIO.write(xc_SeqRecord, output_name,'fasta')
               
                '''SPOL AMPLICONS'''
            else:
                sp_fasta_input.append(sequence.id)
                for row in range(len(genotype_results)):
                    csv_seq_list.append(genotype_results.iloc[row]['Seq ID']) # list of IDs in the csv
                    
                    """ Selection of refseq according to genotype"""
                    if name in genotype_results.iloc[row]['Seq ID']:
                        genotype = (genotype_results.iloc[row]['Genotype'])
                        if genotype[0] in spol_refseqs:
                            spol_refseq = SeqIO.read(f'refseqs/{spol_refseqs[genotype[0]]}.fas','fasta')
                            spol_SeqRecord = []
                            spol_SeqRecord.append(spol_refseq)
                            sequence.id = str(key)
                            spol_SeqRecord.append(sequence)
                            output_name = os.path.join(path_to_tmp, f'{name}_spol_refseq.fas')
                            SeqIO.write(spol_SeqRecord, output_name,'fasta') 

        for fasta in glob(os.path.join(path_to_tmp, '*refseq.fas')):    
            #clustalw_exe = "/usr/bin/clustalw" # Windows users
            #clustalw_cline = ClustalwCommandline(clustalw_exe, infile=fasta)
            #assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
            clustalw_cline = ClustalwCommandline('clustalw2', infile=fasta)
            stdout, stderr = clustalw_cline()

        """It takes the alignments and trims them to the length of the region
        according to the given coordinates"""

        for aln in glob(os.path.join(path_to_tmp, '*refseq.aln')): 

            ''' Analysis of xcore amplicons '''

            if 'xc' in aln:
                xcore_aln_SeqRecords = AlignIO.read(aln, 'clustal')
                xcore_aln_SeqRecords[1].id = id_map[str(xcore_aln_SeqRecords[1].id)]  
                x_aln = mu.trimming_alignment(xcore_aln_SeqRecords, 1374-1, 1838)
                x_aln[0].id = str(x_aln[0].id).replace('WG','xregion') 
                x_aln[0].description = ''
                x_aln[1].description = 'x_region'
                x_alignment = MultipleSeqAlignment(x_aln) 
                output_name = os.path.join(path_to_trimmed, 
                                            f'{x_aln[1].id}_x_region_raw.fas')
                AlignIO.write(x_alignment, output_name,'fasta')

                precore_aln = mu.trimming_alignment(xcore_aln_SeqRecords, 1814-1, 1900)
                precore_aln[0].id = str(precore_aln[0].id).replace('WG','precore') 
                precore_aln[0].description = ''
                precore_aln[1].description = 'precore'
                precore_alignment = MultipleSeqAlignment(precore_aln)
                output_name = os.path.join(path_to_trimmed, 
                                            f'{precore_aln[1].id}_precore_raw.fas')
                AlignIO.write(precore_alignment, output_name,'fasta')

                """Analysis of mutations; cheking the frame, alignment of the
                sequences by codons, translation and detection of insertion
                keeping the numbering of the refseq and recording of the 
                mutations """

                if precore_aln:
                    sampleid = precore_aln[1].id 
                    region = (precore_aln[0].id).split('_')[-1]
                    refseqid = precore_aln[0].id

                    frame = mu.frameshift_check(precore_aln)
                    report[sampleid][region]['frameshift'] = frame
                    
                    if len(frame) == 0:
                        codonAlnRecord = mu.codon_alignment(precore_aln,
                                                            region, 
                                                            path_to_alns, 
                                                            path_to_tmp)

                        codons, aas = mu.translation_of_codons(codonAlnRecord)

                        insertions,\
                        codonsNumb, aasNumb = mu.refseq_numbering(codons, aas)
                        report[sampleid][region]['insertions'] = insertions
                        
                        mutations = mu.scan_mutations(codonsNumb, aasNumb, 
                                                    refseqid)
                        report[sampleid][region]['mutations'] = mutations   

                if x_aln:
                    sampleid = (x_aln[1].id)
                    region = (x_aln[0].id).split('_')[-1]
                    refseqid = x_aln[0].id

                    """ only final region for this gene"""
                    starting_aa,\
                    starting_nt = mu.first_aa_nt_in_partial_fragment(x_aln)

                    frame = mu.frameshift_check(x_aln, first_codon = starting_nt)
                    report[sampleid][region]['frameshift'] = frame

                    if len(frame) == 0:
                            codonAlnRecord = mu.codon_alignment(x_aln,
                                                                region, 
                                                                path_to_alns, 
                                                                path_to_tmp, 
                                                                starting_nt)

                            codons, aas = mu.translation_of_codons(codonAlnRecord)

                            insertions,\
                            codonsNumb, aasNumb = mu.refseq_numbering(codons, aas, 
                                                                    starting_aa)
                            report[sampleid][region]['insertions'] = insertions

                            mutations = mu.scan_mutations(codonsNumb, aasNumb, 
                                                        refseqid)
                            report[sampleid][region]['mutations'] = mutations         

            ''' Analysis of spol amplicons '''

            if 'spol' in aln: 
                genotype_refseqs = { # 2076 (Whole RT) 1851 (to cover position 269)
                        'AY128092.1_Geno_A_polymerase':{'surface': [1070-1,1750],
                                                            'rt':[1045-1,1851]}, #2076]}, # 1851
                        'AY167089.1_Geno_B_polymerase':{'surface': [1064-1,1744], 
                                                            'rt':[1039-1,1845]}, #2070]}, # 1845
                        'KM999990.1_Geno_C_polymerase':{'surface': [1064-1,1744], 
                                                            'rt':[1039-1,1845]}, #2070]},# 1845
                        'X65257.1_Geno_D_polymerase':{'surface': [1031-1,1711], 
                                                            'rt':[1006-1,1812]}, #2037]},# 1812
                        'X75657.1_Geno_E_polymerase':{'surface': [1061-1,1741], 
                                                            'rt':[1036-1,1842]}, #2067]}, # 1842
                        'KX264497.1_Geno_F_polymerase':{'surface': [1064-1,1744], 
                                                            'rt':[1039-1,1845]}, #2070]},# 1845
                        'AF160501.1_Geno_G_polymerase':{'surface': [1061-1,1741], 
                                                            'rt':[1036-1,1842]}, #2067]} #1842
                        'AY090457.1_Geno_H_polymerase':{'surface': [1064-1,1744], 
                                                            'rt':[1039-1,1845]}, #2070]},# 1845
                        'AF241411.1_Geno_I_polymerase':{'surface': [1064-1,1744], 
                                                            'rt':[1039-1,1845]}, #2070]} #1845
                        'AB486012.1_Geno_J_polymerase':{'surface': [1031-1,1711], 
                                                            'rt':[1006-1,1812]} #2037]} #1812 
                                    }

                spol_aln_SeqRecords = AlignIO.read(aln, 'clustal')
                spol_aln_SeqRecords[1].id = id_map[str(spol_aln_SeqRecords[1].id)]
                
                """Trimming according to the coordinates of each reference sequences"""

                if spol_aln_SeqRecords[0].id in genotype_refseqs:

                    start_surface = genotype_refseqs[spol_aln_SeqRecords[0].id]['surface'][0]
                    end_surface = genotype_refseqs[spol_aln_SeqRecords[0].id]['surface'][1]
                    surface_aln = mu.trimming_alignment(spol_aln_SeqRecords, start_surface, end_surface)
                    surface_aln[0].id = str(surface_aln[0].id).replace('polymerase', 'surface') 
                    surface_aln[0].description = ''
                    surface_aln[1].description = 'surface'
                    surface_alignment = MultipleSeqAlignment(surface_aln)
                    output_name = os.path.join(path_to_trimmed, f'{surface_aln[1].id}_surface_raw.fas')
                    AlignIO.write(surface_alignment, output_name,'fasta')

                    start_rt = genotype_refseqs[spol_aln_SeqRecords[0].id]['rt'][0]
                    end_rt = genotype_refseqs[spol_aln_SeqRecords[0].id]['rt'][1]
                    pol_RT_aln = mu.trimming_alignment(spol_aln_SeqRecords, start_rt, end_rt)
                    pol_RT_aln[0].id = str(pol_RT_aln[0].id).replace('polymerase', 'polymeraseRT')
                    pol_RT_aln[0].description = ''
                    pol_RT_aln[1].description = 'pol RT'
                    polrt_alignment = MultipleSeqAlignment(pol_RT_aln)
                    output_name = os.path.join(path_to_trimmed, f'{pol_RT_aln[1].id}_polymeraseRT_raw.fas')
                    AlignIO.write(polrt_alignment, output_name,'fasta')

                    if surface_aln:
                        sampleid = (surface_aln[1].id)
                        region = (surface_aln[0].id).split('_')[-1]
                        refseqid = surface_aln[0].id
                        
                        frame = mu.frameshift_check(surface_aln)
                        report[sampleid][region]['frameshift'] = frame
                        
                        if len(frame) == 0:
                            codonAlnRecord = mu.codon_alignment(surface_aln,
                                                                region, 
                                                                path_to_alns, 
                                                                path_to_tmp)

                            codons, aas = mu.translation_of_codons(codonAlnRecord)
                            
                            insertions,\
                            codonsNumb, aasNumb = mu.refseq_numbering(codons, aas)
                            report[sampleid][region]['insertions'] = insertions
                            
                            mutations = mu.scan_mutations(codonsNumb, aasNumb, 
                                                        refseqid)
                            report[sampleid][region]['mutations'] = mutations                            
                            
                    if pol_RT_aln:
                        sampleid = (pol_RT_aln[1].id)
                        region = (pol_RT_aln[0].id).split('_')[-1] # polymeraseRT
                        refseqid = pol_RT_aln[0].id
                        
                        frame = mu.frameshift_check(pol_RT_aln)
                        report[sampleid][region]['frameshift'] = frame
                        
                        if len(frame) == 0:
                            codonAlnRecord = mu.codon_alignment(pol_RT_aln,
                                                                region, 
                                                                path_to_alns, 
                                                                path_to_tmp)

                            codons, aas = mu.translation_of_codons(codonAlnRecord)
                            
                            insertions,\
                            codonsNumb, aasNumb = mu.refseq_numbering(codons, aas)
                            report[sampleid][region]['insertions'] = insertions
                            
                            mutations = mu.scan_mutations(codonsNumb, aasNumb, 
                                                        refseqid)
                            report[sampleid][region]['mutations'] = mutations    

        """Create a dataframe to use for openpyxl"""
        
        df_results = pd.DataFrame()
        sortedreport = dict(sorted(report.items()))
        row = {'Notes':''}

        csv_molis_ids = [gr.split('_')[1] for gr in genotype_results['Seq ID']]
        csv_tube_numbers = [gr.split('_')[0] for gr in genotype_results['Seq ID']]
        genotype_results['Molis No.'] = csv_molis_ids 
        genotype_results['Tube No.'] = csv_tube_numbers
        
        for seqid, geneScan in sortedreport.items():
            row['Tube No.'] = seqid.split('_')[0]
            row['Molis No.'] = seqid.split('_')[1]  
            for gene in geneScan:
                if geneScan[gene] != {}: 
                    if len(geneScan[gene]['frameshift']) != 0:
                        info =', '.join(geneScan[gene]['frameshift'])
                        row[gene] = info
                    else:
                        if len(geneScan[gene]['insertions']) == 0 and \
                            len(geneScan[gene]['mutations']) == 0:
                            row[gene] = "No mutations"
                        if len(geneScan[gene]['insertions']) == 0 and \
                            len(geneScan[gene]['mutations']) != 0:
                            info =', '.join(geneScan[gene]['mutations'])
                            row[gene] = info
                        if len(geneScan[gene]['insertions']) != 0 and \
                            len(geneScan[gene]['mutations']) != 0:
                            ins_mut =   geneScan[gene]['insertions'] +\
                                        geneScan[gene]['mutations']
                            info =', '.join(ins_mut)
                            row[gene] = info
                else:
                    row[gene] = 'NA'
            df_results = df_results.append(row, ignore_index = True)
        
        combined_df_results = pd.merge(genotype_results, df_results, 
                                        how='left',
                                        on=['Molis No.','Tube No.']) 
        combined_df_results = combined_df_results.replace({np.nan: 'NA' })

        combined_df_results.rename(columns={'surface':'Surface',
                                            'polymeraseRT': 'Polymerase (RT)', 
                                            'xregion': 'X region', 
                                            'precore': 'Precore' }, inplace= True)


        combined_df_results = combined_df_results[['Tube No.', 'Molis No.', 
                                                    'Genotype', 'Surface',
                                                    'Polymerase (RT)', 
                                                    'X region', 'Precore', 
                                                    'Notes']]

        workbook = Workbook()
        ws = workbook.active
        ws1 = workbook.create_sheet("All", 0)
        ws2 = workbook.create_sheet("A", 1)
        ws3 = workbook.create_sheet("B", 2)
        ws4 = workbook.create_sheet("C", 3)
        ws5 = workbook.create_sheet("D", 4)
        ws6 = workbook.create_sheet("E", 5)
        ws7 = workbook.create_sheet("F", 6)
        ws8 = workbook.create_sheet("G", 7)
        ws9 = workbook.create_sheet("H", 8)
        ws10 = workbook.create_sheet("I", 9)
        ws11 = workbook.create_sheet("J", 10)
        ws12 = workbook.create_sheet("CORE", 11)
        ws13 = workbook.create_sheet("Numbering", 12)
        ws14 = workbook.create_sheet("Extended_IUPAC_aa_code", 13)

        for rowdf in dataframe_to_rows(combined_df_results, index=False, 
                                                            header=True):
            ws1.append(rowdf)

        genotypes = ['A','B','C','D','E','F','G','H','I','J']
        worksheets = [ws2, ws3, ws4, ws5, ws6, ws7, ws8, ws9, ws10, ws11, ws12]
        header = ['Tube No.', 'Molis No.', 'Genotype', 'Surface',
                    'Polymerase (RT)', 'X region', 'Precore', 'Notes']
        for ws in worksheets:
            ws.append(header)
        for rowdf in dataframe_to_rows(combined_df_results, index=False, 
                                                            header=False):
            for index, (gt, ws) in enumerate(zip(genotypes, worksheets)):
                if gt in rowdf[2][0]: # Genotype first letter A,G, use np for NaN issue
                    ws.append(rowdf)
            if rowdf[6]!= 'NA': #CORE
                ws12.append(rowdf)

        # numbering worksheet
        acc_list = []
        target_list = []
        for s in spol_refseqs.values():
            acc_number = s.split('_')[0]
            target = s.split('_')[1:]
            target =(' '.join(target)).replace('Genotype ','' )
            acc_list.append(acc_number)
            target_list.append(target)
            
        acc_list.append(xprecore_refseq.id.split('_')[0])
        xprecore = (' '.join(xprecore_refseq.id.split('_')[1:]).
                        replace('Geno C WG','Xregion precore'))
        target_list.append(xprecore)#.replace('WG', 'xregion precore'))
        data = {'Accession Number': acc_list, 
                'Genotype and region for numbering': target_list}
        df_numbering = pd.DataFrame(data)
        
        for r in dataframe_to_rows(df_numbering, index=False, header=True):
            ws13.append(r)

        # amino acid codes worksheet
        aa_code_info = {'Symbol (combination)':['A','R','N','D','C','Q','E',
                                                'G','H','I','M','L','K','F',
                                                'P','S','T','W','Y','V','*',
                                                'B (N/D)','X','Z (E/Q)',
                                                'J (L/I)','U','O'],
                        'Amino acid':['Alanine', 'Arginine','Asparagine',
                                    'Aspartic acid','Cysteine', 'Glutamine',
                                    'Glutamic acid','Glycine','Histidine',
                                    'Isoleucine','Methionine','Leucine',
                                    'Lysine','Phenylalanine','Proline',
                                    'Serine','Threonine','Tryptophan',
                                    'Tyrosine','Valine','STOP',
                                    'Asparagine or Aspartic acid',
                                    'Unknown or other amino acid',
                                    'Glutamic acid or Glutamine',
                                    'Leucine or Isoleucine',
                                    'Selenocysteine','Pyrrolysine']
                        }
        df_amino_acid_code = pd.DataFrame(aa_code_info)
        
        for r in dataframe_to_rows(df_amino_acid_code, index=False, header=True):
            ws14.append(r)

        excel_file = os.path.join(path_to_dir, 
                                f'Results_{run_name}_hbvma_v{__version__}.xlsx') 
        workbook.save(excel_file)


    shutil.rmtree(path_to_tmp) #Comment it out for debugging if needed 

    """ Generation of log"""

    xc_fasta_input = sorted(xc_fasta_input)
    sp_fasta_input = sorted(sp_fasta_input)
    xc_ids = [n.replace('xc','') for n in xc_fasta_input]
    sp_ids = [n.replace('sp','') for n in sp_fasta_input]
    csv_seq_list = sorted(set(csv_seq_list))
    
    xc_fasta_ids = sorted(set(xc_ids))
    sp_fasta_ids = sorted(set(sp_ids))

    not_in_csv =[]
    sp_fasta_not_provided =[]
    xc_fasta_not_provided =[]
    analysed_sp_fasta = []
    analysed_xc_fasta = []
    
    for spname in sp_fasta_ids:
        if spname in csv_seq_list:
            analysed_sp_fasta.append(spname)
        if spname not in csv_seq_list and spname not in not_in_csv:
            not_in_csv.append(spname)

    for xcname in xc_fasta_ids:
        if xcname in csv_seq_list:
            analysed_xc_fasta.append(xcname)
        if xcname not in csv_seq_list and spname not in not_in_csv:
            not_in_csv.append(xcname)
      
    log_name = os.path.join(path_to_dir,f'{run_name}_hbvma_v{__version__}_log.txt')
    
    with open(log_name, 'w') as log:
        info = f'''
    --------------------------------------------------------------------------
            HBV MUTATIONAL ANALYSIS v{__version__} - RUN {run_name} LOG
    --------------------------------------------------------------------------

    --- INPUT ----------------------------------------------------------------
        
    The CSV file contained the following list of samples (and their genotypes if available):


    {', '.join(csv_seq_list)}


    The following FASTA files were provided for the analysis of 
    XCORE amplicons:


    {', '.join(xc_fasta_input)}


    and 
    
    SPOL amplicons:


    {', '.join(sp_fasta_input)}


    --- OUTPUT ---------------------------------------------------------------

    The folowing samples were ANALYSED for 
    
    X REGION and PRECORE:


    {', '.join(analysed_xc_fasta)} 


    and 
    
    for SURFACE and POLYMERASE (RT domain):


    {', '.join(analysed_sp_fasta)}


    and the results have been recorded in {excel_file}. 


    
    <<< WARNING >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    No record was found in the CSV file for the samples 

    {', '.join(not_in_csv)}

    so they were EXCLUDED from the analysis.

    '''
        log.write(info)

    print(f'''
    ANALYSIS COMPLETE
    Please check the details recorded in the log saved at {log_name} 
    before looking at the results.
    Results are available at {excel_file}.
    The alignments used to generate the excel report can be inspected at {path_to_alns}.
    ''')

except FileNotFoundError:
    print(f'''
    The file Sample_list.csv is not in the directory {path_to_input}.
    Please make sure to save it in the same directory where the FASTA sequences are.
    ''')


except KeyError:
    print(f'''
    Please make sure:

    - the FASTA files have been provided in {path_to_input}
    - the header names in the Sample_list.csv in the directory 
      {path_to_input} are labelled as 'Seq ID' and 'Genotype'.
    - the FASTA sequences are labelled according to the convention
        tube_molisAmplicon (i.e. 01_H12345678sp)
    - the sequence id in the CSV file matches the identifier of the FASTA 
      files (either amplicon for the same sample). The name should be the 
      same after excluding sp/xc from the idenfifier of the FASTA file. 

    Example:        "CSV file"        
                Seq ID              Genotype
                01_H12345678         A1
    
                    "FASTA file"
    (correct)   >01_H12345678sp      
                ATGTGTGGTGAT...
                 
    (wrong)     >01_H12345678_sp      
                ATGTGTGGTGAT...
                    ''')
