__version__ = '1.0'
__date__ = '24/02/2022'
__author__ = 'Juan Ledesma'

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.Data import IUPACData
from itertools import product
import os
import re


def trimming_alignment(BioAln, start_position, end):
    """It takes a Biopython alignment and returns a new alignment trimmed to
    the given starting and ending positions. Gaps (-) introduced in the refseq 
    to get it aligned to the query sequence (insertions) can compromise the 
    final region to trimm. Therefore they are not considered in the numbering
    when trimming is done (if str(BioAln[0].seq[pos]) != '-')"""

    # Precore (Bioaln, 1813, 1900)
    position = start_position # 1813 
    stop_position = 0
    for pos in range(position, len(BioAln[0].seq),1): 
        if position < end:# 1900:
            if str(BioAln[0].seq[pos]) != '-':
                position += 1
                stop_position = pos +1 
    aln = BioAln[:, start_position : stop_position]
    return aln 


def first_aa_nt_in_partial_fragment(aln): 
    """Finds the first nucleotide coding for a amino acid. This is particulary 
    useful for X region, as the amplicon only provide information for the end
    of the gene"""

    if aln: 
        nucleotide = 0
        for nt_pos, (ref_nt, que_nt) in enumerate(zip(aln[0],aln[1])):
            if ref_nt != '-' and que_nt != '-':
                nucleotide = nt_pos 
                break
        aa = 0 
        if nucleotide:
            if (nt_pos)%3 == 0:
                nucleotide = nt_pos
                aa = nt_pos/3
            if (nt_pos+1)%3 == 0:
                nucleotide = nt_pos +1
                aa = (nt_pos+1)/3 
            if (nt_pos+2)%3 == 0:
                nucleotide = nt_pos + 2
                aa = (nt_pos+2)/3  
    return int(aa), int(nucleotide) # the numbering is pythonic, not real one


def frameshift_check(aln, first_codon = 0 ): 
    """ It takes a Biopython alignment and analyses the sequences for frameshifts.
     It translates the query nt sequence at three diferent frames, aligns their
     amino acids to aa reference sequence and finds out if the first one (expected 
     to be in frame after trimming) has the best score. By default it would start
     from the first codon, except for X region"""

    frame_analysis = []
    aaRefSeq = aln[0].seq.ungap().translate() # invariable, it should not have gaps, in-frame...
    ntQueySeq = SeqRecord(seq=aln[1][first_codon:].seq.ungap(), 
                                    id=aln[1].id, name ='')
    scores =[]
    for startpos in range(3):
        #BiopythonWarning: Partial codon, len(sequence) not a multiple of three.
        codons = []
        frame = ntQueySeq.seq[startpos:]
        for n in range(0,len(frame),3):
            codon = str(frame[n:n+3])
            if len(codon) == 3:
                codons.append(codon)
        frame = Seq(''.join(codons)).translate() # frame multiple of 3
        score = pairwise2.align.localxx(frame, aaRefSeq, score_only=True)
        scores.append(score)

    if scores[0] < scores[1] or scores[0] < scores[2]:
        frame_analysis.append('Sequence NOT in frame')

    return frame_analysis


def codon_alignment(aln, region, path_to_alns, path_to_tmp, first_nt = 0): 
    """It takes a Biopython alignment and creates new SeqRecords for amino
    acids and nucleotides using the original record.seq after removing the 
    gaps. As Clustal doesn't deal with ambiguous amino acids when doing 
    protein alignments, the characters BZJUO* are replaced with X. The 
    aa alignment is the used to sort the nucleotides according to the 
    positions of the gaps and amino acids. By default it would start
    from the first first nucleotide coding for an aa, except for X region
    """

    #SeqRecord(seq=Seq(), id='', name ='')
    aa_SeqRecords = [] # aa SeqRecord list
    nt_SeqRecords = [] # nt SeqRecord list
    for record in aln:
        if record.id == aln[0].id:
            aa_record = re.sub(r'[B|Z|J|U|O|*]','X', 
                                str(record.seq.ungap().translate()))
            aa = SeqRecord(seq= Seq(aa_record), id = record.id)
            nt = SeqRecord(seq= record.seq.ungap(), id = record.id)
        else:
            if first_nt == 0:
                aa_record = re.sub(r'[B|Z|J|U|O|*]','X', 
                                    str(record[first_nt:].seq.ungap().translate()))
                aa = SeqRecord(seq= Seq(aa_record), id = record.id)
                nt = SeqRecord(seq= record.seq.ungap(), id = record.id)
            else:
                non_amplified_region = Seq('-'*(first_nt)) # needed for the lenght of aa align
                amplicon = non_amplified_region + record[first_nt:].seq.ungap()
                aa_record = re.sub(r'[B|Z|J|U|O|*]','X', str(amplicon.translate()))
                aa = SeqRecord(seq= Seq(aa_record), id = record.id)
                nt = SeqRecord(seq= record[first_nt:].seq.ungap(), id = record.id)
                # gaps not needed in nt aln,  - from aa aln is translated back as ---
    
        aa_SeqRecords.append(aa)
        nt_SeqRecords.append(nt)

    # the sequences must be same length for the alignment by clustal
    # this may produce a tail of --- in both sequences in the final alignment, 
    # not an issue for the mutation scanning
    if len(aa_SeqRecords[0].seq) < len(aa_SeqRecords[1].seq):
        gap_needed = (len(aa_SeqRecords[1].seq)-len(aa_SeqRecords[0].seq))
        aa_SeqRecords[0].seq = aa_SeqRecords[0].seq + '-'*gap_needed
    if len(aa_SeqRecords[0].seq) > len(aa_SeqRecords[1].seq):
        gap_needed = (len(aa_SeqRecords[0].seq)-len(aa_SeqRecords[1].seq))
        aa_SeqRecords[1].seq = aa_SeqRecords[1].seq  + '-'*gap_needed

    fasta_aa_file = os.path.join(path_to_tmp, f'{aln[1].id}_{region}_aa.fas')
    SeqIO.write(aa_SeqRecords, fasta_aa_file,'fasta')
    
    clustalw_cline = ClustalwCommandline('clustalw2', infile=fasta_aa_file)
    stdout, stderr = clustalw_cline()
    aligned_aa_file = fasta_aa_file.replace('.fas','.aln')
    aa_aln = AlignIO.read(aligned_aa_file,'clustal')
    
    ntAlignRecord = []
    for a in range(len(aa_aln)):
        aligned_nt = []
        pos = 0 
        for n in range(len(aa_aln[0])):
            if aa_aln[a][n] == '-':
                aligned_nt.append('---')
                # gaps do not correspond to a codon in the nt sequence, no increment
            if aa_aln[a][n] != '-':
                codon = str(nt_SeqRecords[a].seq[pos: pos+3])
                aligned_nt.append(codon)
                pos+=3      
        ntAlignRecord.append(SeqRecord(seq=Seq(''.join(aligned_nt)), 
                                                id= nt_SeqRecords[a].id, 
                                                description='') )
    
    alignment = MultipleSeqAlignment(ntAlignRecord)
    output_name = os.path.join(path_to_alns, f'{nt_SeqRecords[1].id}_{region}.fas')
    AlignIO.write(alignment, output_name,'fasta')
    return AlignIO.read(output_name,'fasta')



def translation_of_codons(AlnRecord): 
    """It takes two nucleotide sequences as Biopython SeqRecords and returns two 
    lists of tuples with the codons and the amino acids they code for and the 
    positions. """

    Ref = AlnRecord[0]
    Query = AlnRecord[1]
    codons = [] # list of tuples, [('ATG', 'ATG'),
    aas = []
    pos = 0
    for n in range(0,len(Ref.seq),3):
        pos += 1
        cdref = str(Ref.seq[n:n+3])
        aaref = str(Ref.seq[n:n+3].translate())
        cdquery = str(Query.seq[n:n+3])
        aaquery = str(Query.seq[n:n+3].translate())
        
        if len(cdref) == 3 and len(cdquery) == 3:
            codons.append((cdref, cdquery))
            aas.append((aaref, aaquery, pos))
    # BiopythonWarning: Partial codon, len(sequence) not a multiple of three. 
    # Explicitly trim the sequence or add trailing N before translation. 
    # Codons not multiple of 3 (at the end of the sequences) are excluded to 
    # prevent this potential problem.
   
    return codons, aas



def refseq_numbering(codons, aas, starting_aa=0):
    """ Takes 2 lists of tuples: a first one with the codons for the refseq 
    and the query and a second one with the amino acids for each sequence
    and the positions. It returns a list of insertions in the query sequence 
    and two lists of tuples with codons, amino acids and positions keeping 
    the refseq numbering"""

    codonsNumb = []
    aasNumb = []    
    insertions =[]
    posNumb = starting_aa +1
    # aas[n][0], aas[n][1], aas[n][2]) -> aaref, aaquery, position
    for n in range(starting_aa, len(aas)):
    # insertions 
        if codons[n][0] == '---' and codons[n][1] !='---':
            insertions.append((f'insertion of {aas[n][1]} after {posNumb-1}'))
        else:
            codonsNumb.append((codons[n][0], codons[n][1]))
            aasNumb.append((aas[n][0], aas[n][1], posNumb))#aas[n][2]))
            posNumb +=1 # real positions in Refseq (codons[n][0]) do increment
    return insertions, codonsNumb, aasNumb



def translation_of_ambiguous_dna(seq):
    """returns a string with all the possible amino acids given an ambiguous DNA
    as an input"""
    
    d = IUPACData.ambiguous_dna_values
    codons = []
    for i in product(*[d[j] for j in seq]):
       codons.append("".join(i))
    aa = []
    for codon in codons:
        codonSeq =Seq(codon).translate()
        if codonSeq not in aa: # two dif triplets for the same codon 
            aa.append(str(codonSeq))
    aaSet = '/'.join(aa)
    return aaSet



def scan_mutations(codonsNumb, aasNumb, refseqid): 
    """It takes two list of tuples (codons and amino acids/positions of the
    refseq and the query), compares the amino acids at each position and
    records the variations."""

    mutations = []
    wildtype_aas = {'AY128092.1_Geno_A_surface':{
                                122:['K','R'], 194:['A','V'], 
                                207:['S','N'], 209:['V','L'], 
                                213:['I','L']},
                    'AY167089.1_Geno_B_surface':{
                                44:['G','E'], 101:['Q','H'],
                                122:['K','R'], 126:['T','A'],
                                200:['Y','F'], 213:['M','I']},
                    'KM999990.1_Geno_C_surface':{
                                3:['S','N'], 53:['L','S'], 
                                122:['K','R'], 210:['N','S']},
                    'X65257.1_Geno_D_surface':{
                                118:['T','V','A'], 122:['R','I','K'],
                                125:['T','M'], 127:['P','T'], 
                                128:['A','V'], 207:['S','N','R'], # X DOESNT WORK as single letter
                                224:['V','A']},
                    'X75657.1_Geno_E_surface':{
                                3:['N','G','S'], 49:['L','P'], 
                                57:['T','I'], 59:['N','S'], 
                                184:['A','V']},
                    'KX264497.1_Geno_F_surface':{
                                2:['E','D'], 45:['L','T'], 
                                110:['L','I']},
                    'AY090457.1_Geno_H_surface':{30:['Q','K']},

                    # no modifications for 
                    # AF160501.1_Geno_G, 
                    # AF241411.1_Geno_I
                    # AB486012.1_Geno_J

                    'AY128092.1_Geno_A_polymeraseRT':{ #126:['Y','H'],
                                109:['S','P'], 110:['R','G'], 122:['N','H'],
                                126:['Y'], 129:['M','L'], 131:['N','D'], 
                                153:['W','R'], 163:['V','I'], 217:['L','R'],
                                219:['R','S','A'], 253:['I','V']},
                    'AY167089.1_Geno_B_polymeraseRT':{
                                124:['N','D','H'], 131:['D','N'], 
                                134:['S','N','D'], 187:['L','I'], 238:['H','N'], 
                                253:['V','I']},
                    'KM999990.1_Geno_C_polymeraseRT':{
                                191:['L','I'], 124:['Y','N','H'], 125:['Q','K'],
                                129:['M','A'], 134:['D','N'],139:['N','H'], 
                                149:['K','Q'], 219:['S','A'], 223:['S','A'],
                                224:['I','V'], 238:['N','T','H'] },
                    'X65257.1_Geno_D_polymeraseRT':{
                                91:['L','I'], 114:['R','V'], 122:['F','L','I'],
                                124:['H','N','Y'], 126:['R','H'], 129:['M','L'], 
                                130:['Q','P'], 131:['N','D'], 135:['Y','S'], 
                                145:['L','M'], 149:['Q','K'], 151:['F','Y'], 
                                153:['R','W'], 163:['I','V'], 248:['H','N']},
                    'X75657.1_Geno_E_polymeraseRT':{
                                93:['L','I'], 125:['Q','K'], 130:['P','S'],
                                175:['L','I'], 223:['A','S'], 224:['V','A'],
                                237:['P','S'], 248:['N','H'] },
                    'KX264497.1_Geno_F_polymeraseRT':{
                                118:['T','N'], 123:['D','N'], 129:['L','M'], 
                                137:['S','T'], 149:['Q','K'], 246:['N','T','S'],
                                248:['H','N']},
                    'AY090457.1_Geno_H_polymeraseRT':{
                                122:['Y','N'], 253:['I','V']},
                    'LC360507.1_Geno_C_precore':{
                                5:['P','H'], 11:['S','F'], 13:['S','T']},
                    'LC360507.1_Geno_C_xregion':{
                                101:['S','P'], 102:['A','V'], 116:['L','V'], 
                                118:['K','T'], 119:['D','E'], 144:['S','A'], 
                                146:['A','S'], 147:['P','S'], 152:['P','T']
                                }
                } 
    # codonsNumb = ('ATG', 'ATG') codonsNumb[n][0]-ref, codonsNumb[n][1]-query
    # aasNumb = ('M', 'M', 1) aasNumb[n][0]-ref, aasNumb[n][1]-query, aasNumb[n][2]-positions

    for n in range(len(aasNumb)):
        pos =aasNumb[n][2]
        if aasNumb[n][0] != aasNumb[n][1]:
            if refseqid in wildtype_aas:
                if (pos) in wildtype_aas[refseqid]:
                    if aasNumb[n][1] == 'X':
                        aas = translation_of_ambiguous_dna(codonsNumb[n][1])
                        mutations.append(f'{aasNumb[n][0]}{str(pos)}{aas}')
                    else:
                        if aasNumb[n][1] not in wildtype_aas[refseqid][pos]:
                            mutations.append(f'{aasNumb[n][0]}{str(pos)}{aasNumb[n][1]}')
                elif aasNumb[n][1] == 'X':
                    aas = translation_of_ambiguous_dna(codonsNumb[n][1])
                    mutations.append(f'{aasNumb[n][0]}{str(pos)}{aas}') #ambiguous nt 
                elif aasNumb[n][1] == '-':
                    mutations.append(f'deletion of {aasNumb[n][0]}{str(pos)}')#)deletion')
                else:
                    mutations.append(f'{aasNumb[n][0]}{str(pos)}{aasNumb[n][1]}')
            else: # this  works for Genotypes G, J, I 
                if aasNumb[n][1] == 'X':
                        aas = translation_of_ambiguous_dna(codonsNumb[n][1])
                        mutations.append(f'{aasNumb[n][0]}{str(pos)}{aas}')
                elif aasNumb[n][1] == '-':
                    mutations.append(f'deletion of {aasNumb[n][0]}{str(pos)}')#)deletion')
                else:
                    mutations.append(f'{aasNumb[n][0]}{str(pos)}{aasNumb[n][1]}')
                
        if pos == 269:
               if aasNumb[n][1] != 'L':
                   mutations.append(f'L{str(pos)}{aasNumb[n][1]}')
    return mutations
    


"""
# resistance mutations Polymerase

Tenofovir : S106C, H126Y, D134E, M204I/V, L269I (combination)
Entecavir : I169T
Replication: V173L (IN COMb L180M, M204V/I)
Lamivudine/Emtricitabine: L180M
Lamivudine/Adefovir: A181V/T
Adefovir/Lamivudine/Emtricitabine: A181S
Entecavir: T184G/S/I/L/F/C/A/M (IN COMb L180M, M204V/I)
Tenofovir: A194T (especially when combined with L180M & M204V/I)
Tenofovir: A200V (in the presence of L180M, T184L and M204V)
Entecavir: S202I/G/C (in the presence of M204V/I ± L180M. )
Lamivudine, Emtricitabine, Telbivudine, Entecavir: M204V/I
Adefovir: N236T
Entecavir: M250V/I/L  in the presence of M204V/I & L180M.

No conflicts with Genotypes G,H,I and J

#Surface:
codons 120-150

#Precore
Insertions, deletions, frameshift , M1, stop codons

#Xregion
K130M and V131I

dictionary = {106}
 

"""
