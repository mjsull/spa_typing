import sys
import string
import argparse
import os
from itertools import groupby
import urllib


def revseq(seq):
    transtab = string.maketrans('atcgATCG', 'tagcTAGC')
    seq = seq[::-1]
    seq = seq.translate(transtab)
    return seq


def fasta_dict(fasta_name):
    """
    given a fasta file. yield dict of header, sequence
    """
    seqDict = {}
    with open(fasta_name) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = header.next()[1:].strip()
            seq = "".join(s.strip() for s in faiter.next())
            if header in seqDict:
                sys.exit('FASTA contains multiple entries with the same name')
            else:
                seqDict[header] = seq
    return seqDict

def enrichSeq(seq, fortemp, revtemp):
    index = 0
    found = None
    out1 = []
    out_list = []
    while found != -1:
        found = seq.find(fortemp, index)
        index = found + 1
        if found != -1:
            out1.append(found)
    index = 0
    found = None
    out2 = []
    revtempr = revseq(revtemp)
    while found != -1:
        found = seq.find(revtempr, index)
        index = found + 1
        if found != -1:
            out2.append(found)
    for j in out1:
        for k in out2:
            if k - j > 50 and k - j < 5000:
                enriched_seq = seq[j:k+len(revtemp)]
                out_list.append(enriched_seq)
    index = 0
    found = None
    out1 = []
    fortempr = revseq(fortemp)
    while found != -1:
        found = seq.find(fortempr, index)
        index = found + 1
        if found != -1:
            out1.append(found)
    index = 0
    found = None
    out2 = []
    while found != -1:
        found = seq.find(revtemp, index)
        index = found + 1
        if found != -1:
            out2.append(found)
    for j in out2:
        for k in out1:
            if k - j > 50 and k - j < 5000:
                enriched_seq = seq[j:k+len(fortemp)]
                out_list.append(revseq(enriched_seq))
    return out_list


def getSpaTypes(reps, orders):
    rep_dir = os.path.dirname(os.path.realpath(__file__))
    if reps is None:
        reps = os.path.join(rep_dir, 'sparepeats.fasta')
        if not os.path.exists(reps):
            urllib.urlretrieve('http://spa.ridom.de/dynamic/sparepeats.fasta', reps)
        if not os.path.exists(reps):
            sys.exit('Could not download http://spa.ridom.de/dynamic/sparepeats.fasta, download manually and use -r flag')
    if orders is None:
        orders = os.path.join(rep_dir, 'spatypes.txt')
        if not os.path.exists(orders):
            urllib.urlretrieve('http://spa.ridom.de/dynamic/spatypes.txt', orders)
        if not os.path.exists(orders):
            sys.exit('Could not download http://spa.ridom.de/dynamic/spatypes.txt, download manually and use -p flag')
    seqDict = {}
    letDict = {'58': 'B4', '30': 'O2', '54': 'H3', '42': 'M2', '48': 'V2', '45': 'A3', '43': 'X2',
               '60': 'S2', '61': 'W3', '62': 'U3', '57': 'S', '64': 'X3', '49': 'Y2', '66': 'F4',
               '90': 'I', '68': 'E4', '69': 'C4', '80': 'K4', '52': 'R3', '53': 'G3', '02': 'A',
               '03': 'D2', '26': 'T', '01': 'XX', '06': 'G2', '07': 'U', '04': 'Z', '05': 'C',
               '46': 'Y3', '47': 'Z3', '08': 'X', '09': 'A2', '28': 'R', '29': 'F2', '41': 'U2',
               '14': 'I2', '59': 'T3', '78': 'J4', '51': 'P2', '24': 'Q', '56': 'J2', '25': 'O',
               '39': 'E3', '65': 'S3', '76': 'K3', '75': 'I4', '38': 'F3', '73': 'G4', '72': 'P3',
               '71': 'Q3', '70': 'D4', '20': 'D', '74': 'H4', '21': 'F', '11': 'Y', '10': 'C2',
               '13': 'E', '12': 'G', '15': 'W', '22': 'L', '17': 'M', '16': 'K', '19': 'H', '18': 'H2',
               '31': 'N', '23': 'J', '37': 'D3', '36': 'W2', '35': 'C3', '34': 'B', '33': 'P', '55': 'A4',
               '63': 'V3', '32': 'E2', '44': 'Z2', '50': 'T2'}
    typeDict = {}
    seqLengths = set()

    reps_dict = fasta_dict(reps)
    for i in reps_dict:
        seq = reps_dict[i]
        num = i[1:]
        seqDict[seq.upper()] = num
        seqLengths.add(len(seq))
    with open(orders) as f:
        for line in f:
            st, pattern = line.rstrip().split(',')
            typeDict[pattern] = st
    return seqDict, letDict, typeDict, seqLengths


def findPattern(infile, seqDict, letDict, typeDict, seqLengths):
    qDict = fasta_dict(infile)
    seq_list = []
    for i in qDict:
        enriched_seqs = enrichSeq(qDict[i].upper(), 'TAAAGACGATCCTTCGGTGAG', 'CAGCAGTAGTGCCGTTTGCTT')
        seq_list += enriched_seqs
    if seq_list == []:
        for i in qDict:
            enriched_seqs = enrichSeq(qDict[i].upper(), 'AGACGATCCTTCGGTGAGC', 'GCTTTTGCAATGTCATTTACTG')
            seq_list += enriched_seqs
    if seq_list == []:
        for i in qDict:
            enriched_seqs = enrichSeq(qDict[i].upper(), 'ATAGCGTGATTTTGCGGTT', 'CTAAATATAAATAATGTTGTCACTTGGA')
            seq_list += enriched_seqs
    if seq_list == []:
        for i in qDict:
            enriched_seqs = enrichSeq(qDict[i].upper(), 'CAACGCAATGGTTTCATCCA', 'GCTTTTGCAATGTCATTTACTG')
            seq_list += enriched_seqs
    if seq_list == []:
        return ['no enriched sequence.']
    if len(seq_list) > 1:
        sys.stderr.write(' more than one enriched sequence in ' + infile + '\n')
    rep_list = []
    for i in seq_list:
        index = 0
        adjacent = False
        rep_order = []
        while index <= len(i):
            gotit = False
            for j in seqLengths:
                if i[index:index+j] in seqDict:
                    if adjacent or rep_order == []:
                        rep_order.append(seqDict[i[index:index+j]])
                    else:
                        rep_order.append('xx')
                        rep_order.append(seqDict[i[index:index+j]])
                    index += j
                    gotit = True
                    adjacent = True
                    break
            if not gotit:
                index += 1
                adjacent = False
        rep_list.append(rep_order)
    out_list = []
    for i in rep_list:
        let_out = ''
        for j in i:
            if j in letDict:
                let_out += letDict[j] + '-'
            else:
                let_out += 'xx-'
        let_out = let_out[:-1]
        if '-'.join(i) in typeDict:
            type_out = typeDict['-'.join(i)]
        else:
            type_out = '-'.join(i)
        out_list.append(let_out)
        out_list.append(type_out)
    return out_list

parser = argparse.ArgumentParser(prog='get_spa_type.py', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
get_spa_type.py: Get spa types

Version: 0.1.0
License: GPLv3

USAGE: python get_spa_type.py fasta_file.fa
Prints spa type to stdout

Will download sparepeats.fasta and spatypes.txt to repository directory if files not provided or already in directory.



''', epilog="written by mjsull.")


parser.add_argument('-r', '--repeat_file', action='store', help='List of spa repeats (http://spa.ridom.de/dynamic/sparepeats.fasta)')
parser.add_argument('-o', '--repeat_order_file', action='store', help='List spa types and order of repeats (http://spaserver2.ridom.de/dynamic/spatypes.txt)')
parser.add_argument('-f', '--fasta', action='store', help='Fasta file input.')
parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

args = parser.parse_args()

seqDict, letDict, typeDict, seqLengths = getSpaTypes(args.repeat_file, args.repeat_order_file)
the_out = findPattern(args.fasta, seqDict, letDict, typeDict, seqLengths)
sys.stdout.write('Spa type:\t' + '\t'.join(the_out) + '\n')
