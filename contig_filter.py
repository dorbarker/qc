from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import os
import argparse

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--kmer', type = int)
    
    parser.add_argument('--fasta-dir')
    
    parser.add_argument('--organism', required = True,
                        help = 'Name of the organism the contigs \
                                *should* be from.')

    parser.add_argument('--gc-cutoff', nargs = 2, type = float,
                        required = True, metavar = ('LOW', 'HIGH'), 
                        help = 'GC contents outside this range \
                                will be queried against GenBank')

    parser.add_argument('--fragment', type = int, default = 100,
                        help = 'Size (bp) of contig fragment used to query \
                                GenBank (default = 100)')
    
    parser.add_argument('--hits', type = int, default = 50,
                        help = 'Number top hits to consider')

    return parser.parse_args()

def gc_content(seq):

     return round(sum(1.0 for x in seq if x in ('G','C')) / len(seq), 5)

def scan_contig(contig, kmer_width):

    gc = []

    for start in range(len(contig)):
        
        gc.append( gc_content(contig[start : start + kmer_width]) )

    
    return sum(gc) / len(gc)

def select_checkables(fasta_dir, fragment, low, high):

    checkables = {}

    for fasta in (os.path.join(fastq_dir, x) for x in os.listdir(fasta_dir)):
        
        checkables[fasta] = {}

        with open(fasta, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                
                name = rec.id
                s    = str(rec.seq)[:fragment]
                gc   = scan_contig(s)

                if not low <= gc <= high:

                    checkables[fasta][name] = s

    return checkables

def query_genbank(checkable, organism, hits):

    axe = {}

    for fasta in checkable:
        axe[fasta] = []
        for contig in fasta:
            handle = NCBIWWW.qblast('blastn', 'nt', checkable[fasta][contig],
                     megablast = True, hitlist_size = hits, alignments = hits)

            record = NCBIXML.read(handle)

            cj_in_title = []

            for aln in record.alignments:
                for hsp in aln.hsps:
                    cj_in_title.append(organism.lower() in aln.title.lower())
                            
            if not any(cj_in_title):
                axe[fasta].append(contig)

    return axe

def remove_bad_contigs(axe):

    for fasta in axe:
        towrite = []
        with open(fasta, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                
                name = rec.id

                if name not in axe[fasta]:
                    towrite.append(rec)
        
        with open(fasta, 'w') as j:

            SeqIO.write(towrote, j, 'fasta')


def main():

    args = arguments()

if __name__ == '__main__':
    main()
