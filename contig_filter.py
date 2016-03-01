from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from time import sleep
import os
import argparse

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--fasta-dir', required = True,
                        help = 'Directory containing FASTA files')
    
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
                        help = 'Number top hits to consider (default = 50)')

    return parser.parse_args()

def gc_content(seq):
     '''Returns GC content of seq'''

     return sum(1.0 for x in seq if x in ('G','C')) / len(seq)

def select_checkables(fasta_dir, fragment, low, high):
    '''Checks each contig in each FASTA for GC content.

    Out-of-range GC contents are later checked against GenBank.
    '''

    checkables = {}

    for fasta in (os.path.join(fastq_dir, x) for x in os.listdir(fasta_dir)):
        
        checkables[fasta] = {}

        with open(fasta, 'r') as f:
            
            records = 0.0

            for rec in SeqIO.parse(f, 'fasta'):
                
                records += 1.0

                name = rec.id
                s    = str(rec.seq)
                gc   = gc_content(s)

                if not low <= gc <= high:

                    checkables[fasta][name] = s[:fragment]
            
            # If 3/4 contigs have the wrong GC, then the whole isolate
            # will be discarded at a later step and there's
            # no need to make all those BLAST searches

            if records / len(checkables[fasta]) >= 0.75:
                del checkables[fasta]

    return checkables

def query_genbank(checkable, organism, hits):
    '''Any contigs identified in select_checkables() are queried against
    NCBI GenBank using a fragment from the beginning of the contig. 
    
    If the `organism` CLI argument is not found retrieved hits,
    the contig is flagged for removal
    '''


    axe = {}
    redo = {}

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
            
            sleep(2) # avoid hitting NCBI with too many requests

    return axe

def remove_bad_contigs(axe):
    '''Overwrites FASTAs minus any contigs flagged in query_genbank()'''

    for fasta in axe:
        towrite = []
        with open(fasta, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                
                name = rec.id

                if name not in axe[fasta]:
                    towrite.append(rec)
                else:
                    print "Removing contig {} from FASTA {}".format(name, fasta)
        
        with open(fasta, 'w') as j:

            SeqIO.write(towrite, j, 'fasta')

def filter_contigs(fasta_dir, fragment, gc_cutoff, organism, hits):
    '''All functions except the argument parser.'''
    
    checkable = select_checkables(fasta_dir, fragment, *gc_cutoff)
    
    axe = query_genbank(checkable, organism, hits)
    
    remove_bad_contigs(axe)


def main():

    args = arguments()
    
    filter_contigs(args.fasta_dir, args.fragment, args.gc_cutoff,
                   args.gc_cutoff, args.organism, args.hits)

if __name__ == '__main__':
    main()
