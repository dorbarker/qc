from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from time import sleep
from multiprocessing import Pool
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
    parser.add_argument('--min-contig', type = int, default = 500,
                        help = 'Minimum contig size for inclusion (default = 500 bp)')

    parser.add_argument('--fragment', type = int, default = 100,
                        help = 'Size (bp) of contig fragment used to query \
                                GenBank (default = 100)')
    
    parser.add_argument('--hits', type = int, default = 50,
                        help = 'Number top hits to consider (default = 50)')

    return parser.parse_args()

def gc_content(seq):
     '''Returns GC content of seq'''

     return 100.0 * sum(1.0 for x in seq if x in ('G','C')) / len(seq)

def select_checkables(fasta_dir, fragment, min_contig, low, high):
    '''Checks each contig in each FASTA for GC content.

    Out-of-range GC contents are later checked against GenBank.
    '''

    checkables = {}
    too_short = {}

    for fasta in (os.path.join(fasta_dir, x) for x in os.listdir(fasta_dir)):
        
        checkables[fasta] = {}

        too_short[fasta] = []
        
        with open(fasta, 'r') as f:
            
            records = 0.0
            total_len = 0.0
            bad_len = 0.0

            for rec in SeqIO.parse(f, 'fasta'):
                
                records += 1.0

                name = rec.id
                s    = str(rec.seq)
                gc   = gc_content(s)
                
                if len(s) < min_contig:
                    too_short[fasta].append(name)
                    continue

                total_len += len(s)

                if not low <= gc <= high:

                    checkables[fasta][name] = s[:fragment]
                    bad_len += len(s)

            # If 3/4 total contig length have the wrong GC, then the whole isolate
            # will be discarded at a later step and there's
            # no need to make all those BLAST searches
            
            try:
                if bad_len / total_len >= 0.75:
                    del checkables[fasta]
            except ZeroDivisionError:
                del checkables[fasta] # Wouldn't have been queried anyway
            
    remove_bad_contigsurn checkables

def format_queries(checkables, taxid, hits, organism, delay):

    delay = 0
    
    for fasta in checkables:
        for contig in checkables[fasta]:
            
            yield (checkables[fasta], taxid, hits, organism, delay)

            delay += 1 # so as not to spam GenBank, 
                       # queries are spaced ~1 second apart

def manage_queries(checkables, taxid, hits, organism, delay):

    for q in format_queries(checkables, taxid, hits, organism, delay):
        pass # multiprocessing.Pool and apply_async goes here 

def query(fasta_name, contig_name, fragment, taxid, hits, organism, delay):
   
    # staggers queries to not spam NCBI
    sleep(delay)
     
    txid = '{}[taxid]'.format(taxid) if taxid != None else '(none)'
    
    handle = NCBIWWW.qblast('blastn', 'nt', fragment, 
             megablast = True, entrez_query = txid, hitlist_size = hits,
             alignments = hits)

                            
    record = NCBIXML.read(handle)

    correct_organism = []

    for aln in record.alignments:
        for hsp in aln.hsps:
            correct_organism.append(organism.lower() in aln.title.lower())

    # Return tuple with names here because result order is not guaranteed
    # with apply_async and callback logging
    return fasta_name, contig_name, any(correct_organism)


def query_genbank(checkable, organism, taxid, hits):
    '''Any contigs identified in select_checkables() are queried against
    NCBI GenBank using a fragment from the beginning of the contig. 
    
    If the `organism` CLI argument is not found retrieved hits,
    the contig is flagged for removal
    '''

    axe = {}

    for fasta in checkable:
        axe[fasta] = []
        for contig in checkable[fasta]:
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

def filter_contigs(fasta_dir, fragment, min_contig, gc_cutoff,  organism, hits):
    '''All functions except the argument parser.'''
    
    checkable = select_checkables(fasta_dir, fragment, min_contig, *gc_cutoff)
    
    axe = query_genbank(checkable, organism, hits)
    
    remove_bad_contigs(axe)


def main():

    args = arguments()
    
    filter_contigs(args.fasta_dir, args.fragment, args.min_contig, 
            args.gc_cutoff, args.gc_cutoff, args.organism, args.hits)

if __name__ == '__main__':
    main()
