from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import os
import argparse
from io import StringIO

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--fasta-dir', required = True,
                        help = 'Directory containing FASTA files')
    
    parser.add_argument('--db', required = True, help = 'Path to nt database')

    parser.add_argument('--organism', required = True,
                        help = 'Name of the organism the contigs \
                                *should* be from.')

    parser.add_argument('--gc-cutoff', nargs = 2, type = float,
                        required = True, metavar = ('LOW', 'HIGH'), 
                        help = 'GC contents outside this range \
                                will be queried against nt')

    parser.add_argument('--min-contig', type = int, default = 500,
                        help = 'Minimum contig size for inclusion (default = 500 bp)')

    parser.add_argument('--fragment', type = int, default = 100,
                        help = 'Size (bp) of contig fragment used to query \
                                GenBank (default = 100)')
    
    parser.add_argument('--cores', type = int, default = cpu_count,
                        help = 'Number of CPU cores to use')

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
            
    return checkables, too_short

def format_queries(checkables, organism, db_path):

    for fasta in checkables:
        for contig in checkables[fasta]:
            
            q = checkables[fasta][contig]
            yield (fasta, contig, q, db_path, organism)

def dispatch_queries(checkables, organism, db_path, cores):
   
    queries = format_queries(checkables, organism, db_path)
    
    p = Pool(processes = cores)

    results = [p.apply_async(query, args) for args in queries]

    return [m.get() for m in results] 

def query(fasta, contig, q, db_path, organism):
 
    search = blastn(query = StringIO(q), db = db_path, outfmt = 5)

    stdout, stderr = search()

    record = NCBIXML.read(stdout)

    correct = []

    for aln in record.alignments:
        for hsp in aln.hsps:
            correct.append(organism.lower() in aln.title.lower())

    return (fasta, contig, any(correct))

def remove_bad_contigs(fasta_dir, axe):
    '''Overwrites FASTAs minus any contigs flagged in query()'''

    fastas = (os.path.join(fasta_dir, x) for x in os.listdir(fasta_dir))
    
    for fasta in fastas:
        towrite = []
        with open(fasta, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                
                name = rec.id

                if name not in axe[fasta]:
                    towrite.append(rec)
                else:
                    msg = 'Removed contig {} from FASTA {}'.format(name, fasta)
                    print(msg)
        
        with open(fasta, 'w') as j:

            SeqIO.write(towrite, j, 'fasta')

def make_axe_dict(blast_results, too_short):

    axe = {}

    bad_blast = filter(lambda x: not x[2], blast_results)

    for i in bad_blast:
        try:
            axe[i[0]].add(i[1])
        except KeyError:
            axe[i[0]] = set([i[1]])
    
    for j in too_short:
        
        cur_set = set(too_short[j])

        try:
            axe[j] = axe[j] | cur_set

        except KeyError:
            axe[j] = cur_set 
    
    return axe 

def filter_contigs(fasta_dir, frag, min_contig, gc_cutoff, organism, db, cores):
    '''All functions except the argument parser.'''
    
    checkables, too_short  = select_checkables(fasta_dir, frag, 
                                               min_contig, *gc_cutoff)
   
    blast_results = dispatch_queries(checkables, organism, db, cores)

    axe = make_axe_dict(blast_results, too_short)

    remove_bad_contigs(axe)

def main():

    args = arguments()
    
    filter_contigs(args.fasta_dir, args.fragment, args.min_contig, 
                   args.gc_cutoff, args.gc_cutoff, args.organism, args.hits)


if __name__ == '__main__':
    main()
