#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import contig_filter
import version_control
from multiprocessing import cpu_count

def arguments():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--fasta-parent-dir', required = True,
                        help = 'Parent directory of FASTA groups')

    parser.add_argument('--gc-cutoff', nargs = 2, type = float, metavar = ('LOW', 'HIGH'), required = True,
                        help = 'Lower and upper bounds for GC content')

    parser.add_argument('--size-cutoff', nargs = 2, type = int, metavar = ('LOW', 'HIGH'), required = True,
                        help = 'Lower and upper bounds of genome size (bp)')
    
    parser.add_argument('--db', required = True, help = 'Path to nt database') 
    
    parser.add_argument('--prefix', nargs = '+', help = 'Prefices to consider when removing duplicates')
    
    parser.add_argument('--organism', required = True, help = 'Name of the organism to which \
                                                              the contigs *should belong')

    parser.add_argument('--fragment', type = int, default = 100, help = 'Size (bp) of contig fragment \
                                                                        used to query GenBank \
                                                                        (default = 100)')
    parser.add_argument('--min-contig', type = int, default = 500, help = 'Minimum contig size for \
                                                                          inclusion (default = 500 bp)')

    parser.add_argument('--hits', type = int, default = 50, help = 'Number of BLAST hits to consider \
                                                                 during contig filtering (default = 50)')

    parser.add_argument('--cores', type = int, default = cpu_count(), help = 'CPUs to use')
   
    parser.add_argument('--skip-quast', action = 'store_true', help = 'Skip Quast and proceed to filtering')


    return parser.parse_args()

def run_quast(parent_dir, min_contig, ncpu):
    '''Runs Quast on each subdirectory (presumably full of fastas)
    of the parent directory.
    
    Results are saved in /parent_dir/subdir_quast_report/
    '''
    
    non_empty_dir = lambda d: os.path.isdir(d) and os.listdir(d)
    subdirs = [d for d in os.listdir(parent_dir) if non_empty_dir(d) and 'quast_report' not in d]
    
    quastpath = subprocess.check_output(['locate', '--regex', '/quast.py$']).strip() 
    
    for subdir in subdirs:

        s = os.path.join(parent_dir, subdir)
        
        fastas = [os.path.join(s, f) for f in os.listdir(s)]
        cmd = [quastpath, '--threads', str(ncpu), '--min-contig', min_contig, '--no-plots', '--no-html',
               '--output-dir', os.path.join(parent_dir, subdir + '_quast_report')]
        cmd.extend(fastas)
        
        process = subprocess.call(cmd)

def flag_duds(parent_dir):
    '''Extracts a list of assemblies with no contigs >= 500 bp.'''

    report_dirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if 'quast_report' in d]
    
    duds_dict = {}

    for rd in report_dirs:

        with open(os.path.join(rd, 'quast.log'), 'r') as f:

            duds = [line.strip().split()[2] for line in f if 'WARNING: Skipping' in line]
        
        duds_dict[rd] = duds
                    
    return duds_dict

def gc_filter_genomes(parent_dir, low, high):
    '''Flags genomes with out-of-range GC content'''

    report_dirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if 'quast_report' in d]
    
    bad_gc_dict = {}

    for rd in report_dirs:

        with open(os.path.join(rd, 'report.tsv'), 'r') as f:
            for line in f:
                if 'Assembly' in line:
                    names = line.strip().split('\t')[1:]
                elif 'GC (%)' in line:
                    gcs = map(float, line.strip().split('\t')[1:])
                else:
                    continue
            
            bad_gc_dict[rd] = [name for name, gc in zip(names, gcs) if not low <= gc <= high]

    return bad_gc_dict

def size_filter(parent_dir, low, high):
    '''Flags genomes with out-of-range size
    
    File size is heuristically used to estimate genome size
    '''

    fasta_dirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if 'quast_report' not in d]

    bad_size_dict = {}

    for fd in fasta_dirs:
        
        get_size = lambda d, f: os.path.getsize(os.path.join(d, f))
            
        # add '_quast_report' to simplify file writing later
        bad_size_dict[fd + '_quast_report'] = [os.path.splitext(x)[0] for x in os.listdir(fd)
                              if '.fasta' in x and not low <= get_size(fd, x) <= high]

    return bad_size_dict

def fix_name(name, prefices):
    '''Replace special characters with underscores'''
    
    def strip_prefix(s, p):
        return s[s.startswith(p) and len(p):]

    n = name
    
    for c in '-:;\', )(][\\':
        n = n.replace(c, '_')
    
    if prefices != None: 
        for p in prefices:

            n = strip_prefix(n, p)

    return n

def name_collisions(report_path, prefix, bad):
    '''Fixes naming conventions for a strain group.

    Collisions are handled by adding adding every duplicate
    except the best (measured by n50) to the 'bad' list and returning
    it.
    '''

    def match_names(fixed, original):
        '''Returns a dict of the fixed name as a key,
        and the original names that were associated with
        the fixed name.
        '''
     
        d = {}

        for key, value in zip(fixed, original):
            try:
                d[key].append(value)
            except KeyError:
                d[key] = [value]
        
        return d

    def name_n50_pairs(report_path, bad):
        '''Returns dictionary of strain names and their N50 values'''

        with open(report_path, 'r') as f:

            for line in f:
                if 'Assembly' in line:
                    names = line.strip().split('\t')[1:]
                elif 'N50' in line:
                    n50s = map(int, line.strip().split('\t')[1:])
                else:
                    continue

        return {name: n50 for name, n50 in zip(names, n50s) if name not in bad}
   
    def get_bad_twins(matched_names, names_n50s):
        '''For twinned names, give a list of the versions with
        less than the best N50
        '''

        bad = []
        for fixed in matched_names:

            if len(matched_names[fixed]) == 1:
                continue

            m = max( [(names_n50s[name], name) for name in matched_names[fixed]] )[1]
            
            bad.extend( [name for name in matched_names[fixed] if name != m] )

        return bad
    
    names_n50s = name_n50_pairs(report_path, bad)
    
    original_names = names_n50s.keys()
    
    fixed_names = [fix_name(name, prefix) for name in original_names]
    
    matched_names = match_names(fixed_names, original_names)
    
    bad_original_names = get_bad_twins(matched_names, names_n50s)

    return bad_original_names

def correct_names(parent_dir, prefix, bad):
    '''Handle running of name_collions on each fasta directory.

    Return a dictionary of the 'losing' names in the cases of 
    name collisions.
    '''

    report_dirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if 'quast_report' in d]
    
    bad_name_dict = {}

    for rd in report_dirs:

        report_path = os.path.join(rd, 'report.tsv') 
        bad_name_dict[rd] =  name_collisions(report_path, prefix, bad)
    
    return bad_name_dict

def merge_bad(*args):
    '''Flattens dictionary of bad genomes into a set'''

    bad = set()
     
    for method in args:
        for directory in method:
            for genome in method[directory]:
                bad.add(genome)
    return bad

def axe_bad(parent_dir, all_bad):
    
    fasta_dir = [d for d in os.listdir(parent_dir) if 'quast_report' not in d]

    for fasta_dir in fasta_dirs:

        full_fasta_dir_path = os.path.join(parent_dir, fasta_dir)

        for f in os.listdir(full_fasta_dir_path):

            if f[:f.endswith('.fasta') and f.rfind('.fasta')] in all_bad:
                
                try:
                    os.remove( os.path.join(full_fasta_dir_path, fasta))
                except OSError:
                    print("Encountered error removing {}".format(fasta), file = sys.stderr)

def write_reasons(*args):
    '''Write a list of bad genomes to the relevant
    *_quast_report directory
    '''
    bydir = {}
    for method in args:
        for directory in method:
            try:
                bydir[directory].extend(method[directory])
            except KeyError:
                bydir[directory] = method[directory]
    
    for d in bydir:
        with open(os.path.join(d, 'bad_genomes.txt'), 'w') as f:
            f.write( '\n'.join(set(bydir[d])) )

def gc_filter_contigs(parent_dir, fragment, min_contig, gc_cutoff, organism, hits, db):
    '''Check individual contigs for GC content.

    Contigs with out-of-range GC content are checked against GenBank.

    If none of the top hits are the correct organism, the contig is discarded
    '''
    
    fasta_dirs = [d for d in os.listdir(parent_dir) if 'quast_report' not in d]

    for fasta_dir in fasta_dirs:
        contig_filter.filter_contigs(fasta_dir, fragment, min_contig,
                                    gc_cutoff, organism, hits, db)
        
def main():
    
    args = arguments()
    
    version_control.commit(args.fasta_parent_dir, args)

    gc_filter_contigs(args.fasta_parent_dir, args.fragment, args.min_contig,
            args.gc_cutoff, args.organism, args.hits, args.db)

    if not args.skip_quast:
        run_quast(args.fasta_parent_dir, args.min_contig, args.cores)
    
    duds = flag_duds(args.fasta_parent_dir)
    
    bad_gc = gc_filter_genomes(args.fasta_parent_dir, *args.gc_cutoff)
    
    bad_size = size_filter(args.fasta_parent_dir, *args.size_cutoff)    
    
    bad_so_far = merge_bad(duds, bad_gc, bad_size)
    
    bad_names = correct_names(args.fasta_parent_dir, args.prefix, bad_so_far)    
    
    all_bad = merge_bad(duds, bad_gc, bad_size, bad_names)
    
    write_reasons(duds, bad_gc, bad_size, bad_names)
    
    axe_bad(parent_dir, all_bad)
     
    version_control.commit(args.fasta_parent_dir, args)

if __name__ == '__main__':
    main()

