 #    
 #    mapsplice_segments.py	
 #    MapSplice
 #
 #    Copyright (C) 2010 University of Kentucky and
 #                       Kai Wang
 #
 #    Authors: Kai Wang
 #
 #    This program is free software: you can redistribute it and/or modify
 #    it under the terms of the GNU General Public License as published by
 #    the Free Software Foundation, either version 3 of the License, or
 #    (at your option) any later version.
 #
 #    This program is distributed in the hope that it will be useful,
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #    GNU General Public License for more details.
 #
 #    You should have received a copy of the GNU General Public License
 #    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #
import sys
import getopt
import subprocess
import errno
import os, glob
import tempfile
import warnings
import shutil
import math  
import locale

from datetime import datetime, date, time

use_message = '''
Usage:
    mapsplice [Inputs] [options]

Inputs:
    -1/--end_1                     <string>    [ end 1 reads or single end reads          ]
    -2/--end_2                     <string>    [ end 2 reads                              ]
    -c/--chromosome-dir            <string>    [ input chromosomes directory              ]
    -x                             <string>    [ path and prefix of bowtie index          ]
Options:
    -p/--threads                   <int>       [ default: 1                               ]
    -o/--output                    <string>    [ default: ./mapsplice_out                 ]
    --bam                                      [ generate bam output                      ]
    --keep-tmp                                 [ keep the intermediate files              ]
    -l/--seglen                    <int>       [ default: 25                              ]
    --min-map-len                  <int>       [ default: 0                               ]
    --min-len                      <int>       [ default: 25                              ]
    -m/--splice-mis                <int>       [ default: 1                               ]
    --max-append-mis               <int>       [ default: 3                               ]
    -k/max-hits                    <int>       [ default: 4                               ]
    -i/--min-intron                <int>       [ default: 50                              ]
    -I/--max-intron                <int>       [ default: 200,000                         ]
    --qual-scale                   <string>    [ phred64(default) or phred33 or solexa64  ]
    --del                          <int>       [ default: 6                               ]
    --ins                          <int>       [ default: 6                               ]
    --non-canonical                            [ output noncanonical junction             ]
    --fusion | --fusion-non-canonical          [ output fusion junction                   ]
    -h/--help                                  [ print the usage message                  ]
    -v/--version                               [ print the version of MapSplice           ]

'''

ver_message = '''
MapSplice v2.0.1.6
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg
        
class Ver(Exception):
    def __init__(self, msg):
        self.msg = msg
	
output_dir = "mapsplice_out/"
logging_dir = output_dir + "logs/"
canon_in_intron_dir = output_dir + "canonical/"
canon_exceed_intron_dir = output_dir + "canonical_exceed/"
noncanon_in_intron_dir = output_dir + "noncanonical/"
noncanon_exceed_intron_dir = output_dir + "noncanonical_exceed/"
fusion_dir = output_dir + "fusion/"

tophat_dir = output_dir + "tophat/"

temp_dir = output_dir + "tmp/"
original_dir = temp_dir + "original/"
filteredbest_dir = temp_dir + "best/"
comparison_dir = temp_dir + "comparison/"

synthetic_dir = temp_dir + "synthetic/"

pairend_dir = temp_dir + "pairend/"

remap_dir = temp_dir + "remap/"

remap_regions_dir = temp_dir + "remap_regions/"

hmer_dir = temp_dir + "single_anchored_middle/"
hole_dir = temp_dir + "double_anchored/"
head_dir = temp_dir + "single_anchored_head/"
tail_dir = temp_dir + "single_anchored_tail/"
bwtout_dir = temp_dir + "bwtout/"
sam_dir = temp_dir + "sam/"
fusion_dir = temp_dir + "fusion/"
cluster_dir = temp_dir + "cluster/"
cluster_result_dir = cluster_dir + "result/"
cluster_data_dir = cluster_dir + "data/"
cluster_data_parsedPER_dir = cluster_data_dir + "parsedPER/"

fusion_data_dir = fusion_dir + "data/"
fusion_data_single_dir = fusion_data_dir + "single/"
fusion_data_PER_dir = fusion_data_dir + "PER/"
fusion_result_dir = fusion_dir + "result/"
fusion_result_PER_prob_dir = fusion_result_dir + "PER_prob/"
fusion_result_fusionRead_dir = fusion_result_dir + "fusionRead/"
fusion_result_junction_support_dir = fusion_result_dir + "junction_support/"

formated_chrom_dir = temp_dir + "formated_chrom/"

formated_reads_dir = temp_dir + "formated_reads/"

filter_repeats_dir = temp_dir + "filter_repeats/"

bin_dir = sys.path[0] + "/"

rerun_all = 1;
DEBUG = 0;
max_chromosome = 1000;

fail_str = "\t[FAILED]\n"

class Params:
    
    def __init__(self): 
	
	self.min_anchor_length = 8
	self.seed_length = 10
        self.splice_mismatches = 1
	self.segment_mismatches = 1
	self.FASTA_file_extension = "fa"
	self.read_file_suffix = "txt"
        self.min_intron_length = 10
        self.max_intron_length = 200000
        self.island_extension = 0
        self.read_width = 0
        self.rank = 0.0
        self.flank_case = 5
        self.islands_file = ""
        self.read_files_dir = ""
        self.chromosome_files_dir = ""
        self.all_chromosomes_file = ""
	self.repeat_regioins = ""
	self.gene_regions = ""
        
        self.bwt_idx_prefix = ""
        self.bowtie_threads = 1
	self.max_hits = 4
        self.threshold = 1
        self.boundary = 50

        self.num_anchor = 1
        
        self.unmapped_reads = ""
        self.sam_formatted = ""
	self.bam_formatted = ""
	
        self.sam_formatted_25 = ""
	self.bwt_map_25 = ""
	
        self.pileup_file = ""
	
	self.synthetic_mappedreads = ""
	
	self.tophat_mappedreads = ""
	
	self.pairend = ""
        
        self.gamma = 0.1
        self.delta = 0.1
	
	#self.num_seg = 4
	self.seg_len = 25

	self.fix_hole_file = "" 
	
	self.format_flag = ""
	
	self.chrom_size_file = ""
	
	self.bam_file = ""
	
	self.extend_bits = 3;
	
	self.total_fusion_mismatch = 2;
	
	self.total_mismatch = 2;
	
	self.append_mismatch = 2;
	
	self.remap_mismatch = 2;
	
	self.skip_bwt = 0;
	
	self.prefix_match = 1;
	
	self.fullrunning = 0;
	
	self.collect_stat = 0;
	
	self.rm_temp = 1;
	
	self.format_reads = 0;
	
	self.format_chromos = 0;
	
	self.do_fusion = 0;
	
	self.do_cluster = 0;
	
	self.search_whole_chromo = 0;
	
	self.map_segment_directly = 0;
	
	self.run_mapper = 0;
	
	self.max_insert = 3;
	
	self.min_missed_seg = 0;
	
	self.min_map_len = 0;
	
	self.fusion_flank_case = 5;
	
	self.do_annot_gene = 0;
	
	self.annot_gene_file = ""
	
	self.filter_fusion_by_repeat = 0;
	
	self.chromosome_blat_idx = "";
	
	#global rerun_all
	
	#global DEBUG
    
    def parse_cfgfile(self, cfg_file):
	fh = open(cfg_file,"r")
	igot = fh.readlines()
	
	for line in igot:
	    if line[0] == '\n' or line[0] == '#':
		comments = line;
	    else:
		line = line.rstrip();
		command = line.split(' = ');
		
		if command[0] == 'reads_file':
		    self.read_files_dir = command[1];
		if command[0] == 'sam_file':
		    self.sam_formatted = command[1];
		if command[0] == 'chromosome_files_dir':
		    self.chromosome_files_dir = command[1];
		if command[0] == 'Bowtieidx':
		    self.bwt_idx_prefix = command[1];
		if command[0] == 'avoid_regions':
		    self.repeat_regioins = command[1];
		if command[0] == 'interested_regions':
		    self.gene_regions = command[1];
		if command[0] == 'output_dir':
		    global output_dir
		    global logging_dir
		    global canon_in_intron_dir
		    global canon_exceed_intron_dir
		    global noncanon_in_intron_dir
		    global noncanon_exceed_intron_dir
		    global fusion_dir
		    global temp_dir
		    global synthetic_dir
		    global pairend_dir
		    global original_dir
		    global filteredbest_dir
		    global comparison_dir
		    global tophat_dir
		    global remap_dir
		    global remap_regions_dir
		    global hmer_dir
		    global hole_dir
		    global head_dir
		    global tail_dir
		    global bwtout_dir
		    global sam_dir
		    global formated_chrom_dir
		    global formated_reads_dir
		    global fusion_dir
		    global fusion_data_dir
		    global fusion_data_single_dir
		    global fusion_data_PER_dir
		    global fusion_result_dir
		    global fusion_result_PER_prob_dir
		    global fusion_result_junction_support_dir
		    global cluster_dir
		    global cluster_result_dir
		    global cluster_data_dir
		    global cluster_data_parsedPER_dir
		    global fusion_result_fusionRead_dir
		   
		    output_dir = command[1] + "/"
		    logging_dir = output_dir + "logs/"
		    canon_in_intron_dir = output_dir + "canonical/"
		    canon_exceed_intron_dir = output_dir + "canonical_exceed/"
		    noncanon_in_intron_dir = output_dir + "noncanonical/"
		    noncanon_exceed_intron_dir = output_dir + "noncanonical_exceed/"
		    fusion_dir = output_dir + "fusion/"
		    temp_dir = output_dir + "tmp/"
		    
		    synthetic_dir = temp_dir + "synthetic/"
		    
		    pairend_dir = temp_dir + "pairend/"
		    
		    original_dir = temp_dir + "original/"
		    filteredbest_dir = temp_dir + "best/"
		    comparison_dir = temp_dir + "comparison/"
		    
		    tophat_dir = temp_dir + "tophat/"
		    
		    remap_dir = temp_dir + "remap/"
		    
		    remap_regions_dir = temp_dir + "remap_regions/"
		    
		    hmer_dir = temp_dir + "single_anchored_middle/"
		    hole_dir = temp_dir + "double_anchored/"
		    head_dir = temp_dir + "single_anchored_head/"
		    tail_dir = temp_dir + "single_anchored_tail/"
		    bwtout_dir = temp_dir + "bwtout/"
		    
		    sam_dir = temp_dir + "sam/"
		    
		    fusion_dir = temp_dir + "fusion/"
		    
		    cluster_dir = temp_dir + "cluster/"
		    
		    cluster_result_dir = cluster_dir + "result/"
		    cluster_data_dir = cluster_dir + "data/"
		    cluster_data_parsedPER_dir = cluster_data_dir + "parsedPER/"
		    
		    fusion_data_dir = fusion_dir + "data/"
		    fusion_data_single_dir = fusion_data_dir + "single/"
		    fusion_data_PER_dir = fusion_data_dir + "PER/"
		    fusion_result_dir = fusion_dir + "result/"
		    fusion_result_PER_prob_dir = fusion_result_dir + "PER_prob/"
		    fusion_result_junction_support_dir = fusion_result_dir + "junction_support/"
				    
		    formated_chrom_dir = temp_dir + "formated_chrom/"
    
		    formated_reads_dir = temp_dir + "formated_reads/"
		    
		    fusion_result_fusionRead_dir = fusion_result_dir + "fusionRead/"
		   
		if command[0] == 'reads_format':
		    if command[1] == 'FASTA':
			self.format_flag = '-f';
		    elif command[1] == 'FASTQ':
			self.format_flag = '-q';		
		elif command[0] == 'segment_mismatches':
		    self.segment_mismatches = int(command[1]);
		elif command[0] == 'segment_length':
		    self.seg_len = int(command[1]);
		elif command[0] == 'read_length':
		    self.read_width = int(command[1]);
		elif command[0] == 'paired_end':
		    if command[1] == 'yes':
			self.pairend = '1';
		    elif command[1] == 'no':
			self.pairend = '';
		elif command[0] == 'junction_type':
		    if command[1] == 'non-canonical':
			self.flank_case = 0;
		    elif command[1] == 'semi-canonical':
			self.flank_case = 1;
		    elif command[1] == 'canonical':
			self.flank_case = 5;
		elif command[0] == 'fusion_junction_type':
		    if command[1] == 'non-canonical':
			self.fusion_flank_case = 0;
		    elif command[1] == 'semi-canonical':
			self.fusion_flank_case = 1;
		    elif command[1] == 'canonical':
			self.fusion_flank_case = 5;		
		elif command[0] == 'full_running':
		    if command[1] == 'yes':
			self.fullrunning = 1;
		    elif command[1] == 'no':
			self.fullrunning = 0;
		elif command[0] == 'anchor_length':
		    self.min_anchor_length = int(command[1]);
		elif command[0] == 'remove_temp_files':
		    if command[1] == 'yes':
			self.rm_temp = 1;
		    elif command[1] == 'no':
			self.rm_temp = 0;
		elif command[0] == 'remap_mismatches':
		    self.remap_mismatch = int(command[1]);
		elif command[0] == 'splice_mismatches':
		    self.splice_mismatches = int(command[1]);
		elif command[0] == 'min_intron_length':
		    self.min_intron_length = int(command[1]);
		elif command[0] == 'max_intron_length':
		    self.max_intron_length = int(command[1]);
		elif command[0] == 'threads':
		    self.bowtie_threads = int(command[1]);
		elif command[0] == 'search_whole_chromosome':
		    if command[1] == 'yes':
			self.search_whole_chromo = 1;
		    elif command[1] == 'no':
			self.search_whole_chromo = 0;
		elif command[0] == 'map_segment_directly':
		    if command[1] == 'yes':
			self.map_segment_directly = 1;
		    elif command[1] == 'no':
			self.map_segment_directly = 0;
		elif command[0] == 'run_MapPER':
		    if command[1] == 'yes':
			self.run_mapper = 1;
		    elif command[1] == 'no':
			self.run_mapper = 0;
		elif command[0] == 'do_fusion':
		    if command[1] == 'yes':
			self.do_fusion = 1;
		    elif command[1] == 'no':
			self.do_fusion = 0;
		elif command[0] == 'do_cluster':
		    if command[1] == 'yes':
			self.do_cluster = 1;
		    elif command[1] == 'no':
			self.do_cluster = 0;
		elif command[0] == 'max_insert':	
		    self.max_insert = int(command[1]);
		elif command[0] == 'BAM':	
		    self.bam_file = command[1];
		elif command[0] == 'min_missed_seg':	
		    self.min_missed_seg = int(command[1]);
		elif command[0] == 'min_map_len':	
		    self.min_map_len = int(command[1]);
		elif command[0] == 'max_hits':	
		    self.max_hits = int(command[1]);
		elif command[0] == 'do_annot_gene':	    
		    self.do_annot_gene = 1;
		    self.annot_gene_file = (command[1]);
		elif command[0] == 'filter_fusion_by_repeat':   
		    self.chromosome_blat_idx = (command[1]);
		    self.filter_fusion_by_repeat = 1;
		
        
def get_version():
    return "2.0"

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def prepare_output_dir():
    #global output_dir
    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)
	
    #if os.path.exists(canon_in_intron_dir):
        #pass
    #else:        
        #os.mkdir(canon_in_intron_dir)

    #if os.path.exists(canon_exceed_intron_dir):
        #pass
    #else:        
        #os.mkdir(canon_exceed_intron_dir)
	
    #if os.path.exists(noncanon_in_intron_dir):
        #pass
    #else:        
        #os.mkdir(noncanon_in_intron_dir)
	
    #if os.path.exists(noncanon_exceed_intron_dir):
        #pass
    #else:        
        #os.mkdir(noncanon_exceed_intron_dir)	
	
    #if os.path.exists(fusion_dir):
        #pass
    #else:        
        #os.mkdir(fusion_dir)	
	
    if os.path.exists(temp_dir):
        pass
    else:        
        os.mkdir(temp_dir)
	
    #if os.path.exists(synthetic_dir):
        #pass
    #else:        
        #os.mkdir(synthetic_dir)
	
    if os.path.exists(original_dir):
        pass
    else:        
        os.mkdir(original_dir)
	
    if os.path.exists(filteredbest_dir):
        pass
    else:        
        os.mkdir(filteredbest_dir)
	
    #if os.path.exists(comparison_dir):
        #pass
    #else:        
        #os.mkdir(comparison_dir)
	
    #if os.path.exists(tophat_dir):
        #pass
    #else:        
        #os.mkdir(tophat_dir)
	
    if os.path.exists(remap_dir):
	pass
    else:        
	os.mkdir(remap_dir)
	
    #if os.path.exists(remap_regions_dir):
	#pass
    #else:        
	#os.mkdir(remap_regions_dir)
	
    #if os.path.exists(hmer_dir):
	#pass
    #else:        
	#os.mkdir(hmer_dir)
	
    #if os.path.exists(hole_dir):
	#pass
    #else:        
	#os.mkdir(hole_dir)
	
    #if os.path.exists(head_dir):
	#pass
    #else:        
	#os.mkdir(head_dir)
	
    #if os.path.exists(tail_dir):
	#pass
    #else:        
	#os.mkdir(tail_dir)
	
    if os.path.exists(fusion_dir):
	pass
    else:        
	os.mkdir(fusion_dir)
    
    if os.path.exists(cluster_dir):
	pass
    else:        
	os.mkdir(cluster_dir)
	
    if os.path.exists(cluster_result_dir):
	pass
    else:        
	os.mkdir(cluster_result_dir)
	
    if os.path.exists(cluster_data_dir):
	pass
    else:        
	os.mkdir(cluster_data_dir)	

    if os.path.exists(cluster_data_parsedPER_dir):
	pass
    else:        
	os.mkdir(cluster_data_parsedPER_dir)	
	

    #if os.path.exists(bwtout_dir):
	#pass
    #else:        
	#os.mkdir(bwtout_dir)
	
    #if os.path.exists(sam_dir):
	#pass
    #else:        
	#os.mkdir(sam_dir)
	
    #if os.path.exists(formated_chrom_dir):
	#pass
    #else:        
	#os.mkdir(formated_chrom_dir)
	
    #if os.path.exists(formated_reads_dir):
	#pass
    #else:        
	#os.mkdir(formated_reads_dir)
    
    if os.path.exists(fusion_data_dir):
	pass
    else:        
	os.mkdir(fusion_data_dir)
	
    if os.path.exists(fusion_data_single_dir):
	pass
    else:        
	os.mkdir(fusion_data_single_dir)
	
    if os.path.exists(fusion_data_PER_dir):
	pass
    else:        
	os.mkdir(fusion_data_PER_dir)
	
    if os.path.exists(fusion_result_dir):
	pass
    else:        
	os.mkdir(fusion_result_dir)
	
    if os.path.exists(fusion_result_fusionRead_dir):
	pass
    else:        
	os.mkdir(fusion_result_fusionRead_dir)

    if os.path.exists(fusion_result_PER_prob_dir):
	pass
    else:        
	os.mkdir(fusion_result_PER_prob_dir)
	
    if os.path.exists(fusion_result_junction_support_dir):
	pass
    else:        
	os.mkdir(fusion_result_junction_support_dir)

    if os.path.exists(filter_repeats_dir):
	pass
    else:
	os.mkdir(filter_repeats_dir)
     
    #if os.path.exists(pairend_dir):
        #pass
    #else:        
        #os.mkdir(pairend_dir)

def read_sam(path_name):
    
    sequences = []
    flag = 0;
    if os.path.isdir(path_name) == False:
        sequences = sequences + [path_name]
        return sequences
    for infile in glob.glob( os.path.join(path_name, '*.sam') ):
        if flag == 1 :
            sequences = sequences + [infile]
        else:
            sequences = sequences + [infile]
            flag = 1;
    return sequences

def read_chroms(path_name):
    
    sequences = []
    flag = 0;
    if os.path.isdir(path_name) == False:
        sequences = sequences + [path_name]
        return sequences
    for infile in glob.glob( os.path.join(path_name, '*.fa') ):
        if flag == 1 :
            sequences = sequences + [infile]
        else:
            sequences = sequences + [infile]
            flag = 1;
    return sequences

def read_dir_by_suffix(path_name, suffix):
    
    suffix = '*.' + suffix;
    
    sequences = []
    flag = 0;
    if os.path.isdir(path_name) == False:
	sequences = sequences + [path_name]
        return sequences
    for infile in glob.glob( os.path.join(path_name, suffix) ):
        if flag == 1 :
            sequences = sequences + [infile]
        else:
            sequences = sequences + [infile]
            flag = 1;
    return sequences

def read_sequence(path_name):
    
    sequences = ''
    flag = 0;
    if os.path.isdir(path_name) == False:
        return path_name
    for infile in glob.glob( os.path.join(path_name, '*.fa') ):
        if flag == 1 :
            sequences += ',' + infile
        else:
            sequences = infile
            flag = 1;
    return sequences

def read_sequence_by_suffix(path_name, suffix):
    
    suffix = '*.' + suffix;
    
    sequences = ''
    flag = 0;
    
    chromosome_count = 0;
    
    if os.path.isdir(path_name) == False:
        return path_name
    for infile in glob.glob( os.path.join(path_name, suffix) ):
        if flag == 1 :
            sequences += ',' + infile
	    
	    chromosome_count = chromosome_count + 1;
        else:
            sequences = infile
            flag = 1;
    
    if chromosome_count > max_chromosome:
	
	all_chromos_path = read_dir_by_suffix(chromosome_files_dir, suffix)
	cat_files(all_chromos_path, formated_chrom_dir + "combined_chromosomes.fa");
	
	return formated_chrom_dir + "combined_chromosomes.fa";
	#max_chromosome
	#print >> sys.stderr, "Error: too many chromosomes %s, try combine them to build index" % str(chromosome_count)
	
    return sequences

def check_islands(islands_file):
    print >> sys.stderr, "[%s] Checking for islands file" % right_now()

    if os.path.exists(islands_file):
        return islands_file
    else:
        print >> sys.stderr, "Error: Could not find islands file %s" % (islands_file)
        exit(1)

def check_mapsplice():
    print >> sys.stderr, "[%s] Checking for mapsplice exec file" % right_now()

    if os.path.exists("mapsplice") and \
       os.path.exists("mapsplice_128"):
        return
    else:
        os.system("make")
        if os.path.exists("mapsplice") and \
           os.path.exists("mapsplice_128"):
            return
        else:
            print >> sys.stderr, "Error: Could not find mapsplice exec files"
            exit(1)

def check_reads_files(reads_files_or_dir):
    print >> sys.stderr, "[%s] Checking for reads files or directory" + reads_files_or_dir % right_now()

    if os.path.exists(reads_files_or_dir):
        return reads_files_or_dir
    else:
        print >> sys.stderr, "Error: Could not find reads files or directory " + reads_files_or_dir
        exit(1)

def check_file_existence(reads_files_or_dir):
    print >> sys.stderr, "[%s] Checking for files or directory" % right_now()

    if os.path.exists(reads_files_or_dir):
        return reads_files_or_dir
    else:
        print >> sys.stderr, "Error: Could not find files or directory " + reads_files_or_dir
        exit(1)
	
def check_chromo_files(chromosomes_files_or_dir):
    print >> sys.stderr, "[%s] Checking for chromosomes files or directory" % right_now()

    if os.path.exists(chromosomes_files_or_dir):

	mapsplice_log = open(logging_dir + "checking_chromosome_file.log", "w")
	
	all_chromos_path = read_dir_by_suffix(chromosomes_files_or_dir, "fa")
	
	merge_sam_cmd = [bin_dir + "check_input_files"]  
    
	for sam_file in all_chromos_path:
	    merge_sam_cmd = merge_sam_cmd + [sam_file]
	    
	try:    
	    retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
	   
	    if retcode == 2 or retcode == 3:
		return retcode
	    
	    if retcode == 0:
		print >> sys.stderr, "[%s] Checking for chromosomes files or directory passed" % right_now()	    
		
	except OSError, o:
	    if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
		print >> sys.stderr, fail_str, "Error: check_input_files not found on this system"
	    exit(1)
	
        return 0
    
    else:
        print >> sys.stderr, "Error: Could not find chromosomes files or directory " + chromosomes_files_or_dir
        exit(1)
    
def check_allchromo_file(all_chromosomes_file):
    print >> sys.stderr, "[%s] Checking for all chromosomes file" % right_now()

    if os.path.exists(all_chromosomes_file):
        return all_chromosomes_file
    else:
        print >> sys.stderr, "Error: Could not find all chromosomes file " + chromosomes_files_or_dir
        exit(1)

def call_mapsplice(islands_gff, 
                   seed_length,
                   read_width,
                   min_anchor_len,
                   splice_mismatches,
                   min_intron_length,
                   max_intron_length,
                   islands_extension,
                   flank_case,
                   rank,
                   FASTA_file_ext,
                   reads_files_dir,
                   chromosomes_file_dir,
                   num_anchor):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning spliced reads" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(logging_dir + "mapsplice.log", "w")
    
    juncfile = output_dir + "junctions.txt"
    if read_width > 64:
        splice_cmd = [bin_dir + "mapsplice_segments",
                      "-v", output_dir,          # Output dir
                      "-n", str(min_anchor_len), # Anchor length
                      "-m", str(splice_mismatches), # Mismatches allowed in extension
                      "-x", str(max_intron_length), # Maxmimum intron length
                      "-i", str(min_intron_length), # Minimum intron length
                      "-h", str(seed_length), # Seed size for reads
                      "-w", str(read_width), # read width for reads
                      "-p", str(islands_extension), # islands extension
                      "-b", str(350000000), # block size for reading chromosome
                      "-s", FASTA_file_ext, # FASTA file extension
                      "-t", islands_gff, # island location
                      "-u", reads_files_dir, # reads file or directory
                      "-c", chromosomes_file_dir, # chromosomes file or directory
                      "-y", str(1), # nothing important
                      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
                      "-R", str(rank)   #only output rank >= rank
                      ]
    else:
        splice_cmd = [bin_dir + "mapsplice",
                      "-v", output_dir,          # Output dir
                      "-n", str(min_anchor_len), # Anchor length
                      "-m", str(splice_mismatches), # Mismatches allowed in extension
                      "-x", str(max_intron_length), # Maxmimum intron length
                      "-i", str(min_intron_length), # Minimum intron length
                      "-h", str(seed_length), # Seed size for reads
                      "-w", str(read_width), # read width for reads
                      "-p", str(islands_extension), # islands extension
                      "-b", str(350000000), # block size for reading chromosome
                      "-s", FASTA_file_ext, # FASTA file extension
                      "-t", islands_gff, # island location
                      "-u", reads_files_dir, # reads file or directory
                      "-c", chromosomes_file_dir, # chromosomes file or directory
                      "-y", str(1), # nothing important
                      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
                      "-A", str(num_anchor), #number of anchors
                      "-R", str(rank)   #only output rank >= rank
                      ]
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
            exit(1)
        #else:
        #    sortjunc_cmd = [bin_dir + "sortjunc/sortjunc", output_dir + "/" + juncfile]
        #    subprocess.call(sortjunc_cmd)                
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MapSplice not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    
    return juncfile

def call_mapsplice_segment2(islands_gff, 
                   seed_length,
                   read_width,
                   min_anchor_len,
                   splice_mismatches,
                   min_intron_length,
                   max_intron_length,
                   islands_extension,
                   flank_case,
                   FASTA_file_ext,
                   chromosomes_file_dir,
                   num_anchor,
		   num_seg,
		   seg_len,
		   bwtout_25,
		   cur_ouput_dir,
		   extend_bits,
		   tmp_dir):
    
    print >> sys.stderr, "[%s] Aligning spliced reads" % right_now()
    
    splice_cmd = ""
 
    test_file =  hmer_dir + "chrX.txt.fixed"
    if os.path.exists(test_file) and \
	rerun_all == 0:
	return test_file
    if read_width > 36:
	splice_cmd = [bin_dir + "mapsplice_segments",
		      "-v", cur_ouput_dir,          # Output dir
		      "-n", str(min_anchor_len), # Anchor length
		      "-m", str(splice_mismatches), # Mismatches allowed in extension
		      "-x", str(max_intron_length), # Maxmimum intron length
		      "-i", str(min_intron_length), # Minimum intron length
		      "-h", str(seed_length), # Seed size for reads
		      "-w", str(read_width), # read width for reads
		      "-p", str(islands_extension), # islands extension
		      "-b", str(350000000), # block size for reading chromosome
		      #"-s", FASTA_file_ext, # FASTA file extension
		      "-t", islands_gff, # island location
		      "-c", chromosomes_file_dir, # chromosomes file or directory
		      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
		      #"-A", str(num_anchor), #number of anchors
		      "-bwtout", bwtout_25, #bwt segment output
		      "-G", str(num_seg), #number of segments
		      "-L", str(seg_len), #segment length
		      #"-E", str(extend_bits), #extend bits when fix hole and fusion
		      "-tmp", tmp_dir #tmp dir
		      ]

    else:
        splice_cmd = [bin_dir + "mapsplice",
                      "-v", cur_ouput_dir,          # Output dir
                      "-n", str(min_anchor_len), # Anchor length
                      "-m", str(splice_mismatches), # Mismatches allowed in extension
                      "-x", str(max_intron_length), # Maxmimum intron length
                      "-i", str(min_intron_length), # Minimum intron length
                      "-h", str(seed_length), # Seed size for reads
                      "-w", str(read_width), # read width for reads
                      "-p", str(islands_extension), # islands extension
                      "-b", str(350000000), # block size for reading chromosome
                      "-s", FASTA_file_ext, # FASTA file extension
                      "-t", islands_gff, # island location
                      "-u", reads_files_dir, # reads file or directory
                      "-c", chromosomes_file_dir, # chromosomes file or directory
                      "-y", str(1), # nothing important
                      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
                      "-A", str(num_anchor), #number of anchors
                      "-R", str(rank)   #only output rank >= rank
                      ]
	
    print >> sys.stderr, "[%s] " % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
            exit(1)
        #else:
        #    sortjunc_cmd = [bin_dir + "sortjunc/sortjunc", output_dir + "/" + juncfile]
        #    subprocess.call(sortjunc_cmd)                
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MapSplice not found on this system"
        exit(1)

def call_mapsplice_segment(islands_gff, 
                   seed_length,
                   read_width,
                   min_anchor_len,
                   splice_mismatches,
                   min_intron_length,
                   max_intron_length,
                   islands_extension,
                   flank_case,
                   rank,
                   FASTA_file_ext,
                   reads_files_dir,
                   chromosomes_file_dir,
                   num_anchor,
		   num_seg,
		   seg_len,
		   bwtout_25,
		   fix_hole_file,
		   cur_ouput_dir,
		   fqreads,
		   extend_bits,
		   total_mismatch,
		   total_fusion_mismatch,
		   append_mismatch,
		   prefix_match,
		   threads,
		   search_whole_chromo,
		   max_insert):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning spliced reads" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(logging_dir + "mapsplice.log", "w")
    #stdout=mapsplice_log
    
    juncfile = original_dir + "canonical_junctions.txt"
    
    fqreads = "";
    
    test_file =  hmer_dir + "chrX.txt.fixed"
    if os.path.exists(test_file) and \
	   rerun_all == 0:
	return test_file
    if read_width > 0:
	
	if fqreads == "":
	    splice_cmd = [bin_dir + "mapsplice_segments",
			  "-v", cur_ouput_dir,          # Output dir
			  "-n", str(min_anchor_len), # Anchor length
			  "-m", str(splice_mismatches), # Mismatches allowed in extension
			  "-x", str(max_intron_length), # Maxmimum intron length
			  "-i", str(min_intron_length), # Minimum intron length
			  "-h", str(seed_length), # Seed size for reads
			  "-w", str(read_width), # read width for reads
			  "-p", str(islands_extension), # islands extension
			  "-b", str(350000000), # block size for reading chromosome
			  "-s", FASTA_file_ext, # FASTA file extension
			  #"-t", islands_gff, # island location
			  #"-u", reads_files_dir, # reads file or directory
			  "-c", chromosomes_file_dir, # chromosomes file or directory
			  "-y", str(prefix_match), # nothing important
			  "-f", str(flank_case), #if is 1, only output flank string that is not case 0
			  "-A", str(num_anchor), #number of anchors
			  #"-R", str(rank),   #only output rank >= rank
			  "-bwtout", bwtout_25, #bwt segment output
			  "-G", str(num_seg), #number of segments
			  "-L", str(seg_len), #segment length
			  #"-H", fix_hole_file, #fix hole file
			  "-E", str(extend_bits), #extend bits when fix hole and fusion
			  #"-M", str(total_mismatch), #total mismatch on splice reads
			  #"-FM", str(total_fusion_mismatch), #total fusion mismatch
			  #"-P", str(append_mismatch), #append mismatch
			  "-threads", str(threads), #number of threads used
			  "-max_insertion", str(max_insert), #maximal insert
			  "-tmp", temp_dir #tmp dir
			  ]
	else:
	    splice_cmd = [bin_dir + "mapsplice_segments",
			  "-v", cur_ouput_dir,          # Output dir
			  "-n", str(min_anchor_len), # Anchor length
			  "-m", str(splice_mismatches), # Mismatches allowed in extension
			  "-x", str(max_intron_length), # Maxmimum intron length
			  "-i", str(min_intron_length), # Minimum intron length
			  "-h", str(seed_length), # Seed size for reads
			  "-w", str(read_width), # read width for reads
			  "-p", str(islands_extension), # islands extension
			  "-b", str(350000000), # block size for reading chromosome
			  "-s", FASTA_file_ext, # FASTA file extension
			  #"-t", islands_gff, # island location
			  #"-u", reads_files_dir, # reads file or directory
			  "-c", chromosomes_file_dir, # chromosomes file or directory
			  "-y", str(1), # nothing important
			  "-f", str(flank_case), #if is 1, only output flank string that is not case 0
			  "-A", str(num_anchor), #number of anchors
			  #"-R", str(rank),   #only output rank >= rank
			  "-bwtout", bwtout_25, #bwt segment output
			  "-G", str(num_seg), #number of segments
			  "-L", str(seg_len), #segment length
			  #"-H", fix_hole_file, #fix hole file
			  "-E", str(extend_bits), #extend bits when fix hole and fusion
			  #"-M", str(total_mismatch), #total mismatch on splice reads
			  #"-FM", str(total_fusion_mismatch), #total fusion mismatch
			  #"-P", str(append_mismatch), #append mismatch
			  "-threads", str(threads), #number of threads used
			  "-tmp", temp_dir #tmp dir
			  #"-FQ", fqreads #fastq reads
			  ]
    else:
        splice_cmd = [bin_dir + "mapsplice",
                      "-v", cur_ouput_dir,          # Output dir
                      "-n", str(min_anchor_len), # Anchor length
                      "-m", str(splice_mismatches), # Mismatches allowed in extension
                      "-x", str(max_intron_length), # Maxmimum intron length
                      "-i", str(min_intron_length), # Minimum intron length
                      "-h", str(seed_length), # Seed size for reads
                      "-w", str(read_width), # read width for reads
                      "-p", str(islands_extension), # islands extension
                      "-b", str(350000000), # block size for reading chromosome
                      "-s", FASTA_file_ext, # FASTA file extension
                      "-t", islands_gff, # island location
                      "-u", reads_files_dir, # reads file or directory
                      "-c", chromosomes_file_dir, # chromosomes file or directory
                      "-y", str(1), # nothing important
                      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
                      "-A", str(num_anchor), #number of anchors
                      "-R", str(rank)   #only output rank >= rank
                      ]
	
    if search_whole_chromo == 0:
	splice_cmd = splice_cmd + ["-t"] + [islands_gff];
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] separate unique and multiple mapped reads" % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MapSplice not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    
    return juncfile





def call_mapsplice_multithreads(mis_match, 
                                read_format,
                                max_hits,
                                threads,
                                seg_len,
                                chromosome_size,
                                split_index_path,
                                ref_seq_path,
                                flankcase,
                                max_insertion,
                                max_deletion,
                                min_missed_seg,
                                min_map_len,
                                mapsplice_out,
                                bowtie_index,
                                input_reads_1,
                                input_reads_2,
                                unmapped_reads,
                                bowtie_output,
                                juncdb_index,
                                no_original_junc,
                                do_fusion_and_is_paired,
                                max_intron_length,
                                append_mismatch,
                                min_read_len,
                                min_intron,
                                qual_scale,
                                log_file):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Run MapSplice multi-thread" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(log_file, "w")
    
    if os.path.exists(mapsplice_out + ".1") and \
	   rerun_all == 0:
	return mapsplice_out
 
    unmapped_option = "--output_unmapped";
    
    if do_fusion_and_is_paired:
	unmapped_option = "--output_unmapped_pe"
#noremap
#./bowtie -v 1 -f -k 40 -m 40 -p 8 --seg_len 25 --chrom_tab [chrom size file] 
#--split_index_path [split index path] --ref_seq_path [reference sequence path] 
#--mapsplice_out [output_file] /bioinfo/projects/zeng/mapsplice_multi/humanchridxallnohap 
#-1 /bioinfo/projects/kai/ground_truth/paired_reads.fa.1 -2 /bioinfo/projects/kai/ground_truth/paired_reads.fa.2 
#../bowtie.bwtout >../log.txt

#remap
#./bowtie -v 1 -f -k 40 -m 40 -p 8 --seg_len 25 --chrom_tab [chrom size file] 
#--split_index_path [split index path] --ref_seq_path [reference sequence path] 
#--juncdb_index [synthetic index basename] --mapsplice_out [output_file] 
#--optimize_repeats --output_unmapped /bioinfo/projects/zeng/mapsplice_multi/humanchridxallnohap 
#-1 /bioinfo/projects/kai/ground_truth/paired_reads.fa.1 -2 /bioinfo/projects/kai/ground_truth/paired_reads.fa.2 
#../bowtie.bwtout >../log.txt

    if min_map_len > 0:
	if juncdb_index == "":
	    if flankcase >= 5:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  "--min_map_len", str(min_map_len),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  ##"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  ##"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  "-1", input_reads_1, # Maxmimum intron length
		                  "-2", input_reads_2, # Maxmimum intron length
		                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  "--min_map_len", str(min_map_len),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  ##"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  ##"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  
		                  bowtie_index,
		                  input_reads_1, # Maxmimum intron length
		                  bowtie_output]
		    
	    else:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--double_anchor_noncanon",
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  "--min_map_len", str(min_map_len),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  #"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  ##"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  "-1", input_reads_1, # Maxmimum intron length
		                  "-2", input_reads_2, # Maxmimum intron length
		                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--double_anchor_noncanon",
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  "--min_map_len", str(min_map_len),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  #"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  ##"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  input_reads_1, # Maxmimum intron length
		                  bowtie_output]
	else:
	    if flankcase >= 5:
		if input_reads_2 != "":
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              ##"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              ##"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		else:
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              ##"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              ##"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]
	    else:
		if input_reads_2 != "":
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		else:
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              "--min_map_len", str(min_map_len),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]
		    
    else:
	if juncdb_index == "":
	    if flankcase >= 5:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  #"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  #"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  "-1", input_reads_1, # Maxmimum intron length
		                  "-2", input_reads_2, # Maxmimum intron length
		                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  #"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  #"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  input_reads_1, # Maxmimum intron length
		                  bowtie_output]
		    
	    else:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--double_anchor_noncanon",
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  #"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  #"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  "-1", input_reads_1, # Maxmimum intron length
		                  "-2", input_reads_2, # Maxmimum intron length
		                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "-v", str(mis_match),          # Output dir
		                  "--splice_only",
		                  read_format, # Anchor length
		                  "--double_anchor_noncanon",
		                  "--max_ins", str(max_insertion),
		                  "--max_del", str(max_deletion),
		                  #"--max_missing_seg", str(min_missed_seg),
		                  #"--max_repeat_all_hits", str(5000),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
		                  "-k", str(max_hits), # Mismatches allowed in extension
		                  "-m", str(max_hits), # Maxmimum intron length
		                  "-p", str(threads), # Minimum intron length
		                  "--seg_len", str(seg_len), # Seed size for reads
		                  "--chrom_tab", chromosome_size, # read width for reads
		                  #"--split_index_path", split_index_path, # islands extension
		                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
		                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                  bowtie_index,
		                  input_reads_1, # Maxmimum intron length
		                  bowtie_output]
	else:
	    if flankcase >= 5:
		if input_reads_2 != "":
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		else:
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]

	    else:
		if input_reads_2 != "":
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              "-1", input_reads_1, # Maxmimum intron length
			              "-2", input_reads_2, # Maxmimum intron length
			              bowtie_output]
		else:
		    if no_original_junc == False:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--juncdb_index", juncdb_index, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]
		    else:
			splice_cmd = [bin_dir + "mapsplice_multi_thread",
			              "--qual-scale", qual_scale,
			              "--min_intron", str(min_intron),
			              "--max_intron_single", str(max_intron_length),
			              "--min_len", str(min_read_len),
			              "-v", str(mis_match),          # Output dir
			              #"--splice_only",
			              read_format, # Anchor length
			              "--double_anchor_noncanon",
			              "--max_ins", str(max_insertion),
			              "--max_del", str(max_deletion),
			              #"--max_missing_seg", str(min_missed_seg),
			              #"--max_repeat_all_hits", str(5000),
			              "--max_append_mis", str( append_mismatch),
			              "--max_intron_double", str(max_intron_length),
			              "-k", str(max_hits), # Mismatches allowed in extension
			              "-m", str(max_hits), # Maxmimum intron length
			              "-p", str(threads), # Minimum intron length
			              "--seg_len", str(seg_len), # Seed size for reads
			              "--chrom_tab", chromosome_size, # read width for reads
			              #"--split_index_path", split_index_path, # islands extension
			              "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			              "--optimize_repeats",
			              unmapped_option, unmapped_reads,
			              "--mapsplice_out", mapsplice_out, # FASTA file extension
			              bowtie_index,
			              input_reads_1, # Maxmimum intron length
			              bowtie_output]


	    
 
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
	
	print >> sys.stderr, "[%s] return code " % str(retcode)
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Run MapSplice multi-thread failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_multi_thread not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    
    duration = finish_time - start_time

    return mapsplice_out



def call_mapsplice_multithreads_remap(mis_match, 
                                read_format,
                                max_hits,
                                threads,
                                seg_len,
                                chromosome_size,
                                split_index_path,
                                ref_seq_path,
                                flankcase,
                                max_insertion,
                                max_deletion,
                                min_missed_seg,
                                min_map_len,
                                mapsplice_out,
                                bowtie_index,
                                input_reads_1,
                                input_reads_2,
                                unmapped_reads,
                                bowtie_output,
                                juncdb_index,
                                max_intron_length,
                                append_mismatch,
                                log_file):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Run MapSplice multi-thread" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(log_file, "w")
    
    if os.path.exists(mapsplice_out + ".1") and \
	   rerun_all == 0:
	return mapsplice_out
#noremap
#./bowtie -v 1 -f -k 40 -m 40 -p 8 --seg_len 25 --chrom_tab [chrom size file] 
#--split_index_path [split index path] --ref_seq_path [reference sequence path] 
#--mapsplice_out [output_file] /bioinfo/projects/zeng/mapsplice_multi/humanchridxallnohap 
#-1 /bioinfo/projects/kai/ground_truth/paired_reads.fa.1 -2 /bioinfo/projects/kai/ground_truth/paired_reads.fa.2 
#../bowtie.bwtout >../log.txt

#remap
#./bowtie -v 1 -f -k 40 -m 40 -p 8 --seg_len 25 --chrom_tab [chrom size file] 
#--split_index_path [split index path] --ref_seq_path [reference sequence path] 
#--juncdb_index [synthetic index basename] --mapsplice_out [output_file] 
#--optimize_repeats --output_unmapped /bioinfo/projects/zeng/mapsplice_multi/humanchridxallnohap 
#-1 /bioinfo/projects/kai/ground_truth/paired_reads.fa.1 -2 /bioinfo/projects/kai/ground_truth/paired_reads.fa.2 
#../bowtie.bwtout >../log.txt

    if min_map_len > 0:
	if juncdb_index == "":
	    if flankcase >= 5:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  ##"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  ##"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  ##"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  ##"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
		                          
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
		    
	    else:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  ##"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  ##"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
	else:
	    if flankcase >= 5:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  ##"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  ##"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
	    else:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  "--min_map_len", str(min_map_len),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
    else:
	if juncdb_index == "":
	    if flankcase >= 5:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
		    
	    else:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          "--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
	else:
	    if flankcase >= 5:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]
	    else:
		if input_reads_2 != "":
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  "-1", input_reads_1, # Maxmimum intron length
			                  "-2", input_reads_2, # Maxmimum intron length
			                  bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
			                  "-v", str(mis_match),          # Output dir
		                          #"--splice_only",
			                  read_format, # Anchor length
			                  "--double_anchor_noncanon",
			                  "--max_ins", str(max_insertion),
			                  "--max_del", str(max_deletion),
			                  #"--max_missing_seg", str(min_missed_seg),
			                  #"--max_repeat_all_hits", str(5000),
			                  "--max_append_mis", str( append_mismatch),
			                  "--max_intron_double", str(max_intron_length),
			                  "-k", str(max_hits), # Mismatches allowed in extension
			                  "-m", str(max_hits), # Maxmimum intron length
			                  "-p", str(threads), # Minimum intron length
			                  "--seg_len", str(seg_len), # Seed size for reads
			                  "--chrom_tab", chromosome_size, # read width for reads
			                  #"--split_index_path", split_index_path, # islands extension
			                  "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			                  "--juncdb_index", juncdb_index, # block size for reading chromosome
			                  "--optimize_repeats",
			                  "--output_unmapped", unmapped_reads,
			                  "--mapsplice_out", mapsplice_out, # FASTA file extension
			                  bowtie_index,
			                  input_reads_1, # Maxmimum intron length
			                  bowtie_output]

	    
 
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
	
	print >> sys.stderr, "[%s] return code " % str(retcode)
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Run MapSplice multi-thread failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_multi_thread not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    
    duration = finish_time - start_time

    return mapsplice_out

def call_mapsplice_multithreads_fusion(mis_match, 
                                       read_format,
                                       max_hits,
                                       threads,
                                       seg_len,
                                       chromosome_size,
                                       split_index_path,
                                       ref_seq_path,
                                       flankcase,
                                       max_insertion,
                                       max_deletion,
                                       min_missed_seg,
                                       min_map_len,
                                       mapsplice_out,
                                       bowtie_index,
                                       input_reads_1,
                                       input_reads_2,
                                       unmapped_reads,
                                       bowtie_output,
                                       juncdb_index,
                                       fusion_juncdb_index,
                                       cluster,
                                       fusion,
                                       no_original_junc,
                                       no_original_fusion_junc,
                                       max_intron_length,
                                       append_mismatch,
                                       min_read_len,
                                       min_intron,
                                       qual_scale,
                                       log_file):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Run MapSplice multi-thread fusion" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(log_file, "w")
    
    if os.path.exists(fusion) and \
	   rerun_all == 0:
	return fusion
    
#./bowtie -v 1 -q -k 40 -m 40 -p 8 --seg_len 25 --max_ins 3 --max_del 10 --double_anchor_noncanon 
#--chrom_tab [chrom size file] --split_index_path [split index path] --ref_seq_path [reference sequence path] 
#--juncdb_index [synthetic index basename] --mapsplice_out [normal_alignment_output_file] 
#--cluster [fusion cluster file] --fusion [fusion alignment output file] bowtie_index 
#-1 unmapped_read1.fa -2 unmapped_read2.fa bowtie_output_file.bwtout >../log.txt

    if input_reads_2 != "":
	if fusion_juncdb_index == "":
	    if flankcase >= 5:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          #"--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		    
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	    else:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	else:
	    if flankcase >= 5:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		    
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	    else:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          "--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          "-1", input_reads_1, # Maxmimum intron length
			          "-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]

    else:
	if fusion_juncdb_index == "":
	    if flankcase >= 5:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          #"--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		    
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	    else:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--optimize_repeats",
			          #"--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	else:
	    if flankcase >= 5:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		    
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,	                      
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	    else:
		if no_original_junc == False:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--juncdb_index", juncdb_index, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
		else:
		    splice_cmd = [bin_dir + "mapsplice_multi_thread",
		                  "--qual-scale", qual_scale,
		                  "--min_intron", str(min_intron),
		                  "--max_intron_single", str(max_intron_length),
		                  "--min_len", str(min_read_len),
		                  "--max_append_mis", str( append_mismatch),
		                  "--max_intron_double", str(max_intron_length),
			          "-v", str(mis_match),          # Output dir
			          read_format, # Anchor length
			          "--double_anchor_noncanon",
		                  "--fusion_double_anchor_noncanon",
		                  "--fusion_single_anchor_noncanon",
			          "--max_ins", str(max_insertion),
			          "--max_del", str(max_deletion),
			          #"--min_map_len", str(min_map_len),
			          #"--max_missing_seg", str(min_missed_seg),
			          "-k", str(max_hits), # Mismatches allowed in extension
			          "-m", str(max_hits), # Maxmimum intron length
			          "-p", str(threads), # Minimum intron length
			          "--seg_len", str(seg_len), # Seed size for reads
			          "--chrom_tab", chromosome_size, # read width for reads
			          #"--split_index_path", split_index_path, # islands extension
			          "--ref_seq_path", ref_seq_path, # block size for reading chromosome
			          "--fusiondb_index", fusion_juncdb_index,
			          "--optimize_repeats",
			          "--output_unmapped", unmapped_reads,
			          "--mapsplice_out", mapsplice_out, # FASTA file extension
			          #"--cluster", cluster, # FASTA file extension
			          "--fusion", fusion, # FASTA file extension
			          bowtie_index,
			          input_reads_1, # Maxmimum intron length
			          #"-2", input_reads_2, # Maxmimum intron length
			          bowtie_output]
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Run MapSplice multi-thread fusion failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_multi_thread not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    
    duration = finish_time - start_time

    return fusion

def run_alignment_handler(junction_file, 
                   is_paired,
                   max_mate_dist,
                   max_hits,
                   filtered_alignment_base,
                   max_read_length,
                   chromosome_dir,
                   filter_flag,
                   min_insertion,
                   max_deletion,
                   min_anchor,
                   min_junction_anchor,
                   min_mismatch,
                   add_soft_clip,
		   mate_dist_sd,
		   max_anchor_diff,
		   chromosome_size_file,
		   encompassing_fusion_region_extension,
		   input_sam_file,
                   log_file):
    
    start_time = datetime.now()
    
    print >> sys.stderr, "[%s] Run alignment handler" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    filtered_alignment_file = input_sam_file + filtered_alignment_base;
    
    if os.path.exists(filtered_alignment_file) and \
	   rerun_all == 0:
	return filtered_alignment_file

    mapsplice_log = open(log_file, "w")
    
    splice_cmd = [bin_dir + "alignment_handler",
                  junction_file,
                  str(is_paired), # Anchor length
                  str(max_mate_dist), # Mismatches allowed in extension
                  str(max_hits), # Maxmimum intron length
                  filtered_alignment_base, # Minimum intron length
                  str(max_read_length), # Seed size for reads
                  chromosome_dir, # read width for reads
                  str(filter_flag), # islands extension
                  str(min_insertion), # block size for reading chromosome
                  str(max_deletion), # nothing important
                  str(min_anchor), #if is 1, only output flank string that is not case 0
                  str(min_junction_anchor), #number of anchors
                  str(min_mismatch), #number of segments
                  str(add_soft_clip), #segment length
                  str(mate_dist_sd), #segment length
                  str(max_anchor_diff), #segment length
                  chromosome_size_file, #segment length
                  str(encompassing_fusion_region_extension), 
                  input_sam_file, #tmp dir
                  number_of_threads
                  ]

	
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % splice_cmd
    
    try:
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: alignment_handler not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    
    duration = finish_time - start_time
  
    return filtered_alignment_file


def run_alignment_handler_multi(junction_file, 
                   is_paired,
                   max_mate_dist,
                   max_hits,
                   filtered_alignment_base,
                   filtered_alignment_append,
                   max_read_length,
                   chromosome_dir,
                   filter_flag,
                   min_insertion,
                   max_deletion,
                   min_anchor,
                   min_junction_anchor,
                   min_mismatch,
                   add_soft_clip,
		   mate_dist_sd,
                   intron_dist_sd,
		   max_anchor_diff,
		   chromosome_size_file,
		   encompassing_fusion_region_extension,
                   number_of_threads,
                   min_coverage,
                   fragment_length,
                   fragment_length_sd,
                   avearge_fragment_length,
                   boundary,
                   min_isoform_length,
                   min_encompass_count,
                   min_entropy,
		   input_sam_file,                   
                   log_file,
                   err_file):
    
    start_time = datetime.now()
    
    print >> sys.stderr, "[%s] Run alignment handler" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    filtered_alignment_file = input_sam_file + filtered_alignment_append;
    
    if os.path.exists(filtered_alignment_base + ".fil.junc") and \
	   rerun_all == 0:
	return filtered_alignment_file

    mapsplice_log = open(log_file, "w");
    
    err_log = open(err_file, "w");
    
    splice_cmd = [bin_dir + "alignment_handler_multi",
                  junction_file,
                  str(is_paired), # Anchor length
                  str(max_mate_dist), # Mismatches allowed in extension
                  str(max_hits), # Maxmimum intron length
                  filtered_alignment_base, # Minimum intron length
                  str(max_read_length), # Seed size for reads
                  chromosome_dir, # read width for reads
                  str(filter_flag), # islands extension
                  str(min_insertion), # block size for reading chromosome
                  str(max_deletion), # nothing important
                  str(min_anchor), #if is 1, only output flank string that is not case 0
                  str(min_junction_anchor), #number of anchors
                  str(min_mismatch), #number of segments
                  str(add_soft_clip), #segment length
                  str(mate_dist_sd), #segment length
                  str(max_anchor_diff), #segment length
                  chromosome_size_file, #segment length
                  str(encompassing_fusion_region_extension), 
                  str(number_of_threads),
                  str(intron_dist_sd),
                  str(min_coverage),
                  str(fragment_length),
                  str(fragment_length_sd),
                  str(avearge_fragment_length),
                  str(boundary),
                  str(min_isoform_length),
                  str(min_encompass_count),
                  str(min_entropy),
                  input_sam_file #tmp dir
                  ]

	
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % splice_cmd
    
    try:
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log, stderr=err_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter alignment failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: alignment_handler not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    
    duration = finish_time - start_time
  
    return filtered_alignment_file

def filterbyrepeats(in_fusion,
                    in_normal,
                    out_fusion,
                    chromosomes_dir,
                    seg_len,
                    max_repeat_length,
                    max_anchor_diff,
                    min_seg_len,
                    min_entropy,
                    threads,
                    bowidx,
                    chromosome_size_file,
                    logging_dir):
    
    start_time = datetime.now()
    
    print >> sys.stderr, "[%s] Filter fusion by repeats" % start_time.strftime("%c")
    
    if seg_len < 25:
	seg_len = 25
    #splice_cmd = ""
    
    #filtered_alignment_file = input_sam_file + filtered_alignment_append;
    
    if os.path.exists(out_fusion) and \
	   rerun_all == 0:
	return out_fusion

    ##############match fusion to normal junction#################
    
    match_junc_log = open(logging_dir + "matchfusion2normal.log", "w");
    
    match_junc_cmd = [bin_dir + "matchfusion2normal",
                      in_normal,
                      in_fusion,
                      filter_repeats_dir + "fusions_matched_normal",
                      str(10)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % match_junc_cmd
	
    try:
        retcode = subprocess.call(match_junc_cmd, stdout=match_junc_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: match fusion to normal failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: matchfusion2normal not found on this system"
        exit(1)
	
    ##############load fusion chromosome seq##################
    
    load_fusion_chrom_seq_std_log = open(logging_dir + "load_fusion_chrom_seq_std.log", "w");
    
    load_fusion_chrom_seq_std_cmd = [bin_dir + "load_fusion_chrom_seq_std",
                      chromosomes_dir,
                      str(seg_len),
                      str(max_anchor_diff),
                      str(min_seg_len),
                      filter_repeats_dir + "fusions_matched_normal.matched_junc",
                      filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq",
                      str(1),
                      str(min_entropy)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % load_fusion_chrom_seq_std_cmd
    
    try:
        retcode = subprocess.call(load_fusion_chrom_seq_std_cmd, stdout=load_fusion_chrom_seq_std_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: add chromosome sequence failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: load_fusion_chrom_seq_std not found on this system"
        exit(1)
	
	
    ##############align anchor reads by bowtie##################
    
    anchor_bowtie_log = open(logging_dir + "anchor_bowtie.log", "w");
    
    anchor_bowtie_err = open(logging_dir + "anchor_bowtie.err", "w");
    
    readslist = filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq.1.fa" + "," + filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq.2.fa";
    
    anchor_bowtie_cmd = [bin_dir + "bowtie",
                         "-S",
                         "--sam-nohead",
                         "-f",
                         "--threads", str(threads),
                         "-X", str(50000),
                         "--un", filter_repeats_dir + "unmapped_segments",
                         '-a', #"-k", str(40000),
                         "-m", str(40001),
                         "-v", str(3),
                         "--max", filter_repeats_dir + "repeat_segments",
                         bowidx,
                         readslist,
                         filter_repeats_dir + "repeats.sam"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % anchor_bowtie_cmd
	
    try:
        retcode = subprocess.call(anchor_bowtie_cmd, stdout=anchor_bowtie_log, stderr=anchor_bowtie_err)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: map repeat reads failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
        exit(1)
	
	
    ##############sort repeat alignments##################
    
    sort_repeat_log = open(logging_dir + "sort_repeat.log", "w");
    
    sort_repeat_cmd = ["sort",
                       "-k1,1n",
                       "-t:",
                       "-o", filter_repeats_dir + "repeats_sorted.sam",
                       filter_repeats_dir + "repeats.sam"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % sort_repeat_cmd
	
    try:
        retcode = subprocess.call(sort_repeat_cmd, stdout=sort_repeat_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort repeat reads failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    
	
    ##############load fusion chromosome long seq##################
    
    load_fusion_chrom_seq_std_long_log = open(logging_dir + "load_fusion_chrom_seq_std_long_seq.log", "w");
    
    load_fusion_chrom_seq_std_long_cmd = [bin_dir + "load_fusion_chrom_seq_std",
                      chromosomes_dir,
                      str(max_repeat_length),
                      str(max_anchor_diff),
                      str(min_seg_len),
                      filter_repeats_dir + "fusions_matched_normal.matched_junc",
                      filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq_long",
                      str(1),
                      str(min_entropy),
                      str(1)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % load_fusion_chrom_seq_std_long_cmd
    
    try:
        retcode = subprocess.call(load_fusion_chrom_seq_std_long_cmd, stdout=load_fusion_chrom_seq_std_long_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: add chromosome long sequence failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: load_fusion_chrom_seq_std not found on this system"
        exit(1)
	
	
	
    ##############align long fusion junction sequence by bowtie##################
    
    long_anchor_bowtie_log = open(logging_dir + "long_fusion_seq_bowtie.log", "w");
    
    long_anchor_bowtie_err = open(logging_dir + "long_fusion_seq_bowtie.err", "w");
    
    readslonglist = filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq_long.1.fa" + "," + filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq_long.2.fa";
    
    long_anchor_bowtie_cmd = [bin_dir + "bowtie",
                         "-f",
                         "--threads", str(threads),
                         "--un", filter_repeats_dir + "unmapped_segments_long",
                         "-k", str(1),
                         "-m", str(1),
                         "-v", str(2),
                         "--max", filter_repeats_dir + "repeat_segments_long",
                         bowidx,
                         readslonglist,
                         filter_repeats_dir + "repeats_long.sam"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % long_anchor_bowtie_cmd
	
    try:
        retcode = subprocess.call(long_anchor_bowtie_cmd, stdout=long_anchor_bowtie_log, stderr=long_anchor_bowtie_err)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: map long fusion sequence failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie not found on this system"
        exit(1)
	
	
	
    ##############pair repeat reads##################
    
    mapsplice_log = open(logging_dir + "pair_repeat_reads.log", "w");
    
    err_log = open(logging_dir + "pair_repeat_reads.err", "w");
    
    splice_cmd = [bin_dir + "alignment_handler_multi",
                  "",
                  str(1), # Anchor length
                  str(50000), # Mismatches allowed in extension
                  str(4000), # Maxmimum intron length
                  filter_repeats_dir + "_filtered_normal_alignments_fix_pair", # Minimum intron length
                  str(seg_len), # Seed size for reads
                  "", # read width for reads
                  str(2348), # islands extension
                  str(6), # block size for reading chromosome
                  str(10), # nothing important
                  str(0), #if is 1, only output flank string that is not case 0
                  str(10), #number of anchors
                  str(5), #number of segments
                  str(1), #segment length
                  str(100), #segment length
                  str(50), #segment length
                  chromosome_size_file, #segment length
                  str(50000), 
                  str(threads),
                  str(500),
                  str(2),
                  str(400),
                  str(100),
                  str(225),
                  str(36),
                  str(100),
                  str(1),
                  str(0),
                  filter_repeats_dir + "repeats_sorted.sam" #tmp dir
                  ]

	
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % splice_cmd
    
    try:
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log, stderr=err_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pair repeat reads failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: alignment_handler not found on this system"
        exit(1)
       
    ##############sort repeat alignments##################
    
    FilterFusionByNormalPaired_log = open(logging_dir + "FilterFusionByNormalPaired.log", "w");
    
    FilterFusionByNormalPaired_cmd = [bin_dir + "FilterFusionByNormalPaired",
                                      filter_repeats_dir + "repeats_sorted.sam_filtered_normal_alignments_fix_pair.bothunspliced.paired",
                                      filter_repeats_dir + "fusions_matched_normal.matched_junc.add_chrom_seq",
                                      out_fusion,
                                      filter_repeats_dir + "repeat_segments",
                                      filter_repeats_dir + "repeat_segments_long",
                                      filter_repeats_dir + "fusion_repeats_filtered.txt"]
    
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % FilterFusionByNormalPaired_cmd
	
    try:
        retcode = subprocess.call(FilterFusionByNormalPaired_cmd, stdout=FilterFusionByNormalPaired_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filter fusion by normalPaired failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterFusionByNormalPaired not found on this system"
        exit(1)
	
    finish_time = datetime.now()
    
    duration = finish_time - start_time
  
    return filter_repeats_dir + "repeats_sorted.sam_filtered_normal_alignments_fix_pair"


def filterbyannotation(in_fusion,
                       out_fusion,
                       out_fusion_not_well_annotated,
                       out_fusion_combined,
                       circular_RNAs,
                       gtf_file,
                       logging_dir):
    
    start_time = datetime.now()
    
    print >> sys.stderr, "[%s] Filter fusion by annotation" % start_time.strftime("%c")
    
    #splice_cmd = ""
    
    #filtered_alignment_file = input_sam_file + filtered_alignment_append;
    
    if os.path.exists(out_fusion) and \
	   rerun_all == 0:
	return out_fusion

    ##############convert gene gtf file to annotation format#################
    
    gtf2genetab_log = open(logging_dir + "gtf2genetab.log", "w");
    
    gtf2genetab_cmd = [bin_dir + "gtf2genetab",
                       gtf_file,
                       temp_dir + "gene_exons.txt"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % gtf2genetab_cmd
	
    try:
        retcode = subprocess.call(gtf2genetab_cmd, stdout=gtf2genetab_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert gene gtf file to annotation failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: gtf2genetab not found on this system"
        exit(1)
	
    ##############annotate fusion junctions##################
    
    annotate_fusion_log = open(logging_dir + "annotate_fusion.log", "w");
    
    annotate_fusion_cmd = [bin_dir + "search_fusion_gene",
                           "-g", temp_dir + "gene_exons.txt",
                           "-f", in_fusion,
                           "-o", filter_repeats_dir + "filter_repeats.txt.annot"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % annotate_fusion_cmd
    
    try:
        retcode = subprocess.call(annotate_fusion_cmd, stdout=annotate_fusion_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: annotate fusion junction failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: search_fusion_gene not found on this system"
        exit(1)

    ##############add fusion strand consistent##################
    
    add_fusion_strand_consistent_log = open(logging_dir + "add_fusion_strand_consistent.log", "w");
    
    add_fusion_strand_consistent_err = open(logging_dir + "add_fusion_strand_consistent.err", "w");

    add_fusion_strand_consistent_cmd = [bin_dir + "AddFusionStrandConsistent",
                                        filter_repeats_dir + "filter_repeats.txt.annot",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % add_fusion_strand_consistent_cmd
	
    try:
        retcode = subprocess.call(add_fusion_strand_consistent_cmd, stdout=add_fusion_strand_consistent_log, stderr=add_fusion_strand_consistent_err)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: add fusion strand consistent failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: AddFusionStrandConsistent not found on this system"
        exit(1)
	

    ##############separate matched strand fusion from fusion junction##################
    
    SeparateNormalFromFusionJunc_log = open(logging_dir + "SeparateNormalFromFusionJunc.log", "w");
    
    SeparateNormalFromFusionJunc_cmd = [bin_dir + "SeparateNormalFromFusionJunc",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.match",
                                        out_fusion,#output_dir + "fusion.from.fusion",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.match.normal",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.match.circularRNAs"]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % SeparateNormalFromFusionJunc_cmd
	
    try:
        retcode = subprocess.call(SeparateNormalFromFusionJunc_cmd, stdout=SeparateNormalFromFusionJunc_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate matched strand fusion from fusion junction failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SeparateNormalFromFusionJunc not found on this system"
        exit(1)
	
	
    ##############separate normal junction and circular RNAs from fusion junction##################
    
    SeparateNormalFromFusionJunc_log = open(logging_dir + "SeparateNormalFromFusionJunc.log", "w");
    
    SeparateNormalFromFusionJunc_cmd = [bin_dir + "SeparateNormalFromFusionJunc",
                                        filter_repeats_dir + "filter_repeats.txt.annot.strand.notmatch",
                                        out_fusion_not_well_annotated,#output_dir + "fusion.from.fusion",
                                        filter_repeats_dir + "normal.from.fusion",
                                        circular_RNAs]
    
    
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % SeparateNormalFromFusionJunc_cmd
	
    try:
        retcode = subprocess.call(SeparateNormalFromFusionJunc_cmd, stdout=SeparateNormalFromFusionJunc_log)     
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate normal junction and circular RNAs from fusion junction failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SeparateNormalFromFusionJunc not found on this system"
        exit(1)
	

    ###############merge not well annotated and well annotated fusions
    
    sort_annotated_log = open(logging_dir + "sort_annotated.log", "w");
    
    sort_repeat_cmd = ["sort",
                       "-k1,1",
                       "-k2,2n",
                       "-k3,3n",
                       "-t:",
                       "-o", filter_repeats_dir + "combine_well_annotated_and_not_well_annodated.txt",
                       out_fusion_not_well_annotated,
                       out_fusion]
    
    if DEBUG == 1:
	    print >> sys.stderr, "[%s]" % sort_repeat_cmd
    
    try:
	retcode = subprocess.call(sort_repeat_cmd, stdout=sort_annotated_log)     
	
	if retcode != 0:
	    print >> sys.stderr, fail_str, "Error: sort annotated failed"
	    exit(1)
	   
    except OSError, o:
	if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
	    print >> sys.stderr, fail_str, "Error: sort not found on this system"
	exit(1)
		
    ###############separate fusion to inter chromosome and long range##################
	
    sepMPSfusion_log = open(logging_dir + "sepMPSfusion.log", "w");
    
    sepMPSfusion_cmd = [bin_dir + "sepMPSfusion",
                       filter_repeats_dir + "combine_well_annotated_and_not_well_annodated.txt",
                       "1000000"]
    
    if DEBUG == 1:
	    print >> sys.stderr, "[%s]" % sepMPSfusion_cmd
	    
    try:
	retcode = subprocess.call(sepMPSfusion_cmd, stdout=sepMPSfusion_log)     
	
	if retcode != 0:
	    print >> sys.stderr, fail_str, "Error: sepMPSfusion failed"
	    exit(1)
	   
    except OSError, o:
	if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
	    print >> sys.stderr, fail_str, "Error: sepMPSfusion not found on this system"
	exit(1)       

    ###############merge not well annotated and well annotated fusions
	
    out_fusion_combined_log = open(out_fusion_combined, "w");
    
    cat_cmd = ["cat",
               filter_repeats_dir + "combine_well_annotated_and_not_well_annodated.txt.longrange",
               filter_repeats_dir + "combine_well_annotated_and_not_well_annodated.txt.interchr"]
    
    if DEBUG == 1:
	    print >> sys.stderr, "[%s]" % cat_cmd
	    
    try:
	retcode = subprocess.call(cat_cmd, stdout=out_fusion_combined_log)     
	
	if retcode != 0:
	    print >> sys.stderr, fail_str, "Error: cat failed"
	    exit(1)
	   
    except OSError, o:
	if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
	    print >> sys.stderr, fail_str, "Error: cat not found on this system"
	exit(1)     
	
    finish_time = datetime.now()
    
    duration = finish_time - start_time
  
    return output_dir + "fusion.from.fusion"

def call_mapsplice_segment_fusion(islands_gff, 
                   seed_length,
                   read_width,
                   min_anchor_len,
                   splice_mismatches,
                   min_intron_length,
                   max_intron_length,
                   islands_extension,
                   flank_case,
                   rank,
                   FASTA_file_ext,
                   reads_files_dir,
                   chromosomes_file_dir,
                   num_anchor,
		   num_seg,
		   seg_len,
		   bwtout_25,
		   fix_hole_file,
		   cur_ouput_dir,
		   fqreads,
		   extend_bits,
		   total_mismatch,
		   total_fusion_mismatch,
		   append_mismatch,
		   prefix_match,
		   threads,
		   fusion_junction
		   ):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning fusion spliced reads" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(logging_dir + "mapsplice_fusion.log", "w")
    
    juncfile = original_dir + "canonical_junctions.txt"
    
    fqreads = "";
    
    test_file =  fusion_junction
    if os.path.exists(test_file) and \
	   rerun_all == 0:
	return test_file
    if read_width > 0:
	
	if fqreads == "":
	    splice_cmd = [bin_dir + "mapsplice_segments",
			  "-v", cur_ouput_dir,          # Output dir
			  "-n", str(min_anchor_len), # Anchor length
			  "-m", str(splice_mismatches), # Mismatches allowed in extension
			  "-x", str(max_intron_length), # Maxmimum intron length
			  "-i", str(min_intron_length), # Minimum intron length
			  "-h", str(seed_length), # Seed size for reads
			  "-w", str(read_width), # read width for reads
			  "-p", str(islands_extension), # islands extension
			  "-b", str(350000000), # block size for reading chromosome
			  "-s", FASTA_file_ext, # FASTA file extension
			  "-t", islands_gff, # island location
			  #"-u", reads_files_dir, # reads file or directory
			  "-c", chromosomes_file_dir, # chromosomes file or directory
			  "-y", str(prefix_match), # nothing important
			  "-f", str(flank_case), #if is 1, only output flank string that is not case 0
			  "-A", str(num_anchor), #number of anchors
			  #"-R", str(rank),   #only output rank >= rank
			  "-bwtout", bwtout_25, #bwt segment output
			  "-G", str(num_seg), #number of segments
			  "-L", str(seg_len), #segment length
			  #"-H", fix_hole_file, #fix hole file
			  "-E", str(extend_bits), #extend bits when fix hole and fusion
			  #"-M", str(total_mismatch), #total mismatch on splice reads
			  #"-FM", str(total_fusion_mismatch), #total fusion mismatch
			  #"-P", str(append_mismatch), #append mismatch
			  "-threads", str(threads), #number of threads used
			  "-tmp", temp_dir, #tmp dir
			  "-fusion", fusion_dir, #fusion dir
			  "-fusionjunc", fusion_junction #fusion junction
			  ]
	else:
	    splice_cmd = [bin_dir + "mapsplice_segments",
			  "-v", cur_ouput_dir,          # Output dir
			  "-n", str(min_anchor_len), # Anchor length
			  "-m", str(splice_mismatches), # Mismatches allowed in extension
			  "-x", str(max_intron_length), # Maxmimum intron length
			  "-i", str(min_intron_length), # Minimum intron length
			  "-h", str(seed_length), # Seed size for reads
			  "-w", str(read_width), # read width for reads
			  "-p", str(islands_extension), # islands extension
			  "-b", str(350000000), # block size for reading chromosome
			  "-s", FASTA_file_ext, # FASTA file extension
			  "-t", islands_gff, # island location
			  #"-u", reads_files_dir, # reads file or directory
			  "-c", chromosomes_file_dir, # chromosomes file or directory
			  "-y", str(1), # nothing important
			  "-f", str(flank_case), #if is 1, only output flank string that is not case 0
			  "-A", str(num_anchor), #number of anchors
			  #"-R", str(rank),   #only output rank >= rank
			  "-bwtout", bwtout_25, #bwt segment output
			  "-G", str(num_seg), #number of segments
			  "-L", str(seg_len), #segment length
			  #"-H", fix_hole_file, #fix hole file
			  "-E", str(extend_bits), #extend bits when fix hole and fusion
			  #"-M", str(total_mismatch), #total mismatch on splice reads
			  #"-FM", str(total_fusion_mismatch), #total fusion mismatch
			  #"-P", str(append_mismatch), #append mismatch
			  "-threads", str(threads), #number of threads used
			  "-tmp", temp_dir, #tmp dir
			  "-fusion", fusion_dir,  #fusion dir
			  "-fusionjunc", fusion_junction #fusion junction
			  #"-FQ", fqreads #fastq reads
			  ]
    else:
        splice_cmd = [bin_dir + "mapsplice",
                      "-v", cur_ouput_dir,          # Output dir
                      "-n", str(min_anchor_len), # Anchor length
                      "-m", str(splice_mismatches), # Mismatches allowed in extension
                      "-x", str(max_intron_length), # Maxmimum intron length
                      "-i", str(min_intron_length), # Minimum intron length
                      "-h", str(seed_length), # Seed size for reads
                      "-w", str(read_width), # read width for reads
                      "-p", str(islands_extension), # islands extension
                      "-b", str(350000000), # block size for reading chromosome
                      "-s", FASTA_file_ext, # FASTA file extension
                      "-t", islands_gff, # island location
                      "-u", reads_files_dir, # reads file or directory
                      "-c", chromosomes_file_dir, # chromosomes file or directory
                      "-y", str(1), # nothing important
                      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
                      "-A", str(num_anchor), #number of anchors
                      "-R", str(rank)   #only output rank >= rank
                      ]
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] separate unique and multiple mapped reads" % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MapSplice not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    
    return juncfile

def call_mapsplice_segment_fusion_single_anchor(islands_gff, 
                   seed_length,
                   read_width,
                   min_anchor_len,
                   splice_mismatches,
                   min_intron_length,
                   max_intron_length,
                   islands_extension,
                   flank_case,
                   rank,
                   FASTA_file_ext,
                   reads_files_dir,
                   chromosomes_file_dir,
                   num_anchor,
		   num_seg,
		   seg_len,
		   bwtout_25,
		   fix_hole_file,
		   cur_ouput_dir,
		   fqreads,
		   extend_bits,
		   total_mismatch,
		   total_fusion_mismatch,
		   append_mismatch,
		   prefix_match,
		   threads,
		   fusion_junction,
		   fusion_single_anchor,
		   ):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning fusion spliced reads" % start_time.strftime("%c")
    
    splice_cmd = ""
    
    mapsplice_log = open(logging_dir + "mapsplice_fusion.log", "w")
    
    juncfile = original_dir + "canonical_junctions.txt"
    
    fqreads = "";
    
    test_file =  fusion_junction
    if os.path.exists(test_file) and \
	   rerun_all == 0:
	return test_file
    if read_width > 0:
	
	if fqreads == "":
	    splice_cmd = [bin_dir + "mapsplice_segments",
			  "-v", cur_ouput_dir,          # Output dir
			  "-n", str(min_anchor_len), # Anchor length
			  "-m", str(splice_mismatches), # Mismatches allowed in extension
			  "-x", str(max_intron_length), # Maxmimum intron length
			  "-i", str(min_intron_length), # Minimum intron length
			  "-h", str(seed_length), # Seed size for reads
			  "-w", str(read_width), # read width for reads
			  "-p", str(islands_extension), # islands extension
			  "-b", str(350000000), # block size for reading chromosome
			  "-s", FASTA_file_ext, # FASTA file extension
			  "-t", islands_gff, # island location
			  #"-u", reads_files_dir, # reads file or directory
			  "-c", chromosomes_file_dir, # chromosomes file or directory
			  "-y", str(prefix_match), # nothing important
			  "-f", str(flank_case), #if is 1, only output flank string that is not case 0
			  "-A", str(num_anchor), #number of anchors
			  #"-R", str(rank),   #only output rank >= rank
			  "-bwtout", bwtout_25, #bwt segment output
			  "-G", str(num_seg), #number of segments
			  "-L", str(seg_len), #segment length
			  #"-H", fix_hole_file, #fix hole file
			  "-E", str(extend_bits), #extend bits when fix hole and fusion
			  #"-M", str(total_mismatch), #total mismatch on splice reads
			  #"-FM", str(total_fusion_mismatch), #total fusion mismatch
			  #"-P", str(append_mismatch), #append mismatch
			  "-threads", str(threads), #number of threads used
			  "-tmp", temp_dir, #tmp dir
			  "-fusion", fusion_dir, #fusion dir
			  "-fusionjunc", fusion_junction, #fusion junction
			  "-fusionsingleanchor", fusion_single_anchor #fusion junction
			  ]
	else:
	    splice_cmd = [bin_dir + "mapsplice_segments",
			  "-v", cur_ouput_dir,          # Output dir
			  "-n", str(min_anchor_len), # Anchor length
			  "-m", str(splice_mismatches), # Mismatches allowed in extension
			  "-x", str(max_intron_length), # Maxmimum intron length
			  "-i", str(min_intron_length), # Minimum intron length
			  "-h", str(seed_length), # Seed size for reads
			  "-w", str(read_width), # read width for reads
			  "-p", str(islands_extension), # islands extension
			  "-b", str(350000000), # block size for reading chromosome
			  "-s", FASTA_file_ext, # FASTA file extension
			  "-t", islands_gff, # island location
			  #"-u", reads_files_dir, # reads file or directory
			  "-c", chromosomes_file_dir, # chromosomes file or directory
			  "-y", str(1), # nothing important
			  "-f", str(flank_case), #if is 1, only output flank string that is not case 0
			  "-A", str(num_anchor), #number of anchors
			  #"-R", str(rank),   #only output rank >= rank
			  "-bwtout", bwtout_25, #bwt segment output
			  "-G", str(num_seg), #number of segments
			  "-L", str(seg_len), #segment length
			  #"-H", fix_hole_file, #fix hole file
			  "-E", str(extend_bits), #extend bits when fix hole and fusion
			  #"-M", str(total_mismatch), #total mismatch on splice reads
			  #"-FM", str(total_fusion_mismatch), #total fusion mismatch
			  #"-P", str(append_mismatch), #append mismatch
			  "-threads", str(threads), #number of threads used
			  "-tmp", temp_dir, #tmp dir
			  "-fusion", fusion_dir,  #fusion dir
			  "-fusionjunc", fusion_junction, #fusion junction
			  #"-FQ", fqreads #fastq reads
			  "-fusionsingleanchor", fusion_single_anchor #fusion junction
			  ]
    else:
        splice_cmd = [bin_dir + "mapsplice",
                      "-v", cur_ouput_dir,          # Output dir
                      "-n", str(min_anchor_len), # Anchor length
                      "-m", str(splice_mismatches), # Mismatches allowed in extension
                      "-x", str(max_intron_length), # Maxmimum intron length
                      "-i", str(min_intron_length), # Minimum intron length
                      "-h", str(seed_length), # Seed size for reads
                      "-w", str(read_width), # read width for reads
                      "-p", str(islands_extension), # islands extension
                      "-b", str(350000000), # block size for reading chromosome
                      "-s", FASTA_file_ext, # FASTA file extension
                      "-t", islands_gff, # island location
                      "-u", reads_files_dir, # reads file or directory
                      "-c", chromosomes_file_dir, # chromosomes file or directory
                      "-y", str(1), # nothing important
                      "-f", str(flank_case), #if is 1, only output flank string that is not case 0
                      "-A", str(num_anchor), #number of anchors
                      "-R", str(rank)   #only output rank >= rank
                      ]
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] separate unique and multiple mapped reads" % splice_cmd
    
    try:    
        retcode = subprocess.call(splice_cmd, stdout=mapsplice_log)#
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Spliced read alignment failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MapSplice not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    
    return juncfile

def extract_maxlen(log_file):
    fh = open(log_file,"r")
    igot = fh.readlines()
    
    maxlen = 0;
    
    read_format = "";
    
    qual_scale = "";
    
    for line in igot:
	if line[0] == '\n' or line[0] == '#':
	    comments = line;
	else:
	    line = line.rstrip();
	    command = line.split(' ');
	    
	    if command[0] == "maxlen":
		maxlen = int(command[1]);
		
	    if command[0] == "read_format":
		read_format = command[1];
		
	    if command[0] == "qual_scale":
		qual_scale = command[1];
		
	
    if maxlen != 0 and read_format != "":
	return (maxlen, read_format, qual_scale);
    else:
	print >> sys.stderr, fail_str, "Error: No max read length or read_format found";
    
    exit(1);

def separatedmultipleunique(mps_unique_mapreads, mps_multiple_mapreads, cur_output_dir):
    print >> sys.stderr, "[%s] separate unique and multiple mapped reads" % (right_now())

    separatemapreads_cmd = [bin_dir + "separatemapreads", 
		        cur_output_dir,
			"fixed_fixhole_divided_reads.txt",
			"fixed_hole_divided_reads.txt",
			"mapreads_divided_reads.txt",
			"fix_head_tail_divided_reads.txt",
			mps_unique_mapreads,
			mps_multiple_mapreads,
                        "1"]
       
    if os.path.exists(cur_output_dir + mps_unique_mapreads) and \
           os.path.exists(cur_output_dir + mps_multiple_mapreads) and \
	   rerun_all == 0:
	return mps_unique_mapreads
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separated mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separatemapreads not found on this system"
        exit(1)
    return (mps_unique_mapreads, mps_multiple_mapreads);

def merge_sort_sam(merged_sorted_by_tagname, merged_sorted_by_chromooffset, bwtmapped_sam, bwt_25bp_mapped, stat_file):
    print >> sys.stderr, "[%s] merge sort mapped reads" % (right_now())

    bwt_25bp_mapped_mishandt_matched = bwt_25bp_mapped + ".mistailandhead.matched";
    
    bwt_25bp_mapped_mishort_matched = bwt_25bp_mapped + ".mistailorhead.matched";
    
    bwt_25bp_mapped_allmapped = bwt_25bp_mapped + ".allmapped";
    
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    separatemapreads_cmd = [bin_dir + "merge_sort_sam", 
			temp_dir + "splicedreads_remdupdivided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_divided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_f0_divided_reads.txt",
			temp_dir + "fixed_fixhole_f0_divided_reads.txt",
			temp_dir + "fixed_hole_f0_divided_reads.txt",
			temp_dir + "fix_head_tail_noncanon_divided_reads.txt",
			temp_dir + "mapreads_noncanon_divided_reads.txt",
			bwtmapped_sam,
			bwt_25bp_mapped_mishandt_matched,
			bwt_25bp_mapped_mishort_matched,
			bwt_25bp_mapped_allmapped,			
			merged_sorted_by_tagname,
			merged_sorted_by_chromooffset,
			stat_file]
       
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort mapped reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sort_sam not found on this system"
        exit(1)
    return (merged_sorted_by_tagname, merged_sorted_by_chromooffset);

def merge_sort_sam_unspliced(merged_sorted_by_tagname, merged_sorted_by_chromooffset, bwtmapped_sam, bwt_25bp_mapped, stat_file):
    print >> sys.stderr, "[%s] merge sort mapped reads unspliced" % (right_now())

    bwt_25bp_mapped_mishandt_matched = bwt_25bp_mapped + ".mistailandhead.matched";
    
    bwt_25bp_mapped_mishort_matched = bwt_25bp_mapped + ".mistailorhead.matched";
    
    bwt_25bp_mapped_allmapped = bwt_25bp_mapped + ".allmapped";
    
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	
	temp_fs = open(merged_sorted_by_tagname, "w")
	temp_fs.close()
	
	temp_fs = open(bwtmapped_sam, "w")
	temp_fs.close()
	
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    separatemapreads_cmd = [bin_dir + "merge_sort_sam", 
			bwtmapped_sam,
			bwt_25bp_mapped_mishandt_matched,
			bwt_25bp_mapped_mishort_matched,
			bwt_25bp_mapped_allmapped,			
			merged_sorted_by_tagname,
			merged_sorted_by_chromooffset,
			stat_file]
       
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	
	temp_fs = open(merged_sorted_by_tagname, "w")
	temp_fs.close()
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort mapped reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sort_sam not found on this system"
        exit(1)
	
    temp_fs = open(merged_sorted_by_tagname, "w")
    temp_fs.close()
    
    temp_fs = open(bwtmapped_sam, "w")
    temp_fs.close()
    
    return (merged_sorted_by_tagname, merged_sorted_by_chromooffset);

def merge_sort_sam_canon(merged_sorted_by_tagname, merged_sorted_by_chromooffset, stat_file):
    print >> sys.stderr, "[%s] merge sort mapped reads canon" % (right_now())
    
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	
	temp_fs = open(merged_sorted_by_tagname, "w")
	temp_fs.close()
	
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    separatemapreads_cmd = [bin_dir + "merge_sort_sam", 
			temp_dir + "splicedreads_remdupdivided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_divided_reads.txt",		
			merged_sorted_by_tagname,
			merged_sorted_by_chromooffset,
			stat_file]
       
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort mapped reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sort_sam not found on this system"
        exit(1)
	
    temp_fs = open(merged_sorted_by_tagname, "w")
    temp_fs.close()
    return (merged_sorted_by_tagname, merged_sorted_by_chromooffset);

def merge_sort_sam_noncanon(merged_sorted_by_tagname, merged_sorted_by_chromooffset, stat_file):
    print >> sys.stderr, "[%s] merge sort mapped reads noncanon" % (right_now())
    
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	
	temp_fs = open(merged_sorted_by_tagname, "w")
	temp_fs.close()
	
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    separatemapreads_cmd = [bin_dir + "merge_sort_sam", 
			temp_dir + "fixed_fixhole_exceed_f0_divided_reads.txt",
			temp_dir + "fixed_fixhole_f0_divided_reads.txt",
			temp_dir + "fixed_hole_f0_divided_reads.txt",
			temp_dir + "fix_head_tail_noncanon_divided_reads.txt",
			temp_dir + "mapreads_noncanon_divided_reads.txt",			
			merged_sorted_by_tagname,
			merged_sorted_by_chromooffset,
			stat_file]
       
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort mapped reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sort_sam not found on this system"
        exit(1)
	
    temp_fs = open(merged_sorted_by_tagname, "w")
    temp_fs.close()
    return (merged_sorted_by_tagname, merged_sorted_by_chromooffset);

def pairendmappedreads_sep(bwtmapped_sam, bwt_25bp_mapped, fusion_sam, pairend, single, cur_output_dir, maxhits):
    print >> sys.stderr, "[%s] pairendmappedreads_sep reads" % (right_now())

    bwt_25bp_mapped_mishandt_matched = bwt_25bp_mapped + ".mistailandhead.matched";
    
    bwt_25bp_mapped_mishort_matched = bwt_25bp_mapped + ".mistailorhead.matched";
    
    bwt_25bp_mapped_allmapped = bwt_25bp_mapped + ".allmapped";
    
    resortmapreads_cmd = [bin_dir + "pairend", 
			temp_dir + "splicedreads_remdupdivided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_divided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_f0_divided_reads.txt",
			temp_dir + "fixed_fixhole_f0_divided_reads.txt",
			temp_dir + "fixed_hole_f0_divided_reads.txt",
			temp_dir + "fix_head_tail_noncanon_divided_reads.txt",
			temp_dir + "mapreads_noncanon_divided_reads.txt",
			bwtmapped_sam,
			bwt_25bp_mapped_mishandt_matched,
			bwt_25bp_mapped_mishort_matched,
			bwt_25bp_mapped_allmapped,
			fusion_sam,			
			output_dir_paired,
			str(maxhits),
                        output_dir_single]  
       
    if os.path.exists(output_dir_paired) and \
           os.path.exists(output_dir_single) and \
	   rerun_all == 0:
	return (output_dir_paired, output_dir_single)

    try:    
        retcode = subprocess.call(resortmapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pairendmappedreads_sep failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pairend not found on this system"
        exit(1)
    return (output_dir_paired, output_dir_single);

def compare2syntheticsam(mps_all_sam, synthe_sam, compared_out, comprange, all_junc, compared_offset_out, 
			 filter_sam, maxhits, uniquemapped, stat_file, anchor_width):
    print >> sys.stderr, "[%s] compare to synthetic sam" % (right_now())
    
    entrpy_weight = 0.097718;
    pqlen_weight = 0.66478;
    ave_mis_weight = -0.21077;

    separatemapreads_cmd = [bin_dir + "compare2synthesam", 
			mps_all_sam,
			synthe_sam,
			compared_out,
			str(comprange),
			all_junc,
			filter_sam,
			str(maxhits),
			uniquemapped,
			compared_offset_out,
			stat_file,
			str(anchor_width),
			str(entrpy_weight),
			str(pqlen_weight),
			str(ave_mis_weight)]
       
    if os.path.exists(filter_sam) and \
           os.path.exists(uniquemapped) and\
	   os.path.exists(stat_file) and \
	   rerun_all == 0:
	
	return (filter_sam, uniquemapped)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: compare to synthetic sam failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: compare2synthesam not found on this system"
        exit(1)
    return (filter_sam, uniquemapped);

def compare2sam(mps_all_sam, comp_sam, compared_out, comprange, all_junc, compared_offset_out, 
			 filter_sam, maxhits, uniquemapped, stat_file, anchor_width):
    print >> sys.stderr, "[%s] compare to another sam" % (right_now())
    
    separatemapreads_cmd = [bin_dir + "Compare2Sam", 
			mps_all_sam,
			comp_sam,
			compared_out,
			str(comprange),
			all_junc,
			filter_sam,
			str(maxhits),
			uniquemapped,
			compared_offset_out,
			stat_file,
			str(anchor_width)]
       
    if os.path.exists(filter_sam) and \
           os.path.exists(uniquemapped) and \
	   rerun_all == 0:
	return (filter_sam, uniquemapped)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: compare to another sam failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Compare2Sam not found on this system"
        exit(1)
    return (filter_sam, uniquemapped);

def sep_spliced_unspliced(merged_sam, unspliced_sam, spliced_sam):
    print >> sys.stderr, "[%s] separate sam by spliced unspliced" % (right_now())
    print >> sys.stdout, "[%s] separate sam by spliced unspliced" % (merged_sam)
    
    separatemapreads_cmd = [bin_dir + "separate_spliced_unspliced", 
			merged_sam,
			unspliced_sam,
			spliced_sam]
       
    if os.path.exists(unspliced_sam) and \
           os.path.exists(spliced_sam) and \
	   rerun_all == 0:
	
	temp_fs = open(unspliced_sam, "w")
	temp_fs.close()
	
	temp_fs = open(merged_sam, "w")
	temp_fs.close()
	
	return (unspliced_sam, spliced_sam)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate sam by spliced unspliced failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separate_spliced_unspliced not found on this system"
        exit(1)
	
    temp_fs = open(unspliced_sam, "w")
    temp_fs.close()
    
    temp_fs = open(merged_sam, "w")
    temp_fs.close()
	
    return (unspliced_sam, spliced_sam);

def filtersambyanchor(spliced_sam, anchor_width):
    print >> sys.stderr, "[%s] filter sam by anchor" % (right_now())
    print >> sys.stdout, "[%s] filter sam by anchor" % (spliced_sam)
    
    separatemapreads_cmd = [bin_dir + "filterbyanchor", 
			spliced_sam,
			str(anchor_width)]
       
    #if os.path.exists(spliced_sam + ".longanchor") and \
           #os.path.exists(spliced_sam + ".shortanchor"):
	#return (spliced_sam + ".longanchor", spliced_sam + ".shortanchor")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter sam by anchor failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyanchor not found on this system"
        exit(1)
	
    temp_fs = open(spliced_sam + ".longanchor", "w")
    temp_fs.close()
    
    return (spliced_sam + ".longanchor", spliced_sam + ".shortanchor");

def filterbyintronlenhmer(spliced_sam, intron_len, seglen):
    print >> sys.stderr, "[%s] filter by intronlen hmer" % (right_now())
    print >> sys.stdout, "[%s] filter by intronlen hmer" % (spliced_sam)
    
    separatemapreads_cmd = [bin_dir + "filterbyintronlenhmer", 
			spliced_sam,
			str(intron_len),
			str(seglen)]
       
    #if os.path.exists(spliced_sam + ".longanchor") and \
           #os.path.exists(spliced_sam + ".shortanchor"):
	#return (spliced_sam + ".longanchor", spliced_sam + ".shortanchor")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by intronlen hmer failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyintronlenhmer not found on this system"
        exit(1)
	
    temp_fs = open(spliced_sam + ".inintronlenhmer", "w")
    temp_fs.close()
    
    return (spliced_sam + ".inintronlenhmer", spliced_sam + ".notinintronlenhmer");

def filterbysmallexon(spliced_sam, seglen):
    print >> sys.stderr, "[%s] filter by small exon" % (right_now())
    print >> sys.stdout, "[%s] filter by small exon" % (spliced_sam)
    
    separatemapreads_cmd = [bin_dir + "filterbysmallexon", 
			spliced_sam,
			str(seglen)]
       
    #if os.path.exists(spliced_sam + ".longanchor") and \
           #os.path.exists(spliced_sam + ".shortanchor"):
	#return (spliced_sam + ".longanchor", spliced_sam + ".shortanchor")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by small exon failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbysmallexon not found on this system"
        exit(1)
	
    temp_fs = open(spliced_sam + ".insmallexon", "w")
    temp_fs.close()
    
    return (spliced_sam + ".insmallexon", spliced_sam + ".notinsmallexon");

def filternotinunmapped(spliced_sam, unmapped_reads, fa_fq):
    print >> sys.stderr, "[%s] filter spliced sam not in unmapped reads" % (right_now())
    print >> sys.stdout, "[%s] filter spliced sam not in unmapped reads" % (unmapped_reads)
    
    separatemapreads_cmd = [bin_dir + "filternotinunmapped", 
			spliced_sam,
			unmapped_reads,
			str(fa_fq)]
       
    #if os.path.exists(spliced_sam + ".notinunmapped") and \
           #os.path.exists(spliced_sam + ".inunmapped"):
	#return (spliced_sam + ".notinunmapped", spliced_sam + ".inunmapped")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter spliced sam not in unmapped reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filternotinunmapped not found on this system"
        exit(1)
	
    temp_fs = open(spliced_sam + ".inunmapped", "w")
    temp_fs.close()
    
    return (spliced_sam + ".notinunmapped", spliced_sam + ".inunmapped");

def filterbyrepeated_reads(spliced_sam, repeated_reads, fa_fq):
    print >> sys.stderr, "[%s] filter by repeated reads" % (right_now())
    print >> sys.stdout, "[%s] filter by repeated reads" % (repeated_reads)
    
    separatemapreads_cmd = [bin_dir + "filterbyrepeated_reads", 
			spliced_sam,
			repeated_reads,
			str(fa_fq)]
       
    #if os.path.exists(spliced_sam + ".notinunmapped") and \
           #os.path.exists(spliced_sam + ".inunmapped"):
	#return (spliced_sam + ".notinunmapped", spliced_sam + ".inunmapped")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by repeated reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyrepeated_reads not found on this system"
        exit(1)
	
    temp_fs = open(spliced_sam + ".notinrepeat", "w")
    temp_fs.close()
    
    return (spliced_sam + ".inrepeat", spliced_sam + ".notinrepeat");

def filterbyunsplicedmapped(spliced_sam, unspliced_mapped_sam):
    print >> sys.stderr, "[%s] filter by unsplicedmapped" % (right_now())
    print >> sys.stdout, "[%s] filter by unsplicedmapped" % (unspliced_mapped_sam)
    
    separatemapreads_cmd = [bin_dir + "filterbyunsplicedmapped", 
			spliced_sam,
			unspliced_mapped_sam]
       
    #if os.path.exists(spliced_sam + ".notinunmapped") and \
           #os.path.exists(spliced_sam + ".inunmapped"):
	#return (spliced_sam + ".notinunmapped", spliced_sam + ".inunmapped")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by unsplicedmapped failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyunsplicedmapped not found on this system"
        exit(1)
	
    temp_fs = open(spliced_sam + ".notinunsplicedmapped", "w")
    temp_fs.close()
    
    return (spliced_sam + ".notinunsplicedmapped", spliced_sam + ".inunsplicedmapped");

def filtersambyisland(spliced_sam, island_file):
    print >> sys.stderr, "[%s] filter sam by island" % (right_now())
    print >> sys.stdout, "[%s] filter sam by island" % (island_file)
    
    separatemapreads_cmd = [bin_dir + "filterbyisland", 
			spliced_sam,
			island_file]
       
    #if os.path.exists(spliced_sam + ".notinunmapped") and \
           #os.path.exists(spliced_sam + ".inunmapped"):
	#return (spliced_sam + ".notinunmapped", spliced_sam + ".inunmapped")
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter sam by island failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyisland not found on this system"
        exit(1)
    
    temp_fs = open(spliced_sam + ".inisland", "w")
    temp_fs.close()
    
    return (spliced_sam + ".inisland", spliced_sam + ".notinisland");

def FilterMultipleMapped(mps_all_sam, all_junc, filter_sam, maxhits, uniquemapped, stat_file, log_file):
    #print >> sys.stderr, "[%s] Filter multiple mapped" % (right_now())
    
    entrpy_weight = 0.097718;
    pqlen_weight = 0.66478;
    ave_mis_weight = -0.21077;
    
    separatemapreads_cmd = [bin_dir + "FilterMultipleMappedByRead", 
			mps_all_sam,
			all_junc,
			filter_sam,
			str(maxhits),
			uniquemapped,
			stat_file,
			str(entrpy_weight),
			str(pqlen_weight),
			str(ave_mis_weight)]
       
    if os.path.exists(filter_sam) and \
           os.path.exists(uniquemapped) and \
	   rerun_all == 0:
	return (filter_sam, uniquemapped)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    try:    
        retcode = subprocess.call(separatemapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filter multiple mapped failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterMultipleMapped not found on this system"
        exit(1)
    return (filter_sam, uniquemapped);

def merge_sort_sam2(merged_sorted_by_tagname, merged_sorted_by_chromooffset, mapped1, mapped2, stat_file):
    #print >> sys.stderr, "[%s] merge sort 2 mapped reads" % (right_now())
    
    separatemapreads_cmd = [bin_dir + "merge_sort_sam", 
			mapped1,
			mapped2,
			merged_sorted_by_tagname,
			merged_sorted_by_chromooffset,
			stat_file]
       
    if os.path.exists(merged_sorted_by_tagname) and \
           os.path.exists(merged_sorted_by_chromooffset) and \
	   rerun_all == 0:
	
	temp_fs = open(merged_sorted_by_tagname, "w")
	temp_fs.close()
	
	#temp_fs = open(mapped2, "w")
	#temp_fs.close()
	
	#temp_fs = open(mapped1, "w")
	#temp_fs.close()
    
	return (merged_sorted_by_tagname, merged_sorted_by_chromooffset)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort mapped reads failed"
            exit(1)
 
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sort_sam not found on this system"
        exit(1)
    
    temp_fs = open(merged_sorted_by_tagname, "w")
    temp_fs.close()
    
    #temp_fs = open(mapped2, "w")
    #temp_fs.close()
    
    #temp_fs = open(mapped1, "w")
    #temp_fs.close()
    
    return (merged_sorted_by_tagname, merged_sorted_by_chromooffset);

def pairsam(allsam, pairend, single, cur_output_dir, maxhits, log_file):
    print >> sys.stderr, "[%s] pairend mapped reads" % (right_now())
    
    output_dir_paired = cur_output_dir + pairend;

    output_dir_single = cur_output_dir + single;
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(output_dir_paired) and \
           os.path.exists(output_dir_single) and \
	   rerun_all == 0:
	return (output_dir_paired, output_dir_single)

    resortmapreads_cmd = [bin_dir + "pairend", 
		        allsam,
			str(maxhits),
			output_dir_paired,
                        output_dir_single]    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pairend mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pairend not found on this system"
        exit(1)
    return (output_dir_paired, output_dir_single);

def parseCluster(merge_pair_sam, output_dir, log_file):
    print >> sys.stderr, "[%s] parse cluster regions" % (right_now())
    
    resortmapreads_cmd = ["sed", 
                        "-e",
			"s/^/1~/",
                        merge_pair_sam]
    
    sed_log = open(fusion_dir + "sed.fusion.sam", "w")
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=sed_log)
	
	if retcode != 0:
		print >> sys.stderr, fail_str, "Error: sed failed"
		exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sed not found on this system"
        exit(1)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #if os.path.exists(output_file) and \
	   #rerun_all == 0:
	#return output_file

    resortmapreads_cmd = [bin_dir + "parseCluster", 
		        fusion_dir + "sed.fusion.sam",
			output_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: parse cluster regions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: parseCluster not found on this system"
        exit(1)

def cluster(cluster_dir, log_file):
    print >> sys.stderr, "[%s] cluster regions" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(cluster_dir + "result/cluster.txt") and \
	   rerun_all == 0:
	return cluster_dir + "result/cluster.txt"

    resortmapreads_cmd = [bin_dir + "cluster", 
		        cluster_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: cluster failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cluster not found on this system"
        exit(1)
    #return output_file;

def SepSam(sam_file, 
           unspliced,
           spliced,
           small_del,
           small_ins,
           clipped,
           unmapped,
           head,
           log_file):
    
    print >> sys.stderr, "[%s] Separate sam file" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #if os.path.exists(output_file) and \
	   #rerun_all == 0:
	#return output_file

    resortmapreads_cmd = [bin_dir + "SepSam", 
		        sam_file,
                        unspliced,
                        spliced,
                        small_del,
                        small_ins,
                        clipped,
                        unmapped,
                        head]
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Separate sam file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SepSam not found on this system"
        exit(1)
	
def SepSamMappedUnmapped(sam_file, 
                         head,
                         mapped,
                         fusion,
                         unmapped,
                         log_file):
    
    print >> sys.stderr, "[%s] Separate sam file to header mapped and unmapped" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #if os.path.exists(output_file) and \
	   #rerun_all == 0:
	#return output_file

    resortmapreads_cmd = [bin_dir + "SepSamUnmapped", 
                          mapped,
                          fusion,
                          unmapped,
                          head,
                          sam_file]
	  
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    if os.path.exists(mapped) and \
       os.path.exists(fusion) and \
       os.path.exists(unmapped) and \
       os.path.exists(head) and \
               rerun_all == 0:
	return mapped

    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Separate sam file to header mapped and unmapped"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SepSamUnmapped not found on this system"
        exit(1)

def recover_fusion_alignments_order(fusion_original, 
                                    fusion_converted,
                                    log_file):
    
    print >> sys.stderr, "[%s] Convert fusion two lines alignment into one line" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #if os.path.exists(output_file) and \
	   #rerun_all == 0:
	#return output_file

    resortmapreads_cmd = [bin_dir + "recover_fusion_alignments_order", 
                          fusion_original,
                          fusion_converted]
	  
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    if os.path.exists(fusion_converted) and \
               rerun_all == 0:
	    return fusion_converted    

    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Convert fusion two lines alignment into one line"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: recover_fusion_alignments_order not found on this system"
        exit(1)
	
def find_mate_sam_fq(unmapped_sam_file, input_reads, unmapped_fq, log_file):
    print >> sys.stderr, "[%s] find unmapped sam" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #if os.path.exists(output_file) and \
	   #rerun_all == 0:
	#return output_file

    input_read_files = input_reads.split(',')	
    
    resortmapreads_cmd = [bin_dir + "find_mate_sam_fq", 
		        unmapped_sam_file,
                        unmapped_fq]
    
    for input_read_file in input_read_files:
	resortmapreads_cmd = resortmapreads_cmd + [input_read_file]
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: cluster failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cluster not found on this system"
        exit(1)
	
def ReadRegions(region_file, single_sam, reads_file, out_path, format_flag, chrom_dir, output_dir_file, extend_len, log_file):
    print >> sys.stderr, "[%s] read regions" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(output_dir_file) and \
	   rerun_all == 0:
	return output_dir_file
    
    format_int = 0;
    
    if format_flag == "-q":
	format_int = 1;	

    resortmapreads_cmd = [bin_dir + "ReadRegions", 
		        region_file,
			single_sam,
			reads_file,
			out_path,
			str(format_int),
			chrom_dir,
			output_dir_file,
			str(extend_len)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: read regions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: ReadRegions not found on this system"
        exit(1)
    return output_dir_file;

def FilterFusionCandidatesByClusterRegions(region_file, single_sam, fusion_candidate_file, extend_len, seg_num, log_file):
    print >> sys.stderr, "[%s] filter fusion candidate by regions" % (right_now())

#/stage/wk0571/syc/FilterFusionCandidatesByClusterRegions/FilterFusionCandidatesByClusterRegions/FilterFusionCandidatesByClusterRegions 
#../../cluster/result/cluster.txt ../../remap/remapped_single.sam sorted_combine_head_tail.txt 10000 2    

    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(fusion_candidate_file + ".remained.formated") and \
	   rerun_all == 0:
	return fusion_candidate_file + ".remained.formated"

    resortmapreads_cmd = [bin_dir + "FilterFusionCandidatesByClusterRegions", 
		        region_file,
			single_sam,
			fusion_candidate_file,
			str(extend_len),
			str(seg_num)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter fusion candidate by regions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterFusionCandidatesByClusterRegions not found on this system"
        exit(1)
    return fusion_candidate_file + ".remained.formated";

def pairing_fusion_normal_aligned(region_file, single_sam, fusion_candidate_file, extend_len, log_file):
    print >> sys.stderr, "[%s] pairing fusion normal aligned" % (right_now())

#/stage/wk0571/syc/FilterFusionCandidatesByClusterRegions/FilterFusionCandidatesByClusterRegions/FilterFusionCandidatesByClusterRegions 
#../../cluster/result/cluster.txt ../../remap/remapped_single.sam sorted_combine_head_tail.txt 10000 2    

    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(fusion_candidate_file + ".paired") and \
       os.path.exists(fusion_candidate_file + ".fusion_single") and \
       os.path.exists(fusion_candidate_file + ".normal_single") and \
	   rerun_all == 0:
	return fusion_candidate_file + ".paired"

    resortmapreads_cmd = [bin_dir + "pairing_fusion_normal_aligned", 
		        region_file,
			single_sam,
			fusion_candidate_file,
			str(extend_len)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pairing fusion normal aligned failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pairing_fusion_normal_aligned not found on this system"
        exit(1)
    return fusion_candidate_file + ".paired";

def load_fusion_single_anchor_chrom_seq(fusion_candidate_file, fusion_candidate_file_loaded, chrom_dir, log_file):
    print >> sys.stderr, "[%s] load_fusion_single_anchor_chrom_seq" % (right_now())

    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(fusion_candidate_file_loaded) and \
	   rerun_all == 0:
	return fusion_candidate_file_loaded

    resortmapreads_cmd = [bin_dir + "load_fusion_single_anchor_chrom_seq", 
		       fusion_candidate_file,
			fusion_candidate_file_loaded,
			chrom_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: load_fusion_single_anchor_chrom_seq failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: load_fusion_single_anchor_chrom_seq not found on this system"
        exit(1)
    return fusion_candidate_file_loaded;

def generate_bash_file_and_run(director_file, bash_file, abs_path, arguments, bin_dir, log_file):
    print >> sys.stderr, "[%s] generate bash file and run" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #if os.path.exists(bash_file) and \
	   #rerun_all == 0:
	#return bash_file
    quote_bin_dir = "\"" + bin_dir + "\"";
    
    resortmapreads_cmd = [bin_dir + "generate_bash_file_and_run", 
		        director_file,
			bash_file,
			abs_path,
			arguments,
			quote_bin_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: generate bash file and run failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: generate_bash_file_and_run not found on this system"
        exit(1)
	
    assign_exec_cmd = ["chmod",
		       "700",
		       bash_file]
    
    print >> sys.stderr, "[%s] " % assign_exec_cmd
    
    try:    
        retcode = subprocess.call(assign_exec_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: assign executable bash file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: chmod not found on this system"
        exit(1)
    
    print >> sys.stderr, "[%s] execute bash file" % (right_now())
    
    try:    
        #retcode = subprocess.call(bash_file)
	os.popen(bash_file)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: execute bash file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bash file not found on this system"
        exit(1)
	
    return bash_file;

def convert_to_abs_offset(director_file, sam_output_file, fusion_output_file, abs_path, log_file):
    print >> sys.stderr, "[%s] convert to abs offset" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if os.path.exists(sam_output_file) and \
       os.path.exists(fusion_output_file) and \
	   rerun_all == 0:
	return (sam_output_file, fusion_output_file)
    
    resortmapreads_cmd = [bin_dir + "convert_to_abs_offset", 
		        director_file,
			sam_output_file,
			fusion_output_file,
			abs_path]
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert to abs offset failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: generate_bash_file not found on this system"
        exit(1)
	
    return (sam_output_file, fusion_output_file);

def parsePER(PERsam, log_file):
    print >> sys.stderr, "[%s] parse PER sam file" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    test_file = fusion_data_PER_dir + "chrXchrX.txt"
    
    if os.path.exists(test_file) and \
           rerun_all == 0:
	return

    resortmapreads_cmd = [bin_dir + "parsePER", 
		        PERsam,
			fusion_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: parse PER sam file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: parsePER not found on this system"
        exit(1)

def parseSingle(PERsam, log_file):
    print >> sys.stderr, "[%s] parse Single sam file" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    test_file = fusion_data_single_dir + "chrX.txt"
    
    if os.path.exists(test_file) and \
           rerun_all == 0:
	return

    resortmapreads_cmd = [bin_dir + "parseSingle", 
		        PERsam,
			fusion_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: parse Single sam file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: parseSingle not found on this system"
        exit(1)

def PERall(log_file):
    print >> sys.stderr, "[%s] Map PER reads" % (right_now())
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    test_file = fusion_result_PER_prob_dir + "PERprob_chrXchrX.sam"
    
    if os.path.exists(test_file) and \
           rerun_all == 0:
	return

    resortmapreads_cmd = [bin_dir + "MapPERall", 
		        "all",
			bin_dir,
			fusion_dir]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % resortmapreads_cmd
    
    try:    
        retcode = subprocess.call(resortmapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Map PER reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MapPERall not found on this system"
        exit(1)
	
def runPER(PERsam, Singlesam):
    print >> sys.stderr, "[%s] processing PER reads" % (right_now())

    parsePER(PERsam, logging_dir + "PERsam.log")
    
    parseSingle(Singlesam, logging_dir + "Singlesam.log")
    
    PERall(logging_dir + "PERall.log");
	
def pairendmappedreads(allsam, fusion_sam, pairend, single, cur_output_dir, maxhits):
    print >> sys.stderr, "[%s] pairend mapped reads" % (right_now())
    
    output_dir_paired = cur_output_dir + pairend;

    output_dir_single = cur_output_dir + single;
    
    if os.path.exists(output_dir_paired) and \
           os.path.exists(output_dir_single) and \
	   rerun_all == 0:
	return (output_dir_paired, output_dir_single)

    resortmapreads_cmd = [bin_dir + "pairend", 
		        allsam,
			fusion_sam,
			output_dir_paired,
			str(maxhits),
                        output_dir_single]    
    try:    
        retcode = subprocess.call(resortmapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pairend mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pairend not found on this system"
        exit(1)
    return (output_dir_paired, output_dir_single);

def pairendmappedreads3(sam1, sam2, sam3, pairend, single, cur_output_dir, maxhits):
    print >> sys.stderr, "[%s] pairend mapped reads" % (right_now())
    
    output_dir_paired = cur_output_dir + pairend;

    output_dir_single = cur_output_dir + single;
    
    if os.path.exists(output_dir_paired) and \
           os.path.exists(output_dir_single) and \
	   rerun_all == 0:
	return (output_dir_paired, output_dir_single)

    resortmapreads_cmd = [bin_dir + "pairend", 
		        sam1,
			sam2,
			sam3,
			output_dir_paired,
			str(maxhits),
                        output_dir_single]    
    try:    
        retcode = subprocess.call(resortmapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: pairend mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pairend not found on this system"
        exit(1)
    return (output_dir_paired, output_dir_single);

def separatedmultipleunique_non_canon(mps_unique_mapreads, mps_multiple_mapreads, cur_output_dir):
    print >> sys.stderr, "[%s] separate unique and multiple mapped reads" % (right_now())

    separatemapreads_cmd = [bin_dir + "separatemapreads", 
		        cur_output_dir,
			"fixed_fixhole_f0_divided_reads.txt",
			"fixed_hole_f0_divided_reads.txt",
			"mapreads_noncanon_divided_reads.txt",
			"fix_head_tail_noncanon_divided_reads.txt",
			mps_unique_mapreads,
			mps_multiple_mapreads,
                        "1"]
       
    if os.path.exists(cur_output_dir + mps_unique_mapreads) and \
           os.path.exists(cur_output_dir + mps_multiple_mapreads) and \
	   rerun_all == 0:
	return mps_unique_mapreads
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separated mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separatemapreads not found on this system"
        exit(1)
    return (mps_unique_mapreads, mps_multiple_mapreads);

def separatedmultipleunique1(combinefile, unique_mapreads, multiple_mapreads, cur_output_dir):
    print >> sys.stderr, "[%s] separate combined to unique and multiple mapped reads" % (right_now())

    separatemapreads_cmd = [bin_dir + "separatemapreads", 
		        cur_output_dir,
			combinefile,
			unique_mapreads,
			multiple_mapreads,
                        "1"]
       
    if os.path.exists(cur_output_dir + unique_mapreads) and \
           os.path.exists(cur_output_dir + multiple_mapreads) and \
	   rerun_all == 0:
	return unique_mapreads
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separated mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separatemapreads not found on this system"
        exit(1)
    return (unique_mapreads, multiple_mapreads);

def separatedmultipleuniquefusion(combinefile, unique_mapreads, multiple_mapreads, stat_file):
    print >> sys.stderr, "[%s] separate combined fusion to unique and multiple mapped reads" % (right_now())

    separatemapreads_cmd = [bin_dir + "separateuniquefusion", 
			combinefile,
			unique_mapreads,
			multiple_mapreads,
			stat_file]
       
    if os.path.exists(unique_mapreads) and \
           os.path.exists(multiple_mapreads) and \
	rerun_all == 0:
	return unique_mapreads
    
    mapsplice_log = open(logging_dir + "separatedmultipleuniquefusion.log", "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(separatemapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separated mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separatemapreads not found on this system"
        exit(1)
    return (unique_mapreads, multiple_mapreads);

def separateuniquefusion_newfmt(combinefile, unique_mapreads, multiple_mapreads, stat_file):
    print >> sys.stderr, "[%s] separate combined new fmt fusion to unique and multiple mapped reads" % (right_now())

    separatemapreads_cmd = [bin_dir + "separateuniquefusion_newfmt", 
			combinefile,
			unique_mapreads,
			multiple_mapreads,
			stat_file]
       
    if os.path.exists(unique_mapreads) and \
           os.path.exists(multiple_mapreads) and \
	rerun_all == 0:
	return unique_mapreads
    
    mapsplice_log = open(logging_dir + "separateuniquefusion_newfmt.log", "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(separatemapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate combined new fmt fusion to unique and multiple mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separateuniquefusion_newfmt not found on this system"
        exit(1)
    return (unique_mapreads, multiple_mapreads);

def separate_canon_noncanon(combinejunc, canon_junc, noncanon_junc):
    print >> sys.stderr, "[%s] separate combined junction to canon and noncanon" % (right_now())

    separatemapreads_cmd = [bin_dir + "separate_canon_noncanon", 
			combinejunc,
			canon_junc,
			noncanon_junc]
       
    if os.path.exists(canon_junc) and \
           os.path.exists(noncanon_junc) and \
	   rerun_all == 0:
	return (canon_junc, noncanon_junc)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separate combined junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separate_canon_noncanon not found on this system"
        exit(1)
    return (canon_junc, noncanon_junc);

def count_canon_noncanon(combinejunc, canon, noncanon, log_file):
    #print >> sys.stderr, "[%s] count combined junction to canon and noncanon" % (right_now())
    
    separatemapreads_cmd = [bin_dir + "count_canon_noncanon", 
			combinejunc, canon, noncanon]
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(separatemapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: count combined junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: count_canon_noncanon not found on this system"
        exit(1)
	
    return (canon, noncanon)
	
def filterjuncbyminmis(combinejunc, min_mismatch, in_min_mis, not_in_min_mis):
    #print >> sys.stderr, "[%s] filter junction by min mismatch" % (right_now())
    
    separatemapreads_cmd = [bin_dir + "filterbyjuncminmis", 
			combinejunc, str(min_mismatch), in_min_mis, not_in_min_mis]
       
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junction by min mismatch failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyjuncminmis not found on this system"
        exit(1)
	
    return (in_min_mis, not_in_min_mis)

def filterjuncbysmalldeletion(combinejunc, min_intron, in_small_deletion, not_in_small_deletion, log_file):
    #print >> sys.stderr, "[%s] filter junction by small deletion" % (right_now())
    
    separatemapreads_cmd = [bin_dir + "filterjuncbysmalldeletion", 
			combinejunc, str(min_intron), in_small_deletion, not_in_small_deletion]
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    try:    
        retcode = subprocess.call(separatemapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junction by small deletion failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterjuncbysmalldeletion not found on this system"
        exit(1)
	
    return (in_small_deletion, not_in_small_deletion)

def filterfusionjuncbyminmis(combinejunc, min_mis, in_min_mis, not_in_min_mis, log_file):
    print >> sys.stderr, "[%s] filter fusion junction by min mismatch" % (right_now())
    
    separatemapreads_cmd = [bin_dir + "filterfusionjuncbyminmis", 
			combinejunc, str(min_mis), in_min_mis, not_in_min_mis]
    
    #mapsplice_log = open(logging_dir + "filterfusionjuncbyminmis.log", "w")
    ##stdout=mapsplice_log
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(separatemapreads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error:  filter fusion junction by min mismatch failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterfusionjuncbyminmis not found on this system"
        exit(1)
	
    return (in_min_mis, not_in_min_mis)

def filtermappedreads(samfile, notfiltered, filtered, objreads, opt):
    print >> sys.stderr, "[%s] filter mapped reads" % (right_now())

    separatemapreads_cmd = [bin_dir + "filtermappedreads",
			    opt,
			    objreads,
			    notfiltered,
			    filtered,
			    samfile]
       
    if os.path.exists(filtered) and \
           os.path.exists(notfiltered) and \
	   rerun_all == 0:
	return (filtered, notfiltered)
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filtermappedreads not found on this system"
        exit(1)
    return (filtered, notfiltered);

def separatedfullymapped(bwt_mapped_sam, bwt_25bp_mapped, cur_output_dir):
    print >> sys.stderr, "[%s] separate fully mapped unique and multiple mapped reads" % (right_now())

    bwt_25bp_mapped_mishandt_matched = bwt_25bp_mapped + ".mistailandhead.matched";
    
    bwt_25bp_mapped_mishort_matched = bwt_25bp_mapped + ".mistailorhead.matched";
    bwt_25bp_mapped_allmapped = bwt_25bp_mapped + ".allmapped";
    
    if os.path.exists(cur_output_dir + "fully_unique.txt") and \
       os.path.exists(cur_output_dir + "fully_multiple.txt") and \
	   rerun_all == 0:
	return "fully_unique.txt"
    separatemapreads_cmd = [bin_dir + "separatemapreads", 
		        cur_output_dir,
			bwt_mapped_sam,
			bwt_25bp_mapped_mishandt_matched,
			bwt_25bp_mapped_mishort_matched,
			bwt_25bp_mapped_allmapped,
			"fully_unique.txt",
			"fully_multiple.txt",
                        "0"]    
    try:    
        retcode = subprocess.call(separatemapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: separated mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: separatemapreads not found on this system"
        exit(1)
    return ("fully_unique.txt", "fully_multiple.txt");

def countline(countfile, log_file):
    #print >> sys.stderr, "[%s] count line" % (right_now())

    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    countline_cmd = [bin_dir + "countline", 
		        countfile]    
    try:    
        retcode = subprocess.call(countline_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: count line failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: countline not found on this system"
        exit(1)
	
def write_current_stats(cur_stats):
    #print >> sys.stderr, "[%s] count line" % (right_now())

    mapsplice_log = open(logging_dir + "current_stats.log", "w")

    print >> mapsplice_log, cur_stats

    
def filtermultiplemapped(mps_multiple_mapreads, mps_filtered_multiple, junctions, maxhits, cur_output_dir):
    print >> sys.stderr, "[%s] filter multiple mapped reads" % (right_now())

    if os.path.exists(output_dir + mps_filtered_multiple) and \
	   rerun_all == 0:
	return mps_filtered_multiple
    
    filtermapreads_cmd = [bin_dir + "filtermapreads", 
		        cur_output_dir,
			mps_multiple_mapreads,
			str(maxhits),
			junctions,
			mps_filtered_multiple]    
    try:    
        retcode = subprocess.call(filtermapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter multiple mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filtermapreads not found on this system"
        exit(1)
    return mps_filtered_multiple;

def sam2junc_all(converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, min_anchor_length):
    print >> sys.stderr, "[%s] convert sam file to junctions" % (right_now())
    
    if os.path.exists(converted_junc) and \
	   rerun_all == 0:
	return converted_junc
    
   

    sam2junc_cmd = [bin_dir + "newsam2junc", 
		        converted_junc,
			str(read_width),
			chromosomes_file_dir,
			str(min_intron_length),str(max_intron_length),
			str(min_anchor_length),
			temp_dir + "splicedreads_remdupdivided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_divided_reads.txt",
			temp_dir + "fixed_fixhole_exceed_f0_divided_reads.txt",
			temp_dir + "fixed_fixhole_f0_divided_reads.txt",
			temp_dir + "fixed_hole_f0_divided_reads.txt",
			temp_dir + "fix_head_tail_noncanon_divided_reads.txt",
			temp_dir + "mapreads_noncanon_divided_reads.txt"]    
    try:    
        retcode = subprocess.call(sam2junc_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2junc not found on this system"
        exit(1)
    return converted_junc;

def sam2junc4(mapped1, mapped2, mapped3, mapped4, converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, min_anchor_length):
    print >> sys.stderr, "[%s] convert sam file to junctions" % (right_now())
    
    if os.path.exists(converted_junc) and \
	   rerun_all == 0:
	return converted_junc

    sam2junc_cmd = [bin_dir + "newsam2junc", 
		        converted_junc,
			str(read_width),
			chromosomes_file_dir,
			str(min_intron_length),str(max_intron_length),
			str(min_anchor_length),
			mapped1, 
			mapped2,
			mapped3,
			mapped4]    
    try:    
        retcode = subprocess.call(sam2junc_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2junc not found on this system"
        exit(1)
    return converted_junc;

def sam2junc(mps_unique_mapped, mps_multiple_mapped, converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, cur_output_dir, min_anchor_length):
    print >> sys.stderr, "[%s] convert sam file to junctions" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_mps_multiple_mapped = cur_output_dir + mps_multiple_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	return output_dir_converted_junc

    sam2junc_cmd = [bin_dir + "newsam2junc", 
		        output_dir_converted_junc,
			str(read_width),
			chromosomes_file_dir,
			str(min_intron_length),str(max_intron_length),
			str(min_anchor_length),
			output_dir_mps_unique_mapped, 
			output_dir_mps_multiple_mapped]    
    try:    
        retcode = subprocess.call(sam2junc_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2junc not found on this system"
        exit(1)
    return output_dir_converted_junc;

def sam2junc1(mps_unique_mapped, converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, cur_output_dir, min_anchor_length, log_file):
    print >> sys.stderr, "[%s] convert sam file to junctions" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;

    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	
	if os.path.exists(output_dir_converted_junc + ".sepdir"):
	    shutil.rmtree(output_dir_converted_junc + ".sepdir")

	return output_dir_converted_junc
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    sam2junc_cmd = [bin_dir + "newsam2junc", 
		        output_dir_converted_junc,
			str(read_width),
			chromosomes_file_dir,
			str(min_intron_length),str(max_intron_length),
			str(min_anchor_length),
			output_dir_mps_unique_mapped]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2junc not found on this system"
        exit(1)
	
    if os.path.exists(output_dir_converted_junc + ".sepdir"):
	    shutil.rmtree(output_dir_converted_junc + ".sepdir")
	    
    return output_dir_converted_junc;

def junc2bed(junction_file, junction_bed, log_file):
    print >> sys.stderr, "[%s] convert junction file to bed format" % (right_now())

    if os.path.exists(junction_bed) and \
	   rerun_all == 0:
	return junction_bed
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    sam2junc_cmd = [bin_dir + "junc2bed", 
		        junction_file,
			junction_bed]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert junction file to bed format failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: junc2bed not found on this system"
        exit(1)
    return junction_bed;

def junc2bed2(junction_file1, junction_file2, junction_bed, log_file):
    print >> sys.stderr, "[%s] convert junction file to bed format" % (right_now())

    if os.path.exists(junction_bed) and \
	   rerun_all == 0:
	return junction_bed
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    sam2junc_cmd = [bin_dir + "junc2bed", 
		        junction_file1,
			junction_file2,
			junction_bed]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert junction file to bed format failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: junc2bed not found on this system"
        exit(1)
    return junction_bed;

def sam2juncarray(sam_file_array, converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, cur_output_dir, min_anchor_length, log_file):
    print >> sys.stderr, "[%s] convert sam files to junctions" % (right_now())

    #output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    sam2junc_cmd = [bin_dir + "newsam2junc", 
		    output_dir_converted_junc,
		    str(read_width),
		    chromosomes_file_dir,
		    str(min_intron_length),str(max_intron_length),
		    str(min_anchor_length)]
    
    for sam_file in sam_file_array:
	sam2junc_cmd = sam2junc_cmd + [sam_file]
	
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
   
    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	if os.path.exists(output_dir_converted_junc + ".sepdir"):
	    shutil.rmtree(output_dir_converted_junc + ".sepdir")
	    
	return output_dir_converted_junc

    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % sam2junc_cmd
    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2junc not found on this system"
        exit(1)
	
    if os.path.exists(output_dir_converted_junc + ".sepdir"):
	    shutil.rmtree(output_dir_converted_junc + ".sepdir")
	    
    return output_dir_converted_junc;

def sam2juncarray_paired(sam_file_array, converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, cur_output_dir, min_anchor_length, paired, log_file):
    print >> sys.stderr, "[%s] convert array paired sam file to junctions" % (right_now())

    #output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    sam2junc_cmd = [bin_dir + "newsam2junc_paired", 
		    output_dir_converted_junc,
		    str(read_width),
		    chromosomes_file_dir,
		    str(min_intron_length),str(max_intron_length),
		    str(min_anchor_length),
		    paired]
    
    for sam_file in sam_file_array:
	sam2junc_cmd = sam2junc_cmd + [sam_file]
	
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	if os.path.exists(output_dir_converted_junc + ".sepdir"):
	    shutil.rmtree(output_dir_converted_junc + ".sepdir")
	    
	return output_dir_converted_junc

    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert paired sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: newsam2junc_paired not found on this system"
        exit(1)
	
    if os.path.exists(output_dir_converted_junc + ".sepdir"):
	    shutil.rmtree(output_dir_converted_junc + ".sepdir")
	    
    return output_dir_converted_junc;

def bed2bigbed(junction_bed, junction_bigbed, all_chromosomes_file):
    print >> sys.stderr, "[%s] convert junction bed to big bed format" % (right_now())

    if os.path.exists(junction_bigbed) and \
	   rerun_all == 0:
	return junction_bigbed
    
    rm_1stline_cmd = [bin_dir + "remove_firstline", 
		        junction_bed,
			junction_bed + ".rm_1st_line"]    
    try:    
        retcode = subprocess.call(rm_1stline_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: remove first line failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: remove_firstline not found on this system"
        exit(1)
	
    #bedToBigBed input.bed chrom.sizes myBigBed.bb
 	
    bedToBigBed_cmd = [bin_dir + "bedToBigBed", 
		        junction_bed + ".rm_1st_line",
			all_chromosomes_file + ".chrom.sizes",
			junction_bigbed]    
    try:    
        retcode = subprocess.call(bedToBigBed_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: bed to big bed failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bedToBigBed not found on this system"
        exit(1)
	
    #add_1st_line_cmd = ["./add_firstline", 
		        #junction_bigbed + ".notrack",
			#junction_bigbed,
			#"track type=bigBed name=\"big bed junctions\" description=\"Mapsplice junctions\" bigDataUrl=http://bioinfo.cs.uky.edu/projects/kai/mapsplice/best_syn_junction.bb"]    
    #try:    
        #retcode = subprocess.call(add_1st_line_cmd)
       
        ## cvg_islands returned an error 
        #if retcode != 0:
            #print >> sys.stderr, fail_str, "Error: add 1st line failed"
            #exit(1)
    ## cvg_islands not found
    #except OSError, o:
        #if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            #print >> sys.stderr, fail_str, "Error: add_firstline not found on this system"
        #exit(1)
    temp_fs = open(junction_bed + ".rm_1st_line", "w")
    temp_fs.close()

    return junction_bigbed;

def wig2bigwig(coverage_wig, coverage_big_wig, all_chromosomes_file):
    print >> sys.stderr, "[%s] convert coverage wig to big wig format" % (right_now())

    if os.path.exists(coverage_big_wig) and \
	   rerun_all == 0:
	return coverage_big_wig
    
    rm_1stline_cmd = [bin_dir + "remove_firstline", 
		        coverage_wig,
			coverage_wig + ".rm_1st_line"]    
    try:    
        retcode = subprocess.call(rm_1stline_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: remove first line failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: remove_firstline not found on this system"
        exit(1)
	
    #bedToBigBed input.bed chrom.sizes myBigBed.bb
 	
    wigToBigWig_cmd = [bin_dir + "wigToBigWig", 
		        coverage_wig + ".rm_1st_line",
			all_chromosomes_file + ".chrom.sizes",
			coverage_big_wig]    
    try:    
        retcode = subprocess.call(wigToBigWig_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: wig to big wig failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: wigToBigWig not found on this system"
        exit(1)
	
    #add_1st_line_cmd = ["./add_firstline", 
		        #coverage_big_wig + ".notrack",
			#coverage_big_wig,
			#"track type=bigWig name=\"big wig coverage\" description=\"Mapsplice coverage\" bigDataUrl=http://bioinfo.cs.uky.edu/projects/kai/mapsplice/coverage.bw"]    
    #try:    
        #retcode = subprocess.call(add_1st_line_cmd)
       
        ## cvg_islands returned an error 
        #if retcode != 0:
            #print >> sys.stderr, fail_str, "Error: add 1st line failed"
            #exit(1)
    ## cvg_islands not found
    #except OSError, o:
        #if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            #print >> sys.stderr, fail_str, "Error: add_firstline not found on this system"
        #exit(1)

    temp_fs = open(coverage_wig + ".rm_1st_line", "w")
    temp_fs.close()
    
    return coverage_big_wig;

def filter_junc_by_min_mis_lpq(junction_file, remained_junc, filtered_out_junc, min_mismatch, min_lpq, log_file):
    print >> sys.stderr, "[%s] filter junction by min mis and min lpq" % (right_now())

    if os.path.exists(filtered_out_junc) and \
       os.path.exists(remained_junc) and \
	   rerun_all == 0:
	return (remained_junc, filtered_out_junc)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if min_lpq > 0.5:
        min_lpq = 0.5
    
    sam2junc_cmd = [bin_dir + "filter_1hits",
		        junction_file,
			str(min_mismatch),
			str(min_lpq),
			remained_junc,
			filtered_out_junc]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % sam2junc_cmd
    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junction by min mis and min lpq failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filter_1hits not found on this system"
        exit(1)
    return (remained_junc, filtered_out_junc);

def filteroriginalfusion(junction_file, remained_junc, filtered_out_junc, min_mismatch, min_lpq, log_file):
    print >> sys.stderr, "[%s] filter original fusion junction" % (right_now())

    if os.path.exists(filtered_out_junc) and \
       os.path.exists(remained_junc) and \
	   rerun_all == 0:
	return (remained_junc, filtered_out_junc)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if min_lpq > 0.5:
        min_lpq = 0.5
    
    sam2junc_cmd = [bin_dir + "filteroriginalfusion",
		        junction_file,
			str(min_mismatch),
			str(min_lpq),
			remained_junc,
			filtered_out_junc]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % sam2junc_cmd
    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
	if retcode == 100:
            print >> sys.stderr, "Waring: No original fusions found, skip build index step"
	    #exit(1)
            return ("", "");
	
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junction by min mis and min lpq failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filter_1hits not found on this system"
        exit(1)
    return (remained_junc, filtered_out_junc);

def filterjuncbyROCarguNoncanon(junction_file, remained_junc, filtered_out_junc, entropy_weight, lpq_weight, ave_mis_weight, min_score, log_file):
    print >> sys.stderr, "[%s] filter junc by ROC argu noncanonical" % (right_now())

    #if os.path.exists(filtered_out_junc) and \
       #os.path.exists(remained_junc) and \
	#rerun_all == 0:
	#return (remained_junc, filtered_out_junc)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    sam2junc_cmd = [bin_dir + "filterjuncbyROCarguNonCanonical",
		        junction_file,
			str(entropy_weight),
			str(lpq_weight),
			str(ave_mis_weight),
			'0',
			'0',
			str(min_score),
			'5',
			remained_junc,
			filtered_out_junc]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
	if retcode == 100:
            print >> sys.stderr, "Waring: No original junctions found, skip build index step"
	    #exit(1)
            return ("", "");
	    
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junc by ROC argu failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterjuncbyROCargu not found on this system"
        exit(1)
    return (remained_junc, filtered_out_junc);

def filterjuncbyROCargu(junction_file, remained_junc, filtered_out_junc, entropy_weight, lpq_weight, ave_mis_weight, min_score, log_file):
    print >> sys.stderr, "[%s] filter junc by ROC argu" % (right_now())

    if os.path.exists(filtered_out_junc) and \
       os.path.exists(remained_junc) and \
	   rerun_all == 0:
	return (remained_junc, filtered_out_junc)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    sam2junc_cmd = [bin_dir + "filterjuncbyROCargu",
		        junction_file,
			str(entropy_weight),
			str(lpq_weight),
			str(ave_mis_weight),
			str(min_score),
			remained_junc,
			filtered_out_junc]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter junc by ROC argu failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterjuncbyROCargu not found on this system"
        exit(1)
    return (remained_junc, filtered_out_junc);

def SamHandler_sam2junc1(mps_unique_mapped, converted_junc, chromosomes_file_dir, 
	     read_width, min_intron_length, max_intron_length, cur_output_dir):
    print >> sys.stderr, "[%s] SamHandler convert sam file to junctions" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;

    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	return output_dir_converted_junc
    
    sam2junc_cmd = [bin_dir + "SamHandlerSam2junc", 
		        output_dir_converted_junc,
			str(read_width),
			chromosomes_file_dir,
			str(min_intron_length),str(max_intron_length),
			output_dir_mps_unique_mapped]    
    try:    
        retcode = subprocess.call(sam2junc_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SamHandlerSam2junc not found on this system"
        exit(1)
    return output_dir_converted_junc;

def fusionsam2junc1(mps_unique_mapped, converted_junc, read_width, cur_output_dir, log_file):
    print >> sys.stderr, "[%s] convert fusion sam file to junctions" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	return output_dir_converted_junc
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "fusionsam2junc", 
		        output_dir_converted_junc,
			str(read_width),
			output_dir_mps_unique_mapped]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fusionsam2junc not found on this system"
        exit(1)
    return output_dir_converted_junc;



def fusionsam2juncplus1(mps_unique_mapped, converted_junc, read_width, cur_output_dir, log_file):
    print >> sys.stderr, "[%s] convert fusion sam file to junctions" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	return output_dir_converted_junc
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "fusionsam2juncplus1", 
		        output_dir_converted_junc,
			str(read_width),
			output_dir_mps_unique_mapped]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fusionsam2junc not found on this system"
        exit(1)
    return output_dir_converted_junc;

def fusionsam2junc_filteranchor(mps_unique_mapped, converted_junc, read_width, cur_output_dir, min_anchor, log_file):
    print >> sys.stderr, "[%s] convert fusion sam file to junctions and filter by anchor" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	return output_dir_converted_junc
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "fusionsam2junc_filteranchor", 
		        output_dir_converted_junc,
			str(read_width),
			str(min_anchor),
			output_dir_mps_unique_mapped]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam file to junctions and filter by anchor failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fusionsam2junc_filteranchor not found on this system"
        exit(1)
    return output_dir_converted_junc;

def fusionsam2junc_filteranchor_newfmt(mps_unique_mapped, converted_junc, 
                                       read_width, cur_output_dir, min_anchor, 
                                       chromosomes_file_dir, log_file):
    print >> sys.stderr, "[%s] convert new fusion sam file to junctions and filter by anchor" % (right_now())

    output_dir_mps_unique_mapped = cur_output_dir + mps_unique_mapped;
    
    output_dir_converted_junc = cur_output_dir + converted_junc;
    
    if os.path.exists(output_dir_converted_junc) and \
	   rerun_all == 0:
	return output_dir_converted_junc
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "fusionsam2junc_filteranchor_newfmt", 
		        output_dir_converted_junc,
			str(read_width),
			str(min_anchor),
                        chromosomes_file_dir,
			output_dir_mps_unique_mapped]    
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % sam2junc_cmd
	
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert new sam file to junctions and filter by anchor failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fusionsam2junc_filteranchor_newfmt not found on this system"
        exit(1)
    return output_dir_converted_junc;

def extract_fusion_chr_seq(in_junc, converted_junc, chr_dir, extract_len, log_file):
    print >> sys.stderr, "[%s] extract fusion chromosome sequence" % (right_now())

    if os.path.exists(converted_junc) and \
	   rerun_all == 0:
	return converted_junc
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "extract_fusion_chr_seq", 
		        in_junc,
			converted_junc,
			chr_dir,
			str(extract_len)]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: extract fusion chromosome sequence"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: extract_fusion_chr_seq not found on this system"
        exit(1)
    return converted_junc;

def filter_fusion_by_repeat(fusion_junc, fusion_junc_out, fusion_junc_append_out, filtered_fusion, chrom_blat):
    
    print >> sys.stderr, "[%s] filter fusoin by repeat" % (right_now())

    if os.path.exists(filtered_fusion) and \
	   rerun_all == 0:
	return annot_gene
    
    blatJunct_log = open(logging_dir + "blatJunct.log", "w")
    
    filterBlat_log = open(logging_dir + "filterBlat.log", "w")
    
    filtered_fusion_log = open(filtered_fusion, "w")
    
    blatJunct_cmd = ["python", bin_dir + "blatJunct.py", fusion_junc, chrom_blat];
    
    filterBlat_cmd = ["python", bin_dir + "filterBlat.py", fusion_junc, fusion_junc_out];
    
    printJunct_cmd = ["python", bin_dir + "printJunct.py", fusion_junc_append_out];

    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % blatJunct_cmd
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % filterBlat_cmd
	
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % printJunct_cmd
	
    try:    
        retcode = subprocess.call(blatJunct_cmd, stdout=blatJunct_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter fusoin by repeat:blatJunct"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: blatJunct.py not found on this system"
        exit(1)
	
    try:    
        retcode = subprocess.call(filterBlat_cmd, stdout=filterBlat_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter fusoin by repeat:filterBlat"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterBlat.py not found on this system"
        exit(1)
	
    try:    
        retcode = subprocess.call(printJunct_cmd, stdout=filtered_fusion_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter fusoin by repeat:printJunct"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: printJunct.py not found on this system"
        exit(1)
	
    return filtered_fusion;



    
    
def annot_gene(fusion_junc, normal_junc, know_gene, annot_gene, log_file):
    print >> sys.stderr, "[%s] annotate gene from normal junction and fusion junction" % (right_now())

    if os.path.exists(annot_gene) and \
	   rerun_all == 0:
	return annot_gene
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "search_fusion_gene_new", 
		        '-g', know_gene,
			'-f', fusion_junc,
                        '-n', normal_junc,
                        '-o', annot_gene,
                        '-n_header']    
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % sam2junc_cmd
	
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: annotate gene from normal junction and fusion junction"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: search_fusion_gene_new not found on this system"
        exit(1)
    return annot_gene;

def fusionsam2maf(out_maf_file, chrom_size_file, fusion_sam_file, log_file):
    print >> sys.stderr, "[%s] convert fusion sam file to maf file" % (right_now())
   
    if os.path.exists(out_maf_file) and \
	   rerun_all == 0:
	return out_maf_file
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "fusionsam2maf", 
		        out_maf_file,
			chrom_size_file,
			fusion_sam_file]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert fusion sam file to maf file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fusionsam2maf not found on this system"
        exit(1)
    return out_maf_file;

def filterfusionjuncbystartend(normal_junc, fusion_junc, out_put_root, range, log_file):
    print >> sys.stderr, "[%s] filter fusion sam file to by start end" % (right_now())
    
    if os.path.exists(out_put_root + ".matched") and \
	   rerun_all == 0:
	return out_put_root + ".matched"
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log

    sam2junc_cmd = [bin_dir + "MatchStartEndFusion", 
		        normal_junc,
                      fusion_junc,
                      out_put_root,
			str(range)]    
    try:    
        retcode = subprocess.call(sam2junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter fusion sam file to by start end failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: MatchStartEndFusion not found on this system"
        exit(1)
    return out_put_root + ".matched";

def resortmappedreads(mps_unique_mapreads, mps_multiple_mapreads, bwtmapped_sam, bwt_25bp_mapped, resorted_offset, resorted_tagname, cur_output_dir):
    print >> sys.stderr, "[%s] resort mapped reads by offset and tagname" % (right_now())

    #output_dir_mps_unique_mapreads = cur_output_dir + mps_unique_mapreads;
    
    #output_dir_mps_multiple_mapreads = cur_output_dir + mps_multiple_mapreads;
    
    bwt_25bp_mapped_mishandt_matched = bwt_25bp_mapped + ".mistailandhead.matched";
    
    bwt_25bp_mapped_mishort_matched = bwt_25bp_mapped + ".mistailorhead.matched";
    bwt_25bp_mapped_allmapped = bwt_25bp_mapped + ".allmapped";
    
    output_dir_resorted_offset = cur_output_dir + resorted_offset;

    output_dir_resorted_tagname = cur_output_dir + resorted_tagname;
    
    if os.path.exists(output_dir_resorted_offset) and \
           os.path.exists(output_dir_resorted_tagname) and \
	   rerun_all == 0:
	return (output_dir_resorted_offset, output_dir_resorted_tagname)

    resortmapreads_cmd = [bin_dir + "resortmapreads", 
		        mps_unique_mapreads,
			mps_multiple_mapreads,
			bwt_25bp_mapped_mishandt_matched,
			bwt_25bp_mapped_mishort_matched,
			bwt_25bp_mapped_allmapped,
			bwtmapped_sam,
			output_dir_resorted_offset,
                        output_dir_resorted_tagname]    
    try:    
        retcode = subprocess.call(resortmapreads_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: resort mapped reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: resortmapreads not found on this system"
        exit(1)
    return (output_dir_resorted_offset, output_dir_resorted_tagname);



def call_rmap(syn_junctions, 
              seed_length,
              read_width,
              unmapped_reads,
              splice_mismatches,
              anchor_width,
              rank):
    start_time = datetime.now()
    print >> sys.stderr, "[%s] Aligning synthetic junction with rmap" % start_time.strftime("%c")
    
    rmap_cmd = ""
    
    juncfile_rmapped = output_dir + "junctions_rmapped.txt"

    rmap_cmd = [bin_dir + "rmap",
                "-L", str(read_width - anchor_width), # synthetic length
                "-m", str(splice_mismatches), # Mismatches allowed in extension
                "-h", str(seed_length), # Seed size for reads
                "-w", str(read_width), # read width for reads
                "-r", juncfile_rmapped, # output junction file
                "-u", unmapped_reads, # reads file or directory
                "-c", syn_junctions, # synthetic junctions
                "-R", str(rank), #rank for filter
                "whatever"
                ]
    try:    
        retcode = subprocess.call(rmap_cmd)
        
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Aligning synthetic junction with rmap failed"
            exit(1)
           
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: rmap not found on this system"
        exit(1)
       
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    
    return  juncfile_rmapped


def bsb_build_chromosome_index(user_splice_fasta, split_index_dir, splitsize, 
                               leftoverlapsize, rightoverlapsize, log_file):
    print >> sys.stderr, "[%s] Indexing bsb4 chromosome sequences" % (right_now())

    bowtie_build_log = open(log_file, "w")

    bowtie_build_path = bin_dir + "bowtie-build";
    
    bowtie_build_cmd = [bin_dir + "bsb4", 
                        bowtie_build_path,
                        str(splitsize),
                        str(leftoverlapsize),
                        str(rightoverlapsize),
                        split_index_dir,
                        user_splice_fasta]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] bsb4 commandline" % (bowtie_build_cmd)
    try:    
        retcode = subprocess.call(bowtie_build_cmd, 
                                  stdout=bowtie_build_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Splice bsb4 sequence indexing failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bsb4 not found on this system"
        exit(1)
    return user_splice_fasta

def build_chromosome_index(user_splice_fasta, user_splices_out_prefix):
    print >> sys.stderr, "[%s] Indexing chromosome sequences" % (right_now())
    bowtie_build_log = open(logging_dir + "bowtie_build.log", "w")
    
    bowtie_build_cmd = [bin_dir + "bowtie-build", 
                        user_splice_fasta,
                        user_splices_out_prefix]
    if DEBUG == 1:
	print >> sys.stderr, "[%s] bowtie-build commandline" % (bowtie_build_cmd)
    try:    
        retcode = subprocess.call(bowtie_build_cmd, 
                                  stdout=bowtie_build_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Splice sequence indexing failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie-build not found on this system"
        exit(1)
    return user_splices_out_prefix


def check_split_index(split_index_dir, chromo_dir, rerun, file_extension, 
                      splitsize, leftoverlapsize, rightoverlapsize):
    print >> sys.stderr, "[%s] Checking for Bowtie split index files" % right_now()
    
    all_chrom_paths = read_dir_by_suffix(chromo_dir, file_extension);
    
    for a_chrom_path in all_chrom_paths:

	chrom_file = os.path.basename(a_chrom_path);
	
	idx_prefix = split_index_dir + chrom_file + ".1";
	
	idx_fwd_1 = idx_prefix + ".1.ebwt"
	idx_fwd_2 = idx_prefix + ".2.ebwt"
	idx_fwd_3 = idx_prefix + ".3.ebwt"
	idx_fwd_4 = idx_prefix + ".4.ebwt"
	idx_rev_1 = idx_prefix + ".rev.1.ebwt"
	idx_rev_2 = idx_prefix + ".rev.2.ebwt"
	
	if os.path.exists(idx_fwd_3) and \
	   os.path.exists(idx_fwd_4) and \
	       rerun == 0:
	    continue 
	else:
	    
	    chromo_sequences = a_chrom_path;
	    
	    log_file = logging_dir + chrom_file + "_bowtie_build.log";
	    
	    bsb_build_chromosome_index(chromo_sequences, split_index_dir, 
	                                            splitsize, leftoverlapsize, rightoverlapsize, log_file);
	    
	    idx_fwd_1 = idx_prefix + ".1.ebwt"
	    idx_fwd_2 = idx_prefix + ".2.ebwt"
	    idx_fwd_3 = idx_prefix + ".3.ebwt"
	    idx_fwd_4 = idx_prefix + ".4.ebwt"
	    idx_rev_1 = idx_prefix + ".rev.1.ebwt"
	    idx_rev_2 = idx_prefix + ".rev.2.ebwt"
	    
	    #print >> sys.stderr, idx_fwd_1
	    #print >> sys.stderr, idx_fwd_2
	    #print >> sys.stderr, idx_fwd_3
	    #print >> sys.stderr, idx_fwd_4
	    #print >> sys.stderr, idx_rev_1
	    #print >> sys.stderr, idx_rev_2
	    
	    if os.path.exists(idx_fwd_3) and \
	       os.path.exists(idx_fwd_4):
		continue 
	    else:
		print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
		exit(1)
	    
def check_bowtie_index_combined(idx_prefix, combined_chromosome, rerun, file_extension):
    print >> sys.stderr, "[%s] Checking for Bowtie index files" % right_now()
    
    idx_fwd_1 = idx_prefix + ".1.ebwt"
    idx_fwd_2 = idx_prefix + ".2.ebwt"
    idx_fwd_3 = idx_prefix + ".3.ebwt"
    idx_fwd_4 = idx_prefix + ".4.ebwt"
    idx_rev_1 = idx_prefix + ".rev.1.ebwt"
    idx_rev_2 = idx_prefix + ".rev.2.ebwt"
    
    if os.path.exists(idx_fwd_1) and \
       os.path.exists(idx_fwd_2) and \
       os.path.exists(idx_fwd_3) and \
       os.path.exists(idx_fwd_4) and \
       os.path.exists(idx_rev_1) and \
       os.path.exists(idx_rev_2) and \
	   rerun == 0:
        return 
    else:
        idx_prefix = build_chromosome_index(combined_chromosome, idx_prefix);
        idx_fwd_1 = idx_prefix + ".1.ebwt"
        idx_fwd_2 = idx_prefix + ".2.ebwt"
        idx_fwd_3 = idx_prefix + ".3.ebwt"
        idx_fwd_4 = idx_prefix + ".4.ebwt"
        idx_rev_1 = idx_prefix + ".rev.1.ebwt"
        idx_rev_2 = idx_prefix + ".rev.2.ebwt"
        
        if os.path.exists(idx_fwd_1) and \
           os.path.exists(idx_fwd_2) and \
           os.path.exists(idx_fwd_3) and \
           os.path.exists(idx_fwd_4) and \
           os.path.exists(idx_rev_1) and \
           os.path.exists(idx_rev_2):
            return 
        else:
            print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
            exit(1)

def check_bowtie_index(idx_prefix, chromo_dir, rerun, file_extension):
    print >> sys.stderr, "[%s] Checking for Bowtie index files" % right_now()
    
    idx_fwd_1 = idx_prefix + ".1.ebwt"
    idx_fwd_2 = idx_prefix + ".2.ebwt"
    idx_fwd_3 = idx_prefix + ".3.ebwt"
    idx_fwd_4 = idx_prefix + ".4.ebwt"
    idx_rev_1 = idx_prefix + ".rev.1.ebwt"
    idx_rev_2 = idx_prefix + ".rev.2.ebwt"
    
    if os.path.exists(idx_fwd_1) and \
       os.path.exists(idx_fwd_2) and \
       os.path.exists(idx_fwd_3) and \
       os.path.exists(idx_fwd_4) and \
       os.path.exists(idx_rev_1) and \
       os.path.exists(idx_rev_2) and \
	   rerun == 0:
        return 
    else:
        chromo_sequences = read_sequence_by_suffix(chromo_dir, file_extension)
        idx_prefix = build_chromosome_index(chromo_sequences, idx_prefix);
        idx_fwd_1 = idx_prefix + ".1.ebwt"
        idx_fwd_2 = idx_prefix + ".2.ebwt"
        idx_fwd_3 = idx_prefix + ".3.ebwt"
        idx_fwd_4 = idx_prefix + ".4.ebwt"
        idx_rev_1 = idx_prefix + ".rev.1.ebwt"
        idx_rev_2 = idx_prefix + ".rev.2.ebwt"
        
        if os.path.exists(idx_fwd_1) and \
           os.path.exists(idx_fwd_2) and \
           os.path.exists(idx_fwd_3) and \
           os.path.exists(idx_fwd_4) and \
           os.path.exists(idx_rev_1) and \
           os.path.exists(idx_rev_2):
            return 
        else:
            print >> sys.stderr, "Error: Could not find Bowtie index files " + idx_prefix + ".*"
            exit(1)
	    
def bowtie(bwt_idx_prefix, 
           reads_list,
           reads_format, 
           mapped_reads,
           unmapped_reads, 
           bowtie_threads, 
           seed_length,
           max_hits,
	   unmapped_repeat_fasta_name,
	   mismatch,
	   log_file):
    start_time = datetime.now()
    bwt_idx_name = bwt_idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Mapping reads against %s with Bowtie" % (start_time.strftime("%c"), bwt_idx_name)
    
    bwt_map = temp_dir + mapped_reads

    unmapped_reads_fasta_name = temp_dir + unmapped_reads
    
    if os.path.exists(unmapped_reads_fasta_name) and \
	   os.path.exists(bwt_map) and \
	   rerun_all == 0:
	return (bwt_map, unmapped_reads_fasta_name)
    #os.path.exists(unmapped_repeat_fasta_name) and \
    bwt_log = open(log_file, "w")

    # Launch Bowtie
    try:    
        bowtie_cmd = [bin_dir + "bowtie"]
	
	unampped_format = "--un"
	
	repeat_format = "--max"
	
	if reads_format == "-q":
	    unampped_format = "--un"
	    repeat_format = "--max"
	    
	if mismatch > 3:
	    mismatch = 3
     
        bowtie_cmd += [reads_format,
                       "--threads", str(bowtie_threads),
                       unampped_format, unmapped_reads_fasta_name,
                       "-k", str(max_hits),
                       #"-m", str(max_hits + 1),
		       "-v", str(mismatch),
	               "--best",
                       repeat_format, unmapped_repeat_fasta_name,
                       bwt_idx_prefix, 
                       reads_list, 
                       bwt_map]   

	if DEBUG == 1:
	    print >> sys.stderr, "[%s] bowtie commandline" % (bowtie_cmd)
        ret = subprocess.call(bowtie_cmd, 
                              stderr=bwt_log)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute Bowtie"
            exit(1)
            
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to put it in the directory?"
	    exit(1)

    # Success    
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return (bwt_map, unmapped_reads_fasta_name)

def bowtie2(bwt_idx_prefix, 
           reads_list,
           reads_format, 
           mapped_reads,
           unmapped_reads, 
           bowtie_threads, 
           seed_length,
           max_hits,
	   unmapped_repeat_fasta_name,
	   mismatch,
	   log_file):
    start_time = datetime.now()
    bwt_idx_name = bwt_idx_prefix.split('/')[-1]
    print >> sys.stderr, "[%s] Mapping reads against %s with Bowtie" % (start_time.strftime("%c"), bwt_idx_name)
    
    bwt_map = temp_dir + mapped_reads

    unmapped_reads_fasta_name = temp_dir + unmapped_reads
   
    
    if os.path.exists(unmapped_reads_fasta_name) and \
	   os.path.exists(bwt_map) and \
	   rerun_all == 0:
	return (bwt_map, unmapped_reads_fasta_name)
    #os.path.exists(unmapped_repeat_fasta_name) and \
    bwt_log = open(log_file, "w")

    # Launch Bowtie
    try:    
        bowtie_cmd = [bin_dir + "bowtie"]
	
	unampped_format = "--unfa"
	
	repeat_format = "--maxfa"
	
	if reads_format == "-q":
	    unampped_format = "--unfq"
	    repeat_format = "--maxfq"
	    
	if mismatch > 3:
	    mismatch = 3
     
        bowtie_cmd += [reads_format,
                       "--threads", str(bowtie_threads),
                       unampped_format, unmapped_reads_fasta_name,
                       "-k", str(0),
                       "-m", str(max_hits + 1),
		       "-v", str(mismatch),
                       repeat_format, unmapped_repeat_fasta_name,
                       bwt_idx_prefix, 
                       reads_list, 
                       bwt_map]   

	#print >> sys.stderr, "[%s] bowtie commandline" % (bowtie_cmd)
        ret = subprocess.call(bowtie_cmd, 
                              stderr=bwt_log)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute Bowtie"
            exit(1)
            
    # Bowtie not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: Bowtie not found on this system.  Did you forget to put it in the directory?"
	    exit(1)

    # Success    
    finish_time = datetime.now()
    duration = finish_time - start_time
    #print >> sys.stderr, "\t\t\t[%s elapsed]" %  formatTD(duration)
    return (bwt_map, unmapped_reads_fasta_name)

def bowtie2sam(bowtie_mapped, misinfo, sam_formatted, log_file):
    print >> sys.stderr, "[%s] Converting bowtie mapped to SAM format" % (right_now())

    if os.path.exists(sam_formatted) and \
	   rerun_all == 0:
	return sam_formatted
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    #bowtie2sam_cmd = ["bowtie2sam", 
                        #bowtie_mapped,
                        #sam_formatted]  
    bowtie2sam_cmd = [bin_dir + "bowtie2sam", 
                        bowtie_mapped,
                        sam_formatted,
			str(misinfo)]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert bowtie mapped to SAM failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: bowtie2sam not found on this system"
        exit(1)
    return sam_formatted

def FilterBWTByRegions(bowtie_mapped, regions, bowtie_mapped_in_regions, bowtie_mapped_notin_regions, log_file):
    print >> sys.stderr, "[%s] Remove bwt mapped by regions" % (right_now())

    if os.path.exists(bowtie_mapped_in_regions) and \
       os.path.exists(bowtie_mapped_notin_regions) and \
       rerun_all == 0:
	return (bowtie_mapped_in_regions, bowtie_mapped_notin_regions)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log    
    
    bowtie2sam_cmd = [bin_dir + "FilterBWTByRegions", 
                        regions,
                        bowtie_mapped,
			bowtie_mapped_notin_regions,
			bowtie_mapped_in_regions]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % bowtie2sam_cmd
	
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Remove bwt mapped by regions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterBWTByRegions not found on this system"
        exit(1)
	
    return (bowtie_mapped_in_regions, bowtie_mapped_notin_regions)

def fa_fq_oneline(divided_reads, sorted_divide_reads, format_flag):
    print >> sys.stderr, "[%s] fa_fq_oneline divide reads" % (right_now())
 
    if os.path.exists(sorted_divide_reads) and \
	   rerun_all == 0:
	return sorted_divide_reads
    
    if format_flag == "-f":
	line_per_reads = 2
    elif format_flag == "-q":
	line_per_reads = 4
     
    bowtie2sam_cmd = [bin_dir + "fa_fq_oneline", 
                        divided_reads,
                        sorted_divide_reads,
			str(line_per_reads)]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: fa_fq_oneline divide reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: fa_fq_oneline not found on this system"
        exit(1)
    return sorted_divide_reads

def cat_files(files_tobe_cat, cated_file):
    print >> sys.stderr, "[%s] combine to file %s" % (right_now(), cated_file)
   
    if os.path.exists(cated_file) and \
	   rerun_all == 0:
	return cated_file
	
    bowtie2sam_cmd = ["cat"];

    for sam_file in files_tobe_cat:
	bowtie2sam_cmd = bowtie2sam_cmd + [sam_file]
    
    cated_file_fs = open(cated_file, "w")
    
    #bowtie2sam_cmd = bowtie2sam_cmd + [">"] + [cated_file]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] combine to file" % bowtie2sam_cmd
	
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=cated_file_fs)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: combine to final sam"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cat not found on this system"
        exit(1)
	
    return cated_file

def remove_files(files_tobe_remove):
    print >> sys.stderr, "[%s] remove files %s"  % (right_now(), files_tobe_remove)
   
    if DEBUG > 0:
	return;
    
    bowtie2sam_cmd = ["rm"];

    for sam_file in files_tobe_remove:
	bowtie2sam_cmd = bowtie2sam_cmd + [sam_file]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] remove files" % bowtie2sam_cmd
	
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: remove files failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: rm not found on this system"
        exit(1)

def copy_file(files_tobe_copy, copied_file):
    print >> sys.stderr, "[%s] copy file %s"  % (right_now(), copied_file)
   
    if os.path.exists(copied_file) and \
	   rerun_all == 0:
	return copied_file
	
    bowtie2sam_cmd = ["cp",
                      files_tobe_copy,
                      copied_file];

    if DEBUG == 1:
	print >> sys.stderr, "[%s] copy file" % bowtie2sam_cmd
	
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: copy file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cp not found on this system"
        exit(1)
	
    return copied_file


def merge_sort_segment_bwt1(tobe_sorted_1, tobe_sorted_2, merged_sorted):
    print >> sys.stderr, "[%s] merge sort files" % (right_now())
   
    if os.path.exists(merged_sorted) and \
	   rerun_all == 0:
	#temp_fs = open(segment_bwt, "w")
	#temp_fs.close()
	return merged_sorted
    
    bowtie2sam_cmd = ["sort",
		      "-t~",
		      "-k1,1n",
		      "-S", "3500000",
		      "-o", merged_sorted, 
		      "-T", temp_dir,
		      tobe_sorted_1,
		      tobe_sorted_2]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort files failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    #temp_fs = open(segment_bwt, "w")
    #temp_fs.close()
    return merged_sorted

def merge_sort_segment_bwt2(tobe_sorted_1, merged_sorted):
    print >> sys.stderr, "[%s] merge sort files by chromosome" % (right_now())
   
    if os.path.exists(merged_sorted) and \
	   rerun_all == 0:
	#temp_fs = open(segment_bwt, "w")
	#temp_fs.close()
	return merged_sorted

#sort -k2,2 -k4,4n -o sorted_combine_head_tail.txt.remained.formated.sorted 
#sorted_combine_head_tail.txt.remained.formated 
    
    bowtie2sam_cmd = ["sort",
		      "-k2,2",
		      "-k4,4n",
		      "-S", "3500000",
		      "-o", merged_sorted, 
		      "-T", temp_dir,
		      tobe_sorted_1]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge sort files by chromosome failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    #temp_fs = open(segment_bwt, "w")
    #temp_fs.close()
    return merged_sorted

def sort_segment_bwt1(segment_bwt, segment_bwt_sorted):
    print >> sys.stderr, "[%s] sort segmentbwt" % (right_now())
   
    if os.path.exists(segment_bwt_sorted) and \
	   rerun_all == 0:
	temp_fs = open(segment_bwt, "w")
	temp_fs.close()
	return segment_bwt_sorted
    
    bowtie2sam_cmd = ["sort",
		      "-t~",
		      "-k1,1n",
		      "-S", "3500000",
		      "-o", segment_bwt_sorted, 
		      "-T", temp_dir,
		      segment_bwt]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort segmentbwt failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    temp_fs = open(segment_bwt, "w")
    temp_fs.close()
    return segment_bwt_sorted

def sort_segment_bwt(segment_bwt, segment_bwt_sorted):
    print >> sys.stderr, "[%s] sort segmentbwt" % (right_now())
   
    if os.path.exists(segment_bwt_sorted) and \
	   rerun_all == 0:
	temp_fs = open(segment_bwt, "w")
	temp_fs.close()
	return segment_bwt_sorted
    
    bowtie2sam_cmd = ["sort",
		      "-t~",
		      "-k2,2n",
		      "-S", "3500000",
		      "-o", segment_bwt_sorted, 
		      "-T", temp_dir,
		      segment_bwt]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort segmentbwt failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    temp_fs = open(segment_bwt, "w")
    temp_fs.close()
    return segment_bwt_sorted

def sort_by_name(segment_bwt, segment_bwt_sorted):
    print >> sys.stderr, "[%s] sort by name" % (right_now())
   
    if os.path.exists(segment_bwt_sorted) and \
	   rerun_all == 0:
	temp_fs = open(segment_bwt, "w")
	temp_fs.close()
	return segment_bwt_sorted
    
    bowtie2sam_cmd = ["sort",
		      "-t~",
		      "-k2,2",
		      "-S", "3500000",
		      "-o", segment_bwt_sorted, 
		      "-T", temp_dir,
		      segment_bwt]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort by name failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    temp_fs = open(segment_bwt, "w")
    temp_fs.close()
    return segment_bwt_sorted

def sort_by_name_c(files_tobe_sort, sorted_file):
    print >> sys.stderr, "[%s] sort by name" % (right_now())
   
    if os.path.exists(sorted_file) and \
	   rerun_all == 0:
	return sorted_file
    
    env = os.environ.copy()
    
    env.update({"LC_ALL": "C"})
    
    bowtie2sam_cmd = ["sort",
		      "-k1,1",
		      "-S", "3500000",
		      "-o", sorted_file, 
		      "-T", temp_dir]
    
    for sam_file in files_tobe_sort:
	bowtie2sam_cmd = bowtie2sam_cmd + [sam_file]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % bowtie2sam_cmd
	
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, env = env)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort by name failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    return sorted_file

def sort_by_name1(segment_bwt, segment_bwt_sorted):
    print >> sys.stderr, "[%s] sort by name" % (right_now())
   
    if os.path.exists(segment_bwt_sorted) and \
	   rerun_all == 0:
	temp_fs = open(segment_bwt, "w")
	temp_fs.close()
	return segment_bwt_sorted
    
    bowtie2sam_cmd = ["sort",
		      "-k1,1",
		      "-S", "3500000",
		      "-o", segment_bwt_sorted, 
		      "-T", temp_dir,
		      segment_bwt]    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: sort by name failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sort not found on this system"
        exit(1)
	
    temp_fs = open(segment_bwt, "w")
    temp_fs.close()
    return segment_bwt_sorted

def mapsplice_search(min_intron, max_intron, seg_len, ext_len, sorted_read_file, format_flag, chrom_size_file, 
		     whole_reads_bwt, island_ext, island_gap, island_output, segment_bwt, pairend, max_insert):
    print >> sys.stderr, "[%s] mapsplice_search" % (right_now())
 
    if os.path.exists(island_output) and \
	   rerun_all == 0:
	return island_output
    
    mapsplice_log = open(logging_dir + "mapsplice_search.log", "w")
    #stdout=mapsplice_log
    #./mapsplice_search -i 0 -I 50000 -l 25 -n 4 -e 3 -ie 45 -g 50 -r divided_reads_oneline.fa.sorted -m unspliced_map.bwtout -c chrom.txt 
    #-o1 ./hole/ -o2 ./hmer/ -o3 ./head/ -o4 ./tail/ -o5 ./bwtout/ -o6 island.txt -fa unspliced_map_25.bwtout
    read_format = "";
    
    if format_flag == "-f":
	read_format = "-fa"
    elif format_flag == "-q":
	read_format = "-fq"
    
    bowtie2sam_cmd = [bin_dir + "mapsplice_search", 
		        "--max_insertion", str(max_insert),
                        "-i", str(min_intron),
                        "-I", str(max_intron),
			"-l", str(seg_len),
			"-e", str(ext_len),
			"-ie", str(island_ext),
			"-g", str(island_gap),
			"-r", str(sorted_read_file),
			"-m", str(whole_reads_bwt),
			"-c", str(chrom_size_file),
			"-o1", hole_dir,
			"-o2", hmer_dir,
			"-o3", head_dir,
			"-o4", tail_dir,
			"-o5", island_output,
			read_format];
    
    if pairend != "":
	bowtie2sam_cmd = bowtie2sam_cmd + ["-pe"];
	
    bowtie2sam_cmd = bowtie2sam_cmd + [segment_bwt]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] mapsplice_search" % bowtie2sam_cmd
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: fmapsplice_search failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_search not found on this system"
        exit(1)
    return island_output

def mapsplice_report(seg_len, chrom_size_file, ext_bit, sorted_read, format_flag, unspliced_mapped_segment, 
		     max_intron, fusion_file, single_sam_file, paired_sam_file, min_seg, min_fusion_seg, do_fusion, paired):
    print >> sys.stderr, "[%s] mapsplice_report" % (right_now())
 
    #./mapsplice_report -l 25 -n 4 -c chrom.txt -i1 bwtout/ -i2 hole/ -i3 hmer_sort/ -i4 head_sort/ -i5 tail_sort/ -o sam/

    if os.path.exists(single_sam_file) and \
	   rerun_all == 0:
	return single_sam_file
    
    mapsplice_log = open(logging_dir + "mapsplice_report.log", "w")
    #stdout=mapsplice_log
    
    formatfq = "-fq";
    if format_flag == "-f":
	formatfq = "-fa";	
    
    bowtie2sam_cmd = [bin_dir + "mapsplice_report", 
			"-l", str(seg_len),
			"-c", chrom_size_file,
			"-e", str(ext_bit),
			"-r", sorted_read,
			"-s", str(min_seg),
			"-I", str(max_intron),
			formatfq,
			"-i1", unspliced_mapped_segment,
			"-i2", hole_dir,
			"-i3", hmer_dir,
			"-i4", head_dir,
			"-i5", tail_dir,
			"-o", single_sam_file]
    
    if paired_sam_file != "":
	bowtie2sam_cmd = bowtie2sam_cmd + ["-o2", paired_sam_file] + ["-o3", paired_sam_file + ".filtered"];
    
    if paired != "":
	bowtie2sam_cmd = bowtie2sam_cmd + ["-pe"];
	
    if do_fusion == 1: 
	bowtie2sam_cmd = bowtie2sam_cmd + ["-f", fusion_dir]# + ["-v", str(min_fusion_seg)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] mapsplice_report" % bowtie2sam_cmd
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: mapsplice_report failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: mapsplice_report not found on this system"
        exit(1)
	
    return single_sam_file

def RemDup(syn_mapped, syn_mapped_remdup_unique_spliced, syn_mapped_remdup_multiple_spliced, stat_file, log_file): #syn_mapped_remdup_unspliced, 
    print >> sys.stderr, "[%s] Remove duplication of sam" % (right_now())

    if os.path.exists(syn_mapped_remdup_unique_spliced) and \
       os.path.exists(syn_mapped_remdup_multiple_spliced) and \
 	   rerun_all == 0:
	#temp_fs = open(syn_mapped, "w")
	#temp_fs.close()
	return (syn_mapped_remdup_unique_spliced, syn_mapped_remdup_multiple_spliced);

    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    unionjunc_cmd = [bin_dir + "RemDup"] + [syn_mapped] + [syn_mapped_remdup_unique_spliced] #[syn_mapped_remdup_unspliced] + 
    
    unionjunc_cmd = unionjunc_cmd + [syn_mapped_remdup_multiple_spliced] + [stat_file];
    
    try:    
        retcode = subprocess.call(unionjunc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Remove duplication failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: RemDup not found on this system"
        exit(1)
	
    #temp_fs = open(syn_mapped, "w")
    #temp_fs.close()

    return (syn_mapped_remdup_unique_spliced, syn_mapped_remdup_multiple_spliced);

def SepSplicedUnSpliced(combined_sam, spliced_sam, unspliced_sam, stat_file, log_file): #syn_mapped_remdup_unspliced, 
    print >> sys.stderr, "[%s] Separate Spliced UnSpliced sam" % (right_now())

    if os.path.exists(spliced_sam) and \
       os.path.exists(unspliced_sam) and \
 	   rerun_all == 0:
	return (spliced_sam, unspliced_sam);

    mapsplice_log = open(log_file, "w")
    
    unionjunc_cmd = [bin_dir + "SepSplicedUnspliced"] + [combined_sam] + [spliced_sam] + [unspliced_sam] + [stat_file];#[syn_mapped_remdup_unspliced] + 

    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % unionjunc_cmd
    try:    
        retcode = subprocess.call(unionjunc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Separate Spliced UnSpliced sam failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: SepSplicedUnspliced not found on this system"
        exit(1)

    return (spliced_sam, unspliced_sam);

def FilterByParing(tobepairedsam, pairedsam, fusionpairedsam, singlesam, filteredsam, max_dist, stat_file, log_file):
    print >> sys.stderr, "[%s] filter by paring" % (right_now())

    if os.path.exists(pairedsam) and \
       os.path.exists(filteredsam) and \
       rerun_all == 0:
	#temp_fs = open(syn_mapped, "w")
	#temp_fs.close()
	return (pairedsam, filteredsam);

    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    unionjunc_cmd = [bin_dir + "pairing"] + [tobepairedsam] + [pairedsam] + [fusionpairedsam] + [singlesam] + [filteredsam];
    
    unionjunc_cmd = unionjunc_cmd + [stat_file] + [str(max_dist)];
    
    try:    
        retcode = subprocess.call(unionjunc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by paring failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pairing not found on this system"
        exit(1)
	
    #temp_fs = open(syn_mapped, "w")
    #temp_fs.close()

    return (pairedsam, filteredsam);

def FilterReadsByCanonNoncanon(syn_mapped_remdup, junc_file, filtered_canon, filtered_noncanon, filtered_noncanon_canon, filtered_ins, stat_file, log_file):
    print >> sys.stderr, "[%s] Filter Reads By Canon Noncanon" % (right_now())

    if os.path.exists(filtered_canon) and \
           os.path.exists(filtered_noncanon) and \
           os.path.exists(filtered_noncanon_canon) and \
	   rerun_all == 0:
	
	#temp_fs = open(syn_mapped_remdup, "w")
	#temp_fs.close()
	
	return (filtered_canon, filtered_noncanon, filtered_noncanon_canon)
    
    unionjunc_cmd = [bin_dir + "FilterReadsByCanonNoncanonByReads",
		     syn_mapped_remdup, 
		     junc_file,
		     filtered_canon,
		     filtered_noncanon,
		     filtered_noncanon_canon,
		     filtered_ins,
		     stat_file]
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(unionjunc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filter Reads By Canon Noncanon failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterReadsByCanonNoncanon not found on this system"
        exit(1)
	
    #temp_fs = open(syn_mapped_remdup, "w")
    #temp_fs.close()

    return (filtered_canon, filtered_noncanon, filtered_noncanon_canon);

def juncdb_bwt2sam(bowtie_mapped, misinfo, sam_formatted, log_file):
    print >> sys.stderr, "[%s] Converting juncdb bowtie mapped to SAM format" % (right_now())
    
    if os.path.exists(sam_formatted) and \
	   rerun_all == 0:
	return sam_formatted
    #bowtie2sam_cmd = ["bowtie2sam", 
                        #bowtie_mapped,
                        #sam_formatted]  
    bowtie2sam_cmd = [bin_dir + "juncdb_bwt2sam", 
                        bowtie_mapped,
                        sam_formatted,
			str(misinfo)] 
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert juncdb bowtie mapped to SAM failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: juncdb_bwt2sam not found on this system"
        exit(1)
    
    return sam_formatted

def juncdb_fusion_bwt2sam(bowtie_mapped, misinfo, sam_formatted, syn_len, log_file):
    print >> sys.stderr, "[%s] Converting fusion juncdb bowtie mapped to SAM format" % (right_now())
    
    if os.path.exists(sam_formatted) and \
	   rerun_all == 0:
	return sam_formatted
    #bowtie2sam_cmd = ["bowtie2sam", 
                        #bowtie_mapped,
                        #sam_formatted]  
    bowtie2sam_cmd = [bin_dir + "juncdb_fusion_bwt2sam", 
                        bowtie_mapped,
                        sam_formatted,
                        str(syn_len),
			str(misinfo)] 
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] Converting fusion juncdb bowtie mapped to SAM format" % bowtie2sam_cmd
    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert fusion juncdb bowtie mapped to SAM failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: juncdb_fusion_bwt2sam not found on this system"
        exit(1)
	
    return sam_formatted

def replace_fusion_seq(bowtie_mapped, bowtie_mapped_replaced, synthe_fusion_file, log_file):
    print >> sys.stderr, "[%s] Replace fusion read sequence with syn chromosome sequence" % (right_now())
    
    if os.path.exists(bowtie_mapped_replaced) and \
	   rerun_all == 0:
	return bowtie_mapped_replaced

    bowtie2sam_cmd = [bin_dir + "replace_fusion_seq", 
                        bowtie_mapped,
                        bowtie_mapped_replaced,
                        synthe_fusion_file] 
    
    mapsplice_log = open(log_file, "w")
    
    try:    
        retcode = subprocess.call(bowtie2sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Replace fusion read sequence with syn chromosome sequence failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: replace_fusion_seq not found on this system"
        exit(1)
	
    return bowtie_mapped_replaced

def dividereads(reads_files, seg_len, format_flag):
    print >> sys.stderr, "[%s] divide reads" % (right_now())

    mapsplice_log = open(logging_dir + "dividereads.log", "w")
    #stdout=mapsplice_log
    
    divided_reads = temp_dir + "divided_reads.fa"
    
    divided_reads_trunc = temp_dir + "divided_reads_trunc.fa"
    
    if os.path.exists(divided_reads) and \
	   rerun_all == 0:
	
	#temp_fs = open(reads_file, "w")
	#temp_fs.close()    
	return divided_reads
    
    read_files_array = reads_files.split(',')
    
    dividereads_cmd = [bin_dir + "dividereads"]
    
    for read_file in read_files_array:
	dividereads_cmd = dividereads_cmd + [read_file];
    
    format_int = 1;
    if format_flag == "-q":
	format_int = 0;

    dividereads_cmd = dividereads_cmd + [divided_reads_trunc] + [divided_reads] + [str(format_int)] + [str(seg_len)]   
    if DEBUG == 1:
	print >> sys.stderr, "[%s] divide reads" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: divide reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: dividereads not found on this system"
        exit(1)
	
    #temp_fs = open(reads_file, "w")
    #temp_fs.close()
    return divided_reads

def dividereadsindexed(reads_files, seg_len, format_flag):
    print >> sys.stderr, "[%s] divide reads" % (right_now())

    mapsplice_log = open(logging_dir + "dividereads.log", "w")
    #stdout=mapsplice_log
    
    divided_reads = temp_dir + "divided_reads.fa"
    
    divided_reads_trunc = temp_dir + "divided_reads_trunc.fa"
    
    if os.path.exists(divided_reads) and \
	   rerun_all == 0:
	
	#temp_fs = open(reads_file, "w")
	#temp_fs.close()    
	return divided_reads
    
    read_files_array = reads_files.split(',')
    
    dividereads_cmd = [bin_dir + "dividereadsindexed"]
    
    for read_file in read_files_array:
	dividereads_cmd = dividereads_cmd + [read_file];
    
    format_int = 1;
    if format_flag == "-q":
	format_int = 0;

    dividereads_cmd = dividereads_cmd + [divided_reads_trunc] + [divided_reads] + [str(format_int)] + [str(seg_len)]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] divide reads" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: divide reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: dividereads not found on this system"
        exit(1)
	
    #temp_fs = open(reads_file, "w")
    #temp_fs.close()
    return divided_reads

def remove_index_sam(indexed_sam_file, remove_indexed_sam_file):
    print >> sys.stderr, "[%s] remove index of sam" % (right_now())

    mapsplice_log = open(logging_dir + "remove_index_sam.log", "w")
    #stdout=mapsplice_log
    
    if os.path.exists(remove_indexed_sam_file) and \
	   rerun_all == 0:
	
	#temp_fs = open(reads_file, "w")
	#temp_fs.close()    
	return remove_indexed_sam_file
    
    dividereads_cmd = [bin_dir + "remove_index_sam"]
    
    dividereads_cmd = dividereads_cmd + [indexed_sam_file] + [remove_indexed_sam_file];
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] remove index of sam" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: remove index of sam"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: remove_index_sam not found on this system"
        exit(1)
	
    #temp_fs = open(reads_file, "w")
    #temp_fs.close()
    return remove_indexed_sam_file

def merge_paired_end_reads(reads_files, format_flag, paired, merged_reads):
    print >> sys.stderr, "[%s] merge paired end reads" % (right_now())

    mapsplice_log = open(logging_dir + "merge_paired_end_reads.log", "w")
    #stdout=mapsplice_log
     
    if os.path.exists(merged_reads) and \
       rerun_all == 0:
	
	#temp_fs = open(reads_file, "w")
	#temp_fs.close()    
	return merged_reads
    
    read_files_array = reads_files.split(',')
  
    dividereads_cmd = [bin_dir + "merge_paired_end_reads"]
    
    for read_file in read_files_array:
	dividereads_cmd = dividereads_cmd + [read_file];
    
    format_int = 1;
    if format_flag == "-q":
	format_int = 0;
	
    ispaired = 0;
    
    if paired != "":
	ispaired = 1

    dividereads_cmd = dividereads_cmd + [str(ispaired)] + [str(format_int)] + [merged_reads] 
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] merge paired end reads" % dividereads_cmd
	
    #print >> sys.stderr, "[%s] divide reads" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge paired end reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_paired_end_reads not found on this system"
        exit(1)
	
    #temp_fs = open(reads_file, "w")
    #temp_fs.close()
    return merged_reads

def merge_paired_end_reads_rmshort(reads_files, format_flag, paired, seg_len, merged_reads):
    print >> sys.stderr, "[%s] merge paired end reads remove short" % (right_now())

    mapsplice_log = open(logging_dir + "merge_paired_end_reads_rmshort.log", "a")
    #stdout=mapsplice_log
     
    if os.path.exists(merged_reads) and \
       rerun_all == 0:
	
	#temp_fs = open(reads_file, "w")
	#temp_fs.close()    
	return merged_reads
    
    read_files_array = reads_files.split(',')
  
    dividereads_cmd = [bin_dir + "merge_paired_end_reads_rmshort"]
    
    for read_file in read_files_array:
	dividereads_cmd = dividereads_cmd + [read_file];
    
    format_int = 1;
    if format_flag == "-q":
	format_int = 0;
	
    ispaired = 0;
    
    if paired != "":
	ispaired = 1

    dividereads_cmd = dividereads_cmd + [str(seg_len)] + [str(ispaired)] + [str(format_int)] + [merged_reads] 

    if DEBUG == 1:
	print >> sys.stderr, "[%s] remove index of sam" % dividereads_cmd
    #print >> sys.stderr, "[%s] divide reads" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge paired end reads remove short failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_paired_end_reads_rmshort not found on this system"
        exit(1)
	
    #temp_fs = open(reads_file, "w")
    #temp_fs.close()
    return merged_reads

def check_reads_format(reads_files, min_read_len, paired):
    print >> sys.stderr, "[%s] check reads format" % (right_now())

    mapsplice_log = open(logging_dir + "check_reads_format.log", "a")
    #stdout=mapsplice_log
         
    read_files_array = reads_files.split(',')
  
    dividereads_cmd = [bin_dir + "check_reads_format"]
    
    for read_file in read_files_array:
	dividereads_cmd = dividereads_cmd + [read_file];
    
    #format_int = 1;
    #if format_flag == "-q":
	#format_int = 0;
	
    ispaired = 0;
    
    if paired != "":
	ispaired = 1

    dividereads_cmd = dividereads_cmd + [str(0)] + [str(ispaired)] + [str(min_read_len)] 

    if DEBUG == 1:
	print >> sys.stderr, "[%s] check reads format" % dividereads_cmd
    #print >> sys.stderr, "[%s] divide reads" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: check reads format failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: check_reads_format not found on this system"
        exit(1)
	

def index_reads(reads_files, format_flag, merged_reads):
    print >> sys.stderr, "[%s] index reads" % (right_now())

    mapsplice_log = open(logging_dir + "index_reads.log", "w")
    #stdout=mapsplice_log
     
    if os.path.exists(merged_reads) and \
       rerun_all == 0:
	
	#temp_fs = open(reads_file, "w")
	#temp_fs.close()    
	return merged_reads
    
    read_files_array = reads_files.split(',')
    
    #if len(read_files_array) == 1:
	#return reads_files
    
    dividereads_cmd = [bin_dir + "index_reads"]
    
    for read_file in read_files_array:
	dividereads_cmd = dividereads_cmd + [read_file];
    
    format_int = 1;
    if format_flag == "-q":
	format_int = 0;

    dividereads_cmd = dividereads_cmd + [str(format_int)] + [merged_reads] 
    #print >> sys.stderr, "[%s] divide reads" % dividereads_cmd
    try:    
        retcode = subprocess.call(dividereads_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: index failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: index_reads not found on this system"
        exit(1)
	
    #temp_fs = open(reads_file, "w")
    #temp_fs.close()
    return merged_reads

def parsepileup(pileup_file, threshold, boundary):
    print >> sys.stderr, "[%s] Parsing pileup file" % (right_now())

    islands_file = temp_dir + "islands.gff"
    
    if os.path.exists(islands_file) and \
	   rerun_all == 0:
	temp_fs = open(pileup_file, "w")
	temp_fs.close()
	return islands_file
    
    parsepileup_cmd = [bin_dir + "parsepileup", pileup_file, islands_file, 
                      str(threshold), str(boundary), temp_dir]            
    try:    
        retcode = subprocess.call(parsepileup_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Parsing pileup file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: parsepileup not found on this system"
        exit(1)
	
    temp_fs = open(pileup_file, "w")
    temp_fs.close()
    
    return islands_file

#sort_bowtie(char* bowtiefile, char* sortedbowtie, int maximal_intron, const char* readsfile, 
#            string chromdir, int num_seg, int seg_len, char* merged_bowtie)
def merge_comb_bowtie(bowtiefile, sortedbowtie, maximal_intron, readsfile,
		      chromdir, num_seg, seg_len, merged_bowtie):
    print >> sys.stderr, "[%s] Merging combine bowtie file" % (right_now())
    
    mergecombbowtie_cmd = [bin_dir + "merge_comb_bowtie", bowtiefile, sortedbowtie, str(maximal_intron), 
			   readsfile, chromdir, str(num_seg), str(seg_len), merged_bowtie]
    
    if os.path.exists(merged_bowtie) and \
           os.path.exists(sortedbowtie) and \
	   rerun_all == 0:
	
	temp_fs = open(bowtiefile + ".notcombined", "w")
	temp_fs.close()
	
	return bowtiefile
    try:    
        retcode = subprocess.call(mergecombbowtie_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Merging combine bowtie file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_comb_bowtie not found on this system"
        exit(1)
	
    temp_fs = open(bowtiefile + ".notcombined", "w")
    temp_fs.close()
    
    return bowtiefile

def merge_sam(bowtie_mapped_sam, bowtie_mapped_sam_25):
    print >> sys.stderr, "[%s] Merging bowtie sam file" % (right_now())

    #pileup_file = bowtie_mapped_sam + ".pileup"
    #pileup_fs = open(pileup_file, "w")
   
    merged_bam = temp_dir + "merged.sam"
    
    if os.path.exists(merged_bam) and \
	   rerun_all == 0:
	temp_fs = open(bowtie_mapped_sam, "w")
	temp_fs.close()
	temp_fs = open(bowtie_mapped_sam_25, "w")
	temp_fs.close()
	return merged_bam
    
    merge_sam_cmd = [bin_dir + "merge_sam", bowtie_mapped_sam, bowtie_mapped_sam_25, merged_bam]  
    try:    
        retcode = subprocess.call(merge_sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge bowtie failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sam not found on this system"
        exit(1)
	
    temp_fs = open(bowtie_mapped_sam, "w")
    temp_fs.close()
    temp_fs = open(bowtie_mapped_sam_25, "w")
    temp_fs.close()
    return merged_bam

def merge_chromo_sams(all_sams_path, merged_sam, log_file):
    print >> sys.stderr, "[%s] Merging all sam files" % (right_now())
    
    if os.path.exists(merged_sam) and \
	   rerun_all == 0:
	return merged_sam
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    merge_sam_cmd = [bin_dir + "merge_sam"]  
    
    for sam_file in all_sams_path:
	merge_sam_cmd = merge_sam_cmd + [sam_file]
    
    merge_sam_cmd = merge_sam_cmd + [merged_sam]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] Merging all sam files" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Merging all sam files failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sam not found on this system"
        exit(1)
	
    return merged_sam

def read_chromo_sizes(all_sams_path, chromo_size_file, chrom_names, chrom_head, chromo_fai_file):
    print >> sys.stderr, "[%s] reads all chromo sizes" % (right_now())
    
    if os.path.exists(chromo_size_file) and \
       os.path.exists(chrom_names) and \
	   rerun_all == 0:
	return chromo_size_file
    
    mapsplice_log = open(logging_dir + "read_chromo_sizes.log", "w")
    #stdout=mapsplice_log
    
    merge_sam_cmd = [bin_dir + "read_chromo_size"]  
    
    for sam_file in all_sams_path:
	merge_sam_cmd = merge_sam_cmd + [sam_file]
    
    merge_sam_cmd = merge_sam_cmd + [chromo_fai_file] + [chrom_head] + [chromo_size_file] + [chrom_names]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] reads all chromo sizes" % merge_sam_cmd
	
    #print >> sys.stderr, "[%s] Merging all sams file" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: reads all chromo sizes failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: read_chromo_size not found on this system"
        exit(1)
	
    return chromo_size_file


def format_reads_func(reads_file, formated_reads_file):
    print >> sys.stderr, "[%s] format reads" % (right_now())
    
    mapsplice_log = open(logging_dir + "format_reads.log", "w")
    #stdout=mapsplice_log
   
    merge_sam_cmd = [bin_dir + "remove_blankspace_perline"]  
    
    merge_sam_cmd = merge_sam_cmd + [reads_file]
    
    merge_sam_cmd = merge_sam_cmd + [formated_reads_file]
    
    #print >> sys.stderr, "[%s] Merging all sams file" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: format reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: remove_blankspace_perline not found on this system"
        exit(1)
	
    return formated_reads_file

def sam2fq(sam_formatted, converted_read_file, flag, log_file):
    print >> sys.stderr, "[%s] convert sam to fq reads" % (right_now())
    
    mapsplice_log = open(log_file, "w")
   
    merge_sam_cmd = [bin_dir + "sam2fq"]  

    merge_sam_cmd = merge_sam_cmd + [converted_read_file];
    
    merge_sam_cmd = merge_sam_cmd + [flag];
    
    read_files_array = sam_formatted.split(',')

    for read_file in read_files_array:
	merge_sam_cmd = merge_sam_cmd + [read_file];

    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % merge_sam_cmd
	
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam to fq reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2fq not found on this system"
        exit(1)
	
    return converted_read_file

def fqfa2sam(samfile, reads_files, flag, log_file):
    print >> sys.stderr, "[%s] convert unmapped reads to sam" % (right_now())
    
    mapsplice_log = open(log_file, "w")
   
    merge_sam_cmd = [bin_dir + "reads2unmappedsam"]  

    merge_sam_cmd = merge_sam_cmd + [samfile];
    
    merge_sam_cmd = merge_sam_cmd + [str(flag)];
    
    for read_file in reads_files:
	merge_sam_cmd = merge_sam_cmd + [read_file];

    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % merge_sam_cmd
	
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert unmapped reads to sam failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: reads2unmappedsam not found on this system"
        exit(1)
	
    return samfile

def gtf2exons(genegtf, annotation, log_file):
    print >> sys.stderr, "[%s] convert gene gtf to annotation format" % (right_now())
    
    mapsplice_log = open(log_file, "w")
   
    merge_sam_cmd = [bin_dir + "gtf2genetab"]  

    merge_sam_cmd = merge_sam_cmd + [genegtf] + [annotation];
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % merge_sam_cmd
	
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert gene gtf to annotation format failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: gtf2genetab not found on this system"
        exit(1)
	
    return annotation

def sam2fq_array(sam_formatted, converted_read_file, flag, log_file):
    print >> sys.stderr, "[%s] convert sam to fq reads" % (right_now())
    
    mapsplice_log = open(log_file, "w")
   
    merge_sam_cmd = [bin_dir + "sam2fq"]  

    merge_sam_cmd = merge_sam_cmd + [converted_read_file];
    
    merge_sam_cmd = merge_sam_cmd + [flag];
    
    read_files_array = sam_formatted.split(',')

    for read_file in read_files_array:
	merge_sam_cmd = merge_sam_cmd + [read_file];

    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % merge_sam_cmd
	
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert sam to fq reads failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: sam2fq not found on this system"
        exit(1)
	
    return converted_read_file

def bam2sam(bam_file, sam_file):
    print >> sys.stderr, "[%s] convert bam to sam" % (right_now())
    
    if os.path.exists(sam_file) and \
       rerun_all == 0:
	return sam_file        

    mapsplice_log = open(sam_file, "w")
   
    merge_sam_cmd = [bin_dir + "samtools"]  

    merge_sam_cmd = merge_sam_cmd + ["view"] + ["-h"];

    merge_sam_cmd = merge_sam_cmd + [bam_file];
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s] " % merge_sam_cmd
	
    
	
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert bam to sam failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: samtools not found on this system"
        exit(1)
	
    return sam_file

def format_chromos_func(all_chroms_path):
    print >> sys.stderr, "[%s] format chromosome files" % (right_now())
    
    #print >> sys.stderr, "[%s] " % all_chroms_path
   
    merge_sam_cmd = [bin_dir + "remove_blankspace"]  
    
    mapsplice_log = open(logging_dir + "format_chromos.log", "w")
    #stdout=mapsplice_log
    
    for sam_file in all_chroms_path:
	merge_sam_cmd = merge_sam_cmd + [sam_file]
    
    merge_sam_cmd = merge_sam_cmd + [formated_chrom_dir]
    
    #print >> sys.stderr, "[%s] format chromosomes" % merge_sam_cmd
    
    #print >> sys.stderr, "[%s] Merging all sams file" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: format chromosomes failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: remove_blankspace not found on this system"
        exit(1)
	
    return formated_chrom_dir

def merge_sam_cut_alters(unspliced_sam, spliced_sam, merged_sam):
    print >> sys.stderr, "[%s] Merging bowtie sam file cut alters" % (right_now())

    if os.path.exists(merged_sam) and \
	   rerun_all == 0:
	return merged_sam
    
    merge_sam_cmd = [bin_dir + "merge_sam_cutalters", unspliced_sam, spliced_sam, merged_sam]  
    try:    
        retcode = subprocess.call(merge_sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: merge bowtie cut alters failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: merge_sam_cutalters not found on this system"
        exit(1)
	
    return merged_sam

def pileup2wig(pileupfile, wigfile, allchromos_fai):
    print >> sys.stderr, "[%s] convert pileup to wig file" % (right_now())

    if os.path.exists(wigfile) and \
	   rerun_all == 0:
	temp_fs = open(pileupfile, "w")
	temp_fs.close()
	return wigfile
    
    merge_sam_cmd = [bin_dir + "pileup2wig", pileupfile, wigfile, allchromos_fai]  
    try:    
        retcode = subprocess.call(merge_sam_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert pileup to wig file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: pileup2wig not found on this system"
        exit(1)
    
    temp_fs = open(pileupfile, "w")
    temp_fs.close()
    
    return wigfile

def sam2pileup(bowtie_mapped_sam, allchromos_file):
    print >> sys.stderr, "[%s] Converting bowtie sam file to pileup format" % (right_now())

    pileup_file = bowtie_mapped_sam + ".pileup"
    
    bowtie_mapped_bam = bowtie_mapped_sam + ".bam"
    
    bowtie_mapped_bam_sorted = bowtie_mapped_bam + ".sorted"
    
    bowtie_mapped_bam_sorted_bam = bowtie_mapped_bam_sorted + ".bam"
    
    allchromos_file_fai = allchromos_file + ".fai"
    
    if os.path.exists(pileup_file) and \
	   rerun_all == 0:
	temp_fs = open(bowtie_mapped_sam, "w")
	temp_fs.close()
	if os.path.exists(bowtie_mapped_bam):
	    temp_fs = open(bowtie_mapped_bam, "rw")
	    temp_fs.close()
	if os.path.exists(bowtie_mapped_bam_sorted_bam):
	    temp_fs = open(bowtie_mapped_bam_sorted_bam, "rw")
	    temp_fs.close()
	return pileup_file
    
    pileup_fs = open(pileup_file, "w")
    

    
    sam2pileup_cmd = ""
    if bowtie_mapped_sam.endswith(".sam"):
	sam2faidx_cmd = [bin_dir + "samtools", "faidx", allchromos_file]
	
	try:    
	    retcode = subprocess.call(sam2faidx_cmd)
	    
	    if retcode != 0:
		print >> sys.stderr, fail_str, "Error: faidx allchromos_file failed"
		exit(1)
		
	except OSError, o:
	    if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
		print >> sys.stderr, fail_str, "Error: samtools not found on this system"
	    exit(1)
	    
	
	
    	sam2import_cmd = [bin_dir + "samtools", "import", allchromos_file_fai,
                      	bowtie_mapped_sam, bowtie_mapped_bam]
	
	try:    
	    retcode = subprocess.call(sam2import_cmd)
       
	    if retcode != 0:
		print >> sys.stderr, fail_str, "Error: import sam to bam failed"
		exit(1)
	
	except OSError, o:
	    if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
		print >> sys.stderr, fail_str, "Error: samtools not found on this system"
	    exit(1)
	    
	temp_fs = open(bowtie_mapped_sam, "w")
	temp_fs.close()

	sam2sort_cmd = [bin_dir + "samtools", "sort", bowtie_mapped_bam, bowtie_mapped_bam_sorted]
	
	try:    
	    retcode = subprocess.call(sam2sort_cmd)
       
	    if retcode != 0:
		print >> sys.stderr, fail_str, "Error: sort bam failed"
		exit(1)

	except OSError, o:
	    if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
		print >> sys.stderr, fail_str, "Error: samtools not found on this system"
	    exit(1)
	    
	temp_fs = open(bowtie_mapped_bam, "w")
	temp_fs.close()
	
    elif bowtie_mapped_sam.endswith(".bam"):
	bowtie_mapped_bam_sorted_bam = bowtie_mapped_sam
	
    
    sam2pileup_cmd = [bin_dir + "samtools", "pileup",
                     "-f", allchromos_file,
                     bowtie_mapped_bam_sorted_bam, ">", pileup_file]
    try:    
        retcode = subprocess.call(sam2pileup_cmd, stdout=pileup_fs)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: convert bowtie sam to pileup failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: samtools not found on this system"
        exit(1)
	
    temp_fs = open(bowtie_mapped_bam_sorted_bam, "w")
    temp_fs.close()
    
    return pileup_file


def sam2bam(bowtie_mapped_sam, bam_file, allchromos_file_fai):
    print >> sys.stderr, "[%s] Converting bowtie sam file to bam format" % (right_now())

    bowtie_mapped_bam = bam_file
    
    #allchromos_file_fai = allchromos_file;#allchromos_file + ".fai"
    
    if os.path.exists(bowtie_mapped_bam) and \
	   rerun_all == 0:
    
	return bowtie_mapped_bam
    
    sam2pileup_cmd = ""
    if bowtie_mapped_sam.endswith(".sam"):
	#sam2faidx_cmd = [bin_dir + "samtools", "faidx", allchromos_file]
	
	#if DEBUG == 1:
	    #print >> sys.stderr, "[%s] " % sam2faidx_cmd
	
	#try:    
	    #retcode = subprocess.call(sam2faidx_cmd)
	    
	    #if retcode != 0:
		#print >> sys.stderr, fail_str, "Error: faidx allchromos_file failed"
		#exit(1)
		
	#except OSError, o:
	    #if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
		#print >> sys.stderr, fail_str, "Error: samtools not found on this system"
	    #exit(1)
	
    	#sam2import_cmd = [bin_dir + "samtools", "import", allchromos_file_fai,
                      	#bowtie_mapped_sam, bowtie_mapped_bam]
	
	sam2import_cmd = [bin_dir + "samtools", "view", "-S", "-b", "-o",
                      	bowtie_mapped_bam, bowtie_mapped_sam]
	
	if DEBUG == 1:
	    print >> sys.stderr, "[%s] " % sam2import_cmd
	    
	try:    
	    retcode = subprocess.call(sam2import_cmd)
       
	    if retcode != 0:
		print >> sys.stderr, fail_str, "Error: import sam to bam failed"
		exit(1)
	
	except OSError, o:
	    if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
		print >> sys.stderr, fail_str, "Error: samtools not found on this system"
	    exit(1)
	
    elif bowtie_mapped_sam.endswith(".bam"):
	bowtie_mapped_bam = bowtie_mapped_sam
	    
    return bowtie_mapped_bam


def synthetic_junc(junc_file, chromosome_files_dir, read_width, anchor_width):
    print >> sys.stderr, "[%s] Synthetic junctions" % (right_now())
    
    chromosome_files_dir = chromosome_files_dir + "/";
    
    synjunc_log = open(logging_dir + "synjunc.log", "w")
    
    syn_junction = output_dir + "syn_junctions.txt";
    
    syn_junc_cmd = [bin_dir + "synjunc", junc_file, 
                        syn_junction, str(read_width - anchor_width), chromosome_files_dir]            
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=synjunc_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Synthetic junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: synjunc not found on this system"
        exit(1)
    return syn_junction

def syn_fusion_junc_seq(junc_file, syn_junc_file, chromosome_files_dir, synlen, log_file):
    print >> sys.stderr, "[%s] Synthetic fusion junctions sequence" % (right_now())
    
    if os.path.exists(syn_junc_file) and \
	   rerun_all == 0:
	return syn_junc_file
    
    chromosome_files_dir = chromosome_files_dir + "/";
    
    synjunc_log = open(log_file, "w");
    
    syn_junction = syn_junc_file;
    
    syn_junc_cmd = [bin_dir + "syn_fusion_junc_seq", junc_file, 
                        syn_junc_file, chromosome_files_dir, str(synlen)]            
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=synjunc_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Synthetic fusion junctions sequence"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: syn_fusion_junc_seq not found on this system"
        exit(1)
    return syn_junction

def FusionBWA2FusionSam(fusion_BWT_file, fusion_SAM_file, synlen, log_file):
    print >> sys.stderr, "[%s] Convert fusion bwt to sam" % (right_now())
    
    if os.path.exists(fusion_SAM_file) and \
	   rerun_all == 0:
	return fusion_SAM_file
    
    synjunc_log = open(log_file, "w");

    syn_junc_cmd = [bin_dir + "FusionBWA2FusionSam", fusion_BWT_file, 
                        fusion_SAM_file, str(synlen)]            
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=synjunc_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Convert fusion bwt to sam"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FusionBWA2FusionSam not found on this system"
        exit(1)
    return fusion_SAM_file

def FusionBWA2FusionSamNew(fusion_BWT_file, fusion_SAM_file, synlen, log_file):
    print >> sys.stderr, "[%s] Convert fusion bwt to sam" % (right_now())
    
    if os.path.exists(fusion_SAM_file) and \
	   rerun_all == 0:
	return fusion_SAM_file
    
    synjunc_log = open(log_file, "w");

    syn_junc_cmd = [bin_dir + "FusionBWA2FusionSam", fusion_BWT_file, 
                        fusion_SAM_file, str(synlen)]            
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=synjunc_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Convert fusion bwt to sam"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FusionBWA2FusionSam not found on this system"
        exit(1)
    return fusion_SAM_file

def FilterSamByJunc(all_sams_path, junc_bed, remained_sam, filtered_sam, log_file):
    print >> sys.stderr, "[%s] Filter Sam By junction" % (right_now())
    
    if os.path.exists(remained_sam) and \
       os.path.exists(filtered_sam) and \
	   rerun_all == 0:
	return (remained_sam, filtered_sam)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    merge_sam_cmd = [bin_dir + "FilterSamByJunc"] + [junc_bed] + [remained_sam] + [filtered_sam]
   
    for sam_file in all_sams_path:
	merge_sam_cmd = merge_sam_cmd + [sam_file]
   
    #print >> sys.stderr, "[%s] Merging all sams file" % merge_sam_cmd
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filter Sam By junction failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: FilterSamByJunc not found on this system"
        exit(1)
	
    return (remained_sam, filtered_sam)

def AddTagsToSam(sam_file, added_tags_sam_file, paired, reads_file, unmapped_sam, stat_file, format_flag, max_insert, junction, chrom_size, log_file):
    print >> sys.stderr, "[%s] Add tags to sam file" % (right_now())
    
    if os.path.exists(added_tags_sam_file) and \
       os.path.exists(unmapped_sam) and \
	   rerun_all == 0:
	return (added_tags_sam_file, unmapped_sam)
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    merge_sam_cmd = [bin_dir + "AddTagsToSam"] + [sam_file] + [added_tags_sam_file] + [stat_file]

    ispaired = 0;
    
    if paired != "":
	ispaired = 1;
	
    cmd_format = "fa"
    
    if format_flag == '-q':
	cmd_format = "fq"
	
    merge_sam_cmd = merge_sam_cmd + [str(ispaired)] + [reads_file] + [unmapped_sam] + [cmd_format] + [str(max_insert)] + [junction] + [chrom_size]
   
    if DEBUG == 1:
	print >> sys.stderr, "[%s] AddTagsToSam" % merge_sam_cmd
	
    try:    
        retcode = subprocess.call(merge_sam_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Add tags to sam file failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: AddTagsToSam not found on this system"
        exit(1)
	
    return (added_tags_sam_file, unmapped_sam)

def junc_db(junc_file, chromosome_files_dir, min_anchor_width, max_anchor, max_threshold, synthetic_file, log_file):
    print >> sys.stderr, "[%s] Synthetic junctions sequence" % (right_now())
    
    if os.path.exists(synthetic_file) and \
	   rerun_all == 0:
	return synthetic_file
    
    chromosome_files_dir = chromosome_files_dir + "/";
    
    syn_junc_cmd = [bin_dir + "junc_db", str(min_anchor_width), str(max_anchor), str(max_threshold),
                        junc_file, chromosome_files_dir, synthetic_file]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % syn_junc_cmd
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Synthetic junctions sequence failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: junc_db not found on this system"
        exit(1)
    return synthetic_file

def combine_junc(junc_files, combined_junc, log_file):
    print >> sys.stderr, "[%s] combine junctions" % (right_now())
    
    if os.path.exists(combined_junc) and \
	   rerun_all == 0:
	return combined_junc
    
    syn_junc_cmd = [bin_dir + "comb_junc"];
    
    splitted_reads2 = junc_files.split(',');
    
    for read_file in splitted_reads2:
	syn_junc_cmd = syn_junc_cmd + [read_file];
		
    syn_junc_cmd = syn_junc_cmd + [combined_junc];
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % syn_junc_cmd
    
    mapsplice_log = open(log_file, "w")
    
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: combine junctions failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: comb_junc not found on this system"
        exit(1)
    return combined_junc

def fusion_junc_db(junc_file, fusion_junc_file, chromosome_files_dir, min_anchor_width, 
                   max_anchor, max_threshold_each, max_threshold_total, synthetic_file, log_file):
    print >> sys.stderr, "[%s] Synthetic fusion junctions sequence" % (right_now())
    
    if os.path.exists(synthetic_file) and \
	   rerun_all == 0:
	return synthetic_file
    
    chromosome_files_dir = chromosome_files_dir + "/";
    
    syn_junc_cmd = [bin_dir + "junc_db_fusion", str(min_anchor_width), str(max_anchor), str(max_threshold_each), str(max_threshold_total),
                        junc_file, "1", fusion_junc_file, "1", chromosome_files_dir, synthetic_file]
    
    if DEBUG == 1:
	print >> sys.stderr, "[%s]" % syn_junc_cmd
    
    mapsplice_log = open(log_file, "w")
    #stdout=mapsplice_log
    
    try:    
        retcode = subprocess.call(syn_junc_cmd, stdout=mapsplice_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Synthetic fusion junctions sequence failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: junc_db_fusion not found on this system"
        exit(1)
    return synthetic_file

def compare_junc(f1_junc, f2_junc, compare_range, cur_out_dir):
    print >> sys.stderr, "[%s] compare juncitons" % (right_now())
   
    basename = os.path.basename(f1_junc)
    filename = os.path.splitext(basename)
    f1_basename = filename[0]
    
    basename = os.path.basename(f2_junc)
    filename = os.path.splitext(basename)
    f2_basename = filename[0]
    
    if os.path.exists(cur_out_dir):
	pass
    else:        
	os.mkdir(cur_out_dir)
    
    in_f1_and_in_f2 = cur_out_dir + "in_(" + f1_basename + ")_in_(" + f2_basename + ").txt";
    in_f1_not_in_f2 = cur_out_dir + "in_(" + f1_basename + ")_NOTin_(" + f2_basename + ").txt";
    in_f2_not_in_f1 = cur_out_dir + "in_(" + f2_basename + ")_NOTin_(" + f1_basename + ").txt";
    
    compare_junc_cmd = ["mono", bin_dir + "compare_junc.exe", "-h1", "-h2", "-o", cur_out_dir, "-s", "-m", str(compare_range), 
			      f1_junc, f2_junc]            
    try:    
        retcode = subprocess.call(compare_junc_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: compare juncitons failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error:compare_junc.exe not found on this system"
        exit(1)
    return (in_f1_not_in_f2, in_f2_not_in_f1, in_f1_and_in_f2)

def filter_by_expression_level(junc_file, gamma, delta):
    print >> sys.stderr, "[%s] Filtering junction by expression level format" % (right_now())

    filtered_junction = output_dir + "junctions_filtered.txt";
    
    fel_log = open(logging_dir + "filter_byexpresslevel.log", "w")
    
    filter_50bp_cmd = [bin_dir + "filterjuncbyexpresslevel", junc_file, 
                       filtered_junction, str(gamma), str(delta), output_dir]            
    try:    
        retcode = subprocess.call(filter_50bp_cmd, stdout=fel_log)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: Filtering junction by expression level failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterbyexpresslevel not found on this system"
        exit(1)
    return filtered_junction

def filterbyanchor(junction, min_anchor, read_width, filtered, notfiltered):
    print >> sys.stderr, "[%s] filter by anchor length" % (right_now())
    
    if os.path.exists(filtered) and \
           os.path.exists(notfiltered)  and \
	   rerun_all == 0:
	return (filtered, notfiltered)
        
    filterbyintronlen_cmd = [bin_dir + "filterjuncbyanchor", junction, filtered, notfiltered, str(min_anchor), str(read_width)]            
    try:    
        retcode = subprocess.call(filterbyintronlen_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by anchor length failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filterjuncbyanchor not found on this system"
        exit(1)
	
    return (filtered, notfiltered)

def filternotinisland(island, junction, notinisland, inisland):
    print >> sys.stderr, "[%s] filter by island" % (right_now())
    
    if os.path.exists(notinisland) and \
           os.path.exists(inisland)  and \
	   rerun_all == 0:
	return (notinisland, inisland)
    
    filternotinisland_cmd = [bin_dir + "filternotinisland", island, junction, 
			     notinisland, inisland]
    try:    
        retcode = subprocess.call(filternotinisland_cmd)
       
        if retcode != 0:
            print >> sys.stderr, fail_str, "Error: filter by island failed"
            exit(1)
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: filternotinislands not found on this system"
        exit(1)
	
    return (notinisland, inisland)

def process_extsam(samfile, outputjuncname, cur_out_dir, min_intron_length, max_intron_length, min_anchor, chromosome_files_dir, read_width):
    print >> sys.stderr, "[%s] process external sam file" % (right_now())
    
    sam2junc1(samfile, cur_out_dir + outputjuncname, 
	 chromosome_files_dir, read_width, 1, 100000000, "", 1)
    
    sam2junc1(samfile, cur_out_dir + outputjuncname + ".in_intron", 
	 chromosome_files_dir, read_width, min_intron_length, max_intron_length, "", 1)
    
    separate_canon_noncanon(cur_out_dir + outputjuncname + ".in_intron",
			    cur_out_dir + outputjuncname + ".in_intron" + ".canon",
			    cur_out_dir + outputjuncname + ".in_intron" + ".noncanon")
    
    sam2junc1(samfile, cur_out_dir + outputjuncname + ".exceed_intron", 
	 chromosome_files_dir, read_width,max_intron_length+1, 100000000, "", 1)
    
    separate_canon_noncanon(cur_out_dir + outputjuncname + ".exceed_intron",
			    cur_out_dir + outputjuncname + ".exceed_intron" + ".canon",
			    cur_out_dir + outputjuncname + ".exceed_intron" + ".noncanon")
    
    sam2junc1(samfile, cur_out_dir + outputjuncname + ".in_intron.anchor", 
	 chromosome_files_dir, read_width, min_intron_length, max_intron_length, "", min_anchor)
    
    separate_canon_noncanon(cur_out_dir + outputjuncname + ".in_intron.anchor",
			    cur_out_dir + outputjuncname + ".in_intron.anchor" + ".canon",
			    cur_out_dir + outputjuncname + ".in_intron.anchor" + ".noncanon")
    
    sam2junc1(samfile, cur_out_dir + outputjuncname + ".exceed_intron.anchor", 
	 chromosome_files_dir, read_width,max_intron_length+1, 100000000, "", min_anchor)
    
    separate_canon_noncanon(cur_out_dir + outputjuncname + ".exceed_intron.anchor",
			    cur_out_dir + outputjuncname + ".exceed_intron.anchor" + ".canon",
			    cur_out_dir + outputjuncname + ".exceed_intron.anchor" + ".noncanon")
    
    return (cur_out_dir + outputjuncname + ".in_intron" + ".canon",
	    cur_out_dir + outputjuncname + ".in_intron" + ".noncanon",
	    cur_out_dir + outputjuncname + ".exceed_intron" + ".canon",
	    cur_out_dir + outputjuncname + ".exceed_intron" + ".noncanon",
	    cur_out_dir + outputjuncname + ".in_intron.anchor" + ".canon",
	    cur_out_dir + outputjuncname + ".in_intron.anchor" + ".noncanon",
	    cur_out_dir + outputjuncname + ".exceed_intron.anchor" + ".canon",
	    cur_out_dir + outputjuncname + ".exceed_intron.anchor" + ".noncanon")
	    
def print_arguments(argu_file):
    
    
    print >> sys.stderr, "print argugment"
    
    argu_log = open(argu_file, "w")
    
    print >> argu_log, "min_anchor_length=[%s]" % (min_anchor_length)
    
    print >> argu_log, "seed_length=[%s]" % (seed_length)

    print >> argu_log, "splice_mismatches=[%s]" % (splice_mismatches)

    print >> argu_log, "segment_mismatches=[%s]" % (segment_mismatches)

    print >> argu_log, "FASTA_file_extension=[%s]" % (segment_mismatches)

    print >> argu_log, "read_file_suffix=[%s]" % (read_file_suffix)

    print >> argu_log, "min_intron_length=[%s]" % (min_intron_length)

    print >> argu_log, "max_intron_length=[%s]" % (max_intron_length)

    print >> argu_log, "island_extension=[%s]" % (max_intron_length)
    
    print >> argu_log, "read_width=[%s]" % (read_width)

    print >> argu_log, "rank=[%s]" % (rank)

    print >> argu_log, "flank_case=[%s]" % (flank_case)

    print >> argu_log, "fusion_flank_case=[%s]" % (fusion_flank_case)

    print >> argu_log, "islands_file=[%s]" % (islands_file)

    print >> argu_log, "read_files_dir=[%s]" % (read_files_dir)

    print >> argu_log, "chromosome_files_dir=[%s]" % (chromosome_files_dir)
    
    print >> argu_log, "all_chromosomes_file=[%s]" % (all_chromosomes_file)

    print >> argu_log, "repeat_regioins=[%s]" % (repeat_regioins)
    
    print >> argu_log, "gene_regions=[%s]" % (gene_regions)
    
    print >> argu_log, "bwt_idx_prefix=[%s]" % (bwt_idx_prefix)
    
    print >> argu_log, "bowtie_threads=[%s]" % (bowtie_threads)
    
    print >> argu_log, "max_hits=[%s]" % (max_hits)
    
    print >> argu_log, "threshold=[%s]" % (threshold)
    
    print >> argu_log, "boundary=[%s]" % (boundary)

    print >> argu_log, "num_anchor=[%s]" % (num_anchor)
    
    print >> argu_log, "unmapped_reads=[%s]" % (unmapped_reads)

    print >> argu_log, "sam_formatted=[%s]" % (sam_formatted)

    print >> argu_log, "bam_file=[%s]" % (bam_file)
    
    print >> argu_log, "sam_formatted_25=[%s]" % (sam_formatted_25)
 
    print >> argu_log, "bwt_map_25=[%s]" % (bwt_map_25)
    
    print >> argu_log, "pileup_file=[%s]" % (pileup_file)
    
    print >> argu_log, "synthetic_mappedreads=[%s]" % (synthetic_mappedreads)
    
    print >> argu_log, "tophat_mappedreads=[%s]" % (tophat_mappedreads)
    
    print >> argu_log, "pairend=[%s]" % (pairend)
    
    print >> argu_log, "gamma=[%s]" % (gamma)
    
    print >> argu_log, "delta=[%s]" % (delta)
    
    #print >> argu_log, "num_seg=[%s]" % (num_seg)

    print >> argu_log, "seg_len=[%s]" % (seg_len)

    print >> argu_log, "fix_hole_file=[%s]" % (fix_hole_file)
    
    print >> argu_log, "format_flag=[%s]" % (format_flag)
    
    print >> argu_log, "chrom_size_file=[%s]" % (chrom_size_file)
    
    print >> argu_log, "extend_bits=[%s]" % (extend_bits)
    
    print >> argu_log, "total_fusion_mismatch=[%s]" % (total_fusion_mismatch)
    
    print >> argu_log, "total_mismatch=[%s]" % (total_mismatch)
    
    print >> argu_log, "append_mismatch=[%s]" % (append_mismatch)
    
    print >> argu_log, "remap_mismatch=[%s]" % (remap_mismatch)
    
    print >> argu_log, "skip_bwt=[%s]" % (skip_bwt)
    
    print >> argu_log, "prefix_match=[%s]" % (prefix_match)
    
    print >> argu_log, "fullrunning=[%s]" % (fullrunning)
    
    print >> argu_log, "collect_stat=[%s]" % (collect_stat)
    
    print >> argu_log, "rm_temp=[%s]" % (rm_temp)
    
    print >> argu_log, "format_reads=[%s]" % (format_reads)
    
    print >> argu_log, "format_chromos=[%s]" % (format_chromos)
    
    print >> argu_log, "do_fusion=[%s]" % (do_fusion)
    
    print >> argu_log, "do_cluster=[%s]" % (do_cluster)
    
    print >> argu_log, "search_whole_chromo=[%s]" % (search_whole_chromo)
    
    print >> argu_log, "map_segment_directly=[%s]" % (map_segment_directly)
    
    print >> argu_log, "run_mapper=[%s]" % (run_mapper)
    
    print >> argu_log, "max_insert=[%s]" % (max_insert)
    
    print >> argu_log, "min_missed_seg=[%s]" % (min_missed_seg)
    
    print >> argu_log, "do_annot_gene=[%s]" % (do_annot_gene)
    
    print >> argu_log, "do_annot_gene=[%s]" % (annot_gene_file)
    
    print >> argu_log, "do_filter_fusion_by_repeat=[%s]" % (do_filter_fusion_by_repeat)
    
def main(argv=None):
    
    #test = "1,2,3";
    
    #test_split = test.split(',');
    
    #combined_reads = "";
    #for ix in range(len(test_split)):
	#combined_reads = combined_reads + test_split[ix]
	
	#if ix < len(test_split) - 1:
	    #combined_reads = combined_reads + ',';
	
	
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hvo:s:n:m:i:x:w:S:R:B:p:t:u:c:f:a:d:g:T:D:U:M:P:N:X:G:L:H:E:Y:I:Q:A:F:C:e:K:r:O:1:2:l:k:", 
                                        ["version",
                                         "help",  
                                         "seed-length=",
                                         "splice-mis=",
					 "segment-mismatches=",
					 "read-file-suffix=",
                                         "min-anchor-length=",
                                         "min-intron=",
                                         "max-intron=",
                                         "extend-exons=",
                                         "read-width=", 
                                         "non-canonical",
					 "semi-canonical",
                                         "fusion-non-canonical",
					 "fusion-semi-canonical",
                                         "delta=",
                                         "gamma=",
                                         "threshold=",
                                         "boundary=",
                                         "Rank=",
                                         "bowtie-index=",
					 "threads=",
                                         "islands-file=",
                                         "reads-file=",
                                         "chromosome-dir=",
                                         "all-chromosomes-files=",
                                         "FASTA-files=",
                                         "unmapped-reads=",
                                         "Bam-file=",
					 "full-running",
					 "num-anchor=",
                                         "pileup-file=",
					 "numseg=",
					 "seglen=",
					 "fixholefile=",
					 "synthetic=",
					 "tophat=",
					 "pairend",
					 "fastq=",
					 "extend-bits=",
					 "total-mismatch=",
					 "total-fusion-mismatch=",
                                         "output-dir=",
					 "chrom-size=",
					 "skip-bowtie=",
					 "prefix-match=",
					 "collect-stat",
					 "max-hits=",
					 "keep-tmp",
					 "not-rerun-all",
					 "format-chromos",
					 "format-reads",
					 "fusion",
					 "cluster",
					 "DEBUG",
					 "run-MapPER",
	                                 "split-index=",
					 "search-whole-chromosome",
					 "map-segments-directly",
					 "remap-mismatches=",
					 "avoid-regions=",
					 "config=",
					 "interested-regions=",
					 "ins=",
	                                 "del=",
	                                 "bam",
					 "min-missed-seg=",
	                                 "min-map-len=",
	                                 "min-len=",
	                                 "min-entropy=",
	                                 "min-fusion-entropy=",
	                                 "qual-scale=",
	                                 "annotgene=",
	                                 "filter-fusion-by-repeat=",
	                                 "end_1=",
	                                 "end_2=",
	                                 "gene-gtf=",
	                                 "Sorted-bam",
					 "max-append-mis="])
        except getopt.error, msg:
            raise Usage(msg)
        
        min_anchor_length = 8
	seed_length = 10
        splice_mismatches = 1
	segment_mismatches = 1
	FASTA_file_extension = "fa"
	read_file_suffix = "txt"
        min_intron_length = 50
        max_intron_length = 200000
        island_extension = 0
        read_width = 0
	min_read_len = 25
        rank = 0.0
        flank_case = 5
        fusion_flank_case = 5
        islands_file = ""
        read_files_dir = ""
        chromosome_files_dir = ""
        all_chromosomes_file = ""
	repeat_regioins = ""
	gene_regions = ""
	qual_scale = "phred64"
        
        bwt_idx_prefix = ""
        bowtie_threads = 1
	max_hits = 4
        threshold = 1
        boundary = 50

        num_anchor = 1
        
        unmapped_reads = ""
        sam_formatted = ""
	bam_file = ""
	
        sam_formatted_25 = ""
	bwt_map_25 = ""
	
        pileup_file = ""
	
	synthetic_mappedreads = ""
	
	tophat_mappedreads = ""
	
	pairend = ""
        
        gamma = 0.1
        delta = 0.1
	
	num_seg = 4
	
	seg_len = 25

	fix_hole_file = "" 
	
	format_flag = ""
	
	chrom_size_file = ""
	
	extend_bits = 3;
	
	total_fusion_mismatch = 2;
	
	total_mismatch = 2;
	
	append_mismatch = 3;
	
	remap_mismatch = 2;
	
	skip_bwt = 0;
	
	prefix_match = 1;
	
	fullrunning = 0;
	
	collect_stat = 0;
	
	rm_temp = 1;
	
	format_reads = 0;
	
	format_chromos = 0;
	
	do_fusion = 0;
	
	do_cluster = 0;
	
	search_whole_chromo = 0;
	
	map_segment_directly = 0;
	
	run_mapper = 0;
	
	max_insert = 6;
	
	max_delete = 6;
	
	min_missed_seg = 0;
	
	min_map_len = 0;
	
	do_annot_gene = 0;
	
	annot_gene_file = "";
	
	do_filter_fusion_by_repeat = 0;
	
	chrom_blat = ""
	
	junction_tobe_syn = "";
	
	input_reads_1 = "";
	
	input_reads_2 = "";
	
	split_index = "";
	
	min_entropy = -0.0001;
	
	gene_gtf_file = "";
	
	min_entropy_repeats = 1;
	
	dec_sam_file = "";
	
	issortedbam = 0;
	
	global rerun_all
	
	global DEBUG
	
	if len(args) == 1:
	    params = Params()
	    
	    params.parse_cfgfile(args[0]);
	    
	    min_anchor_length = params.min_anchor_length
	    seed_length = params.seed_length
	    splice_mismatches = params.splice_mismatches
	    segment_mismatches = params.segment_mismatches
	    FASTA_file_extension = params.FASTA_file_extension
	    read_file_suffix = params.read_file_suffix
	    min_intron_length = params.min_intron_length
	    max_intron_length = params.max_intron_length
	    island_extension = params.island_extension
	    read_width = params.read_width
	    rank = params.rank
	    flank_case = params.flank_case
	    fusion_flank_case = params.fusion_flank_case
	    islands_file = params.islands_file
	    read_files_dir = params.read_files_dir
	    chromosome_files_dir = params.chromosome_files_dir
	    all_chromosomes_file = params.all_chromosomes_file
	    repeat_regioins = params.repeat_regioins
	    gene_regions = params.gene_regions	    
	    bwt_idx_prefix = params.bwt_idx_prefix
	    bowtie_threads = params.bowtie_threads
	    max_hits = params.max_hits
	    threshold = params.threshold
	    boundary = params.boundary    
	    num_anchor = params.num_anchor	    
	    unmapped_reads = params.unmapped_reads
	    sam_formatted = params.sam_formatted
	    sam_formatted_25 = params.sam_formatted_25
	    bwt_map_25 = params.bwt_map_25	    
	    pileup_file = params.pileup_file	    
	    synthetic_mappedreads = params.synthetic_mappedreads
	    tophat_mappedreads = params.tophat_mappedreads	    
	    pairend = params.pairend	    
	    gamma = params.gamma
	    delta = params.delta	    
	    #num_seg = params.num_seg
	    seg_len = params.seg_len    
	    fix_hole_file = params.fix_hole_file	    
	    format_flag = params.format_flag	    
	    chrom_size_file = params.chrom_size_file	    
	    extend_bits = params.extend_bits	    
	    total_fusion_mismatch = params.total_fusion_mismatch	    
	    total_mismatch = params.total_mismatch	    
	    append_mismatch = params.append_mismatch	    
	    remap_mismatch = params.remap_mismatch	    
	    skip_bwt = params.skip_bwt	    
	    prefix_match = params.prefix_match	    
	    fullrunning = params.fullrunning	    
	    collect_stat = params.collect_stat	    
	    rm_temp = params.rm_temp	    
	    format_reads = params.format_reads	    
	    format_chromos = params.format_chromos	    
	    do_fusion = params.do_fusion	    
	    do_cluster = params.do_cluster	    
	    search_whole_chromo = params.search_whole_chromo	    
	    map_segment_directly = params.map_segment_directly	    
	    run_mapper = params.run_mapper
	    max_insert = params.max_insert
	    bam_file = params.bam_file
	    min_missed_seg = params.min_missed_seg;
	    min_map_len = params.min_map_len;
	    do_annot_gene = params.do_annot_gene;
	    annot_gene_file = params.annot_gene_file;
	    do_filter_fusion_by_repeat = params.filter_fusion_by_repeat;
	    chrom_blat = params.chromosome_blat_idx
                              
        # option processing
	
	if len(opts) == 0:
	    raise Usage(use_message)
        for option, value in opts:
            if option in ("-o", "--output-dir"):
                global output_dir
                global logging_dir
		global canon_in_intron_dir
		global canon_exceed_intron_dir
		global noncanon_in_intron_dir
		global noncanon_exceed_intron_dir
		global fusion_dir
		global temp_dir
		global synthetic_dir
		global pairend_dir
		global original_dir
		global filteredbest_dir
		global comparison_dir
		global tophat_dir
		global remap_dir
		global remap_regions_dir
		global hmer_dir
		global hole_dir
		global head_dir
		global tail_dir
		global bwtout_dir
		global sam_dir
		global formated_chrom_dir
		global formated_reads_dir
		global fusion_dir
		global fusion_data_dir
		global fusion_data_single_dir
		global fusion_data_PER_dir
		global fusion_result_dir
		global fusion_result_PER_prob_dir
		global fusion_result_junction_support_dir
		global cluster_dir
		global cluster_result_dir
		global cluster_data_dir
		global cluster_data_parsedPER_dir
		global fusion_result_fusionRead_dir
		global filter_repeats_dir
		
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
		canon_in_intron_dir = output_dir + "canonical/"
		canon_exceed_intron_dir = output_dir + "canonical_exceed/"
		noncanon_in_intron_dir = output_dir + "noncanonical/"
		noncanon_exceed_intron_dir = output_dir + "noncanonical_exceed/"
		fusion_dir = output_dir + "fusion/"
		temp_dir = output_dir + "tmp/"
		
		synthetic_dir = temp_dir + "synthetic/"
		
		pairend_dir = temp_dir + "pairend/"
		
		original_dir = temp_dir + "original/"
		filteredbest_dir = temp_dir + "best/"
		comparison_dir = temp_dir + "comparison/"
		
		tophat_dir = temp_dir + "tophat/"
		
		remap_dir = temp_dir + "remap/"
		
		remap_regions_dir = temp_dir + "remap_regions/"
		
		hmer_dir = temp_dir + "single_anchored_middle/"
		hole_dir = temp_dir + "double_anchored/"
		head_dir = temp_dir + "single_anchored_head/"
		tail_dir = temp_dir + "single_anchored_tail/"
		bwtout_dir = temp_dir + "bwtout/"
		
		sam_dir = temp_dir + "sam/"
		
		fusion_dir = temp_dir + "fusion/"
		
		cluster_dir = temp_dir + "cluster/"
		
		cluster_result_dir = cluster_dir + "result/"
		cluster_data_dir = cluster_dir + "data/"
		cluster_data_parsedPER_dir = cluster_data_dir + "parsedPER/"
		
		fusion_data_dir = fusion_dir + "data/"
		fusion_data_single_dir = fusion_data_dir + "single/"
		fusion_data_PER_dir = fusion_data_dir + "PER/"
		fusion_result_dir = fusion_dir + "result/"
		fusion_result_PER_prob_dir = fusion_result_dir + "PER_prob/"
		fusion_result_junction_support_dir = fusion_result_dir + "junction_support/"
		fusion_result_fusionRead_dir = fusion_result_dir + "fusionRead/"
				
		formated_chrom_dir = temp_dir + "formated_chrom/"

		formated_reads_dir = temp_dir + "formated_reads/"
		
		filter_repeats_dir = temp_dir + "filter_repeats/"
		
            if option in ("-v", "--version"):
                raise Usage(ver_message) #% (get_version());
                #exit(0);
            if option in ("-h", "--help"):
                raise Usage(use_message)
            #if option in ("-f", "--config"):
                #a = 0;
            #if option in ("-w", "--read-width"):
                #read_width = int(value)
                #if not read_width >= 1:
                    #print >> sys.stderr, "Error: arg to --read-width must be greater than or equal to 1"
                    #exit(1)
            #if option in ("-n", "--min-anchor"):
                #min_anchor_length = int(value)
                #if min_anchor_length < 4:
                    #print >> sys.stderr, "Error: arg to --min-anchor-len must be greater than 3"
                    #exit(1)
            #if option in ("-N", "--num-anchor"):
                #num_anchor = int(value)
                #if num_anchor <= 0:
                    #print >> sys.stderr, "Error: arg to --num-anchor must be greater than 0"
                    #exit(1)
	    if option == "--ins":
                max_insert = int(value)
                if max_insert <= 0:
                    print >> sys.stderr, "Error: arg to --ins must be greater than 0"
                    exit(1)
	    if option == "--del":
                max_delete = int(value)
                if max_delete <= 0:
                    print >> sys.stderr, "Error: arg to --del must be greater than 0"
                    exit(1)	    
	    if option == "--max-append-mis":
                append_mismatch = int(value)
                if append_mismatch < 0:
                    print >> sys.stderr, "Error: arg to --max-append-mis must be greater or equal to 0"
                    exit(1)
            #if option in ("-s", "--read-file-suffix"):
                #read_file_suffix = value
            if option in ("-m", "--splice-mis"):
                splice_mismatches = int(value)
                if not splice_mismatches >= 0:
                    print >> sys.stderr, "Error: arg to --splice-mis must be greater than or equal to 0"
                    exit(1)
	    #if option in ("-e", "--extend-bits"):
                #extend_bits = int(value)
                #if not extend_bits >= 0:
                    #print >> sys.stderr, "Error: arg to --extend-bits must be greater than or equal to 0"
                    #exit(1)
	    #if option in ("-C", "--total-mismatch"):
                #total_mismatch = int(value)
                #if not total_mismatch >= 0:
                    #print >> sys.stderr, "Error: arg to --total-mismatch must be greater than or equal to 0"
                    #exit(1)
	    #if option in ("-F", "--total-fusion-mismatch"):
                #total_fusion_mismatch = int(value)
                #if not total_fusion_mismatch >= 0:
                    #print >> sys.stderr, "Error: arg to --total-fusion-mismatch must be greater than or equal to 0"
                    #exit(1)
	    #if option in ("-E", "--segment-mismatches"):
                #segment_mismatches = int(value)
                #if not segment_mismatches >= 0:
                    #print >> sys.stderr, "Error: arg to --segment-mismatches must be greater than or equal to 0"
                    #exit(1)
            if option in ("-i", "--min-intron"):
                min_intron_length = int(value)
                if min_intron_length <= 0:
                    print >> sys.stderr, "Error: arg to --min-intron must be greater than 0"
                    exit(1)                
            if option in ("-I", "--max-intron"):
                max_intron_length = int(value)
                if max_intron_length <= 0:
                    print >> sys.stderr, "Error: arg to --max-intron must be greater than 0"
                    exit(1)
	    if option in ("-p", "--threads"):
                bowtie_threads = int(value)
                if bowtie_threads <= 0:
                    print >> sys.stderr, "Error: arg to --threads must be greater than 0"
                    exit(1)
            #if option in ("-T", "--interested-regions"):
                #gene_regions = value
	    #if option in ("-G", "--numseg"):
                #num_seg = int(value)
		#
	    if option in ("-l", "--seglen"):
                seg_len = int(value)
	    #if option in ("-O"):
                #junction_tobe_syn = value
		##
	    #if option in ("-H", "--fixholefile"):
                #fix_hole_file = value
		#
	    
            #if option in ("-D", "--boundary"):
                #boundary = int(value)
                #if boundary <= 0:
                    #print >> sys.stderr, "Error: arg to --boundary must be greater than 0"
                    #exit(1)
		
            #if option in ("-p", "--extend-exons"):
                #island_extension = int(value)
                #if island_extension < 0:
                    #print >> sys.stderr, "Error: arg to --extend-exons must be at least 0"
                    #exit(1)
            if option == "--non-canonical":
                flank_case = 0
	    #if option == "--semi-canonical":
		#flank_case = 1
		
            if option == "--fusion-non-canonical":
                fusion_flank_case = 0
		do_fusion = 1
		
	    #if option == "--fusion-semi-canonical":
		#fusion_flank_case = 1
            #if option in ("-R", "--remap-mismatches"):
                #remap_mismatch =int(value)
                #if remap_mismatch < 0 or remap_mismatch > 3:
                    #print >> sys.stderr, "Error: arg to --remap-mismatch must be in [0, 3]"
                    #exit(1)
	    
            #if option in ("-d", "--delta"):
                #delta =float(value)
                #if delta < 0.0000001:
                    #print >> sys.stderr, "Error: arg to --delta must be at least 0.0000001"
                    #exit(1)
            #if option in ("-g", "--gamma"):
                #gamma =float(value)
                #if gamma < 0.0000001:
                    #print >> sys.stderr, "Error: arg to --gamma must be at least 0.0000001"
                    #exit(1)
            #if option in ("-S", "--FASTA-files-ext"):
                #FASTA_file_extension = value
            #if option in ("-a", "--all-chromosomes-files"):
                #all_chromosomes_file = value
            if option in ("-x", "--bowtie-index"):
                bwt_idx_prefix = value
            #if option in ("-u", "--reads-file"):
                #read_files_dir = value
	    #if option in ("-t", "--avoid-regions"):
                #repeat_regioins = value
            #if option in ("-U", "--unmapped-reads"):
                #unmapped_reads = value
            
            #if option in ("-P", "--pileup-file"):
                #pileup_file = value
	    #if option in ("-Y", "--synthetic"):
                #synthetic_mappedreads = value
	    #if option in ("-A", "--tophat"):
                #tophat_mappedreads = value
	    #if option in ("-Q", "--reads-format"):
		#if value == "fq":
		    #format_flag = "-q"
		#elif value == "fa":
		    #format_flag = "-f"
	    
	    #if option in ("-g","--skip-bowtie"):
		#skip_bwt = int(value)
	    #if option in ("-d","--prefix-match"):
		#prefix_match = int(value)
            if option in ("-c", "--chromosome-dir"):
                chromosome_files_dir = value
		chromosome_files_dir = chromosome_files_dir + "/"
	    #if option == "--collect-stat":
                #collect_stat = 1
	    #if option == "--chrom-size":
                #chrom_size_file = value
	    if option == "--bam":
                bam_file = "1" 
	    if option == "--Sorted-bam":
		issortedbam = 1;
	    
	    #if option == "--full-running":
                #fullrunning = 1
		#if fullrunning < 0:
                    #print >> sys.stderr, "Error: arg to --full-running must be greater than or equal to 0"
                    #exit(1)    	    		    
	    if option in ("-k", "--max-hits"):
		max_hits = int(value);
	    #if option == "--min-missed-seg":
		#min_missed_seg = int(value);
	    if option == "--min-map-len":
		min_map_len = int(value);
	    if option == "--min-len":
		min_read_len = int(value);
	    if option == "--gene-gtf":
		gene_gtf_file = value;
	    if option == "--Bam-file":
		dec_sam_file = value	    
	    if option == "--min-entropy":
		min_entropy = float(value);
	    if option == "--min-fusion-entropy":
		min_entropy_repeats = float(value);
	    if option == "--qual-scale":
		if value == "phred64":
		    qual_scale = value;
		elif value == "phred33":
		    qual_scale = value;
		elif value == "solexa64":
		    qual_scale = value;
		else:
		    print >> sys.stderr, "qual scale can only be phred64 or phred33 or solexa64"
                    exit(1)		
	    if option == "--keep-tmp":
		rm_temp = 0;
	    if option == "--not-rerun-all":
		rerun_all = 0
	    #if option == "--format-chromos":
		#format_chromos = 1
	    #if option == "--cluster":
		#do_cluster = 1
	    if option == "--fusion":
		do_fusion = 1
	    #if option == "--search-whole-chromosome":
		#search_whole_chromo = 1
	    #if option == "--map-segments-directly":
		#map_segment_directly = 1
	    if option == "--run-MapPER":
		run_mapper = 1
	    if option == "--DEBUG":
		DEBUG = 1
	    #if option == ("--pairend"):
                #pairend = "1"
	    #if option == ("--annotgene"):
		#do_annot_gene = 1;
		#annot_gene_file = value;
	    #if option == ("--filter-fusion-by-repeat"):
		#do_filter_fusion_by_repeat = 1;
		#chrom_blat = value;
	    if option in ("-1", "--end_1"):
                input_reads_1 = value
	    if option in ("-2", "--end_2"):
                input_reads_2 = value;
		
	    #if option == ("--split-index"):
		#split_index = value;
		
		#print >> sys.stderr, "chrom_blat=[%s]" % (chrom_blat)

	start_time = datetime.now()
	
        prepare_output_dir()
	
	#print_arguments("argu_log");
	print_argu = 1
	if print_argu > 0:
	    #print >> sys.stderr, "print argugment"
	
	    argu_file = logging_dir + "argu_log";
	    
	    argu_log = open(argu_file, "w")
	    
	    
	    print >> argu_log, "output_dir=[%s]" % (output_dir)
	    
	    print >> argu_log, "min_anchor_length=[%s]" % (min_anchor_length)
	    
	    print >> argu_log, "seed_length=[%s]" % (seed_length)
	
	    print >> argu_log, "splice_mismatches=[%s]" % (splice_mismatches)
	
	    print >> argu_log, "segment_mismatches=[%s]" % (segment_mismatches)
	
	    print >> argu_log, "FASTA_file_extension=[%s]" % (FASTA_file_extension)
	
	    print >> argu_log, "read_file_suffix=[%s]" % (read_file_suffix)
	
	    print >> argu_log, "min_intron_length=[%s]" % (min_intron_length)
	
	    print >> argu_log, "max_intron_length=[%s]" % (max_intron_length)
	
	    print >> argu_log, "island_extension=[%s]" % (island_extension)
	    
	    print >> argu_log, "read_width=[%s]" % (read_width)
	
	    print >> argu_log, "rank=[%s]" % (rank)
	
	    print >> argu_log, "flank_case=[%s]" % (flank_case)
	
	    print >> argu_log, "fusion_flank_case=[%s]" % (fusion_flank_case)
	
	    print >> argu_log, "islands_file=[%s]" % (islands_file)
	
	    print >> argu_log, "read_files_dir=[%s]" % (read_files_dir)
	
	    print >> argu_log, "chromosome_files_dir=[%s]" % (chromosome_files_dir)
	    
	    print >> argu_log, "all_chromosomes_file=[%s]" % (all_chromosomes_file)
	
	    print >> argu_log, "repeat_regioins=[%s]" % (repeat_regioins)
	    
	    print >> argu_log, "gene_regions=[%s]" % (gene_regions)
	    
	    print >> argu_log, "bwt_idx_prefix=[%s]" % (bwt_idx_prefix)
	    
	    print >> argu_log, "bowtie_threads=[%s]" % (bowtie_threads)
	    
	    print >> argu_log, "max_hits=[%s]" % (max_hits)
	    
	    print >> argu_log, "threshold=[%s]" % (threshold)
	    
	    print >> argu_log, "boundary=[%s]" % (boundary)
	
	    print >> argu_log, "num_anchor=[%s]" % (num_anchor)
	    
	    print >> argu_log, "unmapped_reads=[%s]" % (unmapped_reads)
	
	    print >> argu_log, "sam_formatted=[%s]" % (sam_formatted)
	
	    print >> argu_log, "bam_file=[%s]" % (bam_file)
	    
	    print >> argu_log, "sam_formatted_25=[%s]" % (sam_formatted_25)
	 
	    print >> argu_log, "bwt_map_25=[%s]" % (bwt_map_25)
	    
	    print >> argu_log, "pileup_file=[%s]" % (pileup_file)
	    
	    print >> argu_log, "synthetic_mappedreads=[%s]" % (synthetic_mappedreads)
	    
	    print >> argu_log, "tophat_mappedreads=[%s]" % (tophat_mappedreads)
	    
	    print >> argu_log, "pairend=[%s]" % (pairend)
	    
	    print >> argu_log, "gamma=[%s]" % (gamma)
	    
	    print >> argu_log, "delta=[%s]" % (delta)
	    
	    #print >> argu_log, "num_seg=[%s]" % (num_seg)
	
	    print >> argu_log, "seg_len=[%s]" % (seg_len)
	
	    print >> argu_log, "fix_hole_file=[%s]" % (fix_hole_file)
	    
	    print >> argu_log, "format_flag=[%s]" % (format_flag)
	    
	    print >> argu_log, "chrom_size_file=[%s]" % (chrom_size_file)
	    
	    print >> argu_log, "extend_bits=[%s]" % (extend_bits)
	    
	    print >> argu_log, "total_fusion_mismatch=[%s]" % (total_fusion_mismatch)
	    
	    print >> argu_log, "total_mismatch=[%s]" % (total_mismatch)
	    
	    print >> argu_log, "append_mismatch=[%s]" % (append_mismatch)
	    
	    print >> argu_log, "remap_mismatch=[%s]" % (remap_mismatch)
	    
	    print >> argu_log, "skip_bwt=[%s]" % (skip_bwt)
	    
	    print >> argu_log, "prefix_match=[%s]" % (prefix_match)
	    
	    print >> argu_log, "fullrunning=[%s]" % (fullrunning)
	    
	    print >> argu_log, "collect_stat=[%s]" % (collect_stat)
	    
	    print >> argu_log, "rm_temp=[%s]" % (rm_temp)
	    
	    print >> argu_log, "format_reads=[%s]" % (format_reads)
	    
	    print >> argu_log, "format_chromos=[%s]" % (format_chromos)
	    
	    print >> argu_log, "do_fusion=[%s]" % (do_fusion)
	    
	    print >> argu_log, "do_cluster=[%s]" % (do_cluster)
	    
	    print >> argu_log, "search_whole_chromo=[%s]" % (search_whole_chromo)
	    
	    print >> argu_log, "map_segment_directly=[%s]" % (map_segment_directly)
	    
	    print >> argu_log, "run_mapper=[%s]" % (run_mapper)
	    
	    print >> argu_log, "max_insert=[%s]" % (max_insert)
	    
	    print >> argu_log, "max_delete=[%s]" % (max_delete)
	    
	    print >> argu_log, "min_missed_seg=[%s]" % (min_missed_seg)
	    
	    print >> argu_log, "min_map_len=[%s]" % (min_map_len)
	     
	    print >> argu_log, "do_annot_gene=[%s]" % (do_annot_gene)
	    
	    print >> argu_log, "do_annot_gene=[%s]" % (annot_gene_file)
	    
	    print >> argu_log, "do_filter_fusion_by_repeat=[%s]" % (do_filter_fusion_by_repeat)
	    
	    print >> argu_log, "do_filter_fusion_by_repeat=[%s]" % (chrom_blat)
	    
	    print >> argu_log, "junction_tobe_syn=[%s]" % (junction_tobe_syn)
	    
	    print >> argu_log, "rerun_all=[%s]" % (rerun_all)
	    
	    print >> argu_log, "DEBUG=[%s]" % (DEBUG)
	    
	    print >> argu_log, "Split_index=[%s]" % (split_index)
	    
	    print >> argu_log, "input_reads_1=[%s]" % (input_reads_1)
	    
	    print >> argu_log, "input_reads_2=[%s]" % (input_reads_2)
	    
	    argu_log.close();
		
	#if read_width == 0:
	    #print >> sys.stderr, "read width must be specified"
	    #exit(0)
	    
	#if format_flag == "":
	    #print >> sys.stderr, "read format must be specified"
	    #exit(0)
	    
	#if bwt_idx_prefix == "":
	    #print >> sys.stderr, "bowtie index must be specified"
	    #exit(0)

	if bwt_idx_prefix == "":
	    print >> sys.stderr, "bowtie index must be specified"
	    exit(0)	    

        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning Mapsplice run (v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------" 

	global bin_dir
	
	print >> sys.stderr, "bin directory: [%s] " % bin_dir
	
	if input_reads_1 != "" and input_reads_2 != "":
	    pairend = "1";
	    
	is_paired = 0;
	
	if pairend != "":
	    is_paired = 1;
	
	#write_current_stats("Beginning MapSplice run");
	
	#no_original_junction = False;
        ## Validate all the input files, check all prereqs before committing 
        ## to the run
	
	check_file_existence(chromosome_files_dir);
	
	#check_file_existence(split_index);

	#splitted_reads1 = input_reads_1.split(',');
	
	#for read_file in splitted_reads1:
	    #check_file_existence(read_file);
	
	#splitted_reads2 = [];
	
	#if input_reads_2 != "":
	    #splitted_reads2 = input_reads_2.split(',');
	    
	    #for read_file in splitted_reads2:
		#check_file_existence(read_file);
	#####Read chromosome size
	
	##check_bowtie_index(bwt_idx_prefix, all_chromosomes_file, 0, FASTA_file_extension)
	
	#check_bowtie_index(bwt_idx_prefix, chromosome_files_dir, 0, FASTA_file_extension)
	
	##check_split_index(split_index, chromosome_files_dir, 0, FASTA_file_extension, 5000000, 4, 5);
	
	all_chromos_path = read_dir_by_suffix(chromosome_files_dir, FASTA_file_extension);
	
	chromo_fai_file = temp_dir + "chromo.fai"
	
	chromo_size_file = read_chromo_sizes(all_chromos_path, 
	                                     temp_dir + "chrom_sizes", 
	                                     fusion_dir + "chrName.txt", 
	                                     temp_dir + "chrom_head",
	                                     chromo_fai_file);
	
	is_fastq = 0;
	
	unmapped_reads_sams = [];
	
	####Filter alignments
	
	max_mate_dist = 50000;
	
	filtered_alignment_base = remap_dir + "_filtered_normal_alignments";
	
	filtered_alignment_append = "_filtered_normal_alignments";
	
	
	
	filter_flag = 12 + 32;
	
	min_insertion = max_insert;
	
	max_deletion = max_delete;
	
	min_anchor = 0;
	
	min_junction_anchor = 10;
	
	min_mismatch = 5;
	
	add_soft_clip = 1;
	
	mate_dist_sd = 100;
	
	max_anchor_diff = 50;
	
	intron_dist_sd = 500;
	
	chromosome_size_file = chromo_size_file;
	
	encompassing_fusion_region_extension = 50000;
	
	#input_sam_file = remapped_sam;
	
	min_coverage = 0;
	
	if is_paired > 0 and do_fusion > 0:
	    filter_flag = filter_flag + 256;
	
	fragment_length = 400;

	fragment_length_sd = 100;

	avearge_fragment_length = 225;

	boundary = 36;

	min_isoform_length = read_width / 2;

	min_encompass_count = 1;
	
	    
	####Filter alignments with fusion
	
	bam2sam(dec_sam_file, fusion_dir + "old_alignment.sam");
	
	input_sam_file = fusion_dir + "old_alignment.sam";

	read_files_dir = sam2fq(input_sam_file, temp_dir + "sam_reads.fq", "0", logging_dir + "sam2fq.log")
		    
	format_flag = '-q';

	check_reads_format(read_files_dir, min_read_len, pairend);
	
	(maxlen, format_flag, qual_scale_detected) = extract_maxlen(logging_dir + "check_reads_format.log");
	
	read_width = maxlen;
	
	max_read_length = read_width;
	
	is_paired = 1;	

	if issortedbam > 0:
	    SepSamMappedUnmapped(fusion_dir + "old_alignment.sam", 
		                 fusion_dir + "dec_7_sam.head",
		                 fusion_dir + "dec_7_sam.mapped.normal",
		                 fusion_dir + "dec_7_sam.mapped.fusion",
		                 fusion_dir + "dec_7_sam.unmapped",
		                 logging_dir + "SepSamMappedUnmapped.log")
	    
	    sort_by_name_c([fusion_dir + "dec_7_sam.mapped.fusion"], fusion_dir + "dec_7_sam.mapped.fusion.sorted");
	    
	    recover_fusion_alignments_order(fusion_dir + "dec_7_sam.mapped.fusion.sorted", 
		                            fusion_dir + "dec_7_sam.mapped.fusion.converted",
		                            logging_dir + "recover_fusion_alignments_order.log");
	
	    files_tobe_sort = [fusion_dir + "dec_7_sam.mapped.fusion.converted"] + [fusion_dir + "dec_7_sam.mapped.normal"];
	    
	    sort_by_name_c(files_tobe_sort, fusion_dir + "converted_alignments.sorted.txt");
	    
	    input_sam_file = fusion_dir + "converted_alignments.sorted.txt";
	    
	else:
	    SepSamMappedUnmapped(fusion_dir + "old_alignment.sam", 
	                         fusion_dir + "dec_7_sam.head",
	                         fusion_dir + "dec_7_sam.mapped.normal.fusion",
	                         fusion_dir + "dec_7_sam.mapped.normal.fusion",
	                         fusion_dir + "dec_7_sam.unmapped",
	                         logging_dir + "SepSamMappedUnmapped.log")	
	    
	    input_sam_file = fusion_dir + "dec_7_sam.mapped.normal.fusion";
  
	filtered_fusion_alignment_base = fusion_dir + "_filtered_fusion_alignments";
	    
	filtered_fusion_alignment_append = "_filtered_fusion_alignments";
	
	#
	
	filter_flag = 12 + 32 + 128 + 1024;
	
	fragment_length = 400;

	fragment_length_sd = 100;

	avearge_fragment_length = 225;

	boundary = 36;

	min_isoform_length = read_width / 2;

	min_encompass_count = 1;
	
	if maxlen > 75:
	    min_encompass_count = 0;
		
	run_alignment_handler_multi("", 
                                    is_paired,
                                    max_mate_dist,
                                    max_hits * 10,
                                    filtered_fusion_alignment_base,
                                    filtered_fusion_alignment_append,
                                    max_read_length,
                                    chromosome_files_dir,
                                    filter_flag,
                                    min_insertion,
                                    max_deletion,
                                    min_anchor,
                                    min_junction_anchor,
                                    min_mismatch,
                                    add_soft_clip,
                                    mate_dist_sd,
                                    intron_dist_sd,
                                    max_anchor_diff,
                                    chromosome_size_file,
                                    encompassing_fusion_region_extension,
                                    bowtie_threads,
                                    min_coverage,
                                    fragment_length,
                                    fragment_length_sd,
                                    avearge_fragment_length,
                                    boundary,
                                    min_isoform_length,
                                    min_encompass_count,
                                    min_entropy,
                                    input_sam_file,
                                    logging_dir + "alignment_handler_fusion.log",
                                    logging_dir + "alignment_handler_fusion.err")
	    
    
	fusion_alignments = input_sam_file + "_filtered_fusion_alignments"
	
	unmapped_reads_sams = [fusion_dir + "dec_7_sam.unmapped"] + [fusion_dir + "_filtered_fusion_alignments.unmapped"]
	
	final_alignments = [fusion_alignments];
	
	final_alignments = final_alignments + [input_sam_file + "_filtered_fusion_alignments.bothunspliced.paired"]
	
	#if do_fusion > 0 and is_paired > 0:
	    #if no_original_fusion == False:
		#final_alignments = [fusion_alignments] + paired_alignments;
	    #else:
		#final_alignments = paired_alignments + [single_alignments] + [fusion_paired_alignments];
	#elif is_paired > 0 and do_fusion == 0:
	    #final_alignments = paired_alignments + [single_alignments] + [fusion_paired_alignments];
	#elif is_paired == 0 and do_fusion == 0:
	    #final_alignments = [single_alignments];
	#elif is_paired == 0 and do_fusion > 0:
	    #if no_original_fusion == False:
		#final_alignments = [fusion_alignments] + [single_alignments];
	    #else:
		#final_alignments = [single_alignments];
	
	final_alignments = final_alignments + unmapped_reads_sams;
	
	#final_alignments_sorted = sort_by_name_c(final_alignments, temp_dir + "final_alignments.sam");
	
	add_head = [temp_dir + "chrom_head"] + final_alignments;
	
	if bam_file != "":
	    cat_files(add_head, temp_dir + "final_alignments_headed.sam");
	    
	    sam2bam(temp_dir + "final_alignments_headed.sam", output_dir + "alignments.bam", chromo_fai_file);
	else:
	    cat_files(add_head, output_dir + "alignments.sam");
        
	#copy_file(filtered_alignment_base + ".fil.junc", output_dir + "junctions.txt");
	
	#copy_file(filtered_alignment_base + ".fil.junc.del", output_dir + "deletions.txt");
	
	#copy_file(filtered_alignment_base + ".fil.junc.ins", output_dir + "insertions.txt");
	
	#copy_file(filtered_alignment_base + ".stat", output_dir + "stats.txt");
	
	if do_fusion > 0:
	    copy_file(filtered_fusion_alignment_base + ".stat", output_dir + "fusion_stat.txt");
	    
	    if is_paired > 0:
		copy_file(filtered_fusion_alignment_base + ".ori.junc.fusion.encompass.no.compass", output_dir + "fusions_raw.txt");
	    else:
		copy_file(filtered_fusion_alignment_base + ".ori.junc.fusion.encompass.no.compass", output_dir + "fusions_raw.txt");

	    #if gene_gtf_file != "":
		#gtf2exons(gene_gtf_file, temp_dir + "gene_exons.txt", logging_dir + "gtf_to_exons.log");
	    
	    min_anchor_repeat = 15;
	    
	    if maxlen / 4 < min_anchor_repeat:
		min_anchor_repeat = 10;
		
	    max_anchor_diff = maxlen / 2;
	    
	    filterbyrepeats(output_dir + "fusions_raw.txt",
	                    output_dir + "junctions.txt",
	                    output_dir + "fusions_candidates.txt",
	                    chromosome_files_dir,
	                    seg_len,
	                    100,
	                    max_anchor_diff,
	                    min_anchor_repeat,
	                    min_entropy_repeats,
	                    bowtie_threads,
	                    bwt_idx_prefix,
	                    chromosome_size_file,
	                    logging_dir)
	    
	    if gene_gtf_file != "":
		filterbyannotation(output_dir + "fusions_candidates.txt",
		                   output_dir + "fusions_well_annotated.txt",
		                   output_dir + "fusions_not_well_annotated.txt",
		                   output_dir + "fusions_combined_annotated_long_intron_and_interchr.txt",
		                   output_dir + "circular_RNAs.txt",
		                   gene_gtf_file,
		                   logging_dir)
        
	if rm_temp and os.path.exists(temp_dir):
	    shutil.rmtree(temp_dir, True)
	    
        finish_time = datetime.now()
        duration = finish_time - start_time
	
	print >> sys.stderr
        print >> sys.stderr, "[%s] Finishing Mapsplice run (time used: %s)" % (right_now(), duration)
        print >> sys.stderr, "-----------------------------------------------" 
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        return 2


if __name__ == "__main__":
    sys.exit(main())
