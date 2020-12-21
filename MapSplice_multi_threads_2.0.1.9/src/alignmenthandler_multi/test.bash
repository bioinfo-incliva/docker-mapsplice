binpath=/bioinfo/projects/kai/
chrompath=/bioinfo/projects/zeng/mapsplice_multi/hg19/
bin=/bioinfo/projects/kai/code/comp_pack/
MPS_SAM1_BASENAME=MPS_ori.SAM
MPS_SAM1_OUTPUT_DIR=./alignment_comparison_result/

#LC_ALL=C sort -k1,1 -t~ -S 10000000 -T ./ -o remap_alignment.sam /bioinfo/projects/zeng/rgasp/R1/R1_alignment_annotation_only.sam

/bioinfo/projects/kai/code/AlignmentHandlerMulti_8_4/test_handler "" 1 50000 40 _filtered_normal_alignments 76 /bioinfo/projects/zeng/mapsplice_multi/hg19/ 12 3 10 0 10 5 1 20 50 /bioinfo/projects/zeng/mapsplice_multi/human_chrom_sizes 50000 12 500 remap_alignment.sam >log.txt

