all: collectstats SetUnmappedBitFlag sam2fq recover_fusion_alignments_order SepSamUnmapped swap_dRanger_and_MPS_matched comp_fusiondb_offset sepdRangerfusion sepMPSfusion generate_fusiongene_convert_coordinate_trim generate_fusiongene_convert_coordinate_trim_dRanger junc2bed SeparateNormalFromFusionJunc search_fusion_gene AddFusionStrandConsistent gtf2genetab FilterFusionByNormalPaired load_fusion_chrom_seq_std matchfusion2normal filteroriginalfusion Convert2FusionAlignment generate_combined_sequence DNA2StdRegion search_unmapped_reads reads2unmappedsam junc_db_fusion bsb4 comb_junc mapsplice_multi_thread check_reads_format samtools alignmenthandler_multi parseCluster cluster junc_db filterjuncbyROCarguNonCanonical filter_1hits newsam2junc fusionsam2junc_filteranchor_newfmt read_chromo_size find_mate_sam_fq alignmenthandler SepSam

OPTFLAGS = -O3

#LDFLAGS = -static

CFLAGS += $(OPTFLAGS)

#CFLAGS += $(LDFLAGS)

collectstats:
	g++ $(CFLAGS) -o bin/collectstats src/collectstats/collectstats.cpp 2>log

SetUnmappedBitFlag:
	g++ $(CFLAGS) -o bin/SetUnmappedBitFlag src/SetUnmappedBitFlag/SetUnmappedBitFlag.cpp 2>log

sam2fq:
	g++ $(CFLAGS) -o bin/sam2fq src/sam2fq/sam2fq.cpp 2>log

recover_fusion_alignments_order:
	g++ $(CFLAGS) -o bin/recover_fusion_alignments_order src/recover_fusion_alignments_order/recover_fusion_alignments_order.cpp src/recover_fusion_alignments_order/JunctionSeed.cpp src/recover_fusion_alignments_order/SamRec.cpp src/recover_fusion_alignments_order/sharedlib.cpp src/recover_fusion_alignments_order/SpliceWay.cpp 2>log

SepSamUnmapped:
	g++ $(CFLAGS) -o bin/SepSamUnmapped src/SepSamUnmapped/SepSamUnmapped.cpp 2>log

swap_dRanger_and_MPS_matched:
	g++ $(CFLAGS) -o bin/swap_dRanger_and_MPS_matched src/swap_dRanger_and_MPS_matched/swap_dRanger_and_MPS_matched.cpp 2>log

comp_fusiondb_offset:
	g++ $(CFLAGS) -o bin/comp_fusiondb_offset src/comp_fusiondb_offset/comp_fusiondb_offset.cpp 2>log

sepdRangerfusion:
	g++ $(CFLAGS) -o bin/sepdRangerfusion src/sepdRangerfusion/sepdRangerfusion.cpp 2>log

sepMPSfusion:
	g++ $(CFLAGS) -o bin/sepMPSfusion src/sepMPSfusion/sepMPSfusion.cpp 2>log

generate_fusiongene_convert_coordinate_trim:

	g++ -O3 -o bin/generate_fusiongene_convert_coordinate_trim src/generate_fusiongene_convert_coordinate/generate_fusiongene_convert_coordinate_trim.cpp src/generate_fusiongene_convert_coordinate/JunctionSeed.cpp src/generate_fusiongene_convert_coordinate/SamRec.cpp src/generate_fusiongene_convert_coordinate/sharedlib.cpp src/generate_fusiongene_convert_coordinate/SpliceWay.cpp 2>log

generate_fusiongene_convert_coordinate_trim_dRanger:

	g++ -O3 -o bin/generate_fusiongene_convert_coordinate_trim_dRanger src/generate_fusiongene_convert_coordinate/generate_fusiongene_convert_coordinate_trim_dRanger.cpp src/generate_fusiongene_convert_coordinate/JunctionSeed.cpp src/generate_fusiongene_convert_coordinate/SamRec.cpp src/generate_fusiongene_convert_coordinate/sharedlib.cpp src/generate_fusiongene_convert_coordinate/SpliceWay.cpp 2>log

junc2bed:
	g++ $(CFLAGS) -o bin/junc2bed src/junc2bed/junc2bed.cpp 2>log

SeparateNormalFromFusionJunc:
	g++ $(CFLAGS) -o bin/SeparateNormalFromFusionJunc src/SeparateNormalFromFusionJunc/SeparateNormalFromFusionJunc.cpp 2>log

search_fusion_gene:
	g++ $(CFLAGS) -o bin/search_fusion_gene src/search_fusion_gene/search_fusion_gene.cpp 2>log

AddFusionStrandConsistent:
	g++ $(CFLAGS) -o bin/AddFusionStrandConsistent src/AddFusionStrandConsistent/AddFusionStrandConsistent.cpp 2>log

gtf2genetab:
	g++ $(CFLAGS) -o bin/gtf2genetab src/gtf2genetab/gtf2genetab.cpp 2>log

bsb4:
	g++ $(CFLAGS) -o bin/bsb4 src/bsb4/bsb4.cpp 2>log

samtools:
	cd ./samtools-0.1.9;make

	cp ./samtools-0.1.9/samtools ./bin/

mapsplice_multi_thread:
	cd ./src/MapSplice;make

	cp ./src/MapSplice/bowtie ./bin/mapsplice_multi_thread

reads2unmappedsam:
	g++ $(CFLAGS) -o bin/reads2unmappedsam src/reads2unmappedsam/reads2unmappedsam.cpp 2>log

DNA2StdRegion:
	g++ $(CFLAGS) -o bin/DNA2StdRegion src/DNA2StdRegion/DNA2StdRegion.cpp 2>log

generate_combined_sequence:
	g++ $(CFLAGS) -o bin/generate_combined_sequence src/generate_combined_sequence/generate_combined_sequence.cpp 2>log

Convert2FusionAlignment:
	g++ $(CFLAGS) -o bin/Convert2FusionAlignment src/Convert2FusionAlignment/Convert2FusionAlignment.cpp src/Convert2FusionAlignment/sharedlib.cpp src/Convert2FusionAlignment/SamRec.cpp src/Convert2FusionAlignment/JunctionSeed.cpp 2>log

check_reads_format:
	g++ $(CFLAGS) -o bin/check_reads_format src/check_reads_format/check_reads_format.cpp 2>log

search_unmapped_reads:
	g++ $(CFLAGS) -o bin/search_unmapped_reads src/search_unmapped_reads/search_unmapped_reads.cpp 2>log

parseCluster:
	g++ $(CFLAGS) -o bin/parseCluster src/cluster/preparseFusion.cpp 2>log
	
cluster:
	g++ $(CFLAGS) -o bin/cluster src/cluster/PE_match_junction.cpp 2>log
	
junc_db:
	g++ $(CFLAGS) -o bin/junc_db src/alignmenthandler/junc_db.cpp src/alignmenthandler/JunctionHandler.cpp src/alignmenthandler/SamRec.cpp src/alignmenthandler/sharedlib.cpp src/alignmenthandler/JunctionSeed.cpp 2>log

junc_db_fusion:
	g++ $(CFLAGS) -o bin/junc_db_fusion src/junction_database_fusion/junction_seq_construction_fusion.cpp 2>log

filterjuncbyROCarguNonCanonical:
	g++ $(CFLAGS) -o bin/filterjuncbyROCarguNonCanonical src/filterjuncbyROCarguNonCanonical/filterjuncbyROCarguNonCanonical.cpp 2>log

filter_1hits:
	g++ $(CFLAGS) -o bin/filter_1hits src/filter_1hits/filter_1hits.cpp 2>log
	
newsam2junc:
	g++ $(CFLAGS) -o bin/newsam2junc src/newsam2junc/newsam2junc.cpp 2>log

fusionsam2junc_filteranchor_newfmt:
	g++ $(CFLAGS) -o bin/fusionsam2junc_filteranchor_newfmt src/fusionsam2junc_filteranchor_newfmt/fusionsam2junc_filteranchor_newfmt.cpp 2>log

read_chromo_size:
	g++ $(CFLAGS) -o bin/read_chromo_size src/read_chromo_size/read_chromo_size.cpp 2>log

find_mate_sam_fq:
	g++ $(CFLAGS) -o bin/find_mate_sam_fq src/find_mate_sam_fq/find_mate_sam_fq.cpp 2>log

alignmenthandler:
	g++ $(CFLAGS) -o bin/alignment_handler src/alignmenthandler/test_handler.cpp src/alignmenthandler/JunctionHandler.cpp src/alignmenthandler/SamRec.cpp src/alignmenthandler/sharedlib.cpp src/alignmenthandler/AlignmentHandler.cpp src/alignmenthandler/JunctionSeed.cpp src/alignmenthandler/FusionSamRec.cpp src/alignmenthandler/UnionExpressedRegions.cpp src/alignmenthandler/disjointset.cpp 2>log

alignmenthandler_multi:
	g++ $(CFLAGS) -o bin/alignment_handler_multi src/alignmenthandler_multi/test_handler.cpp src/alignmenthandler_multi/JunctionHandler.cpp src/alignmenthandler_multi/SamRec.cpp src/alignmenthandler_multi/sharedlib.cpp src/alignmenthandler_multi/AlignmentHandler.cpp src/alignmenthandler_multi/JunctionSeed.cpp src/alignmenthandler_multi/SpliceWay.cpp src/alignmenthandler_multi/FusionSamRec.cpp src/alignmenthandler_multi/UnionExpressedRegions.cpp src/alignmenthandler_multi/disjointset.cpp -lpthread 2>log

comb_junc:
	g++ $(CFLAGS) -o bin/comb_junc src/alignmenthandler_multi/comb_junc.cpp src/alignmenthandler_multi/JunctionHandler.cpp src/alignmenthandler_multi/SpliceWay.cpp src/alignmenthandler_multi/SamRec.cpp src/alignmenthandler_multi/sharedlib.cpp src/alignmenthandler_multi/JunctionSeed.cpp 2>log

SepSam:
	g++ $(CFLAGS) -o bin/SepSam src/SepSam/SepSam.cpp 2>log

filteroriginalfusion:
	g++ $(CFLAGS) -o bin/filteroriginalfusion src/filteroriginalfusion/filteroriginalfusion.cpp 2>log

matchfusion2normal:
	g++ $(CFLAGS) -o bin/matchfusion2normal src/matchfusion2normal/matchfusion2normal.cpp 2>log

load_fusion_chrom_seq_std:
	g++ $(CFLAGS) -o bin/load_fusion_chrom_seq_std src/load_fusion_chrom_seq_std/load_fusion_chrom_seq_std.cpp 2>log

FilterFusionByNormalPaired:
	g++ $(CFLAGS) -o bin/FilterFusionByNormalPaired src/FilterFusionByNormalPaired/FilterFusionByNormalPaired.cpp 2>log

clean:
	rm bin/comp_fusiondb_offset bin/SeparateNormalFromFusionJunc bin/collectstats bin/SetUnmappedBitFlag bin/sam2fq bin/recover_fusion_alignments_order bin/SepSamUnmapped bin/swap_dRanger_and_MPS_matched bin/sepdRangerfusion bin/sepMPSfusion bin/generate_fusiongene_convert_coordinate_trim_dRanger bin/generate_fusiongene_convert_coordinate_trim bin/junc2bed bin/search_fusion_gene bin/AddFusionStrandConsistent bin/gtf2genetab bin/FilterFusionByNormalPaired bin/load_fusion_chrom_seq_std bin/matchfusion2normal bin/filteroriginalfusion bin/DNA2StdRegion bin/generate_combined_sequence bin/Convert2FusionAlignment bin/search_unmapped_reads bin/reads2unmappedsam bin/junc_db_fusion bin/bsb4 bin/fusionsam2junc_filteranchor_newfmt bin/comb_junc bin/mapsplice_multi_thread bin/check_reads_format bin/samtools bin/alignment_handler_multi bin/parseCluster bin/cluster bin/junc_db bin/filterjuncbyROCarguNonCanonical bin/filter_1hits bin/newsam2junc bin/read_chromo_size bin/find_mate_sam_fq bin/alignment_handler bin/SepSam