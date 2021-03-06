#include "AlignmentHandler.h"

bool comp_sam_rec_maplen(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.mappedlen > rhs.mappedlen;
}

bool comp_min_anchor(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.min_anchor > rhs.min_anchor;
}

bool comp_mate_dist(size_t lhs, size_t rhs)
{
	//cout << mate_dist_sd << endl;
	if (lhs >= rhs && lhs - rhs <= mate_dist_sd)
		return false;

	if (rhs > lhs  && rhs - lhs <= mate_dist_sd)
		return false;

	return lhs < rhs;
}

bool comp_intron_dist(size_t lhs, size_t rhs)
{
	//cout << mate_dist_sd << endl;
	if (lhs >= rhs && lhs - rhs <= intron_dist_sd)
		return false;

	if (rhs > lhs  && rhs - lhs <= intron_dist_sd)
		return false;

	return lhs < rhs;
}

bool better_than(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (comp_mate_dist(lhs.mate_dist, rhs.mate_dist) && comp_mate_dist(rhs.mate_dist, lhs.mate_dist))
	{
		if (lhs.total_mismatch == rhs.total_mismatch)
		{
			if (lhs.mappedlen == rhs.mappedlen)
			{
				if (lhs.intron_size == rhs.intron_size)
					return false;
				else 
					return lhs.intron_size < rhs.intron_size;
			}
			else
				return lhs.mappedlen > rhs.mappedlen;
		}			
		else
		{
			return lhs.total_mismatch < rhs.total_mismatch;
		}
	}
	else
		return lhs.mate_dist < rhs.mate_dist;
}

bool comp_dist(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.mappedlen == rhs.mappedlen)
	{
		if (lhs.total_mismatch == rhs.total_mismatch)
		{
			if (!comp_mate_dist(lhs.mate_dist, rhs.mate_dist) && !comp_mate_dist(rhs.mate_dist, lhs.mate_dist))
			{
				//if (lhs.total_ave_mismatch == rhs.total_ave_mismatch)
				//{
					if (!comp_intron_dist(lhs.intron_size, rhs.intron_size) && !comp_intron_dist(rhs.intron_size, lhs.intron_size))
					{
						if (lhs.total_hits == rhs.total_hits)
						{
							if (lhs.total_filter_score == rhs.total_filter_score)
							{
								if (lhs.contiglen == rhs.contiglen)
									return lhs.min_anchor > rhs.min_anchor;
								else;
								return lhs.contiglen > rhs.contiglen;
							}
							else
								return lhs.total_filter_score > rhs.total_filter_score;
						}
						else
							return lhs.total_hits > rhs.total_hits;
					}
					else
						return lhs.intron_size < rhs.intron_size;
				//}
				//else
				//	return lhs.total_ave_mismatch < rhs.total_ave_mismatch;

			}
			else
				return lhs.mate_dist < rhs.mate_dist;
		}			
		else
		{
			return lhs.total_mismatch < rhs.total_mismatch;
		}
	}
	else
		return lhs.mappedlen > rhs.mappedlen;
}

bool comp_dist_0(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.mappedlen == rhs.mappedlen)
	{
		if (lhs.total_mismatch == rhs.total_mismatch)
		{
			if (lhs.mate_dist == rhs.mate_dist)
			{
				//if (lhs.total_ave_mismatch == rhs.total_ave_mismatch)
				//{
					if (lhs.intron_size == rhs.intron_size)
					{
						if (lhs.total_hits == rhs.total_hits)
						{
							if (lhs.total_filter_score == rhs.total_filter_score)
							{
								if (lhs.contiglen == rhs.contiglen)
								{
									if (lhs.min_anchor == rhs.min_anchor)
										return lhs.total_confid > rhs.total_confid;
									else
										return lhs.min_anchor > rhs.min_anchor;
								}
								else
									return lhs.contiglen > rhs.contiglen;
							}
							else
								return lhs.total_filter_score > rhs.total_filter_score;
						}
						else
							return lhs.total_hits > rhs.total_hits;

						//if (lhs.total_filter_score == rhs.total_filter_score)
						//{
						//	if (lhs.contiglen == rhs.contiglen)
						//		return lhs.min_anchor > rhs.min_anchor;
						//	else;
						//	return lhs.contiglen > rhs.contiglen;
						//}
						//else
						//	return lhs.total_filter_score > rhs.total_filter_score;
					}
					else
						return lhs.intron_size < rhs.intron_size;
				//}
				//else
				//	return lhs.total_ave_mismatch < rhs.total_ave_mismatch;

			}
			else
				return lhs.mate_dist < rhs.mate_dist;
		}			
		else
		{
			return lhs.total_mismatch < rhs.total_mismatch;
		}
	}
	else
		return lhs.mappedlen > rhs.mappedlen;
}

bool comp_dist_1(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.intron_size == rhs.intron_size)
	{
		if (lhs.total_hits == rhs.total_hits)
		{
			if (lhs.total_filter_score == rhs.total_filter_score)
			{
				if (lhs.contiglen == rhs.contiglen)
				{
					if (lhs.min_anchor == rhs.min_anchor)
						return lhs.total_confid > rhs.total_confid;
					else
						return lhs.min_anchor > rhs.min_anchor;
				}
				else
					return lhs.contiglen > rhs.contiglen;
			}
			else
				return lhs.total_filter_score > rhs.total_filter_score;
		}
		else
			return lhs.total_hits > rhs.total_hits;
		//if (lhs.total_filter_score == rhs.total_filter_score)
		//{
		//	if (lhs.contiglen == rhs.contiglen)
		//		return lhs.min_anchor > rhs.min_anchor;
		//	else;
		//	return lhs.contiglen > rhs.contiglen;
		//}
		//else
		//	return lhs.total_filter_score > rhs.total_filter_score;
	}
	else
		return lhs.intron_size < rhs.intron_size;
}

bool comp_dist_2(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.total_hits == rhs.total_hits)
	{
		if (lhs.total_filter_score == rhs.total_filter_score)
		{
			if (lhs.contiglen == rhs.contiglen)
			{
				if (lhs.min_anchor == rhs.min_anchor)
					return lhs.total_confid > rhs.total_confid;
				else
					return lhs.min_anchor > rhs.min_anchor;
			}
			else
				return lhs.contiglen > rhs.contiglen;
		}
		else
			return lhs.total_filter_score > rhs.total_filter_score;
	}
	else
		return lhs.total_hits > rhs.total_hits;
	/*if (lhs.total_filter_score == rhs.total_filter_score)
	{
		if (lhs.contiglen == rhs.contiglen)
			return lhs.min_anchor > rhs.min_anchor;
		else;
		return lhs.contiglen > rhs.contiglen;
	}
	else
		return lhs.total_filter_score > rhs.total_filter_score;*/
}

bool comp_mis(const SamRec* lhs, const SamRec* rhs)
{
	return lhs->mis_match < rhs->mis_match;
}

bool comp_dup(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->chrom_name == rhs->chrom_name)
	{
		if (lhs->start == rhs->start)
		{
			if (lhs->splice_way == rhs->splice_way)
			{
				if (lhs->strand_t == rhs->strand_t)
				{
					if (lhs->chrom_name2 == rhs->chrom_name2)
					{
						if (lhs->start2 == rhs->start2)
						{
							if (lhs->splice_way2 == rhs->splice_way2)
								return lhs->strand_t2  < rhs->strand_t2;
							else
								return lhs->splice_way2 < rhs->splice_way2;
						}
						else
							return lhs->start2 < rhs->start2;
					}
					else
						return lhs->chrom_name2 < rhs->chrom_name2;
				}
				else
					return lhs->strand_t  < rhs->strand_t;
			}
			else
				return lhs->splice_way < rhs->splice_way;
		}
		else
			return lhs->start < rhs->start;
	}
	else
		return lhs->chrom_name < rhs->chrom_name;
}

bool
comp_filterscore(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->isexonic && !rhs->isexonic)
		return true;
	else if (!lhs->isexonic && rhs->isexonic)
		return false;

	//if (lhs->mappedlen == rhs->mappedlen)
	//{
	if (lhs->filter_score == rhs->filter_score)
	{
		if (lhs->ave_intron_len == rhs->ave_intron_len)
		{
			if (lhs->mappedlen == rhs->mappedlen)
			{
				if (lhs->ave_junc_mis == rhs->ave_junc_mis)
					return lhs->m_contig_len > rhs->m_contig_len;
				else
					return lhs->ave_junc_mis < rhs->ave_junc_mis;
			}
			else
				return lhs->mappedlen > rhs->mappedlen;
		}
		else
		{
			return lhs->ave_intron_len < rhs->ave_intron_len;
		}
	}
	else
		return lhs->filter_score > rhs->filter_score;
	//}
	//else
	//	return lhs->mappedlen > rhs->mappedlen;
}

bool comp_fusion_filterscore(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->isexonic && !lhs->is_fusion && !rhs->isexonic)
		return true;
	else if (!lhs->isexonic && rhs->isexonic && !rhs->is_fusion)
		return false;

	if (lhs->mappedlen == rhs->mappedlen)
	{
		if (lhs->filter_score == rhs->filter_score)
		{
			if (lhs->ave_junc_mis == rhs->ave_junc_mis)
				return lhs->m_contig_len > rhs->m_contig_len;
			else
				return lhs->ave_junc_mis < rhs->ave_junc_mis;
		}
		else
			return lhs->filter_score > rhs->filter_score;
	}
	else
		return lhs->mappedlen > rhs->mappedlen;
	//}
	//else
	//	return lhs->mappedlen > rhs->mappedlen;
}

bool check_overlap(vector<pair<size_t, int> >& spliceway_vec_1, vector<pair<size_t, int> >& spliceway_vec_2)
{
	vector<pair<size_t, int> >::iterator sp_vec_iter1, sp_vec_iter2;

	for (sp_vec_iter1 = spliceway_vec_1.begin(); sp_vec_iter1 != spliceway_vec_1.end()/* - 1*/; ++sp_vec_iter1)
	{
		//if (sp_vec_iter1->second > 0/* && (sp_vec_iter1 + 1)->second > 0*/)
		{
			for (sp_vec_iter2 = spliceway_vec_2.begin(); sp_vec_iter2 != spliceway_vec_2.end() - 1; ++sp_vec_iter2)
			{
				if (sp_vec_iter2->second > 0 && (sp_vec_iter2 + 1)->second > 0)
				{
					if (sp_vec_iter1->first > sp_vec_iter2->first + sp_vec_iter2->second && 
						sp_vec_iter1->first < (sp_vec_iter2 + 1)->first)
						return true;


					//if (sp_vec_iter1->second > 0/* && (sp_vec_iter1 + 1)->second > 0*/)
					if (sp_vec_iter1->second > 0 &&
						sp_vec_iter1->first + sp_vec_iter1->second > sp_vec_iter2->first + sp_vec_iter2->second && 
						sp_vec_iter1->first + sp_vec_iter1->second - 1 < (sp_vec_iter2 + 1)->first)
						return true;
				}
			}
		}
	}

	return false;
}

struct EstPairingStruct {
	pair<vector<SamRec*>, vector<SamRec*> >* sam_rec_pe_ptr;

	vector<PairedSamRec>* fusion_paired_reads_ptr;

	int do_filter;

	EstPairingStruct(pair<vector<SamRec*>, vector<SamRec*> >* _sam_rec_pe_ptr, vector<PairedSamRec>* _fusion_paired_reads_ptr, int _do_filter)
		: sam_rec_pe_ptr(_sam_rec_pe_ptr), fusion_paired_reads_ptr(_fusion_paired_reads_ptr), do_filter(_do_filter) {}

};

struct FilterPairingStruct {

	vector<PairedSamRec>* pairedsamrec_vec;

	pair<vector<SamRec*>, vector<SamRec*> >* sam_rec_pe_ptr;

	int do_filter;

	FilterPairingStruct(vector<PairedSamRec>* _pairedsamrec_vec, pair<vector<SamRec*>, vector<SamRec*> >* _sam_rec_pe_ptr, int _do_filter)
		: pairedsamrec_vec(_pairedsamrec_vec), sam_rec_pe_ptr(_sam_rec_pe_ptr), do_filter(_do_filter) {}
};

struct FindFusionJuncRegionStruct {
	FusionJuncRegion* cur_fusion_region;

	PairedSamRec* cur_paired_sam_rec;

	vector<FusionJuncRegion>* sorted_fusion_regions;

	FindFusionJuncRegionStruct(FusionJuncRegion* _cur_fusion_region, PairedSamRec* _cur_paired_sam_rec, vector<FusionJuncRegion>* _sorted_fusion_regions)
		: cur_fusion_region(_cur_fusion_region), cur_paired_sam_rec(_cur_paired_sam_rec), sorted_fusion_regions(_sorted_fusion_regions)
	{
	}
};
size_t AlignmentHandler::m_cur_read_alignments_index = 0;// =  0;

size_t AlignmentHandler::m_cur_line_index = 0;

vector<string> AlignmentHandler::m_stored_lines(buf_size);

vector<SamRec> AlignmentHandler::m_stored_sam_recs(buf_size);

vector<pair<vector<SamRec>::iterator, vector<SamRec>::iterator> > AlignmentHandler::m_same_reads_alignments;

vector<pair<vector<SamRec*>, vector<SamRec*> > > AlignmentHandler::m_same_reads_alignments_paired;

vector<bool> AlignmentHandler::m_same_reads_alignments_paired_is_unspliced;

int AlignmentHandler::m_min_ins_static;

int AlignmentHandler::m_do_filter;

size_t AlignmentHandler::m_max_pair_dist;

JunctionHandler AlignmentHandler::m_junction_handler;

JunctionHandler AlignmentHandler::m_junction_handler_filtered;

JunctionHandler AlignmentHandler::m_junction_handler_intermediate;

JunctionHandler AlignmentHandler::m_junction_handler_annotated;

JunctionHandler* AlignmentHandler::m_junction_handler_prev_ptr;

double AlignmentHandler::m_entrpy_weight;

double AlignmentHandler::m_pqlen_weight;

double AlignmentHandler::m_ave_mis_weight;

bool AlignmentHandler::m_add_S;

size_t AlignmentHandler::m_max_hits;

size_t AlignmentHandler::m_spliced;
size_t AlignmentHandler::m_unspliced;
size_t AlignmentHandler::m_insertion;
size_t AlignmentHandler::m_deletion;
size_t AlignmentHandler::m_unique;
size_t AlignmentHandler::m_multiple;
size_t AlignmentHandler::m_canoical;
size_t AlignmentHandler::m_semi_canonical;
size_t AlignmentHandler::m_non_canonical;
size_t AlignmentHandler::m_paired;
size_t AlignmentHandler::m_single;
size_t AlignmentHandler::m_unmapped;
size_t AlignmentHandler::m_clipped;
size_t AlignmentHandler::m_fusion;

string AlignmentHandler::m_prev_tag;

bool AlignmentHandler::m_is_paired;

AlignmentHandler::AlignmentHandler(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, string filtered_alignment_file,
	size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor, char* chr_sz_file,
	size_t min_fusion_coverage, double entrpy_weight, double pqlen_weight, double ave_mis_weight) : /* m_sam_file_handler(alignment_file.c_str(), min_ins), */
	m_junction_file(junction_file), m_sam_files(alignment_files), 
	 m_filtered_alignment_file(filtered_alignment_file),
	/*m_ofs_filtered_alignment(m_filtered_alignment_file.c_str()),*/ m_max_read_width(max_read_width), m_chrom_dir(chrom_dir), 
	m_min_ins(min_ins), m_max_del(max_del),
	m_min_anchor(min_anchor), m_min_mismatch(min_mismatch), 
	m_min_junc_anchor(min_junc_anchor), m_chrom_size_file(chr_sz_file), m_min_fusion_coverage(min_fusion_coverage)
	//m_cur_read_alignments_index(0), m_cur_line_index(0)//, m_stored_lines(buf_size), m_stored_sam_recs(buf_size)
{
	string to_mapper = m_filtered_alignment_file; to_mapper.append(".tomapper");

	m_do_filter = do_filter;

	m_min_ins_static = m_min_ins;

	m_max_pair_dist = max_pair_dist;

	m_entrpy_weight = entrpy_weight;
	
	m_pqlen_weight = pqlen_weight;
	
	m_ave_mis_weight = ave_mis_weight;

	m_add_S = add_S;

	m_max_hits = max_hits;

	m_is_paired = is_paired; 

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0;
	
	m_multiple = 0;
	
	m_canoical = 0;
	
	m_semi_canonical = 0;
	
	m_non_canonical = 0;
	
	m_paired = 0;
	
	m_single = 0;

	m_unmapped = 0;
	
	m_fusion = 0;

	m_clipped = 0;
	
	//m_ofs_to_mapper.open(to_mapper.c_str());
}

bool AlignmentHandler::Init(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, 
	string filtered_alignment_file, size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor,
	double entrpy_weight, double pqlen_weight, double ave_mis_weight) 
{
	Clear();

	m_filtered_alignment_file = filtered_alignment_file;

	//m_ofs_filtered_alignment.open(m_filtered_alignment_file.c_str());

	string to_mapper = m_filtered_alignment_file; to_mapper.append(".tomapper");

	//m_ofs_to_mapper.open(to_mapper.c_str());

	m_entrpy_weight = entrpy_weight;

	m_pqlen_weight = pqlen_weight;

	m_ave_mis_weight = ave_mis_weight;

	m_is_paired = is_paired;

	m_add_S = add_S;

	m_max_pair_dist = max_pair_dist;

	m_max_hits = max_hits;

	//m_sam_file_handler.Init(alignment_file.c_str(), min_ins);

	m_sam_files = alignment_files;

	m_max_read_width = max_read_width;

	m_chrom_dir = chrom_dir;

	m_do_filter = do_filter;

	m_min_ins = min_ins;

	m_max_del = max_del;

	m_min_anchor = min_anchor;

	m_min_junc_anchor = min_junc_anchor;

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0;
	
	m_multiple = 0;
	
	m_canoical = 0;

	m_semi_canonical = 0;
	
	m_non_canonical = 0;
	
	m_paired = 0;
	
	m_single = 0;

	m_unmapped = 0;

	m_fusion = 0;

	m_clipped = 0;

	return true;
}

bool AlignmentHandler::Clear()
{
	m_is_paired = false;

	m_add_S = false;

	m_max_pair_dist = 0;

	m_max_pair_dist = 0;

	m_entrpy_weight = 0;

	m_pqlen_weight = 0;

	m_ave_mis_weight = 0;

	m_sam_files.clear();

	m_sam_file_handler.Clear();

	m_junction_handler.Clear();

	m_junc_hash_ptr = 0;

	m_sam_rec_pe.first.clear();

	m_sam_rec_pe.second.clear();

	m_sam_rec_pe_ptr.first.clear();

	m_sam_rec_pe_ptr.second.clear();

	m_filtered_alignment_file.clear();

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_bothunspliced.close();

	m_ofs_onespliced.close();

	m_chrom_dir.clear();

	m_do_filter = 0;

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0; 
	
	m_multiple = 0;
	
	m_canoical = 0;

	m_semi_canonical = 0;
	
	m_non_canonical = 0;

	

	return true;
}

bool AlignmentHandler::ClearStats()
{
	m_unmapped = 0;

	m_fusion = 0;

	m_clipped = 0;

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0; 
	
	m_multiple = 0;
	
	m_canoical = 0;

	m_semi_canonical = 0;
	
	m_non_canonical = 0;

	m_paired = 0;

	m_single = 0;

	return true;
}

void 
AlignmentHandler::FilterByMinAnchor(vector<SamRec>& read_sam)
{
	sort(read_sam.begin(), read_sam.end(), comp_min_anchor);

	vector<SamRec >::iterator sam_rec_iter1;

	for (sam_rec_iter1 = read_sam.begin(); sam_rec_iter1 != read_sam.end(); ++sam_rec_iter1)
	{
		if (sam_rec_iter1->min_anchor < m_min_anchor)
		{
			read_sam.resize(sam_rec_iter1 - read_sam.begin());

			break;
		}
	}
}

void 
AlignmentHandler::FilterByMapLen(vector<SamRec>& read_sam)
{
	sort(read_sam.begin(), read_sam.end(), comp_sam_rec_maplen);

	vector<SamRec >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin() + 1; sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		if (sam_rec_iter->mappedlen != (sam_rec_iter - 1)->mappedlen)
		{
			read_sam.resize(sam_rec_iter - read_sam.begin());

			break;
		}
	}
}

void
AlignmentHandler::MarkCanonNoncanonByReads(vector<SamRec>& read_sam, JunctionHandler* junc_hash_ptr)
{
	if (read_sam.empty())
		return;

	//FilterByMinAnchor(read_sam);

	//if (read_sam.empty())
	//	return;

	//FilterByMapLen(read_sam);

	vector<SamRec >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin(); sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		//cout << (*sam_rec_iter)->tostring()<<endl;
		if ((sam_rec_iter)->spliceway_vec.size() <= 1)
			continue;

		vector<pair<size_t, int> >::iterator oft_mpl_iter;

		double sumscore = 0;

		double sum_intron_len = 0;

		double junc_anchor_len = 0;

		double sum_ave_mis = 0;

		double pair_rate = 0;

		//cout << "calculating "<<endl;

		for (oft_mpl_iter = (sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
		{
			if ((oft_mpl_iter + 1)->second < 0)
			{
				//sam_rec_iter->is_insert = true;

				continue;
			}

			size_t comb_offset;

			if (oft_mpl_iter->second < 0)
			{
				comb_offset = ((oft_mpl_iter->first - 1) << THIRTY_TWO) + oft_mpl_iter->first - 1;
			}
			else
				comb_offset = ((oft_mpl_iter->first + oft_mpl_iter->second - 1) << THIRTY_TWO) + (oft_mpl_iter + 1)->first;

			CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = junc_hash_ptr->m_junc_hash.find(sam_rec_iter->chrom_name);

			if (chrom_junc_hash_iter == junc_hash_ptr->m_junc_hash.end())
			{
				//cout << "corresponding junction not found" << endl;

				//cout << sam_rec_iter->tostring((sam_rec_iter)->spliceway_vec.size(), oft_mpl_iter - (sam_rec_iter)->spliceway_vec.begin() + 1)<<endl;

				//cout << sam_rec_iter->chrom_name << '\t' << oft_mpl_iter->first + oft_mpl_iter->second - 1 << '\t' << (oft_mpl_iter + 1)->first << endl;

				sam_rec_iter->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				continue;
			}

			JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

			JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(comb_offset);

			if (junc_hash_comb_iter == junc_hash_comb.end())
			{
				//cout << "corresponding junction not found" << endl;

				//cout << sam_rec_iter->tostring((sam_rec_iter)->spliceway_vec.size(), oft_mpl_iter - (sam_rec_iter)->spliceway_vec.begin() + 1)<<endl;

				//cout << sam_rec_iter->chrom_name << '\t' << oft_mpl_iter->first + oft_mpl_iter->second - 1 << '\t' << (oft_mpl_iter + 1)->first << endl;

				sam_rec_iter->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				continue;
			}

			sam_rec_iter->filter_type= (junc_hash_comb_iter->second.m_filtered_type);

			//sam_rec_iter->corresponding_juncs.push_back(&junc_hash_comb_iter->second);
			
			//(sam_rec_iter)->junc_id.push_back(junc_hash_comb_iter->second.m_juncname);

			sumscore += (junc_hash_comb_iter->second.m_entropy * m_entrpy_weight) + (junc_hash_comb_iter->second.m_lpq * m_pqlen_weight) 
				+ (junc_hash_comb_iter->second.m_ave_mismatch * m_ave_mis_weight);

			sum_intron_len += junc_hash_comb_iter->second.m_intronlen;

			junc_anchor_len += junc_hash_comb_iter->second.m_max_prefix_len + junc_hash_comb_iter->second.m_max_suffix_len;

			sum_ave_mis += junc_hash_comb_iter->second.m_ave_mismatch;

			pair_rate += (double)junc_hash_comb_iter->second.m_paired_count;
			
			if (junc_hash_comb_iter->second.m_flankcase >= 5)
			{
				(sam_rec_iter)->canon_count++;
				(sam_rec_iter)->iscanonical = true;
			}
			else if (junc_hash_comb_iter->second.m_flankcase >= 1)
			{
				(sam_rec_iter)->canon_count++;

				(sam_rec_iter)->noncanon_count++;

				(sam_rec_iter)->issemicanonical = true;
			}
			else
			{
				(sam_rec_iter)->noncanon_count++;

				(sam_rec_iter)->isnoncanoical = true;
			}
		}

		sumscore = sumscore / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		sum_intron_len =  sum_intron_len / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		junc_anchor_len = junc_anchor_len / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		sum_ave_mis = sum_ave_mis / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		pair_rate = pair_rate / double ((sam_rec_iter)->spliceway_vec.size() - 1);
		//cout << "assigning "<<endl;
		(sam_rec_iter)->filter_score = sumscore;

		(sam_rec_iter)->junc_anchor_len = junc_anchor_len;

		(sam_rec_iter)->ave_intron_len = sum_intron_len;

		(sam_rec_iter)->ave_junc_mis = sum_ave_mis;

		(sam_rec_iter)->pair_rate = pair_rate;

		(sam_rec_iter)->canon_rate = (double)(sam_rec_iter)->canon_count / (double)((sam_rec_iter)->canon_count + (sam_rec_iter)->noncanon_count);

		//++count;
	}
}


void*
AlignmentHandler::MarkCanonNoncanonByReadsStatic(void * str)
{
	vector<SamRec*>& read_sam = *((vector<SamRec*>*)str);

	JunctionHandler* junc_hash_ptr = m_junction_handler_prev_ptr;

	if (read_sam.empty())
		return 0;

	//FilterByMinAnchor(read_sam);

	//if (read_sam.empty())
	//	return;

	//FilterByMapLen(read_sam);

	vector<SamRec* >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin(); sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		//cout << (*sam_rec_iter)->tostring()<<endl;

		if ((*sam_rec_iter)->is_fusion)
		{
			//if ((*sam_rec_iter)->tag_name.find("R_49939/") != string::npos)
			//	cout <<"mark canonical"<<endl;
			double sumscore = 0;

			double sum_intron_len = 0;

			double junc_anchor_len = 0;

			double junc_hits = 0;

			double sum_ave_mis = 0;

			double pair_rate = 0;

			size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

			string ins_str;

			prefixst = (*sam_rec_iter)->fusion_prefix_st;

			prefixend = (*sam_rec_iter)->fusion_prefix_end;

			prefixlen = (*sam_rec_iter)->fusion_prefix_len;//spliceway_vec.back().second;

			suffixst = (*sam_rec_iter)->fusion_suffix_st;

			suffixlen = (*sam_rec_iter)->fusion_suffix_len;
			
			//spliceway_vec2.front().second;

			//if (prefixend == 98658330 || suffixst == 98658330)
			//	cout << samrec.tostring(0, 0) << endl;

			combined_offset = (prefixend << THIRTY_TWO) + suffixst;

			string chrom_name;

			chrom_name = (*sam_rec_iter)->chrom_name + "~" + (*sam_rec_iter)->chrom_name2;

			CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = junc_hash_ptr->m_junc_hash.find(chrom_name);

			if (chrom_junc_hash_iter == junc_hash_ptr->m_junc_hash.end())
			{
				(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				continue;
			}

			JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

			JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

			if (junc_hash_comb_iter == junc_hash_comb.end())
			{
				(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				continue;
			}

			(*sam_rec_iter)->filter_type = (junc_hash_comb_iter->second.m_filtered_type);

			//(*sam_rec_iter)->corresponding_juncs.push_back(&junc_hash_comb_iter->second);

			sumscore += (junc_hash_comb_iter->second.m_entropy * m_entrpy_weight) + (junc_hash_comb_iter->second.m_lpq * m_pqlen_weight) 
				+ (junc_hash_comb_iter->second.m_ave_mismatch * m_ave_mis_weight);

			sum_intron_len += junc_hash_comb_iter->second.m_intronlen;

			junc_anchor_len += junc_hash_comb_iter->second.m_max_prefix_len + junc_hash_comb_iter->second.m_max_suffix_len;

			junc_hits += junc_hash_comb_iter->second.m_positive_count + junc_hash_comb_iter->second.m_negative_count;

			sum_ave_mis += junc_hash_comb_iter->second.m_ave_mismatch;

			pair_rate += (double)junc_hash_comb_iter->second.m_paired_count;

			if (junc_hash_comb_iter->second.m_flankcase >= 5)
			{
				(*sam_rec_iter)->canon_count++;
				(*sam_rec_iter)->iscanonical = true;
			}
			else if (junc_hash_comb_iter->second.m_flankcase >= 1)
			{
				(*sam_rec_iter)->canon_count++;

				(*sam_rec_iter)->noncanon_count++;

				(*sam_rec_iter)->issemicanonical = true;
			}
			else
			{
				(*sam_rec_iter)->noncanon_count++;

				(*sam_rec_iter)->isnoncanoical = true;
			}

			/*switch (junc_hash_comb_iter->second.m_flankcase) 
			{
			case 0:
				++(*sam_rec_iter)->forward_count;
				break;
			case 1:
				++(*sam_rec_iter)->forward_count;
				break;
			case 2:
				++(*sam_rec_iter)->reverse_count;
				break;
			case 3:
				++(*sam_rec_iter)->reverse_count;
				break;
			case 4:
				++(*sam_rec_iter)->forward_count;
				break;
			case 5:
				++(*sam_rec_iter)->forward_count;
				break;
			case 6:
				++(*sam_rec_iter)->reverse_count;
				break;
			default:
				cout << "none specified case: "<<junc_hash_comb_iter->second.m_flankcase << endl;
				break;
			}*/

			//sumscore = sumscore / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			//sum_intron_len =  sum_intron_len / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			//junc_anchor_len = junc_anchor_len / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			//sum_ave_mis = sum_ave_mis / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			//junc_hits = junc_hits / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			//pair_rate = pair_rate / double ((*sam_rec_iter)->spliceway_vec.size() - 1);
			//cout << "assigning "<<endl;
			(*sam_rec_iter)->filter_score = sumscore;

			(*sam_rec_iter)->junc_anchor_len = junc_anchor_len;

			(*sam_rec_iter)->ave_intron_len = 0;

			(*sam_rec_iter)->ave_junc_mis = sum_ave_mis;

			(*sam_rec_iter)->junc_hits = junc_hits;

			(*sam_rec_iter)->pair_rate = pair_rate;

			(*sam_rec_iter)->canon_rate = (double)(*sam_rec_iter)->canon_count / (double)((*sam_rec_iter)->canon_count + (*sam_rec_iter)->noncanon_count);

			if ((*sam_rec_iter)->tag_name.find("R_49939/") != string::npos)
			{
				cout <<"mark canonical"<<endl;

				cout <<"junc_hash_comb_iter->second.to_normal_junction(0):"<< junc_hash_comb_iter->second.to_normal_junction(0)<<endl;

				cout <<"(*sam_rec_iter)->filter_score:"<<(*sam_rec_iter)->filter_score<<endl;

				cout <<"(*sam_rec_iter)->junc_anchor_len:"<<(*sam_rec_iter)->junc_anchor_len<<endl;

				cout <<"(*sam_rec_iter)->ave_intron_len:"<<(*sam_rec_iter)->ave_intron_len<<endl;

				cout <<"(*sam_rec_iter)->ave_junc_mis:"<<(*sam_rec_iter)->ave_junc_mis<<endl;

				cout <<"(*sam_rec_iter)->junc_hits:"<<(*sam_rec_iter)->junc_hits<<endl;

				cout <<"(*sam_rec_iter)->pair_rate:"<<(*sam_rec_iter)->pair_rate<<endl;

				cout <<"(*sam_rec_iter)->canon_rate:"<<(*sam_rec_iter)->canon_rate<<endl;
				

			}

		}
		else
		{
			if ((*sam_rec_iter)->spliceway_vec.size() <= 1 || (*sam_rec_iter)->is_fusion)
				continue;

			////debug
			//if ((*sam_rec_iter)->tag_name.find("seq.1000169/") != string::npos)
			//	cout << "find seq.1000169/" << endl;

			vector<pair<size_t, int> >::iterator oft_mpl_iter;

			double sumscore = 0;

			double sum_intron_len = 0;

			double junc_anchor_len = 0;

			double junc_hits = 0;

			double sum_ave_mis = 0;

			double pair_rate = 0;

			//cout << "calculating "<<endl;

			for (oft_mpl_iter = (*sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (*sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
			{
				if ((oft_mpl_iter + 1)->second < 0)
				{
					//sam_rec_iter->is_insert = true;

					continue;
				}

				size_t comb_offset;

				if (oft_mpl_iter->second < 0)
				{
					comb_offset = ((oft_mpl_iter->first - 1) << THIRTY_TWO) + oft_mpl_iter->first - 1;
				}
				else
					comb_offset = ((oft_mpl_iter->first + oft_mpl_iter->second - 1) << THIRTY_TWO) + (oft_mpl_iter + 1)->first;

				CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = junc_hash_ptr->m_junc_hash.find((*sam_rec_iter)->chrom_name);

				if (chrom_junc_hash_iter == junc_hash_ptr->m_junc_hash.end())
				{
					//cout << "corresponding junction not found" << endl;

					//cout << sam_rec_iter->tostring((sam_rec_iter)->spliceway_vec.size(), oft_mpl_iter - (sam_rec_iter)->spliceway_vec.begin() + 1)<<endl;

					//cout << sam_rec_iter->chrom_name << '\t' << oft_mpl_iter->first + oft_mpl_iter->second - 1 << '\t' << (oft_mpl_iter + 1)->first << endl;

					(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

					continue;
				}

				JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

				JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(comb_offset);

				if (junc_hash_comb_iter == junc_hash_comb.end())
				{
					//cout << "corresponding junction not found" << endl;

					//cout << sam_rec_iter->tostring((sam_rec_iter)->spliceway_vec.size(), oft_mpl_iter - (sam_rec_iter)->spliceway_vec.begin() + 1)<<endl;

					//cout << sam_rec_iter->chrom_name << '\t' << oft_mpl_iter->first + oft_mpl_iter->second - 1 << '\t' << (oft_mpl_iter + 1)->first << endl;

					(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

					continue;
				}

				//if ((*sam_rec_iter)->tag_name.find("seq.1000169/") != string::npos || (*sam_rec_iter)->tag_name.find("seq.10004370/") != string::npos ||
				//	(*sam_rec_iter)->tag_name.find("seq.10012918/") != string::npos)
				//{
				//	cout << (*sam_rec_iter)->tostring(0, 0) << endl;

				//	cout << junc_hash_comb_iter->second.to_normal_junction(0) << endl;

				//	cout << "m_ave_mismatch\t"<< junc_hash_comb_iter->second.m_ave_mismatch << endl;

				//	cout << "hits\t"<< junc_hash_comb_iter->second.m_positive_count + junc_hash_comb_iter->second.m_negative_count << endl;
				//}

				(*sam_rec_iter)->filter_type = (junc_hash_comb_iter->second.m_filtered_type);

				//(*sam_rec_iter)->corresponding_juncs.push_back(&junc_hash_comb_iter->second);

				//(sam_rec_iter)->junc_id.push_back(junc_hash_comb_iter->second.m_juncname);

				sumscore += (junc_hash_comb_iter->second.m_entropy * m_entrpy_weight) + (junc_hash_comb_iter->second.m_lpq * m_pqlen_weight) 
					+ (junc_hash_comb_iter->second.m_ave_mismatch * m_ave_mis_weight);

				sum_intron_len += junc_hash_comb_iter->second.m_intronlen;

				junc_anchor_len += junc_hash_comb_iter->second.m_max_prefix_len + junc_hash_comb_iter->second.m_max_suffix_len;

				junc_hits += junc_hash_comb_iter->second.m_positive_count + junc_hash_comb_iter->second.m_negative_count;

				sum_ave_mis += junc_hash_comb_iter->second.m_ave_mismatch;

				pair_rate += (double)junc_hash_comb_iter->second.m_paired_count;

				if (junc_hash_comb_iter->second.m_flankcase >= 5)
				{
					(*sam_rec_iter)->canon_count++;
					(*sam_rec_iter)->iscanonical = true;
				}
				else if (junc_hash_comb_iter->second.m_flankcase >= 1)
				{
					(*sam_rec_iter)->canon_count++;

					(*sam_rec_iter)->noncanon_count++;

					(*sam_rec_iter)->issemicanonical = true;
				}
				else
				{
					(*sam_rec_iter)->noncanon_count++;

					(*sam_rec_iter)->isnoncanoical = true;
				}

				//switch (junc_hash_comb_iter->second.m_flankcase) 
				//{
				//case 0:
				//	++(*sam_rec_iter)->forward_count;
				//	break;
				//case 1:
				//	++(*sam_rec_iter)->forward_count;
				//	break;
				//case 2:
				//	++(*sam_rec_iter)->reverse_count;
				//	break;
				//case 3:
				//	++(*sam_rec_iter)->reverse_count;
				//	break;
				//case 4:
				//	++(*sam_rec_iter)->forward_count;
				//	break;
				//case 5:
				//	++(*sam_rec_iter)->forward_count;
				//	break;
				//case 6:
				//	++(*sam_rec_iter)->reverse_count;
				//	break;
				//default:
				//	cout << "none specified case: "<<junc_hash_comb_iter->second.m_flankcase << endl;
				//	break;
				//}
			}

			sumscore = sumscore / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			sum_intron_len =  sum_intron_len / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			junc_anchor_len = junc_anchor_len / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			sum_ave_mis = sum_ave_mis / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			junc_hits = junc_hits / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			pair_rate = pair_rate / double ((*sam_rec_iter)->spliceway_vec.size() - 1);
			//cout << "assigning "<<endl;
			(*sam_rec_iter)->filter_score = sumscore;

			(*sam_rec_iter)->junc_anchor_len = junc_anchor_len;

			(*sam_rec_iter)->ave_intron_len = sum_intron_len;

			(*sam_rec_iter)->ave_junc_mis = sum_ave_mis;

			(*sam_rec_iter)->junc_hits = junc_hits;

			(*sam_rec_iter)->pair_rate = pair_rate;

			(*sam_rec_iter)->canon_rate = (double)(*sam_rec_iter)->canon_count / (double)((*sam_rec_iter)->canon_count + (*sam_rec_iter)->noncanon_count);
		}

		//++count;
	}

	return 0;
}

void AlignmentHandler::LoadSamVec2SamPtr(vector<SamRec>& sam_rec, vector<SamRec*>& sam_rec_ptr)
{
	sam_rec_ptr.clear();

	vector<SamRec>::iterator sam_rec_iter;

	for (sam_rec_iter = sam_rec.begin(); sam_rec_iter != sam_rec.end(); ++sam_rec_iter)
	{
		sam_rec_ptr.push_back(&(*sam_rec_iter));
	}

	if (sam_rec_ptr.size() > 1)
		RemoveDup(sam_rec_ptr);
}

void AlignmentHandler::ProcessPEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	//cout << 1 << endl;
	//if (m_do_filter & 4)
	//{
	//	MarkCanonNoncanonByReads(m_sam_rec_pe.first, junc_handler_filtered);

	//	MarkCanonNoncanonByReads(m_sam_rec_pe.second, junc_handler_filtered);
	//}
	//cout << "start filtering" << endl;
	//SetSamRecUnpaired(m_sam_rec_pe);

	//cout << 2 << endl;

	//LoadSamVec2SamPtr(m_sam_rec_pe.first, m_sam_rec_pe_ptr.first);

	//LoadSamVec2SamPtr(m_sam_rec_pe.second, m_sam_rec_pe_ptr.second);
	
	//1

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best

	//cout << "filtered junc" << endl;

	//cout << 3 << endl;

	//if (m_do_filter & 4)// do junction filtering clipping
	//{
	//	FilterByFilteredJunction(m_sam_rec_pe_ptr.first, junc_handler_filtered);

	//	FilterByFilteredJunction(m_sam_rec_pe_ptr.second, junc_handler_filtered);
	//}

	//cout << 4 << endl;

	//FindSamVecContigLen(m_sam_rec_pe_ptr.first);

	//FindSamVecContigLen(m_sam_rec_pe_ptr.second);

	//2

	//bool paired = false;

	//cout << 5 << endl;

	//vector<PairedSamRec> fusion_paired_reads_ptr;

	//if (m_sam_rec_pe_ptr.first.size() && m_sam_rec_pe_ptr.second.size())
	//{
	//	//2. filter by pairing 

	//	paired = EstablishPairing(m_sam_rec_pe_ptr, fusion_paired_reads_ptr);
	//}

	//FindFusionJuncRegionVec(junc_handler_filtered, fusion_paired_reads_ptr);


	//3. filter multiple

	//cout << 6 << endl;

	//if ((m_do_filter & 8) && paired == false) //do filter multiple unpaired
	//{
	//	FilterSingleMulti(m_sam_rec_pe_ptr.first);

	//	FilterSingleMulti(m_sam_rec_pe_ptr.second);
	//}

	////cout << "convert 2 junc" << endl;

	////cout << 7 << endl;

	//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);

	//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.second);

	//CollectStats(m_sam_rec_pe_ptr.first);

	//CollectStats(m_sam_rec_pe_ptr.second);

	////cout << 8 << endl;

	////cout << "write alignment" << endl;

	//if (paired == false)
	//{
	//	SetFusionBitInfo(m_sam_rec_pe_ptr);

	//	if (m_sam_rec_pe_ptr.first.size() && m_sam_rec_pe_ptr.second.size())
	//	{
	//		WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_fusion_paired);

	//		WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_fusion_paired);
	//	}
	//	else
	//	{
	//		if (m_sam_rec_pe_ptr.first.size())
	//		{
	//			for (size_t i = 0; i < m_sam_rec_pe_ptr.first.size(); ++i)
	//			{
	//				if (m_sam_rec_pe_ptr.first[i]->is_fusion)
	//					m_sam_rec_pe_ptr.first[i]->strand_t |= MATE_UNMAPPED;
	//			}

	//			WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
	//		}
	//		else if (m_sam_rec_pe_ptr.second.size())
	//		{
	//			for (size_t i = 0; i < m_sam_rec_pe_ptr.second.size(); ++i)
	//			{
	//				if (m_sam_rec_pe_ptr.second[i]->is_fusion)
	//					m_sam_rec_pe_ptr.second[i]->strand_t |= MATE_UNMAPPED;
	//			}

	//			WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_single);
	//		}
	//		else
	//			;
	//	}
	//}

	//SetBitInfo(m_sam_rec_pe_ptr);

	//WriteAlignment(m_sam_rec_pe_ptr.first);

	//WriteAlignment(m_sam_rec_pe_ptr.second);

	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			cout << "parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_base_name;

			cout << "group same read alignments"<<endl;

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_base_name)
				{
					//cout << prev_tag << endl;

					//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

					//cout << (sam_iter - 1)->tostring(0, 0) << endl;

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_base_name;

					previous_iter = sam_iter;
				}
			}

			//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			cout << "process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadFilterJuncMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			
			cout << "output to files"<<endl;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				//cout << "SamRecVec2Junc"<<endl;

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				//cout << "CollectStatsStatic"<<endl;

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

				//cout << "WriteAlignment"<<endl;

				WriteAlignment(m_same_reads_alignments_paired_iter->first);

				WriteAlignment(m_same_reads_alignments_paired_iter->second);
				//CollectStats(m_same_reads_alignments_paired_iter->first);

				//CollectStats(m_same_reads_alignments_paired_iter->second);
			}

			//cout << "junction size:" << junc_handler->m_total_junc << endl;

			//cout << "junction size:" << junc_handler->m_total_encompass << endl;

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			cout << "remained lines:"<< index << endl;

			cout << "finish"<<endl;

			//keep the last read's alignments to be processed with next block of alignments

			//m_stored_sam_recs.clear();

			//m_stored_lines.clear();
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{
		//cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_base_name;

		cout << "group same read alignments"<<endl;

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_base_name)
			{
				//cout << prev_tag << endl;

				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_base_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		cout << "process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadFilterJuncMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

			WriteAlignment(m_same_reads_alignments_paired_iter->first);

			WriteAlignment(m_same_reads_alignments_paired_iter->second);
			//CollectStats(m_same_reads_alignments_paired_iter->first);

			//CollectStats(m_same_reads_alignments_paired_iter->second);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		cout << "finish"<<endl;

		//keep the last read's alignments to be processed with next block of alignments

		//m_stored_sam_recs.clear();

		//m_stored_lines.clear();
	}

	cout << "processed_line:"<< processed_line << endl;

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void AlignmentHandler::ProcessPERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			//cout << m_stored_lines.front() << endl;
			//cout << m_stored_lines.back() << endl;
			cout << "parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			//if (m_prev_tag.empty())
			m_prev_tag = sam_iter->tag_base_name;

			cout << "group same read alignments"<<endl;

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_base_name)
				{
					//cout << m_prev_tag << endl;

					//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

					//cout << (sam_iter - 1)->tostring(0, 0) << endl;

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_base_name;

					previous_iter = sam_iter;
				}
			}

				//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

			cout << "process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				//cout << "sam2junction"<<endl;

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

				//cout << "write alignment" << endl;

				WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
				//CollectStats(m_same_reads_alignments_paired_iter->first);

				//CollectStats(m_same_reads_alignments_paired_iter->second);
			}


			//for (sam_iter = m_stored_sam_recs.begin(); sam_iter != m_stored_sam_recs.end(); ++sam_iter)
			//{
			//	junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);
			//}

			//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.second);

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			m_same_reads_alignments_paired_is_unspliced.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				//cout << processed_lines << endl;
				//cout << processed_lines + i << endl;
				//cout << m_stored_lines.size() << endl;
				//cout << m_stored_lines[processed_lines + i] << endl;
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			cout << "remain line:"<< index << endl;

			cout << "finish"<<endl;
			//keep the last read's alignments to be processed with next block of alignments

			//m_stored_sam_recs.clear();

			//m_stored_lines.clear();
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}


	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	cout << "process the end block of file" << endl;

	if (m_stored_sam_recs.size())
	{
		cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		//if (m_prev_tag.empty())
		m_prev_tag = sam_iter->tag_base_name;

		cout << "group same read alignments"<<endl;

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_base_name)
			{
				//cout << prev_tag << endl;

				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_base_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

		cout << "process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

			WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
			//CollectStats(m_same_reads_alignments_paired_iter->first);

			//CollectStats(m_same_reads_alignments_paired_iter->second);
		}


		//for (sam_iter = m_stored_sam_recs.begin(); sam_iter != m_stored_sam_recs.end(); ++sam_iter)
		//{
		//	junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);
		//}

		//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.second);

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		m_same_reads_alignments_paired_is_unspliced.clear();
	}

	cout << "processed_line:"<< processed_line << endl;

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired_is_unspliced.clear();

	m_same_reads_alignments_paired.clear();

	m_prev_tag = "";

	index = 0;
}

void AlignmentHandler::ProcessPEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			cout << "parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_base_name;

			cout << "group same read alignments"<<endl;

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_base_name)
				{
					//cout << prev_tag << endl;

					//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

					//cout << (sam_iter - 1)->tostring(0, 0) << endl;

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_base_name;

					previous_iter = sam_iter;
				}
			}

			//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			cout << "process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadSelectBestMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				//cout << "SamRecVec2Junc"<<endl;

				//if (m_same_reads_alignments_paired_iter->first.size())
				//	cout << m_same_reads_alignments_paired_iter->first.front()->tostring(0, 0) << endl;

				//if (m_same_reads_alignments_paired_iter->second.size())
				//	cout << m_same_reads_alignments_paired_iter->second.front()->tostring(0, 0) << endl;

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				//cout << "CollectStatsStatic"<<endl;

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

				//cout << "WriteAlignment"<<endl;

				//WriteAlignment(m_same_reads_alignments_paired_iter->first);

				//WriteAlignment(m_same_reads_alignments_paired_iter->second);
				//CollectStats(m_same_reads_alignments_paired_iter->first);

				//CollectStats(m_same_reads_alignments_paired_iter->second);
			}

			//cout << "junction size:" << junc_handler->m_total_junc << endl;

			//cout << "junction size:" << junc_handler->m_total_encompass << endl;

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			cout << "remained lines:"<< index << endl;

			cout << "finish"<<endl;

			//keep the last read's alignments to be processed with next block of alignments

			//m_stored_sam_recs.clear();

			//m_stored_lines.clear();
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{
		//cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_base_name;

		cout << "group same read alignments"<<endl;

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_base_name)
			{
				//cout << prev_tag << endl;

				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_base_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		cout << "process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadSelectBestMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

			//WriteAlignment(m_same_reads_alignments_paired_iter->first);

			//WriteAlignment(m_same_reads_alignments_paired_iter->second);
			//CollectStats(m_same_reads_alignments_paired_iter->first);

			//CollectStats(m_same_reads_alignments_paired_iter->second);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		cout << "finish"<<endl;

		//keep the last read's alignments to be processed with next block of alignments

		//m_stored_sam_recs.clear();

		//m_stored_lines.clear();
	}

	cout << "processed_line:"<< processed_line << endl;

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void AlignmentHandler::ProcessSEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			cout << "parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_name;

			cout << "group same read alignments"<<endl;

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_name)
				{
					//cout << prev_tag << endl;

					//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

					//cout << (sam_iter - 1)->tostring(0, 0) << endl;

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_name;

					previous_iter = sam_iter;
				}
			}

			//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			cout << "process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadFilterJuncMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				//junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				WriteAlignment(m_same_reads_alignments_paired_iter->first);

				//WriteAlignment(m_same_reads_alignments_paired_iter->second);
				//CollectStats(m_same_reads_alignments_paired_iter->first);

				//CollectStats(m_same_reads_alignments_paired_iter->second);
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			//cout << "remained lines:"<< index << endl;

			cout << "finish"<<endl;

			//keep the last read's alignments to be processed with next block of alignments

			//m_stored_sam_recs.clear();

			//m_stored_lines.clear();
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{
		//cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_name;

		cout << "group same read alignments"<<endl;

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_name)
			{
				//cout << prev_tag << endl;

				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		cout << "process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadFilterJuncMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));
			//junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			WriteAlignment(m_same_reads_alignments_paired_iter->first);

			//WriteAlignment(m_same_reads_alignments_paired_iter->second);
			//CollectStats(m_same_reads_alignments_paired_iter->first);

			//CollectStats(m_same_reads_alignments_paired_iter->second);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		cout << "finish"<<endl;

		//keep the last read's alignments to be processed with next block of alignments

		//m_stored_sam_recs.clear();

		//m_stored_lines.clear();
	}

	cout << "processed_line:"<< processed_line << endl;

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void AlignmentHandler::ProcessSERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	//LoadSamVec2SamPtr(m_sam_rec_pe.first, m_sam_rec_pe_ptr.first);

	////2

	//vector<SamRec*>::iterator sam_rec_iter;

	//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);

	//CollectStats(m_sam_rec_pe_ptr.first);


	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			//cout << m_stored_lines.front() << endl;
			//cout << m_stored_lines.back() << endl;
			cout << "parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			//if (m_prev_tag.empty())
			m_prev_tag = sam_iter->tag_name;

			cout << "group same read alignments"<<endl;

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_name)
				{
					//cout << prev_tag << endl;

					//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

					//cout << (sam_iter - 1)->tostring(0, 0) << endl;

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_name;

					previous_iter = sam_iter;
				}
			}

			cout << "m_stored_sam_recs size:" << m_stored_sam_recs.size() << endl;

			cout << "last same read index:" << previous_iter - m_stored_sam_recs.begin() << endl;

			//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

			cout << "process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
				//junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				//CollectStats(m_same_reads_alignments_paired_iter->first);

				//CollectStats(m_same_reads_alignments_paired_iter->second);
			}


			//for (sam_iter = m_stored_sam_recs.begin(); sam_iter != m_stored_sam_recs.end(); ++sam_iter)
			//{
			//	junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);
			//}

			//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.second);

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			m_same_reads_alignments_paired_is_unspliced.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				//cout << processed_lines << endl;
				//cout << processed_lines + i << endl;
				//cout << m_stored_lines.size() << endl;
				//cout << m_stored_lines[processed_lines + i] << endl;
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			cout << "remain line:"<< index << endl;

			cout << "finish"<<endl;
			//keep the last read's alignments to be processed with next block of alignments

			//m_stored_sam_recs.clear();

			//m_stored_lines.clear();
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}


	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	cout << "process the end block of file" << endl;

	if (m_stored_sam_recs.size())
	{
		cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		//if (m_prev_tag.empty())
		m_prev_tag = sam_iter->tag_name;

		cout << "group same read alignments"<<endl;

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_name)
			{
				//cout << prev_tag << endl;

				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

		cout << "process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
			//junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			//CollectStats(m_same_reads_alignments_paired_iter->first);

			//CollectStats(m_same_reads_alignments_paired_iter->second);
		}


		//for (sam_iter = m_stored_sam_recs.begin(); sam_iter != m_stored_sam_recs.end(); ++sam_iter)
		//{
		//	junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);
		//}

		//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.second);

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		m_same_reads_alignments_paired_is_unspliced.clear();
	}

	cout << "processed_line:"<< processed_line << endl;

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	m_same_reads_alignments_paired_is_unspliced.clear();

	m_prev_tag = "";

	index = 0;
}

void AlignmentHandler::ProcessSEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			cout << "parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait parse"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_name;

			cout << "group same read alignments"<<endl;

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_name)
				{
					//cout << prev_tag << endl;

					//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

					//cout << (sam_iter - 1)->tostring(0, 0) << endl;

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_name;

					previous_iter = sam_iter;
				}
			}

			//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			cout << "process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadSelectBestMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			cout << "wait process same read alignments"<<endl;

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				//junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				//WriteAlignment(m_same_reads_alignments_paired_iter->first);

				//WriteAlignment(m_same_reads_alignments_paired_iter->second);
				//CollectStats(m_same_reads_alignments_paired_iter->first);

				//CollectStats(m_same_reads_alignments_paired_iter->second);
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			//cout << "remained lines:"<< index << endl;

			cout << "finish"<<endl;

			//keep the last read's alignments to be processed with next block of alignments

			//m_stored_sam_recs.clear();

			//m_stored_lines.clear();
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{
		//cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait parse"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_name;

		cout << "group same read alignments"<<endl;

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_name)
			{
				//cout << prev_tag << endl;

				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		cout << "process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadSelectBestMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		cout << "wait process same read alignments"<<endl;

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));
			//junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			//WriteAlignment(m_same_reads_alignments_paired_iter->first);

			//WriteAlignment(m_same_reads_alignments_paired_iter->second);
			//CollectStats(m_same_reads_alignments_paired_iter->first);

			//CollectStats(m_same_reads_alignments_paired_iter->second);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		cout << "finish"<<endl;

		//keep the last read's alignments to be processed with next block of alignments

		//m_stored_sam_recs.clear();

		//m_stored_lines.clear();
	}

	cout << "processed_line:"<< processed_line << endl;

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void* AlignmentHandler::ParseLines(void * str)
{
	while(true) 
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_line_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_stored_lines.size())
		{
			break;
		}
		//cout << "cur read index:" << read_index << endl;

		//cout << m_stored_lines[read_index] << endl;

		//if (read_index == 999999)
		//	break;

		m_stored_sam_recs[read_index].Set(m_stored_lines[read_index], 3);

		//pthread_mutex_unlock(&inc_num_threads);
	}

	return 0;
}

void* AlignmentHandler::ProcessPEReadMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;
		//pthread_mutex_lock(&inc_num_threads);

		//cout << "lock"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			//pthread_mutex_unlock(&inc_num_threads);
			break;
		}

		
		//cout << "locked"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur read indx:"<<read_index << endl;

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur size:"<< cur_same_read_alignments.second - cur_same_read_alignments.first + 1 << endl;

		//cout << "last index:"<< cur_same_read_alignments.second - m_same_reads_alignments.begin()->first << endl;

		//cout << cur_same_read_alignments.second->tostring(0, 0) << endl;

		bool are_both_ends_unspliced = true;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			//cout << sam_rec_iter->tostring(0, 0) << endl;
			//cout << sam_rec_iter->isexonic << endl;
			//cout << sam_rec_iter->is_fusion << endl;
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if ((sam_rec_iter->isexonic == false && sam_rec_iter->is_fusion == false) || (sam_rec_iter->is_fusion == true))
					are_both_ends_unspliced = false;

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}
			}
			else
			{
				are_both_ends_unspliced = false;

				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					cout << "new format fusion alignment not together:"<<endl <<sam_rec_iter->tostring(0, 0 ) << endl;

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}

				++sam_rec_iter;
			}
		}

		int tmp_do_filter = m_do_filter;

		if (are_both_ends_unspliced)
			tmp_do_filter = tmp_do_filter | 16;

		//cout << "remove duplication 1"<<endl;

		if (cur_same_read_alignments_paired.first.size())
			RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);

		//cout << "remove duplication 2"<<endl;
		if (cur_same_read_alignments_paired.second.size())
			RemoveDupStatic((void*)&cur_same_read_alignments_paired.second);

		SetSamRecUnpaired(cur_same_read_alignments_paired);

		vector<PairedSamRec> fusion_paired_reads;

		//cout << " EstablishPairingStatic " << endl;

		bool paired = false;

		EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads, tmp_do_filter);

		paired = EstablishPairingStatic((void *) &est_pairing);

		if (/*(m_do_filter & 8) &&*/ paired == false && are_both_ends_unspliced) //do filter multiple unpaired
		{
			if ((m_do_filter & 256) == 0)
			{
				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
			}
		}

		if (are_both_ends_unspliced)
		{
			m_same_reads_alignments_paired_is_unspliced[read_index] = true;
			//Set bits
		}

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
				//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_fusion_paired);

				//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_fusion_paired);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.first[i]->is_fusion)
							cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
				}
				else if (cur_same_read_alignments_paired.second.size())
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.second[i]->is_fusion)
							cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_single);
				}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));

		//cout << "are_both_ends_unspliced:"<<are_both_ends_unspliced << endl;
		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first), are_both_ends_unspliced);

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.second), are_both_ends_unspliced);
		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair finished"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);
		//pthread_mutex_unlock(&inc_num_threads);

	}  

	return 0;
}

void* AlignmentHandler::ProcessSEReadMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;
		//pthread_mutex_lock(&inc_num_threads);

		//cout << "lock"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			//pthread_mutex_unlock(&inc_num_threads);
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		bool are_both_ends_unspliced = true;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if (sam_rec_iter->isexonic == false && sam_rec_iter->is_fusion == false)
					are_both_ends_unspliced = false;
			}

			cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
		}

		//cout << cur_same_read_alignments_paired.first.size()<< endl;

		//cout << cur_same_read_alignments_paired.first.front()->tostring(0, 0)<< endl;

		//cout << "remove duplicated "<< endl;

		int tmp_do_filter = m_do_filter;

		if (are_both_ends_unspliced)
			tmp_do_filter = tmp_do_filter | 16;

		if (cur_same_read_alignments_paired.first.size() > 0)
			RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);
	
		//cout << cur_same_read_alignments_paired.first.size()<< endl;

		//cout << cur_same_read_alignments_paired.first.front()->tostring(0, 0)<< endl;

		bool paired = false;

		//cout << "filter single" << endl;

		if (paired == false && are_both_ends_unspliced) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

			//FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		//cout << cur_same_read_alignments_paired.first.size()<< endl;

		//cout << cur_same_read_alignments_paired.first.front()->tostring(0, 0)<< endl;

		//cout << "set both_ends_unspliced" << endl;

		if (are_both_ends_unspliced)
		{
			m_same_reads_alignments_paired_is_unspliced[read_index] = true;
			//Set bits
		}

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
				}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);
		}

		//SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));
	
		//cout << "GenerateAlignment" << endl;

		//cout << cur_same_read_alignments_paired.first.front()->tostring(0, 0)<< endl;

		//cout << "is fusion:" << cur_same_read_alignments_paired.first.front()->is_fusion << endl;

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first), are_both_ends_unspliced);

		//cout << "put in pair finished"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

	}  

	return 0;
}

void* AlignmentHandler::ProcessPEReadFilterJuncMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;
		//pthread_mutex_lock(&inc_num_threads);

		//cout << "lock"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			//pthread_mutex_unlock(&inc_num_threads);
			break;
		}

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "locked"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur read indx:"<<read_index << endl;

		//cout << m_same_reads_alignments[read_index].first->tag_base_name<<endl;

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur size:"<< cur_same_read_alignments.second - cur_same_read_alignments.first + 1 << endl;

		//cout << "last index:"<< cur_same_read_alignments.second - m_same_reads_alignments.begin()->first << endl;

		//cout << cur_same_read_alignments.second->tostring(0, 0) << endl;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}

				++sam_rec_iter;
			}
		}

		//cout << "mark canon noncanon"<<endl;

		//cout << cur_same_read_alignments_paired.first.size() << endl;

		//cout << cur_same_read_alignments_paired.second.size() << endl;

		SetSamRecUnpaired(cur_same_read_alignments_paired);

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.second));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		bool paired = false;

		vector<PairedSamRec> fusion_paired_reads;

		//cout << "est pairing"<<endl;

		//cout << cur_same_read_alignments_paired.first.size() << endl;

		//cout << cur_same_read_alignments_paired.second.size() << endl;

		EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads, m_do_filter);

		paired = EstablishPairingStatic((void *) &est_pairing);

		//cout << "is_paired:" << paired << endl;

		//cout << "FindFusionJuncRegionVecStatic"<<endl;

		//cout << cur_same_read_alignments_paired.first.size() << endl;

		//cout << cur_same_read_alignments_paired.second.size() << endl;

		//if (paired == false)
		//	FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);

		//cout << "filter single"<<endl;

		//cout << cur_same_read_alignments_paired.first.size() << endl;

		//cout << cur_same_read_alignments_paired.second.size() << endl;

		if (/*(m_do_filter & 8) &&*/ paired == false) //do filter multiple unpaired
		{
			//RemoveFusionUnPairedStatic((void*)&(cur_same_read_alignments_paired.first));

			//RemoveFusionUnPairedStatic((void*)&(cur_same_read_alignments_paired.second));

			if ((m_do_filter & 256) == 0)
			{
				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
			}

			EstablishFusionPairingStatic((void *) &est_pairing);

			FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);
		}

		//CollectStatsStatic((void*)&(cur_same_read_alignments_paired.first));

		//CollectStatsStatic((void*)&(cur_same_read_alignments_paired.second));

		//cout << "set bit"<<endl;

		//cout << cur_same_read_alignments_paired.first.size() << endl;

		//cout << cur_same_read_alignments_paired.second.size() << endl;

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
				//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_fusion_paired);

				//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_fusion_paired);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.first[i]->is_fusion)
							cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
				}
				else if (cur_same_read_alignments_paired.second.size())
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.second[i]->is_fusion)
							cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_single);
				}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));
	
		//cout << "GenerateAlignment"<<endl;

		//cout << cur_same_read_alignments_paired.first.size() << endl;

		//cout << cur_same_read_alignments_paired.second.size() << endl;

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first));

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.second));
		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair finished"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);
		//pthread_mutex_unlock(&inc_num_threads);

	}  

	return 0;
}

void* AlignmentHandler::ProcessPEReadSelectBestMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;
		//pthread_mutex_lock(&inc_num_threads);

		//cout << "lock"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			//pthread_mutex_unlock(&inc_num_threads);
			break;
		}

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "locked"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur read indx:"<<read_index << endl;

		//cout << m_same_reads_alignments[read_index].first->tag_base_name<<endl;

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur size:"<< cur_same_read_alignments.second - cur_same_read_alignments.first + 1 << endl;

		//cout << "last index:"<< cur_same_read_alignments.second - m_same_reads_alignments.begin()->first << endl;

		//cout << cur_same_read_alignments.second->tostring(0, 0) << endl;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}

				++sam_rec_iter;
			}
		}

		//cout << "mark canon noncanon"<<endl;

		SetSamRecUnpaired(cur_same_read_alignments_paired);

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.second));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		bool paired = false;

		vector<PairedSamRec> fusion_paired_reads;

		//cout << "est pairing"<<endl;

		EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads, m_do_filter);

		paired = EstablishPairingStatic((void *) &est_pairing);

		//FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);

		//cout << "filter single"<<endl;

		if (/*(m_do_filter & 8) &&*/ paired == false) //do filter multiple unpaired
		{
			//RemoveFusionUnPairedStatic((void*)&(cur_same_read_alignments_paired.first));

			//RemoveFusionUnPairedStatic((void*)&(cur_same_read_alignments_paired.second));

			if ((m_do_filter & 256) == 0)
			{
				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
			}

			//EstablishFusionPairingStatic((void *) &est_pairing);

			//FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);
		}

		//CollectStatsStatic((void*)&(cur_same_read_alignments_paired.first));

		//CollectStatsStatic((void*)&(cur_same_read_alignments_paired.second));

		//cout << "set bit"<<endl;

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
				//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_fusion_paired);

				//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_fusion_paired);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.first[i]->is_fusion)
							cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
				}
				else if (cur_same_read_alignments_paired.second.size())
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.second[i]->is_fusion)
							cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_single);
				}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));
	
		//GenerateAlignment((void *)&(cur_same_read_alignments_paired.first));

		//GenerateAlignment((void *)&(cur_same_read_alignments_paired.second));
		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair finished"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);
		//pthread_mutex_unlock(&inc_num_threads);

	}  

	return 0;
}

void* AlignmentHandler::ProcessSEReadFilterJuncMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;
		//pthread_mutex_lock(&inc_num_threads);

		//cout << "lock"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			//pthread_mutex_unlock(&inc_num_threads);
			break;
		}

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "locked"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur read indx:"<<read_index << endl;

		//cout << m_same_reads_alignments[read_index].first->tag_base_name<<endl;

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur size:"<< cur_same_read_alignments.second - cur_same_read_alignments.first + 1 << endl;

		//cout << "last index:"<< cur_same_read_alignments.second - m_same_reads_alignments.begin()->first << endl;

		//cout << cur_same_read_alignments.second->tostring(0, 0) << endl;

		//for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		//{
		//	//if (sam_rec_iter->end_id == 1)
		//	//{
		//	cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
		//	//}
		//	//else if (sam_rec_iter->end_id == 2)
		//	//{
		//	//	cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
		//	//}
		//}

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));

				++sam_rec_iter;
			}
		}

		RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);
		//cout << "mark canon noncanon"<<endl;

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			//MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.second));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));

			//FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		bool paired = false;

		//vector<PairedSamRec> fusion_paired_reads;

		////cout << "est pairing"<<endl;

		//EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads);

		//paired = EstablishPairingStatic((void *) &est_pairing);

		//FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);

		//cout << "filter single"<<endl;

		if (/*(m_do_filter & 8) &&*/ paired == false) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

			//FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		

		//CollectStatsStatic((void*)&(cur_same_read_alignments_paired.second));

		//cout << "set bit"<<endl;

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				//SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
				//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_fusion_paired);

				//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_fusion_paired);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					//for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					//{
					//	if (!cur_same_read_alignments_paired.first[i]->is_fusion)
					//		cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					//}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
				}
				//else if (cur_same_read_alignments_paired.second.size())
				//{
				//	for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
				//	{
				//		if (cur_same_read_alignments_paired.second[i]->is_fusion)
				//			cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
				//	}

				//	SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
				//	//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_single);
				//}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			//SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		//SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));
	
		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first));

		//GenerateAlignment((void *)&(cur_same_read_alignments_paired.second));
		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair finished"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);
		//pthread_mutex_unlock(&inc_num_threads);

	}  

	return 0;
}

void* AlignmentHandler::ProcessSEReadSelectBestMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;
		//pthread_mutex_lock(&inc_num_threads);

		//cout << "lock"<<endl;

		//pthread_mutex_unlock(&inc_num_threads);

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			//pthread_mutex_unlock(&inc_num_threads);
			break;
		}

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "locked"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur read indx:"<<read_index << endl;

		//cout << m_same_reads_alignments[read_index].first->tag_base_name<<endl;

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);

		//cout << "cur size:"<< cur_same_read_alignments.second - cur_same_read_alignments.first + 1 << endl;

		//cout << "last index:"<< cur_same_read_alignments.second - m_same_reads_alignments.begin()->first << endl;

		//cout << cur_same_read_alignments.second->tostring(0, 0) << endl;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));

				++sam_rec_iter;
			}
		}

		RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);
		//cout << "mark canon noncanon"<<endl;

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			//MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.second));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));

			//FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		bool paired = false;

		//vector<PairedSamRec> fusion_paired_reads;

		////cout << "est pairing"<<endl;

		//EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads);

		//paired = EstablishPairingStatic((void *) &est_pairing);

		//FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);

		//cout << "filter single"<<endl;

		if (/*(m_do_filter & 8) &&*/ paired == false) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

			//FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		

		//CollectStatsStatic((void*)&(cur_same_read_alignments_paired.second));

		//cout << "set bit"<<endl;

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				//SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
				//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_fusion_paired);

				//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_fusion_paired);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					//for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					//{
					//	if (!cur_same_read_alignments_paired.first[i]->is_fusion)
					//		cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					//}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
					//WriteAlignment(m_sam_rec_pe_ptr.first, m_ofs_single);
				}
				//else if (cur_same_read_alignments_paired.second.size())
				//{
				//	for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
				//	{
				//		if (cur_same_read_alignments_paired.second[i]->is_fusion)
				//			cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
				//	}

				//	SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
				//	//WriteAlignment(m_sam_rec_pe_ptr.second, m_ofs_single);
				//}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			//SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		//SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));
	
		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first));

		//GenerateAlignment((void *)&(cur_same_read_alignments_paired.second));
		//pthread_mutex_lock(&inc_num_threads);
		//cout << "put in pair finished"<<endl;
		//pthread_mutex_unlock(&inc_num_threads);
		//pthread_mutex_unlock(&inc_num_threads);

	}  

	return 0;
}

void AlignmentHandler::ProcessFile(string alignment_file, JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered, int fileter_by_junc)
{
	m_junction_handler_prev_ptr = junc_handler_filtered;

	
	junc_handler->ClearJunction();

	if (fileter_by_junc == 0)
	{
		if (m_do_filter & 32)
		{
			if (!m_ofs_unspliced_fusion_paired.is_open() && !m_ofs_unspliced_single.is_open() && !m_ofs_unspliced_paired.is_open() && !m_ofs_onespliced.is_open())
			{
				string bothunspliced_fusion_paired = alignment_file; bothunspliced_fusion_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_fusion_paired.append(".bothunspliced.fusion_paired");

				m_ofs_unspliced_fusion_paired.open(bothunspliced_fusion_paired.c_str());

				string bothunspliced_single = alignment_file; bothunspliced_single.append(basename2(m_filtered_alignment_file)); bothunspliced_single.append(".bothunspliced.single");

				m_ofs_unspliced_single.open(bothunspliced_single.c_str());

				string bothunspliced_paired = alignment_file; bothunspliced_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_paired.append(".bothunspliced.paired");

				m_ofs_unspliced_paired.open(bothunspliced_paired.c_str());

				string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

				m_ofs_onespliced.open(onespliced.c_str());

				m_ifs_alignments.open(alignment_file.c_str());
			}
		}
		else
		{
			m_ofs_unspliced_fusion_paired.close();

			m_ofs_unspliced_single.close();

			m_ofs_unspliced_paired.close();

			m_ofs_onespliced.close();

			string bothunspliced_fusion_paired = alignment_file; bothunspliced_fusion_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_fusion_paired.append(".bothunspliced.fusion_paired");

			m_ofs_unspliced_fusion_paired.open(bothunspliced_fusion_paired.c_str());

			string bothunspliced_single = alignment_file; bothunspliced_single.append(basename2(m_filtered_alignment_file)); bothunspliced_single.append(".bothunspliced.single");

			m_ofs_unspliced_single.open(bothunspliced_single.c_str());

			string bothunspliced_paired = alignment_file; bothunspliced_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_paired.append(".bothunspliced.paired");

			m_ofs_unspliced_paired.open(bothunspliced_paired.c_str());

			string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

			m_ofs_onespliced.open(onespliced.c_str());

			m_ifs_alignments.open(alignment_file.c_str());
		}
	}
	else
	{
		if (m_do_filter & 32)
		{
			if (!m_ofs_filtered_alignment.is_open() && 
				!m_ofs_fusion_paired.is_open() && 
				!m_ofs_single.is_open() && 
				!m_ofs_fusion_std.is_open() && 
				!m_ofs_paired.is_open())
			{
				string to_mapper_str = alignment_file; to_mapper_str.append(basename2(m_filtered_alignment_file));

				m_ofs_filtered_alignment.open(to_mapper_str.c_str());

				string fusion_paired_str = alignment_file; fusion_paired_str.append(basename2(m_filtered_alignment_file)); fusion_paired_str.append(".fusion_paired");

				m_ofs_fusion_paired.open(fusion_paired_str.c_str());

				string single_str = alignment_file; single_str.append(basename2(m_filtered_alignment_file)); single_str.append(".single");

				m_ofs_single.open(single_str.c_str());

				string to_fusion_str = alignment_file; to_fusion_str.append(basename2(m_filtered_alignment_file)); to_fusion_str.append(".fusion");

				m_ofs_fusion_std.open(to_fusion_str.c_str());

				string paired_str = alignment_file; paired_str.append(basename2(m_filtered_alignment_file)); paired_str.append(".paired");

				m_ofs_paired.open(paired_str.c_str());

				string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

				m_ifs_alignments.open(onespliced.c_str());
			}
			else
				return;
		}
		else
		{
			m_ofs_filtered_alignment.close();

			m_ofs_fusion_std.close();

			m_ofs_fusion_paired.close();

			m_ofs_single.close();

			m_ofs_paired.close();

			string to_mapper_str = alignment_file; to_mapper_str.append(basename2(m_filtered_alignment_file));

			m_ofs_filtered_alignment.open(to_mapper_str.c_str());

			string fusion_paired_str = alignment_file; fusion_paired_str.append(basename2(m_filtered_alignment_file)); fusion_paired_str.append(".fusion_paired");

			m_ofs_fusion_paired.open(fusion_paired_str.c_str());

			string single_str = alignment_file; single_str.append(basename2(m_filtered_alignment_file)); single_str.append(".single");

			m_ofs_single.open(single_str.c_str());

			string to_fusion_str = alignment_file; to_fusion_str.append(basename2(m_filtered_alignment_file)); to_fusion_str.append(".fusion");

			m_ofs_fusion_std.open(to_fusion_str.c_str());

			string paired_str = alignment_file; paired_str.append(basename2(m_filtered_alignment_file)); paired_str.append(".paired");

			m_ofs_paired.open(paired_str.c_str());

			string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

			m_ifs_alignments.open(onespliced.c_str());
		}
	}
	

	//parallel

	if (m_is_paired)
	{
		//#ifdef LINUX
		//vector<pthread_t> iThreadIds(threads_number);

		//vector<int> return_values(threads_number);
		//#endif

		//while (!m_ifs_alignments.eof())
		//{
		//	if (index == buf_size - 1)
		//	{
		//		cout << "parse"<<endl;

		//		for (size_t i = 0; i < threads_number; ++i)
		//		{
		//			#ifdef LINUX
		//			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
		//			#endif
		//		}

		//		for (size_t i = 0; i < threads_number; ++i)
		//		{
		//			#ifdef LINUX
		//			if (return_values[i] < 0)
		//				cout << "error creating thread:"<< errno<<endl;
		//			#endif
		//		}

		//		cout << "wait parse"<<endl;

		//		for (size_t i = 0; i < threads_number; ++i)
		//		{
		//			#ifdef LINUX
		//			pthread_join (iThreadIds[i], NULL);
		//			#endif
		//		}

		//		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		//		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		//		string prev_tag = sam_iter->tag_base_name;

		//		cout << "group same read alignments"<<endl;

		//		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		//		{
		//			if (prev_tag != sam_iter->tag_base_name)
		//			{
		//				//cout << prev_tag << endl;

		//				//cout << previous_iter - m_stored_sam_recs.begin() + 1<< '\t' << sam_iter - 1 - m_stored_sam_recs.begin() + 1 << endl;

		//				//cout << (sam_iter - 1)->tostring(0, 0) << endl;

		//				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

		//				prev_tag = sam_iter->tag_base_name;

		//				previous_iter = sam_iter;
		//			}
		//		}

		//		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		//		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		//		cout << "process same read alignments"<<endl;

		//		for (size_t i = 0; i < threads_number; ++i)
		//		{
		//			#ifdef LINUX
		//			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadMulti, NULL);
		//			#endif
		//		}

		//		for (size_t i = 0; i < threads_number; ++i)
		//		{
		//			#ifdef LINUX
		//			if (return_values[i] < 0)
		//				cout << "error creating thread:"<< errno<<endl;
		//			#endif
		//		}

		//		cout << "wait process same read alignments"<<endl;

		//		for (size_t i = 0; i < threads_number; ++i)
		//		{
		//			#ifdef LINUX
		//			pthread_join (iThreadIds[i], NULL);
		//			#endif
		//		}

		//		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		//		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
		//			 m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
		//			 ++m_same_reads_alignments_paired_iter)
		//		{
		//			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

		//			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

		//			CollectStats(m_same_reads_alignments_paired_iter->first);

		//			CollectStats(m_same_reads_alignments_paired_iter->second);
		//		}

		//		
		//		//for (sam_iter = m_stored_sam_recs.begin(); sam_iter != m_stored_sam_recs.end(); ++sam_iter)
		//		//{
		//		//	junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.first);
		//		//}

		//		//junc_handler->SamRecVec2Junc(m_sam_rec_pe_ptr.second);

		//		m_cur_read_alignments_index = 0;

		//		m_cur_line_index = 0;

		//		m_same_reads_alignments.clear();

		//		m_same_reads_alignments_paired.clear();

		//		index = 0;

		//		cout << "finish"<<endl;
		//		//keep the last read's alignments to be processed with next block of alignments

		//		//m_stored_sam_recs.clear();

		//		//m_stored_lines.clear();
		//	}

		//	getline(m_ifs_alignments, m_stored_lines[index++]);
		//}
		//loop
		//read next tag
		//size_t read_id;

		//size_t sam_count = m_sam_file_handler.ReadNextTagAlignPE(m_sam_rec_pe, read_id);

		//while (true)
		//{
		//	// if it is end of file then stop
		//	if (sam_count == 0)
		//		break;

		//	//process each read
		if (fileter_by_junc == 2)
			ProcessPEReadSelectBest(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc ==	1)
			ProcessPEReadFilterJunc(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc == 0)
			ProcessPERead(junc_handler, junc_handler_filtered);

		//	//read next tag
		//	sam_count = m_sam_file_handler.ReadNextTagAlignPE(m_sam_rec_pe, read_id);
		//}
	}
	else//SE reads
	{
				//process each read
			/*if (fileter_by_junc == 2)
				ProcessPEReadSelectBest(junc_handler, junc_handler_filtered);
			else */
		if (fileter_by_junc == 2)
			ProcessSEReadSelectBest(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc ==	1)
			ProcessSEReadFilterJunc(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc == 0)
			ProcessSERead(junc_handler, junc_handler_filtered);
	}

	m_ifs_alignments.close();
	//else//SE reads
	//{
	//	//loop
	//	//read next tag
	//	size_t read_id;

	//	size_t sam_count = m_sam_file_handler.ReadNextTagAlign(m_sam_rec_pe.first, read_id);

	//	while (true)
	//	{
	//		// if it is end of file then stop
	//		if (sam_count == 0)
	//			break;

	//		//process each read
	//		/*if (fileter_by_junc == 2)
	//			ProcessPEReadSelectBest(junc_handler, junc_handler_filtered);
	//		else */
	//		if (fileter_by_junc ==	1)
	//			ProcessSEReadFilterJunc(junc_handler, junc_handler_filtered);
	//		else if (fileter_by_junc == 0)
	//			ProcessSERead(junc_handler, junc_handler_filtered);

	//		//read next tag
	//		sam_count = m_sam_file_handler.ReadNextTagAlign(m_sam_rec_pe.first, read_id);
	//	}
	//}
}


void 
AlignmentHandler::FilterAlignment()
{
	if (!m_junction_file.empty())
	{
		vector<string> annotated_junction;

		annotated_junction.push_back(m_junction_file);

		cout << "read annotated junctions"<< endl;

		m_junction_handler_annotated.ReadJunction(annotated_junction);
	}

	m_junction_handler.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch,  m_min_junc_anchor, m_min_fusion_coverage);

	m_junction_handler_filtered.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor, m_min_fusion_coverage);

	m_junction_handler_intermediate.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor, m_min_fusion_coverage);

	cout << "Allocate memory for Hits" << endl;
	
	m_junction_handler.ReadChromSize(m_chrom_size_file);
	//AllocateHitsMemory();

	JunctionHandler* junction_handler_ptr = &m_junction_handler;

	JunctionHandler* junction_handler_filtered_ptr = &m_junction_handler_filtered;

	JunctionHandler* junction_handler_intermediate_ptr = &m_junction_handler_intermediate;

	string stat_file = m_filtered_alignment_file; stat_file.append(".stat");

	vector<string>::iterator sam_file_iter;

	clock_t t1, t2;

	time_t m1, m2;

	double m_any_time;

	t1=clock();

	m1 = time(NULL);

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_paired.close();

	m_ofs_unspliced_fusion_paired.close();

	m_ofs_unspliced_single.close();

	m_ofs_unspliced_paired.close();

	m_ofs_onespliced.close();

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_ptr, 0);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	////cout << "Convert hits to expressed regions " << endl;

	////Hits2ExpressedRegions();

	//cout << "load flank string 1

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first load flank string: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	WriteStats(stat_file, "original alignment stats");

	ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "original junction stats");

	string ori_junction_file = m_filtered_alignment_file; ori_junction_file.append(".ori.junc");

	string ori_junction_ins_file = m_filtered_alignment_file; ori_junction_ins_file.append(".ori.junc.ins");

	string ori_junction_del_file = m_filtered_alignment_file; ori_junction_del_file.append(".ori.junc.del");

	string ori_junction_fusion_file = m_filtered_alignment_file; ori_junction_fusion_file.append(".ori.junc.fusion");

	string ori_junction_filter_file = m_filtered_alignment_file; ori_junction_filter_file.append(".ori.junc.filter");

	m_junction_handler.WriteJunction(ori_junction_file, ori_junction_ins_file, ori_junction_del_file, ori_junction_fusion_file, ori_junction_filter_file);

	cout << "mark fitered junction "<< endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best
	//16 
	//32 if set, open one output file for all input sam files 
	//   if not, open an output file for each input sam file
	//64 if set, output all alignments to a file for 2nd and 3rd steps process
	//   if not set, output both unspliced to a file, the others to another file for 2nd and 3rd step process
	//128 if 64 is not set
	         // if 128 is not set, output all both unspliced to a file
	         // if 128 is set, output both unspliced and normal paired to a file, and single and fusion paired to another file for 2nd and 3rd step process	  
	//256 if set, not filter single reads if 128 is not set, for normal filtering with fusion step and paired end reads
	//    if not set, filter single end reads

	if (m_do_filter & 4)
	{
		m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);

		m_junction_handler.ClearStats();

		m_junction_handler.CollectStats();

		m_junction_handler.WriteStats(stat_file, "original filtered junction stats");

		string ori_fil_junction_file = m_filtered_alignment_file; ori_fil_junction_file.append(".ori.fil.junc");

		string ori_fil_junction_ins_file = m_filtered_alignment_file; ori_fil_junction_ins_file.append(".ori.fil.junc.ins");

		string ori_fil_junction_del_file = m_filtered_alignment_file; ori_fil_junction_del_file.append(".ori.fil.junc.del");

		string ori_fil_junction_fusion_file = m_filtered_alignment_file; ori_fil_junction_fusion_file.append(".ori.fil.junc.fusion");

		string ori_fil_junction_filter_file = m_filtered_alignment_file; ori_fil_junction_filter_file.append(".ori.fil.junc.filter");

		m_junction_handler.WriteJunction(ori_fil_junction_file, ori_fil_junction_ins_file, ori_fil_junction_del_file, ori_fil_junction_fusion_file, ori_fil_junction_filter_file);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "filtering: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	/***********add intermediate filtering here*****************/

	t1=clock();

	m1 = time(NULL);

	if (m_do_filter & 8)//do select best
		m_do_filter = m_do_filter | 16;

	cerr << "Second filtering"<< endl;

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_paired.close();

	m_ofs_unspliced_fusion_paired.close();

	m_ofs_unspliced_single.close();

	m_ofs_unspliced_paired.close();

	m_ofs_onespliced.close();

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_ptr, 2);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	WriteStats(stat_file, "intermediate alignment stats");

	ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "intermediate junction stats");

	string int_junction_file = m_filtered_alignment_file; int_junction_file.append(".inter.junc");

	string int_junction_ins_file = m_filtered_alignment_file; int_junction_ins_file.append(".inter.junc.ins");

	string int_junction_del_file = m_filtered_alignment_file; int_junction_del_file.append(".inter.junc.del");

	string int_junction_fusion_file = m_filtered_alignment_file; int_junction_fusion_file.append(".inter.junc.fusion");

	string int_junction_filter_file = m_filtered_alignment_file; int_junction_filter_file.append(".inter.junc.filter");

	m_junction_handler.WriteJunction(int_junction_file, int_junction_ins_file, int_junction_del_file, int_junction_fusion_file, int_junction_filter_file);

	cout << "mark intermediate fitered junction "<< endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best
	//16 
	//32 if set, open one output file for all input sam files 
	//   if not, open an output file for each input sam file
	//64 if set, output all alignments to a file for 2nd and 3rd steps process
	//   if not set, output both unspliced to a file, the others to another file for 2nd and 3rd step process
	//128 if 64 is not set
	         // if 128 is not set, output all both unspliced to a file
	         // if 128 is set, output both unspliced and normal paired to a file, and single and fusion paired to another file for 2nd and 3rd step process	

	if (m_do_filter & 4)
	{
		m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);

		m_junction_handler.ClearStats();

		m_junction_handler.CollectStats();

		m_junction_handler.WriteStats(stat_file, "intermediate filtered junction stats");

		string int_fil_junction_file = m_filtered_alignment_file; int_fil_junction_file.append(".inter.fil.junc");

		string int_fil_junction_ins_file = m_filtered_alignment_file; int_fil_junction_ins_file.append(".inter.fil.junc.ins");

		string int_fil_junction_del_file = m_filtered_alignment_file; int_fil_junction_del_file.append(".inter.fil.junc.del");

		string int_fil_junction_fusion_file = m_filtered_alignment_file; int_fil_junction_fusion_file.append(".inter.fil.junc.fusion");

		string int_fil_junction_filter_file = m_filtered_alignment_file; int_fil_junction_filter_file.append(".inter.fil.junc.filter");

		m_junction_handler.WriteJunction(int_fil_junction_file, int_fil_junction_ins_file, int_fil_junction_del_file, int_fil_junction_fusion_file, int_fil_junction_filter_file);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "filtering: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;


	/****************************/

	//m_junction_handler.Clear();

	//cout << "Union expression regions " << endl;

	t1=clock();

	m1 = time(NULL);

	//UnionSets(&m_junction_handler, &m_expressed_regions);

	m_junction_handler.LoadFusionJuncToSortVec();

	m_junction_handler.SortFusionJunc();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process prepare: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best
	//16 
	//32 if set, open one output file for all input sam files 
	//   if not, open an output file for each input sam file
	//64 if set, output all alignments to a file for 2nd and 3rd steps process
	//   if not set, output both unspliced to a file, the others to another file for 2nd and 3rd step process
	//128 if 64 is not set
	         // if 128 is not set, output all both unspliced to a file
	         // if 128 is set, output both unspliced and normal paired to a file, and single and fusion paired to another file for 2nd and 3rd step process	

	t1=clock();

	m1 = time(NULL);

	if (m_do_filter & 8)//do select best
		m_do_filter = m_do_filter | 16;

	cerr << "third filtering"<< endl;

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_paired.close();

	m_ofs_unspliced_fusion_paired.close();

	m_ofs_unspliced_single.close();

	m_ofs_unspliced_paired.close();

	m_ofs_onespliced.close();

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_ptr, 1);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	WriteStats(stat_file, "filtered alignment stats");

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "filtered junction");

	string fil_junction_file = m_filtered_alignment_file; fil_junction_file.append(".fil.junc");

	string fil_junction_ins_file = m_filtered_alignment_file; fil_junction_ins_file.append(".fil.junc.ins");

	string fil_junction_del_file = m_filtered_alignment_file; fil_junction_del_file.append(".fil.junc.del");

	string fil_junction_fusion_file = m_filtered_alignment_file; fil_junction_fusion_file.append(".fil.junc.fusion");

	string fil_junction_filter_file = m_filtered_alignment_file; fil_junction_filter_file.append(".fil.junc.filter");

	m_junction_handler.WriteJunction(fil_junction_file, fil_junction_ins_file, fil_junction_del_file, fil_junction_fusion_file, fil_junction_filter_file);
	
	string ori_junction_fusion_encompass_file = m_filtered_alignment_file; ori_junction_fusion_encompass_file.append(".ori.junc.fusion.encompass");

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third write junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	//fusion encompassing 

	cout << "gene fuion struct" << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.GenerateFusionStruct();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "generate structure: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	cout << "write fusion struct" << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.WriteFusionJunctionWithEncompassingAlignments(ori_junction_fusion_encompass_file, junction_handler_ptr);

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "write fusion junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;
}

void
AlignmentHandler::FilterAlignmentFiltered()
{
	if (!m_junction_file.empty())
	{
		vector<string> annotated_junction;

		annotated_junction.push_back(m_junction_file);

		cout << "read annotated junctions"<< endl;

		m_junction_handler_annotated.ReadJunction(annotated_junction);
	}

	m_junction_handler.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch,  m_min_junc_anchor, m_min_fusion_coverage);

	m_junction_handler_filtered.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor, m_min_fusion_coverage);

	//m_junction_handler_intermediate.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor);

	cout << "Allocate memory for Hits" << endl;

	//AllocateHitsMemory();

	JunctionHandler* junction_handler_ptr = &m_junction_handler;

	JunctionHandler* junction_handler_filtered_ptr = &m_junction_handler_filtered;

	JunctionHandler* junction_handler_intermediate_ptr = &m_junction_handler_intermediate;

	string stat_file = m_filtered_alignment_file; stat_file.append(".stat");

	vector<string>::iterator sam_file_iter;

	clock_t t1, t2;

	time_t m1, m2;

	double m_any_time;

	t1=clock();

	m1 = time(NULL);

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_filtered_ptr, 0);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	////cout << "Convert hits to expressed regions " << endl;

	////Hits2ExpressedRegions();

	//cout << "load flank string 1

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first load flank string: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	WriteStats(stat_file, "original alignment stats");

	ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "original junction stats");

	string ori_junction_file = m_filtered_alignment_file; ori_junction_file.append(".ori.junc");

	string ori_junction_ins_file = m_filtered_alignment_file; ori_junction_ins_file.append(".ori.junc.ins");

	string ori_junction_del_file = m_filtered_alignment_file; ori_junction_del_file.append(".ori.junc.del");

	string ori_junction_fusion_file = m_filtered_alignment_file; ori_junction_fusion_file.append(".ori.junc.fusion");

	string ori_junction_filter_file = m_filtered_alignment_file; ori_junction_filter_file.append(".ori.junc.filter");

	m_junction_handler.WriteJunction(ori_junction_file, ori_junction_ins_file, ori_junction_del_file, ori_junction_fusion_file, ori_junction_filter_file);

	cout << "mark fitered junction "<< endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best

	if (m_do_filter & 4)
	{
		m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);

		m_junction_handler.ClearStats();

		m_junction_handler.CollectStats();

		m_junction_handler.WriteStats(stat_file, "original filtered junction stats");

		string ori_fil_junction_file = m_filtered_alignment_file; ori_fil_junction_file.append(".ori.fil.junc");

		string ori_fil_junction_ins_file = m_filtered_alignment_file; ori_fil_junction_ins_file.append(".ori.fil.junc.ins");

		string ori_fil_junction_del_file = m_filtered_alignment_file; ori_fil_junction_del_file.append(".ori.fil.junc.del");

		string ori_fil_junction_fusion_file = m_filtered_alignment_file; ori_fil_junction_fusion_file.append(".ori.fil.junc.fusion");

		string ori_fil_junction_filter_file = m_filtered_alignment_file; ori_fil_junction_filter_file.append(".ori.fil.junc.filter");

		m_junction_handler.WriteJunction(ori_fil_junction_file, ori_fil_junction_ins_file, ori_fil_junction_del_file, ori_fil_junction_fusion_file, ori_fil_junction_filter_file);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "filtering: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	//if (m_do_filter & 8)//do select best
	//	m_do_filter = m_do_filter | 16;

	/***********add intermediate filtering here*****************/

	//t1=clock();

	//m1 = time(NULL);

	//cerr << "Second filtering"<< endl;

	//for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	//{
	//	ProcessFile(*sam_file_iter, junction_handler_intermediate_ptr, junction_handler_ptr, 2);
	//}

	//t2=clock();

	//m2 = time(NULL);

	//m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	//cerr << "second process: " << m_any_time << endl;

	//cerr << "seconds passed " << m2 - m1 << endl;

	//t1=clock();

	//m1 = time(NULL);

	//m_junction_handler_intermediate.LoadFlankString();

	//t2=clock();

	//m2 = time(NULL);

	//m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	//cerr << "second process: " << m_any_time << endl;

	//cerr << "seconds passed " << m2 - m1 << endl;



	//t1=clock();

	//m1 = time(NULL);

	//WriteStats(stat_file, "intermediate alignment stats");

	//ClearStats();

	//m_junction_handler_intermediate.CollectStats();

	//m_junction_handler_intermediate.WriteStats(stat_file, "intermediate junction stats");

	//string int_junction_file = m_filtered_alignment_file; int_junction_file.append(".inter.junc");

	//string int_junction_ins_file = m_filtered_alignment_file; int_junction_ins_file.append(".inter.junc.ins");

	//string int_junction_del_file = m_filtered_alignment_file; int_junction_del_file.append(".inter.junc.del");

	//string int_junction_fusion_file = m_filtered_alignment_file; int_junction_fusion_file.append(".inter.junc.fusion");

	//string int_junction_filter_file = m_filtered_alignment_file; int_junction_filter_file.append(".inter.junc.filter");

	//m_junction_handler_intermediate.WriteJunction(int_junction_file, int_junction_ins_file, int_junction_del_file, int_junction_fusion_file, int_junction_filter_file);

	//cout << "mark intermediate fitered junction "<< endl;

	////do_filter
	////1 pair
	////2 pairing filter
	////4 filter junc
	////8 select best

	//if (m_do_filter & 4)
	//{
	//	m_junction_handler_intermediate.MarkFiltered(m_is_paired, m_junction_handler_annotated);

	//	m_junction_handler_intermediate.ClearStats();

	//	m_junction_handler_intermediate.CollectStats();

	//	m_junction_handler_intermediate.WriteStats(stat_file, "intermediate filtered junction stats");

	//	string int_fil_junction_file = m_filtered_alignment_file; int_fil_junction_file.append(".inter.fil.junc");

	//	string int_fil_junction_ins_file = m_filtered_alignment_file; int_fil_junction_ins_file.append(".inter.fil.junc.ins");

	//	string int_fil_junction_del_file = m_filtered_alignment_file; int_fil_junction_del_file.append(".inter.fil.junc.del");

	//	string int_fil_junction_fusion_file = m_filtered_alignment_file; int_fil_junction_fusion_file.append(".inter.fil.junc.fusion");

	//	string int_fil_junction_filter_file = m_filtered_alignment_file; int_fil_junction_filter_file.append(".inter.fil.junc.filter");

	//	m_junction_handler_intermediate.WriteJunction(int_fil_junction_file, int_fil_junction_ins_file, int_fil_junction_del_file, int_fil_junction_fusion_file, int_fil_junction_filter_file);
	//}

	//t2=clock();

	//m2 = time(NULL);

	//m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	//cerr << "filtering: " << m_any_time << endl;

	//cerr << "seconds passed " << m2 - m1 << endl;


	/****************************/

	//m_junction_handler.Clear();

	//cout << "Union expression regions " << endl;

	t1=clock();

	m1 = time(NULL);

	//UnionSets(&m_junction_handler, &m_expressed_regions);

	m_junction_handler.LoadFusionJuncToSortVec();

	m_junction_handler.SortFusionJunc();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process prepare: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	////do_filter
	////1 pair
	////2 pairing filter
	////4 filter junc
	////8 select best

	t1=clock();

	m1 = time(NULL);

	if (m_do_filter & 8)//do select best
		m_do_filter = m_do_filter | 16;

	cerr << "third filtering"<< endl;

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_filtered_ptr, junction_handler_ptr, 1);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler_filtered.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	WriteStats(stat_file, "filtered alignment stats");

	t1=clock();

	m1 = time(NULL);

	m_junction_handler_filtered.CollectStats();

	m_junction_handler_filtered.WriteStats(stat_file, "filtered junction");

	string fil_junction_file = m_filtered_alignment_file; fil_junction_file.append(".fil.junc");

	string fil_junction_ins_file = m_filtered_alignment_file; fil_junction_ins_file.append(".fil.junc.ins");

	string fil_junction_del_file = m_filtered_alignment_file; fil_junction_del_file.append(".fil.junc.del");

	string fil_junction_fusion_file = m_filtered_alignment_file; fil_junction_fusion_file.append(".fil.junc.fusion");

	string fil_junction_filter_file = m_filtered_alignment_file; fil_junction_filter_file.append(".fil.junc.filter");

	m_junction_handler_filtered.WriteJunction(fil_junction_file, fil_junction_ins_file, fil_junction_del_file, fil_junction_fusion_file, fil_junction_filter_file);
	
	string ori_junction_fusion_encompass_file = m_filtered_alignment_file; ori_junction_fusion_encompass_file.append(".ori.junc.fusion.encompass");

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third write junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	//fusion encompassing 

	cout << "gene fuion struct" << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.GenerateFusionStruct();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "generate structure: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	cout << "write fusion struct" << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.WriteFusionJunctionWithEncompassingAlignments(ori_junction_fusion_encompass_file, junction_handler_filtered_ptr);

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "write fusion junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;
}

void 
AlignmentHandler::WriteAlignment(vector<SamRec*>& sam_rec_ptr)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion)
			m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;
		else
			m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;

		if ((*sam_rec_ptr_iter)->is_fusion)
			m_ofs_fusion_std << (*sam_rec_ptr_iter)->output_line<<endl;

		if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
			m_ofs_single << (*sam_rec_ptr_iter)->output_line<<endl; 
		else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
			m_ofs_fusion_paired << (*sam_rec_ptr_iter)->output_line<<endl;
		else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
		{
			m_ofs_paired << (*sam_rec_ptr_iter)->output_line<<endl;
			//cerr << "paired"<<endl;
		}
		else
			cerr << (*sam_rec_ptr_iter)->paired_type << endl;
	}

}

void 
AlignmentHandler::WriteAlignmentUnspliced(vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator pair_samvec_iter, bool is_unspliced/*;vector<SamRec*>& sam_rec_ptr*/)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	//bool is_unspliced = true;

	//for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end() && is_unspliced; ++sam_rec_ptr_iter)
	//{
	//	if ((*sam_rec_ptr_iter)->isexonic == false && (*sam_rec_ptr_iter)->is_fusion == false)
	//		is_unspliced = false;
	//}

	//for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end() && is_unspliced; ++sam_rec_ptr_iter)
	//{
	//	if ((*sam_rec_ptr_iter)->isexonic == false && (*sam_rec_ptr_iter)->is_fusion == false)
	//		is_unspliced = false;
	//}

	if (is_unspliced && ((m_do_filter & 64) == 0))
	{
		if ((m_do_filter & 128) == 0)
		{
			for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end() && is_unspliced; ++sam_rec_ptr_iter)
			{
				//	if ((*sam_rec_ptr_iter)->is_fusion)
				//	m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;
				//else
				//	m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;

				//if ((*sam_rec_ptr_iter)->is_fusion)
				//	m_ofs_fusion_std << (*sam_rec_ptr_iter)->output_line<<endl;

				if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_unspliced_single << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_unspliced_fusion_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
					//cerr << "paired"<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
				//if 
				//m_ofs_bothunspliced << (*sam_rec_ptr_iter)->tostring(0, 0) << endl;
			}

			for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end() && is_unspliced; ++sam_rec_ptr_iter)
			{
				if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_unspliced_single << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_unspliced_fusion_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
					//cerr << "paired"<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
				//m_ofs_bothunspliced << (*sam_rec_ptr_iter)->tostring(0, 0) << endl;
			}
		}
		else
		{
			for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end() && is_unspliced; ++sam_rec_ptr_iter)
			{
				//	if ((*sam_rec_ptr_iter)->is_fusion)
				//	m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;
				//else
				//	m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;

				//if ((*sam_rec_ptr_iter)->is_fusion)
				//	m_ofs_fusion_std << (*sam_rec_ptr_iter)->output_line<<endl;

				if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
					//cerr << "paired"<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
				//if 
				//m_ofs_bothunspliced << (*sam_rec_ptr_iter)->tostring(0, 0) << endl;
			}

			for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end() && is_unspliced; ++sam_rec_ptr_iter)
			{
				if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
					//cerr << "paired"<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
				//m_ofs_bothunspliced << (*sam_rec_ptr_iter)->tostring(0, 0) << endl;
			}
		}
	}
	else
	{
		for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end(); ++sam_rec_ptr_iter)
		{
			m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line << endl;
		}

		for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end(); ++sam_rec_ptr_iter)
		{
			m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line << endl;
		}
	}
}

void
AlignmentHandler::SetBitInfo(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr)
{
	for (size_t i = 0; i < sam_rec_pe_ptr.first.size(); ++i)
	{
		sam_rec_pe_ptr.first[i]->strand_t |= IS_PAIRED;

		if (!sam_rec_pe_ptr.first[i]->is_fusion)
		{
			sam_rec_pe_ptr.first[i]->strand_t |= IS_FIRST_END;
		}
		else
			sam_rec_pe_ptr.first[i]->strand_t2 |= IS_PAIRED;
	}

	for (size_t i = 0; i < sam_rec_pe_ptr.second.size(); ++i)
	{
		sam_rec_pe_ptr.second[i]->strand_t |= IS_PAIRED;

		if (!sam_rec_pe_ptr.second[i]->is_fusion)
		{
			sam_rec_pe_ptr.second[i]->strand_t |= IS_SECOND_END;
		}
		else
			sam_rec_pe_ptr.second[i]->strand_t2 |= IS_PAIRED;
	}
}

void* AlignmentHandler::SetBitInfoStatic(void * str)
{
	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *((pair<vector<SamRec*>, vector<SamRec*> >*)str);

	for (size_t i = 0; i < sam_rec_pe_ptr.first.size(); ++i)
	{
		sam_rec_pe_ptr.first[i]->strand_t |= IS_PAIRED;

		if (!sam_rec_pe_ptr.first[i]->is_fusion)
		{
			sam_rec_pe_ptr.first[i]->strand_t |= IS_FIRST_END;
		}
		else
			sam_rec_pe_ptr.first[i]->strand_t2 |= IS_PAIRED;
	}

	for (size_t i = 0; i < sam_rec_pe_ptr.second.size(); ++i)
	{
		sam_rec_pe_ptr.second[i]->strand_t |= IS_PAIRED;

		if (!sam_rec_pe_ptr.second[i]->is_fusion)
		{
			sam_rec_pe_ptr.second[i]->strand_t |= IS_SECOND_END;
		}
		else
			sam_rec_pe_ptr.second[i]->strand_t2 |= IS_PAIRED;
	}
}

void
AlignmentHandler::SetFusionBitInfo(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr)
{
	for (size_t i = 0; i < sam_rec_pe_ptr.first.size(); ++i)
	{
		
		if (sam_rec_pe_ptr.first[i]->is_fusion)
		{
			sam_rec_pe_ptr.first[i]->setfusionbit();
		}
	}

	for (size_t i = 0; i < sam_rec_pe_ptr.second.size(); ++i)
	{
		//sam_rec_pe_ptr.second[i]->strand_t |= IS_PAIRED;

		if (sam_rec_pe_ptr.second[i]->is_fusion)
		{
			sam_rec_pe_ptr.second[i]->setfusionbit();
		}
	}
}

void* AlignmentHandler::SetFusionBitInfoStatic(void * str)
{
	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *((pair<vector<SamRec*>, vector<SamRec*> >*)str);

	for (size_t i = 0; i < sam_rec_pe_ptr.first.size(); ++i)
	{
		if (sam_rec_pe_ptr.first[i]->is_fusion)
		{
			sam_rec_pe_ptr.first[i]->setfusionbit();
		}
	}

	for (size_t i = 0; i < sam_rec_pe_ptr.second.size(); ++i)
	{
		if (sam_rec_pe_ptr.second[i]->is_fusion)
		{
			sam_rec_pe_ptr.second[i]->setfusionbit();
		}
	}

	return 0;
}

void
AlignmentHandler::WriteAlignment(vector<SamRec*>& sam_rec_ptr, ofstream& cur_ofs)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		
		if ((*sam_rec_ptr_iter)->is_fusion)
			cur_ofs << (*sam_rec_ptr_iter)->tostandfusion()<<endl;
		else
			cur_ofs << (*sam_rec_ptr_iter)->tostring(sam_rec_ptr.size(), sam_rec_ptr_iter - sam_rec_ptr.begin() + 1)<<endl;
	}
}

void AlignmentHandler::SetPairedType(vector<SamRec*>& sam_rec_ptr, PAIRED_TYPE paired_type)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion)
			(*sam_rec_ptr_iter)->paired_type = paired_type;
		else
			(*sam_rec_ptr_iter)->paired_type = paired_type;
	}
}

void* AlignmentHandler::GenerateAlignment(void * str, bool fusion_std)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>* )str);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	if (sam_rec_ptr.empty())
		return 0;

	size_t is_primary_count = 0;

	//cout << "GenerateAlignment 1" << endl;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		//cout << "g 0" << endl;
		if ((*sam_rec_ptr_iter)->is_fusion)
		{
			//cout << "g 1" << endl;
			if (fusion_std)
			{
				//cout << "g 2" << endl;
				(*sam_rec_ptr_iter)->output_line = (*sam_rec_ptr_iter)->tostandfusion();
			}
			else
			{
				//cout << "g 3" << endl;
				(*sam_rec_ptr_iter)->output_line = (*sam_rec_ptr_iter)->tostring(0, 0);
			}
		}
		else
		{
			//cout << "g 4" << endl;

			(*sam_rec_ptr_iter)->output_line = (*sam_rec_ptr_iter)->tostring(sam_rec_ptr.size(), sam_rec_ptr_iter - sam_rec_ptr.begin() + 1);
		}

		//cout << "g 5" << endl;
		//cout << (*sam_rec_ptr_iter)->output_line << endl;

		if (((*sam_rec_ptr_iter)->strand_t & IS_PRIMARY) == 0)
			++is_primary_count;

		if ((*sam_rec_ptr_iter)->is_fusion)
		{
			if (((*sam_rec_ptr_iter)->strand_t2 & IS_PRIMARY) == 0)
				++is_primary_count;
		}
	}

	//cout << "GenerateAlignment 2" << endl;

	if (is_primary_count != 1 && fusion_std)
	{
		cout << "is_primary_count:"<< is_primary_count << endl;

		cout <<sam_rec_ptr.front()->tag_name<<endl;
	}

	return 0;
}

void
AlignmentHandler::FilterCanonNonCanonByReads(vector<SamRec>& read_sam, vector<SamRec* >& filtered_read_sam_ptr)
{
	int count=0, tag_count = 0, unspliced = 0, filtered_canon = 0, filtered_noncanon = 0, filtered_canon_noncanon = 0;

	/*MarkCanonNoncanonByReads(read_sam);*/

	vector<SamRec>::iterator samrec_iter;

	bool canon = false, noncanon = false, is_insert = false, not_insert = false;

	for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
	{
		if (samrec_iter->canon_count)
			canon =true;

		if (samrec_iter->noncanon_count)
			noncanon =true;

		if (samrec_iter->issmallins)
			is_insert = true;
		else
			not_insert = true;
	}

	int insert_model;

	if (is_insert && not_insert)
		insert_model = 0;
	else if (!is_insert && not_insert)
		insert_model = 1;
	else if (is_insert && !not_insert)
		insert_model = -1;
	else
	{
		insert_model = -2;

		cout << "not possible"<<endl;
	}

	if ((canon && !noncanon) || (!canon && noncanon) || (!canon && !noncanon))
	{
		if (insert_model == 0)
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (samrec_iter->issmallins == false)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
					//ofs_filtered_canon << samrec_iter->tostring()<<endl;
				}					
				else 
					/*ofs_fitlered_ins << samrec_iter->tostring()<<endl*/;
			}
		}
		else
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				filtered_read_sam_ptr.push_back(&(*samrec_iter));
				//ofs_filtered_canon << samrec_iter->tostring()<<endl;
			}
		}

		++filtered_canon;
	}
	else if (canon && noncanon)
	{
		double max_canon_rate = 0;

		bool noncan_can = false, noncan = false, unspliced = false;

		if (insert_model == 0)
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (max_canon_rate < samrec_iter->canon_rate && samrec_iter->issmallins == false)
					max_canon_rate = samrec_iter->canon_rate;
			}

			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (samrec_iter->issmallins)
				{
					/*ofs_fitlered_ins << samrec_iter->tostring()<<endl*/;
				}
				else if (max_canon_rate == samrec_iter->canon_rate)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
					/*ofs_filtered_canon << samrec_iter->tostring()<<endl*/;
					//++count;
				}
				else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
				{
					/*ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl*/;

					noncan_can = true;
				}
				else if (samrec_iter->noncanon_count)
				{
					/*ofs_filtered_noncanon << samrec_iter->tostring()<<endl*/;

					noncan = true;
				}
				else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
					/*unspliced_ofs << samrec_iter->tostring()<<endl*/;

					unspliced = true;
				}

			}
		}
		else
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (max_canon_rate < samrec_iter->canon_rate)
					max_canon_rate = samrec_iter->canon_rate;
			}

			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{

				if (max_canon_rate == samrec_iter->canon_rate)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
					//ofs_filtered_canon << samrec_iter->tostring()<<endl;
					//++count;
				}
				else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
				{
					/*ofs_fitlered_noncanon_canon << samrec_iter->tostring()<<endl*/;

					noncan_can = true;
				}
				else if (samrec_iter->noncanon_count)
				{
					/*ofs_filtered_noncanon << samrec_iter->tostring()<<endl*/;

					noncan = true;
				}
				else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
				{
					//cout << "should not be here one reads"<<endl;
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
					//unspliced_ofs << samrec_iter->tostring()<<endl;

					unspliced = true;
				}
			}
		}
	}
}



void
AlignmentHandler::FilterPairedSamRec(vector<PairedSamRec>& pairedsamrec_vec, pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr)
{
	if (pairedsamrec_vec.size() == 0)
		return;

	vector<PairedSamRec>::iterator psm_iter;

	sort(pairedsamrec_vec.begin(), pairedsamrec_vec.end(), comp_dist);

	vector<PairedSamRec>::iterator vps_iter;

	if (m_do_filter & 16)
	{
		sam_rec_pe_ptr.first.clear();

		sam_rec_pe_ptr.second.clear();
	}

	vps_iter = pairedsamrec_vec.begin();

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

		vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;//static_cast<long> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.first->mate_match = "=";
		else
			vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

		if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
	}

	if (vps_iter->paired_sam_rec.first->is_fusion)
	{
		vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

		vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

		vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

		vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
	}

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

		vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist; //static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.second->mate_match = "=";
		else
			vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

		if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

	}

	if (vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

		vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

		vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

		vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
	}

	//output

	if (m_do_filter & 16)
	{
		sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

		sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
	}

	for (vps_iter = pairedsamrec_vec.begin() + 1; vps_iter != pairedsamrec_vec.end(); ++vps_iter)
	{
		if ((vps_iter - 1)->paired_sam_rec.first->tag_name.find("seq.10001997/") != string::npos)
		{
			cout << (vps_iter - 1)->paired_sam_rec.first->tostring(1,1) << endl << (vps_iter - 1)->paired_sam_rec.second->tostring(1,1) << endl;

			cout << (vps_iter)->paired_sam_rec.first->tostring(1,1) << endl << (vps_iter)->paired_sam_rec.second->tostring(1,1) << endl;

			cout <<"mappedlen" << endl<< (vps_iter - 1)->mappedlen << endl <<  (vps_iter)->mappedlen << endl;

			cout <<"total_mismatch" << endl<< (vps_iter - 1)->total_mismatch << endl <<  (vps_iter)->total_mismatch << endl;

			cout <<"mate_dist" << endl<< (vps_iter - 1)->mate_dist << endl <<  (vps_iter)->mate_dist << endl;

			cout <<"total_ave_mismatch" << endl<< (vps_iter - 1)->total_ave_mismatch << endl <<  (vps_iter)->total_ave_mismatch << endl;

			cout <<"intron_size" << endl<< (vps_iter - 1)->intron_size << endl <<  (vps_iter)->intron_size << endl;

			cout <<"total_pairing_rate" << endl<< (vps_iter - 1)->total_pairing_rate << endl <<  (vps_iter)->total_pairing_rate << endl;

			cout <<"total_filter_score" << endl<< (vps_iter - 1)->total_filter_score << endl <<  (vps_iter)->total_filter_score << endl;

			cout <<"total_filter_score" << endl<< (vps_iter - 1)->total_hits << endl <<  (vps_iter)->total_hits << endl;

			cout <<"mappedlen" << endl<< (vps_iter - 1)->total_anchor_len << endl <<  (vps_iter)->total_anchor_len << endl;
			
		}

		if ((vps_iter - 1)->mappedlen == (vps_iter)->mappedlen && !comp_mate_dist(((vps_iter - 1)->mate_dist), ((vps_iter)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (vps_iter - 1)->total_mismatch == (vps_iter)->total_mismatch 
			//&& (vps_iter - 1)->total_pairing_rate == (vps_iter)->total_pairing_rate 
			//&& (vps_iter - 1)->total_anchor_len == (vps_iter)->total_anchor_len
			&& (vps_iter - 1)->total_ave_mismatch == (vps_iter)->total_ave_mismatch
			//&& (vps_iter - 1)->total_filter_score == (vps_iter)->total_filter_score
			&& !comp_intron_dist((vps_iter - 1)->intron_size, (vps_iter)->intron_size)
			&& (vps_iter - 1)->contiglen == (vps_iter)->contiglen
			&& (vps_iter - 1)->min_anchor == (vps_iter)->min_anchor
			//&& (vps_iter - 1)->total_pairing_rate == (vps_iter)->total_pairing_rate
			/*&& (vps_iter - 1)->intron_size == (vps_iter)->intron_size*/)
		{

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

				vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;//static_cast<int> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.first->mate_match = "=";
				else
					vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

				if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.first->is_fusion)
			{
				vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

				vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

				vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
			}

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{

				vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

				vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist;//static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.second->mate_match = "=";
				else
					vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

				if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

				vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

				vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
			}

			if (m_do_filter & 16)
			{
				sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

				sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
			}
		}
		else
			break;
	}
}

void* 
AlignmentHandler::FilterPairedSamRecStatic(void * str)
{
	FilterPairingStruct* filter_pairing = (FilterPairingStruct*)str;

	vector<PairedSamRec>& pairedsamrec_vec = *(filter_pairing->pairedsamrec_vec);
	
	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(filter_pairing->sam_rec_pe_ptr);

	int do_filter = filter_pairing->do_filter;

	//cout << "test pairedsamrec_vec:" << pairedsamrec_vec.size() << endl;
	//EstPairingStruct* est_pairing_ptr = (EstPairingStruct*) str;

	//vector<PairedSamRec> paired_reads_ptr;

	//pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(est_pairing_ptr->sam_rec_pe_ptr);

	//vector<PairedSamRec>& fusion_paired_reads_ptr = *(est_pairing_ptr->fusion_paired_reads_ptr);

	if (pairedsamrec_vec.size() == 0)
		return 0;

	vector<PairedSamRec>::iterator psm_iter;

	sort(pairedsamrec_vec.begin(), pairedsamrec_vec.end(), comp_dist_0);

	vector<PairedSamRec>::iterator same_vps_iter;

	for (same_vps_iter = pairedsamrec_vec.begin() + 1; same_vps_iter != pairedsamrec_vec.end(); ++same_vps_iter)
	{
		if ((same_vps_iter - 1)->mappedlen == (same_vps_iter)->mappedlen && !comp_mate_dist(((same_vps_iter - 1)->mate_dist), ((same_vps_iter)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (same_vps_iter - 1)->total_mismatch == (same_vps_iter)->total_mismatch 
			//&& (vps_iter - 1)->total_pairing_rate == (vps_iter)->total_pairing_rate 
			//&& (vps_iter - 1)->total_anchor_len == (vps_iter)->total_anchor_len
			/*&& (same_vps_iter - 1)->total_ave_mismatch == (same_vps_iter)->total_ave_mismatch*/)
		{
		}
		else
			break;
	}

	if (same_vps_iter - pairedsamrec_vec.begin() > 1)
		sort(pairedsamrec_vec.begin(), same_vps_iter, comp_dist_1);

	vector<PairedSamRec>::iterator same_vps_iter2;

	for (same_vps_iter2 = pairedsamrec_vec.begin() + 1; same_vps_iter2 != pairedsamrec_vec.end(); ++same_vps_iter2)
	{
		if ((same_vps_iter2 - 1)->mappedlen == (same_vps_iter2)->mappedlen && !comp_mate_dist(((same_vps_iter2 - 1)->mate_dist), ((same_vps_iter2)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& !comp_intron_dist(((same_vps_iter2 - 1)->intron_size), ((same_vps_iter2)->intron_size))
			&& (same_vps_iter2 - 1)->total_mismatch == (same_vps_iter2)->total_mismatch 
			//&& (vps_iter - 1)->total_pairing_rate == (vps_iter)->total_pairing_rate 
			//&& (vps_iter - 1)->total_anchor_len == (vps_iter)->total_anchor_len
			/*&& (same_vps_iter2 - 1)->total_ave_mismatch == (same_vps_iter2)->total_ave_mismatch*/)
		{
		}
		else
			break;
	}

	if (pairedsamrec_vec.begin()->intron_size && same_vps_iter2 - pairedsamrec_vec.begin() > 1)
		sort(pairedsamrec_vec.begin(), same_vps_iter2, comp_dist_2);

	vector<PairedSamRec>::iterator vps_iter;

	if (do_filter & 16)
	{
		sam_rec_pe_ptr.first.clear();

		sam_rec_pe_ptr.second.clear();
	}

	vps_iter = pairedsamrec_vec.begin();

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

		vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;//static_cast<long> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.first->mate_match = "=";
		else
			vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

		if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
	}

	if (vps_iter->paired_sam_rec.first->is_fusion)
	{
		vps_iter->paired_sam_rec.first->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

		vps_iter->paired_sam_rec.first->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

		vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

		//cout << "vps_iter->paired_sam_rec.first->fusion_prefix_len:"<< vps_iter->paired_sam_rec.first->fusion_prefix_len<<endl;

		vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

		vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

		//cout << "vps_iter->paired_sam_rec.first->fusion_suffix_len:"<< vps_iter->paired_sam_rec.first->fusion_suffix_len<<endl;

		vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
	}

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

		vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist; //static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.second->mate_match = "=";
		else
			vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

		if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

	}

	if (vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

		vps_iter->paired_sam_rec.second->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

		vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

		//cout << "vps_iter->paired_sam_rec.second->fusion_prefix_len:"<< vps_iter->paired_sam_rec.second->fusion_prefix_len<<endl;

		vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

		vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

		//cout << "vps_iter->paired_sam_rec.second->fusion_suffix_len:"<< vps_iter->paired_sam_rec.second->fusion_suffix_len<<endl;

		vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
	}

	//cout << vps_iter->paired_sam_rec.first->tag_name << endl;

	//cout << vps_iter->paired_sam_rec.second->tag_name << endl;

	//cout << vps_iter->paired_sam_rec.first->strand_t << endl;

	//cout << vps_iter->paired_sam_rec.second->strand_t << endl;


	if (vps_iter->paired_sam_rec.first->is_fusion)
	{
		if (vps_iter->paired_sam_rec.first->fusion_prefix_paired)
			vps_iter->paired_sam_rec.first->strand_t -= IS_PRIMARY;
		else if (vps_iter->paired_sam_rec.first->fusion_suffix_paired)
			vps_iter->paired_sam_rec.first->strand_t2 -= IS_PRIMARY;
	}
	else
		vps_iter->paired_sam_rec.first->strand_t -= IS_PRIMARY;

	if (vps_iter->paired_sam_rec.second->is_fusion)
	{
		if (vps_iter->paired_sam_rec.second->fusion_prefix_paired)
			vps_iter->paired_sam_rec.second->strand_t -= IS_PRIMARY;
		else if (vps_iter->paired_sam_rec.second->fusion_suffix_paired)
			vps_iter->paired_sam_rec.second->strand_t2 -= IS_PRIMARY;
	}
	else
		vps_iter->paired_sam_rec.second->strand_t -= IS_PRIMARY;

	//cout << vps_iter->paired_sam_rec.first->strand_t << endl;

	//cout << vps_iter->paired_sam_rec.second->strand_t << endl;
	
	//output

	if (do_filter & 16)
	{
		sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

		sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
	}

	for (vps_iter = pairedsamrec_vec.begin() + 1; vps_iter != pairedsamrec_vec.end(); ++vps_iter)
	{
		if ((vps_iter - 1)->paired_sam_rec.first->tag_name.find("R_15755/") != string::npos || (vps_iter - 1)->paired_sam_rec.first->tag_name.find("R_15755/") != string::npos ||
			(vps_iter - 1)->paired_sam_rec.first->tag_name.find("ILLUMINA-2F52BD_0038:2:114:18746:16774#0/") != string::npos || (vps_iter - 1)->paired_sam_rec.first->tag_name.find("seq.10012918/") != string::npos)
		{
			cout << (vps_iter - 1)->paired_sam_rec.first->tostring(1,1) << endl << (vps_iter - 1)->paired_sam_rec.second->tostring(1,1) << endl;

			cout << (vps_iter)->paired_sam_rec.first->tostring(1,1) << endl << (vps_iter)->paired_sam_rec.second->tostring(1,1) << endl;

			cout <<"mappedlen" << endl<< (vps_iter - 1)->mappedlen << endl <<  (vps_iter)->mappedlen << endl;

			cout <<"total_mismatch" << endl<< (vps_iter - 1)->total_mismatch << endl <<  (vps_iter)->total_mismatch << endl;

			cout <<"mate_dist" << endl<< (vps_iter - 1)->mate_dist << endl <<  (vps_iter)->mate_dist << endl;

			cout <<"total_ave_mismatch" << endl<< (vps_iter - 1)->total_ave_mismatch << endl <<  (vps_iter)->total_ave_mismatch << endl;

			cout <<"intron_size" << endl<< (vps_iter - 1)->intron_size << endl <<  (vps_iter)->intron_size << endl;

			cout <<"total_pairing_rate" << endl<< (vps_iter - 1)->total_pairing_rate << endl <<  (vps_iter)->total_pairing_rate << endl;

			cout <<"total_filter_score" << endl<< (vps_iter - 1)->total_filter_score << endl <<  (vps_iter)->total_filter_score << endl;

			cout <<"total_hits" << endl<< (vps_iter - 1)->total_hits << endl <<  (vps_iter)->total_hits << endl;

			cout <<"mappedlen" << endl<< (vps_iter - 1)->total_anchor_len << endl <<  (vps_iter)->total_anchor_len << endl;

			cout <<"min_anchor" << endl<< (vps_iter - 1)->min_anchor << endl <<  (vps_iter)->min_anchor << endl;

			cout <<"contiglen" << endl<< (vps_iter - 1)->contiglen << endl <<  (vps_iter)->contiglen << endl;
			
		}

		if ((vps_iter - 1)->mappedlen == (vps_iter)->mappedlen && !comp_mate_dist(((vps_iter - 1)->mate_dist), ((vps_iter)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (vps_iter - 1)->total_mismatch == (vps_iter)->total_mismatch 
			//&& (vps_iter - 1)->total_pairing_rate == (vps_iter)->total_pairing_rate 
			//&& (vps_iter - 1)->total_anchor_len == (vps_iter)->total_anchor_len
			//&& (vps_iter - 1)->total_ave_mismatch == (vps_iter)->total_ave_mismatch
			&& (vps_iter - 1)->intron_size == (vps_iter)->intron_size
			&& (vps_iter - 1)->total_hits == (vps_iter)->total_hits
			&& (vps_iter - 1)->total_filter_score == (vps_iter)->total_filter_score
			//&& !comp_intron_dist((vps_iter - 1)->intron_size, (vps_iter)->intron_size)
			&& (vps_iter - 1)->contiglen == (vps_iter)->contiglen
			&& (vps_iter - 1)->min_anchor == (vps_iter)->min_anchor
			//&& (vps_iter - 1)->total_pairing_rate == (vps_iter)->total_pairing_rate
			/*&& (vps_iter - 1)->intron_size == (vps_iter)->intron_size*/)
		{

			if (vps_iter->paired_sam_rec.first->mate_match != "*" || vps_iter->paired_sam_rec.second->mate_match != "*")
				continue;

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

				vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;//static_cast<int> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.first->mate_match = "=";
				else
					vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

				if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.first->is_fusion)
			{
				vps_iter->paired_sam_rec.first->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

				vps_iter->paired_sam_rec.first->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

				vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

				vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

				vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
			}

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{

				vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

				vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist;//static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.second->mate_match = "=";
				else
					vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

				if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.second->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

				vps_iter->paired_sam_rec.second->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

				vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

				vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

				vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
			}

			if (do_filter & 16)
			{
				sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

				sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
			}
		}
		else
			break;
	}

	return 0;
}

bool
AlignmentHandler::EstablishPairing(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr, vector<PairedSamRec>& fusion_paired_reads_ptr)
{
	//chrom pairid strand offset line
	//hash_map<string, hash_map<size_t, hash_map<unsigned short, map<size_t, vector<SamRec* > > > > > mapped_reads;

	//map<int, vector<pair<char, string> > > paired_reads;

	vector<PairedSamRec> paired_reads_ptr;

	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	//bool end1_mapped = false, end2_mapped = false;

	vector<SamRec*>::iterator sam_rec_ptr_iter1, sam_rec_ptr_iter2;

	for (sam_rec_ptr_iter1 = sam_rec_pe_ptr.first.begin(); sam_rec_ptr_iter1 != sam_rec_pe_ptr.first.end(); ++sam_rec_ptr_iter1)
	{
		for (sam_rec_ptr_iter2 = sam_rec_pe_ptr.second.begin(); sam_rec_ptr_iter2 != sam_rec_pe_ptr.second.end(); ++sam_rec_ptr_iter2)
		{
			long s_1_st, s_1_end, s_2_st, s_2_end;

			if ((*sam_rec_ptr_iter1)->is_fusion)
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->fusion_prefix_st);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->fusion_suffix_end);
			}
			else
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->start);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->end);
			}


			if ((*sam_rec_ptr_iter2)->is_fusion)
			{
				s_2_st = static_cast <long>( (*sam_rec_ptr_iter2)->fusion_prefix_st);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->fusion_suffix_end);
			}
			else
			{
				s_2_st = static_cast <long>((*sam_rec_ptr_iter2)->start);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->end);
			}

			long mate_dist1 = 0; //= (long(s_1_st - s_2_end));

			long mate_dist2 = 0; //= (long(s_1_end - s_2_st));

			bool crossed = false;

			bool not_crossed = false;

			bool is_fusion_paired = false;

			size_t strand_1 = 0, strand_2 = 0;

			size_t prefix_len = 0, suffix_len = 0;

			vector<SpliceWay> left_splice_ways;

			vector<SpliceWay> right_splice_ways;

			string chr1 = "", chr2 = "";

			if (!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec) &&
				!check_overlap((*sam_rec_ptr_iter2)->spliceway_vec, (*sam_rec_ptr_iter1)->spliceway_vec)
				//((*sam_rec_ptr_iter1)->end < (*sam_rec_ptr_iter2)->start/* && (*sam_rec_ptr_iter1)->end >= (*sam_rec_ptr_iter2)->end*/ ||
				//(*sam_rec_ptr_iter2)->end < (*sam_rec_ptr_iter1)->start/* && (*sam_rec_ptr_iter2)->end >= (*sam_rec_ptr_iter1)->end*/ )
				)
			{
				//((*sam_rec_ptr_iter1)->start <= (*sam_rec_ptr_iter2)->start && (*sam_rec_ptr_iter1)->end >= (*sam_rec_ptr_iter2)->end ) ||
				//((*sam_rec_ptr_iter2)->start <= (*sam_rec_ptr_iter1)->start && (*sam_rec_ptr_iter2)->end >= (*sam_rec_ptr_iter1)->end ))

				not_crossed = true;

				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				if (m_do_filter & 16 && (((strand_1 ^ strand_2) & IS_REVERSE) == 0) ||
					(abs(mate_dist1) >= m_max_pair_dist && abs(mate_dist2) >= m_max_pair_dist))
				{
					fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
						(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
						chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));
				}
			}
			else if (m_do_filter & 16 &&
					!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
					(*sam_rec_ptr_iter1)->chrom_name != (*sam_rec_ptr_iter2)->chrom_name/* ||
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec)*/)
			{
				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
					(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
					chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));
			}
			else if ((*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion)
			{
				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->end && abs(long(s_1_st - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					is_fusion_paired = true;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1 + abs(mate_dist1) + (*sam_rec_ptr_iter2)->spliceway_vec.back().second;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end;
					
				}

				else if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1 + abs(mate_dist1) + (*sam_rec_ptr_iter2)->spliceway_vec.front().second;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2;

					is_fusion_paired = true;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start;
				}

				else if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->start && abs(long(s_1_end - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1;

					is_fusion_paired = true;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2 + abs(mate_dist2) + (*sam_rec_ptr_iter2)->spliceway_vec.front().second;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->start2 - (*sam_rec_ptr_iter2)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start2;

				}

				else if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2 + abs(mate_dist2) + (*sam_rec_ptr_iter2)->spliceway_vec.back().second;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));



					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->end2 - (*sam_rec_ptr_iter2)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end2;
				}

			}
			else if (!(*sam_rec_ptr_iter1)->is_fusion && (*sam_rec_ptr_iter2)->is_fusion)
			{
				if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->spliceway_vec.back().second + abs(mate_dist2) + (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;

				}

				else if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->spliceway_vec.front().second + abs(mate_dist1) + (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));


					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;
				}

				else if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;

					prefix_len = (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2 + (*sam_rec_ptr_iter1)->spliceway_vec.front().second + abs(mate_dist1);

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));


					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start2;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter2)->fusion_suffix_end;
				}

				else if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;

					prefix_len = (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2 + (*sam_rec_ptr_iter1)->spliceway_vec.back().second + abs(mate_dist2);

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end2;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter2)->fusion_suffix_end;
				}
			}
			else if ((*sam_rec_ptr_iter1)->is_fusion && (*sam_rec_ptr_iter2)->is_fusion)
			{
				continue;
				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
					(*sam_rec_ptr_iter1)->strand1 == '+' && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->fusion_prefix_st)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
					(*sam_rec_ptr_iter1)->strand1 == '-' && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->fusion_prefix_st)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
					(*sam_rec_ptr_iter1)->strand2 == '+' && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->fusion_suffix_end)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
					(*sam_rec_ptr_iter1)->strand2 == '-' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->fusion_suffix_end)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;
				}

				//if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2 && 
				//	(*sam_rec_ptr_iter1)->strand1 == '+' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
				//	(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->fusion_prefix_st)
				//	not_crossed = true;

				//if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				//	(*sam_rec_ptr_iter1)->strand1 == '-' && (*sam_rec_ptr_iter2)->strand1 == '+' &&
				//	(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->fusion_prefix_st)
				//	not_crossed = true;

				//if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
				//	(*sam_rec_ptr_iter1)->strand2 == '+' && (*sam_rec_ptr_iter2)->strand2 == '-' &&
				//	(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->fusion_suffix_end)
				//	not_crossed = true;

				//if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
				//	(*sam_rec_ptr_iter1)->strand2 == '-' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
				//	(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->fusion_suffix_end)
				//	not_crossed = true;
			}

			if (!not_crossed)
				crossed = true;

			if (/*!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name &&*/
				(((strand_1 ^ strand_2) & IS_REVERSE) || is_fusion_paired)&&
				!crossed &&
				(abs(mate_dist1) < m_max_pair_dist ||
				 abs(mate_dist2) < m_max_pair_dist) /*&&
				 (*sam_rec_ptr_iter1)->filter_type.find(NOT_FILTERED) != (*sam_rec_ptr_iter1)->filter_type.end() && 
				 (*sam_rec_ptr_iter2)->filter_type.find(NOT_FILTERED) != (*sam_rec_ptr_iter2)->filter_type.end()*/)
			{
				paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
					(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
					chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));

				//m_ofs_filtered_alignment << (*sam_rec_ptr_iter1)->tostring(sam_rec_pe_ptr.first.size(), sam_rec_ptr_iter1 - sam_rec_pe_ptr.first.begin() + 1)<<endl;
				//m_ofs_filtered_alignment << (*sam_rec_ptr_iter2)->tostring(sam_rec_pe_ptr.second.size(), sam_rec_ptr_iter2 - sam_rec_pe_ptr.second.begin() + 1)<<endl;

				paired = true;
			}

		}
	}

	FilterPairedSamRec(paired_reads_ptr, sam_rec_pe_ptr);

	//m_ofs_filtered_alignment << endl;
	return paired;

}



bool AlignmentHandler::EstablishPairingStatic(void * str)
{
	//chrom pairid strand offset line
	//hash_map<string, hash_map<size_t, hash_map<unsigned short, map<size_t, vector<SamRec* > > > > > mapped_reads;

	//map<int, vector<pair<char, string> > > paired_reads;

	EstPairingStruct* est_pairing_ptr = (EstPairingStruct*) str;

	vector<PairedSamRec> paired_reads_ptr;

	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(est_pairing_ptr->sam_rec_pe_ptr);

	vector<PairedSamRec>& fusion_paired_reads_ptr = *(est_pairing_ptr->fusion_paired_reads_ptr);

	int do_filter = est_pairing_ptr->do_filter;
	
	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	//bool end1_mapped = false, end2_mapped = false;

	vector<SamRec*>::iterator sam_rec_ptr_iter1, sam_rec_ptr_iter2;

	for (sam_rec_ptr_iter1 = sam_rec_pe_ptr.first.begin(); sam_rec_ptr_iter1 != sam_rec_pe_ptr.first.end(); ++sam_rec_ptr_iter1)
	{
		for (sam_rec_ptr_iter2 = sam_rec_pe_ptr.second.begin(); sam_rec_ptr_iter2 != sam_rec_pe_ptr.second.end(); ++sam_rec_ptr_iter2)
		{
			long s_1_st, s_1_end, s_2_st, s_2_end;

			if ((*sam_rec_ptr_iter1)->is_fusion)
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->fusion_prefix_st);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->fusion_suffix_end);
			}
			else
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->start);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->end);
			}


			if ((*sam_rec_ptr_iter2)->is_fusion)
			{
				s_2_st = static_cast <long>( (*sam_rec_ptr_iter2)->fusion_prefix_st);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->fusion_suffix_end);
			}
			else
			{
				s_2_st = static_cast <long>((*sam_rec_ptr_iter2)->start);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->end);
			}

			long mate_dist1 = 0; //= (long(s_1_st - s_2_end));

			long mate_dist2 = 0; //= (long(s_1_end - s_2_st));

			bool crossed = false;

			bool not_crossed = false;

			bool is_fusion_paired = false;

			bool fusion_prefix_matched = false, fusion_suffix_matched = false;

			size_t strand_1 = 0, strand_2 = 0;

			size_t prefix_len = 0, suffix_len = 0;

			vector<SpliceWay> left_splice_ways;

			vector<SpliceWay> right_splice_ways;

			string chr1 = "", chr2 = "";

			if (!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec) && 
				!check_overlap((*sam_rec_ptr_iter2)->spliceway_vec, (*sam_rec_ptr_iter1)->spliceway_vec)
				//((*sam_rec_ptr_iter1)->end < (*sam_rec_ptr_iter2)->start/* && (*sam_rec_ptr_iter1)->end >= (*sam_rec_ptr_iter2)->end*/ ||
				//(*sam_rec_ptr_iter2)->end < (*sam_rec_ptr_iter1)->start/* && (*sam_rec_ptr_iter2)->end >= (*sam_rec_ptr_iter1)->end*/ )
				)
			{
				//((*sam_rec_ptr_iter1)->start <= (*sam_rec_ptr_iter2)->start && (*sam_rec_ptr_iter1)->end >= (*sam_rec_ptr_iter2)->end ) ||
				//((*sam_rec_ptr_iter2)->start <= (*sam_rec_ptr_iter1)->start && (*sam_rec_ptr_iter2)->end >= (*sam_rec_ptr_iter1)->end ))

				not_crossed = true;

				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				//if (m_do_filter & 16 && (((strand_1 ^ strand_2) & IS_REVERSE) == 0) ||
				//	(abs(mate_dist1) >= m_max_pair_dist && abs(mate_dist2) >= m_max_pair_dist))
				//{
				//	//cout << "fusion paired"<<endl;

				//	//cout << (*sam_rec_ptr_iter1)->tostring(0, 0) << endl;

				//	//cout << (*sam_rec_ptr_iter2)->tostring(0, 0) << endl;

				//	fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
				//		(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
				//		chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));
				//}
			}
			else if (m_do_filter & 16 &&
					!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
					(*sam_rec_ptr_iter1)->chrom_name != (*sam_rec_ptr_iter2)->chrom_name/* ||
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec)*/)
			{
				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				//cout << "fusion paired"<<endl;

				//cout << (*sam_rec_ptr_iter1)->tostring(0, 0) << endl;

				//cout << (*sam_rec_ptr_iter2)->tostring(0, 0) << endl;

				//fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
				//	(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
				//	chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));
			}
			else if ((*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion)
			{
				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->end && abs(long(s_1_st - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1 + abs(mate_dist1) + (*sam_rec_ptr_iter2)->spliceway_vec.back().second;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end;
					
				}

				else if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1 + abs(mate_dist1) + (*sam_rec_ptr_iter2)->spliceway_vec.front().second;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start;
				}

				else if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->start && abs(long(s_1_end - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2 + abs(mate_dist2) + (*sam_rec_ptr_iter2)->spliceway_vec.front().second;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->start2 - (*sam_rec_ptr_iter2)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start2;

				}

				else if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2 + abs(mate_dist2) + (*sam_rec_ptr_iter2)->spliceway_vec.back().second;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->end2 - (*sam_rec_ptr_iter2)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end2;
				}

			}
			else if (!(*sam_rec_ptr_iter1)->is_fusion && (*sam_rec_ptr_iter2)->is_fusion)
			{
				if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->spliceway_vec.back().second + abs(mate_dist2) + (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;

				}

				else if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->spliceway_vec.front().second + abs(mate_dist1) + (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));


					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;
				}

				else if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;

					prefix_len = (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2 + (*sam_rec_ptr_iter1)->spliceway_vec.front().second + abs(mate_dist1);

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));


					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start2;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter2)->fusion_suffix_end;
				}

				else if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist)
				{
					not_crossed = true;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;

					prefix_len = (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2 + (*sam_rec_ptr_iter1)->spliceway_vec.back().second + abs(mate_dist2);

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end2;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					//(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter2)->fusion_suffix_end;
				}
			}
			else if ((*sam_rec_ptr_iter1)->is_fusion && (*sam_rec_ptr_iter2)->is_fusion)
			{
				continue;
				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
					(*sam_rec_ptr_iter1)->strand1 == '+' && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->fusion_prefix_st)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
					(*sam_rec_ptr_iter1)->strand1 == '-' && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->fusion_prefix_st)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
					(*sam_rec_ptr_iter1)->strand2 == '+' && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->fusion_suffix_end)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
					(*sam_rec_ptr_iter1)->strand2 == '-' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->fusion_suffix_end)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;
				}

				//if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2 && 
				//	(*sam_rec_ptr_iter1)->strand1 == '+' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
				//	(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->fusion_prefix_st)
				//	not_crossed = true;

				//if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				//	(*sam_rec_ptr_iter1)->strand1 == '-' && (*sam_rec_ptr_iter2)->strand1 == '+' &&
				//	(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->fusion_prefix_st)
				//	not_crossed = true;

				//if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
				//	(*sam_rec_ptr_iter1)->strand2 == '+' && (*sam_rec_ptr_iter2)->strand2 == '-' &&
				//	(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->fusion_suffix_end)
				//	not_crossed = true;

				//if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
				//	(*sam_rec_ptr_iter1)->strand2 == '-' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
				//	(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->fusion_suffix_end)
				//	not_crossed = true;
			}

			if (!not_crossed)
				crossed = true;

			if (/*!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name &&*/
				(((strand_1 ^ strand_2) & IS_REVERSE) || is_fusion_paired)&&
				!crossed &&
				(abs(mate_dist1) < m_max_pair_dist ||
				 abs(mate_dist2) < m_max_pair_dist) /*&&
				 (*sam_rec_ptr_iter1)->filter_type.find(NOT_FILTERED) != (*sam_rec_ptr_iter1)->filter_type.end() && 
				 (*sam_rec_ptr_iter2)->filter_type.find(NOT_FILTERED) != (*sam_rec_ptr_iter2)->filter_type.end()*/)
			{
				paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
					(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
					chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, 
					fusion_prefix_matched, fusion_suffix_matched));

				//m_ofs_filtered_alignment << (*sam_rec_ptr_iter1)->tostring(sam_rec_pe_ptr.first.size(), sam_rec_ptr_iter1 - sam_rec_pe_ptr.first.begin() + 1)<<endl;
				//m_ofs_filtered_alignment << (*sam_rec_ptr_iter2)->tostring(sam_rec_pe_ptr.second.size(), sam_rec_ptr_iter2 - sam_rec_pe_ptr.second.begin() + 1)<<endl;

				paired = true;
			}

		}
	}

	FilterPairingStruct filter_pairing(&paired_reads_ptr, &sam_rec_pe_ptr, do_filter);

	FilterPairedSamRecStatic( (void *) (&filter_pairing));
	//FilterPairedSamRec(paired_reads_ptr, sam_rec_pe_ptr);

	//m_ofs_filtered_alignment << endl;
	//return paired;

	//cout << "EstablishPairing:"<<paired << endl;

	return paired;

}

bool AlignmentHandler::EstablishFusionPairingStatic(void * str)
{
	//chrom pairid strand offset line
	//hash_map<string, hash_map<size_t, hash_map<unsigned short, map<size_t, vector<SamRec* > > > > > mapped_reads;

	//map<int, vector<pair<char, string> > > paired_reads;

	EstPairingStruct* est_pairing_ptr = (EstPairingStruct*) str;

	vector<PairedSamRec> paired_reads_ptr;

	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(est_pairing_ptr->sam_rec_pe_ptr);

	vector<PairedSamRec>& fusion_paired_reads_ptr = *(est_pairing_ptr->fusion_paired_reads_ptr);

	int do_filter = est_pairing_ptr->do_filter;
	
	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	//bool end1_mapped = false, end2_mapped = false;

	vector<SamRec*>::iterator sam_rec_ptr_iter1, sam_rec_ptr_iter2;

	for (sam_rec_ptr_iter1 = sam_rec_pe_ptr.first.begin(); sam_rec_ptr_iter1 != sam_rec_pe_ptr.first.end(); ++sam_rec_ptr_iter1)
	{
		for (sam_rec_ptr_iter2 = sam_rec_pe_ptr.second.begin(); sam_rec_ptr_iter2 != sam_rec_pe_ptr.second.end(); ++sam_rec_ptr_iter2)
		{
			long s_1_st, s_1_end, s_2_st, s_2_end;

			if ((*sam_rec_ptr_iter1)->is_fusion)
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->fusion_prefix_st);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->fusion_suffix_end);
			}
			else
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->start);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->end);
			}


			if ((*sam_rec_ptr_iter2)->is_fusion)
			{
				s_2_st = static_cast <long>( (*sam_rec_ptr_iter2)->fusion_prefix_st);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->fusion_suffix_end);
			}
			else
			{
				s_2_st = static_cast <long>((*sam_rec_ptr_iter2)->start);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->end);
			}

			long mate_dist1 = 0; //= (long(s_1_st - s_2_end));

			long mate_dist2 = 0; //= (long(s_1_end - s_2_st));

			bool crossed = false;

			bool not_crossed = false;

			bool is_fusion_paired = false;

			bool fusion_prefix_matched = false, fusion_suffix_matched = false;

			size_t strand_1 = 0, strand_2 = 0;

			size_t prefix_len = 0, suffix_len = 0;

			vector<SpliceWay> left_splice_ways;

			vector<SpliceWay> right_splice_ways;

			string chr1 = "", chr2 = "";

			if (!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec) && 
				!check_overlap((*sam_rec_ptr_iter2)->spliceway_vec, (*sam_rec_ptr_iter1)->spliceway_vec)
				//((*sam_rec_ptr_iter1)->end < (*sam_rec_ptr_iter2)->start/* && (*sam_rec_ptr_iter1)->end >= (*sam_rec_ptr_iter2)->end*/ ||
				//(*sam_rec_ptr_iter2)->end < (*sam_rec_ptr_iter1)->start/* && (*sam_rec_ptr_iter2)->end >= (*sam_rec_ptr_iter1)->end*/ )
				)
			{
				//((*sam_rec_ptr_iter1)->start <= (*sam_rec_ptr_iter2)->start && (*sam_rec_ptr_iter1)->end >= (*sam_rec_ptr_iter2)->end ) ||
				//((*sam_rec_ptr_iter2)->start <= (*sam_rec_ptr_iter1)->start && (*sam_rec_ptr_iter2)->end >= (*sam_rec_ptr_iter1)->end ))

				not_crossed = true;

				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				if (m_do_filter & 16 && (((strand_1 ^ strand_2) & IS_REVERSE) == 0) ||
					(abs(mate_dist1) >= m_max_pair_dist && abs(mate_dist2) >= m_max_pair_dist))
				{
					//cout << "fusion paired"<<endl;

					//cout << (*sam_rec_ptr_iter1)->tostring(0, 0) << endl;

					//cout << (*sam_rec_ptr_iter2)->tostring(0, 0) << endl;

					fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
						(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
						chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));
				}
			}
			else if (m_do_filter & 16 &&
					!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
					(*sam_rec_ptr_iter1)->chrom_name != (*sam_rec_ptr_iter2)->chrom_name/* ||
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec)*/)
			{
				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				//cout << "fusion paired"<<endl;

				//cout << (*sam_rec_ptr_iter1)->tostring(0, 0) << endl;

				//cout << (*sam_rec_ptr_iter2)->tostring(0, 0) << endl;

				fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
					(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
					chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false));
			}
			

		}
	}

	//FilterPairingStruct filter_pairing(&paired_reads_ptr, &sam_rec_pe_ptr, do_filter);

	//FilterPairedSamRecStatic( (void *) (&filter_pairing));
	////FilterPairedSamRec(paired_reads_ptr, sam_rec_pe_ptr);

	////m_ofs_filtered_alignment << endl;
	////return paired;

	////cout << "EstablishPairing:"<<paired << endl;

	return paired;

}

void 
AlignmentHandler::FilterByMisMatch(vector<SamRec*>& sam_rec_ptr)
{
	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_mis);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->mis_match != (*(sam_rec_ptr_iter-1))->mis_match)
		{
			sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());
			break;
		}
	}
}

void* AlignmentHandler::FilterByMisMatchStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_mis);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->mis_match != (*(sam_rec_ptr_iter-1))->mis_match)
		{
			sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());
			break;
		}
	}
}

void 
AlignmentHandler::RemoveDup(vector<SamRec*>& sam_rec_ptr)
{
	vector<SamRec*> filtered_sam_rec_ptr;

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_dup);

	vector<SamRec*>::iterator sam_rec_iter;

	filtered_sam_rec_ptr.push_back(*sam_rec_ptr.begin());

	for (sam_rec_iter = sam_rec_ptr.begin() + 1; sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->chrom_name != (*(sam_rec_iter - 1))->chrom_name ||
			(*sam_rec_iter)->start != (*(sam_rec_iter - 1))->start || 
			(*sam_rec_iter)->ori_splice_way != (*(sam_rec_iter - 1))->ori_splice_way || 
			(*sam_rec_iter)->chrom_name2 != (*(sam_rec_iter - 1))->chrom_name2 ||
			(*sam_rec_iter)->start2 != (*(sam_rec_iter - 1))->start2 || 
			(*sam_rec_iter)->ori_splice_way2 != (*(sam_rec_iter - 1))->ori_splice_way2)
			filtered_sam_rec_ptr.push_back(*sam_rec_iter);
		else
		{
			//cout <<"duplication "<<endl;
			//cout << (*sam_rec_iter)->tostring(0,0)<<endl;
			//cout << (*(sam_rec_iter-1))->tostring(0,0)<<endl;
			//getchar();
		}
	}

	sam_rec_ptr = filtered_sam_rec_ptr;
}

void* AlignmentHandler::RemoveDupStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	vector<SamRec*> filtered_sam_rec_ptr;

	//cout << "sort"<<endl;

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_dup);

	vector<SamRec*>::iterator sam_rec_iter;

	filtered_sam_rec_ptr.push_back(*sam_rec_ptr.begin());

	//cout << "remove"<<endl;

	//cout << "sam_rec_ptr size:" << sam_rec_ptr.size()<< endl;

	for (sam_rec_iter = sam_rec_ptr.begin() + 1; sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->chrom_name != (*(sam_rec_iter - 1))->chrom_name ||
			(*sam_rec_iter)->start != (*(sam_rec_iter - 1))->start || 
			(*sam_rec_iter)->ori_splice_way != (*(sam_rec_iter - 1))->ori_splice_way ||
			(*sam_rec_iter)->chrom_name2 != (*(sam_rec_iter - 1))->chrom_name2 ||
			(*sam_rec_iter)->start2 != (*(sam_rec_iter - 1))->start2 || 
			(*sam_rec_iter)->ori_splice_way2 != (*(sam_rec_iter - 1))->ori_splice_way2)
			filtered_sam_rec_ptr.push_back(*sam_rec_iter);
		else
		{
			//cout <<"duplication "<<endl;
			//cout << (*sam_rec_iter)->tostring(0,0)<<endl;
			//cout << (*(sam_rec_iter-1))->tostring(0,0)<<endl;
			//getchar();
		}
	}

	sam_rec_ptr = filtered_sam_rec_ptr;

	//cout << "finish"<<endl;
}

void 
AlignmentHandler::AllocateHitsMemory()
{
	ifstream ifs(m_chrom_size_file);

	if (ifs.is_open())
	{
		vector<bool> dumy;

		while (!ifs.eof() )
		{
			string line;

			getline(ifs, line);

			if (line.empty()/* || line[0] == '@'*/)
				continue;

			char chrom_name[1000];

			size_t chrom_size;

			sscanf(line.c_str(), "%s\t%llu", chrom_name, &chrom_size);

			m_mapped_pos.insert(make_pair(chrom_name, dumy));

			//cout << line << endl;
			//cout << "resize " << endl;
			m_mapped_pos[chrom_name].resize(chrom_size, false);
			//cout << "resize finished" << endl;
		}
	}
	else
	{
		cerr <<"Can't open file: " << m_chrom_size_file << endl;
		exit(0);
	}
}

void
AlignmentHandler::Hits2ExpressedRegions()
{
	hash_map<string, vector<bool> >::iterator mapped_pos_iter;
	
	for (mapped_pos_iter = m_mapped_pos.begin(); mapped_pos_iter != m_mapped_pos.end(); ++mapped_pos_iter)
	{
		vector<pair<size_t, size_t> >& cur_islands = m_expressed_regions[mapped_pos_iter->first];

		vector<bool>& cur_hits = mapped_pos_iter->second;

		bool inisland = false;

		size_t tmpstart = 0, tmpend = 0;

		for (size_t i = 0; i < cur_hits.size(); ++i)
		{
			if (cur_hits[i] && !inisland)
			{
				inisland = true;
				tmpstart = i;
			}
			else if(!cur_hits[i] && inisland)
			{
				tmpend = i;

				inisland = false;

				if (cur_islands.size() && (tmpstart - cur_islands[cur_islands.size() - 1].second) <= boundary)
					cur_islands[cur_islands.size() - 1].second = tmpend/* + 45*/;
				else
				{
					//if ((int)tmpstart - 45 < 0)
					//	cur_islands.push_back(make_pair(1, tmpend + 45));
					//else
					cur_islands.push_back(make_pair(tmpstart/* - 45*/, tmpend/* + 45*/));
				}
			}
		}

		m_disjointset_map[mapped_pos_iter->first].Resize(cur_islands.size());
	}

	string islands_file = m_filtered_alignment_file; islands_file.append(".islands");

	ofstream test_islands(islands_file.c_str());

	hash_map<string, vector<pair<size_t, size_t> > >::iterator chrom_islands_iter;

	for (chrom_islands_iter = m_expressed_regions.begin(); chrom_islands_iter != m_expressed_regions.end(); ++chrom_islands_iter)
	{
		vector<pair<size_t, size_t> >::iterator island_iter;

		for (island_iter = chrom_islands_iter->second.begin(); island_iter != chrom_islands_iter->second.end(); ++island_iter)
		{
			test_islands << chrom_islands_iter->first << '\t' << island_iter->first << '\t' << island_iter->second << '\t' <<island_iter->second - island_iter->first<< endl;
		}
	}

	//hash_map<string, vector<pair<size_t, size_t> > > m_expressed_regions;
}

bool
AlignmentHandler::FindRegion(pair<size_t, size_t>& cur_region, vector<pair<size_t, size_t> >& sorted_regions, size_t& find_region_idx)
{
	vector<pair<size_t, size_t> >::iterator express_region_iter;

	express_region_iter = lower_bound(sorted_regions.begin(), sorted_regions.end(), cur_region, compare_pair_region);

	if (express_region_iter != sorted_regions.begin())
		--express_region_iter;

	bool find_region = false;

	vector<pair<size_t, size_t> >::iterator cur_express_region_iter = express_region_iter;

	while (true)
	{
		if (cur_express_region_iter == sorted_regions.end() || cur_express_region_iter->first > cur_region.first)
		{
			break;
		}
		else if (cur_express_region_iter->first  <= cur_region.first + 3 && cur_express_region_iter->second + 3 >= cur_region.first)
		{
			find_region = true;
			break;
		}

		++cur_express_region_iter;
	}

	if (find_region)
		find_region_idx = cur_express_region_iter - sorted_regions.begin();

	return find_region;
}

void
AlignmentHandler::UnionSets(JunctionHandler* junc_handler, hash_map<string, vector<pair<size_t, size_t> > >* expressed_regions)
{
	junc_handler->LoadJuncToSortVec();

	junc_handler->SortJunc();

	hash_map<string, vector<pair<size_t, size_t> > >::iterator expressed_regions_chrom_iter;

	hash_map<size_t, UnionExpressedRegions > dumy;

	for (expressed_regions_chrom_iter = expressed_regions->begin(); expressed_regions_chrom_iter != expressed_regions->end(); ++expressed_regions_chrom_iter)
	{
		if (junc_handler->m_junc_sort.find(expressed_regions_chrom_iter->first) == junc_handler->m_junc_sort.end())
			continue;

		vector<JunctionSeed*>& sorted_juncs = (junc_handler->m_junc_sort.find(expressed_regions_chrom_iter->first))->second;

		vector<pair<size_t, size_t> >& sorted_regions = expressed_regions_chrom_iter->second;

		DisjointSet& cur_disjointSet = (m_disjointset_map.find(expressed_regions_chrom_iter->first))->second;

		vector<JunctionSeed*>::iterator junc_seed_iter;

		for (junc_seed_iter = sorted_juncs.begin(); junc_seed_iter != sorted_juncs.end(); ++junc_seed_iter)
		{
			pair<size_t, size_t> cur_junc_region_start = make_pair((*junc_seed_iter)->m_start, (*junc_seed_iter)->m_start);

			pair<size_t, size_t> cur_junc_region_end = make_pair((*junc_seed_iter)->m_end, (*junc_seed_iter)->m_end);

			size_t start_region_idx = -1, end_region_idx = -1;

			bool find_start_region = FindRegion(cur_junc_region_start, sorted_regions, start_region_idx);

			bool find_end_region = FindRegion(cur_junc_region_end, sorted_regions, end_region_idx);

			if (find_start_region && find_end_region && start_region_idx != end_region_idx)
			{
				cur_disjointSet.Union(start_region_idx, end_region_idx);
			}
			else if (start_region_idx != end_region_idx)
			{
				cout << "region not found: "<<expressed_regions_chrom_iter->first << '\t' << (*junc_seed_iter)->m_start << '\t' << (*junc_seed_iter)->m_end<< endl;//cout error information
				continue;
			}
		}

		m_unioned_expressed_regions_map.insert(make_pair(expressed_regions_chrom_iter->first, dumy));

		hash_map<size_t, UnionExpressedRegions >& cur_unioned_expressed_regions = m_unioned_expressed_regions_map.find(expressed_regions_chrom_iter->first)->second;

		//vector<pair<size_t, size_t> >::iterator sorted_regions_iter;
		for (size_t i = 0; i < sorted_regions.size(); ++i)
		{
			size_t source_id = cur_disjointSet.FindSet(i);

			cur_unioned_expressed_regions[source_id].m_unioned_regions_ids.push_back(i);

			cur_unioned_expressed_regions[source_id].m_unioned_regions.push_back(sorted_regions[i]);

			cur_unioned_expressed_regions[source_id].m_unioned_region_length += sorted_regions[i].second - sorted_regions[i].first + 1;
		}
	}

	string expressed_regions_file = m_filtered_alignment_file; expressed_regions_file.append(".expressed_regions");

	ofstream expressed_regions_ofs(expressed_regions_file.c_str());

	hash_map<string, hash_map<size_t, UnionExpressedRegions > >::iterator chrom_expressed_regions_iter;

	hash_map<size_t, UnionExpressedRegions >::iterator source_expressed_regions_iter;

	for (chrom_expressed_regions_iter = m_unioned_expressed_regions_map.begin(); chrom_expressed_regions_iter != m_unioned_expressed_regions_map.end(); ++chrom_expressed_regions_iter)
	{
		for (source_expressed_regions_iter = chrom_expressed_regions_iter->second.begin(); source_expressed_regions_iter != chrom_expressed_regions_iter->second.end(); ++source_expressed_regions_iter)
		{
			expressed_regions_ofs << chrom_expressed_regions_iter->first << '\t' << source_expressed_regions_iter->first << '\t' << source_expressed_regions_iter->second.m_unioned_region_length << '\t';
			
			for (size_t i = 0; i < source_expressed_regions_iter->second.m_unioned_regions.size(); ++i)
				expressed_regions_ofs<<source_expressed_regions_iter->second.m_unioned_regions[i].first<<'\t'
				<<source_expressed_regions_iter->second.m_unioned_regions[i].second<< '\t';

			expressed_regions_ofs << endl;
		}
	}
}

void
AlignmentHandler::FindFusionJuncRegionVec(JunctionHandler* junc_handler, vector<PairedSamRec>& fusion_paired_reads_ptr)
{
	//junc_handler->LoadFusionJuncToSortVec();

	//junc_handler->SortFusionJunc();

	vector<PairedSamRec>::iterator paired_sam_rec_iter;

	for (paired_sam_rec_iter = fusion_paired_reads_ptr.begin(); paired_sam_rec_iter != fusion_paired_reads_ptr.end(); ++paired_sam_rec_iter)
	{
		FusionJuncRegion cur_fusion_junc_region;
		
		string chrom_name;

		if (paired_sam_rec_iter->paired_sam_rec.first->chrom_name < paired_sam_rec_iter->paired_sam_rec.second->chrom_name ||
		    (paired_sam_rec_iter->paired_sam_rec.first->chrom_name == paired_sam_rec_iter->paired_sam_rec.second->chrom_name &&
			 paired_sam_rec_iter->paired_sam_rec.first->start <= paired_sam_rec_iter->paired_sam_rec.second->start)
			 )
		{
			chrom_name = paired_sam_rec_iter->paired_sam_rec.first->chrom_name + "~" + paired_sam_rec_iter->paired_sam_rec.second->chrom_name;

			cur_fusion_junc_region.set(paired_sam_rec_iter->paired_sam_rec.first->start, paired_sam_rec_iter->paired_sam_rec.first->end,
				paired_sam_rec_iter->paired_sam_rec.second->start, paired_sam_rec_iter->paired_sam_rec.second->end);
		}
		else
		{
			chrom_name = paired_sam_rec_iter->paired_sam_rec.second->chrom_name + "~" + paired_sam_rec_iter->paired_sam_rec.first->chrom_name;

			cur_fusion_junc_region.set(paired_sam_rec_iter->paired_sam_rec.second->start, paired_sam_rec_iter->paired_sam_rec.second->end,
				paired_sam_rec_iter->paired_sam_rec.first->start, paired_sam_rec_iter->paired_sam_rec.first->end);
		}

		if (junc_handler->m_fuson_junc_sort.find(chrom_name) != junc_handler->m_fuson_junc_sort.end())
		{
			vector<FusionJuncRegion>& cur_chrom_fusion_junc_regions = (junc_handler->m_fuson_junc_sort.find(chrom_name))->second;

			FindFusionJuncRegion(cur_fusion_junc_region, *paired_sam_rec_iter, cur_chrom_fusion_junc_regions);
		}
	}
}

void*
AlignmentHandler::FindFusionJuncRegionVecStatic(void * str)
{
	JunctionHandler* junc_handler = m_junction_handler_prev_ptr;
	
	vector<PairedSamRec>& fusion_paired_reads_ptr = *((vector<PairedSamRec>*)str);
	//junc_handler->LoadFusionJuncToSortVec();

	//cout << "fusion_paired_reads_ptr:"<<fusion_paired_reads_ptr.size()<<endl;


	//junc_handler->SortFusionJunc();

	vector<PairedSamRec>::iterator paired_sam_rec_iter;

	for (paired_sam_rec_iter = fusion_paired_reads_ptr.begin(); paired_sam_rec_iter != fusion_paired_reads_ptr.end(); ++paired_sam_rec_iter)
	{
		FusionJuncRegion cur_fusion_junc_region;
		
		string chrom_name;

		if (paired_sam_rec_iter->paired_sam_rec.first->chrom_name < paired_sam_rec_iter->paired_sam_rec.second->chrom_name ||
		    (paired_sam_rec_iter->paired_sam_rec.first->chrom_name == paired_sam_rec_iter->paired_sam_rec.second->chrom_name &&
			 paired_sam_rec_iter->paired_sam_rec.first->start <= paired_sam_rec_iter->paired_sam_rec.second->start)
			 )
		{
			chrom_name = paired_sam_rec_iter->paired_sam_rec.first->chrom_name + "~" + paired_sam_rec_iter->paired_sam_rec.second->chrom_name;

			cur_fusion_junc_region.set(paired_sam_rec_iter->paired_sam_rec.first->start, paired_sam_rec_iter->paired_sam_rec.first->end,
				paired_sam_rec_iter->paired_sam_rec.second->start, paired_sam_rec_iter->paired_sam_rec.second->end);
		}
		else
		{
			chrom_name = paired_sam_rec_iter->paired_sam_rec.second->chrom_name + "~" + paired_sam_rec_iter->paired_sam_rec.first->chrom_name;

			cur_fusion_junc_region.set(paired_sam_rec_iter->paired_sam_rec.second->start, paired_sam_rec_iter->paired_sam_rec.second->end,
				paired_sam_rec_iter->paired_sam_rec.first->start, paired_sam_rec_iter->paired_sam_rec.first->end);
		}

		if (junc_handler->m_fuson_junc_sort.find(chrom_name) != junc_handler->m_fuson_junc_sort.end())
		{
			vector<FusionJuncRegion>& cur_chrom_fusion_junc_regions = (junc_handler->m_fuson_junc_sort.find(chrom_name))->second;

			FindFusionJuncRegionStruct ffr(&cur_fusion_junc_region, &(*paired_sam_rec_iter), &cur_chrom_fusion_junc_regions);

			FindFusionJuncRegionStatic((void*)&ffr);
		}
	}
	
	return 0;
}
void 
AlignmentHandler::FindFusionJuncRegion(FusionJuncRegion& cur_fusion_region, PairedSamRec& cur_paired_sam_rec, vector<FusionJuncRegion>& sorted_fusion_regions)
{
	vector<FusionJuncRegion>::iterator fusion_region_iter;

	fusion_region_iter = lower_bound(sorted_fusion_regions.begin(), sorted_fusion_regions.end(), cur_fusion_region, comp_fusion_junc_sort);

	if (fusion_region_iter != sorted_fusion_regions.begin())
		--fusion_region_iter;

	while(fusion_region_iter != sorted_fusion_regions.begin() && fusion_region_iter->m_doner_end > cur_fusion_region.m_doner_end)
		--fusion_region_iter;

	bool find_region = false;

	size_t regions_count = 0;

	vector<FusionJuncRegion>::iterator cur_fusion_region_iter = fusion_region_iter;

	while (true)
	{
		if (cur_fusion_region_iter == sorted_fusion_regions.end() || cur_fusion_region_iter->m_doner_st > cur_fusion_region.m_doner_end)
		{
			break;
		}
		else if (cur_fusion_region_iter->m_doner_st <= cur_fusion_region.m_doner_st &&
			cur_fusion_region_iter->m_doner_end >= cur_fusion_region.m_doner_end &&
			cur_fusion_region_iter->m_acceptor_st <= cur_fusion_region.m_acceptor_st &&
			cur_fusion_region_iter->m_acceptor_end >= cur_fusion_region.m_acceptor_end)
		{
			++regions_count;

			//fragment 1
			cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;

			if (cur_paired_sam_rec.paired_sam_rec.second->strand_t & IS_REVERSE)
				cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

			cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_FIRST_END;

			if (cur_paired_sam_rec.paired_sam_rec.first->chrom_name == cur_paired_sam_rec.paired_sam_rec.second->chrom_name)
				cur_paired_sam_rec.paired_sam_rec.first->mate_match = "=";
			else
				cur_paired_sam_rec.paired_sam_rec.first->mate_match = cur_paired_sam_rec.paired_sam_rec.second->chrom_name;

			cur_paired_sam_rec.paired_sam_rec.first->mate_offset = cur_paired_sam_rec.paired_sam_rec.second->start;

			cur_paired_sam_rec.paired_sam_rec.first->mate_diff = static_cast<long> (cur_paired_sam_rec.paired_sam_rec.first->end) - static_cast<long> (cur_paired_sam_rec.paired_sam_rec.second->start);

			//fragment 2
			cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

			if (cur_paired_sam_rec.paired_sam_rec.first->strand_t & IS_REVERSE)
				cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

			cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_FIRST_END;

			if (cur_paired_sam_rec.paired_sam_rec.first->chrom_name == cur_paired_sam_rec.paired_sam_rec.second->chrom_name)
				cur_paired_sam_rec.paired_sam_rec.second->mate_match = "=";
			else
				cur_paired_sam_rec.paired_sam_rec.second->mate_match = cur_paired_sam_rec.paired_sam_rec.first->chrom_name;

			cur_paired_sam_rec.paired_sam_rec.second->mate_offset = cur_paired_sam_rec.paired_sam_rec.first->start;

			cur_paired_sam_rec.paired_sam_rec.second->mate_diff = static_cast<long> (cur_paired_sam_rec.paired_sam_rec.second->end) - static_cast<long> (cur_paired_sam_rec.paired_sam_rec.first->start);

			//cur_fusion_region_iter->m_junc_seed_ptr->m_fusion_encompassing_reads.push_back(make_pair(*(cur_paired_sam_rec.paired_sam_rec.first), *(cur_paired_sam_rec.paired_sam_rec.second)));

			//cur_fusion_region_iter->m_junc_seed_ptr->m_fusion_encompassing_reads_doner.push_back(*(cur_paired_sam_rec.paired_sam_rec.first));//cur_paired_sam_rec.paired_sam_rec.first->tostring(0,0));

			//cur_fusion_region_iter->m_junc_seed_ptr->m_fusion_encompassing_reads_acceptor.push_back(*(cur_paired_sam_rec.paired_sam_rec.second));//cur_paired_sam_rec.paired_sam_rec.second->tostring(0,0));
			//find_region = true;
			//break;
		}

		++cur_fusion_region_iter;
	}

	//if (find_region)
	//	find_region_idx = cur_express_region_iter - sorted_regions.begin();

	//return find_region;
}

void*
AlignmentHandler::FindFusionJuncRegionStatic(void * str)
{
	FindFusionJuncRegionStruct* ffrs = (FindFusionJuncRegionStruct*) str;

	FusionJuncRegion& cur_fusion_region = *(ffrs->cur_fusion_region);

	PairedSamRec& cur_paired_sam_rec = *(ffrs->cur_paired_sam_rec);

	vector<FusionJuncRegion>& sorted_fusion_regions = *(ffrs->sorted_fusion_regions);

	vector<FusionJuncRegion>::iterator fusion_region_iter;

	fusion_region_iter = lower_bound(sorted_fusion_regions.begin(), sorted_fusion_regions.end(), cur_fusion_region, comp_fusion_junc_sort);

	if (fusion_region_iter != sorted_fusion_regions.begin())
		--fusion_region_iter;

	while(fusion_region_iter != sorted_fusion_regions.begin() && fusion_region_iter->m_doner_end > cur_fusion_region.m_doner_end)
		--fusion_region_iter;

	bool find_region = false;

	size_t regions_count = 0;

	vector<FusionJuncRegion>::iterator cur_fusion_region_iter = fusion_region_iter;

	while (true)
	{
		if (cur_fusion_region_iter == sorted_fusion_regions.end() || cur_fusion_region_iter->m_doner_st > cur_fusion_region.m_doner_end)
		{
			break;
		}
		else if (cur_fusion_region_iter->m_doner_st <= cur_fusion_region.m_doner_st &&
			cur_fusion_region_iter->m_doner_end >= cur_fusion_region.m_doner_end &&
			cur_fusion_region_iter->m_acceptor_st <= cur_fusion_region.m_acceptor_st &&
			cur_fusion_region_iter->m_acceptor_end >= cur_fusion_region.m_acceptor_end)
		{
			++regions_count;

			//fragment 1
			cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;

			if (cur_paired_sam_rec.paired_sam_rec.second->strand_t & IS_REVERSE)
				cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

			cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_FIRST_END;

			if (cur_paired_sam_rec.paired_sam_rec.first->chrom_name == cur_paired_sam_rec.paired_sam_rec.second->chrom_name)
				cur_paired_sam_rec.paired_sam_rec.first->mate_match = "=";
			else
				cur_paired_sam_rec.paired_sam_rec.first->mate_match = cur_paired_sam_rec.paired_sam_rec.second->chrom_name;

			cur_paired_sam_rec.paired_sam_rec.first->mate_offset = cur_paired_sam_rec.paired_sam_rec.second->start;

			cur_paired_sam_rec.paired_sam_rec.first->mate_diff = static_cast<long> (cur_paired_sam_rec.paired_sam_rec.first->end) - static_cast<long> (cur_paired_sam_rec.paired_sam_rec.second->start);

			//fragment 2
			cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

			if (cur_paired_sam_rec.paired_sam_rec.first->strand_t & IS_REVERSE)
				cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

			cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_SECOND_END;//IS_FIRST_END;

			if (cur_paired_sam_rec.paired_sam_rec.first->chrom_name == cur_paired_sam_rec.paired_sam_rec.second->chrom_name)
				cur_paired_sam_rec.paired_sam_rec.second->mate_match = "=";
			else
				cur_paired_sam_rec.paired_sam_rec.second->mate_match = cur_paired_sam_rec.paired_sam_rec.first->chrom_name;

			cur_paired_sam_rec.paired_sam_rec.second->mate_offset = cur_paired_sam_rec.paired_sam_rec.first->start;

			cur_paired_sam_rec.paired_sam_rec.second->mate_diff = static_cast<long> (cur_paired_sam_rec.paired_sam_rec.second->end) - static_cast<long> (cur_paired_sam_rec.paired_sam_rec.first->start);

			//if ( cur_paired_sam_rec.paired_sam_rec.second->tag_base_name == "UNC9-SN296_125:8:1:11416:65749/")
			//{
			//	cout << "UNC9-SN296_125:8:1:11416:65749/ is being encompmassed"<<endl;

			//	cout << cur_paired_sam_rec.paired_sam_rec.second->tostring(0, 0)<<endl;
			//}

			//if ( cur_paired_sam_rec.paired_sam_rec.first->tag_base_name == "UNC9-SN296_125:8:1:11416:65749/")
			//{
			//	cout << "UNC9-SN296_125:8:1:11416:65749/ is being encompmassed"<<endl;

			//	cout << cur_paired_sam_rec.paired_sam_rec.first->tostring(0, 0)<<endl;
			//}

			//add encompass read to junction 
			if (cur_fusion_region.m_doner_st == cur_paired_sam_rec.paired_sam_rec.first->start &&
				cur_fusion_region.m_doner_end == cur_paired_sam_rec.paired_sam_rec.first->end &&
				cur_fusion_region.m_acceptor_st == cur_paired_sam_rec.paired_sam_rec.second->start &&
				cur_fusion_region.m_acceptor_end == cur_paired_sam_rec.paired_sam_rec.second->end)
			{
				//left_splice_ways.push_back(SpliceWayTrue(*v_iter));

				//left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

				cur_paired_sam_rec.paired_sam_rec.first->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, true));

				cur_paired_sam_rec.paired_sam_rec.second->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, false));
			}
			else if (cur_fusion_region.m_doner_st == cur_paired_sam_rec.paired_sam_rec.second->start &&
					 cur_fusion_region.m_doner_end == cur_paired_sam_rec.paired_sam_rec.second->end &&
					 cur_fusion_region.m_acceptor_st == cur_paired_sam_rec.paired_sam_rec.first->start &&
					 cur_fusion_region.m_acceptor_end == cur_paired_sam_rec.paired_sam_rec.first->end)
			{
				cur_paired_sam_rec.paired_sam_rec.first->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, false));

				cur_paired_sam_rec.paired_sam_rec.second->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, true));

				//SpliceWayTrue::SpliceWayTrue(string* cnp, size_t st, vector<pair<size_t, int> >* svp) : chrom_name(*cnp), start(st), spliceway_vec(*svp)

				//left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));
			}
			else
			{
				cerr <<"no matched alignments?"<<endl;
			}
			

			//pthread_mutex_lock(&fusion_lock);
		
			//cur_fusion_region_iter->m_junc_seed_ptr->m_fusion_encompassing_reads.push_back(cur_paired_sam_rec.paired_sam_rec.first->tostring(0,0));

			//cur_fusion_region_iter->m_junc_seed_ptr->m_fusion_encompassing_reads.push_back(cur_paired_sam_rec.paired_sam_rec.second->tostring(0,0));

			//pthread_mutex_unlock(&fusion_lock);
			//find_region = true;
			//break;
		}

		++cur_fusion_region_iter;
	}

	//if (find_region)
	//	find_region_idx = cur_express_region_iter - sorted_regions.begin();

	//return find_region;
}

void 
AlignmentHandler::SetHits(vector<SamRec*>& samrecs)
{
	vector<SamRec*>::iterator sam_rec_vec_iter;

	for (sam_rec_vec_iter = samrecs.begin(); sam_rec_vec_iter != samrecs.end(); ++sam_rec_vec_iter)
	{
		SamRec2Hits(*sam_rec_vec_iter);
	}
}

void
AlignmentHandler::SamRec2Hits(SamRec* sam_rec_ptr)
{
	vector<pair<size_t, int> >::iterator splice_way_iter;

	hash_map<string, vector<bool> >::iterator mapped_pos_iter;
	
	mapped_pos_iter = m_mapped_pos.find(sam_rec_ptr->chrom_name);

	if (mapped_pos_iter == m_mapped_pos.end())
		return;

	vector<bool>& mapped_pos_ref = mapped_pos_iter->second;

	//cout << sam_rec_ptr->tostring(1,1) << endl;

	for (splice_way_iter = sam_rec_ptr->spliceway_vec.begin(); splice_way_iter != sam_rec_ptr->spliceway_vec.end(); ++splice_way_iter)
	{
		if (splice_way_iter->second > 0)
		{
			for (size_t offset = splice_way_iter->first; offset < splice_way_iter->first + splice_way_iter->second; ++offset)
			{
				mapped_pos_ref[offset] = true;
			}
		}
	}
}

void
AlignmentHandler::FindSamVecContigLen(vector<SamRec*>& sam_rec_ptr)
{
	vector<SamRec*>::iterator sam_rec_iter;

	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		hash_map<string, vector<pair<size_t, size_t> > >::iterator chrom_expressed_regions_iter;

		chrom_expressed_regions_iter = m_expressed_regions.find((*sam_rec_iter)->chrom_name);

		if (chrom_expressed_regions_iter == m_expressed_regions.end())
			continue;

		vector<pair<size_t, size_t> >& sorted_regions = chrom_expressed_regions_iter->second;

		if (m_unioned_expressed_regions_map.find((*sam_rec_iter)->chrom_name) == m_unioned_expressed_regions_map.end())
			continue;

		hash_map<size_t, UnionExpressedRegions >& cur_unioned_expressed_regions = m_unioned_expressed_regions_map.find((*sam_rec_iter)->chrom_name)->second;

		pair<size_t, size_t> cur_junc_region_start = make_pair((*sam_rec_iter)->start, (*sam_rec_iter)->start);
		
		size_t start_region_idx = -1;

		bool find_start_region = FindRegion(cur_junc_region_start, sorted_regions, start_region_idx);

		DisjointSet& cur_disjointSet = (m_disjointset_map.find((*sam_rec_iter)->chrom_name))->second;

		if (find_start_region)
		{
			size_t source_idx = cur_disjointSet.FindSet(start_region_idx);

			(*sam_rec_iter)->m_contig_len = cur_unioned_expressed_regions[source_idx].m_unioned_region_length;
		}

		//if (((*sam_rec_iter)->isspliced || (*sam_rec_iter)->issmallins || (*sam_rec_iter)->is_fusion) && real_sam_count)
		//	SamRec2Junc(**sam_rec_iter, real_sam_count, cur_read_junc);
	}
}

//bool
//AlignmentHandler::FindSamContigLen(vector<SamRec*>& sam_rec_ptr)
//{
//}

void
AlignmentHandler::FilterByFilteredJunction(vector<SamRec*>& sam_rec_ptr, JunctionHandler* junc_handler)
{
	vector<SamRec*> filtered_sam_rec_ptr;

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		vector<JunctionSeed*>::iterator junc_iter;

		bool filtered = false;

		//junc_handler->FindCorrespondingJunction(*sam_rec_ptr_iter);

		//if ((*sam_rec_ptr_iter)->issmallins)
		//	filtered = true;

		//for (junc_iter = (*sam_rec_ptr_iter)->corresponding_juncs.begin(); junc_iter != (*sam_rec_ptr_iter)->corresponding_juncs.end(); ++junc_iter)
		//{
		//	if ((*junc_iter)->m_filtered_type != NOT_FILTERED)
		//	{
		//		filtered = true;		
		//		break;
		//	}
		//}

		if ((*sam_rec_ptr_iter)->filter_type != NOT_FILTERED)
			filtered = true;

		bool trimmed = false;

		if (filtered)
		{
			//if ((*sam_rec_ptr_iter)->filter_type.find(FILTERED_BY_SMALL_ANCHOR) == (*sam_rec_ptr_iter)->filter_type.end())
			//	cout << "corresponding junction:\t"<<(*junc_iter)->to_normal_junction(0)<< endl;
			//cout << "original junction alignment:\t"<<(*sam_rec_ptr_iter)->tostring(0,0) << endl;
			if (!(*sam_rec_ptr_iter)->is_fusion)
				trimmed = (*sam_rec_ptr_iter)->clip_by_small_anchor(m_add_S/**junc_iter*/);
			else
				trimmed = false;//(*sam_rec_ptr_iter)->clip_fusion_alignment(m_add_S/**junc_iter*/);

			//if (trimmed)
			//	cout << "clipped alignment:\t"<<(*sam_rec_ptr_iter)->tostring(0,0) << endl;
			//else 
			//	cout << "alignment filtered"<<endl;
		}

		if (!filtered || trimmed)
			filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);

	}

	sam_rec_ptr = filtered_sam_rec_ptr;

	if (sam_rec_ptr.size() > 1)
		RemoveDup(sam_rec_ptr);
}

void* AlignmentHandler::FilterByFilteredJunctionStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);;
	
	JunctionHandler* junc_handler = m_junction_handler_prev_ptr;

	vector<SamRec*> filtered_sam_rec_ptr;

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		vector<JunctionSeed*>::iterator junc_iter;

		bool filtered = false;

		//junc_handler->FindCorrespondingJunction(*sam_rec_ptr_iter);

		//if ((*sam_rec_ptr_iter)->issmallins)
		//	filtered = true;

		//for (junc_iter = (*sam_rec_ptr_iter)->corresponding_juncs.begin(); junc_iter != (*sam_rec_ptr_iter)->corresponding_juncs.end(); ++junc_iter)
		//{
		//	if ((*junc_iter)->m_filtered_type != NOT_FILTERED)
		//	{
		//		//cout << (*sam_rec_ptr_iter)->tostring(0, 0) << endl;

		//		//cout << "corresponding junction is filtered" << endl;
		//		//
		//		//cerr << "corresponding junction is filtered" << endl;

		//		filtered = true;

		//		break;
		//	}
		//}

		if ((*sam_rec_ptr_iter)->filter_type != NOT_FILTERED) 
		{
			//cout << "corresponding junction is filtered by small anchor" << endl;
			//cerr << "corresponding junction is filtered by small anchor" << endl;
			filtered = true;
		}

		bool trimmed = false;

		if (filtered)
		{
			//if ((*sam_rec_ptr_iter)->filter_type.find(FILTERED_BY_SMALL_ANCHOR) == (*sam_rec_ptr_iter)->filter_type.end())
			//	cout << "corresponding junction:\t"<<(*junc_iter)->to_normal_junction(0)<< endl;
			//cout << "original junction alignment:\t"<<(*sam_rec_ptr_iter)->tostring(0,0) << endl;
			if (!(*sam_rec_ptr_iter)->is_fusion)
				trimmed = (*sam_rec_ptr_iter)->clip_by_small_anchor(m_add_S/**junc_iter*/);
			else
			{
				trimmed = false;
				//(*sam_rec_ptr_iter)->clip_fusion_alignment(m_add_S/**junc_iter*/);
			}

			//if (trimmed)
			//	cout << "clipped alignment:\t"<<(*sam_rec_ptr_iter)->tostring(0,0) << endl;
			//else 
			//	cout << "alignment filtered"<<endl;
		}

		if (!filtered || trimmed)
			filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);

	}

	sam_rec_ptr = filtered_sam_rec_ptr;

	if (sam_rec_ptr.size() > 1)
		RemoveDupStatic((void*)&sam_rec_ptr);
}
void
AlignmentHandler::FilterSingleMulti(vector<SamRec*>& sam_rec_ptr)
{
	if (sam_rec_ptr.size() <= 1)
	{
		return;
	}

	if (sam_rec_ptr.size() > m_max_hits)
	{
		sam_rec_ptr.clear();
		return;
	}

	FilterByMisMatch(sam_rec_ptr);

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_filterscore);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*(sam_rec_ptr_iter - 1))->isexonic && (*(sam_rec_ptr_iter))->isexonic &&
			(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen)
			;
		else if (!(*(sam_rec_ptr_iter - 1))->isexonic && !(*(sam_rec_ptr_iter))->isexonic && 
			(*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
			(*(sam_rec_ptr_iter - 1))->ave_intron_len == (*(sam_rec_ptr_iter))->ave_intron_len &&
			(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
			(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
			&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
			)
			;
		else
		{
			sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

			break;
		}
	}
}

void* AlignmentHandler::FilterSingleMultiStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	if (sam_rec_ptr.size() == 1)
	{
		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;
	}

	if (sam_rec_ptr.size() <= 1)
	{
		return 0;
	}

	if (sam_rec_ptr.size() > m_max_hits)
	{
		sam_rec_ptr.clear();
		return 0;
	}

	bool is_fusion = false, is_not_fusion = true;

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		//if ((*sam_rec_ptr_iter)->tag_name.find("R_49939/") != string::npos )
		//{
		//	cout << (*sam_rec_ptr_iter)->tostring(1,1) << endl ;

		//	cout <<"filter_score" << endl<< (*sam_rec_ptr_iter)->filter_score << endl;

		//	cout <<"mappedlen" << endl<< (*sam_rec_ptr_iter)->mappedlen << endl;

		//	cout <<"ave_junc_mis" << endl<< (*sam_rec_ptr_iter)->ave_junc_mis << endl;

		//	cout <<"m_contig_len" << endl<< (*sam_rec_ptr_iter)->m_contig_len << endl;
		//}
		if ((*sam_rec_ptr_iter)->is_fusion == true)
			is_fusion = true;
		else
			is_not_fusion = false;
	}

	if (is_fusion && is_not_fusion == false)
	{
		FilterByMisMatchStatic((void*)&sam_rec_ptr);

		sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_fusion_filterscore);

		vector<SamRec*>::iterator sam_rec_ptr_iter;

		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;

		if (!(m_do_filter & 16))
		{
			return 0;
		}

		for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
		{
			if ((*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
				(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
				&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
				)
				;
			else
			{
				sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

				break;
			}
		}

		return 0;
	}
	else if (is_fusion && is_not_fusion)
	{
		FilterByMisMatchStatic((void*)&sam_rec_ptr);

		sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_fusion_filterscore);

		vector<SamRec*>::iterator sam_rec_ptr_iter;

		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;

		if (!(m_do_filter & 16))
		{
			return 0;
		}

		for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
		{
			if ((*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
				(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
				&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
				)
				;
			else
			{
				sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

				break;
			}
		}

		return 0;
	}
	else 
	{
		FilterByMisMatchStatic((void*)&sam_rec_ptr);

		sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_filterscore);

		vector<SamRec*>::iterator sam_rec_ptr_iter;

		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;

		if (!(m_do_filter & 16))
		{
			return 0;
		}

		for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
		{
			if ((*(sam_rec_ptr_iter - 1))->isexonic && (*(sam_rec_ptr_iter))->isexonic &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen)
				;
			else if (!(*(sam_rec_ptr_iter - 1))->isexonic && !(*(sam_rec_ptr_iter))->isexonic && 
				(*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
				(*(sam_rec_ptr_iter - 1))->ave_intron_len == (*(sam_rec_ptr_iter))->ave_intron_len &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
				(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
				&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
				)
				;
			else
			{
				sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

				break;
			}
		}
	}

	return 0;
}

void* AlignmentHandler::RemoveFusionUnPairedStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	if (!(m_do_filter & 16))
	{
		return 0;
	}

	bool is_fusion = false, is_not_fusion = false;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion == true)
			is_fusion = true;
		else
			is_not_fusion = false;
	}

	if (is_fusion && is_not_fusion)
	{
		cout << "is fusion and is not fusion:"<<(*sam_rec_ptr_iter)->tag_name<<endl;
	}
	else if (is_fusion)
	{
		sam_rec_ptr.clear();
	}
}

void AlignmentHandler::SetSamRecUnpaired(pair<vector<SamRec>, vector<SamRec> >& m_sam_rec_pe)
{
	vector<SamRec>::iterator sam_rec_iter;

	for (sam_rec_iter = m_sam_rec_pe.first.begin(); sam_rec_iter != m_sam_rec_pe.first.end(); ++sam_rec_iter)
	{
		sam_rec_iter->set_unmapped();
	}

	for (sam_rec_iter = m_sam_rec_pe.second.begin(); sam_rec_iter != m_sam_rec_pe.second.end(); ++sam_rec_iter)
	{
		sam_rec_iter->set_unmapped();
	}
}

void AlignmentHandler::SetSamRecUnpaired(pair<vector<SamRec*>, vector<SamRec*> >& m_sam_rec_pe_ptr)
{
	vector<SamRec*>::iterator sam_rec_iter;

	for (sam_rec_iter = m_sam_rec_pe_ptr.first.begin(); sam_rec_iter != m_sam_rec_pe_ptr.first.end(); ++sam_rec_iter)
	{
		(*sam_rec_iter)->set_unmapped();
	}

	for (sam_rec_iter = m_sam_rec_pe_ptr.second.begin(); sam_rec_iter != m_sam_rec_pe_ptr.second.end(); ++sam_rec_iter)
	{
		(*sam_rec_iter)->set_unmapped();
	}
}

bool AlignmentHandler::MergeStdFusion(vector<SamRec>::iterator sam_rec_ptr1, vector<SamRec>::iterator sam_rec_ptr2)
{
	//fusion 1
	sam_rec_ptr1->splice_way2 = sam_rec_ptr2->splice_way;

	sam_rec_ptr1->ori_splice_way2 = sam_rec_ptr2->ori_splice_way;

	sam_rec_ptr1->chrom_name2 = sam_rec_ptr2->chrom_name;

	sam_rec_ptr1->start2 = sam_rec_ptr2->start;

	sam_rec_ptr1->end2 = sam_rec_ptr2->end;

	sam_rec_ptr1->strand_t2 = sam_rec_ptr2->strand_t;

	sam_rec_ptr1->mappedlen1 = sam_rec_ptr1->mappedlen;

	sam_rec_ptr1->mappedlen2 = sam_rec_ptr2->mappedlen;

	sam_rec_ptr1->mappedlen = sam_rec_ptr1->mappedlen1 + sam_rec_ptr1->mappedlen2;

	sam_rec_ptr1->fusion_prefix_len = sam_rec_ptr1->mappedlen1;

	sam_rec_ptr1->fusion_suffix_len = sam_rec_ptr1->mappedlen2;

	//cout << "sam_rec_ptr1->fusion_prefix_len:"<<sam_rec_ptr1->fusion_prefix_len<<endl;

	//cout << "sam_rec_ptr1->fusion_suffix_len:"<<sam_rec_ptr1->fusion_suffix_len<<endl;

	sam_rec_ptr1->intron_size = sam_rec_ptr1->intron_size + sam_rec_ptr2->intron_size;

	sam_rec_ptr1->spliceway_vec2 = sam_rec_ptr2->spliceway_vec;

	sam_rec_ptr1->is_fusion = true;

	if (sam_rec_ptr1->strand_t & IS_REVERSE)
	{
		sam_rec_ptr1->strand1 = '-';
	}
	else
	{
		sam_rec_ptr1->strand1 = '+';
	}

	if (sam_rec_ptr1->strand_t2 & IS_REVERSE)
	{
		sam_rec_ptr1->strand2 = '-';
	}
	else
	{
		sam_rec_ptr1->strand2 = '+';
	}

	if (sam_rec_ptr1->strand1 == '+')
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->start;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->end;
	}
	else
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->end;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->start;
	}

	if (sam_rec_ptr1->strand2 == '+')
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->start2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->end2;
	}
	else
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->end2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->start2;
	}

	if (sam_rec_ptr1->min_anchor > sam_rec_ptr2->min_anchor)
	{
		sam_rec_ptr1->min_anchor = sam_rec_ptr2->min_anchor;
	}

	if (sam_rec_ptr1->max_anchor < sam_rec_ptr2->max_anchor)
	{
		sam_rec_ptr1->max_anchor = sam_rec_ptr2->max_anchor;
	}

	if (sam_rec_ptr1->chrom_name < sam_rec_ptr1->chrom_name2 || (sam_rec_ptr1->chrom_name == sam_rec_ptr1->chrom_name2 && sam_rec_ptr1->start <= sam_rec_ptr1->start2))
	{
	}
	else
	{
		sam_rec_ptr1->Swap();
	}
}


void AlignmentHandler::CollectStats(const vector<SamRec*>& sam_rec_ptr)
{
	if (sam_rec_ptr.empty())
		return;
	//size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_non_canonical, m_paired, m_single;

	vector<SamRec*>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;
	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((*sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((*sam_rec_iter)->isexonic)
			isexonic = true;

		if ((*sam_rec_iter)->isspliced)
			isspliced = true;

		if ((*sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((*sam_rec_iter)->issmallins)
			issmallins = true;

		if ((*sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((*sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((*sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((*sam_rec_iter)->isclipped)
			isclipped = true;

		if ((*sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}

	if (isexonic)
		++m_unspliced;

	if (isspliced)
		++m_spliced;

	if (issmalldel)
		++m_deletion;

	if (issmallins)
		++m_insertion;

	if (ispaired)
		++m_paired;
	else if (!isunmapped)
		++m_single;

	if (iscanon)
		++m_canoical;

	if (issemicanon)
		++m_semi_canonical;

	if (isnoncanon)
		++m_non_canonical;

	if (isclipped)
		++m_clipped;

	if (isfusion)
		++m_fusion;

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
		++m_multiple;
	else if (sam_rec_ptr.size() == 1)
		++m_unique;
}

void* AlignmentHandler::CollectStatsStatic(void * str)
{
	const vector<SamRec*>& sam_rec_ptr = (*(vector<SamRec*>*)str);

	if (sam_rec_ptr.empty())
		return 0;
	//size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_non_canonical, m_paired, m_single;

	vector<SamRec*>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;
	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((*sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((*sam_rec_iter)->isexonic)
			isexonic = true;

		if ((*sam_rec_iter)->isspliced)
			isspliced = true;

		if ((*sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((*sam_rec_iter)->issmallins)
			issmallins = true;

		if ((*sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((*sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((*sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((*sam_rec_iter)->isclipped)
			isclipped = true;

		if ((*sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}

	if (isexonic)
		++m_unspliced;

	if (isspliced)
		++m_spliced;

	if (issmalldel)
		++m_deletion;

	if (issmallins)
		++m_insertion;

	if (ispaired)
		++m_paired;
	else if (!isunmapped)
		++m_single;

	if (iscanon)
		++m_canoical;

	if (issemicanon)
		++m_semi_canonical;

	if (isnoncanon)
		++m_non_canonical;

	if (isclipped)
		++m_clipped;

	if (isfusion)
		++m_fusion;

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
		++m_multiple;
	else if (sam_rec_ptr.size() == 1)
		++m_unique;

	return 0;
}

void AlignmentHandler::CollectStats(const vector<SamRec>& sam_rec_ptr)
{
	if (sam_rec_ptr.empty())
		return;
	//size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_non_canonical, m_paired, m_single;

	vector<SamRec>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;
	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((sam_rec_iter)->isexonic)
			isexonic = true;

		if ((sam_rec_iter)->isspliced)
			isspliced = true;

		if ((sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((sam_rec_iter)->issmallins)
			issmallins = true;

		if ((sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((sam_rec_iter)->isclipped)
			isclipped = true;

		if ((sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}
	

	if (isexonic)
		++m_unspliced;

	if (isspliced)
		++m_spliced;

	if (issmalldel)
		++m_deletion;

	if (issmallins)
		++m_insertion;

	if (ispaired)
		++m_paired;
	else if (!isunmapped)
		++m_single;

	if (iscanon)
		++m_canoical;

	if (issemicanon)
		++m_semi_canonical;

	if (isnoncanon)
		++m_non_canonical;

	if (isclipped)
		++m_clipped;

	if (isfusion)
		++m_fusion;

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
		++m_multiple;
	else if (sam_rec_ptr.size() == 1)
		++m_unique;
}

void AlignmentHandler::WriteStats(string stat_file, string headline)
{
	cout << "write stat"<<endl;

	ofstream ofs(stat_file.c_str(), ofstream::app);

	ofs << headline << endl;

	ofs <<"unmapped\t"<<m_unmapped<<endl;

	ofs <<"unspliced\t"<<m_unspliced<<endl;

	ofs <<"spliced\t"<<m_spliced<<endl;

	ofs <<"deletion\t"<<m_deletion<<endl;

	ofs <<"insertion\t"<<m_insertion<<endl;

	ofs <<"paired\t"<<m_paired<<endl;

	ofs <<"single\t"<<m_single<<endl;

	ofs <<"canoical\t"<<m_canoical<<endl;

	ofs <<"semi_canonical\t"<<m_semi_canonical<<endl;

	ofs <<"non_canonical\t"<<m_non_canonical<<endl;

	ofs <<"multiple\t"<<m_multiple<<endl;

	ofs <<"unique\t"<<m_unique<<endl;

	ofs <<"clipped\t"<<m_clipped<<endl;

	ofs <<"fusion\t"<<m_fusion<<endl;

	//m_junction_handler.WriteStats(stat_file);
}

