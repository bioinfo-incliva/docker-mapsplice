#include "AlignmentHandler.h"
//#include "JunctionAccumulator.h"

int main(int argc, char** argv)
{
	ReadNextTagAlignHandler<FusionSamRec> test_fusion_sam;

	pair<vector<FusionSamRec>, vector<FusionSamRec> > v_fusion_rec_pe;

	string alignment_file = "fusion.remapped.multiple";

	test_fusion_sam.Init(alignment_file.c_str(), 3);

	size_t read_id;

	size_t sam_count = test_fusion_sam.ReadNextTagAlignPE(v_fusion_rec_pe, read_id);

}