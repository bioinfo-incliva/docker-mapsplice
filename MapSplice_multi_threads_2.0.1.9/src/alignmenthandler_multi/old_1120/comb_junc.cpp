#include "JunctionHandler.h"

int main(int argc, char** argv)
{
	vector<string> junction_file1;

	for (int i = 1; i < argc - 1; ++i)
		junction_file1.push_back(argv[i]);

	JunctionHandler junc_handler1, junc_handler2;

	junc_handler1.ReadJunction(junction_file1);

	string junction_file = argv[argc - 1];

	string junction_ins_file = argv[argc - 1];  junction_ins_file.append(".ins");

	string junction_del_file = argv[argc - 1]; junction_del_file.append(".del");

	string junction_fusion_file = argv[argc - 1]; junction_fusion_file.append(".fusion");

	string junction_filtered_file = argv[argc - 1]; junction_filtered_file.append(".filtered");

	junc_handler1.WriteJunction(junction_file, junction_ins_file, junction_del_file, junction_fusion_file, junction_filtered_file);

}