#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string.h>
using namespace std;

void
readchrom(const char* filename, string& longseq)
{
	cout << " read chrom: " << filename << endl;
	size_t size;  

	ifstream longfile(filename);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	if (longfile.is_open())
	{
		string skipline;
		getline(longfile,skipline);

		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);

			if (line.empty())
				continue;
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";

	cout <<"chrom size:"<< longseq.size() << endl;
}


char
complement(int i) {
	static const int b2c_size = 20;
	static const char b2c[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
	};
	static const char b2cl[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
	};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

struct exon {

	size_t start, end;

	string chr;

	string gene;

	string transcript;

	int strand;

	exon() {}

	exon(size_t _start, size_t _end, string _chr, string _gene, string _transcript, int _strand) : start(_start), end(_end), chr(_chr), gene(_gene), transcript(_transcript), strand(_strand)
	{
	}

	string tostring(int id)
	{
		char buf[10000];

		sprintf(buf, "%s\t%d\t%s\t%s\tE%d\t%llu\t%llu", chr.c_str(), strand, gene.c_str(), transcript.c_str(), id, start, end);

		return buf;
	}
};

bool comp_exon(const exon& lhs, const exon& rhs)
{
	if (lhs.chr == rhs.chr)
	{
		if (lhs.start == rhs.start)
		{
			return lhs.end < rhs.end;
		}
		else
		{
			return lhs.start < rhs.start;
		}
	}
	else
	{
		return lhs.chr < rhs.chr;
	}
}

struct transcript {

	string transcript_id;

	char transcript_strand;

	map<int, exon> exons;

	vector<pair<int, int> > exon_region;

	transcript (string& _transcript_id, char _transcript_strand) : transcript_id(_transcript_id), transcript_strand(_transcript_strand)
	{
	}
};

struct gene {

	string chrname;

	string geneid;

	string protein_type;

	char gene_strand;

	map<string, transcript> transcripts;

	
	gene (string _chrname, string _geneid, string _protein_type, char _gene_strand) : chrname(_chrname), geneid(_geneid), protein_type(_protein_type), gene_strand(_gene_strand)
	{
	}
};


void read_gtf_file(string gtf_file, string genetabout)
{
	ifstream ifs_gtf(gtf_file.c_str());

	ofstream ofs(genetabout.c_str());

	vector<exon> exons;

	if (ifs_gtf.is_open())
	{
		while(!ifs_gtf.eof())
		{
			string line;

			getline(ifs_gtf, line);

			if (line.empty())
				continue;

			char chr[1000], protein_type[1000], exontype[1000];

			size_t start, end;

			char skip1[1000], strand, skip2[1000];

			char gene_id[1000], geneidname[1000]; 

			char transcript_id[1000], transcriptidname[1000]; 

			char exon_id[1000], exonidname[1000]; 

			char gname[1000], gene_name[1000]; 

			char tname[1000], transcript_name[1000]; 

			char protein_id[1000], protein_id_name[1000];

			//2       protein_coding  CDS     74900584        74900607        .       +       0        
			//gene_id "ENSG00000135622"; transcript_id "ENST00000420077"; exon_number "5"; 
			//gene_name "SEMA4F"; transcript_name "SEMA4F-002"; protein_id "ENSP00000416490";

			sscanf(line.c_str(), "%s\t%s\t%s\t%llu\t%llu\t%s\t%c\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chr, protein_type, exontype,
				&start, &end, skip1, &strand, skip2, gene_id, geneidname, transcript_id, transcriptidname, exon_id, exonidname,
				gname, gene_name, tname, transcript_name, protein_id, protein_id_name);

			if (string(gname) == "gene_biotype")
				continue;

			if (string(exontype) != "exon")
				continue;

			string protein_type_str = protein_type;

			if ((protein_type_str == "protein_coding" || 
				protein_type_str == "processed_transcript" || 
				protein_type_str == "IG_C_gene" || 
				protein_type_str == "IG_D_gene" || 
				protein_type_str == "IG_J_gene" || 
				protein_type_str == "IG_V_gene"))
				;
			else
				continue;

			string genename = gene_name; 
			
			genename = genename.substr(1, genename.length() - 3);

			string transcriptname = transcript_name;

			transcriptname = transcriptname.substr(1, transcriptname.length() - 3);

			int strandint = 0;

			if (strand == '+')
				strandint = 1;

			exons.push_back(exon(start, end, chr, genename, transcriptname, strandint));

			//12	protein_coding	exon	124114716	124114832	.	-	.	 gene_id "ENSG00000111361"; transcript_id "ENST00000228958"; exon_number "4"; gene_name "EIF2B1"; transcript_name "EIF2B1-003";

		}
	}
	else
	{
		cout <<"can't open:" << gtf_file<<endl;
	}

	cout << "sort"<<endl;

	sort(exons.begin(), exons.end(), comp_exon);

	vector<exon>::iterator exon_iter;

	int id = 0;

	cout << "write"<<endl;

	for (exon_iter = exons.begin(); exon_iter != exons.end(); ++exon_iter)
	{
		ofs << exon_iter->tostring(++id) << endl;
	}
}

int main(int argc, char** argv)
{
	cout << "gtf genetab_out"<<endl;

	read_gtf_file(argv[1], argv[2]);

}