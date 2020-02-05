#include <Rcpp.h>
#include <fstream>
#include <api/BamReader.h>
using namespace Rcpp;
using namespace std;
using namespace BamTools;

int sort(string clipped_filename, unsigned threads = 2){
	string cmd_pre = "sort -g -k1 -k2 " + clipped_filename;
	string cmd_post = " | uniq > " + clipped_filename + "-sorted";
	string cmd_sort = cmd_pre + cmd_post;
	string cmd_sort_parallel = cmd_pre + " --parallel=" + to_string(threads) + cmd_post;

  cout << "Sorting " << clipped_filename << "\n";
	int value = system(cmd_sort_parallel.c_str());
	if (value != 0) {
		cout << "Sorting with single thread" << "\n";
		value = system(cmd_sort.c_str());
	}

	return value;
}
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector get_refName(string bam_filename) {
	BamReader reader;
	if ( !reader.Open(bam_filename)) {
		cerr << "Could not open input BAM file\n";
	}

	auto& refData = reader.GetReferenceData();
	int n = refData.size();
	NumericVector refID(n);
  CharacterVector refName(n);

	for (int i = 0; i < n; i++) {
		refID[i] = i;
		refName[i] = refData[i].RefName;
	}
	refID.names() = refName;
	return refID;
}

// [[Rcpp::export]]
void clip(string bam_filename, string out_prefix, int minimum_read_length = 5, unsigned threads = 2) {
	cout << "Clipping Reads\n";

	BamReader reader;
	if ( !reader.Open(bam_filename)) {
		cerr << "Could not open input BAM file\n";
		//return false;
	}

	BamAlignment al;
	string f_clipped_filename = out_prefix + "-f";
	string r_clipped_filename = out_prefix + "-r";
	ofstream f_clipped(f_clipped_filename);
	cout << "Writing " << f_clipped_filename << "\n";
	ofstream r_clipped(r_clipped_filename);
	cout << "Writing " << r_clipped_filename << "\n";

	while(reader.GetNextAlignment(al)) {
		if (al.IsDuplicate()) {
			continue;
		}

		string refID = to_string(al.RefID);

		auto cigars = al.CigarData;

		if (cigars.size() == 0) {
			continue;
		}

		auto& cigars_front_type = cigars.front().Type;
		auto& cigars_back_type = cigars.back().Type;

		int cigars_front_length = cigars.front().Length;
		int cigars_back_length = cigars.back().Length;

		if (('S' == cigars_front_type || 'H' == cigars_front_type)
				&& cigars_front_length >= minimum_read_length) {

			if ('H' == cigars_front_type){
				cigars_front_length = 0;
			}
			if ('H' == cigars_back_type) {
				cigars_back_length = 0;
			}

			int clipped_pos = al.Position;
			string clipped_seq;
			clipped_seq = al.QueryBases.substr(0, cigars_front_length);

			f_clipped << refID << "\t" << clipped_pos << "\t" << clipped_seq << "\n";
		}

		if (('S' == cigars_back_type || 'H' == cigars_back_type)
				&& cigars_back_length >= minimum_read_length) {

			if ('S' != cigars_front_type){
				cigars_front_length = 0;
			}
			if ('S' != cigars_back_type) {
				cigars_back_length = 0;
			}

			int delta = 0;
			int n_del = 0;
			int n_in = 0;

			for (auto cigar : al.CigarData) {
				if ('D' == cigar.Type || 'N' == cigar.Type) {
					n_del += cigar.Length;
				} else if ('I' == cigar.Type) {
					n_in += cigar.Length;
				}
			}
			delta = n_del - n_in;

			int clipped_pos = al.Position + al.QueryBases.size() - cigars_front_length - cigars_back_length + delta + 1;
			string clipped_seq;
			clipped_seq = al.QueryBases.substr(al.QueryBases.size() - cigars_back_length);

			r_clipped << refID << "\t" << clipped_pos << "\t" << clipped_seq << "\n";
		}
	}


	cout << "Closing " << f_clipped_filename << "\n";
	f_clipped.close();

	cout << "Closing " << r_clipped_filename << "\n";
	r_clipped.close();

	sort(f_clipped_filename, threads);
	sort(r_clipped_filename, threads);
}
