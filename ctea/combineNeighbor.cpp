#include <Rcpp.h>
#include <fstream>
#include <boost/lexical_cast.hpp>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
void combine(string input_filename, string output_filename, char f_or_r, int64_t minimum_base_gap = 2) {

	ifstream in(input_filename, ifstream::binary);

	string line;

	string prev_refID = "";
	int64_t prev_pos = 0;

	string refID;
	int64_t pos;
	string family;
	string more_family;
	string base;

	string pos_combined = "";
	int64_t pos_longbase = 0;
	string family_longbase = "";
	string more_family_longbase = "";
	string base_combined = "";
	string base_longest = "";
	int64_t count = 0;
	char clv;

	ofstream combined(output_filename);
	cout << "Writing " << output_filename << "\n";
	combined << "refID\tpos\tori\tcnt\tfamily\tlongclip\tmoreFamily\tseq\tmorePos\tmoreSeq\n";

	while(getline(in, line)){

		stringstream linestream(line);
		string value;

		if(!line.empty()){
			getline(linestream, value, '\t');
			refID = value;
			getline(linestream, value, '\t');
			pos = boost::lexical_cast<int64_t>(value);
			getline(linestream, value, '\t');
			family = value;
			getline(linestream, value, '\t');
			more_family = value;
			getline(linestream, value, '\t');
			base = value;

			if (refID == prev_refID
											&& abs(pos - prev_pos) <= minimum_base_gap) {
				if (base.size() > base_longest.size()) {
					base_longest = base;
					pos_longbase = pos;
					family_longbase = family;
					more_family_longbase = more_family;
				}
				pos_combined += "," + to_string(pos);
				base_combined += "," + base;
			}
			else {
				if (base_longest != "") {
					if (base_longest.size() >= 19) {
						clv = '2';
					} else if (base_longest.size() >= 12) {
						clv = '1';
					} else {
						clv = '0';
					}
				  combined << refID << "\t" << pos_longbase << "\t" << f_or_r << "\t" << count << "\t";
					combined << family_longbase << "\t" << clv << "\t" << more_family_longbase << "\t";
					combined << base_longest << "\t" << pos_combined << "\t" << base_combined << "\n";
				}

				count = 0;
				base_longest = base;
				pos_longbase = pos;
				family_longbase = family;
				pos_combined = to_string(pos);
				base_combined = base;
			}

			count++;
			prev_refID = refID;
			prev_pos = pos;
		}
	}

	combined.close();
	cout << "Closing " << output_filename << "\n";

}
