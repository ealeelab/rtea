from libcpp.string cimport string
from libcppfstream cimport ifstream, ofstream
from libcpp cimport bool
from cython.operator cimport dereference as deref
from cython.operator cimport postincrement as postinc
from libcpp.vector cimport vector
from libcpp cimport abs

# Function to combine the entries
cpdef void combine(string input_filename, string output_filename, char f_or_r, int64_t minimum_base_gap=2):
    cdef ifstream in(input_filename, ifstream.binary)
    cdef ofstream combined(output_filename)
    cdef string line
    cdef string prev_refID = ""
    cdef int64_t prev_pos = 0
    cdef string refID
    cdef int64_t pos
    cdef string family
    cdef string more_family
    cdef string base
    cdef string pos_combined = ""
    cdef int64_t pos_longbase = 0
    cdef string family_longbase = ""
    cdef string more_family_longbase = ""
    cdef string base_combined = ""
    cdef string base_longest = ""
    cdef int64_t count = 0
    cdef char clv

    print("Writing " + output_filename)
    combined << "refID\tpos\tori\tcnt\tfamily\tlongclip\tmoreFamily\tseq\tmorePos\tmoreSeq\n"

    while getline(in, line):
        cdef vector[string] values = split_line(line, '\t')
        refID = deref(values.begin())
        pos = stoi(deref(values.begin() + 1))
        family = deref(values.begin() + 2)
        more_family = deref(values.begin() + 3)
        base = deref(values.begin() + 4)

        if refID == prev_refID and abs(pos - prev_pos) <= minimum_base_gap:
            if base.size() > base_longest.size():
                base_longest = base
                pos_longbase = pos
                family_longbase = family
                more_family_longbase = more_family
            pos_combined += "," + str(pos)
            base_combined += "," + base
        else:
            if base_longest:
                if base_longest.size() >= 19:
                    clv = '2'
                elif base_longest.size() >= 12:
                    clv = '1'
                else:
                    clv = '0'
                combined << prev_refID << "\t" << pos_longbase << "\t" << f_or_r << "\t" << count << "\t"
                combined << family_longbase << "\t" << clv << "\t" << more_family_longbase << "\t"
                combined << base_longest << "\t" << pos_combined << "\t" << base_combined << "\n"
            count = 0
            base_longest = base
            pos_longbase = pos
            family_longbase = family
            pos_combined = str(pos)
            base_combined = base
        count += 1
        prev_refID = refID
        prev_pos = pos

    combined.close()
    print("Closing " + output_filename)
