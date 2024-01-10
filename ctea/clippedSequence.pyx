# clippedSequence.pyx
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdlib cimport system
import numpy as np
cimport numpy as np

# Define the sort function
cpdef int sort(string clipped_filename, unsigned threads=2):
    cdef string cmd_pre = "sort -g -k1 -k2 " + clipped_filename
    cdef string cmd_post = " | uniq > " + clipped_filename + "-sorted"
    cdef string cmd_sort = cmd_pre + cmd_post
    cdef string cmd_sort_parallel = cmd_pre + " --parallel=" + str(threads) + cmd_post

    print("Sorting " + clipped_filename)
    cdef int value = system(cmd_sort_parallel.c_str())
    if value != 0:
        print("Sorting with a single thread")
        value = system(cmd_sort.c_str())

    return value

# Define the get_refName function
cpdef np.ndarray[np.int32_t, ndim=1] get_refName(string bam_filename):
    from bamtools.api cimport BamReader
    cdef BamReader reader
    if not reader.Open(bam_filename):
        print("Could not open input BAM file")

    cdef vector[string] refName
    for refData in reader.GetReferenceData():
        refName.push_back(refData.RefName)

    cdef int n = refName.size()
    cdef np.ndarray[np.int32_t, ndim=1] refID = np.arange(n, dtype=np.int32)
    refID.names = refName
    return refID

# Define the clip function
cpdef void clip(string bam_filename, string out_prefix, int minimum_read_length=5, unsigned threads=2):
    print("Clipping Reads")

    from bamtools.api cimport BamReader
    from libcpp.vector cimport vector
    cdef BamReader reader
    if not reader.Open(bam_filename):
        print("Could not open input BAM file")

    cdef BamAlignment al
    cdef string f_clipped_filename = out_prefix + "-f"
    cdef string r_clipped_filename = out_prefix + "-r"
    cdef ofstream f_clipped(f_clipped_filename)
    print("Writing " + f_clipped_filename)
    cdef ofstream r_clipped(r_clipped_filename)
    print("Writing " + r_clipped_filename)

    while reader.GetNextAlignment(al):
        if al.IsDuplicate():
            continue

        cdef string refID = str(al.RefID)
        cdef int cigars_front_length, cigars_back_length
        cdef char cigars_front_type, cigars_back_type

        if al.CigarData.size() == 0:
            continue

        cigars_front_type = al.CigarData.front().Type
        cigars_back_type = al.CigarData.back().Type
        cigars_front_length = al.CigarData.front().Length
        cigars_back_length = al.CigarData.back().Length

        if (cigars_front_type == 'S' or cigars_front_type == 'H') and cigars_front_length >= minimum_read_length:
            if cigars_front_type == 'H':
                cigars_front_length = 0
            if cigars_back_type == 'H':
                cigars_back_length = 0

            cdef int clipped_pos = al.Position
            cdef string clipped_seq
            clipped_seq = al.QueryBases[:cigars_front_length]

            f_clipped << refID << "\t" << clipped_pos << "\t" << clipped_seq << "\n"

        if (cigars_back_type == 'S' or cigars_back_type == 'H') and cigars_back_length >= minimum_read_length:
            if cigars_front_type != 'S':
                cigars_front_length = 0
            if cigars_back_type != 'S':
                cigars_back_length = 0

            cdef int delta = 0
            cdef int n_del = 0
            cdef int n_in = 0

            for cigar in al.CigarData:
                if cigar.Type == 'D' or cigar.Type == 'N':
                    n_del += cigar.Length
                elif cigar.Type == 'I':
                    n_in += cigar.Length
            delta = n_del - n_in

            cdef int clipped_pos = al.Position + len(al.QueryBases) - cigars_front_length - cigars_back_length + delta + 1
            cdef string clipped_seq
            clipped_seq = al.QueryBases[-cigars_back_length:]

            r_clipped << refID << "\t" << clipped_pos << "\t" << clipped_seq << "\n"

    print("Closing " + f_clipped_filename)
    f_clipped.close()

    print("Closing " + r_clipped_filename)
    r_cl
