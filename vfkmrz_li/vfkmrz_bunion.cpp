#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>

#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

using namespace std;


// this program scans its input (fastq text stream) for forward k mers,

// usage:
//    g++ -O3 --std=c++11 -o vfkmrz_bunion vfkmrz_bunion.cpp
//    ./vfkmrz_bunion -k1 </path/to/kmer_list1> -k2 </path/to/kmer_list2>
//
// standard fastq format only for input, otherwise failure is almost guaranteed. 

// global variable declaration starts here
constexpr auto k = 31;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;

// output file path
constexpr auto out_path = "/dev/stdout";

// get time elapsed since when it all began in milliseconds.
long chrono_time() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// number of bits per single nucleotide base
constexpr int bpb = 2;

template <class int_type>
int_type bit_encode(const char c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    }

    assert(false);
}


template <class int_type>
char bit_decode(const int_type bit_code) {
    switch (bit_code) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    }
    assert(false);
}

template <class int_type>
void make_code_dict(int_type* code_dict) {
    code_dict['A'] = bit_encode<int_type>('A');
    code_dict['C'] = bit_encode<int_type>('C');
    code_dict['G'] = bit_encode<int_type>('G');
    code_dict['T'] = bit_encode<int_type>('T');
}

template <class int_type>
int_type seq_encode(const char* buf, int len, const int_type* code_dict, const int_type b_mask) {
    int_type seq_code = 0;
    for (int i=0;  i < len;  ++i) {
        const int_type b_code = code_dict[buf[i]];
        seq_code |= ((b_code & b_mask) << (bpb * (len - i - 1)));
    }
    return seq_code;
}

template <class int_type>
void seq_decode(char* buf, const int len, const int_type seq_code, int_type* code_dict, const int_type b_mask) {
    for (int i=0;  i < len;  ++i) {
        const int_type b_code = (seq_code >> (bpb * (len - i - 1))) & b_mask;
        buf[i] = bit_decode<int_type>(b_code);
    }
}


template <class int_type>
void bit_load(const char* k_path, vector<char>& buffer, vector<int_type>& k_vec, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int cur_pos = 0;
    char seq_buf[k];

    //auto fh = fstream(out_path, ios::out | ios::binary);

    int fd;
    fd = open(k_path, O_RDONLY);

    bool has_wildcard = false;

    while (true) {

        const ssize_t bytes_read = read(fd, window, step_size);

        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
            cerr << "unknown fetal error!" << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0;  i < bytes_read;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
                ++n_lines;
                cur_pos = 0;

                if (has_wildcard) {
                    has_wildcard = false;
                    continue;    
                }

                auto code = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);
                k_vec.push_back(code);
            } else {
                if (c == 'N') {
                    has_wildcard = true;    
                }

                seq_buf[cur_pos++] = c;
            }
        }

        //fh.write(&kmers[0], kmers.size());

        cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }

    close(fd);

    auto timeit = chrono_time();
    cerr << "start in-place sorting and dereplicating for the first kmer list." << endl;
    cerr << "the first kmer list has " << k_vec.size() << " kmers" << endl;
}

template <class int_type>
void vfkmrz_bunion(const char* k1_path, const char* k2_path) {	
    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

    vector<int_type> kdb;
    vector<char> buffer(buffer_size);

    bit_load<int_type>(k1_path, buffer, kdb, code_dict, b_mask);	

    auto timeit = chrono_time();
    sort(kdb.begin(), kdb.end());
    typename vector<int_type>::iterator ip = unique(kdb.begin(), kdb.end());
    kdb.resize(std::distance(kdb.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the first kmer list has " << kdb.size() << " unique kmers" << endl;


    vector<int_type> kqr;
    bit_load<int_type>(k2_path, buffer, kqr, code_dict, b_mask);	
    /*
       for(auto it = kqr.begin(); it != kqr.end(); ++it){
       cout << *it << "\n";    
       }
     */
    timeit = chrono_time();
    sort(kqr.begin(), kqr.end());
    unique(kqr.begin(), kqr.end());
    ip = unique(kqr.begin(), kqr.end());
    kqr.resize(std::distance(kqr.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the second kmer list has " << kqr.size() << " unique kmers" << endl;

    vector<int_type> kmer_union(kdb.size()+kqr.size());

    cerr << "start merging two vectors." << endl;
    set_union(kdb.begin(), kdb.end(), kqr.begin(), kqr.end(), kmer_union.begin());
    unique(kmer_union.begin(), kmer_union.end());
    ip = unique(kmer_union.begin(), kmer_union.end());
    kmer_union.resize(std::distance(kmer_union.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs for merging" << endl;
    cerr << "the kmer union has " << kmer_union.size() << " unique kmers\n";

    char seq_buf[k];
    ofstream fh(out_path, ofstream::out | ofstream::binary);

    for (ip = kmer_union.begin(); ip != kmer_union.end(); ++ip) {
        seq_decode(seq_buf, k, *ip, code_dict, b_mask);    
        fh << seq_buf << "\n";
    }

    fh.close();
}

int main(int argc, char** argv){		

    int k1_i = -1, k2_i = -1;	

    if (argc == 5) {
        for(int i = 0; i < argc; ++i) {
            if (!strcmp(argv[i], "-k1")) {
                k1_i = i + 1;
            }

            if (!strcmp(argv[i], "-k2")) {
                k2_i = i + 1;
            }

            cerr << "arg " << i << ": " << argv[i] << endl;
        }

        if (k1_i == -1 || k2_i == -1) {
            cerr << argv[0] << "takes exactly two arguments (-k1 and -k2)!" << endl;
            exit(EXIT_FAILURE);
        } else if (k1_i > 4 or k2_i > 4) {
            cerr << argv[0] << "takes exactly two arguments (-k1 and -k2)!" << endl;
            exit(EXIT_FAILURE);
        } else {
            if (1 <= k && k <=4){
                vfkmrz_bunion<uint_fast8_t>(argv[k1_i], argv[k2_i]);		
            } else if (5 <= k && k <= 8) {
                vfkmrz_bunion<uint_fast16_t>(argv[k1_i], argv[k2_i]);		
            } else if (8 <= k && k <= 16) {
                vfkmrz_bunion<uint_fast32_t>(argv[k1_i], argv[k2_i]);		
            } else if (17 <= k && k <= 32) {
                vfkmrz_bunion<uint_fast64_t>(argv[k1_i], argv[k2_i]);		
            } else if (33 <= k && k <= 128) {
                vfkmrz_bunion<uint64_t>(argv[k1_i], argv[k2_i]);		
            } else {
                exit(EXIT_FAILURE);
            }
        }
    } else {
        cerr << argv[0] << "takes exactly two arguments (-k1 and -k2)!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
