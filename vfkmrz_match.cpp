#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>

#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

#include "flat_hash_map.hpp"

using namespace std;


// this program scans its input (fastq text stream) for forward k mers,

// usage:
//    g++ -O3 --std=c++11 -o vfkmrz_fastq vfkmrz_fastq.cpp
//    gzip -dc /path/exp.fastq.gz | ./vfkmrz_fastq
// 	  or 
//    zcat /path/exp.fastq.gz | ./vfkmrz_fastq
//    or 
//    cat /path/exp.fastq | ./vfkmrz_fastq
//    the output path can be specified in a variable below
//
// standard fastq format only for input, otherwise failure is almost guaranteed. 

// global variable declaration starts here
constexpr auto k = 31;
constexpr auto r_len = 90;
constexpr auto n_kmer_per_read = 3;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 64 * 1024 * 1024;
constexpr auto buffer_size = 64 * 1024 * 1024;

// maximum lines when reached; also the max memory controller
constexpr auto seg_l = 1000*1000*5;

// maximum lines when reached the program exits; for testing or practical use
constexpr auto max_l = 1000*1000*100;

constexpr auto max_load = 10*1000*1000;

// output file path
constexpr auto out_path = "/dev/stdout";

// gigantic vectorization
void kmer_search(vector<char>& kmers, const char* buf, int end_pos) {    
    for (int i = 0;  i <= end_pos - k;  ++i) {
        for (int j = i;  j < i+k;  ++j)
            kmers.push_back(buf[j]);
            
        kmers.push_back('\n');
    }
}

// get time elapsed since when it all began in milliseconds.
long chrono_time() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

void kmer_match() {	
    auto t_start = chrono_time();
	
    vector<char> buffer(buffer_size);
    char* window = buffer.data();

    uintmax_t n_lines = 0;

	int cur_pos = 0;
    int snp_cur = 0;

	char seq_buf[k];
	char snp_buf[12];

    int field = 0;
    //unordered_map<string, uintmax_t> kdb = {};
    //unordered_map<string, uintmax_t> kdbc = {};

    ska::flat_hash_map<string, uintmax_t> kdb = {};
    ska::flat_hash_map<string, uintmax_t> kdbc = {};
    //auto fh = fstream(out_path, ios::out | ios::binary);

    int fd;
    fd = open("/Users/jasonshi/Documents/zjshi_github/beta/kmerization/kmer31_db.txt", O_RDONLY);

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

                for(int j = snp_cur; j < 12; ++j) {
                    snp_buf[j] = ' ';
                }

                string str(seq_buf);
                string istr(snp_buf);
                stringstream ss(istr);
                uintmax_t snp_id = 0;
                ss >> snp_id;

/*
                for(auto it = snp_buf.begin(); it != snp_buf.end(); ++it){
                    cout << *it << "\n";    
                }
*/
                kdb.insert({str, snp_id});
                kdbc.insert({str, 0});

                cur_pos = snp_cur = field = 0;

                if (n_lines > max_load)
                    break;
            } else if (c == ' ') {
                field = 1;
            } else {
                if (field) {
                    snp_buf[snp_cur++] = c; 
                } else {
                    seq_buf[cur_pos++] = c;
                }
            }
        }

        
        //fh.write(&kmers[0], kmers.size());
        
        cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
        cerr << "the kdb size is " << kdb.size() << "\n";

        if (n_lines > max_load)
            break;
    }
    
/*
    for(auto it = kdb.begin(); it != kdb.end(); ++it){
        cout << it->first << ":" << it->second << "\n";    
    }

    return;
*/
    //fh.close();
	
    t_start = chrono_time();

    n_lines = 0;
    cur_pos = 0;

    int kmer_count = 0;
    
    close(fd);

    int fc;
    fc = open("/Users/jasonshi/Documents/zjshi_github/beta/kmerization/vfkmrz_fastq.out", O_RDONLY);

    int n_kmer_per_read = 3;
    bool has_wildcard = false;
    ska::flat_hash_map<uintmax_t, int> foot_print= {};
    
    while (true) {

        const ssize_t bytes_read = read(fc, window, step_size);
		
        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
        	cerr << "unknown fetal error!" << endl;
			exit(EXIT_FAILURE);
		}

        for (int i = 0;  i < bytes_read;  ++i) {
            char c = window[i];
            if (c == '\n') {
                ++n_lines;
                ++kmer_count;

                cur_pos = 0;

                if (kmer_count > n_kmer_per_read) {
                    kmer_count = 0;
                    foot_print.clear();
                }

                if (has_wildcard) {
                    has_wildcard = false;
                    continue;    
                }

                if (kdb.find(seq_buf) != kdb.end()){
                    if (!foot_print[kdb[seq_buf]]) {
                        ++kdbc[seq_buf];
                        foot_print[kdb[seq_buf]] = 1;    
                    } else {
                        cerr << seq_buf << "\n";    
                    }
                }
                
            } else {
                if (c == 'N') {
                    has_wildcard = true;
                }

                seq_buf[cur_pos++] = c;
            }
        }
        
        if (n_lines > 100*1000*1000)
            break;

        //fh.write(&kmers[0], kmers.size());
        cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }
    
    close(fc); 

    auto fh = fstream(out_path, ios::out | ios::binary);   

    for(auto it = kdbc.begin(); it != kdbc.end(); ++it){
        fh << it->first << "\t" << it->second << "\n";    
    }

    fh.close();

    auto t_time = (chrono_time() - t_start) / 1000;

    cerr << "done!" << endl;
    cerr << "total lines: " << n_lines << endl;
    cerr << "mappng speed: " << (double)n_lines / (double)t_time << " lines/sec" << endl;
}

int main(int argc, char** argv){		
    kmer_match();
    return 0;
}
