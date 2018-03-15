#include <iostream>
#include <fstream>
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

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;

// maximum lines when reached; also the max memory controller
constexpr auto seg_l = 1000*1000*10;

// maximum lines when reached the program exits; for testing or practical use
constexpr auto max_l = 1000*1000*100;

// output file path
constexpr auto out_path = "./vfkmrz_union.out";

// get time elapsed since when it all began in milliseconds.
long chrono_time() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

void vfkmrz_union(const char* k1_path, const char* k2_path) {	
    auto t_start = chrono_time();
	
    vector<char> buffer(buffer_size);
    char* window = buffer.data();

    uintmax_t n_lines = 0;

	int cur_pos = 0;
	char seq_buf[k];
	
    //auto fh = fstream(out_path, ios::out | ios::binary);

    int fd;
    fd = open(k1_path, O_RDONLY);

    vector<string> kdb;
    
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
                string str(seq_buf);
                kdb.push_back(str);
                cur_pos = 0;
            } else {
                seq_buf[cur_pos++] = c;
            }
        }
        
        //fh.write(&kmers[0], kmers.size());
        
        cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }

    close(fd);

    auto timeit = chrono_time();
    cerr << "start in-place sorting and dereplicating for the first kmer list." << endl;
    cerr << "the first kmer list has " << kdb.size() << " kmers" << endl;
    sort(kdb.begin(), kdb.end());
    vector<string>::iterator ip = unique(kdb.begin(), kdb.end());
    kdb.resize(std::distance(kdb.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the first kmer list has " << kdb.size() << " unique kmers" << endl;


    //fh.close();
    n_lines = cur_pos = 0;

    int fc;
    fc = open(k2_path, O_RDONLY);

    vector<string> queries;

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
                cur_pos = 0;

                string str(seq_buf);

                queries.push_back(str);
            } else {
                seq_buf[cur_pos++] = c;
            }
        }
        
        if (n_lines > 100*1000*1000)
            break;
        //fh.write(&kmers[0], kmers.size());
        cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }
    
    close(fc); 

/*
    for(auto it = queries.begin(); it != queries.end(); ++it){
        cout << *it << "\n";    
    }
*/
    timeit = chrono_time();
    cerr << "start in-place sorting and dereplicating." << endl;
    cerr << "the second kmer list has " << queries.size() << " kmers" << endl;
    sort(queries.begin(), queries.end());
    unique(queries.begin(), queries.end());
    ip = unique(queries.begin(), queries.end());
    queries.resize(std::distance(queries.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the second kmer list has " << queries.size() << " unique kmers" << endl;

    vector<string> kmer_union(kdb.size()+queries.size());

    cerr << "start merging two vectors." << endl;
    set_union(kdb.begin(), kdb.end(), queries.begin(), queries.end(), kmer_union.begin());
    unique(kmer_union.begin(), kmer_union.end());
    ip = unique(kmer_union.begin(), kmer_union.end());
    kmer_union.resize(std::distance(kmer_union.begin(), ip));
    cerr << "Done!\n" << "It takes " << (chrono_time() - timeit) / 1000 << " secs for merging" << endl;
    cerr << "the kmer union has " << kmer_union.size() << " unique kmers\n";

    auto t_time = (chrono_time() - t_start) / 1000;

    cerr << "done!" << endl;
    cerr << "total lines: " << n_lines << endl;
    cerr << "mapping speed: " << (double)n_lines / (double)t_time << " lines/sec" << endl;
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

            cout << "arg " << i << ": " << argv[i] << endl;
		}
		
		if (k1_i == -1 || k2_i == -1) {
			cerr << argv[0] << "takes exactly two arguments (-k1 and -k2)!" << endl;
			exit(EXIT_FAILURE);
		} else if (k1_i > 4 or k2_i > 4) {
			cerr << argv[0] << "takes exactly two arguments (-k1 and -k2)!" << endl;
			exit(EXIT_FAILURE);
		} else {
	        vfkmrz_union(argv[k1_i], argv[k2_i]);		
		}
    } else {
        cerr << argv[0] << "takes exactly two arguments (-k1 and -k2)!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
