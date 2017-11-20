#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>

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
constexpr auto offset = 31;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 128 * 1024 * 1024;
constexpr auto buffer_size = 128 * 1024 * 1024;

// maximum lines when reached; also the max memory controller
constexpr auto seg_l = 1000*1000*5;

// maximum lines when reached the program exits; for testing or practical use
constexpr auto max_l = 1000*1000*100;

// output file path
//constexpr auto out_path = "./vfkmrz_fastq.out";
constexpr auto out_path = "/dev/stdout";

// gigantic vectorization
void kmer_search(vector<char>& kmers, const char* buf, int end_pos, int ofs) {    
    int i = 0;

    for (;  i <= end_pos - k;  i+=ofs) {
        for (int j = i;  j < i+k;  ++j)
            kmers.push_back(buf[j]);
            
        kmers.push_back('\n');
    }

    if (i < end_pos) {
        for (int j = end_pos - k;  j < end_pos;  ++j)
            kmers.push_back(buf[j]);
           
        kmers.push_back('\n');
    }
}

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

void vfkmrz_fastq() {	
    auto t_start = chrono_time();
	
    vector<char> buffer(buffer_size);
    char* window = buffer.data();

    uintmax_t n_lines = 0, n_reads = 0, n_bases = 0;

    int l_label = 2;
    uintmax_t n_lines_local = 0;
	
	int cur_pos = 0;
	char seq_buf[r_len];
	
	vector<char> kmers;

    auto fh = fstream(out_path, ios::out | ios::binary);

    while (true) {

        const ssize_t bytes_read = read(fileno(stdin), window, step_size);
		
        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
        	cerr << "unknown fetal error!" << endl;
			exit(EXIT_FAILURE);
		}

        for (int i = 0;  i < bytes_read;  ++i) {
            //char c = toupper(window[i]);
            char c = window[i];
            if (c == '\n') {
            	++n_lines;
                ++n_lines_local;
                
                if (l_label == 3) {
                	// assert (cur_pos == r_len);
                	++n_reads; 
                	kmer_search(kmers, seq_buf, cur_pos, offset);       
                	        	
                	l_label = cur_pos = 0;
                } else
                	++l_label;
                
                if (n_lines > max_l)
                	break; 
            } else {
            	if (l_label == 3) {
            		seq_buf[cur_pos++] = c;
            		++n_bases;
            	}
            }
        }
        
        if (n_lines > max_l) 
        	break;
              
        if (n_lines_local > seg_l){
        	/* iteration is slow!!
    		for (auto it = kmers.begin();  it != kmers.end(); ++it) {
        		//fh<<*it<<'\n';
        		//cout<<*it<<'\n';
    		}
    		*/
    		
    		fh.write(&kmers[0], kmers.size());
    		
        	kmers.clear();
    		n_lines_local = 0;
        	
        	cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    	}
    }
    
    if (kmers.size() != 0) {
    	fh.write(&kmers[0], kmers.size());
    		
    	kmers.clear();
    	n_lines_local = 0;
        	
        cerr << n_lines << " lines were scanned after " << (chrono_time() - t_start) / 1000 << " seconds" << endl;
    }
    
    fh.close();
    
    cerr << "done!" << endl;
    cerr << "total lines: " << n_lines << endl;
    cerr << "total reads: " << n_reads << endl;
	cerr << "total bases: " << n_bases << endl;
}

int main(int argc, char** argv){		
    vfkmrz_fastq();
    return 0;
}
