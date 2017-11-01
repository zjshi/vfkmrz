#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>

#include <unistd.h>
#include <assert.h>

using namespace std;


// this program scans its input (fasta format) for forward k mers,
//
// usage:
//    g++ -O3 --std=c++11 -o vfkmrz_fasta vfkmrz_fasta.cpp
//    cat /path/exp.fasta | ./vfkmrz_fasta
//
// standard fasta format only for input, otherwise failure is almost guaranteed. 

// global variable declaration starts here

// size of kmer
constexpr auto k = 7;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 128 * 1024 * 1024;
constexpr auto buffer_size = 128 * 1024 * 1024 + k - 1;

// maximum lines when reached; also the max memory controller
constexpr auto seg_l = 1000*1000*5;

// maximum lines when reached the program exits; for testing or practical use
constexpr auto max_l = 1000*1000*80;

// output file path
constexpr auto out_path = "./vfkmrz_fasta.out";

// toggle to organize kmer in the output file by each sequence or not; default false (pure kmers)
// the option is need here because fasta file can have sequences of different lengths
constexpr auto organized_by_seq = false;

// global variable declaration ends here

// gigantic vectorization based kmer search
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

void vfkmrz_fasta() {	
    auto t_start = chrono_time(); // time when program starts
	
    vector<char> buffer(buffer_size); // buffer zone for both file reading and sequence windowing 
    char* window = buffer.data();
	
    uintmax_t n_lines = 0, n_reads = 0, n_bases = 0; // variable names are self explaining
    
    bool go = false; // if true; buffer zone starts to take true bases from the beginning
    uintmax_t n_lines_local = 0; // line counter, once reached seg_l, program pause for writing disk
	
	// residue is the leftover sequence due to the break point somewhere in sequence
	// cur_pos is the cursor position that is local to each bufferred sequence fragment, or is its length.
	uintmax_t residue = 0, cur_pos = 0; 
	
	vector<char> kmers;

    auto fh = fstream(out_path, ios::out | ios::binary);

    while (true) {

        const ssize_t bytes_read = read(fileno(stdin), window + residue, step_size);
		
        if (bytes_read == 0) 
        	break;
        	
        if (bytes_read == (ssize_t) -1) {
        	cerr << "unknown fetal error!" << endl;
			exit(EXIT_FAILURE);
		}
		
		cur_pos = 0;
        for (int i = 0;  i < bytes_read + residue;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
            	go = true;

            	++n_lines;
                ++n_lines_local;
                                   	     	        	  
                if (n_lines > max_l)
                	break;
                
            } else if (c == '>'){
            	go = false;
            	++n_reads;
            	
            	if (organized_by_seq) {
            		kmers.push_back('=');
            		kmers.push_back('\n');
            	}
            	
            	if (cur_pos >= k)
            		kmer_search(kmers, window, cur_pos);
				// else kmer can not be found

            	cur_pos = 0;
            } else {
            	if (go) {
            		window[cur_pos++] = c;
            		++n_bases;
            	}
            }
        }
         
        if (n_lines > max_l) {
        	cerr << "reached max line: " << max_l << endl;
        	break;
        }
        
        n_bases = (n_bases-residue);
        
        if (cur_pos < k) {
            residue = cur_pos;
        } else {
        	kmer_search(kmers, window, cur_pos);	
            residue = k - 1;
        }
        
        for (int i = 0;  i < residue;  ++i) {
        	window[i] = window[cur_pos - residue + i];
        }
        
        if (n_lines_local > seg_l) {
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
    vfkmrz_fasta();
    return 0;
}
