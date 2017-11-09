#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstring>

#include <unistd.h>
#include <assert.h>

using namespace std;

constexpr auto step_size = 128 * 1024 * 1024;
constexpr auto buffer_size = 128 * 1024 * 1024;

void fasta_stats() {	
    cerr << "read and count in fasta" << endl;

    vector<char> buffer(buffer_size); // buffer zone for both file reading and sequence windowing 
    char* window = buffer.data();
	
    uintmax_t n_lines = 0, n_reads = 0, n_bases = 0; // variable names are self explaining
    
    bool go = false; // if true; buffer zone starts to take true bases from the beginning
	
    while (true) {

        const ssize_t bytes_read = read(fileno(stdin), window, step_size);
		
        if (bytes_read == 0) {
            break;
        }

        if (bytes_read == (ssize_t) -1) {
            throw runtime_error("unknown fetal error!");
        }

        for (int i = 0;  i < bytes_read;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
            	go = true;
            	++n_lines;
            } else if (c == '>'){
            	go = false;
            	++n_reads;
            } else {
            	if (go) {
            		++n_bases;
            	}
            }
        }
        
    }
    
    cerr << "total lines: " << n_lines << endl;
    cerr << "total reads: " << n_reads << endl;
	cerr << "total bases: " << n_bases << endl;
}

void fastq_stats() {
    cerr << "read and count in fastq" << endl;

    vector<char> buffer(buffer_size); // buffer zone for both file reading and sequence windowing 
    char* window = buffer.data();

    uintmax_t n_lines = 0, n_reads = 0, n_bases = 0;

    int l_label = 2;
    uintmax_t n_lines_local = 0;
	
	int cur_pos = 0;

    while (true) {

        const ssize_t bytes_read = read(fileno(stdin), window, step_size);
		
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
                ++n_lines_local;
                
                if (l_label == 3) {
                	++n_reads; 
                	l_label = cur_pos = 0;
                } else
                	++l_label;
            } else {
            	if (l_label == 3) {
            		++n_bases;
            	}
            }
        }
        
    }
    
    
    cerr << "done!" << endl;
    cerr << "total lines: " << n_lines << endl;
    cerr << "total reads: " << n_reads << endl;
	cerr << "total bases: " << n_bases << endl;     
}

int main(int argc, char** argv){		
	cout << argc << ": " << argv[0] << endl;
	cout << argc << ": " << argv[1] << endl;

    if (argc == 1) {
	    fastq_stats();
    } else if (argc == 2) {
        if (!strcmp(argv[1], "-fa")){
            fasta_stats();
        } else if (!strcmp(argv[1], "-fq")){
            fastq_stats();
        } else {
            cerr << "unknown fetal error!" << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        cerr << "unknown fetal error!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

