#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>

#include <unistd.h>
#include <assert.h>

using namespace std;

constexpr auto r_len = 90;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 128 * 1024 * 1024;
constexpr auto buffer_size = 128 * 1024 * 1024;

// maximum lines when reached; also the max memory controller
constexpr auto seg_l = 1000*1000*10;

// output file path
constexpr auto out_path = "/dev/stdout";


void fastq2fasta() {	
    vector<char> buffer(buffer_size);
    char* window = buffer.data();

    int l_label = 2;
    uintmax_t n_lines = 0, n_lines_local = 0;
	
	vector<char> padded;

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
            char c = toupper(window[i]);
			
            if (c == '\n') {
            	++n_lines;
                ++n_lines_local;
                
                if (l_label == 3) {
					l_label = -1;
                    padded.push_back('\n');
				} else if (l_label == 2) {
                    padded.push_back('\n');
                }

                ++l_label;
            } else {
            	if (l_label == 3) {
                    padded.push_back(c);
            	} else if (l_label == 2) {
					if (c == '@') {
                        padded.push_back('>');
					} else {
                        padded.push_back(c); 
					}
				}
            }
        }
        
        if (n_lines_local > seg_l){
    		fh.write(&padded[0], padded.size());
    		
        	padded.clear();
    		n_lines_local = 0;

			cout << n_lines << " lines were scanned.\n";
    	}
    }
    
    if (padded.size() != 0) {
    	fh.write(&padded[0], padded.size());
    		
    	padded.clear();
    	n_lines_local = 0;

		cout << n_lines << " lines were scanned.\n";
    }
    
    fh.close();
    
    cerr << "done!" << endl;
}

int main(int argc, char** argv){		
    fastq2fasta();
    return 0;
}

