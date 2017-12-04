#!/usr/bin/env python

import sys, os, gzip, argparse, ujson, hashlib
from operator import itemgetter
from collections import defaultdict, Counter
from time import time

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Usage: python vfkmrz_wrapper.py -f <fastq or fasta files>  [--out <out>] [--offset <kmer jump size>] [--max-reads <num of reads will be processed>]""")
    parser.add_argument('-1', type=str, dest='fp', required=True,
        help="""FASTA/FASTQ file containing short sequencing reads.
Multiple files can be accepted but should be comma separated.
Otherwise FASTA/FASTQ containing unpaired reads.
Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)""")
    parser.add_argument('--out', type=str, dest='output', default='/dev/stdout',
        help="""Path to output file (/dev/stdout)""")
    parser.add_argument('--max-reads', dest='max_reads', default=float('inf'), type=int,
        help="""Number of reads to use from input file(s) (use all)""")
    parser.add_argument('--k', dest='kmer_size', default=7, type=int,
        help="""Kmer offset when searching reads (2)""")
    parser.add_argument('--offset', dest='kmer_offset', default=1, type=int,
        help="""Kmer offset when searching reads (2)""")
    parser.add_argument('--compile-overwrite', dest='overwrite', action='store_true', help="""Overwrite existing compilation files""")

    args = vars(parser.parse_args())

    return args

def run_command_no_comm(cmd, env=None):
    import subprocess as sp
    p = None
    if env:
        p = sp.call(cmd, shell=True)
    else:
        p = sp.call(cmd, shell=True)
    
    if p != 0:
        err_msg =  "\nError: the following returned non-zero status: '%s':\n" % cmd
        sys.exit(err_msg)

def run_command(cmd, env=None):
	import subprocess as sp
	if env:
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, env=env)
	else:
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = p.communicate()
	if p.returncode != 0:
		err_msg =  "\nError: the following returned non-zero status: '%s':\n" % cmd
		err_msg += "\n%s" % err
		sys.exit(err_msg)
	else:
		return out, err

def main():
    """ Search reads versus kmer database and assign to alleles

    Loops over paired read files.
    For each file pair, loops over paired reads.
    Search each read pair, against kmer database.
    Classify read based on its distribution of hits to alleles.
    Stop when max_reads or all reads have been searched

    Args:
        alleles:        dic; {allele_id:[count_kmers, count_reads, (genome, group, allele)]}
        kmers:            dic; {kmer:allele_id}
        kmer_offset:    int; offset to use when kmerizing reads
        paths_m1:        str; comma-separated list of paths to mate-pair 1
        paths_m2:        str; comma-separated list of paths to mate-pair 2; None if no mate-pair 2
        max_reads:        int; stop when this many reads have been searched; Inf if no max
    """
    sys.stderr.write("\n[kmerize] kmerizing sequencing reads\n")

    args = parse_args()

    kmer_size = args["kmer_size"]
    offset = args["kmer_offset"]
    max_line = args["max_reads"]
    output = args["output"]

    fpaths = input_paths(args["fp"])

    for fpath in fpaths:
        output = do_vfkmrz(fpath, output, kmer_size, offset, max_line)

def input_paths(fpaths):
    """ Formats input paths to sequencing files

    Args:
        paths_m1:    str; comma-separated list of paths to mate-pair 1
        paths_m2:    str; comma-separated list of paths to mate-pair 2; None if no mate-pair 2

    Returns:
        List of tuples, where each tuple indicates one pair of sequencing files.
        If the reads are not paired-end, the second element of the tuple is None.
        Ex: [(/path/to/file1_1.fq, /path/to/file1_2.fq), (/path/to/file2_1.fq, /path/to/file2_2.fq)]
    """
    paths = []
    paths = fpaths.split(',')

    return paths

def do_vfkmrz(dat_path, output, kmer_size, offset, max_line, overwrite=True):
    vfkmrz_paths, vfkmrz_params = vfkmrz_prepare(dat_path, output, kmer_size, offset, max_line)

    if os.path.isfile(vfkmrz_params["dat_path"]):
        compile_vfkmrz_atom(vfkmrz_paths, vfkmrz_params, overwrite)
        output = run_vfkmrz_atom(vfkmrz_paths, vfkmrz_params)

    return output

def compile_vfkmrz_atom(vfkmrz_paths, vfkmrz_params, overwrite=True):
    k, r_len, kpr = vfkmrz_params["kmer_size"], vfkmrz_params["r_len"], vfkmrz_params["n_kmer_per_read"]
    offset, max_line = vfkmrz_params["kmer_offset"], vfkmrz_params["max_line"] 

    sys.stderr.write("[compiling vfkmrz]\n")
    sys.stderr.write("{}parameteters used for compilation\n".format(" "*4))
    sys.stderr.write("{}kmer size: {}\n".format(" "*8, k))
    sys.stderr.write("{}input reads length: {}\n".format(" "*8, r_len))
    sys.stderr.write("{}num of kmer per read: {}\n".format(" "*8, kpr))
    sys.stderr.write("{}kmer jump step size: {}\n".format(" "*8, offset))
    sys.stderr.write("{}max num of lines to process: {}\n".format(" "*8, max_line))

    para_str = "{}-{}-{}-{}-{}-{}".format(vfkmrz_params["dat_path"], k, r_len, kpr, offset, max_line)
    hash_val = hashlib.md5(para_str).hexdigest()

    vfkmrz_paths["vfkmrz_fastq_src"] = vfkmrz_paths["vfkmrz_src"] + "/" + "vfkmrz_fastq.cpp"
    vfkmrz_paths["vfkmrz_fasta_src"] = vfkmrz_paths["vfkmrz_src"] + "/" + "vfkmrz_fasta.cpp"

    vfkmrz_paths["vfkmrz_fastq_tmp"] = "{}/{}_fastq.cpp".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)
    vfkmrz_paths["vfkmrz_fasta_tmp"] = "{}/{}_fasta.cpp".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)

    vfkmrz_paths["vfkmrz_output_tmp"] = "{}/{}.out".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)
    vfkmrz_paths["vfkmrz_errput_tmp"] = "{}/{}.err".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)

    vfkmrz_paths["vfkmrz_fastq_bin"] = "{}/{}_fastq".format(vfkmrz_paths["vfkmrz_bin"], hash_val)
    vfkmrz_paths["vfkmrz_fasta_bin"] = "{}/{}_fasta".format(vfkmrz_paths["vfkmrz_bin"], hash_val)

    if os.path.isfile(vfkmrz_paths["vfkmrz_fastq_bin"]) and os.path.isfile(vfkmrz_paths["vfkmrz_fasta_bin"]) and (overwrite == False):
        sys.stderr.write("compiled binaries for vfkmrz were found, skip compilation\n")
        sys.stderr.write("or use --compile-overwrite to overwrite)\n")
        pass
    else:
        vfkmrz_fastq_lines = []


        with open(vfkmrz_paths["vfkmrz_fastq_src"], "r") as fh:
            for line in fh:
                line = line.rstrip()
                if "constexpr auto k = 31;" == line:
                    print "--------------------replacing kmer size line"
                    vfkmrz_fastq_lines.append("constexpr auto k = {};".format(k))
                elif "constexpr auto r_len = 90;" == line:
                    vfkmrz_fastq_lines.append("constexpr auto r_len = {};".format(r_len))
                elif "constexpr auto offset = 31;" == line:
                    vfkmrz_fastq_lines.append("constexpr auto offset = {};".format(offset))
                elif "constexpr auto max_l = 1000*1000*80;" == line:
                    if max_line == float('inf'):
                        vfkmrz_fasta_lines.append("constexpr auto max_l = numeric_limits<uintmax_t>::max();") 
                    else:
                        vfkmrz_fasta_lines.append("constexpr auto max_l = {};".format(max_line))
                else:
                    vfkmrz_fastq_lines.append(line)

        with open(vfkmrz_paths["vfkmrz_fastq_tmp"], "w") as fh:
            fh.write("\n".join(vfkmrz_fastq_lines))

        vfkmrz_fasta_lines = []
        with open(vfkmrz_paths["vfkmrz_fasta_src"], "r") as fh:
            for line in fh:
                line = line.rstrip()
                if "constexpr auto k = 31;" == line:
                    vfkmrz_fasta_lines.append("constexpr auto k = {};".format(k))
                elif "constexpr auto offset = 31;" == line:
                    vfkmrz_fasta_lines.append("constexpr auto offset = {};".format(offset))
                elif "constexpr auto max_l = 1000*1000*80;" == line:
                    if max_line == float('inf'):
                        vfkmrz_fasta_lines.append("constexpr auto max_l = numeric_limits<uintmax_t>::max();") 
                    else:
                        vfkmrz_fasta_lines.append("constexpr auto max_l = {};".format(max_line))
                else:
                    vfkmrz_fasta_lines.append(line)

        with open(vfkmrz_paths["vfkmrz_fasta_tmp"], "w") as fh:
            fh.write("\n".join(vfkmrz_fasta_lines))

        command = ""
        command += "g++ -O3 --std=c++11 {} -o {} ".format(vfkmrz_paths["vfkmrz_fastq_tmp"], vfkmrz_paths["vfkmrz_fastq_bin"])
        command += " && "
        command += "g++ -O3 --std=c++11 {} -o {} ".format(vfkmrz_paths["vfkmrz_fasta_tmp"], vfkmrz_paths["vfkmrz_fasta_bin"])

        sys.stderr.write("{}compiling vfkmrz\n".format(" "*4))
        sys.stderr.write("{}{}\n".format(" "*8, command))

        environ = os.environ.copy()
        run_command(command, environ)

        sys.stderr.write("{}vfkmrz compilation done!\n".format(" "*4))
        sys.stderr.write("{}vfkmrz binary paths:\n".format(" "*4))
        sys.stderr.write("{}vfkmrz_fastq: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_fastq_bin"]))
        sys.stderr.write("{}vfkmrz_fasta: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_fasta_bin"]))

    return os.path.isfile(vfkmrz_paths["vfkmrz_fastq_bin"]) and os.path.isfile(vfkmrz_paths["vfkmrz_fasta_bin"])

def run_vfkmrz_atom(vfkmrz_paths, vfkmrz_params):
    sys.stderr.write("[executing vfkmrz]\n")

    dat_path = vfkmrz_params["dat_path"]
    dtype = vfkmrz_params["dat_type"]

    sys.stderr.write("{}parameteters used for execution\n".format(" "*4))
    sys.stderr.write("{}input data path: {}\n".format(" "*8, dat_path))
    sys.stderr.write("{}input file type: {}\n".format(" "*8, dtype))

    command = ""
    if dtype == "fasta" or dtype == "fastq":
        command += "cat "
    elif dtype == "fasta.gz" or dtype == "fastq.gz":
        command += "gzip -dc "

    command += dat_path
    command += " | "

    if dtype == "fasta" or dtype == "fasta.gz":
        command += vfkmrz_paths["vfkmrz_fasta_bin"]
    elif dtype == "fastq" or dtype == "fastq.gz":
        command += vfkmrz_paths["vfkmrz_fastq_bin"]

    command += " > "
    command += vfkmrz_paths["vfkmrz_output"]

    sys.stderr.write("{}executing vfkmrz\n".format(" "*4))
    sys.stderr.write("{}{}\n".format(" "*8, command))

    environ = os.environ.copy()
    run_command_no_comm(command, environ)

    sys.stderr.write("{}vfkmrz execution done!\n".format(" "*4))
    sys.stderr.write("{}temporary output file generated.\n".format(" "*4))
    sys.stderr.write("{}output: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_output"]))

    return vfkmrz_paths["vfkmrz_output"]

def vfkmrz_prepare(dat_path, output, kmer_size, offset, max_line):
    vfkmrz_paths = prepare_vfkmrz_paths()
    vfkmrz_paths["vfkmrz_output"] = output

    vfkmrz_params = prepare_vfkmrz_params(dat_path, kmer_size, offset, max_line)

    return vfkmrz_paths, vfkmrz_params

def prepare_vfkmrz_paths():
    parent_dir, src_dir = search_vfkmrz_src()
    assert os.path.isdir(parent_dir) and os.path.isdir(src_dir), "problematic vfkmrz direcory: {} (parent) and {} (src)".format(parent_dir, src_dir)

    sys.stderr.write("[preparing path tree for vfkmrz]\n")

    bin_dir = parent_dir + "/bin"
    tmp_dir = parent_dir + "/tmp"

    try:
        os.makedirs(bin_dir)
    except:
        pass

    try:
        os.makedirs(tmp_dir)
    except: pass

    assert os.path.isdir(bin_dir) and os.path.isdir(tmp_dir), "problematic vfkmrz direcory: {} (bin) and {} (tmp)".format(bin_dir, tmp_dir)

    vfkmrz_paths = {}

    vfkmrz_paths["vfkmrz_src"] = src_dir
    vfkmrz_paths["vfkmrz_bin"] = bin_dir
    vfkmrz_paths["vfkmrz_tmp"] = tmp_dir

    vfkmrz_paths["vfkmrz_fastq_src"] = "{}/vfkmrz_fastq.cpp".format(src_dir)
    vfkmrz_paths["vfkmrz_fasta_src"] = "{}/vfkmrz_fasta.cpp".format(src_dir)

    sys.stderr.write("{}vfkmrz path tree preparation done!\n".format(" "*4))
    sys.stderr.write("{}vfkmrz paths generated.\n".format(" "*4))
    sys.stderr.write("{}src direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_src"]))
    sys.stderr.write("{}bin direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_bin"]))
    sys.stderr.write("{}tmp direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_tmp"]))

    return vfkmrz_paths

def prepare_vfkmrz_params(dat_path, kmer_size, offset, max_line):
    sys.stderr.write("[preparing parameters for vfkmrz]\n")

    vfkmrz_params = {}

    vfkmrz_params["kmer_size"] = kmer_size

    if offset == -1:
        vfkmrz_params["kmer_offset"] = vfkmrz_params["kmer_size"]
    else:
        vfkmrz_params["kmer_offset"] = offset

    vfkmrz_params["dat_path"] = dat_path
    vfkmrz_params["dat_type"] = detect_file_type(dat_path)
    sys.stderr.write("{}vfkmrz input file type detected: {}\n".format(" "*4, vfkmrz_params["dat_type"]))

    if vfkmrz_params["dat_type"] in ["fastq", "fastq.gz"]:
        vfkmrz_params["r_len"] = detect_read_len(dat_path, vfkmrz_params["dat_type"])
        sys.stderr.write("{}vfkmrz input read length detected: {}\n".format(" "*4, vfkmrz_params["r_len"]))

        kmer_pos = range(0, vfkmrz_params["r_len"]-vfkmrz_params["kmer_size"], vfkmrz_params["kmer_offset"])
        vfkmrz_params["n_kmer_per_read"] = len(kmer_pos)
        if kmer_pos[-1] + vfkmrz_params["kmer_offset"] <= vfkmrz_params["r_len"]:
            vfkmrz_params["n_kmer_per_read"] = vfkmrz_params["n_kmer_per_read"] + 1

        sys.stderr.write("{}vfkmrz n_kmer_per_read inferred: {}\n".format(" "*4, vfkmrz_params["n_kmer_per_read"]))

    vfkmrz_params["max_line"] = max_line

    sys.stderr.write("{}vfkmrz parameter preparation done!\n".format(" "*4))

    return vfkmrz_params

def search_vfkmrz_src():
    parent_dir = ""
    src_dir = ""

    targets = ["vfkmrz_fasta.cpp", "vfkmrz_fastq.cpp", "vfkmrz_bunion.cpp", "vfkmrz_match.cpp", "flat_hash_map.hpp"]

    sys.stderr.write("[looking for vfkmrz]\n")

    for _ in os.environ['VFKMRZPATH'].strip(':').split(':'):
        if os.path.isdir(_):
            sys.stderr.write("{}Searching into path: {}\n".format(" "*4, _))

            for __ in os.listdir(_):
                vfkmrz_contents = []

                if os.path.isdir(_ + "/" + __):
                    sys.stderr.write("{}Searching into path: {}\n".format(" "*8, _ + "/" + __))
                    vfkmrz_contents += os.listdir(_ + "/" + __)

                vfkmrz_contents.append(__)

                is_src = True
                for target in targets:
                    if target not in vfkmrz_contents:
                        is_src = False
                        sys.stderr.write("{}{} not found!\n".format(" "*8, target))
                        break
                    else:
                        sys.stderr.write("{}{} found!\n".format(" "*8, target))

                if is_src == True:
                    src_dir = _ + "/" + __
                    parent_dir = _

            if src_dir == "" or parent_dir == "":
                sys.stderr.write("{}Path not eligible: {}\n".format(" "*8, _))
            else:
                sys.stderr.write("{}Eligible vfkmrz path found: {}\n".format(" "*8, _))
                break

    return parent_dir, src_dir


def detect_file_type(fpath):
    """ Determine if sequence file is FASTA or FASTQ format

    Opens input file and checks the first character.
    '>' indicates fasta and '@' indicates fastq.
    Any other character results in sys.exit

    Args:
        fpath:    path to reads file

    Returns:
        fasta, fastq or gz
    """
    splitted = fpath.split(".")

    dtype = ""

    if splitted[-1] in ["fasta", "fastq"]:
        dtype = splitted[-1]
    elif splitted[-1] == "gz":
        if splitted[-2] in ["fasta", "fastq"]:
            dtype = splitted[-2] + "." + splitted[-1]

    assert dtype in ["fasta", "fastq", "fasta.gz", "fastq.gz"]

    return dtype

def detect_read_len(fpath, dtype):
    r_len = -1

    assert dtype in ["fastq", "fastq.gz"]

    command = ""

    if dtype == "fastq":
        command += "cat "
    elif dtype == "fastq.gz":
        command += "gzip -dc "

    command += fpath
    command += " | head"

    environ = os.environ.copy()
    out, err = run_command(command, environ)

    first_read = out.split("\n")[1]
    r_len = len(first_read)

    return r_len


if __name__ == "__main__":
    main()
