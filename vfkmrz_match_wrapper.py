#!/usr/bin/env python

import sys, os, gzip, argparse, ujson, hashlib
from operator import itemgetter
from collections import defaultdict, Counter
from time import time

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description="""Usage: python vfkmrz_match_wrapper.py --db <db> --in <kmer list> --kpr <num of kmer per read> [--max-reads <num of reads>] [--out <out>]""")
    parser.add_argument('--db', type=str, dest='db', required=True,
        help="""directory to the files of kmer databases""")
    parser.add_argument('--in', type=str, dest='input', required=True,
        help="""path the input kmer list file""")
    parser.add_argument('--kpr', type=str, dest='kpr', required=True,
        help="""number of kmer per read; this argument is needed for genotyping, may not be needed for pure kmer counting""")
    parser.add_argument('--out', type=str, dest='output', default='/dev/stdout',
        help="""Path to output file (/dev/stdout)""")
    parser.add_argument('--max-reads', dest='max_reads', default=float('inf'), type=int,
        help="""Number of reads to use from input file(s) (use all)""")
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
    sys.stderr.write("\n[bunion] bit encoding kmer lists and searching the union set\n")
    
    args = parse_args()

    db_path = args["db"]
    dat_path = args["input"]
    kpr = args["kpr"]

    max_line = args["max_reads"]
    output = args["output"]
    
    overwrite = args["overwrite"]

    k1 = detect_kmer_size(db_path, True)
    k2 = detect_kmer_size(dat_path)

    assert k1 == k2

    kmer_size = k1

    do_vfkmrz_match(dat_path, db_path, kmer_size, kpr, output, max_line, overwrite)

def do_vfkmrz_match(dat_path, db_path, kmer_size, kpr, output, max_line, overwrite=True):
    vfkmrz_paths, vfkmrz_params = vfkmrz_prepare(dat_path, db_path, kmer_size, kpr, output, max_line)

    if os.path.isfile(vfkmrz_params["dat_path"]):
        compile_vfkmrz_atom(vfkmrz_paths, vfkmrz_params, overwrite)
        output = run_vfkmrz_atom(vfkmrz_paths, vfkmrz_params)

    return output

def compile_vfkmrz_atom(vfkmrz_paths, vfkmrz_params, overwrite=True):
    k, kpr = vfkmrz_params["kmer_size"], vfkmrz_params["n_kmer_per_read"]
    max_line, db_path = vfkmrz_params["max_line"], vfkmrz_params["vfkmrz_db_list"]

    sys.stderr.write("[compiling vfkmrz]\n")
    sys.stderr.write("{}parameteters used for compilation\n".format(" "*4))
    sys.stderr.write("{}kmer size: {}\n".format(" "*8, k))
    sys.stderr.write("{}num of kmer per read: {}\n".format(" "*8, kpr))
    sys.stderr.write("{}max num of lines to process: {}\n".format(" "*8, max_line))
    sys.stderr.write("{}path of kmer db: {}\n".format(" "*8, db_path))

    para_str = "{}-{}-{}-{}".format(k, kpr, max_line, db_path)
    hash_val = hashlib.md5(para_str).hexdigest()

    vfkmrz_paths["vfkmrz_match_src"] = vfkmrz_paths["vfkmrz_src"] + "/" + "vfkmrz_match.cpp"
    vfkmrz_paths["vfkmrz_match_tmp"] = "{}/{}_match.cpp".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)
    vfkmrz_paths["vfkmrz_match_bin"] = "{}/{}_match".format(vfkmrz_paths["vfkmrz_bin"], hash_val)

    vfkmrz_paths["vfkmrz_output_tmp"] = "{}/{}.out".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)
    vfkmrz_paths["vfkmrz_errput_tmp"] = "{}/{}.err".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)


    if os.path.isfile(vfkmrz_paths["vfkmrz_match_bin"]) and (overwrite == False):
        sys.stderr.write("compiled binaries for vfkmrz were found, skip compilation\n")
        sys.stderr.write("or use --compile-overwrite to overwrite)\n")
        pass
    else:
        vfkmrz_match_lines = []
        with open(vfkmrz_paths["vfkmrz_match_src"], "r") as fh:
            for line in fh:
                if "constexpr auto k = 31;" == line:
                    vfkmrz_match_lines.append("constexpr auto k = {};".format(k))
                elif "constexpr auto offset = 31;" == line:
                    vfkmrz_match_lines.append("constexpr auto offset = {};".format(offset))
                elif "constexpr auto n_kmer_per_read = 3;" == line:
                    vfkmrz_match_lines.append("constexpr auto n_kmer_per_read = {};".format(kpr))
                elif "constexpr auto max_l = 1000*1000*80;" == line:
                    vfkmrz_match_lines.append("constexpr auto max_l = {};".format(max_line))
                elif "constexpr auto db_path = \"\";" == line:
                    vfkmrz_match_lines.append("constexpr auto max_l = \"{}\";".format(db_path))
                else:
                    vfkmrz_match_lines.append(line)

        with open(vfkmrz_paths["vfkmrz_match_tmp"], "w") as fh:
            fh.write("\n".join(vfkmrz_match_lines))

        command = ""
        command += "g++ -O3 --std=c++14 {} -o {} ".format(vfkmrz_paths["vfkmrz_match_src"], vfkmrz_paths["vfkmrz_match_bin"])

        sys.stderr.write("{}compiling vfkmrz\n".format(" "*4))
        sys.stderr.write("{}{}\n".format(" "*8, command))

        environ = os.environ.copy()
        run_command(command, environ)

        sys.stderr.write("{}vfkmrz compilation done!\n".format(" "*4))
        sys.stderr.write("{}vfkmrz binary paths:\n".format(" "*4))
        sys.stderr.write("{}vfkmrz_match: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_match_bin"]))

    return os.path.isfile(vfkmrz_paths["vfkmrz_match_bin"])

def run_vfkmrz_atom(vfkmrz_paths, vfkmrz_params):
    sys.stderr.write("[executing vfkmrz]\n")

    dat_path = vfkmrz_params["dat_path"]

    sys.stderr.write("{}parameteters used for execution\n".format(" "*4))
    sys.stderr.write("{}input data path: {}\n".format(" "*8, dat_path))


    command = ""
    command += "cat "
    command += dat_path
    command += " | "
    command += vfkmrz_paths["vfkmrz_match_bin"]

    command += " > "
    command += vfkmrz_params["vfkmrz_output"]

    sys.stderr.write("{}executing vfkmrz\n".format(" "*4))
    sys.stderr.write("{}{}\n".format(" "*8, command))

    environ = os.environ.copy()
    run_command_no_comm(command, environ)

    sys.stderr.write("{}vfkmrz execution done!\n".format(" "*4))
    sys.stderr.write("{}temporary output file generated.\n".format(" "*4))
    sys.stderr.write("{}tmp output: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_output_tmp"]))

    return vfkmrz_paths["vfkmrz_output"]

def vfkmrz_prepare(dat_path, db_path, kmer_size, kpr, output, max_line):
    vfkmrz_paths = prepare_vfkmrz_paths()

    vfkmrz_params = prepare_vfkmrz_params(dat_path, db_path, kmer_size, kpr, max_line)
    vfkmrz_params["vfkmrz_output"] = output

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

    vfkmrz_paths["vfkmrz_match_src"] = "{}/vfkmrz_match.cpp".format(src_dir)

    sys.stderr.write("{}vfkmrz path tree preparation done!\n".format(" "*4))
    sys.stderr.write("{}vfkmrz paths generated.\n".format(" "*4))
    sys.stderr.write("{}src direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_src"]))
    sys.stderr.write("{}bin direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_bin"]))
    sys.stderr.write("{}tmp direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_tmp"]))

    return vfkmrz_paths

def prepare_vfkmrz_params(dat_path, db_path, kmer_size, kpr, max_line):
    sys.stderr.write("[preparing parameters for vfkmrz match]\n")

    vfkmrz_params = {}

    vfkmrz_params["dat_path"] = dat_path
    vfkmrz_params["vfkmrz_db_list"] = db_path
    vfkmrz_params["kmer_size"] = kmer_size
    vfkmrz_params["n_kmer_per_read"] = kpr
    vfkmrz_params["max_line"] = max_line
    
    sys.stderr.write("{}vfkmrz n_kmer_per_read inferred: {}\n".format(" "*4, vfkmrz_params["n_kmer_per_read"]))

    sys.stderr.write("{}vfkmrz parameter preparation done!\n".format(" "*4))

    return vfkmrz_params

def search_vfkmrz_src():
    parent_dir = ""
    src_dir = ""

    targets = ["vfkmrz_match.cpp", "flat_hash_map.hpp"]

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


def detect_kmer_size(fpath, db=False):
    k = 0

    command = ""
    command += "head "

    command += fpath

    environ = os.environ.copy()
    out, err = run_command(command, environ)

    if db == True:
        first_line = out.split("\n")[0]
        first_kmer = first_line.split(" ")[0]
        k = len(first_kmer)
    else:
        first_kmer = out.split("\n")[0]
        k = len(first_kmer)

    assert k > 0

    return k

if __name__ == "__main__":
    main()
