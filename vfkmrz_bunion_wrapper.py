#!/usr/bin/env python

import sys, os, gzip, argparse, ujson, hashlib
from operator import itemgetter
from collections import defaultdict, Counter
from time import time

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description="""Usage: python vfkmrz_bunion_wrapper.py -1 <m1> -2 <m2> [--out <out>]""")
    parser.add_argument('-k1', type=str, dest='m1', required=True,
        help="""path to the first kmer list""")
    parser.add_argument('-k2', type=str, dest='m2', required=True,
        help="""path to the second kmer list""")
    parser.add_argument('--out', type=str, dest='output', default='/dev/stdout',
        help="""Path to output file (/dev/stdout)""")
    parser.add_argument('--k', dest='kmer_size', default=0, type=int,
        help="""k size for each kmer. All kmers from both input files should have the same k.""")
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

    dat1_path = args["m1"]
    dat2_path = args["m2"]

    output = args["output"]
    
    kmer_size = args["kmer_size"]
    overwrite = args["overwrite"]

    k1 = detect_kmer_size(dat1_path)
    k2 = detect_kmer_size(dat2_path)

    assert k1 == k2                

    if kmer_size > 0:
        assert k1 == kmer_size
    else:
        kmer_size = k1

    do_vfkmrz_bunion(dat1_path, dat2_path, kmer_size, output, overwrite)


def do_vfkmrz_bunion(dat1_path, dat2_path, kmer_size, output, overwrite=True):
    vfkmrz_paths, vfkmrz_params = vfkmrz_prepare(dat1_path, dat2_path, kmer_size, output)

    if os.path.isfile(vfkmrz_params["dat1_path"]) and os.path.isfile(vfkmrz_params["dat2_path"]):
        compile_vfkmrz_atom(vfkmrz_paths, vfkmrz_params, overwrite)
        output = run_vfkmrz_atom(vfkmrz_paths, vfkmrz_params)

    return output

def compile_vfkmrz_atom(vfkmrz_paths, vfkmrz_params, overwrite=True):
    k = vfkmrz_params["kmer_size"]

    sys.stderr.write("[compiling vfkmrz]\n")
    sys.stderr.write("{}parameteters used for compilation\n".format(" "*4))
    sys.stderr.write("{}kmer size: {}\n".format(" "*8, k))

    para_str = "{}".format(k)
    hash_val = hashlib.md5(para_str).hexdigest()

    vfkmrz_paths["vfkmrz_bunion_src"] = vfkmrz_paths["vfkmrz_src"] + "/" + "vfkmrz_bunion.cpp"
    vfkmrz_paths["vfkmrz_bunion_tmp"] = "{}/{}_fastq.cpp".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)
    vfkmrz_paths["vfkmrz_bunion_bin"] = "{}/{}_bunion".format(vfkmrz_paths["vfkmrz_bin"], hash_val)

    vfkmrz_paths["vfkmrz_output_tmp"] = "{}/{}.out".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)
    vfkmrz_paths["vfkmrz_errput_tmp"] = "{}/{}.err".format(vfkmrz_paths["vfkmrz_tmp"], hash_val)


    if os.path.isfile(vfkmrz_paths["vfkmrz_bunion_bin"]) and (overwrite == False):
        sys.stderr.write("compiled binaries for vfkmrz were found, skip compilation\n")
        sys.stderr.write("or use --compile-overwrite to overwrite)\n")
        pass
    else:
        vfkmrz_bunion_lines = []
        with open(vfkmrz_paths["vfkmrz_bunion_src"], "r") as fh:
            for line in fh:
                line = line.rstrip()
                if "constexpr auto k = 31;" == line:
                    vfkmrz_bunion_lines.append("constexpr auto k = {};".format(k))
                else:
                    vfkmrz_bunion_lines.append(line)

        with open(vfkmrz_paths["vfkmrz_bunion_tmp"], "w") as fh:
            fh.write("\n".join(vfkmrz_bunion_lines))

        command = ""
        command += "g++ -O3 --std=c++11 {} -o {} ".format(vfkmrz_paths["vfkmrz_bunion_tmp"], vfkmrz_paths["vfkmrz_bunion_bin"])

        sys.stderr.write("{}compiling vfkmrz\n".format(" "*4))
        sys.stderr.write("{}{}\n".format(" "*8, command))

        environ = os.environ.copy()
        run_command(command, environ)

        sys.stderr.write("{}vfkmrz compilation done!\n".format(" "*4))
        sys.stderr.write("{}vfkmrz binary paths:\n".format(" "*4))
        sys.stderr.write("{}vfkmrz_fastq: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_bunion_bin"]))

    return os.path.isfile(vfkmrz_paths["vfkmrz_bunion_bin"]) 

def run_vfkmrz_atom(vfkmrz_paths, vfkmrz_params):
    sys.stderr.write("[executing vfkmrz]\n")

    dat1_path = vfkmrz_params["dat1_path"]
    dat2_path = vfkmrz_params["dat2_path"]

    sys.stderr.write("{}parameteters used for execution\n".format(" "*4))
    sys.stderr.write("{}input data path (1): {}\n".format(" "*8, dat1_path))
    sys.stderr.write("{}input data path (2): {}\n".format(" "*8, dat2_path))

    command = ""
    command += vfkmrz_paths["vfkmrz_bunion_bin"]
    command += " -k1 "
    command += vfkmrz_params["dat1_path"]
    command += " -k2 "
    command += vfkmrz_params["dat2_path"]


    command += " > "
    command += vfkmrz_params["vfkmrz_output"]

    sys.stderr.write("{}executing vfkmrz\n".format(" "*4))
    sys.stderr.write("{}{}\n".format(" "*8, command))

    environ = os.environ.copy()
    run_command_no_comm(command, environ)

    sys.stderr.write("{}vfkmrz execution done!\n".format(" "*4))
    sys.stderr.write("{}output file generated.\n".format(" "*4))
    sys.stderr.write("{}output: {}\n".format(" "*8, vfkmrz_params["vfkmrz_output"]))

    return vfkmrz_params["vfkmrz_output"]

def vfkmrz_prepare(dat1_path, dat2_path, kmer_size, output):
    vfkmrz_paths = prepare_vfkmrz_paths()

    vfkmrz_params = prepare_vfkmrz_params(dat1_path, dat2_path, kmer_size)
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

    vfkmrz_paths["vfkmrz_bunion_src"] = "{}/vfkmrz_fastq.cpp".format(src_dir)

    sys.stderr.write("{}vfkmrz path tree preparation done!\n".format(" "*4))
    sys.stderr.write("{}vfkmrz paths generated.\n".format(" "*4))
    sys.stderr.write("{}src direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_src"]))
    sys.stderr.write("{}bin direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_bin"]))
    sys.stderr.write("{}tmp direcory: {}\n".format(" "*8, vfkmrz_paths["vfkmrz_tmp"]))

    return vfkmrz_paths

def prepare_vfkmrz_params(dat1_path, dat2_path, kmer_size):
    sys.stderr.write("[preparing parameters for vfkmrz]\n")

    vfkmrz_params = {}

    vfkmrz_params["dat1_path"] = dat1_path
    vfkmrz_params["dat2_path"] = dat2_path
    vfkmrz_params["kmer_size"] = kmer_size

    sys.stderr.write("{}vfkmrz parameter preparation done!\n".format(" "*4))

    return vfkmrz_params

def search_vfkmrz_src():
    parent_dir = ""
    src_dir = ""

    targets = ["vfkmrz_bunion.cpp"]

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
                    if os.path.isdir(_ + "/" + __):
                        src_dir = _ + "/" + __
                    else:
                        src_dir = _
                    parent_dir = _

            if src_dir == "" or parent_dir == "":
                sys.stderr.write("{}Path not eligible: {}\n".format(" "*8, _))
            else:
                sys.stderr.write("{}Eligible vfkmrz path found: {}\n".format(" "*8, _))
                break

    return parent_dir, src_dir

def detect_kmer_size(fpath):
    k = 0

    command = ""
    command += "head "

    command += fpath

    environ = os.environ.copy()
    out, err = run_command(command, environ)

    first_kmer = out.split("\n")[1]
    k = len(first_kmer)

    assert k > 0

    return k

if __name__ == "__main__":
    main()
