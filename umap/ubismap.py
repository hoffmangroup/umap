from argparse import ArgumentParser
from datetime import datetime
import os
import pandas as pd
import re
import subprocess
from handle_fasta import FastaHandler


class ArgHandler:
    def __init__(self):
        """Makes an argument class

        Uses ArgumentParser and generates some local arguments as well
        for direct use by other classes and functions.
        """
        parser = ArgumentParser(
            description="This scripts is a wrapper that writes qsub "
            "job submission commands for executing other scripts of the "
            "software in order to identify mappability of a given genome "
            "for a range of various read lengths. This script assumes "
            "that you are using a cluster system that accepts the following "
            "parameters: -N, -e, -o, -t, -tc, -cwd, "
            "-b y, -terse. If these options do not exist in your "
            "cluster, specify -write so you can modify the "
            "qsub job submissions manually.")
        parser.add_argument(
            "fasta_path",
            help="Path to the genome fasta file.")
        parser.add_argument(
            "chrsize_path",
            help="Path to a 2-column file where the first column "
            "is the chromosome name and the second column is its size")
        parser.add_argument(
            "out_dir",
            help="Path to directory to create output files/folders")
        parser.add_argument(
            "queue_name",
            help="Queue name for qsub job submission.")
        parser.add_argument(
            "bowtie_path",
            help="Path to bowtie-build executable")
        parser.add_argument(
            "--kmers",
            default=[24, 36, 50, 100],
            nargs="*",
            help="Kmer length for mappability. e.g. 24 36 50 100")
        parser.add_argument(
            "-GenomeReady",
            action="store_true",
            help="If in the 'out_dir' "
            "there already exists a /chrs and /genome subdirectory "
            "where genome directory has a genome.fasta "
            "with bowtie index suffix as 'BisMap_bowtie.ind' or "
            "'Umap_bowtie.ind' if --Bismap is not specified and "
            "the ./chrs directory has indivudal chromosome "
            "FASTA files, specify this option")
        parser.add_argument(
            "-Bismap",
            action="store_true",
            help="Specify -Bismap if double genome indexing is expected. "
            "This would create a genome that is concatenation of forward "
            "and reverse complement. If -C2T or -G2A is expected, "
            "this must be specified")
        parser.add_argument(
            "-C2T",
            action="store_true",
            help="If --Bismap is provided, specify --C2T or --G2A")
        parser.add_argument(
            "-G2A",
            action="store_true",
            help="If --Bismap is provided, specify --C2T or --G2A")
        parser.add_argument(
            "-ExitAfterIndexing",
            action="store_true",
            help="If you only want the index, specify this option")
        parser.add_argument(
            "-SimultaneousJobs",
            default=120,
            type=int,
            help="Number of jobs to run simultaneously")
        parser.add_argument(
            "-var_id",
            default="SGE_TASK_ID",
            help="Environmental variable for accessing job IDs. "
            "By default is set to SGE_TASK_ID assuming a sungrid "
            "engine environment.")
        parser.add_argument(
            "-write_script",
            default="False",
            help="Specify -write <Path to output job submission file> "
            "if instead of direct execution, "
            "you want to save the job submission file.")
        parser.add_argument(
            "-pipe",
            action="store_true",
            help="If -pipe is specified, the software command "
            "will be piped into the qsub command. For example "
            "instead of: 'qsub -q <queuename> -N <jobname> "
            "python <script>.py <arg1> <arg2>', 'echo python "
            "<script>.py <arg1> <arg2> | qsub -q <queuename> "
            "-N <jobname>' will be written/executed.")
        parser.add_argument(
            "-chunk",
            default=1000000,
            type=int,
            help="Length of chromosomal chunks")
        args = parser.parse_args()
        out_dir = args.out_dir
        chrom_dir, dir_kmers, index_dir = make_dir_structure(out_dir)
        print("{} Jobs will be run simultaneously at each step".format(
                args.SimultaneousJobs))
        # Create paths and key parameters
        chrsize_path = "{}/chrsize.tsv".format(out_dir)
        if args.chrsize_path != chrsize_path:
            subprocess.call(["cp", args.chrsize_path, chrsize_path])
        chrsize_df = pd.read_csv(chrsize_path, sep="\t", names=["Chr", "Size"])
        LenChrs = len(pd.unique(chrsize_df["Chr"]))
        conversion = "None"
        genome_path = "{}/genome.fa".format(index_dir)
        # Make sure write is a new file
        if args.write_script != "False":
            with open(args.write_script, "w") as out_link:
                out_link.write(
                    "# Bismap/Umap job submission commands follow.\n")
        source_dir = os.getcwd()
        if args.C2T:
            if args.G2A:
                raise ValueError("Cannot specify both -C2T and -G2A")
            conversion = "C2T"
        elif args.G2A:
            conversion = "G2A"
        if args.Bismap:
            # It is possible that chrx_RC information does not exist
            # in the chrsize.tsv file. The function below would check
            # and add this information if necessary and would update
            # chromosome length
            len_chrs_after_complement = complement_chrsize_file(chrsize_path)
            if len_chrs_after_complement != LenChrs:
                LenChrs = len_chrs_after_complement
                print "Length of chromosomes changed to %d to account for\
                reverse complements" % LenChrs
        idx_path = index_unique_kmer_jobids(chrsize_path, args.chunk)
        self.genome_path = genome_path
        self.chrom_dir = chrom_dir
        self.dir_kmers = dir_kmers
        self.chrsize_path = chrsize_path
        self.out_dir = out_dir
        self.idx_path = idx_path
        self.conversion = conversion
        self.Bismap = args.Bismap
        self.queue_name = args.queue_name
        self.kmers = args.kmers
        self.SimultaneousJobs = args.SimultaneousJobs
        self.ExitAfterIndexing = args.ExitAfterIndexing
        self.var_id = args.var_id
        self.write_script = args.write_script
        self.source_dir = source_dir
        self.GenomeReady = args.GenomeReady
        self.fasta_path = args.fasta_path
        self.bowtie_path = args.bowtie_path
        self.LenChrs = LenChrs
        self.pipe = args.pipe


def pipe_job(job_list):
    '''Assumes the qsub command has the
    qsub ... -cwd -b y <executing command>.
    It will change it to: echo <executing command>
    | qsub ...'''
    for ind_job in range(len(job_list)):
        if job_list[ind_job] == "-cwd":
            if job_list[ind_job + 1] == "-b":
                if job_list[ind_job + 2] == "y":
                    ind_cwd = ind_job
    try:
        qsub_list = job_list[:ind_cwd]
        exec_list = job_list[(ind_cwd + 3):]
        piped_list = [
            ['echo "'] + exec_list + ['" | '] + qsub_list]
        return piped_list
    except:
        raise ValueError(
            "The submitted command to pipe_job must "
            "include -cwd -b y.")


def index_genome(queue_name, converted_path,
                 bowtie_build_path, Bismap, write_path, pipe):
    path_genome = "{}/genome.fa".format(converted_path)
    if Bismap:
        index_suffix = "BisMap_bowtie.ind"
    else:
        index_suffix = "Umap_bowtie.ind"
    index_job_line = ["qsub", "-q", queue_name, "-terse",
                      "-N", "Index-Bowtie", "-o",
                      "{}/index_genome.LOG".format(converted_path),
                      "-e", "{}/index_genome.ERR".format(converted_path),
                      "-cwd", "-b", "y",
                      bowtie_build_path, path_genome,
                      "{}/{}".format(converted_path, index_suffix)]
    if pipe:
        index_job_line = pipe_job(index_job_line)
    if write_path != "False":
        job_num = write_job(index_job_line, write_path)
    else:
        try:
            job_line = subprocess.Popen(index_job_line, stdout=subprocess.PIPE)
            job_num = job_line.stdout.read().rstrip()
        except:
            raise ValueError(
                "Your system is incompatible with direct job "
                "submission. Please specify -write and manually "
                "modify the job submission options.")
    print "job id %s from indexing genome" % job_num
    return index_suffix, job_num


def get_unique_kmer_jobnums(idx_path):
    ind_df = pd.read_csv(idx_path, sep="\t")
    last_ind = ind_df.shape[0] + 1
    return last_ind


def write_job(list_out, write_path):
    write_link = open(write_path, "a")
    out_str = "\nJOBID=$({})\n".format(
        " ".join(list_out))
    write_link.write(out_str)
    write_link.write(
        'IDPARTS=($(echo $JOBID | tr "." "\\n"))\n')
    write_link.write(
        'WAITID=${IDPARTS[0]}\n')
    write_link.write(
        'echo "Submitted JOB ID $WAITID"\n')
    write_link.close()
    return "$WAITID"


def make_unique_kmers(chrom_dir, dir_kmers,
                      kmer_size, queue_name,
                      source_dir, conversion_job_id,
                      idx_path, chrsize_path,
                      write_path, var_id,
                      job_lim, pipe):
    """Generates kmer.gz files of the genome
    """
    # chr_files = os.listdir(chrom_dir)
    # chr_fas = sublist_by_reg(".fa$", chr_files)
    out_dir = "{}/k{}".format(dir_kmers, kmer_size)
    ind_jobs = get_unique_kmer_jobnums(idx_path)
    array_param = "1-{}".format(ind_jobs)
    get_kmers = ["qsub", "-q",
                 queue_name, "-t",
                 array_param, "-N",
                 "Bismap.UniqueKmers",
                 "-terse",
                 "-tc", str(job_lim),
                 "-hold_jid", conversion_job_id,
                 "-cwd", "-b", "y",
                 "-o", "{}/Bismap.UniqueKmers.LOG".format(dir_kmers),
                 "-e", "{}/Bismap.UniqueKmers.ERR".format(dir_kmers),
                 "python", "get_kmers.py",
                 chrsize_path, out_dir,
                 chrom_dir, idx_path,
                 "--var_id", var_id, "--kmer", "k{}".format(kmer_size)]
    if pipe:
        get_kmers = pipe_job(get_kmers)
    if write_path != "False":
        job_num = write_job(get_kmers, write_path)
    else:
        try:
            job_line = subprocess.Popen(get_kmers, stdout=subprocess.PIPE)
        except:
            raise ValueError("Incompatible system? Try -write.")
        job_num = job_line.stdout.read().rstrip()
        print "job id %s from creating unique kmers" % job_num
        job_num = job_num.split(".")[0]
    return job_num, ind_jobs


def run_bowtie(queue_name, dir_kmers, bowtie_path, index_dir,
               index_suffix, kmer_job_id,
               index_job_id, kmer, source_dir,
               Bismap, num_jobs, write_path, var_id,
               job_lim, pipe):
    """Calls bowtie on kmer.gz files
    """
    bowtie_dir = "/".join(bowtie_path.split("/")[:-1])
    kmer_folder = "{}/k{}".format(dir_kmers, kmer)
    array_param = "1-{}".format(num_jobs)
    wait_param = "%s,%s" % (index_job_id, kmer_job_id)
    run_bowtie = ["qsub", "-q",
                  queue_name, "-t",
                  array_param, "-N",
                  "Bismap.RunBowtie",
                  "-terse",
                  "-tc", str(job_lim),
                  "-hold_jid", wait_param,
                  "-cwd", "-b", "y",
                  # "-l", "\"mem_requested={%s}G\"" % (mem_GB),
                  "-o", "{}/Bismap.RunBowtie.LOG".format(dir_kmers),
                  "-e", "{}/Bismap.RunBowtie.ERR".format(dir_kmers),
                  "python", "run_bowtie.py", kmer_folder,
                  bowtie_dir, index_dir, index_suffix,
                  "-var_id", var_id]
    if Bismap:
        run_bowtie.append("-Bismap")
    if pipe:
        run_bowtie = pipe_job(run_bowtie)
    if write_path != "False":
        job_num = write_job(run_bowtie, write_path)
    else:
        try:
            job_line = subprocess.Popen(run_bowtie, stdout=subprocess.PIPE)
        except:
            raise ValueError("Incompatible system? Try -write")
        job_num = job_line.stdout.read()
        print "job id %s from mapping with bowtie" % job_num
        job_num = job_num.split(".")[0]
    return job_num


def unify_bowtie(queue_name, bowtie_job_id, dir_kmers,
                 kmer, chrom_dir, source_dir,
                 LenChrs, write_path, var_id,
                 chrsize_path, job_lim, pipe):
    """Generates qsub command for interpreting bowtie outputs

    Saves collectinve bowtie.gz files of each chromosome into
    a binary unsigned 8-bit integer file for each chromosome
    which is an array with the same length as the chromosome
    """
    kmer_folder = "{}/k{}".format(dir_kmers, kmer)
    len_chrs = LenChrs
    array_param = "1-%d:1" % (len_chrs)
    chr_order_path = chrsize_path
    wait_param = bowtie_job_id
    unify_bw_lst = ["qsub", "-q",
                    queue_name,
                    "-t", array_param,
                    "-N",
                    "Bismap.UnifyBowtie",
                    "-terse",
                    "-tc", str(job_lim),
                    "-hold_jid", wait_param,
                    "-cwd", "-b", "y",
                    "-o", "{}/Bismap.UnifyBowtie.LOG".format(dir_kmers),
                    "-e", "{}/Bismap.UnifyBowtie.ERR".format(dir_kmers),
                    "python", "unify_bowtie.py",
                    kmer_folder, chr_order_path, "-var_id", var_id]
    if pipe:
        unify_bw_lst = pipe_job(unify_bw_lst)
    if write_path != "False":
        job_num = write_job(unify_bw_lst, write_path)
    else:
        try:
            job_line = subprocess.Popen(unify_bw_lst, stdout=subprocess.PIPE)
        except:
            raise ValueError("Incompatible system? try -write")
        job_num = job_line.stdout.read()
        job_num = job_num.split(".")[0]
    TempDir = "%s/TEMPs" % kmer_folder
    make_dir(TempDir)
    move_intermediary_files = ["qsub", "-q",
                               queue_name, "-N",
                               "Moving.Intermediary.Files",
                               "-hold_jid", job_num,
                               "-terse",
                               "-cwd", "-b", "y",
                               "-o", "Bismap.FileMov.LOG",
                               "-e", "Bismap.FileMov.ERR",
                               "mv",
                               "%s/*kmer*" % (kmer_folder),
                               TempDir]
    if write_path != "False":
        job_num = write_job(move_intermediary_files, write_path)
    else:
        job_line = subprocess.Popen(
            move_intermediary_files, stdout=subprocess.PIPE)
    move_intermediary_files_2 = move_intermediary_files[:-2] +\
        ["%s/*bowtie*" % (kmer_folder), TempDir]
    if write_path != "False":
        job_num = write_job(move_intermediary_files_2, write_path)
    else:
        job_line = subprocess.Popen(
            move_intermediary_files_2,
            stdout=subprocess.PIPE)
        job_num = job_line.stdout.read()
    return job_num


def sublist_by_reg(regex, my_list):
    subsetted_list = []
    for item in my_list:
        RE = re.search(regex, item)
        if RE:
            subsetted_list.append(item)
    return subsetted_list


def make_dir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def make_dir_structure(out_dir):
    """Generates subfolders needed for Umap in out_dir
    """
    chrom_dir = "{}/chrs".format(out_dir)
    path_genome = "{}/genome".format(out_dir)
    dir_kmers = "{}/kmers".format(out_dir)
    make_dir(chrom_dir)
    make_dir(path_genome)
    make_dir(dir_kmers)
    return chrom_dir, dir_kmers, path_genome


def combine_files(bowtie_unify_id, source_dir, dir_kmers, queue_name,
                  chrsize_path, write_path, var_id, LenChrs, pipe):
    """Creates the qsub command for merging data of several k-mers
    """
    combine_job_line = ["qsub", "-q",
                        queue_name, "-N",
                        "CombineUmappedFiles",
                        "-terse",
                        "-t", "1-{}".format(LenChrs),
                        "-hold_jid", bowtie_unify_id,
                        "-cwd", "-b", "y",
                        "-o", "{}/Bismap.combine.LOG".format(dir_kmers),
                        "-e", "{}/Bismap.combine.ERR".format(dir_kmers),
                        "python", "combine_umaps.py", dir_kmers,
                        chrsize_path, "-var_id", var_id]
    if pipe:
        combine_job_line = pipe_job(combine_job_line)
    if write_path != "False":
        write_job(combine_job_line, write_path)
    else:
        subprocess.call(combine_job_line)


def complement_chrsize_file(chrsize_path):
    """Complement chromosome size file

    If using Bismap, an extra set of chromosomes
    will be generated with the _RC suffix.
    These need to be added to the chromosome size
    file as well.

    :param chrsize_path: Path to 2-column file of chrom\tsize\n

    :returns: The actual number of chromosomes
    """
    if not os.path.exists(chrsize_path):
        print "%s didn't exist!" % chrsize_path
    ad_lines = []
    len_chrs = 0
    NO_RC = True
    with open(chrsize_path, "r") as chrsize_link:
        for each_line in chrsize_link:
            len_chrs = len_chrs + 1
            if "_RC" in each_line:
                NO_RC = False
            chrsize_list = each_line.rstrip().split("\t")
            chrsize_list[0] = chrsize_list[0] + "_RC"
            ad_lines.append("\t".join(chrsize_list) + "\n")
    if len(ad_lines) > 0 and NO_RC:
        with open(chrsize_path, "a") as chrsize_link:
            for each_line in ad_lines:
                len_chrs = len_chrs + 1
                print "adding %s" % each_line
                chrsize_link.write(each_line)
        print "Added reverse complemented chromosomes"
    else:
        "Reverse complemented chromosomes existed"
    return len_chrs


def process_genome(GenomeReady, genome_path,
                   Bismap, fasta_path,
                   conversion, out_dir, chr_dir,
                   queue_name, bowtie_path,
                   chrsize_path, write_script, pipe):
    """Generates necessary fasta files for Umap and Bismap

    Umap and Bismap require fasta files for each chromosome.
    Bismap requires the fasta files to have a C>T or G>A
    conversion as well. Additionally, Bismap requires a
    reverse complemented set of chromosome where the
    C>T or G>A conversion occurs after reverse complementation.

    :param GenomeRead: Boolean indicating if index of genome exists
    :param genome_path: Path to the output fasta file of the genome
    :param Bismap: Boolean indicating if Bismap option is selected
    :param fasta_path: Path to the already existing input fasta file
    :param conversion: One of None, C2T or G2A
    :param out_dir: Directory for Umap output
    :param chr_dir: Directory for making individual chromosome fastas
    :param bowtie_path: Path to executable bowtie-build binary
    :param chrsize_path: Path to 2-column file of chromosome\tsize\n
    :param write_script: File path for writing the qsub command
        or the string 'False'
    :param pipe: Boolean indicating if the qsub command should be piped

    :returns: index_suffix, index_job_id (qsub wait IDs)
    """
    if GenomeReady:
        print("Assuming the genome and index files exist")
        if Bismap:
            index_suffix = "BisMap_bowtie.ind"
        else:
            index_suffix = "Umap_bowtie.ind"
        index_job_id = "1"
    else:
        print "Started copying/reverse complementing/converting"
        if not os.path.exists(genome_path):
            if Bismap:
                print "Bismap reverse complementation started\
                at %s" % str(datetime.now())
                FastaObj = FastaHandler(fasta_path, genome_path, chrsize_path,
                                        chr_dir, True, conversion)
            else:
                print "Umap genome is being processed"
                FastaObj = FastaHandler(fasta_path, genome_path, chrsize_path,
                                        chr_dir, False, "None")
            FastaObj.handle_fasta()
        elif Bismap:
            print "Assuming that %s/genome/genome.fa includes reverse\
            complemented chromosomes." % out_dir
        print "Indexing the genome started at %s" % str(datetime.now())
        index_suffix, index_job_id = index_genome(
            queue_name,
            "{}/genome".format(out_dir),
            bowtie_path, Bismap, write_script, pipe)
        print "Done with indexing at %s" % str(datetime.now())
    return index_suffix, index_job_id


def index_unique_kmer_jobids(chrsize_path, CHUNK_SIZE=1e6):
    """Generate index file mapping job_ids to chromosomal chunks

    Umap (or Bismap) requires around 200 hours of computation.
    All of the scripts accept -job_id arguments to perform the
    embarrassingly parallel process on different chromosomal
    chunks using a 0-based index that refers to part of the
    genome. This script creates a 4-column index file with
    the following columns:
    Index\tChromosome\tStart\tend\n

    :param chrsize_path: Path to 2-column file of chrom\tsize\n
    :param CHUNK_SIZE: The maximum length of each chromosomal chunk
    :type CHUNK_SIZE: Integer

    :returns: Path to index file
    """
    CHUNK_SIZE = int(CHUNK_SIZE)
    ind_path = "{}/chrsize_index.tsv".format(
        "/".join(chrsize_path.split("/")[:-1]))
    if not os.path.exists(ind_path):
        ind_link = open(ind_path, "w")
        ind_link.write("Index\tChromosome\tStart\tEnd\n")
        start = 1
        ind = 0
        # Identifying chromosome, start, end for job_id
        with open(chrsize_path, "r") as chrsize_link:
            for chrsize_line in chrsize_link:
                chr, len_chr = chrsize_line.rstrip().split("\t")
                end = int(len_chr) + start
                for pos in range(start, end, CHUNK_SIZE):
                    if pos < end:
                        if pos + CHUNK_SIZE - 1 > end:
                            pos_end = end
                        else:
                            pos_end = pos + CHUNK_SIZE - 1
                        ind_link.write(
                            "\t".join(
                                [str(ind), chr, str(pos), str(int(pos_end))]) +
                            "\n")
                        ind = ind + 1
        ind_link.close()
    return ind_path


if __name__ == "__main__":
    args = ArgHandler()
    # genome_path, chrom_dir, dir_kmers, chrsize_path, out_dir = args_list[:5]
    # idx_path, conversion, Bismap, queue_name = args_list[5:9]
    # kmers, SimultaneousJobs, ExitAfterIndexing = args_list[9:12]
    # var_id, write_script, source_dir, GenomeReady = args_list[12:16]
    # fasta_path, bowtie_path, LenChrs, pipe = args_list[16:]
    index_dir = os.path.dirname(args.genome_path)
    conversion_job_id = "1"  # initiate job dependency with an invalid ID
    index_suffix, index_job_id = process_genome(
        args.GenomeReady, args.genome_path,
        args.Bismap, args.fasta_path,
        args.conversion, args.out_dir,
        args.chrom_dir, args.queue_name,
        args.bowtie_path, args.chrsize_path,
        args.write_script, args.pipe)
    if args.ExitAfterIndexing:
        pass
        print "Exiting because -ExitAfterIndexing is specified"
    else:
        for kmer in args.kmers:
            kmer = int(kmer)
            make_dir("{}/k{}".format(args.dir_kmers, kmer))
            kmer_job_id, num_jobs = make_unique_kmers(
                args.chrom_dir, args.dir_kmers, kmer,
                args.queue_name, args.source_dir,
                conversion_job_id,
                args.idx_path, args.chrsize_path,
                args.write_script, args.var_id,
                args.SimultaneousJobs, args.pipe)
            bowtie_job_id = run_bowtie(
                args.queue_name, args.dir_kmers,
                args.bowtie_path,
                index_dir, index_suffix,
                kmer_job_id, index_job_id, kmer,
                args.source_dir, args.Bismap,
                num_jobs, args.write_script, args.var_id,
                args.SimultaneousJobs, args.pipe)
            bowtie_unify_id = unify_bowtie(
                args.queue_name, bowtie_job_id,
                args.dir_kmers,
                kmer, args.chrom_dir,
                args.source_dir, args.LenChrs,
                args.write_script,
                args.var_id, args.chrsize_path,
                args.SimultaneousJobs, args.pipe)
            conversion_job_id = kmer_job_id
            print "Jobs submitted for k%d" % kmer
        combine_files(bowtie_unify_id, args.source_dir,
                      args.dir_kmers,
                      args.queue_name,
                      args.chrsize_path, args.write_script,
                      args.var_id,
                      args.LenChrs, args.pipe)
        conversion_job_id = bowtie_unify_id
        print("Successfully done with creating all jobs")
