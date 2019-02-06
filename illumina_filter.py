#!/usr/bin/env python

"""
This script...
NOTE: the files <barcodes.txt> and <sample_coords.txt> need to be in the working directory, along with the fastq files
and this script.
"""
########################################################################################################################
########################################################################################################################

import argparse
import subprocess
import os
import shutil
import random
import glob


parser = argparse.ArgumentParser(description="NOTE: the files <barcodes.txt> and <sample_coords.txt> need to be in the "
                                             "working directory, along with the fastq files and this script.")
parser.add_argument("read_length", type=int, help="Theoretical maximum read length from sequencing run "
                                                  "(e.g. 150 for 150 bp paired-end reads")
parser.add_argument("min_overlap", type=int, help="minimum overlap value for Pear merge")
parser.add_argument("min_length", type=int, help="minimum length value for Pear merge")
parser.add_argument("max_length", type=int, help="maximum overlap value for Pear merge")
results = parser.parse_args()

########################################################################################################################
########################################################################################################################


# define some functions


def sbatch(job_name, command, nodes=1, ntasks=1, cpus_per_task=2, mem_per_cpu=4, time=59, dep='', part=''):
    """main slurm submission function"""
    if dep != '':
        dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)
    if part != '':
        part = '--partition={}'.format(part)
    log_no = str(random.randint(1, 1000))
    print("\ncpus_per_task: ", cpus_per_task, "\nntasks: ", ntasks, "\nnodes: ", nodes)
    sbatch_command = "sbatch -J {} -o {}.out -e {}.err --mail-user=cjackson1245@gmail.com --mail-type=ALL -t 00:{}:00 \
    --mem-per-cpu={}000 --nodes={} --ntasks={} --cpus-per-task={} --wrap='{}' {} {}".format(
        job_name, job_name + "_" + log_no, job_name + "_" + log_no, time, mem_per_cpu, nodes, ntasks, cpus_per_task,
        command, dep, part)
    sbatch_response = subprocess.getoutput(sbatch_command)
    print(sbatch_response)
    job_id = sbatch_response.split(' ')[-1].strip()
    return job_id


def createfolder(directory):
    """creates a folder if it doesn't already exist"""
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def trimmomatic(dep=""):
    """trims the raw fastq files"""
    command = "module load Trimmomatic/0.38-Java-1.8.0_71 parallel/20170206-GCC-6.2.0; " \
              "srun=\"srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK\"; echo ${{srun}}; " \
              "parallel=\"parallel --delay 1 -j $SLURM_NTASKS --joblog runtask.log --resume\"; echo ${{parallel}}; " \
              "cd ddrad_project/01-trimmed; for i in *; do " \
              "echo ${{srun}} java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 -threads 4 " \
              "../00-raw/${{i}}/${{i}}*R1* ../00-raw/${{i}}/${{i}}*R2* " \
              "./${{i}}/${{i}}_forward_paired.fastq " \
              "./${{i}}/${{i}}_forward_unpaired.fastq " \
              "./${{i}}/${{i}}_reverse_paired.fastq " \
              "./${{i}}/${{i}}_reverse_unpaired.fastq " \
              "ILLUMINACLIP:../../TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 " \
              "SLIDINGWINDOW:4:20 MINLEN:{}; done | ${{parallel}}; cd ../.."\
        .format(trim_min)
    print(command)
    job_id = sbatch("trimmomatic", command, dep=dep, time=59, nodes=1, ntasks=2, cpus_per_task=4)
    return job_id


def fastqc(dep=""):
    """runs fastqc on trimmed files"""
    command = "module load fastqc parallel/20141122-GCC-4.9.2; " \
              "srun=\"srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK\"; " \
              "parallel=\"parallel --delay 1 -j $SLURM_NTASKS --joblog runtask.log --resume\"; " \
              "cd ddrad_project/02-fastqc; for i in *; do " \
              "echo ${srun} fastqc -o ${i} ../01-trimmed/${i}/${i}_forward_paired.fastq " \
              "../01-trimmed/${i}/${i}_reverse_paired.fastq; done | ${parallel}; cd ../.."
    print(command)
    job_id = sbatch("fastqc", command, dep=dep, time=59, ntasks=4, cpus_per_task=1)
    return job_id


def pear(dep=""):
    """merges reads using the program Pear"""
    command = "module load PEAR/0.9.11-GCC-6.2.0 parallel/20170206-GCC-6.2.0; " \
              "srun=\"srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK\" ; " \
              "parallel=\"parallel --delay 1 -j $SLURM_NTASKS --joblog runtask.log --resume\"; " \
              "cd ddrad_project/03-merged_reads; for i in *; do " \
              "echo \"${{srun}} pear -f ../01-trimmed/${{i}}/${{i}}_forward_paired.fastq " \
              "-r ../01-trimmed/${{i}}/${{i}}_reverse_paired.fastq " \
              "-o ${{i}}/${{i}} -v {} -n {} -m {} -j 4 -b 33 -y 2G && cat " \
              "${{i}}/${{i}}.assembled.fastq " \
              "${{i}}/${{i}}.unassembled.forward.fastq > " \
              "${{i}}/${{i}}_merged_unmergedR1.fastq\"; done | ${{parallel}}; cd ../.."\
        .format(results.min_overlap, results.min_length, results.max_length)
    print(command)
    job_id = sbatch('pear', command, dep=dep, time=59, ntasks=2, cpus_per_task=4)
    return job_id


def kraken(dep="", part="", database="standard", standard_db="/data/projects/KRAKEN/standarddb"):
    """removes contaminating reads using the program Kraken"""
    command = "module load Kraken/1.0-intel-2017.u2-Perl-5.24.1 parallel/20170206-GCC-6.2.0; " \
              "srun=\"srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK\" ;" \
              "parallel=\"parallel --delay 1 -j $SLURM_NTASKS --joblog runtask.log --resume\" ;" \
              "cd ddrad_project/04-remove_contaminants; for i in *; do echo \"${{srun}} kraken " \
              "../03-merged_reads/${{i}}/${{i}}_merged_unmergedR1.fastq --unclassified-out ${{i}}/${{i}}_unc.fastq " \
              "--classified-out ${{i}}/${{i}}_classified_{}.fastq --fastq-input --output ${{i}}_output_{} --db {} " \
              "--threads 4\"; done | ${{parallel}}; cd ../..".format(database, database, standard_db)
    print(command)
    job_id = sbatch('kraken', command, dep=dep, ntasks=1, cpus_per_task=4, mem_per_cpu=4)
    return job_id


def demultiplex(dep="", barcode_file="../../barcodes.txt", enzyme="ecoRI"):
    """demultiplexes reads bases on a provided <bacrcodes.txt file> containing barcode sequences"""
    command = "module load Stacks/2.3-spartan_gcc-6.2.0 parallel/20170206-GCC-6.2.0; " \
              "srun=\"srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK\" ; " \
              "parallel=\"parallel --delay 1 -j $SLURM_NTASKS --joblog runtask.log --resume\" ; " \
              "cd ddrad_project/05-demultiplex; for i in *; do echo \"${{srun}} process_radtags " \
              "-f ../04-remove_contaminants/${{i}}/${{i}}_unc.fastq -o ${{i}} -b {} -q -e {} " \
              "-i fastq \"; done | ${{parallel}}; cd ../.. ".format(barcode_file, enzyme)
    print(command)
    job_id = sbatch('demultiplex', command, dep=dep, mem_per_cpu=4, ntasks=2, cpus_per_task=4)
    return job_id


def rename(dep=""):
    """copies and renames fastq files, changing the barcode sequence in the filename to the sample number, according
    to a provided <sample_coords.txt> file"""
    subprocess.run("awk '{ print >> \"\"$1\"_sample_coords.txt\" }' sample_coords.txt", shell=True)
    command = 'module load parallel/20170206-GCC-6.2.0; srun="srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK" ; ' \
              'parallel="parallel --delay 1 -j $SLURM_NTASKS --joblog runtask.log" ; ' \
              'cd ddrad_project/06-rename; ' \
              'for i in *; do ' \
              'echo "while read -r line; do ' \
              r'barcode=\$(echo \$line | cut -d \" \" -f 5) ; ' \
              r'echo \$barcode ;' \
              r'barcode_upper=\${barcode^^} ; ' \
              r'echo \$barcode_upper ;' \
              r'sample=\$(echo \${line} | cut -d \" \" -f 4) ; ' \
              r'echo \$sample; ' \
              r'files_to_move=\$(ls ../05-demultiplex/${i}/*\${barcode_upper}* 2> /dev/null); ' \
              r'echo \$files_to_move; ' \
              r'if [[ \"\${files_to_move}\" != \"\" ]]; then ' \
              r'for name in \${files_to_move}; do ' \
              r'sample_name=\$(echo \$name | sed -r \"s/.*(sample.*)/\\1/\"); ' \
              r'echo \${sample_name}; ' \
              r'cp \"\${name}\" \"${i}/\${sample_name/\${barcode_upper}/\${sample}}\"; done; fi; ' \
              r'done < ../../${i}_sample_coords.txt" ; done | ${parallel} ; cd ../.. '
    print(command)
    job_id = sbatch('rename', command, dep=dep, mem_per_cpu=2, ntasks=1, cpus_per_task=1, time=2)
    return job_id


########################################################################################################################
########################################################################################################################

# do some pre-processing:
file_names = []
for file in glob.glob("*.fastq"):
    ID = file.split("_")[0]
    file_names.append(ID)
unique_IDs = set(file_names)

# create project output folders
folders_to_create = ["ddrad_project/00-raw", "ddrad_project/01-trimmed", "ddrad_project/02-fastqc",
                     "ddrad_project/03-merged_reads", "ddrad_project/04-remove_contaminants",
                     "ddrad_project/05-demultiplex"]
for item in folders_to_create:
    createfolder(item)

for ID in unique_IDs:
    createfolder("ddrad_project/00-raw/{}".format(ID))
    createfolder("ddrad_project/01-trimmed/{}".format(ID))
    createfolder("ddrad_project/02-fastqc/{}".format(ID))
    createfolder("ddrad_project/03-merged_reads/{}".format(ID))
    createfolder("ddrad_project/04-remove_contaminants/{}".format(ID))
    createfolder("ddrad_project/05-demultiplex/{}".format(ID))
    createfolder("ddrad_project/06-rename/{}".format(ID))

for ID in unique_IDs:
    for item in glob.glob("{}*fastq".format(ID)):
        shutil.copy(item, "ddrad_project/00-raw/{}".format(ID))


# calculate minimum read length for Trimmomatic
trim_min = round(results.read_length - (0.05 * results.read_length))
print(trim_min)

# create adapter file for trimmomatic using TruSeq2-PE-2 sequences
with open("TruSeq3-PE-2.fa", "w") as adapter_file:
    adapter_file.write(">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\n"
                       "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PE1_rc\n"
                       "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n>PE2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE2_rc\n"
                       "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC")


########################################################################################################################
########################################################################################################################

# Run the pipeline!

trimmomatic_job_id = trimmomatic()

fastqc_jobid = fastqc()

pear_job_id = pear(dep=trimmomatic_job_id)

kraken_job_id = kraken(dep=pear_job_id, database="minikraken",
                       standard_db="/data/cephfs/punim0569/chris_j/michael/minikraken_20171013_4GB")

demultiplex_job_id = demultiplex(dep=kraken_job_id)

rename_job_id = rename(dep=demultiplex_job_id)

########################################################################################################################
########################################################################################################################
