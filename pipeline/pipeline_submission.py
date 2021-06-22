#!/usr/bin/env python
"""
Submission script for metatranscriptomics pipelie
"""
import os
import inspect
from pathlib import Path
from subprocess import run, PIPE
from typing import Tuple
import argparse

dependency = ' --dependency=afterok:%d'
submit_cmd = "sbatch -A def-dsteinke -J %s -c 32 --mem=0 -t %d-%d:%d:00 -o %s"
prologue = f'--wrap=cp -r %s/%s ${{SLURM_TMPDIR}} && cd ${{SLURM_TMPDIR}} && '
epilogue = 'cp -r {} {}'

def process_run(instance, **kwargs):
    if instance is not None:
        return int(instance.stdout.strip().split()[-1]), instance.returncode == 0
    else:
        phiter = kwargs['y']
        return next(phiter), True


def dprint(*args, **kwargs):
    print(args, '\n')


def init() -> Tuple[Path, Path]:
    """
    Setup some initial paths
    :return: Path to pipeline and to trimming
    """
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    pipename = 'METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE'
    cwd = Path().resolve()
    pipepath = cwd.joinpath(pipename)
    pipepath.mkdir(parents=True, exist_ok=True)
    step1_dir = pipepath.joinpath('Step1_trimming').resolve()
    step1_dir.mkdir(exist_ok=True, parents=True)
    return pipepath, step1_dir


def format_step1(step1_direct: Path, score: int, min_length:int = 25) -> str:
    """
    Format individual calls to be passed as wrap
    :param pref: input of reads
    :param score: Phred score to trim by
    :param outdir: Path to the direct the outputs
    :param min_length: lenght threshold
    :return: String with the command
    """
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    epi = epilogue.format('fastqc_reports Trimming', step1_direct)
    exe_dir = os.environ['EBROOTTRIMMOMATIC']
    trimm_out = '${SLURM_TMPDIR}/Trimming'
    pref = '{= s/_R1.fastq// =}'
    baseout = f'{trimm_out}/{pref}_phred{score}'
    illumina_opt = f'{exe_dir}/adapters/TruSeq3-PE.fa:2:30:10'
    fastqc_output = '${SLURM_TMPDIR}/fastqc_reports'
    line = f'java -jar {exe_dir}/trimmomatic-0.39.jar PE {pref}_R1.fastq ' \
           f'{pref}_R2.fastq ILLUMINACLIP:{illumina_opt} LEADING:{score} ' \
           f'TRAILING:{score} SLIDINGWINDOW:4:{score} MINLEN:{min_length} ' \
           f'-threads 6 -baseout {baseout}.fastq && spades.py -1 ' \
           f'{baseout}_1P.fastq -2 {baseout}_2P.fastq --only-error-correction' \
           f' --disable-gzip-output -o {trimm_out}/error_correction -t 6 && ' \
           f"sed -r -i 's/ BH:.{{2,6}}//g' {trimm_out}/error_correction/corrected/" \
           f'{pref}*.fastq && mv {trimm_out}/error_correction/corrected/*1P.00' \
           f'.0_0.cor.fastq {trimm_out}/{pref}_R1_error_corrected.fastq && ' \
           f'mv {trimm_out}/error_correction/corrected/*2P.00.0_0.cor.fastq' \
           f' {trimm_out}/{pref}_R2_error_corrected.fastq && fastqc ' \
           f'{trimm_out}/{pref}.fastq -o {fastqc_output} -t 6 && {epi}'
    return line


# STEP 1 on all inputs for all Phreds
def step1_trimming(input_path: Path, step1_dir: Path):
    """
    Submit step one by phred, using parallel to run all inputs in batches of 5
    :param inputs: List with path to inputs
    :param outdir: Directory where to place the outputs
    :return:
    """
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    wrap_cmd = (prologue % (input_path, '*')) + 'parallel -j5 %s ::: ${SLURM_TMPDIR}/*.fastq'
    phreds = [5, 10, 15, 20]
    phiter = iter(phreds)
    step1_insts = [run(
        (submit_cmd % (f'Step1_phred{score}', 0, 3, 0, f'step1_{score}_%A.out')
         ).split() + [wrap_cmd % (format_step1(step1_dir, score))], stdout=PIPE
    ) for score in phreds]
    step1_jobIDs, step1_returncodes = zip(*[process_run(x, y=phiter)
                                            for x in step1_insts])
    step1_ok = all(step1_returncodes)
    assert step1_ok
    step1 = zip(phreds, step1_jobIDs)
    return step1, step1_dir


def sortmerna(wrap_line: str, refpath: Path, outdir: Path, step1_res: zip,
              trim_path: Path):
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    phred, jid = step1_res
    subline = submit_cmd  % ('Sortmerna', 0, 2, 30, 'step2_sortmerna_%%A.out')
    subline += dependency % jid
    wrap_line = wrap_line + f'cp -R ${{SLURM_TMPDIR}}/SORTMERNA {outdir}'
    perl = '{= s/1P_error_corrected.fastq// =}'
    line = f'mkdir -p ${{SLURM_TMPDIR}}/SORTMERNA && sortmerna --ref ' \
           f'{refpath}/silva-bac-16s-id90.fasta --ref ' \
           f'{refpath}/silva-arc-16s-id95.fasta --ref {refpath}/silva-' \
           f'euk-18s-id95.fasta --ref {refpath}/silva-euk-28s-id98.fasta' \
           f' --ref {refpath}/silva-arc-23s-id98.fasta --ref {refpath}/' \
           f'silva-bac-23s-id98.fasta --ref {refpath}/rfam-5.8s-database' \
           f'-id98.fasta --ref {refpath}/rfam-5s-database-id98.fasta ' \
           f'--reads {perl}1P_error_corrected.fastq --reads ' \
           f'{perl}2P_error_corrected.fastq --paired_in --out2 -other -fastx' \
           f' 1 -num_alignments 1 -v -workdir ${{SLURM_TMPDIR}}/SORTMERNA ' \
           f'--threads 1:1:3'
    args = subline.split() + [wrap_line % (trim_path, f'*phred{phred}*', 11, line)]
    return args

def barnap(wrap_line: str, refpath: Path, outdir: Path, step1_res: zip,
              trim_path: Path):
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    pass

def rnafilter(wrap_line: str, refpath: Path, outdir: Path, step1_res: zip,
              trim_path: Path):
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    pass

def keiner(wrap_line: str, refpath: Path, outdir: Path, step1_res: zip,
           trim_path: Path):
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    pass


# Step2
def step2_sorting(step1: zip, step1_dir: Path, refpath: Path, outdir: Path):
    """
    Submit step2 by tool. Multiprocess samples within
    :param step1_jobIDs: List with jobids
    :param step1_dir: Path of step1 results
    :return:
    """
    frame = inspect.currentframe()
    print(f'Executing {frame.f_code.co_name} with {frame.f_locals} arguments\n')
    trim_dir = step1_dir.joinpath('Trimming')
    step2_dir = outdir.joinpath('Step2_sorting').resolve()
    step2_dir.mkdir(exist_ok=True, parents=True)
    wrap_cmd = prologue + ' parallel -j%d %s ::: ${SLURM_TMPDIR}/*1P_error_corrected.fastq && '
    jobs = [run(x(wrap_cmd, refpath, step2_dir, step1_res, trim_dir))
            for step1_res in step1 for x in (sortmerna, barnap, rnafilter, keiner)]
    return jobs

if __name__ == "__main__":
    opts = argparse.ArgumentParser()
    opts.add_argument('input_path', help='Path to input files')
    opts.add_argument('reference_path', help='Path to reference folders')
    opts.add_argument('--dry_run', action='store_true', default=False,
                      help='just show the commands')
    options = opts.parse_args()
    if options.dry_run:
        run = dprint
    outdir, step1_dir = init()
    step1_jobIDs, step1_dir = step1_trimming(options.input_path, step1_dir)
    ac = step2_sorting(step1_jobIDs, step1_dir, options.reference_path, outdir)

