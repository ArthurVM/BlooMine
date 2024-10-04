from .utilities import *
from .moi import *

def runBlooMine(args, fq_box, flank_dict, f1f2_fasta, rundir):
    """ Run BlooMine using both flanks
    """

    for fqs in fq_box:
        prefix = fqs[0]
        read_files = fqs[1:]

        mkchdir(f"{rundir}/{prefix}")

        for flank_id, flank_fasta in flank_dict.items():

            mkchdir(f"{rundir}/{prefix}/{flank_id}")

            flank1_tmp, flank2_tmp = flank_fasta
            
            ## check if the MOI report exists
            if not checkRun(prefix, f"{rundir}/{prefix}/{flank_id}"):

                ## check if reads need decompressing
                read_files = [prepareFQ(f, f"{rundir}/{prefix}/{flank_id}") for f in read_files]
                fqs_joined = ",".join(read_files)

                print(
                    f""" 
                    BLOOMINE RUN STARTED:
                    fq ID : {prefix}
                    flank ID : {flank_id}
                    fq list  : {fqs_joined}\n"""
                    )
                
                ## run the first flank
                f1_runline = f"{args.bloomine_exec_bin}\
                -t {flank1_tmp} \
                -q {fqs_joined} \
                --prefix {prefix}.{flank_id}_F1 \
                --kmer {args.kmer} \
                --false_positive {args.false_positive} \
                --FP-sim {args.FP_sim} \
                --SP-error {args.SP_error} \
                --threads {args.threads}"
                
                if args.on_disk:
                    f1_runline += " --on-disk"

                subprocess.run(f1_runline + f" > {prefix}.{flank_id}.log", shell=True)

                if os.path.exists(f"{rundir}/{prefix}/{flank_id}/{prefix}.{flank_id}_F1_BMfiltered.fq"):
                    
                    ## run the second flank
                    f2_runline = f"{args.bloomine_exec_bin}\
                    -t {flank2_tmp} \
                    -q {prefix}.{flank_id}_F1_BMfiltered.fq \
                    --prefix {prefix}.{flank_id} \
                    --kmer {args.kmer} \
                    --false_positive {args.false_positive} \
                    --FP-sim {args.FP_sim} \
                    --SP-error {args.SP_error} \
                    --threads 1"

                    subprocess.run(f2_runline + f" >> {prefix}.{flank_id}.log", shell=True)
                    
                    outfq = f"{prefix}.{flank_id}_BMfiltered.fq"

                    if os.path.exists(outfq):
                        # run MOI
                        runMOI(prefix, outfq, f1f2_fasta)
                    
                    else:
                        print(f"Can't find {rundir}/{prefix}/{flank_id}/{prefix}.{flank_id}_F1_BMfiltered.fq\n")
                    
                else:
                    print(f"Can't find {rundir}/{prefix}/{flank_id}/{prefix}.{flank_id}_BMfiltered.fq\n")
            
            else:
                print(f"Found {rundir}/{prefix}/{flank_id}/{prefix}.report.txt. Skipping to prevent data overwrite.")
                continue
    
    print("Cleaning...")
    subprocess.run(f"rm {rundir}/{prefix}/{flank_id}/tmp*", shell=True)


def runMOI(prefix, fastq, target, min_anchor_kmer=11):
    target_dict = isolate_target(fastq, target, min_anchor_kmer)

    report_file = f"./{prefix}.report.txt"

    report_out(target_dict, fastq, target, report_file)