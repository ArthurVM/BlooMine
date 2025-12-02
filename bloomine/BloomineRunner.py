from .utilities import *
from .moi import *
from .runCommands import *


class BloomineRunner:
    """
    A class to manage a BlooMine run and track progress.
    """


    def __init__(self, args, fqbox, fq_prefix, flank_prefix, flank_fasta, baserundir, min_anchor_kmer=11):

        ## initialise base attributes
        self.args = args
        self.fq_box = fqbox
        self.fqprefix = fq_prefix
        self.runprefix = f"{fq_prefix}.{flank_prefix}"
        self.flank_fasta = flank_fasta
        self.baserundir = baserundir
        self.flank_prefix = flank_prefix

        ## initialise output directory and log files
        self.rundir = f"{self.baserundir}/{self.fqprefix}/{self.flank_prefix}"
        mkchdir(self.rundir)
        self.stdout = open(f"{self.rundir}/stdout.txt", 'w')
        self.stderr = open(f"{self.rundir}/stderr.txt", 'w')

        ## process flanks
        f1_fasta, f2_fasta = self._processFlankFasta()
        self.f1_fasta = f1_fasta
        self.f2_fasta = f2_fasta

        ## initialise BlooMine argument parameters
        self.kmer = int(self.args.kmer)
        self.false_positive = float(self.args.false_positive)
        self.FP_sim = float(self.args.FP_sim)
        self.SP_error = float(self.args.SP_error)
        self.min_anchor_kmer = min_anchor_kmer

        ## initialise output files
        self.F1_fqout_prefix = f"{self.runprefix}_F1"
        self.F2_fqout_prefix = f"{self.runprefix}"
        self.moi_report = f"{self.rundir}/{self.runprefix}.report.txt"
        self.combined_log = f"{self.rundir}/{self.runprefix}_combined_flank_scores.tsv"


    def run(self, p_table):
        
        # print(
        #     f"\tBLOOMINE RUN STARTED:\n\tfq ID    : {self.fqprefix}\n\tflank ID : {self.flank_prefix}\n\tfq list  : {fqjoined}\n\trun dir  : {self.rundir}\n"""
        #     )
        
        flank1_fq_flag, flank2_fq_flag, moi_report_flag = self._checkRun()


        ###### Screen first flank ######
        runcode = "Flank 1"

        if not flank1_fq_flag:

            ## only prepare and decompress FASTQ files if Flank 1 screening is required
            read_files = [ prepareFQ(f, f"{self.baserundir}/{self.fqprefix}") for f in self.fq_box ]
            fqjoined = ",".join(read_files)

            p_table.update(self.fqprefix, self.flank_prefix, runcode, "...", "yellow")   ## update progress table as running
            returncode = self._execBloomine( runcode, self.f1_fasta, fqjoined, self.F1_fqout_prefix, self.args.threads, self.args.on_disk, 1 )

            col = "red" if returncode != 0 else "green"
            message = "Finished" if returncode != 0 else returncode
            p_table.update(self.fqprefix, self.flank_prefix, runcode, message, col)  ## update progress table
        
        else:
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "found", "green")      ## update progress table as found

        f1_out = f"{self.rundir}/{self.F1_fqout_prefix}_BMfiltered.fq"


        ###### Screen second flank ######
        runcode = "Flank 2"

        if not flank2_fq_flag and os.path.exists(f1_out):
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "...", "yellow")   ## update progress table as running
            returncode = self._execBloomine( runcode, self.f2_fasta, f1_out, self.F2_fqout_prefix, 1, False, 2 )
            
            col = "red" if returncode != 0 else "green"
            message = "Finished" if returncode != 0 else returncode
            p_table.update(self.fqprefix, self.flank_prefix, runcode, message, col)  ## update progress table
        
        elif os.path.exists(f1_out):
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "found", "green")      ## update progress table as found


        f2_out = f"{self.rundir}/{self.F2_fqout_prefix}_BMfiltered.fq"

        # Merge flank logs into a single combined TSV for this probe set
        self._merge_flank_logs()


        ###### MOI analysis ######
        runcode = "MOI"
        
        if not self.args.skip_moi:
            if not moi_report_flag and os.path.exists(f2_out):
                p_table.update(self.fqprefix, self.flank_prefix, runcode, "...", "yellow")  ## update progress table as running
                self._runMOI(self.runprefix, f2_out)
                p_table.update(self.fqprefix, self.flank_prefix, runcode, "Finished", "green")     ## update progress table as success

            elif os.path.exists(f2_out):
                p_table.update(self.fqprefix, self.flank_prefix, runcode, "found", "green")      ## update progress table as found


        ###### cleanup ######
        self._cleanRun()

    
    def _execBloomine(self, run_ID, flank, fq, outprefix, threads, on_disk, flank_number):
        """
        Run BlooMine for a single flank screen.
        """
        f1_runline = f"{self.args.bloomine_exec_bin}\
        -t {flank} \
        -q {fq} \
        --prefix {outprefix} \
        --kmer {self.kmer} \
        --false_positive {self.false_positive} \
        --FP-sim {self.FP_sim} \
        --SP-error {self.SP_error} \
        --flank_number {flank_number} \
        --threads {threads}"

        if on_disk:
            f1_runline += " --on-disk"

        returncode, stout, stderr = command(f1_runline, run_ID).run_comm_quiet(1, self.stdout, self.stderr, exit=0)

        return returncode
    

    def _runMOI(self, prefix, fastq):
        """ Run MOI analysis
        """
        target_dict = isolate_target(fastq, self.flank_fasta, self.min_anchor_kmer)

        report_out(target_dict, fastq, self.flank_fasta, self.moi_report)
    

    def _processFlankFasta(self):
        """ Take the flank pair multifasta and split into individual fastas
        """
        flanks = list(SeqIO.parse(self.flank_fasta, "fasta"))

        if len(flanks) != 2:
            raise ValueError(f"Expected 2 flank sequences in {self.flank_fasta}, found {len(flanks)}")

        flank1, flank2 = flanks

        flank1_path = f"{self.baserundir}/tmp.{self.flank_prefix}_f1.fa"
        flank2_path = f"{self.baserundir}/tmp.{self.flank_prefix}_f2.fa"

        with open(flank1_path, "w") as f1out:
            f1out.write(f">{flank1.description}\n{flank1.seq}\n")

        with open(flank2_path, "w") as f2out:
            f2out.write(f">{flank2.description}\n{flank2.seq}\n")

        return flank1_path, flank2_path
    

    def _checkRun(self):
        """ Check the output from previous BlooMine runs
        """

        flank1_fq_flag = False
        flank2_fq_flag = False
        moi_report_flag = False

        if os.path.exists(f"{self.rundir}/{self.F1_fqout_prefix}_BMfiltered.fq"):
            # print(f"\t\tFound {self.F1_fqout_prefix}_BMfiltered.fq.\n\t\tStarting from flank 2 to prevent data overwrite.\n")
            flank1_fq_flag = True
        
        if os.path.exists(f"{self.rundir}/{self.F2_fqout_prefix}_BMfiltered.fq"):
            # print(f"\t\tFound {self.F2_fqout_prefix}_BMfiltered.fq.\n\t\tStarting from MOI analysis to prevent data overwrite.\n")
            flank2_fq_flag = True

        if os.path.exists(self.moi_report):
            # print(f"\t\tFound {self.moi_report}.\n\t\tSkipping to prevent data overwrite.\n")
            moi_report_flag = True

        return flank1_fq_flag, flank2_fq_flag, moi_report_flag


    def _cleanRun(self):
        """ Clean up tmp files
        """
        ## remove tmp files in the run directory
        for filename in os.listdir(self.rundir):
            if filename.startswith("tmp"):
                os.remove(os.path.join(self.rundir, filename))

        self.stdout.close()
        self.stderr.close()


    def _merge_flank_logs(self):
        """Merge flank 1 and flank 2 score logs into a single TSV."""
        f1_log = f"{self.rundir}/{self.F1_fqout_prefix}_flank_scores.tsv"
        f2_log = f"{self.rundir}/{self.F2_fqout_prefix}_flank_scores.tsv"

        if not (os.path.exists(f1_log) and os.path.exists(f2_log)):
            return

        def parse(log_path):
            scores = {}
            threshold = None
            with open(log_path) as fh:
                fh.readline()
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    read_id, flank, score, thr, rc, passed = parts
                    try:
                        flank = int(flank)
                        score_val = int(score)
                        rc_val = int(rc)
                        thr_val = float(thr)
                    except ValueError:
                        continue
                    threshold = thr_val if threshold is None else threshold
                    key = (read_id, rc_val)
                    scores.setdefault(key, {}).setdefault(flank, 0)
                    scores[(read_id, rc_val)][flank] = max(scores[(read_id, rc_val)][flank], score_val)
            return scores, threshold

        f1_scores, f1_thr = parse(f1_log)
        f2_scores, f2_thr = parse(f2_log)

        read_ids = set([rid for rid, _ in f1_scores.keys()] + [rid for rid, _ in f2_scores.keys()])
        with open(self.combined_log, "w") as out:
            out.write("\t".join([
                "read_id",
                "flank_1_score",
                "flank_1_RC_score",
                "flank_2_score",
                "flank_2_RC_score",
                "threshold",
                "pass"
            ]) + "\n")

            threshold = max([t for t in [f1_thr, f2_thr] if t is not None], default=None)

            for rid in sorted(read_ids):
                f1_forward = f1_scores.get((rid, 0), {}).get(1)
                f1_rc = f1_scores.get((rid, 1), {}).get(1)
                f2_forward = f2_scores.get((rid, 0), {}).get(2)
                f2_rc = f2_scores.get((rid, 1), {}).get(2)

                f1_best = max([s for s in [f1_forward, f1_rc] if s is not None], default=None)
                f2_best = max([s for s in [f2_forward, f2_rc] if s is not None], default=None)

                pass_flag = None
                if f1_thr is not None and f2_thr is not None:
                    pass_flag = int((f1_best is not None and f1_best >= f1_thr) and (f2_best is not None and f2_best >= f2_thr))

                out.write("\t".join([
                    rid,
                    "" if f1_forward is None else str(f1_forward),
                    "" if f1_rc is None else str(f1_rc),
                    "" if f2_forward is None else str(f2_forward),
                    "" if f2_rc is None else str(f2_rc),
                    "" if threshold is None else str(threshold),
                    "" if pass_flag is None else str(pass_flag)
                ]) + "\n")
