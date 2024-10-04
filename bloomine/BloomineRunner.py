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
            returncode = self._execBloomine( runcode, self.f1_fasta, fqjoined, self.F1_fqout_prefix, self.args.threads, self.args.on_disk )

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
            returncode = self._execBloomine( runcode, self.f2_fasta, f1_out, self.F2_fqout_prefix, 1, False )
            
            col = "red" if returncode != 0 else "green"
            message = "Finished" if returncode != 0 else returncode
            p_table.update(self.fqprefix, self.flank_prefix, runcode, message, col)  ## update progress table
        
        elif os.path.exists(f1_out):
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "found", "green")      ## update progress table as found


        f2_out = f"{self.rundir}/{self.F2_fqout_prefix}_BMfiltered.fq"


        ###### MOI analysis ######
        runcode = "MOI"

        if not moi_report_flag and os.path.exists(f2_out):
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "...", "yellow")  ## update progress table as running
            self._runMOI(self.runprefix, f2_out)
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "Finished", "green")     ## update progress table as success

        elif os.path.exists(f2_out):
            p_table.update(self.fqprefix, self.flank_prefix, runcode, "found", "green")      ## update progress table as found


        ###### cleanup ######
        self._cleanRun()

    
    def _execBloomine(self, run_ID, flank, fq, outprefix, threads, on_disk):
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
        for rec in SeqIO.parse(self.flank_fasta, "fasta"):
            with open(f"{self.baserundir}/tmp.{rec.id}_f1.fa", "w") as f1out, open(f"{self.baserundir}/tmp.{rec.id}_f2.fa", "w") as f2out:
                f1out.write(f">{rec.id}\n{rec.seq}")
                f2out.write(f">{rec.id}\n{rec.seq}")

        return f"{self.baserundir}/tmp.{rec.id}_f1.fa", f"{self.baserundir}/tmp.{rec.id}_f2.fa"
    

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
        ## Remove tmp files in the run directory
        for filename in os.listdir(self.rundir):
            if filename.startswith("tmp"):
                os.remove(os.path.join(self.rundir, filename))

        self.stdout.close()
        self.stderr.close()