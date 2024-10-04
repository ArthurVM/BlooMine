from .utilities import *
from .moi import *
from .ProgressTable import ProgressTable
from .BloomineRunner import BloomineRunner


class RunManager:
    """
    A class to manage the BlooMine run and track progress.
    """

    def __init__(self, args, fq_box, flank_dict, baserundir):
        self.args = args
        self.fq_box = fq_box
        self.flank_dict = flank_dict
        self.baserundir = baserundir


    def run(self):
        """
        Run BlooMine for all flanks and samples.
        """

        p_table = ProgressTable([fqs[0] for fqs in self.fq_box], list(self.flank_dict.keys()))

        for fqs in self.fq_box:
            self._run_sample(fqs, p_table)

        p_table.close()


    def _run_sample(self, fqs, p_table):
        """
        Run BlooMine for a specific sample.
        """
        fq_prefix = fqs[0]
        read_files = fqs[1:]

        ## make run directory
        mkchdir(f"{self.baserundir}/{fq_prefix}")

        for flank_id, flanks_fasta in self.flank_dict.items():

            bm_run = BloomineRunner(self.args, read_files, fq_prefix, flank_id, flanks_fasta, self.baserundir)
            bm_run.run(p_table)

            ## update progress table
            p_table.update(fq_prefix, flank_id, "F1", "F", "green")
        
        self._cleanRun(fq_prefix)


    def _cleanRun(self, fq_prefix):
        """ Clean up tmp files
        """
        fq_rundir = os.path.join(self.baserundir, fq_prefix)

        ## Remove tmp files in the base directory
        for filename in os.listdir(fq_rundir):
            if filename.startswith("tmp"):
                os.remove(os.path.join(fq_rundir, filename))
