import json
import os
from collections import Counter

from Bio import SeqIO

from .utilities import *
from .moi import *
from .ProgressTable import ProgressTable
from .BloomineRunner import BloomineRunner
from .polyfamily import collect_probe_sets, bin_reads_by_probe, flank_score_records


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

            bm_run = BloomineRunner(self.args, read_files, fq_prefix, flank_id, flanks_fasta, self.baserundir, min_anchor_kmer=self.args.min_kmer)
            bm_run.run(p_table)

            ## update progress table
            p_table.update(fq_prefix, flank_id, "F1", "F", "green")
        
        if getattr(self.args, "polyfamily", False):
            self._run_polyfamily(fq_prefix)

        self._write_flank_score_log(fq_prefix)
        self._cleanRun(fq_prefix)


    def _run_polyfamily(self, fq_prefix):
        """
        Post-process reads across probe sets, assigning each read to the
        highest-scoring probe family and summarising STR variants to JSON.
        """
        sample_dir = os.path.join(self.baserundir, fq_prefix)

        all_reads, _, _ = collect_probe_sets(sample_dir)
        probe_bins = bin_reads_by_probe(all_reads)

        if not probe_bins:
            print(f"No polyfamily-eligible reads found for {fq_prefix}")
            return

        def norm_read_id(rid: str) -> str:
            return rid[1:] if rid.startswith("@") else rid

        ## map probe -> path of BMfiltered reads from flank 2 (final output)
        fastq_paths = {}
        for probe_id in self.flank_dict:
            runprefix = f"{fq_prefix}.{probe_id}"
            fq_path = os.path.join(sample_dir, probe_id, f"{runprefix}_BMfiltered.fq")
            if os.path.exists(fq_path):
                fastq_paths[probe_id] = fq_path

        ## build a lookup of read_id -> SeqRecord per probe for assigned reads
        assigned_ids = set(norm_read_id(read_id) for reads in probe_bins.values() for read_id, *_ in reads)
        read_store: Dict[str, Dict[str, SeqIO.SeqRecord]] = {}

        for probe_id, fq_path in fastq_paths.items():
            probe_reads = {}
            for rec in SeqIO.parse(fq_path, "fastq"):
                rid = norm_read_id(rec.description.strip())
                if rid in assigned_ids and rid not in probe_reads:
                    probe_reads[rid] = rec
            read_store[probe_id] = probe_reads

        probe_json = {}
        tmp_files = []

        for probe_id, reads in probe_bins.items():
            read_records = [read_store.get(probe_id, {}).get(norm_read_id(rid)) for rid, *_ in reads]
            read_records = [r for r in read_records if r is not None]
            if not read_records:
                continue

            tmp_fq = os.path.join(sample_dir, f"tmp.polyfamily.{probe_id}.fq")
            tmp_files.append(tmp_fq)
            SeqIO.write(read_records, tmp_fq, "fastq")

            flank_fasta = self.flank_dict.get(probe_id)
            if flank_fasta and os.path.exists(flank_fasta):
                target_dict = isolate_target(tmp_fq, flank_fasta, self.args.min_kmer)
                counts = Counter(str(seq) for seq in target_dict.values())
                probe_json[probe_id] = [[seq, count] for seq, count in sorted(counts.items(), key=lambda x: x[1], reverse=True)]

        ## write JSON summary if data is collected
        if probe_json:
            outfile = os.path.join(sample_dir, f"{fq_prefix}.polyfamily.json")
            with open(outfile, "w") as fh:
                json.dump(probe_json, fh, indent=2)
            print(f"Wrote polyfamily summary to {outfile}")

        ## clean up temp FASTQs
        for tmp in tmp_files:
            if os.path.exists(tmp):
                os.remove(tmp)

    def _write_flank_score_log(self, fq_prefix):
        """
        Write a log containing per-read flank scores and thresholds for each probe set.
        """
        sample_dir = os.path.join(self.baserundir, fq_prefix)
        outfile = os.path.join(sample_dir, f"{fq_prefix}.flank_scores.log")

        header = [
            "probe_set",
            "read_id",
            "flank_1_score",
            "flank_1_RC_score",
            "flank_2_score",
            "flank_2_RC_score",
            "threshold",
            "pass",
        ]

        with open(outfile, "w") as fh:
            fh.write("\t".join(header) + "\n")

            for probe_id in self.flank_dict:
                flank_dir = os.path.join(sample_dir, probe_id)

                runprefix = f"{fq_prefix}.{probe_id}"
                combined_log = os.path.join(flank_dir, f"{runprefix}_combined_flank_scores.tsv")
                if not os.path.exists(combined_log):
                    continue

                records = flank_score_records(combined_log)

                for read_id, vals in sorted(records.items()):
                    f1_score, f1_rc_score, f2_score, f2_rc_score, thr, pass_flag = vals

                    fh.write(
                        "\t".join([
                            str(probe_id),
                            str(read_id),
                            "" if f1_score is None else str(f1_score),
                            "" if f1_rc_score is None else str(f1_rc_score),
                            "" if f2_score is None else str(f2_score),
                            "" if f2_rc_score is None else str(f2_rc_score),
                            "" if thr is None else str(thr),
                            "" if pass_flag is None else str(pass_flag),
                        ]) + "\n"
                    )

        print(f"Wrote flank score log to {outfile}")


    def _cleanRun(self, fq_prefix):
        """ Clean up tmp files
        """
        fq_rundir = os.path.join(self.baserundir, fq_prefix)

        ## remove tmp files in the base directory
        for filename in os.listdir(fq_rundir):
            if filename.startswith("tmp"):
                os.remove(os.path.join(fq_rundir, filename))
