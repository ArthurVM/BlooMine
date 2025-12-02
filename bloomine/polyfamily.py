import argparse
import json
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

## read_id, sequence (as captured from log), flank1 score, flank2 score, total score
ReadTuple = Tuple[str, str, int, int, int]

## read_id -> (f1, f1_rc, f2, f2_rc, threshold, pass)
FlankRecord = Dict[str, Tuple[Optional[int], Optional[int], Optional[int], Optional[int], Optional[float], Optional[int]]]


def _parse_stdout(log_path: str) -> Tuple[Dict[int, Dict[str, Tuple[str, int]]], Dict[int, float]]:
    """
    Parse a BlooMine stdout log and return flank-specific scores and thresholds.

    Returns
    -------
    tuple
        ({1: {read_id: (seq, score)}, 2: {read_id: (seq, score)}}, {1: thr1, 2: thr2})
    """
    flank_scores: Dict[int, Dict[str, Tuple[str, int]]] = {1: {}, 2: {}}
    min_thresholds: Dict[int, float] = {1: 0.0, 2: 0.0}
    current_flank = None

    with open(log_path, "r") as handle:
        lines = handle.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        flank_match = re.match(r"Flank\s+(\d+)", line)
        if flank_match:
            current_flank = int(flank_match.group(1))
            i += 1
            continue

        if current_flank is None:
            i += 1
            continue

        thr_match = re.search(r"MINIMUM SCORE THRESHOLD SET TO\s*:\s*([0-9.]+)", line)
        if thr_match:
            min_thresholds[current_flank] = float(thr_match.group(1))
            i += 1
            continue

        if line.startswith("@") and "SCORE" not in line:
            read_id = line.lstrip("@").strip()
            seq_line = lines[i + 1].strip() if i + 1 < len(lines) else ""

            score_line = None
            j = i + 2
            while j < len(lines):
                candidate = lines[j]
                if "SCORE" in candidate:
                    score_line = candidate
                    j += 1
                    break
                if candidate.startswith("@") or candidate.strip().startswith("Flank"):
                    break
                j += 1

            score_match = (
                re.search(r"SCORE\s*:\s*(\d+)", score_line) if score_line else None
            )
            if score_match:
                score_val = int(score_match.group(1))
                existing = flank_scores[current_flank].get(read_id)
                if existing is None or score_val > existing[1]:
                    flank_scores[current_flank][read_id] = (seq_line, score_val)

            i = j
            continue

        i += 1

    return flank_scores, min_thresholds


def _combine_flanks(
    flank_scores: Dict[int, Dict[str, Tuple[str, int]]], thresholds: Dict[int, float]
) -> List[ReadTuple]:
    """Merge flank 1 and 2 calls into a list of tuples, filtered by thresholds."""
    combined: List[ReadTuple] = []
    th1 = thresholds.get(1, 0.0)
    th2 = thresholds.get(2, 0.0)
    read_ids = set(flank_scores[1]) & set(flank_scores[2])

    for read_id in sorted(read_ids):
        f1_seq, f1_score = flank_scores[1][read_id]
        f2_seq, f2_score = flank_scores[2][read_id]
        if f1_score < th1 or f2_score < th2:
            continue
        read_seq = f2_seq or f1_seq
        combined.append((read_id, read_seq, f1_score, f2_score, f1_score + f2_score))

    return combined


def collect_probe_sets(root_dir: str) -> Tuple[Dict[str, List[ReadTuple]], Dict[str, List[ReadTuple]], Dict[str, Dict[int, float]]]:
    """Walk the run directory, parse combined flank score TSVs, and return all tuples and maxima."""
    all_reads: Dict[str, List[ReadTuple]] = {}
    max_reads: Dict[str, List[ReadTuple]] = defaultdict(list)
    thresholds_by_probe: Dict[str, Dict[int, float]] = {}

    for dirpath, _, filenames in os.walk(root_dir):
        combined_files = [f for f in filenames if f.endswith("_combined_flank_scores.tsv")]
        if not combined_files:
            continue

        for fname in combined_files:
            log_path = os.path.join(dirpath, fname)
            probe_set = os.path.basename(dirpath)

            records = flank_score_records(log_path)
            tuples: List[ReadTuple] = []
            thresholds = {}

            for read_id, vals in records.items():
                f1, f1_rc, f2, f2_rc, thr, passed = vals

                ## take the best score per flank (forward vs RC)
                best_f1 = max([s for s in (f1, f1_rc) if s is not None], default=None)
                best_f2 = max([s for s in (f2, f2_rc) if s is not None], default=None)

                if thr is not None:
                    thresholds[1] = thr
                    thresholds[2] = thr

                ## skip if explicitly failed or missing flank scores
                if passed is not None and passed == 0:
                    continue
                if best_f1 is None or best_f2 is None:
                    continue

                tuples.append((read_id, "", best_f1, best_f2, best_f1 + best_f2))

            all_reads[probe_set] = tuples
            thresholds_by_probe[probe_set] = thresholds

            if tuples:
                max_sum = max(t[-1] for t in tuples)
                max_reads[probe_set] = [t for t in tuples if t[-1] == max_sum]

    return all_reads, max_reads, thresholds_by_probe


def choose_best_probes(all_reads: Dict[str, List[ReadTuple]]) -> Dict[str, str]:
    """Return a mapping of read_id -> probe with the highest total score."""
    best: Dict[str, Tuple[str, int]] = {}

    for probe_id, tuples in all_reads.items():
        for read_id, _seq, _s1, _s2, total in tuples:
            current = best.get(read_id)
            if current is None or total > current[1] or (total == current[1] and probe_id < current[0]):
                best[read_id] = (probe_id, total)

    return {read_id: probe for read_id, (probe, _score) in best.items()}


def bin_reads_by_probe(all_reads: Dict[str, List[ReadTuple]]) -> Dict[str, List[ReadTuple]]:
    """Filter to only the highest-scoring probe per read and return bins by probe."""
    assignments = choose_best_probes(all_reads)
    bins: Dict[str, List[ReadTuple]] = defaultdict(list)

    for probe_id, tuples in all_reads.items():
        for read in tuples:
            read_id = read[0]
            if assignments.get(read_id) == probe_id:
                bins[probe_id].append(read)

    return bins


def flank_score_records(log_path: str) -> FlankRecord:
    """
    Parse a combined flank score TSV and return per-read flank scores and flags.

    Returns
    -------
    records : {read_id: (f1, f1_rc, f2, f2_rc, threshold, pass_flag)}
    """
    records: FlankRecord = {}

    if not os.path.exists(log_path):
        return records

    with open(log_path, "r") as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            (
                read_id,
                f1_score,
                f1_rc_score,
                f2_score,
                f2_rc_score,
                threshold,
                passed,
            ) = parts + [None] * (7 - len(parts))

            def to_int(val):
                try:
                    return int(val)
                except (TypeError, ValueError):
                    return None

            def to_float(val):
                try:
                    return float(val)
                except (TypeError, ValueError):
                    return None

            records[read_id] = (
                to_int(f1_score),
                to_int(f1_rc_score),
                to_int(f2_score),
                to_int(f2_rc_score),
                to_float(threshold),
                to_int(passed),
            )

    return records
