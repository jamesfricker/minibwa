#!/usr/bin/env python3
"""Human-focused benchmark and QA harness for minibwa.

The default synthetic run is intentionally small enough for CI. Larger
HMF/GIAB runs use the same metric collection by passing a prepared reference,
reads, and optional BED stratifications.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ALT_RE = re.compile(r"(_alt$|_alt\b|GL\d+|KI\d+|HLA-ALT)", re.IGNORECASE)
HLA_RE = re.compile(r"(hla|chr6_gl.*_alt)", re.IGNORECASE)
DECOY_RE = re.compile(r"(decoy|hs37d5|hs38d1)", re.IGNORECASE)
RT_RE = re.compile(r"Real time:\s*([0-9.]+)\s*sec;\s*CPU:\s*([0-9.]+)\s*sec;\s*Peak RSS:\s*([0-9.]+)\s*GB")


@dataclass
class Alignment:
    name: str
    flag: int
    contig: str
    pos: int
    mapq: int

    @property
    def mapped(self) -> bool:
        return self.contig != "*" and not (self.flag & 4)

    @property
    def primary(self) -> bool:
        return not (self.flag & 0x100) and not (self.flag & 0x800)


@dataclass
class BedSet:
    name: str
    intervals: dict[str, list[tuple[int, int]]]


def run(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
    proc = subprocess.run(cmd, text=True, capture_output=True, **kwargs)
    if proc.returncode:
        sys.stderr.write(proc.stderr)
        sys.stderr.write(proc.stdout)
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc


def dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def write_fastq(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def make_synthetic_fixture(out_dir: Path) -> tuple[Path, Path, list[dict[str, object]], list[Path]]:
    rng = random.Random(27)
    chr1 = dna(rng, 1600)
    primary_target = dna(rng, 140)
    chr1 = chr1[:300] + primary_target + chr1[440:]

    hla_shared = dna(rng, 150)
    hla_primary = dna(rng, 500) + hla_shared + dna(rng, 500)
    hla_alt = dna(rng, 500) + hla_shared[:70] + "G" + hla_shared[71:] + dna(rng, 500)

    decoy = dna(rng, 400) + dna(rng, 140) + dna(rng, 400)
    decoy_target = decoy[400:540]

    duplicated = dna(rng, 110)
    chr1 = chr1[:900] + duplicated + chr1[1010:]
    hla_primary = hla_primary[:900] + duplicated + hla_primary[1010:]

    ref = out_dir / "synthetic-human.fa"
    reads = out_dir / "synthetic-human.fq"
    write_fasta(
        ref,
        [
            ("chr1", chr1),
            ("chr6", hla_primary),
            ("chr6_GL000250v2_alt", hla_alt),
            ("hs37d5_decoy", decoy),
        ],
    )
    write_fastq(
        reads,
        [
            ("primary_unique", primary_target),
            ("hla_primary_shared", hla_shared),
            ("hla_alt_near", hla_alt[500:650]),
            ("decoy_unique", decoy_target),
            ("low_mappability", duplicated),
        ],
    )

    (out_dir / "giab-high-confidence.bed").write_text("chr1\t300\t440\n", encoding="ascii")
    (out_dir / "low-mappability.bed").write_text("chr1\t900\t1010\nchr6\t900\t1010\n", encoding="ascii")
    (out_dir / "hla.bed").write_text("chr6\t500\t650\nchr6_GL000250v2_alt\t500\t650\n", encoding="ascii")
    (out_dir / "blacklist.bed").write_text("hs37d5_decoy\t400\t540\n", encoding="ascii")

    expectations = [
        {"read": "primary_unique", "contig": "chr1", "pos": 301, "max_dist": 0, "min_mapq": 30},
        {"read": "hla_primary_shared", "contig": "chr6", "pos": 501, "max_dist": 0, "max_mapq": 10},
        {"read": "hla_alt_near", "contig": "chr6_GL000250v2_alt", "pos": 501, "max_dist": 0},
        {"read": "decoy_unique", "contig": "hs37d5_decoy", "pos": 401, "max_dist": 0, "min_mapq": 30},
        {"read": "low_mappability", "max_mapq": 10},
    ]
    beds = [
        out_dir / "giab-high-confidence.bed",
        out_dir / "low-mappability.bed",
        out_dir / "hla.bed",
        out_dir / "blacklist.bed",
    ]
    return ref, reads, expectations, beds


def parse_sam(path: Path) -> list[Alignment]:
    out: list[Alignment] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            out.append(Alignment(fields[0], int(fields[1]), fields[2], int(fields[3]), int(fields[4])))
    return out


def parse_bed(path: Path) -> BedSet:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise SystemExit(f"bad BED line in {path}: {line}")
            intervals.setdefault(fields[0], []).append((int(fields[1]) + 1, int(fields[2])))
    return BedSet(path.stem, intervals)


def overlaps(aln: Alignment, bed: BedSet) -> bool:
    for start, end in bed.intervals.get(aln.contig, []):
        if start <= aln.pos <= end:
            return True
    return False


def mapq_bins(alignments: list[Alignment]) -> dict[str, int]:
    bins = {"0": 0, "1-9": 0, "10-29": 0, "30-59": 0, "60": 0}
    for aln in alignments:
        if not aln.mapped:
            continue
        if aln.mapq == 0:
            bins["0"] += 1
        elif aln.mapq < 10:
            bins["1-9"] += 1
        elif aln.mapq < 30:
            bins["10-29"] += 1
        elif aln.mapq < 60:
            bins["30-59"] += 1
        else:
            bins["60"] += 1
    return bins


def load_expectations(path: Path | None) -> list[dict[str, object]]:
    if path is None:
        return []
    expectations: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise SystemExit(f"bad expectation line in {path}: {line}")
            item: dict[str, object] = {"read": fields[0], "contig": fields[1], "pos": int(fields[2])}
            if len(fields) > 3 and fields[3]:
                item["max_dist"] = int(fields[3])
            expectations.append(item)
    return expectations


def index_file(prefix: Path) -> Path:
    return Path(str(prefix) + ".mbw")


def validate(
    metrics: dict[str, object], alignments: list[Alignment], expectations: list[dict[str, object]]
) -> list[str]:
    failures: list[str] = []
    primary_by_read = {aln.name: aln for aln in alignments if aln.primary}
    if float(metrics["alignment_rate"]) < 0.99:
        failures.append(f"alignment_rate below 0.99: {metrics['alignment_rate']}")
    for exp in expectations:
        read = str(exp["read"])
        aln = primary_by_read.get(read)
        if aln is None or not aln.mapped:
            failures.append(f"{read}: missing mapped primary alignment")
            continue
        contig = exp.get("contig")
        if contig is not None and aln.contig != contig:
            failures.append(f"{read}: expected contig {contig}, observed {aln.contig}")
        pos = exp.get("pos")
        if pos is not None and abs(aln.pos - int(pos)) > int(exp.get("max_dist", 0)):
            failures.append(f"{read}: expected position {pos}, observed {aln.pos}")
        if "min_mapq" in exp and aln.mapq < int(exp["min_mapq"]):
            failures.append(f"{read}: expected MAPQ >= {exp['min_mapq']}, observed {aln.mapq}")
        if "max_mapq" in exp and aln.mapq > int(exp["max_mapq"]):
            failures.append(f"{read}: expected MAPQ <= {exp['max_mapq']}, observed {aln.mapq}")
    return failures


def collect_metrics(
    alignments: list[Alignment],
    elapsed: float,
    stderr: str,
    beds: list[BedSet],
) -> dict[str, object]:
    primaries = [aln for aln in alignments if aln.primary]
    mapped = [aln for aln in primaries if aln.mapped]
    rt_match = RT_RE.search(stderr)
    metrics: dict[str, object] = {
        "reads": len(primaries),
        "mapped_reads": len(mapped),
        "alignment_rate": round(len(mapped) / len(primaries), 6) if primaries else 0.0,
        "wall_seconds": round(elapsed, 6),
        "minibwa_real_seconds": float(rt_match.group(1)) if rt_match else None,
        "minibwa_cpu_seconds": float(rt_match.group(2)) if rt_match else None,
        "peak_rss_gb": float(rt_match.group(3)) if rt_match else None,
        "mapq": mapq_bins(mapped),
        "alt_primaries": sum(1 for aln in mapped if ALT_RE.search(aln.contig)),
        "hla_primaries": sum(1 for aln in mapped if HLA_RE.search(aln.contig)),
        "decoy_primaries": sum(1 for aln in mapped if DECOY_RE.search(aln.contig)),
        "stratifications": {},
    }
    strat = metrics["stratifications"]
    assert isinstance(strat, dict)
    for bed in beds:
        strat[bed.name] = sum(1 for aln in mapped if overlaps(aln, bed))
    return metrics


def write_tsv(path: Path, metrics: dict[str, object]) -> None:
    mapq = metrics["mapq"]
    strat = metrics["stratifications"]
    assert isinstance(mapq, dict)
    assert isinstance(strat, dict)
    flat = {
        "reads": metrics["reads"],
        "mapped_reads": metrics["mapped_reads"],
        "alignment_rate": metrics["alignment_rate"],
        "wall_seconds": metrics["wall_seconds"],
        "minibwa_real_seconds": metrics["minibwa_real_seconds"],
        "minibwa_cpu_seconds": metrics["minibwa_cpu_seconds"],
        "peak_rss_gb": metrics["peak_rss_gb"],
        "mapq_0": mapq["0"],
        "mapq_1_9": mapq["1-9"],
        "mapq_10_29": mapq["10-29"],
        "mapq_30_59": mapq["30-59"],
        "mapq_60": mapq["60"],
        "alt_primaries": metrics["alt_primaries"],
        "hla_primaries": metrics["hla_primaries"],
        "decoy_primaries": metrics["decoy_primaries"],
    }
    flat.update({f"bed_{name}": count for name, count in strat.items()})
    with path.open("w") as handle:
        handle.write("\t".join(flat) + "\n")
        handle.write("\t".join("" if value is None else str(value) for value in flat.values()) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minibwa", default="./minibwa", help="minibwa executable")
    parser.add_argument("--out-dir", default=".context/human-benchmark", help="output directory")
    parser.add_argument("--reference", help="external FASTA for HMF/GIAB-style runs")
    parser.add_argument("--index-prefix", help="existing or output index prefix")
    parser.add_argument("--reads", help="single-end FASTA/FASTQ reads")
    parser.add_argument("--reads1", help="paired-end FASTA/FASTQ read 1")
    parser.add_argument("--reads2", help="paired-end FASTA/FASTQ read 2")
    parser.add_argument("--bed", action="append", default=[], help="BED stratification; may be repeated")
    parser.add_argument("--expect-primary-tsv", help="read, contig, 1-based position, optional max distance")
    parser.add_argument("--threads", default=os.environ.get("THREADS", "1"), help="mapping threads")
    parser.add_argument("--keep", action="store_true", help="keep an existing output directory")
    parser.add_argument("--no-validate", action="store_true", help="collect metrics without failing thresholds")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    if out_dir.exists() and not args.keep:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    synthetic = args.reference is None
    expectations = load_expectations(Path(args.expect_primary_tsv)) if args.expect_primary_tsv else []
    bed_paths = [Path(path) for path in args.bed]
    if synthetic:
        reference, reads, synthetic_expectations, synthetic_beds = make_synthetic_fixture(out_dir)
        expectations = synthetic_expectations
        bed_paths = synthetic_beds
    else:
        reference = Path(args.reference)
        if args.reads:
            reads = Path(args.reads)
        elif args.reads1 and args.reads2:
            reads = Path(args.reads1)
        else:
            raise SystemExit("external runs require --reads or --reads1/--reads2")

    index_prefix = Path(args.index_prefix) if args.index_prefix else out_dir / "index"
    if not index_file(index_prefix).exists():
        run([args.minibwa, "index", str(reference), str(index_prefix)])

    sam = out_dir / "aln.sam"
    cmd = [args.minibwa, "map", "-t", str(args.threads), "-x", "sr", str(index_prefix), str(reads)]
    if not synthetic and args.reads1 and args.reads2:
        cmd.append(str(Path(args.reads2)))
    start = time.monotonic()
    with sam.open("w") as sam_handle:
        proc = subprocess.run(cmd, text=True, stdout=sam_handle, stderr=subprocess.PIPE)
    elapsed = time.monotonic() - start
    if proc.returncode:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(cmd)}")
    (out_dir / "map.stderr").write_text(proc.stderr, encoding="ascii")

    alignments = parse_sam(sam)
    beds = [parse_bed(path) for path in bed_paths]
    metrics = collect_metrics(alignments, elapsed, proc.stderr, beds)
    failures = [] if args.no_validate else validate(metrics, alignments, expectations)
    metrics["validation_failures"] = failures

    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="ascii")
    write_tsv(out_dir / "metrics.tsv", metrics)
    print(json.dumps(metrics, indent=2, sort_keys=True))

    if failures:
        for failure in failures:
            print(f"validation failure: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
