#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from collections import OrderedDict, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "waterbird_metadata_nanostats_summary.tsv"

METADATA_FILE = ROOT / "metadata-15757625-processed-ok.tsv"
WATERBIRD_FILES = [ROOT / "waterbird1.tsv", ROOT / "waterbird2.tsv", ROOT / "waterbird3.tsv"]

METADATA_FIELDS = [
    "accession",
    "study",
    "object_status",
    "bioproject_accession",
    "biosample_accession",
    "library_id",
    "title",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "metadata_platform",
    "instrument_model",
    "design_description",
    "filetype",
    "fasta_file",
    "assembly",
    "source_filenames",
    "project_id",
    "isolation_source",
    "source_type",
    "host",
    "collection_date",
    "collected_by",
    "locality",
    "preservation",
    "waterbird_platform",
    "sequencing_approach",
    "barcode",
    "run_date",
    "project_directory",
    "comment",
    "if_repeated",
    "specimenID",
]

STATS_KEYS = [
    "number_of_reads",
    "number_of_bases",
    "median_read_length",
    "mean_read_length",
    "read_length_stdev",
    "n50",
    "mean_qual",
    "median_qual",
    "longest_read_with_Q_1",
    "longest_read_with_Q_2",
    "longest_read_with_Q_3",
    "longest_read_with_Q_4",
    "longest_read_with_Q_5",
    "highest_Q_read_with_length_1",
    "highest_Q_read_with_length_2",
    "highest_Q_read_with_length_3",
    "highest_Q_read_with_length_4",
    "highest_Q_read_with_length_5",
    "reads_gt_Q5",
    "reads_gt_Q7",
    "reads_gt_Q10",
    "reads_gt_Q12",
    "reads_gt_Q15",
]

STATS_HEADER = [f"raw_{key}" for key in STATS_KEYS] + [f"filt_{key}" for key in STATS_KEYS]


def merge_value(record: dict[str, str], field: str, value: str) -> None:
    value = (value or "").strip()
    if not value:
        return
    current = record.get(field, "")
    if not current:
        record[field] = value
        return
    if current == value:
        return
    existing = [piece.strip() for piece in current.split(" | ") if piece.strip()]
    if value not in existing:
        record[field] = current + " | " + value


def add_source(record: dict[str, str], source_name: str) -> None:
    merge_value(record, "source_files", source_name)


def load_metadata() -> tuple[OrderedDict[str, dict[str, str]], list[str]]:
    records: OrderedDict[str, dict[str, str]] = OrderedDict()
    order: list[str] = []

    def get_record(sample_id: str) -> dict[str, str]:
        if sample_id not in records:
            records[sample_id] = {"sample_id": sample_id}
            order.append(sample_id)
        return records[sample_id]

    with METADATA_FILE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = (row.get("sample_name") or "").strip()
            if not sample_id:
                continue
            record = get_record(sample_id)
            add_source(record, METADATA_FILE.name)

            merge_value(record, "accession", row.get("accession", ""))
            merge_value(record, "study", row.get("study", ""))
            merge_value(record, "object_status", row.get("object_status", ""))
            merge_value(record, "bioproject_accession", row.get("bioproject_accession", ""))
            merge_value(record, "biosample_accession", row.get("biosample_accession", ""))
            merge_value(record, "library_id", row.get("library_ID", ""))
            merge_value(record, "title", row.get("title", ""))
            merge_value(record, "library_strategy", row.get("library_strategy", ""))
            merge_value(record, "library_source", row.get("library_source", ""))
            merge_value(record, "library_selection", row.get("library_selection", ""))
            merge_value(record, "library_layout", row.get("library_layout", ""))
            merge_value(record, "metadata_platform", row.get("platform", ""))
            merge_value(record, "instrument_model", row.get("instrument_model", ""))
            merge_value(record, "design_description", row.get("design_description", ""))
            merge_value(record, "filetype", row.get("filetype", ""))
            merge_value(record, "fasta_file", row.get("fasta_file", ""))
            merge_value(record, "assembly", row.get("assembly", ""))

            filenames = [row.get(f"filename{i}", "") if i else row.get("filename", "") for i in range(9)]
            filenames = [value.strip() for value in filenames if value and value.strip()]
            if filenames:
                merge_value(record, "source_filenames", "; ".join(filenames))

    return records, order


def merge_waterbird(records: OrderedDict[str, dict[str, str]], order: list[str]) -> None:
    for path in WATERBIRD_FILES:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample_id = (row.get("sampleID") or "").strip()
                if not sample_id:
                    continue
                if sample_id not in records:
                    records[sample_id] = {"sample_id": sample_id}
                    order.append(sample_id)
                record = records[sample_id]
                add_source(record, path.name)

                merge_value(record, "project_id", row.get("projectID", ""))
                merge_value(record, "isolation_source", row.get("isolation_source", ""))
                merge_value(record, "source_type", row.get("source_type", ""))
                merge_value(record, "host", row.get("host", ""))
                merge_value(record, "collection_date", row.get("collection_date", ""))
                merge_value(record, "collected_by", row.get("collected_by", ""))
                merge_value(record, "locality", row.get("locality", ""))
                merge_value(record, "preservation", row.get("preservation", ""))
                merge_value(record, "waterbird_platform", row.get("platform", ""))
                merge_value(record, "sequencing_approach", row.get("sequencing_approach", ""))
                merge_value(record, "barcode", row.get("barcode", ""))
                merge_value(record, "run_date", row.get("run_date", ""))
                merge_value(record, "project_directory", row.get("project_directory", ""))
                merge_value(record, "comment", row.get("comment", ""))
                merge_value(record, "if_repeated", row.get("if_repeated", ""))
                merge_value(record, "specimenID", row.get("specimenID", ""))


def normalize_nanostats_key(raw_key: str) -> str | None:
    key = raw_key.strip()
    mapping = {
        "number_of_reads": "number_of_reads",
        "number_of_bases": "number_of_bases",
        "median_read_length": "median_read_length",
        "mean_read_length": "mean_read_length",
        "read_length_stdev": "read_length_stdev",
        "n50": "n50",
        "mean_qual": "mean_qual",
        "median_qual": "median_qual",
        "longest_read_(with_Q):1": "longest_read_with_Q_1",
        "longest_read_(with_Q):2": "longest_read_with_Q_2",
        "longest_read_(with_Q):3": "longest_read_with_Q_3",
        "longest_read_(with_Q):4": "longest_read_with_Q_4",
        "longest_read_(with_Q):5": "longest_read_with_Q_5",
        "highest_Q_read_(with_length):1": "highest_Q_read_with_length_1",
        "highest_Q_read_(with_length):2": "highest_Q_read_with_length_2",
        "highest_Q_read_(with_length):3": "highest_Q_read_with_length_3",
        "highest_Q_read_(with_length):4": "highest_Q_read_with_length_4",
        "highest_Q_read_(with_length):5": "highest_Q_read_with_length_5",
        "Reads >Q5:": "reads_gt_Q5",
        "Reads >Q7:": "reads_gt_Q7",
        "Reads >Q10:": "reads_gt_Q10",
        "Reads >Q12:": "reads_gt_Q12",
        "Reads >Q15:": "reads_gt_Q15",
    }
    return mapping.get(key)


def parse_nanostats(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            key = normalize_nanostats_key(row[0])
            if key:
                values[key] = row[1].strip()
    return values


def resolve_sample_dirs() -> list[tuple[str, str, Path, Path]]:
    pattern = re.compile(r"^(?P<sample>.+?)_(?P<run>waterbird_\d+)_(?P<kind>raw|filt)$")
    rows: list[tuple[str, str, Path, Path]] = []
    for raw_dir in sorted(ROOT.glob("*_raw")):
        match = pattern.match(raw_dir.name)
        if not match:
            continue
        sample_id = match.group("sample")
        run_directory = match.group("run")
        filt_dir = ROOT / f"{sample_id}_{run_directory}_filt"
        rows.append((sample_id, run_directory, raw_dir, filt_dir))
    return rows


def main() -> None:
    records, order = load_metadata()
    merge_waterbird(records, order)

    fieldnames = ["sample_id", "source_files", *METADATA_FIELDS, "output_run_directory", "raw_dir", "filt_dir", *STATS_HEADER]

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()

        for sample_id, run_directory, raw_dir, filt_dir in resolve_sample_dirs():
            record = records.setdefault(sample_id, {"sample_id": sample_id})
            if sample_id not in order:
                order.append(sample_id)
            row = {name: record.get(name, "") for name in fieldnames}
            row["sample_id"] = sample_id
            row["output_run_directory"] = run_directory
            row["raw_dir"] = raw_dir.name
            row["filt_dir"] = filt_dir.name

            raw_stats = parse_nanostats(raw_dir / "NanoStats.txt")
            filt_stats = parse_nanostats(filt_dir / "NanoStats.txt")
            for key in STATS_KEYS:
                row[f"raw_{key}"] = raw_stats.get(key, "")
                row[f"filt_{key}"] = filt_stats.get(key, "")

            writer.writerow(row)

    print(f"Wrote {OUTPUT.name}")


if __name__ == "__main__":
    main()