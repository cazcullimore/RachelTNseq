#!/usr/bin/env python3
"""
gff_to_ptt.py

Convert a Geneious-exported GFF3 annotation plus matching FASTA into
legacy NCBI PTT files, one PTT per contig/replicon.

Designed specifically to handle Geneious GFF files where annotations
look like:

    contig_4  Geneious  gene  1  978  .  +  .  Name=repB gene
    contig_4  Geneious  CDS   1  978  .  +  0  Name=repB CDS

and also cases where the gene feature has no Name:

    contig_4  Geneious  gene  1054 1173 . - . .
    contig_4  Geneious  CDS   1054 1173 . - 0 Name=hypothetical protein CDS

Important behavior
------------------
* PTT entries are generated from CDS features only.
* A gene feature with identical coordinates/strand is used to recover
  the gene symbol where possible.
* Repeated gene names are NOT merged.
* Each CDS receives a unique generated Synonym/locus ID.
* Geneious suffixes such as " gene" and " CDS" are removed.
* URL-encoded GFF values such as "%2C" are decoded.
* FASTA sequence IDs must normally match GFF column 1.
* GFF ##sequence-region lengths are also read and used as a fallback.

Usage
-----
python3 gff_to_ptt.py \
    -g annotation.gff3 \
    -f genome.fasta \
    -o ptt_out

Output example
--------------
ptt_out/contig_1.ptt
ptt_out/contig_2.ptt
ptt_out/contig_3.ptt
ptt_out/contig_4.ptt
"""

import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote


def open_maybe_gz(path):
    """Open plain-text or gzip-compressed file."""
    path = str(path)

    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")

    return open(path, "rt", encoding="utf-8", errors="replace")


def clean_text(value):
    """
    Clean a value before writing it to a tab-delimited PTT file.
    """
    if value is None:
        return ""

    value = unquote(str(value))
    value = value.replace("\t", " ")
    value = value.replace("\r", " ")
    value = value.replace("\n", " ")
    return value.strip()


def strip_geneious_suffix(value, feature_type=None):
    """
    Remove Geneious feature-type suffixes from Name fields.

    Examples:
        repB gene                -> repB
        repB CDS                 -> repB
        hypothetical protein CDS -> hypothetical protein
        copA ncRNA               -> copA
    """
    value = clean_text(value)

    if not value:
        return ""

    suffixes = []

    if feature_type:
        suffixes.append(" " + feature_type)

    suffixes.extend(
        [
            " gene",
            " CDS",
            " ncRNA",
            " tRNA",
            " rRNA",
            " tmRNA",
        ]
    )

    # Longest suffix first.
    for suffix in sorted(set(suffixes), key=len, reverse=True):
        if value.lower().endswith(suffix.lower()):
            value = value[: -len(suffix)].rstrip()
            break

    return value


def parse_attributes(attr_string):
    """
    Parse GFF3 or simple GTF-style attributes.

    Geneious in this file primarily uses:
        Name=repB gene

    but the parser also accepts conventional GFF attributes such as:
        ID=...
        locus_tag=...
        gene=...
        product=...
        protein_id=...
    """
    attrs = {}

    attr_string = attr_string.strip()

    if not attr_string or attr_string == ".":
        return attrs

    for field in attr_string.split(";"):
        field = field.strip()

        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key.strip()] = clean_text(value)

        else:
            # Basic support for:
            # gene_id "abc"
            match = re.match(r'(\S+)\s+"?([^"]+)"?', field)

            if match:
                attrs[match.group(1)] = clean_text(match.group(2))

    return attrs


def parse_fasta_lengths(fasta_path):
    """
    Return:
        {seqid: sequence_length}

    The seqid is the first whitespace-delimited token after ">".
    """
    lengths = {}

    current_name = None
    current_length = 0

    with open_maybe_gz(fasta_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")

            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_length

                header = line[1:].strip()

                if not header:
                    raise ValueError("Encountered an empty FASTA header.")

                current_name = header.split()[0]
                current_length = 0

            else:
                if current_name is None:
                    if not line.strip():
                        continue

                    raise ValueError(
                        "FASTA sequence data was encountered before the first header."
                    )

                current_length += len(line.strip())

    if current_name is not None:
        lengths[current_name] = current_length

    return lengths


def read_gff(gff_path):
    """
    Read the GFF.

    Returns:
        genes_by_coordinate
        cds_by_seqid
        gff_sequence_lengths
        feature_counts
        attribute_counts

    genes_by_coordinate is keyed by:
        (seqid, start, end, strand)

    This is critical for this Geneious file because gene symbols are NOT
    guaranteed to be unique. Two copies of a gene can legitimately have
    the same Name at different positions.
    """
    genes_by_coordinate = {}
    cds_by_seqid = defaultdict(list)
    gff_sequence_lengths = {}

    feature_counts = Counter()
    attribute_counts = Counter()

    with open_maybe_gz(gff_path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            # Read Geneious sequence-region declarations.
            if line.startswith("##sequence-region"):
                parts = line.split()

                if len(parts) >= 4:
                    seqid = parts[1]

                    try:
                        region_start = int(parts[2])
                        region_end = int(parts[3])
                    except ValueError:
                        continue

                    # Usually 1..length.
                    gff_sequence_lengths[seqid] = (
                        region_end - region_start + 1
                    )

                continue

            if line.startswith("#"):
                continue

            columns = line.split("\t")

            if len(columns) < 9:
                # prints even for normal fasta in end of gff
                # print(
                #     f"WARNING: skipping malformed GFF line {line_number}: "
                #     f"expected 9 columns, found {len(columns)}",
                #     file=sys.stderr,
                # )
                continue

            seqid = columns[0]
            source = columns[1]
            feature_type = columns[2]
            start_text = columns[3]
            end_text = columns[4]
            score = columns[5]
            strand = columns[6]
            phase = columns[7]
            attr_string = columns[8]

            feature_counts[feature_type] += 1

            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                print(
                    f"WARNING: invalid coordinates on GFF line {line_number}",
                    file=sys.stderr,
                )
                continue

            attrs = parse_attributes(attr_string)

            for key in attrs:
                attribute_counts[key] += 1

            coordinate_key = (seqid, start, end, strand)

            # Save gene features so that CDS features can be paired by
            # exact coordinates rather than by gene name.
            if feature_type == "gene":
                raw_name = attrs.get("Name", "")

                gene_name = (
                    attrs.get("gene")
                    or attrs.get("locus_tag")
                    or attrs.get("gene_id")
                    or strip_geneious_suffix(raw_name, "gene")
                )

                genes_by_coordinate[coordinate_key] = {
                    "name": clean_text(gene_name),
                    "attrs": attrs,
                    "line_number": line_number,
                }

                continue

            # PTT is a protein table, so only CDS features become records.
            if feature_type != "CDS":
                continue

            cds_by_seqid[seqid].append(
                {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "phase": phase,
                    "source": source,
                    "attrs": attrs,
                    "line_number": line_number,
                }
            )

    return (
        genes_by_coordinate,
        cds_by_seqid,
        gff_sequence_lengths,
        feature_counts,
        attribute_counts,
    )


def looks_like_gene_symbol(name):
    """
    Conservative heuristic used only when no matching named 'gene'
    feature is available.

    This is intentionally conservative. Descriptions such as
    'hypothetical protein' should NOT become gene symbols.
    """
    if not name:
        return False

    if " " in name:
        return False

    if len(name) > 30:
        return False

    # Permit normal microbial gene-name characters.
    if not re.fullmatch(r"[A-Za-z0-9_.:+-]+", name):
        return False

    return True


def build_ptt_records(genes_by_coordinate, cds_by_seqid):
    """
    Convert CDS records into PTT-ready dictionaries.

    Geneious does not provide locus_tag in this GFF, so Synonym is
    generated as:

        contig_1_000001
        contig_1_000002
        ...

    This avoids collisions when the same gene name occurs multiple times.
    """
    result = {}

    for seqid, cds_list in cds_by_seqid.items():
        cds_list = sorted(
            cds_list,
            key=lambda record: (
                record["start"],
                record["end"],
                record["strand"],
                record["line_number"],
            ),
        )

        records = []

        for index, cds in enumerate(cds_list, 1):
            attrs = cds["attrs"]

            coordinate_key = (
                seqid,
                cds["start"],
                cds["end"],
                cds["strand"],
            )

            paired_gene = genes_by_coordinate.get(coordinate_key)

            cds_name = strip_geneious_suffix(
                attrs.get("Name", ""),
                "CDS",
            )

            # Prefer conventional attributes if they exist.
            conventional_gene = clean_text(
                attrs.get("gene")
                or attrs.get("gene_name")
                or ""
            )

            locus_tag = clean_text(
                attrs.get("locus_tag")
                or attrs.get("gene_id")
                or ""
            )

            protein_id = clean_text(
                attrs.get("protein_id")
                or attrs.get("proteinId")
                or ""
            )

            product = clean_text(attrs.get("product") or "")

            paired_gene_name = ""

            if paired_gene:
                paired_gene_name = clean_text(paired_gene["name"])

            # Determine gene symbol.
            if conventional_gene:
                gene_name = conventional_gene

            elif paired_gene_name:
                gene_name = paired_gene_name

            elif looks_like_gene_symbol(cds_name):
                # If Geneious supplied only the CDS Name and it looks
                # convincingly like a gene symbol, use it.
                gene_name = cds_name

            else:
                gene_name = ""

            # Geneious often puts either a gene symbol OR a product
            # description in Name=... CDS.
            #
            # When it is a descriptive name, preserve it as Product.
            if not product:
                if cds_name and cds_name != gene_name:
                    product = cds_name

                elif cds_name:
                    # No real protein product was supplied. Keep something
                    # informative instead of a completely empty field.
                    product = cds_name

            if not product:
                product = "hypothetical protein"

            # Generate a unique synonym if the GFF lacks locus_tag.
            #
            # Do NOT use gene_name directly because this GFF contains
            # repeated gene symbols at different genomic locations.
            synonym = locus_tag or f"{seqid}_{index:06d}"

            records.append(
                {
                    "start": cds["start"],
                    "end": cds["end"],
                    "strand": cds["strand"],
                    "pid": protein_id,
                    "gene": gene_name,
                    "synonym": synonym,
                    "code": "-",
                    "cog": "-",
                    "product": product,
                    "line_number": cds["line_number"],
                }
            )

        result[seqid] = records

    return result


def safe_filename(value):
    """
    Make a seqid safe for use as a filename without changing the value
    written inside the PTT.
    """
    value = value.replace("/", "_")
    value = value.replace("\\", "_")

    if not value:
        return "replicon"

    return value


def write_ptt(seqid, records, replicon_length, out_dir):
    """
    Write one legacy PTT file.

    Format:
        <description> - 1..<length>
        <number> proteins
        Location Strand Length PID Gene Synonym Code COG Product
    """
    filename = safe_filename(seqid) + ".ptt"
    output_path = Path(out_dir) / filename

    with open(output_path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{seqid} - 1..{replicon_length}\n")
        handle.write(f"{len(records)} proteins\n")
        handle.write(
            "Location\tStrand\tLength\tPID\tGene\t"
            "Synonym\tCode\tCOG\tProduct\n"
        )

        for index, record in enumerate(records, 1):
            start = record["start"]
            end = record["end"]

            nucleotide_length = end - start + 1

            pid = record["pid"] or str(index)
            gene = record["gene"] or "-"
            synonym = record["synonym"] or f"{seqid}_{index:06d}"
            product = record["product"] or "-"

            handle.write(
                f"{start}..{end}\t"
                f"{record['strand']}\t"
                f"{nucleotide_length}\t"
                f"{pid}\t"
                f"{gene}\t"
                f"{synonym}\t"
                f"{record['code']}\t"
                f"{record['cog']}\t"
                f"{product}\n"
            )

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a Geneious GFF3 annotation and matching FASTA "
            "to legacy PTT files."
        )
    )

    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        help="Input GFF3 annotation file (.gff, .gff3, optionally .gz)",
    )

    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="Matching genome FASTA file (.fa/.fasta/.fna, optionally .gz)",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Directory in which PTT files will be written",
    )

    args = parser.parse_args()

    output_directory = Path(args.outdir)
    output_directory.mkdir(parents=True, exist_ok=True)

    print("Reading FASTA...", file=sys.stderr)

    try:
        fasta_lengths = parse_fasta_lengths(args.fasta)
    except Exception as exc:
        sys.exit(f"ERROR reading FASTA: {exc}")

    if not fasta_lengths:
        sys.exit("ERROR: no FASTA sequences were found.")

    print(
        f"Found {len(fasta_lengths)} FASTA sequence(s).",
        file=sys.stderr,
    )

    print("Reading GFF...", file=sys.stderr)

    try:
        (
            genes_by_coordinate,
            cds_by_seqid,
            gff_sequence_lengths,
            feature_counts,
            attribute_counts,
        ) = read_gff(args.gff)

    except Exception as exc:
        sys.exit(f"ERROR reading GFF: {exc}")

    total_cds = sum(len(records) for records in cds_by_seqid.values())

    if total_cds == 0:
        feature_summary = ", ".join(
            f"{name}={count}"
            for name, count in feature_counts.most_common()
        )

        sys.exit(
            "ERROR: no CDS features were found in the GFF.\n"
            f"Feature types seen: {feature_summary or '(none)'}"
        )

    print(
        f"Found {total_cds} CDS feature(s) on "
        f"{len(cds_by_seqid)} replicon(s).",
        file=sys.stderr,
    )

    if attribute_counts:
        print(
            "GFF attribute keys: "
            + ", ".join(
                f"{key}={count}"
                for key, count in attribute_counts.most_common()
            ),
            file=sys.stderr,
        )

    records_by_seqid = build_ptt_records(
        genes_by_coordinate,
        cds_by_seqid,
    )

    written = []
    skipped = []

    for seqid in sorted(
        records_by_seqid,
        key=lambda name: (
            records_by_seqid[name][0]["start"]
            if records_by_seqid[name]
            else 0
        ),
    ):
        records = records_by_seqid[seqid]

        fasta_length = fasta_lengths.get(seqid)
        gff_length = gff_sequence_lengths.get(seqid)

        if fasta_length is not None:
            replicon_length = fasta_length

            if (
                gff_length is not None
                and fasta_length != gff_length
            ):
                print(
                    f"WARNING: {seqid}: FASTA length is {fasta_length:,}, "
                    f"but GFF ##sequence-region length is "
                    f"{gff_length:,}. Using FASTA length.",
                    file=sys.stderr,
                )

        elif gff_length is not None:
            replicon_length = gff_length

            print(
                f"WARNING: '{seqid}' was not found as an exact FASTA "
                f"header token. Using GFF ##sequence-region length "
                f"{gff_length:,} instead.",
                file=sys.stderr,
            )

        else:
            skipped.append(seqid)

            fasta_examples = list(fasta_lengths.keys())[:10]

            print(
                f"WARNING: skipping '{seqid}': no matching FASTA "
                f"header and no ##sequence-region length.\n"
                f"FASTA header examples: {fasta_examples}",
                file=sys.stderr,
            )

            continue

        output_path = write_ptt(
            seqid,
            records,
            replicon_length,
            output_directory,
        )

        written.append(output_path)

        named_count = sum(
            1 for record in records if record["gene"]
        )

        unnamed_count = len(records) - named_count

        print(
            f"Wrote {output_path}: "
            f"{len(records)} proteins "
            f"({named_count} named, {unnamed_count} unnamed)"
        )

    if not written:
        sys.exit(
            "ERROR: no PTT files were written. "
            "Check GFF/FASTA sequence identifiers."
        )

    total_written_records = sum(
        len(records_by_seqid[seqid])
        for seqid in records_by_seqid
        if (output_directory / f"{safe_filename(seqid)}.ptt").exists()
    )

    print()
    print(
        f"Done: wrote {len(written)} PTT file(s) containing "
        f"{total_written_records} CDS records."
    )

    if skipped:
        print(
            "Skipped replicons: " + ", ".join(skipped),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()

