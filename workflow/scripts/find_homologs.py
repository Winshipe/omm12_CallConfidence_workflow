#from __future__ import annotations # doesn't work in a snakemake context b/c script is placed in some wrapper
"""
scripts/find_homologs.py
========================
Snakemake helper script – called by the `find_homologs` rule in annotate.smk.

Purpose
-------
1. Parse the MMseqs2 easy-cluster output TSV (two columns: representative,
   member) to collect every gene ID that appears in a multi-member cluster,
   i.e. genes that have a homolog in at least one other reference genome.

2. For each such homologous gene, look up its genomic coordinates from the
   Prodigal FASTA headers embedded in the nucleotide gene FASTAs (.fna).

3. Write a 6-column BED-style TSV:
       chrom  start  stop  name  score  strand

Snakemake I/O contract (set in the `find_homologs` rule)
---------------------------------------------------------
  snakemake.input[0]   : str        – path to clusters.txt (two-column TSV)
  snakemake.input[1:]  : list[str]  – paths to per-reference .fna FASTAs
  snakemake.output[0]  : str        – path for the output BED-style TSV
"""


from Bio import SeqIO


# ---------------------------------------------------------------------------
# Step 1 – identify homologous gene IDs from the cluster file
# ---------------------------------------------------------------------------

def find_homolog_ids(cluster_path: str) :#-> set[str]:
    """
    Read a two-column MMseqs2 cluster TSV and return the set of all gene IDs
    that belong to a multi-member cluster (i.e. genes with at least one
    homolog).

    The cluster file format is:
        <representative_id>  <member_id>
    A gene is its own representative when representative == member; any row
    where they differ indicates a genuine cross-sequence homology.

    Parameters
    ----------
    cluster_path:
        Path to the *_cluster.tsv produced by ``mmseqs easy-cluster``.

    Returns
    -------
    set[str]
        All gene IDs (both representative and member) that participate in a
        multi-member cluster.
    """
    homologs: set[str] = set()

    with open(cluster_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cluster_rep, member = line.split("\t")
            # A singleton cluster has representative == member; skip those.
            if cluster_rep != member:
                homologs.add(cluster_rep)
                homologs.add(member)

    return homologs


# ---------------------------------------------------------------------------
# Step 2 – extract gene coordinates from Prodigal FASTA headers
# ---------------------------------------------------------------------------

def get_gene_coords(fa_path: str):# -> dict[str, tuple[str, str, str]]:
    """
    Parse a Prodigal-generated nucleotide FASTA (.fna) and return a mapping of
    gene ID → (start, stop, strand).

    Expected Prodigal header format
    --------------------------------
    >contig_23611_1 # 2 # 271 # 1 # ID=1_1;partial=10;start_type=Edge;...
     ^--- name ---^   ^s^ ^e^  ^strand^

    Fields are space-delimited with ' # ' as the separator:
        [0] gene name
        [1] start coordinate (1-based)
        [2] stop  coordinate (1-based, inclusive)
        [3] strand  (1 = forward, -1 = reverse)
        [4] misc attributes (key=value pairs)

    Parameters
    ----------
    fa_path:
        Path to a Prodigal .fna output file.

    Returns
    -------
    dict mapping gene_id → (start, stop, strand_str)
    """
    gene_coords: dict[str, tuple[str, str, str]] = {}

    for record in SeqIO.parse(fa_path, "fasta"):
        # SeqIO sets record.description to the full header line (without '>').
        # Split on ' # ' to get the five Prodigal fields.
        parts = record.description.split(" # ")
        if len(parts) < 4:
            # Skip malformed headers rather than crashing the whole run.
            print("Malformed header in prodigal output:")
            print(record.description)
            continue
        name = parts[0]
        start = parts[1]
        stop = parts[2]
        strand = parts[3]
        gene_coords[name] = (start, stop, strand)

    return gene_coords


# ---------------------------------------------------------------------------
# Step 3 – collect coords from all input FASTAs
# ---------------------------------------------------------------------------

def read_fasta_files(fa_paths: list[str]):# -> dict[str, dict[str, tuple[str, str, str]]]:
    """
    Build a nested dict of { fasta_path: { gene_id: (start, stop, strand) } }
    by calling get_gene_coords on every supplied .fna file.

    Parameters
    ----------
    fa_paths:
        List of Prodigal .fna file paths (one per reference genome).

    Returns
    -------
    dict mapping each fasta path to its gene_coords dict.
    """
    fasta_entries: dict[str, dict[str, tuple[str, str, str]]] = {}

    for fa_path in fa_paths:
        fasta_entries[fa_path] = get_gene_coords(fa_path)

    return fasta_entries


# ---------------------------------------------------------------------------
# Step 4 – write BED-style output
# ---------------------------------------------------------------------------

def write_bed(
    homologs: set[str],
    fasta_entries: dict[str, dict[str, tuple[str, str, str]]],
    out_path: str,
):# -> None:
    """
    For every homologous gene ID, look up its coordinates and write one line
    to the output file in BED-like format:

        chrom  start  stop  name  score  strand

    The contig name (chrom) is reconstructed by stripping the trailing
    Prodigal gene index from the gene ID (e.g. 'contig_23611_1' → 'contig_23611').
    Score is set to '.' (unknown) as MMseqs2 cluster output carries no score.

    Parameters
    ----------
    homologs:
        Set of gene IDs identified as homologs.
    fasta_entries:
        Nested dict of { fasta_path: { gene_id: (start, stop, strand) } }.
    out_path:
        Destination file path for the BED-style TSV.
    """
    template = "{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}\n"

    with open(out_path, "w") as out_fh:
        # Write a header so the file is self-documenting.
        out_fh.write("#chrom\tstart\tstop\tname\tscore\tstrand\n")

        for homolog in sorted(homologs):  # sorted for deterministic output
            for fa_path, coord_dict in fasta_entries.items():
                if homolog not in coord_dict:
                    continue
                start, stop, strand = coord_dict[homolog]
                # Prodigal names genes as <contig>_<index>; drop the index to
                # recover the contig / chromosome name.
                chrom = "_".join(homolog.split("_")[:-1])
                # Convert Prodigal strand integer to BED convention (+/-)
                strand_symbol = "+" if strand == "1" else "-"
                out_fh.write(
                    template.format(
                        chrom=chrom,
                        start=start,
                        stop=stop,
                        name=homolog,
                        score=".",
                        strand=strand_symbol,
                    )
                )


# ---------------------------------------------------------------------------
# Entry point (called by Snakemake's `script:` directive)
# ---------------------------------------------------------------------------

def main():# -> None:
    """
    Orchestrate the three steps using paths supplied by Snakemake:
      - snakemake.input[0]   : cluster TSV
      - snakemake.input[1:]  : list of .fna FASTA files
      - snakemake.output[0]  : output BED-style TSV
    """
    cluster_path: str = snakemake.input[0]          
    fa_paths: list[str] = list(snakemake.input[1:])  
    out_path: str = snakemake.output[0]             

    homologs = find_homolog_ids(cluster_path)
    fasta_entries = read_fasta_files(fa_paths)
    write_bed(homologs, fasta_entries, out_path)


main()
