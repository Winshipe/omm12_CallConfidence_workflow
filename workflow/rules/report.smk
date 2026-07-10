"""
rules/report.smk
================
Generate a per-reference PDF QA/QC report by rendering an R Markdown
template with data produced by the rest of the pipeline.

What the report contains
------------------------
  • A stacked figure (complete_plot) with:
      - Detection rate (% mutations called above quality threshold) per
        genomic window, with one line per scenario.
      - A genomic annotation track showing transposable elements, IS
        elements, transposons, and regions with close cross-reference
        homologs – the features most likely to confound variant calling.
  • Three summary tables:
      - Per-scenario detection statistics.
      - Detection broken down by mutation type (transition / transversion).
      - Annotated region counts and total/mean span.

Inputs consumed
---------------
  assessment TSVs   — results/assessment/{scenario}_assessment.tsv
                       (one per scenario that includes this ref_id)
  challenging TSV   — results/annotation/{ref_id}/{ref_id}.challenging.tsv
                       (cross-reference homolog BED from find_homologs)
  mmseqs hit files  — results/annotation/{ref_id}/{db}_hits
                       (one per configured database)

Outputs
-------
  results/report/{ref_id}_report.pdf

Configuration keys used
-----------------------
  config["scenarios"]               — to find which scenarios include each ref
  config["annotation"]["mmseqs_databases"] — database names and paths
  config["report"]["genome_length"] — (optional) x-axis limit in nt;
                                       defaults to 0 (auto)
  config["report"]["window_size"]   — (optional) bin size for windowed
                                       summaries; defaults to 2000 nt
"""


# ---------------------------------------------------------------------------
# Helper: find all scenarios that include a given ref_id
# ---------------------------------------------------------------------------

def scenarios_for_ref(ref_id):
    """
    Return a list of scenario names (strings) whose contributor list
    includes ref_id with a non-zero mutated_fraction.
    """
    return [
        scenario
        for scenario, contribs in config["scenarios"].items()
        for c in contribs
        if c["ref_id"] == ref_id and c["mutated_fraction"] > 0
    ]


# ---------------------------------------------------------------------------
# Input functions
# ---------------------------------------------------------------------------

def report_assessment_tsvs(wildcards):
    """
    Return the list of assessment TSV paths for every scenario that
    includes this ref_id.
    """
    return [
        f"results/assessment/{scenario}_all_replicates.tsv"
        for scenario in scenarios_for_ref(wildcards.ref_id)
    ]

def report_mapped_reads(wildcards):
    """
    Return the list of assessment TSV paths for every scenario that
    includes this ref_id.
    """
    return [
        f"results/variant_calling/{scenario}/rep1/aligned.sam"
        for scenario in scenarios_for_ref(wildcards.ref_id)
    ]

def report_annotation_hits(wildcards):
    """
    Return the list of mmseqs hit file paths for this ref_id — one file
    per configured database.
    """
    # config["annotation"]["mmseqs_databases"] may be a dict (name→path)
    # or a list of paths depending on config style.
    dbs = config["annotation"]["mmseqs_databases"]
    if isinstance(dbs, dict):
        db_names = list(dbs.keys())
    else:
        # List of paths — derive a name from the basename
        db_names = [os.path.basename(p) for p in dbs]
    return [
        f"results/annotation/{wildcards.ref_id}/{db}_hits"
        for db in db_names
    ]


# ---------------------------------------------------------------------------
# Main report rule
# ---------------------------------------------------------------------------

rule generate_report:
    """
    Render the R Markdown report template for one reference genome and
    produce a self-contained PDF.

    The R Markdown template receives all data paths as named parameters
    so it never relies on hard-coded file locations.  Parameters are
    passed as comma-separated strings because Snakemake wildcards/params
    do not support lists in all contexts.
    """
    input:
        # R Markdown template (shared across all references)
        rmd="workflow/scripts/genome_report.Rmd",
        # Assessment result files (one per relevant scenario)
        assessment=report_assessment_tsvs,
        # Challenging / homologous regions BED-style TSV
        challenging="results/annotation/{ref_id}/{ref_id}.challenging.tsv",
        # mmseqs hit files (one per database)
        db_hits=report_annotation_hits,
        gff="results/annotation/{ref_id}/{ref_id}.gff",
        mut_reads="results/simulated/{ref_id}/rep1/mutated/{ref_id}_.sam",
        unmut_reads="results/simulated/{ref_id}/rep1/unmutated/{ref_id}_.sam",
        mapped_reads=report_mapped_reads
    output:
        html="results/report/{ref_id}_report.html",
    params:
        # Genome length for the x-axis.  Read from config if present,
        # otherwise default to 0 (R Markdown will compute it from the data).
        genome_length=lambda wc: config.get("report", {}).get("genome_length", 0),
        # Genomic window size (nt) used for the windowed detection summaries.
        window_size=lambda wc: config.get("report", {}).get("window_size", 2000),
        # Build the comma-separated strings that R Markdown expects.
        assessment_tsvs=lambda wc, input: ",".join(input.assessment),
        db_hit_files=lambda wc, input: ",".join(input.db_hits),
        db_names=lambda wc: (
            ",".join(config["annotation"]["mmseqs_databases"].keys())
            if isinstance(config["annotation"]["mmseqs_databases"], dict)
            else ",".join(
                os.path.basename(p)
                for p in config["annotation"]["mmseqs_databases"]
            )
        ),
    log:
        "logs/report/{ref_id}_report.log",
    conda:
        "../envs/r_report.yaml"
    shell:
        # Render the R Markdown document from the command line using Rscript.
        #
        # We pass every input path as a named parameter so the template is
        # fully data-agnostic.  The output PDF is written directly to its
        # final location.
        #
        # Flags explained:
        #   --vanilla        : skip .Rprofile / .Renviron to ensure reproducibility
        #   rmarkdown::render: the standard R Markdown rendering function
        #   output_file      : absolute path prevents working-directory issues
        #   params           : named list of parameters injected into the template
        """
        mkdir -p results/report
        Rscript --vanilla -e "
          rmarkdown::render(
            input       = '{input.rmd}',
            output_file = file.path('"""+os.getcwd()+"""', '{output.html}'),
            params      = list(
              scenario        = '{wildcards.scenario}',
              workingdir      = '""" + os.getcwd() + """',
              ref_id          = '{wildcards.ref_id}',
              assessment_tsvs = '{params.assessment_tsvs}',
              challenging_tsv = '{input.challenging}',
              db_hit_files    = '{params.db_hit_files}',
              db_names        = '{params.db_names}',
              genome_length   = {params.genome_length},
              window_size     = {params.window_size},
              unmut_sam       = '{input.unmut_reads}',
              mut_sam         = '{input.mut_reads}',
              mapped_sam      = '{input.mapped_reads}',
              gff             = '{input.gff}' 
            )
          )
        " &> {log}
        """
