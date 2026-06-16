"""
Publication-oriented Snakemake workflow for the Pantoea/Kalamiella LPS analysis.

Run from the directory containing this Snakefile and the helper scripts, for example:

    snakemake --cores 8 \
      --config genome_dir=data/genomes \
               bakta_db=resources/bakta_db \
               eggnog_data_dir=resources/eggnog \
               go_obo=resources/go-basic.obo \
               lps_terms=resources/lps_go_terms.tsv \
               group_map=metadata.tsv

The workflow intentionally uses relative, configurable paths. No machine-specific
paths are required. If a sample table is supplied, it must be tab-delimited with
columns: sample, fasta. Otherwise all *.fa, *.fna, and *.fasta files in genome_dir
are used and sample names are inferred from file names.
"""

import csv
from pathlib import Path


config.setdefault("genome_dir", "data/genomes")
config.setdefault("sample_table", "")
config.setdefault("results_dir", "results/genome_annotations")
config.setdefault("analysis_dir", "results/lps_analysis")
config.setdefault("plots_dir", "results/lps_analysis/plots")
config.setdefault("shared_genes_dir", "results/lps_analysis/lps_shared_genes")
config.setdefault("selection_dir", "results/lps_analysis/selection_curation")
config.setdefault("phase2_dir", "results/lps_analysis/selection_phase2")

config.setdefault("bakta", "bakta")
config.setdefault("bakta_db", "resources/bakta_db")
config.setdefault("emapper", "emapper.py")
config.setdefault("eggnog_data_dir", "resources/eggnog")
config.setdefault("diamond", "diamond")
config.setdefault("threads_per_job", 2)
config.setdefault("eggnog_method", "diamond")
config.setdefault("eggnog_sensmode", "sensitive")

config.setdefault("go_obo", "resources/go-basic.obo")
config.setdefault("lps_terms", "resources/lps_go_terms.tsv")
config.setdefault("goslim_map", "resources/goslim_prok_mapping.tsv")
config.setdefault("group_map", "metadata.tsv")
config.setdefault("metadata", "metadata.tsv")
config.setdefault("reference_isolate", "")


SCRIPT_DIR = Path(workflow.basedir) / "scripts"
GENOME_DIR = Path(config["genome_dir"])
RESULTS_DIR = Path(config["results_dir"])
ANALYSIS_DIR = Path(config["analysis_dir"])
PLOTS_DIR = Path(config["plots_dir"])
SHARED_GENES_DIR = Path(config["shared_genes_dir"])
SELECTION_DIR = Path(config["selection_dir"])
PHASE2_DIR = Path(config["phase2_dir"])
THREADS_PER_JOB = int(config["threads_per_job"])


def load_samples():
    sample_table = str(config.get("sample_table", "")).strip()
    if sample_table:
        with open(sample_table, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"sample", "fasta"}.issubset(reader.fieldnames or []):
                raise ValueError("sample_table must contain tab-delimited columns: sample, fasta")
            rows = [(row["sample"], row["fasta"]) for row in reader if row.get("sample") and row.get("fasta")]
    else:
        exts = ["*.fa", "*.fna", "*.fasta"]
        paths = []
        for ext in exts:
            paths.extend(sorted(GENOME_DIR.glob(ext)))
        rows = [(p.name.rsplit(".", 1)[0], str(p)) for p in sorted(paths)]

    if not rows:
        raise ValueError(
            "No genomes were found. Provide --config sample_table=... or place FASTA files in genome_dir."
        )
    return dict(rows)


FASTA_BY_SAMPLE = load_samples()
SAMPLES = sorted(FASTA_BY_SAMPLE)
TARGET_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]


rule all:
    input:
        ANALYSIS_DIR / "lps_eggnog_hits.tsv",
        ANALYSIS_DIR / "lps_strict_hits.tsv",
        ANALYSIS_DIR / "lps_eggnog_go_terms.tsv",
        ANALYSIS_DIR / "rfb_eggnog_locations.tsv",
        ANALYSIS_DIR / "rfb_operon_summary.tsv",
        ANALYSIS_DIR / "lps_strict_prevalence_diff.tsv",
        ANALYSIS_DIR / "lps_go_category_summary.tsv",
        ANALYSIS_DIR / "lps_ko_category_summary.tsv",
        PLOTS_DIR / "lps_terms_norm.png",
        PLOTS_DIR / "lps_goslim_norm.png",
        PLOTS_DIR / "oantigen_gene_fractions.png",
        expand(str(SELECTION_DIR / "{gene}.ffn"), gene=TARGET_GENES),
        PHASE2_DIR / "diamond_curation_summary.tsv",


rule run_bakta:
    input:
        genome=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
    output:
        faa=RESULTS_DIR / "{sample}" / "bakta_output" / "{sample}.faa",
        ffn=RESULTS_DIR / "{sample}" / "bakta_output" / "{sample}.ffn",
        gff3=RESULTS_DIR / "{sample}" / "bakta_output" / "{sample}.gff3",
        proteins=RESULTS_DIR / "{sample}" / "proteins.faa",
    params:
        outdir=lambda wildcards: RESULTS_DIR / wildcards.sample / "bakta_output",
        prefix=lambda wildcards: wildcards.sample,
        bakta=config["bakta"],
        db=config["bakta_db"],
    threads: THREADS_PER_JOB
    log:
        RESULTS_DIR / "{sample}" / "logs" / "bakta.log",
    shell:
        r"""
        set -euo pipefail
        rm -rf {params.outdir}
        mkdir -p {params.outdir} $(dirname {log})

        {params.bakta} \
          --db {params.db} \
          --output {params.outdir} \
          --prefix {params.prefix} \
          --threads {threads} \
          --force \
          {input.genome} \
          > {log} 2>&1

        cp {output.faa} {output.proteins}
        """


rule eggnog_mapper:
    input:
        faa=RESULTS_DIR / "{sample}" / "proteins.faa",
    output:
        hits=RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.hits",
        seed=RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.seed_orthologs",
        ann=RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations",
    params:
        outdir=lambda wildcards: RESULTS_DIR / wildcards.sample / "eggnog_out",
        prefix=lambda wildcards: wildcards.sample,
        emapper=config["emapper"],
        data_dir=config["eggnog_data_dir"],
        method=config["eggnog_method"],
        sensmode=config["eggnog_sensmode"],
    threads: THREADS_PER_JOB
    log:
        RESULTS_DIR / "{sample}" / "logs" / "eggnog_mapper.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})

        {params.emapper} \
          -i {input.faa} \
          -o {params.prefix} \
          --itype proteins \
          --cpu {threads} \
          --output_dir {params.outdir} \
          --data_dir {params.data_dir} \
          -m {params.method} \
          --sensmode {params.sensmode} \
          --dmnd_iterate no \
          --pfam_realign none \
          > {log} 2>&1
        """


rule extract_lps_hits:
    input:
        expand(str(RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations"), sample=SAMPLES),
    output:
        proteins=ANALYSIS_DIR / "lps_eggnog_proteins.faa",
        nucleotides=ANALYSIS_DIR / "lps_eggnog_nucleotides.ffn",
        hits=ANALYSIS_DIR / "lps_eggnog_hits.tsv",
    shell:
        r"""
        python {SCRIPT_DIR}/extract_lps_from_eggnog.py \
          --root {RESULTS_DIR} \
          --out-proteins {output.proteins} \
          --out-nucleotides {output.nucleotides} \
          --out-hits {output.hits}
        """


rule extract_strict_lps_hits:
    input:
        expand(str(RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations"), sample=SAMPLES),
    output:
        proteins=ANALYSIS_DIR / "lps_strict_proteins.faa",
        nucleotides=ANALYSIS_DIR / "lps_strict_nucleotides.ffn",
        hits=ANALYSIS_DIR / "lps_strict_hits.tsv",
    shell:
        r"""
        python {SCRIPT_DIR}/extract_lps_from_eggnog.py \
          --root {RESULTS_DIR} \
          --strict \
          --out-proteins {output.proteins} \
          --out-nucleotides {output.nucleotides} \
          --out-hits {output.hits}
        """


rule extract_lps_go_terms:
    input:
        expand(str(RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations"), sample=SAMPLES),
        go_obo=config["go_obo"],
    output:
        ANALYSIS_DIR / "lps_eggnog_go_terms.tsv",
    shell:
        r"""
        python {SCRIPT_DIR}/extract_lps_go_terms.py \
          --root {RESULTS_DIR} \
          --obo {input.go_obo} \
          --strict \
          --out {output}
        """


rule map_rfb_from_eggnog:
    input:
        expand(str(RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations"), sample=SAMPLES),
    output:
        hits=ANALYSIS_DIR / "rfb_eggnog_locations.tsv",
        operons=ANALYSIS_DIR / "rfb_operon_sequences.fna",
        proteins=ANALYSIS_DIR / "rfb_proteins.faa",
    shell:
        r"""
        python {SCRIPT_DIR}/find_rfb_from_eggnog.py \
          --root {RESULTS_DIR} \
          --out-hits {output.hits} \
          --out-operons {output.operons} \
          --out-proteins {output.proteins}
        """


rule find_rfb_operon:
    input:
        expand(str(RESULTS_DIR / "{sample}" / "bakta_output" / "{sample}.gff3"), sample=SAMPLES),
    output:
        summary=ANALYSIS_DIR / "rfb_operon_summary.tsv",
        hits=ANALYSIS_DIR / "rfb_operon_hits.tsv",
    shell:
        r"""
        python {SCRIPT_DIR}/find_rfb_operon.py \
          --root {RESULTS_DIR} \
          --out-summary {output.summary} \
          --out-hits {output.hits}
        """


rule extract_shared_lps_genes:
    input:
        ANALYSIS_DIR / "lps_strict_hits.tsv",
    output:
        directory(SHARED_GENES_DIR),
    shell:
        r"""
        python {SCRIPT_DIR}/extract_shared_lps_genes.py \
          --root {RESULTS_DIR} \
          --outdir {output}
        """


rule build_strict_lps_presence:
    input:
        hits=ANALYSIS_DIR / "lps_strict_hits.tsv",
        group_map=config["group_map"],
    output:
        wide=ANALYSIS_DIR / "lps_strict_presence.tsv",
        long=ANALYSIS_DIR / "lps_strict_presence_long.tsv",
        prevalence=ANALYSIS_DIR / "lps_strict_prevalence_diff.tsv",
    run:
        import csv
        from collections import defaultdict

        def norm(name):
            return name.replace("\ufeff", "").replace("\u00a0", " ").strip().lower()

        def pick(header, candidates):
            lookup = {norm(h): h for h in header}
            for candidate in candidates:
                if candidate in lookup:
                    return lookup[candidate]
            return None

        groups = {}
        with open(input.group_map, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            sample_col = pick(reader.fieldnames or [], ["sample", "sample id", "isolate", "genome"])
            group_col = pick(reader.fieldnames or [], ["group", "source", "category"])
            if not sample_col or not group_col:
                raise ValueError("group_map must contain sample/isolate and group/source columns")
            for row in reader:
                sample = row.get(sample_col, "").strip()
                group = row.get(group_col, "").strip()
                if sample and group:
                    groups[sample] = group

        genes = set()
        present = defaultdict(set)
        with open(input.hits, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            isolate_col = pick(reader.fieldnames or [], ["isolate", "sample", "genome"])
            gene_col = pick(reader.fieldnames or [], ["gene", "gene_id", "preferred_name", "preferred"])
            if not isolate_col or not gene_col:
                raise ValueError("lps_strict_hits must contain isolate and gene/preferred_name columns")
            for row in reader:
                isolate = row.get(isolate_col, "").strip()
                gene = row.get(gene_col, "").strip()
                if isolate and gene:
                    genes.add(gene)
                    present[isolate].add(gene)

        samples = sorted(groups)
        genes = sorted(genes)
        Path(output.wide).parent.mkdir(parents=True, exist_ok=True)

        with open(output.wide, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["isolate", "group"] + genes)
            for sample in samples:
                writer.writerow([sample, groups[sample]] + [int(gene in present[sample]) for gene in genes])

        with open(output.long, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["isolate", "group", "gene_id", "present"])
            for sample in samples:
                for gene in genes:
                    writer.writerow([sample, groups[sample], gene, int(gene in present[sample])])

        by_group = defaultdict(list)
        for sample, group in groups.items():
            by_group[group].append(sample)

        with open(output.prevalence, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene", "group", "present", "total", "prevalence"])
            for gene in genes:
                for group in sorted(by_group):
                    total = len(by_group[group])
                    count = sum(1 for sample in by_group[group] if gene in present[sample])
                    prevalence = count / total if total else 0
                    writer.writerow([gene, group, count, total, f"{prevalence:.6f}"])


rule plot_lps_terms:
    input:
        hits=ANALYSIS_DIR / "lps_eggnog_hits.tsv",
        terms=config["lps_terms"],
        group_map=config["group_map"],
    output:
        tsv=ANALYSIS_DIR / "lps_terms_norm.tsv",
        plot=PLOTS_DIR / "lps_terms_norm.png",
    shell:
        r"""
        python {SCRIPT_DIR}/plot_lps_terms.py \
          --lps-hits {input.hits} \
          --annotations-root {RESULTS_DIR} \
          --lps-terms {input.terms} \
          --group-map {input.group_map} \
          --out-tsv {output.tsv} \
          --out-plot {output.plot}
        """


rule plot_lps_goslim:
    input:
        hits=ANALYSIS_DIR / "lps_eggnog_hits.tsv",
        goslim=config["goslim_map"],
        group_map=config["group_map"],
    output:
        tsv=ANALYSIS_DIR / "lps_goslim_norm.tsv",
        plot=PLOTS_DIR / "lps_goslim_norm.png",
    shell:
        r"""
        python {SCRIPT_DIR}/plot_lps_goslim.py \
          --lps-hits {input.hits} \
          --annotations-root {RESULTS_DIR} \
          --goslim-map {input.goslim} \
          --group-map {input.group_map} \
          --out-tsv {output.tsv} \
          --out-plot {output.plot}
        """


rule plot_oantigen_genes:
    input:
        expand(str(RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations"), sample=SAMPLES),
        group_map=config["group_map"],
    output:
        tsv=ANALYSIS_DIR / "oantigen_gene_fractions.tsv",
        plot=PLOTS_DIR / "oantigen_gene_fractions.png",
    shell:
        r"""
        python {SCRIPT_DIR}/plot_oantigen_genes.py \
          --annotations-root {RESULTS_DIR} \
          --group-map {input.group_map} \
          --out-tsv {output.tsv} \
          --out-plot {output.plot}
        """


rule build_lps_go_category_table:
    input:
        terms=config["lps_terms"],
        group_map=config["group_map"],
        annotations=expand(str(RESULTS_DIR / "{sample}" / "eggnog_out" / "{sample}.emapper.annotations"), sample=SAMPLES),
    output:
        mapping=ANALYSIS_DIR / "lps_go_term_category_map.tsv",
        table=ANALYSIS_DIR / "lps_go_category_counts.tsv",
        summary=ANALYSIS_DIR / "lps_go_category_summary.tsv",
    shell:
        r"""
        python {SCRIPT_DIR}/build_lps_go_category_table.py \
          --annotations-root {RESULTS_DIR} \
          --lps-terms {input.terms} \
          --group-map {input.group_map} \
          --mapping-out {output.mapping} \
          --out-table {output.table} \
          --out-summary {output.summary}
        """


rule build_lps_ko_category_table:
    input:
        gene_list=ANALYSIS_DIR / "lps_strict_prevalence_diff.tsv",
        presence_long=ANALYSIS_DIR / "lps_strict_presence_long.tsv",
    output:
        mapping=ANALYSIS_DIR / "lps_ko_category_map.tsv",
        table=ANALYSIS_DIR / "lps_ko_category_counts.tsv",
        summary=ANALYSIS_DIR / "lps_ko_category_summary.tsv",
    shell:
        r"""
        python {SCRIPT_DIR}/build_lps_ko_category_table.py \
          --gene-list {input.gene_list} \
          --presence-long {input.presence_long} \
          --out-mapping {output.mapping} \
          --out-table {output.table} \
          --out-summary {output.summary}
        """


rule curate_selection_orthologs:
    input:
        nucleotides=ANALYSIS_DIR / "lps_eggnog_nucleotides.ffn",
        proteins=ANALYSIS_DIR / "lps_eggnog_proteins.faa",
        rfb_locations=ANALYSIS_DIR / "rfb_eggnog_locations.tsv",
        lps_hits=ANALYSIS_DIR / "lps_strict_hits.tsv",
        groups=config["group_map"],
    output:
        expand(str(SELECTION_DIR / "{gene}.ffn"), gene=TARGET_GENES),
    shell:
        r"""
        python {SCRIPT_DIR}/curate_selection_orthologs.py \
          --root {RESULTS_DIR} \
          --groups {input.groups} \
          --nucleotides {input.nucleotides} \
          --proteins {input.proteins} \
          --rfb-locations {input.rfb_locations} \
          --lps-hits {input.lps_hits} \
          --outdir {SELECTION_DIR}
        """


rule prepare_selection_phase2:
    input:
        expand(str(SELECTION_DIR / "{gene}.ffn"), gene=TARGET_GENES),
    output:
        PHASE2_DIR / "diamond_curation_summary.tsv",
    params:
        reference=lambda wildcards: (
            f"--reference-isolate {config['reference_isolate']}"
            if str(config.get("reference_isolate", "")).strip()
            else ""
        ),
        diamond=config["diamond"],
    shell:
        r"""
        python {SCRIPT_DIR}/prepare_selection_phase2.py \
          --indir {SELECTION_DIR} \
          --outdir {PHASE2_DIR} \
          --diamond {params.diamond} \
          {params.reference}
        """
