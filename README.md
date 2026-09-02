# Genome annotation & KEGG module completeness pipeline

A Nextflow pipeline (`pipeline.nf`) that takes per-genome HMMER (KO) and InterProScan
annotations and produces:

- a combined per-accession annotation table (KO + InterPro member databases + taxonomy)
- a characterized / uncharacterized classification for every protein
- KEGG module completeness per genome (accounting for KO complexes and alternative KOs)
- presence/absence of KOs that aren't part of any KEGG module
- normalized counts of uncharacterized proteins per taxonomic rank
- a synteny (gene-neighborhood) analysis around uncharacterized-protein hits

A set of standalone `tool_*.py` / `tool_r_*.R` scripts (not called from `pipeline.nf`)
run a parallel arCOG-based annotation/classification workflow and turn its output (and
the KEGG module-completeness output) into plots.

## Layout

```
.
├── pipeline.nf          # Nextflow workflow, 13 processes, wraps scripts/pipeline_*.py
├── requirements.txt      # Python deps
├── data/                 # reference tables the pipeline is configured against (see below)
└── scripts/
    ├── pipeline_*.py     # called from pipeline.nf, in workflow order (see table)
    ├── select_kegg_modules.sh  # switch data/kegg_modules_architecture.txt between manual/original
    ├── tool_*.py         # standalone arCOG / plotting-prep tools
    └── r/
        └── tool_r_*.R    # plotting scripts consuming the tool_*.py / pipeline output
```

## Requirements

- Python 3.8+, [Biopython](https://biopython.org/) (`pip install -r requirements.txt`)
- R with `ggplot2`, `dplyr` (boxplot) and `pheatmap` (heatmap) for the plotting scripts
- [Nextflow](https://www.nextflow.io/) to run `pipeline.nf`
- Precomputed per-genome HMMER and InterProScan output files (see "Preparing your
  genome data directory" below)

Tested with Nextflow 24.04.4-17.0.6 and Python 3.11.9.

## Preparing your genome data directory

Before running the pipeline, generate per-genome annotations and lay everything out in
one directory — this becomes the directory you launch Nextflow *from* (see "Running the
Nextflow pipeline" below).

### 1. Run InterProScan per genome

```bash
interproscan.sh -cpu $n -f tsv -iprlookup -i $protein -o $out
```
`$n` = number of CPUs, `$protein` = the genome's protein FASTA, `$out` = where to save
the output.

### 2. Run HMMER against the KOfam profiles (not `kofam_scan`) per genome

Download the profiles and `ko_list` from
[the KofamKOALA site](https://www.genome.jp/tools/kofamkoala/), then, with HMMER 3.4:

```bash
hmmsearch --cpu=$n -o $errlog -T0 --tblout=$tbl $hmm $protein
```
`$n` = number of CPUs, `$protein` = the genome's protein FASTA, `$tbl` = where to save
the output, `$hmm` = the concatenated KOfam profiles file.

### 3. Prepare `taxonomy.tsv`

One line per genome, tab-separated, **no header**:

```
genome_accession    domain    phylum    class    order    family    genus    species    strain
```

### 4. Lay out the directory

Put everything below in one directory, using these exact filenames (`$genome` = the
genome accession as it appears in `taxonomy.tsv`, without the `$`):

| File | Contents |
|---|---|
| `$genome.faa` | genome protein FASTA |
| `$genome.gff` or `$genome_feature_table.txt` | feature table/GFF (used by the synteny step) |
| `$genome_interpro.tsv` | InterProScan output from step 1 |
| `$genome_hmmer.tsv` | HMMER output from step 2 |
| `taxonomy.tsv` | from step 3 |
| `$genome_arcogs.tsv` *(optional)* | only needed for the arCOG scripts — see "Running the standalone arCOG & plotting workflow" below |

## Reference data (`data/`)

| File | Used by | What it is |
|---|---|---|
| `kegg_modules_architecture.txt` | `mapModule`, `analyzeChar`, `prepCutoffs`, `calcModule` | KEGG module definitions (KOs, `+` for complexes, `,` for alternatives), blocks separated by a dashed line |
| `kegg_modules_architecture_original.txt` | — (reference only) | Unmodified copy of the same file, kept for comparison — **currently byte-identical** to `kegg_modules_architecture.txt` (see note below) |
| `manual_modules_architecture.txt` | — (not wired into `pipeline.nf`) | Manually curated extra modules (oxygen/sulfur/nitrogen reductases, archaeal riboflavin/heme/cobalamin variants, etc.) in the same format — see note below |
| `ko00001_level_C.tsv` | `mergeTables` | KO → name → KEGG pathway ("level C") mapping |
| `ko_cutoffs.tsv` | `filterHmmer`, `prepCutoffs` | Per-KO HMMER score/e-value cutoffs (`knum, threshold, score_type, ...`) |
| `all_ko_names.tsv` | `analyzeChar` | KO → short name lookup |
| `kegg_module.txt` | `analyzeChar` | KEGG module ID → module name, one per line |
| `nonmodule_kos_list.txt` | `calcNonmodule` | KOs with no KEGG module affiliation |
| `tigr_to_remove.txt`, `panther_to_remove.txt` | `classifyUnch` | Curated lists of TIGRFAM/PANTHER IDs treated as "uncharacterized" |

### Choosing a KEGG module set

`pipeline.nf` always reads `data/kegg_modules_architecture.txt` — that file is a
*working copy*, and which module set it holds is switched by copying one of the other
two files over it:

```bash
scripts/select_kegg_modules.sh manual     # use manual_modules_architecture.txt
                                           # (the manually curated modules described
                                           # in the publication - oxygen/sulfur/
                                           # nitrogen reductases, archaeal riboflavin/
                                           # heme/cobalamin variants, etc.)

scripts/select_kegg_modules.sh original   # use kegg_modules_architecture_original.txt
                                           # (unmodified KEGG modules only)
```

Run whichever you want before `nextflow run /absolute/path/to/this/repo/pipeline.nf`
(see "Running the Nextflow pipeline" below). `data/kegg_modules_architecture.txt` ships
set to `original`. Note that
`tool_prep_matrices_for_plotting.py`'s `marker_modules` list (used for the marker-module
heatmap) only resolves against the `manual` set — switch to `manual` before that step if
you're using it. That list also includes two entries, `"DsrAB"` and `"Qmo"`, that aren't
module IDs in either file — that lookup will `KeyError` if reached.

## Pipeline steps

Each `scripts/pipeline_*.py` can be run standalone (`python3 scripts/<name>.py --help`
for its arguments) or as part of the Nextflow workflow.

| # | Nextflow process | Script | What it does |
|---|---|---|---|
| 1 | `parseInterpro` | `pipeline_parse_interpro.py` | Parse raw InterProScan TSVs into one per-accession table |
| 2 | `filterHmmer` | `pipeline_filter_hmmer.py` | Filter HMMER hits by KO score/e-value cutoffs, keep best hit per accession |
| 3 | `createAccgen` | `pipeline_prepare_acc_gen.py` | Build accession → genome map from `.faa` headers |
| 4 | `prepCutoffs` | `pipeline_prep_cutoffs.py` | Build per-module table of KO score cutoffs |
| 5 | `mergeTables` | `pipeline_merge_table.py` | Merge InterPro + HMMER + taxonomy into one full table |
| 6 | `summaryAnno` | `pipeline_calc_general_presence_annotations.py` | Presence/absence matrix + counts across annotation DBs |
| 7 | `classifyUnch` | `pipeline_classify_un_characterized.py` | Classify each accession as characterized / uncharacterized |
| 8 | `countUnch` | `pipeline_normalize_uncharacterized_counts.py` | Normalize uncharacterized counts per genome, summarize per taxon |
| 9 | `mapModule` | `pipeline_map_kegg_module.py` | Map characterized hits to KEGG module(s) |
| 10 | `analyzeChar` | `pipeline_analyze_characterized.py` | Per-genome, per-module KO presence/absence ("fullness"), split per module |
| 11 | `calcModule` | `pipeline_calc_kegg_modules_completeness.py` | Module completeness per genome (complexes + alternatives), binned + taxonomy |
| 12 | `calcNonmodule` | `pipeline_calc_presence_nonmodule_kos.py` | Presence/absence of KOs with no module affiliation |
| 13 | `syntenyUnch` | `pipeline_synteny_uncharacterized.py` | Gene-neighborhood ("synteny") analysis around uncharacterized hits |

### Standalone tools (`scripts/tool_*.py`, `scripts/r/tool_r_*.R`)

| Script | What it does |
|---|---|
| `tool_arcog_filter.py` | Best-hit filtering of raw arCOG search output |
| `tool_arcogs_map_info.py` | Annotate filtered arCOG hits with model/function/definition + taxonomy |
| `tool_arcog_count_hits_per_genome.py` | Count arCOG hits per genome, % unclassified per phylum |
| `tool_boxplot_prep.py` | Reshape arCOG/pipeline "uncharacterized %" tables to long format for boxplots |
| `tool_prep_matrices_for_plotting.py` | Presence/partial/absence matrices per taxonomic order, for heatmaps |
| `r/tool_r_boxplot.R` | Boxplot of uncharacterized % per taxon from `tool_boxplot_prep.py` output |
| `r/tool_r_heatmap.R` | Heatmap (e.g. of marker-module presence) from `tool_prep_matrices_for_plotting.py` output |
| `r/tool_r_stacked_barplot.R` | Stacked barplot of Present/Partial/Absent per taxon from the same matrices |

All three R scripts take `<input> <output>` as command-line arguments now (see each
file's header comment for the exact usage), instead of a hardcoded `/path/to/...`.

## Running the Nextflow pipeline

Nextflow distinguishes between where the *pipeline* lives (`projectDir`, this repo) and
where you *run it from* (`launchDir`). `pipeline.nf` reads `taxonomy.tsv` and the
per-genome files from `launchDir`, so:

1. `cd` into the genome data directory you prepared above (the one containing
   `taxonomy.tsv` and the per-genome files)
2. run the pipeline from there, pointing at wherever you cloned/downloaded this repo:

```bash
nextflow run /absolute/path/to/this/repo/pipeline.nf
```

This creates a `work/` directory in your data directory (Nextflow's scratch space —
safe to delete once a run finishes) and publishes results to `./output/`, also inside
your data directory.

### Using custom KO cutoffs

`filterHmmer` and `prepCutoffs` read their KO score/e-value cutoffs from a
Nextflow parameter, `params.ko_cutoffs`, which defaults to the shipped
`data/ko_cutoffs.tsv`. To use your own cutoffs, point the pipeline at your file on the
command line instead of overwriting the shipped default:

```bash
nextflow run /absolute/path/to/this/repo/pipeline.nf --ko_cutoffs /path/to/your_cutoffs.tsv
```

Same tabular format as the shipped file (`KO\tcutoff\tscore_type`, no header). This
also makes it easy to keep several cutoff sets around and switch between runs without
touching anything in the repo.

## Running the standalone arCOG & plotting workflow

These scripts aren't part of `pipeline.nf` — run them by hand, in this order, from the
repo root (or adjust the `scripts/...` paths below to wherever you put them).

### arCOG annotation

1. Download the arCOG files, keeping the directory structure, from
   https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/
2. Search each genome's proteins against the arCOG sequences with DIAMOND:

   ```bash
   diamond blastp -q $protein -d $dbfile -p $n -f 6 $format -k 0 --ultra-sensitive --out $out
   ```
   where `$format` is
   `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp ppos`.
   Save each genome's result as `$genome_arcogs.tsv` next to `taxonomy.tsv`, per the
   directory layout in "Preparing your genome data directory" — `tool_arcog_filter.py`
   expects to find them there.

3. Filter to the single best hit per query:

   ```bash
   python3 scripts/tool_arcog_filter.py /path/to/taxonomy.tsv /path/to/save/filtered_arcog_file.tsv
   ```
   > The original documentation's example passes the raw DIAMOND output file as the
   > first argument. The script actually takes `taxonomy.tsv` as the first argument and
   > derives each genome's `$genome_arcogs.tsv` path from the directory it lives in —
   > pass `taxonomy.tsv`, not a DIAMOND output file, as shown above.

4. Map arCOG model/functional-category/definition and taxonomy onto the filtered hits:

   ```bash
   python3 scripts/tool_arcogs_map_info.py /path/to/arcog/files/directory/ /path/to/filtered_arcog_file.tsv /path/to/taxonomy.tsv /path/to/save/arcog_hit_map.tsv
   ```
   > Called `tool_arcog_map_info.py` (no "s") in the original documentation — the
   > script in this repo is `tool_arcogs_map_info.py`.

5. Count characterized/uncharacterized arCOG hits per genome and summarize per phylum:

   ```bash
   python3 scripts/tool_arcog_count_hits_per_genome.py /path/to/arcog_hit_map.tsv /path/to/taxonomy.tsv /path/to/save/arcog_counts_per_genome.tsv /path/to/save/arcog_counts_per_phylum.tsv
   ```

### Plotting

**Boxplot** of uncharacterized % per taxon — reshape, then plot:

```bash
# from the KEGG pipeline's countUnch output:
python3 scripts/tool_boxplot_prep.py pipeline /path/to/uncharacterized_counts_normalized_per_genome.tsv /path/to/taxonomy.tsv /path/to/save/boxplot_table.tsv

# or from the arCOG output above:
python3 scripts/tool_boxplot_prep.py arcog /path/to/arcog_hit_map.tsv /path/to/taxonomy.tsv /path/to/save/boxplot_table.tsv
```
then `Rscript scripts/r/tool_r_boxplot.R <boxplot_table.tsv> <output>`.

**Stacked barplot / heatmap** of module completeness — reshape, then plot:

```bash
python3 scripts/tool_prep_matrices_for_plotting.py /path/to/modules_presence_per_genome.tsv /path/to/taxonomy.tsv /path/to/module_completeness_percent.tsv
```
then `Rscript scripts/r/tool_r_stacked_barplot.R <input> <output>` for per-module
completeness, or `Rscript scripts/r/tool_r_heatmap.R <input> <output>` for the
marker-module heatmap (this one needs `data/kegg_modules_architecture.txt` set to the
`manual` module set — see "Choosing a KEGG module set" above).

## Bugs fixed in this pass

1. **`pipeline_analyze_characterized.py` and `pipeline_map_kegg_module.py`** parsed
   the KEGG module architecture file with `for line in mf: mr = mf.read().split(...)`.
   The `for` loop consumed the first line before `mf.read()` read (and this code
   parsed) the remainder of the file in one go, so the first line of the
   module-architecture file was silently dropped. Fixed to parse the whole file
   directly, matching `pipeline_prep_cutoffs.py` and
   `pipeline_calc_kegg_modules_completeness.py`, which never had this bug.
2. **`pipeline_calc_kegg_modules_completeness.py`** had a duplicated block (the M00617
   "combination module" fix + the `save_pergenome` write). The first copy capped the
   summed completeness at 1 (100%); an identical second copy immediately after
   recomputed the same sum *without* the cap and overwrote the output again, so the
   *uncapped* value was the one that actually ended up in
   `modules_presence_per_genome.tsv`. Removed the duplicate; the capped value is now
   the one that's saved.
3. **`r/tool_r_heatmap.R`** built its color palette with `length(breaks) - 1`, but the
   breaks vector was assigned to `breaks1` — `breaks` was never defined, so the script
   would fail with `object 'breaks' not found` on every run. Fixed to reference
   `breaks1`.
4. **`tool_prep_matrices_for_plotting.py`**: `count_taxa` was initialized but never
   populated, then indexed with `count_taxa[curr_ord]` later to turn raw counts into
   percentages — a `KeyError` on every run. Restored the missing block that increments
   `count_taxa[phyl]` per genome while parsing the taxonomy file (per your confirmation
   of the intended logic).

## Found but not yet resolved

None outstanding

## What changed in this cleanup pass (Python pipeline/tool scripts)

- `sys.argv[N]` positional reads → `argparse`, so every script now has `--help` and
  named arguments, and can be imported without executing
- Added a module docstring to every script explaining what it does and (for
  pipeline scripts) which Nextflow process calls it
- Wrapped each script's body in `main()` / `if __name__ == "__main__":`
- `open(f, "w+")` → `open(f, "w")` (cosmetic — these files are never read back)
- `pipeline.nf` updated to call `scripts/pipeline_*.py` and read from `data/`
- The R scripts now take CLI arguments, have a usage header comment, and dropped a
  no-op `mutate(Taxon = Taxon)` in the boxplot script
- Added `scripts/select_kegg_modules.sh` to automate the manual/original module-set
  copy step described in the original documentation (see "Choosing a KEGG module set")
- `ko_cutoffs.tsv` path is now a Nextflow param (`params.ko_cutoffs`), overridable with
  `--ko_cutoffs` instead of requiring you to overwrite `data/ko_cutoffs.tsv` by hand
  (see "Using custom KO cutoffs")

## Contact

Questions about the pipeline: mail@val-k.science (subject: "Question protein characterization
pipeline").

## Funding

This project has received funding from the European Research Council (ERC) under the
European Union's Horizon 2020 research and innovation programme (grant agreement no.
803768).

## Note

The original documentation has been updated using Claude Sonnet 5, with a subsequent quality control by the author of the pipeline.

This tool was produced as a part of the publication "Navigating the archaeal frontier: insights and projections from bioinformatic pipelines" by Karavaeva et al., which is available as open-access at https://doi.org/10.3389/fmicb.2024.1433224
