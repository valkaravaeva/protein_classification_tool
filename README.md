#steps to run the tool: # tested with the nextflow version 24.04.4-17.0.6; python version 3.11.9

Biopython has to be installed! It is a dependency!

1. Run Interpro per genome, example:

    interproscan.sh -cpu $n -f tsv -iprlookup -i $protein -o $out #version

    where $n is the number of cpus; $protein is your fasta file; $out is the path to savefile

2. Run Hmmer with kofam models (not kofam scan!) per genome, example:

    Go to https://www.genome.jp/tools/kofamkoala/ and download the profiles and ko_list

    run hmmer #version 3.4

    hmmsearch --cpu=$n -o $errlog -T0 --tblout=$tbl $hmm $protein

    where $n is the number of cpus; $protein is your genome faa file; $tbl is the path to savefile; $hmm is the concatenated profiles file

3. Prepare the taxonomy table without a header with one line per genome; tab-separated; with following columns: genome_accession; domain; phylum; class; order; family; genus; species; strain

4. Put everything into one directory (folder)

    the directory should contain:

    the genome protein fasta file named: $genome.faa
    the genome feature table or gff file (for synteny analysis), named: $genome.gff or $genome_feature_table.txt
    the interpro output: $genome_interpro.tsv
    the hmmer output: $genome_hmmer.tsv
    taxonomy file named: taxonomy.tsv
    (optional if running arcog scripts, see 6.): $genome_arcogs.tsv

    replace the "$genome" in the filenames with your genome accession (remove the $)

5. Install Nextflow, see https://www.nextflow.io/docs/latest/install.html

    Download the pipeline files (clone the Git directory)

    make sure Nextflow is working

    go to the directory where your genome and taxonomy, etc files are

    from that directory run the pipeline via "nextflow run /path/to/the/tool/directory/pipeline.nf" #put absolute path to the directory where the tool/pipeline files are

    if you want to run the pipeline with the manual modules (as described in the publication), go to the directory where the pipeline is and create a copy of the "manual_modules_architecture.txt" file, then rename it to "kegg_modules_architecture.txt" (to use the original modules again, copy the "kegg_modules_architecture_original.txt" and rename it to "kegg_modules_architecture.txt")

    if you want to provide custom KO cutoffs/thresholds, save them as "ko_cutoffs.tsv" in the pipeline directory (tabular format, KO\tcutoff per line, no header)

    the output files will be in the output/ directory, the work/ directory is automatically created by Nextflow, it can be deleted afterwards

6. the Git repository also contains:

6.1. example scripts to filter and analyse arcogs
    6.1.1. download arcog files, keep the directory structure as it is: https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/
    6.1.2. run diamond your genome vs arcogs (find arcog sequences files where you downloaded them)
            example run diamond:
            diamond blastp -q $protein -d $dbfile -p $n -f 6 $format -k 0 --ultra-sensitive --out $out
            where $protein is your genome faa file, $dbfile is the arcog sequence file, $n is # CPUs, $out is your outputfile, and $format is "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp ppos"
    6.1.3. filter the diamond output using tool_arcog_filter.py,
            example: python3 tool_arcog_filter.py /path/to/diamond/arcog/output.out /path/to/save/filtered_arcog_file.tsv
    6.1.4 map arcog information to the hits,
            example: python3 tool_arcog_map_info.py /path/to/arcog/files/directory/ /path/to/filtered_arcog_file.tsv /path/to/taxonomy.tsv /path/to/save/arcog_hit_map.tsv
    6.1.5. analyze the characterized/uncharacterized ratio via arcog, using tool_arcog_count_hits_per_genome.py,
            example: python3 tool_arcog_count_hits_per_genome.py /path/to/arcog_hit_map.tsv /path/to/taxonomy.tsv /path/to/save/arcog_counts_per_genome.tsv /path/to/save/arcog_counts_per_phylum.tsv

6.2. example plotting scripts and the python prep scripts
    6.2.1. for boxplot:
        run tool_boxplot_prep.py,
            example for pipeline: python3 pipeline /path/to/uncharacterized_counts_normalized_per_genome.tsv /path/to/taxonomy.tsv /path/to/save/boxplot_table.tsv
            example for arcogs: python3 arcog /path/to/arcog_hit_map.tsv /path/to/taxonomy.tsv /path/to/save/boxplot_table.tsv
        then use the R script tool_r_boxplot.R, replace the paths to your respective paths
    6.2.2. for stacked plots and heatmaps:
        run tool_prep_matrices_for_plotting.py,
            example: python3 tool_prep_matrices_for_plotting.py /path/to/modules_presence_per_genome.tsv /path/to/taxonomy.tsv /path/to/module_completeness_percent.tsv
            to plot stacked barplot per module: use tool_r_stacked_barplot.R, replace the paths to your respective paths
            to plot heatmap of markers modules: use tool_r_heatmap.R, replace the paths to your respective paths

For questions, email: val.karavaeva@univie.ac.at, with the "Question uncharacterized pipeline" in the topic
