process parseInterpro {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    output:
        path "merged_parsed_interpro.tsv"

    """
        python3 $projectDir/pipeline_parse_interpro.py "$launchDir/taxonomy.tsv" "merged_parsed_interpro.tsv"
    """
}

process filterHmmer {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    output:
        path "merged_filtered_hmmer.tsv"

    """
    python3 $projectDir/pipeline_filter_hmmer.py "$launchDir/taxonomy.tsv" "merged_filtered_hmmer.tsv" "$projectDir/ko_cutoffs.tsv"
    """
}

process createAccgen {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    output:
        path "merged_accgen.tsv"

    """
    python3 $projectDir/pipeline_prepare_acc_gen.py "$launchDir/taxonomy.tsv" "merged_accgen.tsv"
    """
}

process mergeTables {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path interproout
        path hmmerout
        path accgenout
    output:
        path "merged_full_table.tsv"

    """
    python3 $projectDir/pipeline_merge_table.py $hmmerout $interproout $accgenout "$launchDir/taxonomy.tsv" "$projectDir/ko00001_level_C.tsv" "merged_full_table.tsv"
    """
}

process summaryAnno {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path fullout
    output:
        path "total_presence_annotations_matrix.tsv"
        path "total_presence_annotations_counts.tsv"
        path "total_presence_annotations_counts_summary.tsv"

    """
    python3 $projectDir/pipeline_calc_general_presence_annotations.py $fullout "total_presence_annotations_matrix.tsv" "total_presence_annotations_counts.tsv" "total_presence_annotations_counts_summary.tsv"
    """
}

process classifyUnch {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path fullout
        path accgenout
    output:
        path "uncharacterized_proteins.tsv", emit: unch_sv
        path "characterized_proteins.tsv", emit: ch_sv
        path "classify_summary.txt"
    """
    python3 $projectDir/pipeline_classify_un_characterized.py $fullout "$projectDir/tigr_to_remove.txt" "$projectDir/panther_to_remove.txt" $accgenout "$launchDir/taxonomy.tsv" "uncharacterized_proteins.tsv" "characterized_proteins.tsv" "classify_summary.txt"
    """
}

process countUnch {
    publishDir (
        path: "./output/uncharacterized_hits_counts",
        mode: "copy",
        overwrite: true
    )
    input:
        path unchout
    output:
        path("uncharacterized_counts_normalized_per_*.tsv", arity: "1..*")
    """
    python3 $projectDir/pipeline_normalize_uncharacterized_counts.py "$launchDir/taxonomy.tsv" $unchout
    """
}

process mapModule {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path chout
    output:
        path "characterized_proteins_with_kegg_modules.tsv", emit: chsv_mod
        path "characterized_grouped_by_kegg_modules.tsv", emit: ch_group
    """
    python3 $projectDir/pipeline_map_kegg_module.py $chout "$projectDir/kegg_modules_architecture.txt" "characterized_proteins_with_kegg_modules.tsv" "characterized_grouped_by_kegg_modules.tsv"
    """
}

process analyzeChar {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path chgroup
    output:
        path "characterized_hits_module_fullness_presence_absence.tsv", emit: modfull
        path "yes_hit_modules.txt"
        path "no_hit_modules.txt", emit: nohit
        path "hits_per_module", emit: permoddir
        path("hits_per_module/no_hits/characterized_hits_module_fullness_presence_absence_*.tsv", arity: "1..*")
        path("hits_per_module/yes_hits/characterized_hits_module_fullness_presence_absence_*.tsv", arity: "1..*")
    """
    python3 $projectDir/pipeline_analyze_characterized.py $chgroup "$projectDir/kegg_modules_architecture.txt" "$launchDir/taxonomy.tsv" "$projectDir/kegg_module.txt" "$projectDir/all_ko_names.tsv" "characterized_hits_module_fullness_presence_absence.tsv" "no_hit_modules.txt" "yes_hit_modules.txt"
    """
}

process prepCutoffs {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    output:
        path "kegg_cutoffs_per_module.txt"
    """
    python3 $projectDir/pipeline_prep_cutoffs.py "$projectDir/kegg_modules_architecture.txt" "$projectDir/ko_cutoffs.tsv" "kegg_cutoffs_per_module.txt"
    """
}

process calcModule {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path accgenout
        path nohitf
        path modulefull
        path moddir
        path kegg_cutoff
    output:
        path "complex_presence_per_genome.tsv"
        path "alternative_presence_per_genome.tsv"
        path "modules_presence_per_genome.tsv"
        path "module_completeness_percent.tsv"
    """
    python3 $projectDir/pipeline_calc_kegg_modules_completeness.py "$projectDir/kegg_modules_architecture.txt" $accgenout $nohitf $kegg_cutoff "$launchDir/taxonomy.tsv" $modulefull "complex_presence_per_genome.tsv" "alternative_presence_per_genome.tsv" "modules_presence_per_genome.tsv" "module_completeness_percent.tsv" $moddir
    """
}

process calcNonmodule {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path chsvmod
    output:
        path "nonmodule_ko_presence_per_genome.tsv"
        path "nonmodule_ko_presence_percentage.tsv"
    """
    python3 $projectDir/pipeline_calc_presence_nonmodule_kos.py $chsvmod "$launchDir/taxonomy.tsv" "$projectDir/nonmodule_kos_list.txt" "nonmodule_ko_presence_per_genome.tsv" "nonmodule_ko_presence_percentage.tsv"
    """
}

process syntenyUnch {
    publishDir (
        path: "./output",
        mode: "copy",
        overwrite: true
    )
    input:
        path unchout
        path fullout
    output:
        path "synteny_merged_table.tsv"
        path "synteny_indexed_table.tsv"
        path "synteny_neighbors.tsv"
        path "synteny_count_uncharacterized_blocks.tsv"
        path "synteny_count_uncharacterized_blocks_more_than_ten.tsv"
        path "synteny_count_uncharacterized_blocks_per_genome.tsv"
    """
    python3 $projectDir/pipeline_synteny_uncharacterized.py "$launchDir/taxonomy.tsv" $unchout $fullout "synteny_merged_table.tsv" "synteny_indexed_table.tsv" "synteny_neighbors.tsv" "synteny_count_uncharacterized_blocks.tsv" "synteny_count_uncharacterized_blocks_more_than_ten.tsv" "synteny_count_uncharacterized_blocks_per_genome.tsv"
    """
}

workflow {
    parseInterpro()
    filterHmmer()
    createAccgen()
    prepCutoffs()
    mergeTables(parseInterpro.out,filterHmmer.out,createAccgen.out)
    summaryAnno(mergeTables.out)
    classifyUnch(mergeTables.out,createAccgen.out)
    countUnch(classifyUnch.out.unch_sv)
    mapModule(classifyUnch.out.ch_sv)
    analyzeChar(mapModule.out.ch_group)
    calcModule(createAccgen.out,analyzeChar.out.nohit,analyzeChar.out.modfull,analyzeChar.out.permoddir,prepCutoffs.out)
    calcNonmodule(mapModule.out.chsv_mod)
    syntenyUnch(classifyUnch.out.unch_sv,mergeTables.out)
}
