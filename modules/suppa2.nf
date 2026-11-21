/*
 * Alternative splicing events and PSI (SUPPA2)
 */

process SUPPA2 {
    tag "suppa2"
    publishDir "${params.outdir}/suppa2", mode: params.publish_dir_mode, pattern: "{events_*.ioe,psi/*.psi}"

    input:
    path gtf          // Bambu extended_annotations.gtf
    path counts       // Bambu counts_transcript.txt (TXNAME, GENEID, samples...)

    output:
    path "events_*.ioe"              , emit: ioe
    path "psi/*_samples.psi"         , emit: psi_per_sample
    path "psi/*_samples_filtered.psi", emit: psi_filtered
    path "versions.yml"              , emit: versions

    script:
    """
    # SUPPA2 expression file (transcript IDs kept exactly as in the Bambu GTF; bin/suppa2_expression.awk)
    awk -f "\$(command -v suppa2_expression.awk)" ${counts} > expression.tsv
    echo "Expression file: \$(( \$(wc -l < expression.tsv) - 1 )) transcripts; samples: \$(head -n 1 expression.tsv | tr '\\t' ' ')"

    suppa.py generateEvents -i ${gtf} -f ioe -e ${params.suppa2_events} -o events

    mkdir -p psi
    for ioe_file in events_*_strict.ioe; do
        [ -s "\$ioe_file" ] || continue
        prefix=\$(basename "\$ioe_file" _strict.ioe)
        suppa.py psiPerEvent -i "\$ioe_file" -e expression.tsv -o "psi/\${prefix}_samples"

        # Keep events with a PSI value (not nan) in at least 2 samples (bin/psi_min_samples.awk)
        psi_file="psi/\${prefix}_samples.psi"
        if [ -s "\$psi_file" ]; then
            awk -v min=2 -f "\$(command -v psi_min_samples.awk)" "\$psi_file" > "psi/\${prefix}_samples_filtered.psi"
            echo "\$prefix: \$(( \$(wc -l < "\$psi_file") - 1 )) events, \$(( \$(wc -l < "psi/\${prefix}_samples_filtered.psi") - 1 )) with PSI in >= 2 samples"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python3 -c "import importlib.metadata as m; print(m.version('suppa'))" 2>/dev/null || echo 2.4)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p psi
    touch events_SE_strict.ioe psi/events_SE_samples.psi psi/events_SE_samples_filtered.psi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}

/*
 * Differential splicing (SUPPA2 diffSplice), each non-reference condition vs the reference.
 * Runs only with >= 2 replicates in the reference and in the compared condition (decided in main.nf).
 * As in the SUPPA2 workflow, the events of all types are merged and tested together (one empirical
 * background per comparison). PSI values come from the count-based psiPerEvent output (PSI is a
 * within-gene ratio, so counts and CPM give the same PSI); the expression files passed with -e are the
 * Bambu CPM values. bin/suppa2_results.py adds a genome-wide Benjamini-Hochberg padj.
 */

process SUPPA2_DIFFSPLICE {
    tag "diffsplice"
    publishDir "${params.outdir}/suppa2/diffsplice", mode: params.publish_dir_mode, pattern: "*.{dpsi,psivec,tsv}"

    input:
    path ioe            // events_<type>_strict.ioe
    path psi            // events_<type>_samples.psi
    path cpm            // Bambu CPM_transcript.txt (TXNAME, GENEID, samples...)
    path design         // design.csv (sample,condition[,covariates])
    val  reference      // reference condition
    val  conditions     // non-reference conditions to compare with the reference

    output:
    path "*.dpsi"                                 , emit: dpsi, optional: true
    path "*.psivec"                               , emit: psivec, optional: true
    path "suppa2_diffsplice_results.tsv"          , emit: results
    path "suppa2_diffsplice_significant.tsv"      , emit: significant
    path "versions.yml"                           , emit: versions

    script:
    def args = params.suppa2_diff_args ?: ''
    def conds = conditions.collect { c -> "'${c}'" }.join(' ')
    """
    awk -f "\$(command -v suppa2_expression.awk)" ${cpm} > expression_cpm.tsv
    split_cols() {  # split_cols <condition> <matrix> <output>
        awk -v cond="\$1" -f "\$(command -v split_condition.awk)" fileno=1 ${design} fileno=2 "\$2" > "\$3"
    }

    # Merge the event types that have events (a header-only PSI file would make diffSplice fail)
    : > events_all.ioe; : > events_all.psi
    for ioe_file in events_*_strict.ioe; do
        type=\$(basename "\$ioe_file" _strict.ioe); type=\${type#events_}
        psi_file="events_\${type}_samples.psi"
        if [ ! -s "\$psi_file" ] || [ "\$(wc -l < "\$psi_file")" -le 1 ]; then
            echo "\$type: no events - not tested"
            continue
        fi
        [ -s events_all.ioe ] || head -n 1 "\$ioe_file" > events_all.ioe
        [ -s events_all.psi ] || head -n 1 "\$psi_file" > events_all.psi
        tail -n +2 "\$ioe_file" >> events_all.ioe
        tail -n +2 "\$psi_file" >> events_all.psi
        echo "\$type: \$(( \$(wc -l < "\$psi_file") - 1 )) events"
    done

    dpsi_args=""
    if [ -s events_all.psi ]; then
        for cond in '${reference}' ${conds}; do
            mkdir -p "split/\$cond"
            split_cols "\$cond" expression_cpm.tsv "split/\$cond/expression.tpm"
            split_cols "\$cond" events_all.psi "split/\$cond/\$cond.psi"
        done
        for cond in ${conds}; do
            # The PSI files are named after the conditions: SUPPA2 uses the file names in the .dpsi header
            prefix="\${cond}_vs_${reference}"
            suppa.py diffSplice ${args} -i events_all.ioe \\
                -p "split/${reference}/${reference}.psi" "split/\$cond/\$cond.psi" \\
                -e "split/${reference}/expression.tpm" "split/\$cond/expression.tpm" \\
                -o "\$prefix"
            # SUPPA2 names the columns after the PSI file names up to the first '.', which truncates
            # condition names containing '.' (e.g. '2.trt' -> '2_1'): write the full names
            n_ref=\$(head -n 1 "split/${reference}/${reference}.psi" | awk -F '\t' '{ print NF }')
            awk -v ref='${reference}' -v cond="\$cond" -v nref="\$n_ref" 'BEGIN { FS = OFS = "\t" }
                NR == 1 { for (i = 1; i <= NF; i++) \$i = (i <= nref ? ref "_" i : cond "_" (i - nref)) } { print }' \\
                "\$prefix.psivec" > psivec.tmp && mv psivec.tmp "\$prefix.psivec"
            awk -v ref='${reference}' -v cond="\$cond" 'BEGIN { FS = OFS = "\t" }
                NR == 1 { \$1 = ref "-" cond "_dPSI"; \$2 = ref "-" cond "_p-val" } { print }' \\
                "\$prefix.dpsi" > dpsi.tmp && mv dpsi.tmp "\$prefix.dpsi"
            dpsi_args="\$dpsi_args \$prefix=\$prefix.dpsi"
        done
    else
        echo "No splicing events with PSI values: diffSplice not run"
    fi

    # All tested events with a genome-wide padj, and the significant ones (padj < suppa2_alpha,
    # |dPSI| > suppa2_dpsi_cutoff; dPSI = condition - reference)
    suppa2_results.py --alpha ${params.suppa2_alpha} --dpsi ${params.suppa2_dpsi_cutoff} \$dpsi_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python3 -c "import importlib.metadata as m; print(m.version('suppa'))" 2>/dev/null || echo 2.4)
    END_VERSIONS
    """

    stub:
    """
    touch treatment_vs_control.dpsi treatment_vs_control.psivec
    printf 'comparison\\tevent_type\\tevent_id\\tgene_id\\tdPSI\\tpvalue\\tpadj\\n' > suppa2_diffsplice_results.tsv
    cp suppa2_diffsplice_results.tsv suppa2_diffsplice_significant.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
