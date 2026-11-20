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
 * PSI values come from the count-based psiPerEvent output (PSI is a within-gene ratio, so counts and
 * CPM give the same PSI); the expression files passed with -e are the Bambu CPM values.
 */

process SUPPA2_DIFFSPLICE {
    tag "diffsplice"
    publishDir "${params.outdir}/suppa2/diffsplice", mode: params.publish_dir_mode, pattern: "*.{dpsi,psivec,tsv}"

    input:
    path ioe            // events_<type>_strict.ioe
    path psi            // events_<type>_samples.psi
    path cpm            // Bambu CPM_transcript.txt (TXNAME, GENEID, samples...)
    path design         // design.csv (sample,condition[,covariates])
    val reference       // reference condition
    val conditions      // non-reference conditions to compare with the reference

    output:
    path "*.dpsi"                                 , emit: dpsi, optional: true
    path "*.psivec"                               , emit: psivec, optional: true
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

    sig_args=""
    for cond in '${reference}' ${conds}; do
        mkdir -p "split/\$cond"
        split_cols "\$cond" expression_cpm.tsv "split/\$cond/expression.tpm"
    done
    for ioe_file in events_*_strict.ioe; do
        [ -s "\$ioe_file" ] || continue
        type=\$(basename "\$ioe_file" _strict.ioe); type=\${type#events_}
        psi_file="events_\${type}_samples.psi"
        [ -s "\$psi_file" ] || continue
        for cond in '${reference}' ${conds}; do
            split_cols "\$cond" "\$psi_file" "split/\$cond/\$type.psi"
        done
        for cond in ${conds}; do
            # Files are named after the conditions: SUPPA2 uses the PSI file names in the .dpsi header
            mkdir -p "cmp_\$cond"
            ln -sf "../split/${reference}/\$type.psi" "cmp_\$cond/${reference}.psi"
            ln -sf "../split/\$cond/\$type.psi" "cmp_\$cond/\$cond.psi"
            prefix="\${cond}_vs_${reference}_\${type}"
            suppa.py diffSplice ${args} -i "\$ioe_file" \\
                -p "cmp_\$cond/${reference}.psi" "cmp_\$cond/\$cond.psi" \\
                -e "split/${reference}/expression.tpm" "split/\$cond/expression.tpm" \\
                -o "\$prefix"
            sig_args="\$sig_args cmp=\${cond}_vs_${reference} type=\$type \$prefix.dpsi"
        done
    done

    # Combined table of significant events (p < suppa2_alpha, |dPSI| > suppa2_dpsi_cutoff; dPSI = condition - reference)
    awk -v alpha=${params.suppa2_alpha} -v dpsi=${params.suppa2_dpsi_cutoff} -f "\$(command -v suppa2_significant.awk)" \$sig_args < /dev/null \\
        > suppa2_diffsplice_significant.tsv
    echo "Significant events (p < ${params.suppa2_alpha}, |dPSI| > ${params.suppa2_dpsi_cutoff}): \$(( \$(wc -l < suppa2_diffsplice_significant.tsv) - 1 ))"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        suppa: \$(python3 -c "import importlib.metadata as m; print(m.version('suppa'))" 2>/dev/null || echo 2.4)
    END_VERSIONS
    """

    stub:
    """
    touch treatment_vs_control_SE.dpsi treatment_vs_control_SE.psivec
    printf 'comparison\\tevent_type\\tevent_id\\tgene_id\\tdPSI\\tp_value\\n' > suppa2_diffsplice_significant.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: "stub"
    END_VERSIONS
    """
}
