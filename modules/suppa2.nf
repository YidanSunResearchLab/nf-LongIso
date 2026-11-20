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
