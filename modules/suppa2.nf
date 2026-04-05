process SUPPA2 {
    tag "suppa2_quant"
    publishDir "${params.outdir}/suppa2", mode: 'copy'

    input:
    path gtf                     
    path counts                  

    output:
    path "events_*.ioe", emit: ioe
    path "psi/*.psi", emit: psi_per_sample
    path "psi/*_filtered.psi", emit: psi_filtered   

    script:
    """
    rm -rf events_* psi cleaned_counts.txt || true

    # 1. Clean expression file + strip version numbers from IDs
    python3 -c "
import sys

with open('${counts}', 'r') as f_in, open('cleaned_counts.txt', 'w') as f_out:
    header_line = next(f_in).rstrip('\\n').split('\\t')
    sample_names = header_line[2:]
    sample_names = [
        s.replace('.sorted', '')
         .replace('.bam', '')
         .replace('.bai', '')
         .replace('.Aligned', '')
         .strip()
        for s in sample_names
    ]
    
    if not sample_names:
        sys.exit('Error: No sample columns detected in header')
    
    f_out.write('\\t'.join(sample_names) + '\\n')
    
    for line in f_in:
        parts = line.rstrip('\\n').split('\\t')
        n_samples = len(sample_names)
        if len(parts) >= n_samples + 2:
            tx_id = parts[0].split('.')[0]          # <- remove version suffix
            out_row = [tx_id] + parts[2 : 2 + n_samples]
            if len(out_row) == n_samples + 1:
                f_out.write('\\t'.join(out_row) + '\\n')

print('Detected sample names:', ', '.join(sample_names))
print('Number of samples:', len(sample_names))
"

    # Debug: show what we produced
    echo "=== cleaned_counts.txt (first 10 lines) ==="
    head -n 10 cleaned_counts.txt || echo "File is empty!"
    echo ""
    echo "Field counts:"
    awk -F'\\t' '{print \"Line \" NR \": \" NF \" fields\"}' cleaned_counts.txt | head -n 10

    # 2. Generate events
    suppa.py generateEvents -i ${gtf} -f ioe -e SE SS MX RI FL -o events

    # 3. Calculate PSI
    mkdir -p psi

    for ioe_file in events_*_strict.ioe; do
        if [ -s "\$ioe_file" ]; then
            prefix=\$(basename \$ioe_file _strict.ioe)
            echo "Calculating PSI for: \$prefix"

            suppa.py psiPerEvent \
                -i "\$ioe_file" \
                -e cleaned_counts.txt \
                -o "psi/\${prefix}_samples"

                        # 4. Filter: keep events with non-NaN in ≥ 2 samples
            psi_file="psi/\${prefix}_samples.psi"
            filtered_file="psi/\${prefix}_samples_filtered.psi"

            if [ -s "\$psi_file" ]; then
                awk -v samples=\$(awk -F'\\t' 'NR==1 {print NF}' "\$psi_file") '
                NR==1 {print; next}
                {
                    non_nan = 0
                    for (i=2; i<=samples; i++) {          # use the passed variable
                        if (\$i != "nan" && \$i != "NA" && \$i != "") non_nan++
                    }
                    if (non_nan >= 2) print
                }' "\$psi_file" > "\$filtered_file"

                echo "Filtered \$psi_file → \$filtered_file (≥2 non-NaN PSI values)"
                echo "  Original lines: \$(wc -l < "\$psi_file")"
                echo "  Filtered lines: \$(wc -l < "\$filtered_file")"
            else
                echo "No PSI file produced for \$prefix"
            fi
        fi
    done
    """
}