process TE_ANALYSIS {
    tag "TE_analysis"
    publishDir "${params.outdir}/TE", mode: 'copy'
    
    input:
    path iso_gtf
    path te_gtf
    path counts
    path junctions_bed
    path ref_gtf
    
    output:
    path "TE_family_counts.tsv"       , emit: te_counts
    path "TE_instance_counts.tsv"     , emit: te_instances
    path "isoform_TE_exonization.tsv" , emit: te_exonization
    
    script:
    """
    set -euo pipefail

    ############################################################
    # 1. DEFINE FORBIDDEN ZONE (All annotated Exons/UTRs)
    ############################################################
    # Extract valid exons only (require strand, valid coords)
    awk -F'\\t' '
      \$3 == "exon" && NF >= 8 && \$4 >= 1 && \$5 > \$4 && \$7 != "" {
        print \$1 "\\t" (\$4 - 1) "\\t" \$5 "\\t.\\t.\\t" \$7
      }
    ' ${ref_gtf} > ref_exons_raw.bed

    # Debug: fail early if no exons found
    if [ ! -s ref_exons_raw.bed ]; then
      echo "ERROR: No valid exons extracted from ref_gtf = ${ref_gtf}"
      echo "First 20 lines of ref_gtf:"
      head -n 20 ${ref_gtf}
      echo "Check: does the file have \$3 == \"exon\" records? Is it tab-delimited?"
      exit 1
    fi

    # Sort (should now succeed)
    bedtools sort -i ref_exons_raw.bed > ref_exons_sorted.bed || {
      echo "Error sorting ref_exons_raw.bed"
      head -n 10 ref_exons_raw.bed
      exit 1
    }

    ############################################################
    # 2. PREP GENOMIC BEDS
    ############################################################
    # Isoform exons
    awk '\$3=="exon"' ${iso_gtf} | \
      sed -E 's/.*transcript_id "([^"]+)".*/\\1\\t&/' | \
      awk -v OFS="\\t" '{print \$2, \$5-1, \$6, \$1, ".", \$8}' > iso_exons.bed

    # Pure intronic exons
    bedtools intersect -v -a iso_exons.bed -b ref_exons_sorted.bed > iso_pure_intronic_exons.bed

    # TE BED – prevent negative / invalid coordinates
    awk '\$3=="exon" || \$3=="transcript"' ${te_gtf} | \
      sed -E 's/.*(gene_id|transcript_id) "([^"]+)".*/\\2\\t&/' | \
      awk -v OFS="\\t" '{
        if (NF < 8 || \$5 == "" || \$6 == "" || \$8 == "") next;
        start = \$5 - 1;
        if (start < 0) start = 0;
        if (start >= \$6) next;
        print \$2, start, \$6, \$1, ".", \$8
      }' | sort -u > te.bed

    # Safety check: te.bed must not be empty
    if [ ! -s te.bed ]; then
      echo "ERROR: te.bed is empty"
      echo "First 10 lines of te_gtf:"
      head -n 10 ${te_gtf}
      exit 1
    fi

    # Optional: check for any remaining negative starts (using safe awk)
    NEGATIVE_CHECK=\$(awk '\$2 < 0 {print \$0; exit 1}' te.bed || true)
    if [ -n "\$NEGATIVE_CHECK" ]; then
      echo "WARNING: Negative start coordinates still found in te.bed"
      echo "\$NEGATIVE_CHECK" | head -n 5
    fi

    ############################################################
    # 3. INTERSECTIONS
    ############################################################
    bedtools intersect -s -f 0.9 -wa -wb -a iso_pure_intronic_exons.bed -b te.bed > exon_overlap_strict.txt

    awk '\$3=="transcript"' ${iso_gtf} | \
      sed -E 's/.*transcript_id "([^"]+)".*/\\1\\t&/' | \
      awk -v OFS="\\t" '{print \$2, \$5-1, \$6, \$1, ".", \$8}' > isoforms.bed
    bedtools intersect -s -wa -wb -a isoforms.bed -b te.bed > exon_overlap_general.txt

    grep 'exon_number "1"' ${iso_gtf} | \
      sed -E 's/.*transcript_id "([^"]+)".*/\\1\\t&/' | \
      awk -v OFS="\\t" '{print \$2, \$5-1, \$6, \$1, ".", \$8}' > first_exons.bed
    bedtools intersect -s -wa -wb -a first_exons.bed -b te.bed > prom_overlap.txt

    awk '\$3=="exon"' ${iso_gtf} | \
      sed -E 's/.*transcript_id "([^"]+)".*/\\1\\t&/' > all_exons_tmp.txt
    awk -v OFS="\\t" '{
        id=\$1; chr=\$2; start=\$5; end=\$6; strand=\$8;
        if (strand == "+") { if (end > max_end[id]+0) { max_end[id]=end; line[id]=chr"\\t"start-1"\\t"end"\\t"id"\\t.\\t"strand } }
        else { if (start < min_start[id] || min_start[id]=="") { min_start[id]=start; line[id]=chr"\\t"start-1"\\t"end"\\t"id"\\t.\\t"strand } }
    } END { for (i in line) print line[i] }' all_exons_tmp.txt > last_exons.bed
    bedtools intersect -s -wa -wb -a last_exons.bed -b te.bed > polya_overlap.txt

    ############################################################
    # 4. SPLICE JUNCTION VALIDATION (Strict Exonization Only)
    ############################################################
    awk -v OFS="\\t" '{print \$1, \$2-5, \$2+5, \$4}' exon_overlap_strict.txt > left_sites.bed
    awk -v OFS="\\t" '{print \$1, \$3-5, \$3+5, \$4}' exon_overlap_strict.txt > right_sites.bed

    bedtools intersect -u -a left_sites.bed -b ${junctions_bed} > left_ok.bed || touch left_ok.bed
    bedtools intersect -u -a right_sites.bed -b ${junctions_bed} > right_ok.bed || touch right_ok.bed

    awk '{print \$4}' left_ok.bed right_ok.bed | sort | uniq -c | awk '\$1==2 {print \$2}' > validated_ids.txt

    ############################################################
    # 5. JOINING AND SUMMARIES
    ############################################################
    head -n 1 ${counts} | sed 's/\\.sorted//g; s/\\.bam//g' > counts_cleaned.tsv
    tail -n +2 ${counts} >> counts_cleaned.tsv

    SAMPLES=\$(head -n 1 counts_cleaned.tsv | cut -f3-)
    NUM_SAMPLES=\$(head -n 1 counts_cleaned.tsv | awk '{print NF-2}')
    HEADER="transcript_id\\tgene_id\\ttx_chr\\ttx_start\\ttx_end\\tTE_id\\tTE_chr\\tTE_start\\tTE_end\\t\$SAMPLES"

    echo -e "\$HEADER" > isoform_TE_exonization.tsv
    awk -F'\\t' '
        ARGIND==1 { ok[\$1]=1; next }
        ARGIND==2 { data[\$4] = \$1"\\t"\$2"\\t"\$3"\\t"\$10"\\t"\$7"\\t"\$8"\\t"\$9; next }
        ARGIND==3 { if (\$1 in data && \$1 in ok) { printf "%s\\t%s\\t%s", \$1, \$2, data[\$1]; for(i=3; i<=NF; i++) printf "\\t%s", \$i; printf "\\n" } }
    ' validated_ids.txt exon_overlap_strict.txt counts_cleaned.tsv >> isoform_TE_exonization.tsv

    awk -F'\\t' 'NR==FNR { data[\$4]=\$1"\\t"\$2"\\t"\$3"\\t"\$10"\\t"\$7"\\t"\$8"\\t"\$9; next } \
        (\$1 in data) { printf "%s\\t%s\\t%s", \$1, \$2, data[\$1]; for(i=3;i<=NF;i++) printf "\\t%s", \$i; printf "\\n" }' \
        exon_overlap_general.txt counts_cleaned.tsv > all_joined_temp.tsv

    echo -e "TE_Family\\t\$SAMPLES" > TE_family_counts.tsv
    awk -F'\\t' -v n="\$NUM_SAMPLES" '{
        fam=\$6; sub(/_dup[0-9]+/, "", fam); for(i=1; i<=n; i++) sum[fam,i] += \$(i+9); ids[fam]=1
    } END { for (f in ids) { printf "%s", f; for(j=1; j<=n; j++) printf "\\t%s", (sum[f,j]?sum[f,j]:0); printf "\\n" } }' all_joined_temp.tsv >> TE_family_counts.tsv

    echo -e "TE_id\\tTE_chr\\tTE_start\\tTE_end\\t\$SAMPLES" > TE_instance_counts.tsv
    awk -F'\\t' -v n="\$NUM_SAMPLES" '{
        inst=\$6"\\t"\$7"\\t"\$8"\\t"\$9; for(i=1; i<=n; i++) inst_sum[inst,i] += \$(i+9); locs[inst]=1
    } END { for (l in locs) { printf "%s", l; for(j=1; j<=n; j++) printf "\\t%s", (inst_sum[l,j]?inst_sum[l,j]:0); printf "\\n" } }' all_joined_temp.tsv >> TE_instance_counts.tsv

    rm all_joined_temp.tsv
    """
}