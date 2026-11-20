/*
 * Normalised design table (sample,condition) for the R steps.
 *
 * Built by a process (not collectFile) so that the file is cached: a collectFile output is
 * re-created on every launch, which would make -resume re-run every step that reads it.
 */

process MAKE_DESIGN {
    tag "design"

    input:
    val rows            // ["sample,condition", ...]; names are validated by the samplesheet schema

    output:
    path "design.csv", emit: csv

    script:
    def lines = (['sample,condition'] + rows).collect { r -> "'${r}'" }.join(' ')
    """
    printf '%s\\n' ${lines} > design.csv
    """

    stub:
    def lines = (['sample,condition'] + rows).collect { r -> "'${r}'" }.join(' ')
    """
    printf '%s\\n' ${lines} > design.csv
    """
}
