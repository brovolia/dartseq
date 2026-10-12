process BULLSEYE_R_GLM {
    tag "$meta.id"
    label 'process_medium'

    container 'community.wave.seqera.io/library/r-aod_r-base_r-data.table_r-dplyr_pruned:50e53b1dbf9130ec'

    input:
    tuple val(meta), path(cov_tsv), path(mut_tsv), path(observations_tsv)
    val(bullseye_glm_mock)
    path(bullseye_functions_r_file)
    val(bullseye_glm_design)
    val(bullseye_glm_min_cov)
    val(bullseye_glm_link)

    output:
    tuple val(meta), path("${prefix}.glm.tsv"), emit: tsv
    tuple val("${task.process}"), val('bullseye_r'), val("${bullseye_r_version}"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.bullseye"
    bullseye_r_version = bullseye_glm_mock ? 'mock' : 'local'

    if (bullseye_glm_mock) {
        """
        cat > ${prefix}.glm.tsv << 'EOF'
feature	pvalue	lfc
mock_site	1	0
EOF
        """
    } else {
        """
        cat > run_bullseye_glm.R << 'EOF'
        suppressPackageStartupMessages({
            library(data.table)
            library(dplyr)
            library(aod)
        })

        source("${bullseye_functions_r_file}")

        cov <- fread("${cov_tsv}")
        mut <- fread("${mut_tsv}")
        observations <- as.data.frame(fread("${observations_tsv}"))

        if (!"sample" %in% colnames(observations)) {
            stop("observations must contain a 'sample' column")
        }

        # rownames<- only survives on a data.frame - data.table silently drops it,
        # so observations is converted above before this assignment (Bullseye()
        # requires rownames(colData) to match the sample columns below).
        rownames(observations) <- observations\$sample

        # Bullseye() expects data_cov/data_mut with rownames = site ID and one
        # column per sample (matching rownames(colData)) - GATHER_SITES output
        # instead keeps chr/start/end/cluster/strand as columns, so build the site
        # ID from those and drop them before handing the matrices to Bullseye().
        meta_cols <- c("#chr", "start", "end", "cluster", "strand")
        site_id <- paste(cov[["#chr"]], cov\$start, cov\$end, cov\$strand, sep = "_")

        data_cov <- as.data.frame(cov)[, setdiff(colnames(cov), meta_cols), drop = FALSE]
        data_mut <- as.data.frame(mut)[, setdiff(colnames(mut), meta_cols), drop = FALSE]
        rownames(data_cov) <- site_id
        rownames(data_mut) <- site_id

        result <- Bullseye(
            data_cov = data_cov,
            data_mut = data_mut,
            colData = observations,
            design = as.formula("${bullseye_glm_design}"),
            min.cov = ${bullseye_glm_min_cov},
            return.df.only = TRUE,
            link = "${bullseye_glm_link}"
        )

        fwrite(as.data.table(result), file = "${prefix}.glm.tsv", sep = "\t")
        EOF

        Rscript run_bullseye_glm.R
        """
    }
}
