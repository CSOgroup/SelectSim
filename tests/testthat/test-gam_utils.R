test_that("maf2gam produces a binary matrix from luad_maf", {
    data(luad_maf, package = "SelectSim")
    gam <- maf2gam(luad_maf,
                   sample.col = TCGA_maf_schema$column$sample,
                   gene.col   = TCGA_maf_schema$column$gene)
    expect_true(is.matrix(gam))
    expect_true(nrow(gam) > 0)
    expect_true(ncol(gam) > 0)
    expect_true(all(gam %in% c(0L, 1L, TRUE, FALSE, NA)))
})

test_that("filter_maf_truncating retains only truncating variant types", {
    data(luad_maf, package = "SelectSim")
    filtered <- filter_maf_truncating(luad_maf, schema = TCGA_maf_schema)
    expect_true(nrow(filtered) > 0)
    expect_true(all(
        filtered[[TCGA_maf_schema$column$mutation.type]] %in%
        TCGA_maf_schema$mutation.type$truncating
    ))
})

test_that("filter_maf_missense retains only missense variant types", {
    data(luad_maf, package = "SelectSim")
    filtered <- filter_maf_missense(luad_maf, schema = TCGA_maf_schema)
    expect_true(nrow(filtered) > 0)
    expect_true(all(
        filtered[[TCGA_maf_schema$column$mutation.type]] %in%
        TCGA_maf_schema$mutation.type$missense
    ))
})

test_that("filter_maf_sample subsets to requested samples", {
    data(luad_maf, package = "SelectSim")
    sample_col <- TCGA_maf_schema$column$sample
    five_samples <- unique(luad_maf[[sample_col]])[seq_len(5)]
    filtered <- filter_maf_sample(luad_maf, samples = five_samples,
                                  sample.col = sample_col)
    expect_true(nrow(filtered) > 0)
    expect_true(all(filtered[[sample_col]] %in% five_samples))
})

test_that("stat_maf_sample returns counts summing to total mutations", {
    data(luad_maf, package = "SelectSim")
    counts <- stat_maf_sample(luad_maf, column = TCGA_maf_schema$column$sample)
    expect_equal(sum(counts), nrow(luad_maf))
})

test_that("stat_maf_gene returns a named table with all genes", {
    data(luad_maf, package = "SelectSim")
    counts <- stat_maf_gene(luad_maf, column = TCGA_maf_schema$column$gene)
    expect_true(is.table(counts))
    expected_genes <- unique(luad_maf[[TCGA_maf_schema$column$gene]])
    expect_true(all(expected_genes %in% names(counts)))
})
