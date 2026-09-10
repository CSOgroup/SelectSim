make_synthetic_am <- function(n_genes = 20, n_samples = 30, seed = 1) {
  set.seed(seed)
  samples <- paste0("S", seq_len(n_samples))
  genes <- paste0("G", seq_len(n_genes))
  mk_gam <- function(p) {
    matrix(rbinom(n_genes * n_samples, 1, p), nrow = n_genes,
           dimnames = list(genes, samples))
  }
  M <- list(Nonsense = mk_gam(0.2), Missense = mk_gam(0.1))
  mk_full_tmb <- function(gam) {
    csum <- colSums(gam)
    data.frame(sample = names(csum), mutation = csum, row.names = names(csum))
  }
  list(M = M, tmb = list(Nonsense = mk_full_tmb(M$Nonsense), Missense = mk_full_tmb(M$Missense)),
       samples = samples, genes = genes)
}

test_that("factor sample.class/alteration.class no longer breaks get.blocks() naming", {
  fx <- make_synthetic_am()
  sample.class.chr <- setNames(rep(c("A", "B"), length.out = length(fx$samples)), fx$samples)
  sample.class.factor <- factor(sample.class.chr)
  alteration.class <- setNames(rep("MUT", length(fx$genes)), fx$genes)

  al_chr <- new.AL.general(fx, feat.covariates = alteration.class,
                            sample.covariates = sample.class.chr, min.freq = 1)
  al_fac <- new.AL.general(fx, feat.covariates = alteration.class,
                            sample.covariates = sample.class.factor, min.freq = 1)

  blocks_chr <- get.blocks(al_chr)
  blocks_fac <- get.blocks(al_fac)

  # Previously: factor input made which() drop names, leaving blocks unnamed
  # (length(names(...)) == 0), which cascaded into empty template matrices.
  expect_true(all(vapply(blocks_fac$sample.blocks, function(b) length(names(b)) > 0, logical(1))))
  expect_equal(blocks_chr$sample.blocks, blocks_fac$sample.blocks)
  expect_equal(al_chr$tmb$total, al_fac$tmb$total)
})

test_that("selectX() runs to completion with a factor sample.class (regression for #11)", {
  fx <- make_synthetic_am()
  sample.class <- factor(setNames(rep(c("BM", "PCT"), length.out = length(fx$samples)), fx$samples))
  alteration.class <- setNames(rep("MUT", length(fx$genes)), fx$genes)

  expect_no_error({
    result <- selectX(M = fx, sample.class = sample.class,
                       alteration.class = alteration.class,
                       min.freq = 1, n.permut = 10, n.cores = 1, verbose = FALSE)
  })
  expect_true(is.list(result))
})

test_that("partial (sparse) tmb tables are zero-filled by sample identity, not recycled", {
  fx <- make_synthetic_am()
  # Drop most Missense rows, keeping only a handful of samples - mirrors the
  # reporter's tmb tables, which only listed samples with >=1 mutation of that type.
  sparse_tmb <- fx$tmb$Missense[1:5, ]
  fx_sparse <- fx
  fx_sparse$tmb$Missense <- sparse_tmb

  sample.class <- setNames(rep("sample", length(fx$samples)), fx$samples)
  alteration.class <- setNames(rep("MUT", length(fx$genes)), fx$genes)

  expect_no_warning({
    al <- new.AL.general(fx_sparse, feat.covariates = alteration.class,
                          sample.covariates = sample.class, min.freq = 1)
  })

  # Samples absent from the sparse Missense table should contribute mutation = 0,
  # not a positionally-recycled value from an unrelated sample.
  dropped_samples <- setdiff(fx$samples, sparse_tmb$sample)
  expect_true(all(al$tmb$Missense[dropped_samples, "mutation"] == 0))

  # Total should equal Nonsense (full) + Missense (zero-filled for dropped
  # samples, actual value otherwise), matched by sample identity.
  missense_zero_filled <- setNames(rep(0, length(fx$samples)), fx$samples)
  missense_zero_filled[sparse_tmb$sample] <- sparse_tmb$mutation
  expected_total <- fx$tmb$Nonsense[fx$samples, "mutation"] + missense_zero_filled[fx$samples]
  expect_equal(unname(al$tmb$total[fx$samples]), unname(expected_total))
})

test_that("am$tmb with an unknown sample id errors clearly", {
  fx <- make_synthetic_am()
  bad_tmb <- fx$tmb$Nonsense
  bad_tmb$sample[1] <- "NOT_A_REAL_SAMPLE"
  rownames(bad_tmb)[1] <- "NOT_A_REAL_SAMPLE"
  fx$tmb$Nonsense <- bad_tmb

  sample.class <- setNames(rep("sample", length(fx$samples)), fx$samples)
  alteration.class <- setNames(rep("MUT", length(fx$genes)), fx$genes)

  expect_error(
    new.AL.general(fx, feat.covariates = alteration.class,
                    sample.covariates = sample.class, min.freq = 1),
    "not present in am\\$M"
  )
})
