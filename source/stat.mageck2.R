stat.mageck2 <- function (untreated.list, treated.list, namecolumn = 1, fullmatchcolumn = 2, 
          norm.fun = "median", extractpattern = expression("^(.+?)_.+"), 
          mageckfolder = NULL, sort.criteria = "neg", adjust.method = "fdr", 
          filename = NULL, fdr.pval = 0.05) 
{
  non = lapply(c(untreated.list, treated.list), function(x) stopifnot(identical(x[, 
                                                                                  namecolumn], treated.list[[1]][, namecolumn])))
  gene.names = sub(extractpattern, "\\1", treated.list[[1]][, 
                                                            namecolumn], perl = TRUE)
  design.names = treated.list[[1]][, namecolumn]
  untreated.list <- do.call("cbind", lapply(untreated.list, 
                                            FUN = function(x) return(x[, fullmatchcolumn])))
  treated.list <- do.call("cbind", lapply(treated.list, FUN = function(x) return(x[, 
                                                                                   fullmatchcolumn])))
  ncol.untreated = ncol(untreated.list)
  ncol.treated = ncol(treated.list)
  designs = as.character(design.names)
  dataset.combined = data.frame(designs = designs, genes = gene.names, 
                                stringsAsFactors = FALSE)
  for (i in 1:ncol.treated) {
    dataset.combined[, 2 + i] = as.numeric(treated.list[, 
                                                        i])
  }
  for (i in 1:ncol.untreated) {
    dataset.combined[, 2 + ncol.treated + i] = as.numeric(untreated.list[, 
                                                                         i])
  }
  if (is.null(filename)) {
    filename = "mageckanalysisfile"
  }
  if (is.null(mageckfolder)) {
    dirstore = getwd()
  }
  else {
    dirstore = mageckfolder
  }
  dataset.combined.file = paste(dirstore, "/", filename, "_MAGeCK_sgRNA.tab", 
                                sep = "")
  write.table(dataset.combined, file = dataset.combined.file, 
              row.names = FALSE, quote = FALSE, sep="\t")
  treated.samples = 0
  for (i in seq(from = 1, to = (0 + ncol.treated - 1), by = 1)) {
    treated.samples = paste(treated.samples, i, sep = ",")
  }
  untreated.samples = 0 + ncol.treated
  for (i in seq(from = (0 + ncol.treated + 1), to = (0 + ncol.treated + 
                                                     ncol.untreated - 1), by = 1)) {
    untreated.samples = paste(untreated.samples, i, sep = ",")
  }
  mageckstring = paste("mageck test", "--count-table", dataset.combined.file, 
                       "--treatment-id", treated.samples, "--control-id", untreated.samples, 
                       "--norm-method", norm.fun, "--sort-criteria", sort.criteria, 
                       "--gene-test-fdr-threshold ", fdr.pval, "--adjust-method", 
                       adjust.method, "--output-prefix", filename, sep = " ")
  system(mageckstring)
  data.mageck.genes = load.file(paste(filename, "gene_summary.txt", 
                                      sep = "."))
  data.mageck.sgrna = load.file(paste(filename, "sgrna_summary.txt", 
                                      sep = "."))
  dataset.return = data.frame(row.names = data.mageck.genes$id, 
                              genes = data.mageck.genes$id, pos = data.mageck.genes[, 
                                                                                    paste("pos.", adjust.method, sep = "")], rank.pos = data.mageck.genes[, 
                                                                                                                                                          "pos.rank"], neg = data.mageck.genes[, paste("neg.", 
                                                                                                                                                                                                       adjust.method, sep = "")], rank.neg = as.numeric(data.mageck.genes[, 
                                                                                                                                                                                                                                                                          "neg.rank"]), sgrna.neg.good = as.numeric(data.mageck.genes[, 
                                                                                                                                                                                                                                                                                                                                      "neg.goodsgrna"]), sgrna.pos.good = as.numeric(data.mageck.genes[, 
                                                                                                                                                                                                                                                                                                                                                                                                       "pos.goodsgrna"]), stringsAsFactors = FALSE)
  return(list(genes = dataset.return, sgRNA = data.mageck.sgrna))
}
