stat.wilcox2 <- function (untreated.list = list(NULL, NULL), treated.list = list(NULL, 
                                                                                 NULL), namecolumn = 1, fullmatchcolumn = 2, normalize = TRUE, 
                          norm.fun = median, extractpattern = expression("^(.+?)_.+"), 
                          controls = NULL, control.picks = 300, sorting = TRUE) 
{
  non = lapply(c(untreated.list, treated.list), function(x) stopifnot(identical(x[, 
                                                                                  namecolumn], treated.list[[1]][, namecolumn])))
  gene.names = sub(extractpattern, "\\1", treated.list[[1]][, 
                                                            namecolumn], perl = TRUE)
  designs = treated.list[[1]][, namecolumn]
  untreated.list = apply(do.call("cbind", lapply(untreated.list, 
                                                 function(x) {
                                                   x[, fullmatchcolumn] = apply(x, 1, function(v) {
                                                     if (as.numeric(v[fullmatchcolumn]) == 0 || is.na(v[fullmatchcolumn])) {
                                                       return(as.numeric(1))
                                                     }
                                                     else {
                                                       return(as.numeric(v[fullmatchcolumn]))
                                                     }
                                                   })
                                                   if (normalize) {
                                                     return(x[, fullmatchcolumn]/norm.fun(x[, fullmatchcolumn]) + 
                                                              1)
                                                   }
                                                   else {
                                                     return(x[, fullmatchcolumn])
                                                   }
                                                 })), 1, mean)
  treated.list = apply(do.call("cbind", lapply(treated.list, 
                                               function(x) {
                                                 x[, fullmatchcolumn] = apply(x, 1, function(v) {
                                                   if (as.numeric(v[fullmatchcolumn]) == 0 || is.na(v[fullmatchcolumn])) {
                                                     return(as.numeric(1))
                                                   }
                                                   else {
                                                     return(as.numeric(v[fullmatchcolumn]))
                                                   }
                                                 })
                                                 if (normalize) {
                                                   return(x[, fullmatchcolumn]/norm.fun(x[, fullmatchcolumn]) + 
                                                            1)
                                                 }
                                                 else {
                                                   return(x[, fullmatchcolumn])
                                                 }
                                               })), 1, mean)
  dataset.combined <- data.frame(untreated = as.numeric(untreated.list), 
                                 treated = as.numeric(treated.list), genes = gene.names, 
                                 foldchange = as.numeric(treated.list)/as.numeric(untreated.list), 
                                 row.names = designs, stringsAsFactors = FALSE)
  for (i in 1:nrow(dataset.combined)) {
    if (!is.finite(dataset.combined[i, 4])) {
      dataset.combined[i, 4] = 1
    }
    else if (dataset.combined[i, 4] == 0) {
      dataset.combined[i, 4] = 1
    }
  }
  if (!is.null(controls)) {
    control.test = dataset.combined$foldchange[dataset.combined$genes == 
                                                 controls]
  }
  else {
    random.picked = dataset.combined[sample(nrow(dataset.combined), 
                                            control.picks), "foldchange"]
    control.test = random.picked
  }
  pvals = do.call("rbind.data.frame", lapply(split(dataset.combined, 
                                                   f = dataset.combined$genes), FUN = function(x) {
                                                     c(mean(x$untreated), mean(x$treated), mean(x$foldchange), 
                                                       wilcox.test(x$foldchange, control.test, alternative = "t")$p.value)
                                                   }))
  names(pvals) = c("untreated", "treated", "foldchange", "p.value")
  pvals$p.value = p.adjust(pvals$p.value, method = "BH", n = length(pvals$p.value))
  row.names(pvals) <- row.names(pvals) <- unique(dataset.combined$genes)
  if (sorting) {
    return(pvals[order(pvals$p.value), ])
  }
  else {
    return(pvals)
  }
}