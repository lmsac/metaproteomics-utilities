arrange.spectronaut.candidates = function(candidates, comparison, p_adjust = 'fdr') {
  
  if (comparison %in% candidates$`Comparison (group1/group2)`) {
    candidates = subset(candidates, `Comparison (group1/group2)` == comparison)
  }
  else {
    comparison = sub('^([^/]) / ([^/])$', '\\2 / \\1', comparison)
    if (comparison %in% candidates$`Comparison (group1/group2)`) {
      candidates = subset(
        candidates, 
        `Comparison (group1/group2)` == comparison,
        select = c('ProteinGroups', 'AVG Log2 Ratio', 'Pvalue', 'Qvalue', 'Ratio')
      )
      candidates$`AVG Log2 Ratio` = - candidates$`AVG Log2 Ratio`
      candidates$Ratio = 1 / candidates$Ratio
    }
  }
  
  data = data.frame(
    protein = sapply(strsplit(candidates$ProteinGroups, ';'), head, 1),
    fc = candidates$Ratio,
    log2fc = candidates$`AVG Log2 Ratio`,
    pvalue = candidates$Pvalue,
    stringsAsFactors = FALSE
  )
  
  if (is.na(p_adjust) || is.null(p_adjust)) {
    data$`-log10p` = -log10(data$pvalue)
  }
  else if (p_adjust == 'qvalue') {
    data$qvalue = candidates$Qvalue
    data$`-log10p` = -log10(data$qvalue)
  }
  else {
    data$adjusted.pvalue = p.adjust(data$pvalue, method = p_adjust)
    data$`-log10p` = -log10(data$adjusted.pvalue)
  }
  
  data
}


library(readr)

candidates = read_delim(
  'candidates.xls', '\t', 
  escape_double = FALSE, trim_ws = TRUE
)


differential_protein = arrange.spectronaut.candidates(
  candidates, 'P / H', p_adjust = 'bonferroni'
)

differential_protein = plot.volcano.differential_protein(
  differential_protein, p_adjust = TRUE,
  log2FC_threshold = 1, p_threshold = 0.01,
  width = 8, height = 8,
  # filename = 'volcano.svg',
  filename = 'volcano.tiff',
  # hide_text = TRUE,
  return_data = TRUE
)

write.csv(
  differential_protein,
  file = 'differential_protein.csv',
  row.names = FALSE
)
