calculate.fc_pvalue = function (protein_report, run_groups, method = 'wilcox', p_adjust = 'fdr', ...) {
  if (is.function(method)) {
    test_func = method
  }
  else if (method == 'wilcox') {
    test_func = wilcox.test
  }
  else if (method == 't') {
    test_func = t.test
  }
  else {
    stop('invalid test method')
  }
  
  run_names = grep('PG\\.Quantity$', colnames(protein_report), value = TRUE)
  group_runs = lapply(run_groups, function(g) grep(g, run_names, value = TRUE))
  
  lapply(names(group_runs), function(x) message(paste0(x, ': ', head(group_runs[[x]], 1), ' ...', length(group_runs[[x]]), ' runs')))
  
  group_mean = t(apply(protein_report[, run_names], 1, function(x) {
    value = sapply(group_runs, function(r) mean(x[r]))
  }))
  colnames(group_mean) = paste0(
    'mean.', names(group_runs)
  )
  
  fc = group_mean[, 2] / group_mean[, 1]
  
  pvalue = apply(protein_report[, run_names], 1, function(x) {
    value1 = log10(x[group_runs[[1]]])
    value2 = log10(x[group_runs[[2]]])
    
    pvalue = test_func(
      value1, value2, ... = ...
    )$p.value
  })
  
  differential_protein = data.frame(
    protein = protein_report$PG.ProteinAccessions,
    group_mean,
    fc = fc,
    pvalue = pvalue,
    
    stringsAsFactors = FALSE
  )
  
  if (!is.na(p_adjust) && !is.null(p_adjust)) {
    adjusted_pvalue = p.adjust(pvalue, method = p_adjust)
    differential_protein$adjusted.pvalue = adjusted_pvalue
  }
  
  differential_protein
}


library(ggplot2)

plot.volcano.differential_protein = function(differential_protein,
                        log2FC_threshold = 1,
                        p_threshold = 0.05,
                        p_adjust = TRUE,
                        color = c('down' = '#339dff', 'none' = 'darkgrey', 'up' = '#ff3355'),
                        filename = 'volcano_differential_protein.svg',
                        width = 12, height = 3.5, unit = 'cm',
                        hide_text = FALSE,
                        return_data = FALSE) {
  
  data = differential_protein
  data$log2fc = log2(data$fc)
  if (p_adjust) {
    data$`-log10p` = -log10(data$adjusted.pvalue)
  } else {
    data$`-log10p` = -log10(data$pvalue)
  }
  data$regulation = ifelse(
    data$`-log10p` <= -log10(p_threshold), 'none', 
    ifelse(
      data$log2fc > log2FC_threshold, 'up', 
      ifelse(data$log2fc < -log2FC_threshold, 'down', 'none')
    )
  )
  
  pl = ggplot(data = data, aes(
    x = log2fc, y = `-log10p`, 
    color = regulation
  )) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_manual(values = color) +
    scale_x_continuous(
      name = 'Log2 Fold Change',
      limits = c(
        -max(abs(data$log2fc)), 
        max(abs(data$log2fc))
      )
    ) +
    scale_y_continuous(
      name = paste(
        '-Log10', 
        if (!p_adjust) 'P-Value' else 'Adjusted P-Value'
      )
    ) +
    geom_vline(
      xintercept = c(-log2FC_threshold, log2FC_threshold), 
      lty = 4, col = 'black'
    ) +
    geom_hline(
      yintercept = -log10(p_threshold), 
      lty = 4, col = 'black'
    ) +
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_line(), 
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.y = if (hide_text) element_blank() else element_text(color = 'black'),
      axis.title.x = if (hide_text) element_blank() else element_text(color = 'black'),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      axis.text.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      strip.text = element_text(face = 'italic'),
      legend.position = 'none'
    )
  
  
  if (!is.null(filename)) {
    ggsave(
      filename = filename, 
      pl, 
      width = width, height = height, unit = unit,
      bg = 'transparent'
    )
    
    if (return_data) {
      data = data
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = data)
    }
    else {
      pl
    }
  }
}



library(readr)

protein_report = read_csv('protein_report.csv')

run_groups = c(
  'H' = 'H[0-9]+',
  'F' = 'F[0-9]+'
)

differential_protein = calculate.fc_pvalue(
  protein_report, run_groups,
  method = 'wilcox', p_adjust = 'fdr'
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

write.csv(differential_protein, 'differential_protein.csv', row.names = FALSE) 
