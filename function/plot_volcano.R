library(ggplot2)

plot.volcano = function(candidates, comparison,
                        p_adjust = 'fdr',
                        fc_threshold = 2,
                        p_threshold = 0.05,
                        color = c('#339dff', 'darkgrey', '#ff3355'),
                        filename = 'volcano.svg',
                        width = 12, height = 3.5, unit = 'cm',
                        hide.text = FALSE,
                        return_data = FALSE) {
  
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
    p_value = candidates$Pvalue,
    stringsAsFactors = FALSE
  )
  
  if (is.na(p_adjust) || is.null(p_adjust)) {
    data$`-log10p` = -log10(data$p_value)
  }
  else if (p_adjust == 'qvalue') {
    data$q_value = candidates$Qvalue
    data$`-log10p` = -log10(data$q_value)
  }
  else {
    data$adjusted_p_value = p.adjust(data$p_value, method = p_adjust)
    data$`-log10p` = -log10(data$adjusted_p_value)
  }
  
  
  pl = ggplot(data = data, aes(
    x = log2fc, y = `-log10p`, 
    color = ifelse(
      `-log10p` <= -log10(p_threshold), 'none', 
      ifelse(
        log2fc > log2(fc_threshold), 'up', 
        ifelse(log2fc < -log2(fc_threshold), 'down', 'none')
      )
    )
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
        if (is.na(p_adjust) || is.null(p_adjust)) 'P-Value' else 'Adjusted P-Value'
      )
    ) +
    geom_vline(
      xintercept = c(-log2(fc_threshold), log2(fc_threshold)), 
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
      axis.title.y = if (hide.text) element_blank() else element_text(color = 'black'),
      axis.title.x = if (hide.text) element_blank() else element_text(color = 'black'),
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

candidates = read_delim(
  'candidates.xls', '\t', 
  escape_double = FALSE, trim_ws = TRUE
)


fc_pvalue = plot.volcano(
  candidates, 'P / H', p_adjust = 'bonferroni',
  fc_threshold = 2, p_threshold = 0.05,
  width = 8, height = 8,
  filename = 'volcano.tiff',
  # hide.text = TRUE,
  return_data = TRUE
)

write.csv(
  fc_pvalue,
  file = 'volcano.csv',
  row.names = FALSE
)


differential_protein = subset(fc_pvalue, abs(log2fc) > 1 & `-log10p` > -log10(0.05))
  
write.csv(
  differential_protein,
  file = 'differential_protein.csv',
  row.names = FALSE
)
