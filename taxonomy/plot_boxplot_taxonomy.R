library(ggplot2)


plot.boxplot.taxonomy = function(taxonomy_quantity, taxonomy_list,
                                 filename = 'boxplot_taxonomy.svg',
                                 width = 14, height = 4, unit = 'cm',
                                 hide.text = FALSE,
                                 return_data = FALSE) {
  data = do.call(rbind, lapply(taxonomy_list, function(taxon) {
    row = taxonomy_quantity$name == taxon
    columns1 = grep('value1', colnames(taxonomy_quantity), value = TRUE)
    columns2 = grep('value2', colnames(taxonomy_quantity), value = TRUE)
    
    data.frame(
      taxon = taxon,
      group = c(rep('group1', length(columns1)), rep('group2', length(columns2))),
      sample = c(columns1, columns2),
      quantity = as.numeric(taxonomy_quantity[row, c(columns1, columns2)]),
      stringsAsFactors = FALSE
    )
  }))
  
  data$taxon = factor(data$taxon, levels = taxonomy_list)
  
  pl = ggplot(data, aes(x = group, y = log10(quantity), fill = group)) + 
    geom_boxplot(position = position_dodge(), alpha = 0.75) +
    geom_jitter(position = position_jitter(width = 0.2), size = 1) +
    scale_fill_manual(values = c('#339dff', '#ff3355')) +
    scale_y_continuous(name = 'log10 Quantity') + 
    facet_wrap(
      vars(taxon), nrow = 2, 
      scales = 'free',
      labeller = if (hide.text) function(x) '' else  'label_value'
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(),
      axis.line.y = element_line(), 
      axis.title.y = if (hide.text) element_blank() else element_text(color = 'black'),
      axis.ticks.y = element_line(),
      axis.text.y = element_text(color = 'black'),
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
      data
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

taxonomy_quantity = read_csv('taxonomy_quantity.csv')

taxonomy_list = c(
  'Proteobacteria',
  'Porphyromonadaceae',
  'Streptococcaceae',
  'Prevotellaceae',
  'Coriobacteriales',
  'Corynebacteriales',
  'Veillonellaceae',
  'Akkermansiaceae',
  'Piscirickettsiaceae',
  'Phyllobacteriaceae'
)

# taxonomy_list = c(
#   'Chloroflexi',
#   'Anaerolineae',
#   'Barnesiellaceae',
#   'Prolixibacteraceae',
#   'Eggerthellaceae',
#   'Rikenellaceae',
#   'Erysipelotrichaceae',
#   'Aspergillaceae'
# )

plot.boxplot.taxonomy(
  taxonomy_quantity,
  taxonomy_list,
  filename = 'boxplot_taxonomy.svg', 
  # hide.text = TRUE,
  width = 14, height = 4, unit = 'cm'
)
