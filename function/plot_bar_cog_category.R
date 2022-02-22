library(ggplot2)


COG_categories = c(
  A = 'RNA processing and modification',
  B = 'chromatin structure and dynamics',
  C = 'energy production and conversion',
  D = 'cell division and chromosome partitioning',
  E = 'amino acid transport and metabolism',
  F = 'nucleotide transport and metabolism',
  G = 'carbohydrate transport and metabolism', 
  H = 'coenzyme transport and metabolism',
  I = 'lipid transport and metabolism',
  J = 'translation, including ribosomal structure, and biogenesis',
  K = 'transcription',
  L = 'replication, recombination, and repair',
  M = 'cell wall, membrane, and envelope biogenesis', 
  N = 'cell motility',
  O = 'post-translational modification, protein turnover, and chaperones',
  P = 'inorganic ion transport and metabolism',
  Q = 'secondary metabolite biosynthesis, transport, and catabolism',
  R = 'general functional prediction only',
  S = 'function unknown',
  T = 'signal transduction mechanisms',
  U = 'intracellular trafficking, secretion, and vesicular transport',
  V = 'defense mechanisms',
  X = 'mobilome: prophages, transposons',
  Z = 'cytoskeleton'
)


plot.bar.cog_category = function(protein_annotation,
                             filename = 'bar_cog_category.svg', 
                             width = 8, height = 8, unit = 'cm',
                             hide.text = FALSE,
                             return_data = FALSE) {
  data = do.call(rbind, lapply(list(
    list('Up', protein_annotation$log2fc > 0), 
    list('Down', protein_annotation$log2fc < 0)
  ), function(x) {
    count = table(unlist(strsplit(protein_annotation$COG_category[x[[2]]], '')))
    
    data.frame(
      regulation = x[[1]],
      category = names(count),
      count = as.integer(count),
      stringsAsFactors = FALSE
    )
  }))
  
  data = subset(data, category != '-')
  
  pl = ggplot(
    data = data, 
    aes(
      x = factor(category, levels = rev(sort(unique(category)))), 
      y = ifelse(regulation == 'Up', count, -count), 
      fill = regulation
    )
  ) +
    geom_bar(stat = 'identity', alpha = 0.75) +
    geom_text(
      aes(
        label = count,
        hjust = ifelse(
          regulation == 'Up', 
          sqrt(pmax((count / max(count[regulation == 'Up']) - 0.75) / 0.25, 0)),
          1 - sqrt(pmax((count / max(count[regulation != 'Up']) - 0.75) / 0.25, 0))
        )
      ),
      color = 'black', vjust = 'center', size = 3
    ) +
    scale_fill_manual(values = c('#339dff', '#ff3355')) +
    scale_x_discrete(
      name = 'COG category',
      position = 'top',
      labels = function(x) paste0('[', x, '] ', COG_categories[x])
    ) + 
    scale_y_continuous(
      name = 'Frequency'
    ) +
    theme(
      axis.line.x = element_line(), 
      axis.line.y = element_blank(), 
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.x = if (hide.text) element_blank() else element_text(color = 'black'),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(),
      axis.text.y = if (hide.text) element_blank() else element_text(color = 'black'),
      axis.text.x = element_text(color = 'black'),
      legend.position = 'none'
    ) +
    coord_flip()
  
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


emapper_annotations = read_delim(
  'emapper.annotations.tsv', 
  '\t', escape_double = FALSE, trim_ws = TRUE, 
  skip = 4
)

differential_protein = read_csv('differential_protein.csv')


differential_protein_annotation = merge(
  differential_protein, 
  emapper_annotations, 
  by.x = 'protein', by.y = '#query'
)

write.csv(
  differential_protein_annotation, 
  'differential_protein_annotation.csv', 
  row.names = FALSE
)


plot.bar.cog_category(
  differential_protein_annotation,
  filename = 'bar_cog_category.svg', 
  # hide.text = TRUE,
  width = 8, height = 8, unit = 'cm'
)









