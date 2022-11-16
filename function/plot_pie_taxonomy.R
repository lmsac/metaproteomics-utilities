library(ggplot2)


plot.pie.taxonomy = function(protein_annotation,
                             level,
                             parent_taxonomy = NULL, topn = 5,
                             filename = 'pie_taxonomy.svg', 
                             width = 6, height = 6, unit = 'cm',
                             hide_text = FALSE,
                             return_data = FALSE) {
  if (is.null(parent_taxonomy) || is.na(parent_taxonomy[1])) {
    eggNOG_OGs = protein_annotation$eggNOG_OGs
  }
  else {
    eggNOG_OGs = grep(
      parent_taxonomy, 
      protein_annotation$eggNOG_OGs, 
      value = TRUE
    )
  }
  
  taxonomy_count = sort(table(
    sapply(strsplit(eggNOG_OGs, ','), function(s) {
      s = unique(sub('^[^|]+\\|(.*)', '\\1', s))
      if (length(s) < level) {
        'NA'
      }
      else {
        s[level]
      }
    })
  ), decreasing = TRUE)
  taxonomy_count = c(taxonomy_count[names(taxonomy_count) != 'NA'], taxonomy_count['NA'])
  
  if (length(taxonomy_count) > topn) {
    taxonomy_count = c(
      head(taxonomy_count, topn),
      'Other' = sum(taxonomy_count[-(1:topn)])
    )
  }
  names(taxonomy_count)[names(taxonomy_count) == 'NA'] = 'Other'
  
  data = data.frame(
    taxonomy = names(taxonomy_count),
    count = as.integer(taxonomy_count),
    stringsAsFactors = FALSE
  )
  
  pl = ggplot(
    data = data, 
    mapping = aes(x = 1, y = count, fill = taxonomy)
  ) + 
    geom_bar(stat = 'identity', position = 'stack', width = 1, alpha = 0.33) +
    coord_polar(theta = 'y') +
    scale_fill_brewer(
      name = 'Taxonomy',
      labels = local({
        labels = paste0(
          data$taxonomy, ' (', 
          round(data$count / sum(data$count) * 100, 1), '%)'
        )
        names(labels) = data$taxonomy
        labels
      }),
      palette = 'Set1'
    ) +
    theme(
      axis.line.x = element_blank(), 
      axis.line.y = element_blank(), 
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.position = if (hide_text) 'none' else 'right'
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


differential_protein_annotation = read_csv('differential_protein_annotation.csv')


taxonomy_count = plot.pie.taxonomy(
  differential_protein_annotation,
  level = 3, topn = 3,
  filename = 'pie_taxonomy.svg',
  # hide_text = TRUE,
  width = 6, height = 6, unit = 'cm',
  return_data = TRUE
)


lapply(1:2, function(i) {
  plot.pie.taxonomy(
    differential_protein_annotation,
    level = 4, topn = 4,
    parent_taxonomy = taxonomy_count$taxonomy[i],
    filename = paste0('pie_taxonomy_', taxonomy_count$taxonomy[i], '.svg'),
    # hide_text = TRUE,
    width = 6, height = 6, unit = 'cm'
  )
})

