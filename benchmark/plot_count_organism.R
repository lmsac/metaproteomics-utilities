plot.count.organism = function(reports,
                               fill = c('#339dff', '#ff3355', '#65c3ba'),
                               color = c('#194e7f', '#7f192a', '#3c756f'),
                               nrow = 2, ncol = NULL,
                               filename = 'bar_organism.svg', 
                               width = 12, height = 8, unit = 'cm',
                               hide.text = FALSE,
                               return_data = FALSE) {
  count_organism = lapply(reports, function(report) { 
    count = table(sub('^(\\S)\\S+\\s(\\S+).*', '\\1. \\2', report$organism))
    data.frame(
      organism = names(count),
      number = as.integer(count),
      stringsAsFactors = FALSE
    )
  })
  
  data = do.call(rbind, lapply(names(count_organism), function(name) {
    cbind(
      dataset = name,
      count_organism[[name]]
    )
  }))
  
  pl = ggplot(
    data, 
    aes(
      x = if (length(reports) == 1) organism else dataset, 
      y = number, 
      color = dataset, fill = dataset
    )
  ) +
    geom_bar(
      alpha = 0.25,
      stat = 'identity'
    ) +
    geom_text(
      aes(
        label = number,
        hjust = ifelse(number > max(number) * 2/3, 1, 0)
      ),
      color = 'black', vjust = 'center',
      angle = 90
    ) +
    
    scale_fill_manual(values = fill) +
    scale_colour_manual(values = color) +
    
    scale_y_continuous(
      name = '# Identifications'
    ) 
  
  if (length(reports) > 1) {
    pl = pl + 
      facet_wrap(
        vars(organism), 
        nrow = nrow, ncol = ncol, 
        labeller = if (hide.text) function(x) '' else  'label_value'
      )
  }
  
  pl = pl +
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_line(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = if (hide.text) element_blank() else element_text(color = 'black'),
      axis.ticks.y = element_line(),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = if (hide.text || length(reports) > 1) element_blank() 
                    else element_text(color = 'black', angle = 90, hjust = 1),
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
      count_organism
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = count_organism)
    }
    else {
      pl
    }
  }
}