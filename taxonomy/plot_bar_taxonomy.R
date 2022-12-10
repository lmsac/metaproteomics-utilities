library(ggplot2)


plot.barplot.taxonomy = function(
  taxonomy_quantity, 
  run_groups,
  min_percentage = 0.05,
  levels = c('phylum', 'class', 'order', 'family', 'genus'),
  filename = 'barplot_taxonomy_abundance.svg',
  width = 10, height = 8, unit = 'cm',
  colors = c(
    "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
    "#FFFF33", "#A65628","#F781BF", "#1B9E77", "#D95F02","#7570B3", 
    "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#66C2A5", 
    "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", 
    "#B3B3B3", "#1A1A1A", "#67001F"
  ),
  return_data = FALSE
) {
  data = cbind(
    taxonomy_quantity[, c('level', 'name')],
    sapply(run_groups, function(group) {
      rowMeans(taxonomy_quantity[, grep(group, colnames(taxonomy_quantity))])
    })
  )
  
  data = do.call(rbind, lapply(
    setdiff(unique(data$level), 'root'), 
    function(level) {
      data_level = data[data$level == level & data$name != '-', ]
      
      data_level[, -1:-2] = t(apply(
        data_level[, -1:-2], 1, 
        function(x) x / colSums(data[match('root', data$level), -1:-2])
      ))
      
      if (!is.na(min_percentage) && !is.null(min_percentage)) {
        other_row = apply(data_level[, -1:-2], 1, min) < min_percentage * 0.01
        if (sum(other_row) > 0) {
          data_other = data_level[1, ]
          data_other$name = paste0('Other(<', min_percentage, '%)')
          data_other[, -1:-2] = colSums(data_level[other_row, -1:-2])
          data_level = rbind(
            data_level[!other_row, ],
            data_other
          )
        }
      }
      
      data_level
    }
  ))

  plots = lapply(
    levels,
    function(level) {
      df = reshape2::melt(
        data[data$level == level, ], 
        id.vers = c('level', 'name'),
        variable.name = 'group',
        value.name = 'abundance'
      )
      
      df$name = factor(df$name, levels = c(
        grep('^Other', unique(df$name), value = TRUE),
        grep('^Other', unique(df$name), invert = TRUE, value = TRUE)
      ))
      
      ggplot(df) +
        geom_bar(
          aes(x = group, y = abundance * 100, fill = name),
          stat = 'identity'
        ) +
        scale_fill_manual(
          name = stringr::str_to_title(level), 
          values = rep(colors, length.out = length(levels(df$name)))
        ) +
        scale_y_continuous(name = 'Relative Abundance (%)') +
        theme(
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(color = 'black', angle = 90, vjust = 0.5),
          axis.line.y = element_line(), 
          axis.title.y = element_text(color = 'black'),
          axis.ticks.y = element_line(),
          axis.text.y = element_text(color = 'black'),
          legend.key.size = unit(0.75, 'char'),
          legend.text = element_text(size = 8)
        )
    }
  )

  # plots = cowplot::align_plots(plotlist = plots, align='hv')
  # plots = lapply(plots, cowplot::ggdraw)
  names(plots) = levels
  
  if (!is.null(filename)) {
    lapply(names(plots), function(level) {
      ggsave(
        filename = sub(
          '^([^.]*)(.[A-Za-z0-9]+)?', 
          paste0('\\1_', level, '\\2'), 
          filename
        ), 
        plots[[level]], 
        width = width, height = height, unit = unit,
        bg = 'transparent'
      )
    })
      
    if (return_data) {
      data
    }
  }
  else {
    if (return_data) {
      list(plots = plots, data = data)
    }
    else {
      plots
    }
  }
}


run_groups = c(
  'W Sp' = '^value1',
  'Y Sp' = '^value2',
  'B Sp' = '^value3',
  'W Su' = '^value4',
  'Y Su' = '^value5',
  'B Su' = '^value6',
  'W Au' = '^value7',
  'Y Au' = '^value8',
  'B Au' = '^value9'
)

plots = plot.barplot.taxonomy(
  taxonomy_quantity,
  run_groups = run_groups,
  filename = NULL
)

plots = cowplot::align_plots(plotlist = plots, align='hv')
plots = lapply(plots, cowplot::ggdraw)

lapply(names(plots), function(level) {
  ggsave(
    filename = sub(
      '^([^.]*)(.[A-Za-z0-9]+)?', 
      paste0('\\1_', level, '\\2'), 
      'barplot_taxonomy_abundance_quantity.svg'
    ), 
    plots[[level]], 
    width = 9.5, height = 8, unit = 'cm',
    bg = 'transparent'
  )
})


plots = plot.barplot.taxonomy(
  taxonomy_quantity_pepcount_all,
  run_groups = run_groups,
  filename = NULL
)

plots = cowplot::align_plots(plotlist = plots, align='hv')
plots = lapply(plots, cowplot::ggdraw)

lapply(names(plots), function(level) {
  ggsave(
    filename = sub(
      '^([^.]*)(.[A-Za-z0-9]+)?', 
      paste0('\\1_', level, '\\2'), 
      'barplot_taxonomy_abundance_pepcount.svg'
    ), 
    plots[[level]], 
    width = 9.5, height = 8, unit = 'cm',
    bg = 'transparent'
  )
})

plots = plot.barplot.taxonomy(
  taxonomy_abundance_metegenomics,
  run_groups = c(
    'W Sp' = 'Spring_W',
    'Y Sp' = 'Spring_Y',
    'B Sp' = 'Spring_B',
    'W Su' = 'Summer_W',
    'Y Su' = 'Summer_Y',
    'B Su' = 'Summer_B',
    'W Au' = 'Autumn_W',
    'Y Au' = 'Autumn_Y',
    'B Au' = 'Autumn_B'
  ),
  filename = NULL,
  min_percentage = NA
)

plots[1:3] = cowplot::align_plots(plotlist = plots[1:3], align='hv')
plots[4:5] = cowplot::align_plots(plotlist = plots[4:5], align='hv')
plots = lapply(plots, cowplot::ggdraw)

lapply(names(plots)[1:3], function(level) {
  ggsave(
    filename = sub(
      '^([^.]*)(.[A-Za-z0-9]+)?', 
      paste0('\\1_', level, '\\2'), 
      'barplot_taxonomy_abundance_metegenomics.svg'
    ), 
    plots[[level]], 
    width = 9, height = 8, unit = 'cm',
    bg = 'transparent'
  )
})
lapply(names(plots)[4:5], function(level) {
  ggsave(
    filename = sub(
      '^([^.]*)(.[A-Za-z0-9]+)?', 
      paste0('\\1_', level, '\\2'), 
      'barplot_taxonomy_abundance_metegenomics.svg'
    ), 
    plots[[level]], 
    width = 12.5, height = 8, unit = 'cm',
    bg = 'transparent'
  )
})
