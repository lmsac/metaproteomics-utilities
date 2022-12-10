library(ggplot2)


plot.boxplot.taxonomy.alpha = function(
  taxonomy_quantity, 
  run_groups,
  levels = c('phylum', 'class', 'order', 'family', 'genus'),
  filename = 'boxplot_taxonomy_alpha.svg',
  width = 10, height = 8, unit = 'cm',
  colors = c('#339dff', '#ff3355'),
  return_data = FALSE
) {
  group_run = lapply(run_groups, function(group) {
    grep(group, colnames(taxonomy_quantity), value = TRUE)
  })
  
  data = taxonomy_quantity[, c('level', 'name', unlist(group_run))]
  
  data = do.call(rbind, lapply(
    setdiff(unique(data$level), 'root'), 
    function(level) {
      data_level = data[data$level == level & data$name != '-', ]
      
      data_level[, -1:-2] = t(apply(
        data_level[, -1:-2], 1, 
        function(x) x / colSums(data[match('root', data$level), -1:-2])
      ))
      
      other = 1 - colSums(data_level[, -1:-2])
      
      shannon = rowSums(apply(rbind(data_level[, -1:-2], other), 1, function(x) {
        ifelse(x == 0, 0, -x * log(x))
      }))
      
      data.frame(
        level = level,
        group = unlist(lapply(
          names(group_run), 
          function(x) rep(x, length(group_run[[x]]))
        )),
        sample = unlist(group_run),
        shannon = as.numeric(shannon),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  plots = lapply(
    levels,
    function(level) {
      df = data[data$level == level, ]
      df$group = factor(df$group, levels = unique(df$group))
      
      ggplot(df, aes(x = group, y = shannon, fill = group)) + 
        geom_boxplot(position = position_dodge(), alpha = 0.75) +
        geom_jitter(position = position_jitter(width = 0.2), size = 1) +
        scale_fill_manual(values = colors) +
        scale_y_continuous(name = 'Shannon Index') +
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



library(readr)

taxonomy_quantity = read_csv('taxonomy_quantity.csv')

plots = plot.boxplot.taxonomy.alpha(
  taxonomy_quantity,
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
  ),
  filename = NULL,
  colors = c(
    'W Sp' = '#339dff',
    'Y Sp' = '#ff3355',
    'B Sp' = '#65c3ba',
    'W Su' = '#339dff',
    'Y Su' = '#ff3355',
    'B Su' = '#65c3ba',
    'W Au' = '#339dff',
    'Y Au' = '#ff3355',
    'B Au' = '#65c3ba'
  )
)


plots[['genus']] = plots[['genus']] + 
  theme(legend.position = 'none') +
  stat_compare_means(
    label = 'p.format', 
    method = 'wilcox.test', 
    comparisons = list(
      c('W Sp', 'Y Sp'),
      c('W Sp', 'B Sp'),
      c('Y Sp', 'B Sp'),
      c('W Su', 'Y Su'),
      c('W Su', 'B Su'),
      c('Y Su', 'B Su'),
      c('W Au', 'Y Au'),
      c('W Au', 'B Au'),
      c('Y Au', 'B Au')
    )

    # comparisons = list(
    #   c('W Sp', 'W Su'),
    #   c('W Sp', 'W Au'),
    #   c('W Su', 'W Au'),
    #   c('Y Sp', 'Y Su'),
    #   c('Y Sp', 'Y Au'),
    #   c('Y Su', 'Y Au'),
    #   c('B Sp', 'B Su'),
    #   c('B Sp', 'B Au'),
    #   c('B Su', 'B Au')
    # )

    # comparisons = list(
    #   c('W Sp', 'Y Sp'),
    #   c('W Au', 'Y Au'),
    #   c('W Su', 'B Su'),
    #   c('W Sp', 'B Sp'),
    #   c('W Au', 'B Au'),
    #   c('B Sp', 'B Su'),
    #   c('B Sp', 'B Au')
    # )
  )

ggsave(
  filename = sub(
    'boxplot_taxonomy_alpha_genus.svg'
  ), 
  plots[['genus']], 
  width = 9.5, height = 8, unit = 'cm',
  bg = 'transparent'
)
