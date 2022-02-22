library(ggplot2)

calculate.ratio = function(report, run_groups) {
  run_names = grep('^quantity', colnames(report), value = TRUE)
  
  result = report[, !colnames(report) %in% run_names]
  sapply(1:length(run_groups), function(i) {
    index = grep(run_groups[[i]], run_names)
    
    avg = apply(report[, run_names[index]], 1, function(x) {
      mean(x, na.rm = TRUE)
    })
    
    result[[names(run_groups)[i]]] <<- avg
  })
  
  result = subset(result, rowSums(
    is.na(result[names(run_groups)]) | 
      (result[names(run_groups)] == 0)
  ) == 0)
  
  
  apply(combn(length(run_groups), 2), 2, function(idx) {
    result[[paste0(names(run_groups)[idx[2]], '/', names(run_groups)[idx[1]])]] <<- 
      result[[names(run_groups)[idx[2]]]] / result[[names(run_groups)[idx[1]]]]
  })
  
  result
}


plot.fold.change = function(reports, run_groups,
                            theoretical_fold_change,
                            invert_if_reduce = FALSE,
                            ylim = if (invert_if_reduce) c(0.25, 25) else c(0.04, 25),
                            breaks = if (invert_if_reduce) c(0.5, 1, 2, 4, 8, 16) 
                                     else c(0.1, 0.2, 0.5, 1, 2, 5, 10),
                            fill = c('#339dff', '#ff3355', '#65c3ba'),
                            color = c('#194e7f', '#7f192a', '#3c756f'),
                            nrow = 2, ncol = NULL, box_width = 0.12, digits = 2,
                            filename = 'fc_organism.svg',
                            width = 12, height = 8, unit = 'cm',
                            hide.text = FALSE,
                            return_data = FALSE) {
  fold_change = lapply(reports, function(report) {
    ratio = calculate.ratio(report, run_groups)
    ratio$organism = sub('^(\\S)\\S+\\s(\\S+).*', '\\1. \\2', ratio$organism)
    ratio
  })
  
  data = do.call(rbind, lapply(names(fold_change), function(name) {
    data = reshape2::melt(
      subset(
        fold_change[[name]], 
        select = setdiff(colnames(fold_change[[name]]), names(run_groups))
      ),
      measure.vars = grep('/', colnames(fold_change[[name]]), value = TRUE),
      variable.name = 'group', 
      value.name = 'fc'
    )
    data$dataset = factor(name, levels = names(fold_change))
    data
  }))
  
  if (is.data.frame(theoretical_fold_change)) {
    theoretical_fold_change = list(theoretical_fold_change)
  }
  
  if (is.list(theoretical_fold_change) && length(theoretical_fold_change) > 0) {
    plots = lapply(1:length(theoretical_fold_change), function(i) {
      theoretical_data = theoretical_fold_change[[i]]
      
      data = subset(
        data, 
        paste(organism, group) %in% 
          paste(theoretical_data$organism, theoretical_data$group)
      )
      
      if (invert_if_reduce) {
        data$fc = ifelse(
          data$organism %in% theoretical_data$organism[theoretical_data$fc < 1],
          1 / data$fc,
          data$fc
        )
        
        theoretical_data$fc = with(theoretical_data, ifelse(fc < 1, 1 / fc, fc))
      }
      
      pl = ggplot(data, aes(x = dataset, y = fc, fill = dataset, color = dataset)) + 
        geom_violin(alpha = 0.5) +
        stat_summary(
          fun.data = function(x) {
            r = quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
            names(r) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
            r
          },
          geom = 'boxplot',
          fill = 'white',
          mapping = aes(width = box_width)
        ) +
        geom_text(
          data = aggregate(fc ~  dataset + organism, data, median), 
          aes(
            label = format(fc, digits = digits), 
            y = if (invert_if_reduce) min(fc) else fc, 
            hjust = if (invert_if_reduce) 1.05 else ifelse(fc < median(fc), -0.33, 1.33)
          ),
          angle = 90, 
          vjust = if (invert_if_reduce) 1.05 else 1.25
          , color = 'black'
        ) +
        geom_hline(
          data = theoretical_data,
          mapping = aes(yintercept = fc),
          linetype = 'dashed'
        ) +
        scale_y_continuous(
          name = 'Fold Change',
          trans = 'log2',
          breaks = breaks
        ) +
        coord_cartesian(ylim = ylim) +
        facet_wrap(
          vars(organism), 
          nrow = nrow, ncol = ncol, 
          labeller = if (hide.text) function(x) '' else  'label_value'
        ) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = color) +
        theme(
          axis.line.y = element_line(), 
          axis.line.x = element_blank(), 
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = if (hide.text) element_blank() else element_text(color = 'black'),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = 'black'),
          strip.text = element_text(face = 'italic'),
          legend.position = 'none'
        )
      
      
      if (!is.null(filename)) {
        ggsave(
          filename = sub('^([^.]*)(.[A-Za-z0-9]+)?', paste0('\\1_', i, '\\2'), filename), 
          pl, 
          width = width, height = height, unit = unit,
          bg = 'transparent'
        )
      }
      else {
        pl
      }
      
    })
    
    if (!is.null(filename)) {
      if (return_data) {
         data = fold_change
      }
    }
    else {
      if (return_data) {
        list(plots = plots, data = fold_change)
      }
      else if (length(plots) == 1) {
        plots[1]
      }
      else {
        plots
      }
    }
  }
}

