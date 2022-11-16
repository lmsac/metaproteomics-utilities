library(ggplot2)

calculate.ratio = function(report, run_groups, min_frequency = 0.6) {
  run_names = grep('^quantity', colnames(report), value = TRUE)
  
  result = report[, !colnames(report) %in% run_names]
  sapply(1:length(run_groups), function(i) {
    index = grep(run_groups[[i]], sub('^quantity\\.', '', run_names))
    
    avg = apply(report[, run_names[index]], 1, function(x) {
      if (sum(!is.na(x)) >= length(x) * min_frequency) {
        mean(x, na.rm = TRUE)
      } else {
        NA
      }
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


plot.compare.fold.change = function(reports, run_groups,
                            theoretical_fold_change, ...,
                            invert_if_reduce = FALSE,
                            show_error = TRUE, mark_best = TRUE,
                            ylim = if (invert_if_reduce) c(0.25, 25) else c(0.04, 25),
                            breaks = if (invert_if_reduce) c(0.5, 1, 2, 4, 8, 16) else c(0.1, 0.2, 0.5, 1, 2, 5, 10),
                            fill = c('#339dff', '#ff3355', '#65c3ba'),
                            color = c('#194e7f', '#7f192a', '#3c756f'),
                            nrow = 2, ncol = NULL,
                            filename = 'foldchange_organism.svg',
                            width = 15, height = 10, unit = 'cm',
                            hide_text = FALSE,
                            return_data = FALSE) {
  fold_change = lapply(reports, function(report) {
    calculate.ratio(report, run_groups, ... = ...)
  })
  
  data = do.call(rbind, lapply(names(fold_change), function(name) {
    data = reshape2::melt(
      subset(
        fold_change[[name]], 
        select = setdiff(colnames(fold_change[[name]]), names(run_groups))
      ),
      measure.vars = grep('/', colnames(fold_change[[name]]), value = TRUE),
      variable.name = 'group', 
      value.name = 'fc',
      na.rm = TRUE
    )
    data$dataset = factor(name, levels = names(fold_change))
    data
  }))
  
  if (is.data.frame(theoretical_fold_change)) {
    theoretical_fold_change = list(theoretical_fold_change)
  }
  
  if (is.list(theoretical_fold_change) && length(theoretical_fold_change) > 0) {
    plots_data = lapply(1:length(theoretical_fold_change), function(i) {
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
      
      data_median = cbind(
        aggregate(fc ~  dataset + organism, data, median),
        count = aggregate(fc ~  dataset + organism, data, length)$fc
      )
      if (show_error || mark_best) {
        data_median$error = data_median$fc - theoretical_data$fc[match(
          data_median$organism,
          theoretical_data$organism
        )]
        
        if (mark_best) {
          data_median$rank = 0
          lapply(unique(data_median$organism), function(org) {
            data_median$rank[data_median$organism == org] <<-
              rank(abs(data_median$error[data_median$organism == org]))
          })
        }
      }
      
      pl = ggplot(data, aes(x = dataset, y = fc, fill = dataset, color = dataset)) + 
        geom_violin(alpha = 0.5, size = 0) +
        stat_summary(
          fun.data = function(x) {
            r = quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
            names(r) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
            r
          },
          geom = 'boxplot',
          fill = 'white',
          width = 0.15
        ) +
        geom_text(
          data = data_median, 
          aes(
            label = if (show_error) sprintf('%+.2f', error) else sprintf('%.2f', fc), 
            y = max(fc), 
            hjust = -0.05
          ),
          angle = 90, size = 3,
          vjust = 1.5, 
          color = if (mark_best) {
            ifelse(data_median$rank == 1, 'black', 'darkgrey') 
          } else {
            'black'
          }
        ) +
        geom_text(
          data = data_median, 
          aes(
            label = sprintf('italic(n)==%d', count), 
            y = min(fc), 
            hjust = 1.25
          ),
          angle = 90, size = 2,
          vjust = 1.75, 
          color = 'black', parse = TRUE
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
          labeller = if (hide_text) function(x) '' else  'label_value'
        ) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = color) +
        theme(
          axis.line.y = element_line(), 
          axis.line.x = element_blank(), 
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y = if (hide_text) element_blank() else element_text(color = 'black'),
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
        data_median
      }
      else {
        list(pl, data_median)
      }
    })
    
    
    if (!is.null(filename)) {
      if (return_data) {
        data = list(data = fold_change, summary = plots_data)
      }
    }
    else {
      plots = lapply(plots_data, head, 1)
      if (return_data) {
        list(plots = plots, data = fold_change, summary = lapply(plots_data, tail, 1))
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

