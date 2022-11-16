library(ggplot2)

calculate.cv = function(report, run_groups, min_frequency = 1) {
  run_names = grep('^quantity', colnames(report), value = TRUE)
  
  result = report[, !colnames(report) %in% run_names, drop = FALSE]
  sapply(1:length(run_groups), function(i) {
    index = grep(run_groups[[i]], sub('^quantity\\.', '', run_names))
    
    cv = apply(report[, run_names[index]], 1, function(x) {
      if (sum(!is.na(x)) >= length(x) * min_frequency) {
        sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      } else {
        NA
      }
    })
    
    result[[names(run_groups)[i]]] <<- cv
  })
  
  result = subset(result, rowSums(!is.na(result[names(run_groups)])) > 0)
  
  result
}


plot.cv.distribution = function(report, run_groups, ...,
                                        fill = c('#339dff', '#ff3355', '#65c3ba'),
                                        color = c('#194e7f', '#7f192a', '#3c756f'),
                                        filename = 'cv_distribution.svg', 
                                        width = 5, height = 5, unit = 'cm',
                                        hide_text = FALSE,
                                        return_data = FALSE) {
  cv = calculate.cv(
    report,
    run_groups = run_groups, ... = ...
  )
  data = reshape2::melt(
    cv,
    measure.vars = names(run_groups),
    variable.name = 'group',
    value.name = 'cv',
    na.rm = TRUE
  )
  
  data_median = aggregate(cv ~ group, data, median)
  pl = ggplot(
    data, 
    aes(x = cv, color = group, fill = group)
  ) +
    geom_density(alpha = 0.25) +
    geom_vline(
      data = data_median,
      aes(xintercept = cv, color = group),
      linetype = 'dashed'
    ) +
    geom_text(
      data = data_median,
      aes(
        label = sprintf(ifelse(abs(cv) >= 0.01, '%.1f%%', '%.2f%%'), cv * 100), 
        color = group
      ),
      y = Inf, hjust = -0.25, vjust = 2
    ) +
    
    scale_fill_manual(values = fill) +
    scale_colour_manual(values = color) +
    scale_x_continuous(
      name = 'CV', breaks = seq(0, 1, 0.2)
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    facet_grid(rows = vars(group)) +
    theme(
      axis.line.x = element_line(),
      axis.line.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = 'black'),
      strip.background = element_blank(),
      strip.text.x = {
        if (length(unique(data$group)) == 1) element_blank()
        else element_text(color = 'black')
      },
      legend.position = 'none'
    )
  
  
  if (hide_text) {
    pl = pl + theme(
      strip.text.x = element_blank(),
      axis.title.x = element_blank()
    )
  }
  
  if (!is.null(filename)) {
    ggsave(
      filename = filename, 
      pl, 
      width = width, height = height, unit = unit,
      bg = 'transparent'
    )
    
    if (return_data) {
      list(data = cv, summary = data_median)
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = cv, summary = data_median)
    }
    else {
      pl
    }
  }
}


plot.compare.cv.distribution = function(reports, run_groups, ...,
                                fill = c('#339dff', '#ff3355', '#65c3ba'),
                                color = c('#194e7f', '#7f192a', '#3c756f'),
                                filename = 'compare_cv_distribution.svg', 
                                width = 12, height = 5, unit = 'cm',
                                hide_text = FALSE,
                                return_data = FALSE) {
  cv = lapply(reports, function(report) {
    calculate.cv(
      report,
      run_groups = run_groups, ... = ...
    )
  })
  data = do.call(rbind, lapply(names(cv), function(name) {
    data = reshape2::melt(
      cv[[name]],
      measure.vars = names(run_groups),
      variable.name = 'group',
      value.name = 'cv',
      na.rm = TRUE
    )
    data$dataset = factor(name, levels = names(cv))
    data
  }))
  
  data_median = aggregate(cv ~ dataset + group, data, median)
  pl = ggplot(
    data, 
    aes(x = cv, color = dataset, fill = dataset)
  ) +
    geom_density(alpha = 0.25) +
    geom_vline(
      data = data_median,
      aes(xintercept = cv, color = dataset),
      linetype = 'dashed'
    ) +
    geom_text(
      data = data_median,
      aes(
        label = sprintf(ifelse(abs(cv) >= 0.01, '%.1f%%', '%.2f%%'), cv * 100), 
        color = dataset
      ),
      y = Inf, hjust = -0.25, vjust = 2
    ) +
    
    scale_fill_manual(values = fill) +
    scale_colour_manual(values = color) +
    scale_x_continuous(
      name = 'CV', breaks = seq(0, 1, 0.2)
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    facet_grid(rows = vars(dataset), cols = vars(group)) +
    theme(
      axis.line.x = element_line(),
      axis.line.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = 'black'),
      strip.background = element_blank(),
      strip.text.y = {
        if (length(unique(data$dataset)) == 1) element_blank() 
        else element_text(color = 'black')
      },
      strip.text.x = {
        if (length(unique(data$group)) == 1) element_blank()
        else element_text(color = 'black')
      },
      legend.position = 'none'
    )
  
  
  if (hide_text) {
    pl = pl + theme(
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      axis.title.x = element_blank()
    )
  }
  
  if (!is.null(filename)) {
    ggsave(
      filename = filename, 
      pl, 
      width = width, height = height, unit = unit,
      bg = 'transparent'
    )
    
    if (return_data) {
      data = list(data = cv, summary = data_median)
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = cv, summary = data_median)
    }
    else {
      pl
    }
  }
}
