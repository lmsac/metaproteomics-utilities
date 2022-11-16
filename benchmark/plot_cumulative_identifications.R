library(ggplot2)

get.cumulative.identifications = function(report, id_column, run_groups = NULL, ...) {
  quant_col = grep('^quantity\\.', colnames(report), value = TRUE)
  
  ident_names = lapply(quant_col, function(col) {
    unique(report[[id_column]][!is.na(report[[col]])])
  })
  names(ident_names) = sub('^quantity\\.', '', quant_col)
  
  if (!is.null(run_groups)) {
    run_names = lapply(run_groups, function(run_group) {
      grep(run_group, names(ident_names), value = TRUE)
    })
  }
  else {
    run_names = list('All' = names(ident_names))
  }
  
  cumulative_numbers = do.call(rbind, lapply(1:length(run_names), function(j) {
    ident_names = ident_names[run_names[[j]]]
    
    cumulative_ident_names = Reduce(
      function(cum, ident_name) {
        list(
          sparse = union(cum$sparse, ident_name),
          full = intersect(cum$full, ident_name)
        )
      },
      ident_names,
      init = list(
        sparse = ident_names[[1]],
        full = ident_names[[1]]
      ),
      accumulate = TRUE
    )[-1]
    
    cumulative_numbers = do.call(rbind, lapply(1:length(cumulative_ident_names), function(run) {
      data.frame(
        group = names(run_names)[j],
        run = run,
        full = length(cumulative_ident_names[[run]]$full),
        sparse = length(cumulative_ident_names[[run]]$sparse),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  cumulative_numbers
}


plot.cumulative.identifications = function(report, id_column, ..., 
                                           filename = 'cumulative_identifications.svg',
                                           width = 12, height = 8, unit = 'cm',
                                           return_data = FALSE) {
  cumulative_identifications = get.cumulative.identifications(report, id_column, ...)
  
  data = reshape(
    cumulative_identifications, 
    direction = 'long', 
    varying = c('full', 'sparse'), 
    timevar = 'category', 
    times = c('full', 'sparse'), 
    v.names = 'count'
  )
  
  pl = ggplot(
    data,
    mapping = aes(x = as.factor(run), y = count)
  ) +
    geom_bar(
      aes(fill = category),
      stat = 'identity', 
      position = position_dodge(),
      alpha = 0.75, color = NA
    ) +
    geom_text(
      aes(label = count, color = category),
      position = position_dodge(width = 0.9),
      angle = 90, hjust = 'right',
      size = 3
    ) +
    scale_fill_manual(values = c(
      sparse = '#339dff', full = '#65c3ba'
    )) +
    scale_colour_manual(values = c(
      sparse = '#194e7f', full = '#3c756f'
    )) +
    scale_x_discrete(
      name = 'Cumulative Runs'
    ) +
    scale_y_continuous(
      name = paste('#', stringr::str_to_title(id_column)),
      expand = expansion(mult = c(0, 0.075))
    ) +
    theme(
      axis.line.x = element_line(), 
      axis.line.y = element_line(), 
      panel.background = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.title.y = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = element_text(color = 'black'),
      legend.position = 'none'
    )
  
  if (length(unique(data$group)) > 1) {
    pl = pl + 
      facet_grid(cols = vars(group))
  }
  
  if (!is.null(filename)) {
    ggsave(
      filename = filename, 
      pl, 
      width = width, height = height, unit = unit,
      bg = 'transparent'
    )
    if (return_data) {
      cumulative_identifications
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = cumulative_identifications)
    }
    else {
      pl
    }
  }
}