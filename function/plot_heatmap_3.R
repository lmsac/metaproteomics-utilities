library(ggplot2)

plot.heatmap = function(protein_report, protein_annotation,
                        log_transform = FALSE,
                        autoscale = FALSE,
                        filename = 'heatmap.svg', 
                        width = 12, height = 8, unit = 'cm',
                        hide.text = FALSE,
                        return_data = FALSE) {
  protein_quantity = protein_report[
    match(
      protein_annotation$protein, 
      sapply(strsplit(protein_report$PG.ProteinAccessions, ';'), head, 1)
    ),
    grep('Quantity', colnames(protein_report)),
    drop = FALSE
  ]
  if (log_transform) {
    protein_quantity = log10(protein_quantity)
  }
  if (autoscale) {
    protein_quantity = t(apply(
      protein_quantity, 1, 
      function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    ))
  }
  
  data = data.frame(
    protein_annotation,
    protein_quantity,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  df = reshape2::melt(
    data,
    measure.vars = grep('Quantity', colnames(data), value = TRUE),
    variable.name = 'sample',
    value.name = 'quantity'
  )
  
  df$protein = factor(df$protein, levels = rev(unique(df$protein)))
  
  pl = ggplot(df) + 
    geom_tile(aes(x = sample, y = protein, fill = quantity)) +
    scale_fill_gradientn(
      name = ifelse(
        autoscale, 'z-score', 
        ifelse(log_transform, 'log10 quantity', 'quantity')
      ),
      rescaler = function(x, ...) scales::rescale_mid(
        x, ..., 
        mid = mean(df$quantity)
      ),
      values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
      colors = c(
        '#323898', '#4474B6', '#75B0D0', '#AADAE7', '#E4F2F8',
        '#FFFEC1', 
        '#F7DE9A', '#FDAD62', '#F57246', '#D63026', '#A50024'
      ),
      limits = mean(df$quantity) + 
        c(-1, 1) * max(abs(df$quantity - mean(df$quantity)))
    ) +
    scale_x_discrete(
      name = 'Sample'
    ) +
    scale_y_discrete(
      name = 'Protein', 
      position = 'right',
      labels = function (x) {
        idx = match(x, df$protein)
        df$description[idx]
      }
    ) +
    guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = 0.5)) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = if (hide.text) element_blank() else element_text(color = 'black'),
      legend.position = 'bottom',
      legend.key.size = unit(0.4, 'cm'), 
      legend.title = if (hide.text) element_blank() else element_text(size = 7), 
      legend.text = if (hide.text) element_blank() else element_text(size = 6)
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



