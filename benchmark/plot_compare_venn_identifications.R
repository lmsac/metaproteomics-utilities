library(VennDiagram)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

plot.compare.venn.identification = function(reports, id_column, filename,
                                            run_groups = NULL, 
                                            id_leading_entry = TRUE,
                                            min_frequency = 0.6,
                                            fill = c('#339dff', '#ff3355', '#65c3ba'), 
                                            height = 3, width = 3, unit = 'cm', margin = 0.15,
                                            ..., return_data = FALSE) {
  ident_names = lapply(reports, function(report) {
    quant_col = grep('^quantity\\.', colnames(report), value = TRUE)
    
    ident_names = lapply(quant_col, function(col) {
      ident_names = unique(report[[id_column]][!is.na(report[[col]])])
      if (id_leading_entry) {
        ident_names = unique(sapply(strsplit(ident_names, ';'), head, 1))
      }
      ident_names
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
    
    ident_names = lapply(1:length(run_names), function(j) {
      ident_names = ident_names[run_names[[j]]]
      
      run_counts = table(unlist(ident_names))
      names(run_counts)[run_counts >= length(ident_names) * min_frequency]
    })
    
    ident_names = Reduce(intersect, ident_names)
  })
  
  pl = VennDiagram::venn.diagram(
    ident_names,
    fill = array(fill, dim = length(ident_names)), 
    alpha = rep(0.5, length(ident_names)),
    fontfamily = "sans", cat.fontfamily = "sans",
    col = 'white', imagetype = 'svg',
    filename = filename, 
    height = height, width = width, units = unit,
    ext.text = FALSE,
    ...
  )
  
  if (!is.null(filename)) {
    if (return_data) {
      ident_names
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = ident_names)
    }
    else {
      pl
    }
  }
}