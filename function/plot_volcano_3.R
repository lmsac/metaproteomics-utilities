calculate.fc_pvalue = function (protein_report, run_groups) {
  run_names = grep('PG\\.Quantity$', colnames(protein_report), value = TRUE)
  group_runs = lapply(run_groups, function(g) grep(g, run_names, value = TRUE))
  
  group_mean = t(apply(protein_report[, run_names], 1, function(x) {
    value = sapply(group_runs, function(r) mean(x[r]))
  }))
  colnames(group_mean) = paste0(
    'mean.', names(group_runs)
  )
  
  fc = apply(combn(length(group_runs), 2), 2, function(x) {
    group_mean[, x[2]] / group_mean[, x[1]]
  })
  colnames(fc) = paste0(
    'fc.',
    apply(combn(length(group_runs), 2), 2, function(x) {
      paste(names(group_runs)[rev(x)], collapse = '.')
    })
  )
  
  pvalue.KW = matrix(apply(protein_report[, run_names], 1, function(x) {
    value = unlist(lapply(group_runs, function(r) x[r]))
    group = unlist(lapply(names(group_runs), function(g) 
      rep(g, length(group_runs[[g]]))
    ))
    
    pvalueKW = kruskal.test(value, group)$p.value
  }))
  colnames(pvalue.KW) = 'pvalue.KW'
  
  pvalue.W = t(apply(protein_report[, run_names], 1, function(x) {
    value = unlist(lapply(group_runs, function(r) x[r]))
    group = unlist(lapply(names(group_runs), function(g) 
      rep(g, length(group_runs[[g]]))
    ))
    group = factor(group, levels = names(group_runs))
    
    pvalueW = pairwise.wilcox.test(
      value, group, 
      p.adjust.method = 'none'
    )$p.value
    
    pvalueW[lower.tri(pvalueW, diag = TRUE)]
  }))
  colnames(pvalue.W) = paste0(
    'pvalue.',
    apply(combn(length(group_runs), 2), 2, function(x) {
      paste(names(group_runs)[rev(x)], collapse = '.')
    })
  )
  
  pvalue.KW = cbind(
    pvalue.KW,
    adjusted.pvalue.KW = p.adjust(pvalue.KW, method = 'fdr')
  )
  
  pvalue.W = cbind(
    pvalue.W,
    apply(pvalue.W, 2, p.adjust, method = 'fdr')
  )
  colnames(pvalue.W)[(ncol(pvalue.W) / 2 + 1):ncol(pvalue.W)] = paste0(
    'adjusted.', 
    colnames(pvalue.W)[(ncol(pvalue.W) / 2 + 1):ncol(pvalue.W)]
  )
  
  differential_protein = data.frame(
    protein = protein_report$PG.ProteinAccessions,
    group_mean,
    fc,
    pvalue.W,
    pvalue.KW,
    stringsAsFactors = FALSE
  )
}


find.differential_significance = function(
  differential_protein,
  log2FC_threshold = 1,
  pvalue_threshold = 0.05,
  pvalue_threshold.KW = 0.05
) {
  differential_protein$significance = ''
  
  invisible(lapply(
    rev(grep('^fc', colnames(differential_protein), value = TRUE)),
    function (x) {
      sig = differential_protein$adjusted.pvalue.KW <= pvalue_threshold.KW &
        abs(log2(differential_protein[, x])) >= log2FC_threshold & 
        differential_protein[, sub('^fc', 'adjusted.pvalue', x)] <= pvalue_threshold
      group = strsplit(sub('^fc\\.', '', x), '\\.')[[1]][ifelse(
        differential_protein[, x] > 1, 1, 2
      )]
      
      differential_protein$significance <<- ifelse(
        sig & sapply(
          1:length(differential_protein$significance), 
          function (i) !grepl(
            paste0(group[i], '+'), 
            differential_protein$significance[i],
            fixed = TRUE
          )
        ),
        paste0(differential_protein$significance, group, '+'),
        differential_protein$significance
      )
    }
  ))
  
  differential_protein$significance[differential_protein$significance == ''] = 'NS'
  differential_protein
}


library(volcano3D)
library(ggplot2)

build.polar_object = function(
  differential_protein, protein_report,
  run_groups,
  log2FC_threshold = 1,
  pvalue_threshold = 0.05
) {
  run_names = grep('PG.Quantity$', colnames(protein_report), value = TRUE)
  group_runs = lapply(run_groups, function(g) grep(g, run_names, value = TRUE))
  
  sample_data = data.frame(
    ID = unlist(group_runs),
    contrast = unlist(lapply(
      names(group_runs), 
      function (x) rep(x, length(group_runs[[x]]))
    )),
    stringsAsFactors = FALSE
  )
  
  pvalues = differential_protein
  colnames(pvalues) = gsub(
    '\\.', '_', 
    sub(
      'adjusted\\.pvalue', 'padj',
      sub(
        '^(pvalue|adjusted\\.pvalue|fc||mean)\\.(.*)', 
        '\\2_\\1', 
        colnames(pvalues)
      )
    )
  )
  
  pvalues[, grep('_fc$', colnames(pvalues))] = 
    log2(pvalues[, grep('_fc$', colnames(pvalues))])
  colnames(pvalues) = sub('_fc$', '_log2fc', colnames(pvalues))
  
  pvalues[, paste0(c(names(group_runs)[c(length(group_runs), 1)], 'log2fc'), collapse = '_')] =
    -pvalues[, paste0(c(names(group_runs)[c(length(group_runs), 1)], 'log2fc'), collapse = '_')]
  colnames(pvalues) = sub(
    paste0(names(group_runs)[c(length(group_runs), 1)], collapse = '_'),
    paste0(names(group_runs)[c(1, length(group_runs))], collapse = '_'), 
    colnames(pvalues)
  )
  
  polar = polar_coords(
    sampledata = sample_data,
    contrast = 'contrast',
    pvalues = pvalues,
    expression = data.frame(protein_report[, sample_data$ID], check.names = FALSE),
    groups = rev(names(group_runs)),
    p_col_suffix = 'pvalue',
    padj_col_suffix = 'padj',
    fc_col_suffix = 'log2fc',
    multi_group_prefix = 'KW',
    non_sig_name = 'NS',
    significance_cutoff = pvalue_threshold,
    label_column = 'protein',
    fc_cutoff = log2FC_threshold,
    cutoff_criteria = 'padj'
  )
}

plot.polar.differential_protein = function(
  differential_protein, protein_report,
  run_groups,
  polar = NULL,
  log2FC_threshold = 1,
  pvalue_threshold = 0.05,
  color = c('green4', 'gold3', 'red', 'purple', 'blue', 'cyan'),
  filename = 'polar_differential_protein.svg',
  width = 10, height = 10, units = 'cm',
  return_data = FALSE
) {
  if (is.null(polar)) {
    polar = build.polar_object(
      differential_protein, protein_report, run_groups,
      log2FC_threshold,
      pvalue_threshold
    )
  }
  
  plot = radial_ggplot(
    polar = polar,
    colours = color,
    colour_scale = 'discrete',
    marker_size = 1.2,
    marker_alpha = 0.7,
    legend_size = 10
  ) + ggplot2::theme(legend.position = 'bottom')
  
  if (!is.null(filename)) {
    ggsave(
      filename = filename, 
      pl, 
      width = width, height = height, unit = unit,
      bg = 'transparent'
    )
    
    if (return_data) {
      data = polar
    }
  }
  else {
    if (return_data) {
      list(plot = pl, data = polar)
    }
    else {
      pl
    }
  }
}

plot.volcano_pairwise.differential_protein = function(
  differential_protein, protein_report,
  run_groups,
  polar = NULL,
  log2FC_threshold = 1,
  pvalue_threshold = 0.05,
  color = c('green4', 'gold3', 'red', 'purple', 'blue', 'cyan'),
  filename = 'volcano_differential_protein.svg',
  width = 10, height = 10, units = 'cm',
  return_data = FALSE
) {
  if (is.null(polar)) {
    polar = build.polar_object(
      differential_protein, protein_report, run_groups,
      log2FC_threshold,
      pvalue_threshold
    )
  }
  
  plots = volcano_trio(
    polar = polar,
    p_cutoff = pvalue_threshold,
    cutoff_criteria = 'padj',
    fc_cutoff = log2FC_threshold,
    colour_scheme = 'upregulated',
    colours = c('green4', 'gold3', 'red', 'purple', 'blue', 'cyan', "grey60"),
    marker_size = 1,
    marker_alpha = 0.7
  ) 
  
  plots = plots[names(plots) != 'All']
  plots = lapply(plots, function(p) p + theme(legend.position = 'none'))
  
  if (!is.null(filename)) {
    lapply(names(plots), function(name) {
      ggsave(
        filename = sub(
          '^([^.]*)(.[A-Za-z0-9]+)?', 
          paste0('\\1_', name, '\\2'), 
          filename
        ),
        plots[[name]], 
        width = width, height = height, unit = units,
        bg = 'transparent'
      )
    })
    
    if (return_data) {
      data = polar
    }
  }
  else {
    if (return_data) {
      list(plots = plots, data = polar)
    }
    else {
      plots
    }
  }
}



library(readr)

protein_report = read_csv('protein_report.pivot.csv')

run_groups = c(
  'Sp' = 'SP',
  'Su' = 'SU',
  'Au' = 'AU'
)

differential_protein = find.differential_proteins(protein_report, run_groups)
differential_protein = find.differential_significance(differential_protein)

write.csv(differential_protein, 'differential_protein.csv', row.names = FALSE)  


polar = plot.polar.differential_protein(
  differential_protein, protein_report,
  run_groups,
  return_data = TRUE
)

plot.volcano_pairwise.differential_protein(
  polar = polar
)

