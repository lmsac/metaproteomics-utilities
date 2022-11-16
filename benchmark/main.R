source('arrange_report.R')
source('plot_run_identifications.R')
source('plot_cumulative_identifications.R')
source('plot_compare_venn_identifications.R')
source('plot_compare_cv_distribution.R')
source('plot_compare_foldchange.R')

library(readr)


level = 'protein'
# DIA
report_files = list(
  'directDIA' = 'DIA/directDIA/protein_directDIA_report.csv',
  'SN DDALib' = 'DIA/Spectronaut_DDALib/protein_DDALib_report.csv',
  'DIA-NN' = 'DIA/DIA-NN/report.pg_matrix.tsv'
)
# DDA LFQ
report_files = list(
  'PEAKS' = 'DDA/PEAKS/proteins.csv',
  'MaxQuant' = 'DDA/MaxQuant/proteinGroups.txt',
  'FragPipe' = 'DDA/FragPipe/combined_protein.tsv'
)
# TMT
report_files = list(
  'PEAKS' = 'TMT/PEAKS/proteins.csv',
  'MaxQuant' = 'TMT/MaxQuant/proteinGroups.txt',
  # 'FragPipe' = 'TMT/FragPipe/tmt-report/abundance_protein_MD.tsv'
  'FragPipe' = 'TMT/FragPipe/protein.tsv'
)

report_files = list(
  'directDIA' = 'DIA/directDIA/protein_directDIA_report.csv',
  'LFQ-DDA' = 'DDA/PEAKS/proteins.csv',
  'TMT' = 'TMT/PEAKS/proteins.csv'
)


level = 'peptide'
# DIA
report_files = list(
  'directDIA' = 'DIA/directDIA/peptide_directDIA_report.csv',
  'SN DDALib' = 'DIA/Spectronaut_DDALib/peptide_DDALib_report.csv',
  'DIA-NN' = 'DIA/DIA-NN/report.pr_matrix.tsv'
)
# DDA LFQ
report_files = list(
  'PEAKS' = 'DDA/PEAKS/protein-peptides.csv',
  'MaxQuant' = 'DDA/MaxQuant/peptides.txt',
  'FragPipe' = 'DDA/FragPipe/combined_peptide.tsv'
)
# TMT
report_files = list(
  'PEAKS' = 'TMT/PEAKS/protein-peptides.csv',
  'MaxQuant' = 'TMT/MaxQuant/peptides.txt',
  #'FragPipe' = 'TMT/FragPipe/tmt-report/abundance_peptide_MD.tsv'
  'FragPipe' = 'TMT/FragPipe/peptide.tsv'
)

report_files = list(
  'directDIA' = 'DIA/directDIA/peptide_directDIA_report.csv',
  'LFQ-DDA' = 'DDA/PEAKS/protein-peptides.csv',
  'TMT' = 'TMT/PEAKS/protein-peptides.csv'
)



reports = lapply(report_files, function(f) {
  if (endsWith(f, '.csv')) {
    read_csv(f)
  } else if (endsWith(f, '.tsv') || endsWith(f, '.txt')) {
    read_delim(f, '\t', guess_max = 100000)
  } else {
    stop(paste('invalid format:', f))
  }
})

reports = lapply(reports, function(report) {
  if (check.spectronaut.report(report, level)) {
    report = arrange.spectronaut.report(report, level)
  }
  else if (check.diann.report(report, level)) {
    report = arrange.diann.report(report, level)
  }
  else if (check.peaks.report(report, level)) {
    report = arrange.peaks.report(report, level)
  }
  else if (check.maxquant.report(report, level)) {
    report = arrange.maxquant.report(report, level)
  }
  else if (check.fragpipe.report(report, level)) {
    report = arrange.fragpipe.report(report, level)
  }
  else {
    stop('unknown format')
  }
})


# assign runs per sample group
# Spectronaut|PEAKS-LFQ|MaxQuant-LFQ|PEAKS-TMT|MaxQuant-TMT|FragPipe-LFQ|FragPipe-TMT
run_groups = c(
  'S1' = 'DIA-1|Sample [1-3]|1-|(126|127N|127C)|^[1-3]( 1)?$|1_|sample-0[1-3]',
  'S2' = 'DIA-2|Sample [4-6]|2-|(128C|129N|129C)|^[5-7]( 1)?$|2_|sample-0[5-7]',
  'S3' = 'DIA-3|Sample [7-9]|3-|(130N|130C|131)|^([8-9]|10)( 1)?$|3_|sample-(0[8-9]|10)'
)

# TMT10-128N not used
if ('quantity.sample-04' %in% colnames(reports$FragPipe)) {
  reports$FragPipe$`quantity.sample-04` = NULL
}
if ('quantity.4' %in% colnames(reports$MaxQuant)) {
  reports$MaxQuant$`quantity.4` = NULL
}
if ('quantity.4 1' %in% colnames(reports$MaxQuant)) {
  reports$MaxQuant$`quantity.4 1` = NULL
}

run_groups = c(
  'S1' = 'FB1-DIA|Sample [1-3]|1-|(126|127N|127C)|^[1-3]( 1)?$|1_|sample-0[1-3]$',
  'S2' = 'FB2-DIA|Sample [4-6]|2-|(128N|128C|129N)|^[4-6]( 1)?$|2_|sample-0[4-6]$',
  'S3' = 'FB3-DIA|Sample [7-9]|3-|(129C|130N|130C)|^[7-9]( 1)?$|3_|sample-0[7-9]$'
)

# TMT10-131 not used
if ('quantity.sample-10' %in% colnames(reports$FragPipe)) {
  reports$FragPipe$`quantity.sample-10` = NULL
}
if ('quantity.sample-01131N' %in% colnames(reports$FragPipe)) {
  reports$FragPipe$`quantity.sample-01131N` = NULL
}
if ('quantity.10' %in% colnames(reports$MaxQuant)) {
  reports$MaxQuant$`quantity.10` = NULL
}
if ('quantity.10 1' %in% colnames(reports$MaxQuant)) {
  reports$MaxQuant$`quantity.10 1` = NULL
}




# assign organism
fasta_file = '12mixture bacteria.fasta'

fasta = seqinr::read.fasta(fasta_file, seqtype = 'AA', as.string = TRUE)

reports = lapply(reports, function(report) {
  if (!'organism' %in% colnames(report)) {
    pg_col = which(colnames(report) == 'protein')
    report = cbind(
      report[, 1:pg_col, drop = FALSE],
      organism = find.organism(report, fasta),
      report[, -(1:pg_col), drop = FALSE]
    )
  }
  report
})


append.xlsx(
  reports,
  sheets = sapply(names(reports), function(name) paste(name, level)), 
  file = 'report_matrices.xlsx'
)



# identification counts
plot.compare.venn.identification(
  reports, level, run_groups = run_groups,
  filename = paste0('compare_', paste(names(reports), collapse = '_'), '_venn_', level, '.svg')
)

ident_count = local({
  ident_count = lapply(names(reports), function(x) {
    run_ident_count = plot.run.identifications(
      reports[[x]], level, run_groups = run_groups,
      filename = paste0('', x, '_run_identifications_', level, '.svg'),
      return_data = TRUE
    )
    
    cumulative_ident_count = plot.cumulative.identifications(
      reports[[x]], level, run_groups = NULL,
      filename = paste0('', x, '_cumulative_identifications_', level, '.svg'),
      return_data = TRUE
    )
    
    list(individual = run_ident_count, cumulative = cumulative_ident_count)
  })
  names(ident_count) = names(reports)
  ident_count
})

invisible(lapply(names(ident_count), function(name) {
  append.xlsx(
    ident_count[[name]],
    sheets = paste0(name, ' ', level, c('s per run', 's cumulative')), 
    file = 'identification_count.xlsx'
  )
}))



# quantification CV
invisible(lapply(1:length(reports), function(i) {
  plot.cv.distribution(
    reports[[i]], run_groups = run_groups,
    filename = paste0('', names(reports)[i], '_cv_distribution_', level, '.svg'),
    fill = rep(c('#339dff', '#ff3355', '#65c3ba')[i], length(run_groups)),
    color = rep(c('#194e7f', '#7f192a', '#3c756f')[i], length(run_groups))
  )
}))

cv_distribution = plot.compare.cv.distribution(
  reports, run_groups = run_groups,
  filename = paste0('compare_', paste(names(reports), collapse = '_'), '_cv_distribution_', level, '.svg'),
  return_data = TRUE
)

append.xlsx(
  cv_distribution$summary,
  sheets = paste(level, 'CV', 'summary'), 
  file = 'cv_distribution.xlsx'
)
append.xlsx(
  cv_distribution$data,
  sheets = sapply(names(cv_distribution$data), function(name) paste(name, level, 'CV')), 
  file = 'cv_distribution.xlsx'
)



# fold change per organism
# 12 mix
theoretical_fold_change = list(
  'S2/S1' = list(
    "K. pneumoniae"    = 1/5,
    "B. fragilis"      = 5/1,
    "P. aeruginosa"    = 2/1,
    "M. morganii"      = 1/5,
    "E. casseliflavus" = 5/2,
    "C. butyricum"     = 2/1,
    "E. faecalis"      = 2/1,
    "E. asburiae"      = 1/5,
    "C. freundii"      = 5/2,
    "L. acidophilus"   = 2/5,
    "E. coli"          = 3/2,
    "K. aerogenes"     = 2/1
  ),
  'S3/S1' = list(
    "K. pneumoniae"    = 2/5,
    "B. fragilis"      = 2/1,
    "P. aeruginosa"    = 5/1,
    "M. morganii"      = 2/5,
    "E. casseliflavus" = 1/2,
    "C. butyricum"     = 5/1,
    "E. faecalis"      = 5/1,
    "E. asburiae"      = 2/5,
    "C. freundii"      = 1/2,
    "L. acidophilus"   = 1/5,
    "E. coli"          = 1/2,
    "K. aerogenes"     = 5/1
  )
)

# 6 mix
theoretical_fold_change = list(
  'S2/S1' = list(
    "P. aeruginosa"    = 2/1,
    "M. morganii"      = 1/5,
    "C. butyricum"     = 5/2,
    "E. faecalis"      = 5/1,
    "E. asburiae"      = 2/5,
    "L. acidophilus"   = 1/2
  ),
  'S3/S1' = list(
    "P. aeruginosa"    = 5/1,
    "M. morganii"      = 2/5,
    "C. butyricum"     = 1/2,
    "E. faecalis"      = 2/1,
    "E. asburiae"      = 1/5,
    "L. acidophilus"   = 5/2
  )
)

theoretical_fold_change = lapply(names(theoretical_fold_change), function(group) {
  data.frame(
    group = group,
    organism = names(theoretical_fold_change[[group]]),
    fc = as.numeric(theoretical_fold_change[[group]]),
    stringsAsFactors = FALSE
  )
})

theoretical_fold_change[[3]] = data.frame(
  group = 'S3/S2',
  organism = theoretical_fold_change[[1]]$organism,
  fc = theoretical_fold_change[[2]]$fc / theoretical_fold_change[[1]]$fc,
  stringsAsFactors = FALSE
)


fold_change = plot.compare.fold.change(
  reports, run_groups = run_groups, 
  theoretical_fold_change = theoretical_fold_change, min_frequency = 0.6,
  filename = paste0('compare_', paste(names(reports), collapse = '_'), '_foldchange_organism_', level, '.svg'),
  return_data = TRUE
)

append.xlsx(
  do.call(rbind, fold_change$summary),
  sheets = paste(level, 'FC', 'summary'), 
  file = 'fold_change_organism.xlsx'
)
append.xlsx(
  fold_change$data,
  sheets = sapply(names(fold_change$data), function(name) paste(name, level, 'FC')), 
  file = 'fold_change_organism.xlsx'
)


