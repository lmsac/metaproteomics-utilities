
append.xlsx = function(reports, file, sheets = names(reports)) {
  library(openxlsx)
  
  if (file.exists(file)) {
    wb = loadWorkbook(file)
  }
  else {
    wb = createWorkbook()
  }
  
  if (length(reports) > 0) {
    lapply(1:length(reports), function(i) {
      sheet = sheets[i]
      addWorksheet(wb, sheet)
      writeData(wb, sheet, reports[[i]])
    })
  }
  
  saveWorkbook(wb, file = file, overwrite = TRUE)
}



library(readr)

level = 'protein'
report_files = list(
  'directDIA' = 'DIA/directDIA/protein_report.csv',
  'LFQ-DDA' = 'DDA/PEAKS/proteins.csv',
  'TMT' = 'TMT/PEAKS/proteins.csv'
)
# report_files = list(
#   'DDALib' = 'DIA/DDALib/protein_report.csv',
#   'directDIA' = 'DIA/directDIA/protein_report.csv'
# )
# report_files = list(
#   'All' = 'DIA/directDIA/protein_report.csv',
#   'Ribosomal' = 'DIA/directDIA_ribosomal/protein_report.csv'
# )


# level = 'peptide'
# report_files = list(
#   'directDIA' = 'DIA/directDIA/peptide_report.csv',
#   'LFQ-DDA' = 'DDA/PEAKS/protein-peptides.csv',
#   'TMT' = 'TMT/PEAKS/protein-peptides.csv'
# )
# report_files = list(
#   'DDALib' = 'DIA/DDALib/peptide_report.csv',
#   'directDIA' = 'DIA/directDIA/peptide_report.csv'
# )


reports = lapply(report_files, read_csv)

reports = lapply(reports, function(report) {
  if (any(c('PEP.StrippedSequence', 'PG.ProteinAccessions') %in% colnames(report))) {
    report = arrange.spectromine.report(report)
  }
  else if (any(c('Peptide', 'Accession') %in% colnames(report))) {
    report = arrange.peaks.report(report)
  }
  else {
    stop('unknown format')
  }
})






plot.venn(
  reports, 
  id = level, 
  filename = paste0('venn_', level, '.svg'), 
  hide.text = TRUE
)

append.xlsx(
  reports,
  sheets = sapply(names(reports), function(name) paste(name, level)), 
  file = 'report.xlsx'
)


run_groups = list(
  'S1' = '\\.[1-3]$',
  'S2' = '\\.[4-6]$',
  'S3' = '\\.[7-9]$'
)

cv = plot.cv.distribution(
  reports, 
  run_groups = run_groups, 
  filename = paste0('density_cv_', level, '.svg'), 
  hide.text = TRUE,
  return_data = TRUE
)

append.xlsx(
  cv,
  sheets = sapply(names(cv), function(name) paste(name, level, 'CV')), 
  file = 'cv.xlsx'
)



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


# theoretical_fold_change = list(
#   'S2/S1' = list(
#      "B. fragilis"      = 5/1,
#      "P. aeruginosa"    = 2/1,
#      "E. casseliflavus" = 5/2
#   ),
#   'S3/S1' = list(
#      "K. pneumoniae"    = 2/5,
#      "M. morganii"      = 2/5,
#      "C. butyricum"     = 5/1,
#      "E. faecalis"      = 5/1,
#      "E. asburiae"      = 2/5,
#      "C. freundii"      = 1/2,
#      "E. coli"          = 1/2
#   ),
#   'S3/S2' = list(
#     "L. acidophilus"   = 1/2,
#     "K. aerogenes"     = 5/2
#   )
# )


# theoretical_fold_change = list(
#   'S2/S1' = list(
#     "P. aeruginosa"    = 1/2,
#     "M. morganii"      = 5/1,
#     "C. butyricum"     = 2/5,
#     "E. faecalis"      = 1/5,
#     "E. asburiae"      = 5/2,
#     "L. acidophilus"   = 2/1
#   ),
#   'S3/S1' = list(
#     "P. aeruginosa"    = 5/2,
#     "M. morganii"      = 2/1,
#     "C. butyricum"     = 1/5,
#     "E. faecalis"      = 2/5,
#     "E. asburiae"      = 1/2,
#     "L. acidophilus"   = 5/1
#   ),
#   'S3/S2' = list(
#     "P. aeruginosa"    = 5/1,
#     "M. morganii"      = 2/5,
#     "C. butyricum"     = 1/2,
#     "E. faecalis"      = 2/1,
#     "E. asburiae"      = 1/5,
#     "L. acidophilus"   = 5/2
#   )
# )



theoretical_fold_change = lapply(names(theoretical_fold_change), function(group) {
  data.frame(
    group = group,
    organism = names(theoretical_fold_change[[group]]),
    fc = as.numeric(theoretical_fold_change[[group]]),
    stringsAsFactors = FALSE
  )
})

# theoretical_fold_change[[3]] = data.frame(
#   group = 'S3/S2',
#   organism = theoretical_fold_change[[1]]$organism,
#   fc = theoretical_fold_change[[2]]$fc / theoretical_fold_change[[1]]$fc,
#   stringsAsFactors = FALSE
# )
  
# theoretical_fold_change = do.call(rbind, theoretical_fold_change)


fold_change = plot.fold.change(
  reports, 
  run_groups = run_groups, 
  theoretical_fold_change = theoretical_fold_change,
  # box_width = 0.20, digits = 3, nrow = 2,
  filename = paste0('fc_organsim_', level, '.svg'), 
  hide.text = TRUE,
  return_data = TRUE
)

append.xlsx(
  fold_change,
  sheets = sapply(names(fold_change), function(name) paste(name, level, 'FC')), 
  file = 'fc.xlsx'
)


count_organsim = plot.count.organism(
  reports, 
  filename = paste0('bar_organsim_', level, '.svg'), 
  hide.text = TRUE,
  return_data = TRUE
)
