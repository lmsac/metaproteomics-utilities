arrange.peaks.report = function(report) {
  if ('Peptide' %in% colnames(report)) {
    peptide = gsub(
      '\\([+\\-][0-9.]+\\)', '',
      gsub('(^[A-Z]\\.)|(\\.[A-Z]$)', '', report$Peptide)
    )
    
    report = subset(report, !duplicated(peptide))
    peptide = peptide[!duplicated(peptide)]
    protein = report$`Protein Accession`
    
    analyte = data.frame(
      peptide = peptide,
      protein = protein,
      stringsAsFactors = FALSE
    )
    
    run_name_pattern = 'Sample [0-9]+'
  }
  else {
    report = subset(report, !duplicated(`Protein Group`))
    
    protein = report$Accession
    organism = sub('^.* OS=([^=]+) [A-Z]+=.*', '\\1', report$Description)
    
    analyte = data.frame(
      protein = protein,
      organism = organism,
      stringsAsFactors = FALSE
    )
    
    run_name_pattern = 'Sample [0-9]+ Area \\(top-3 peptides\\)'
  }
  
  analyte$protein = ifelse(
    grepl('^(sp|tr)\\|([^|]+)\\|(.*)', analyte$protein),
    sub('^(sp|tr)\\|([^|]+)\\|(.*)', '\\2', analyte$protein),
    analyte$protein
  )
  
  run_name = grep(
    run_name_pattern, 
    colnames(report), value = TRUE
  )
  if (length(run_name) == 0) {
    run_name = grep(
      'Intensity TMT',
      colnames(report), value = TRUE
    )
  }
  quantity = subset(report, select = run_name)
  lapply(run_name, function(name) {
    if (is.character(quantity[[name]])) {
      quantity[[name]] <<- as.numeric(quantity[[name]])
    }
  })
  colnames(quantity) = paste('quantity', 1:ncol(quantity))
    
  data.frame(
    analyte,
    quantity,
    stringsAsFactors = FALSE
  )
}

arrange.spectromine.report = function(report) {
  if ('PEP.StrippedSequence' %in% colnames(report)) {
    peptide = report$PEP.StrippedSequence
    protein = sapply(strsplit(report$PG.ProteinAccessions, ';'), head, 1)
    
    analyte = data.frame(
      peptide = peptide,
      protein = protein,
      stringsAsFactors = FALSE
    )
  }
  else {
    protein = sapply(strsplit(report$PG.ProteinAccessions, ';'), head, 1)
    organism = sapply(strsplit(report$PG.Organisms, ';'), head, 1)
    
    analyte = data.frame(
      protein = protein,
      organism = organism,
      stringsAsFactors = FALSE
    )
  }
  
  run_name = grep(
    '\\.(PEP|PG)\\.TMT', 
    colnames(report), value = TRUE
  )
  is_tmt = length(run_name) > 0
  if (length(run_name) == 0) {
    run_name = grep(
      '\\.(PEP|PG)\\.Label-Free Quant', 
      colnames(report), value = TRUE
    )
  } 
  if (length(run_name) == 0) {
    run_name = grep(
      '\\.(PEP|PG)\\.Quantity', 
      colnames(report), value = TRUE
    )
  }
    
  quantity = subset(report, select = run_name)
  lapply(run_name, function(name) {
    if (is.character(quantity[[name]])) {
      quantity[[name]] <<- as.numeric(quantity[[name]])
    }
  })
  
  if (is_tmt) {
    quantity = sapply(unique(stringr::str_extract(run_name, '\\.(PEP|PG)\\.TMT.*')), function(label) {
      intensity = rowSums(
        subset(quantity, select = grep(label, run_name, fixed = TRUE, value = TRUE)),
        na.rm = TRUE
      )
      intensity[intensity == 0] = NA
      intensity
    })
  }
  
  colnames(quantity) = paste('quantity', 1:ncol(quantity))
  
  result = data.frame(
    analyte,
    quantity,
    stringsAsFactors = FALSE
  )
  result = subset(result, rowSums(!is.na(quantity)) > 0)
  result
}
