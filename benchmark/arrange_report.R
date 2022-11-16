check.spectronaut.report = function(report, level) {
  if (level == 'protein') {
    any(colnames(report) == 'PG.ProteinAccessions') &&
      any(grepl('\\.PG\\.Quantity$', colnames(report)))
  } else if (level == 'peptide') {
    any(colnames(report) == 'PEP.StrippedSequence') &&
      any(grepl('\\.PEP\\.Quantity$', colnames(report)))
  } else {
    FALSE
  }
}

arrange.spectronaut.report = function(report, level) {
  data = NULL
  
  pg_col = which(colnames(report) == 'PG.ProteinAccessions')
  if (any(pg_col)) {
    protein = report[[pg_col]]
    protein = gsub('(^|;)([a-z]{2})\\|([^|;]+)\\|([^;]*)', '\\1\\3', protein)
    if (is.null(data)) {
      data = data.frame(protein = protein, stringsAsFactors = FALSE)
    } else {
      data$protein = protein
    }
  }
  
  if (level == 'peptide') {
    pep_col = which(colnames(report) == 'PEP.StrippedSequence')
    peptide = report[[pep_col]]
    if (is.null(data)) {
      data = data.frame(peptide = peptide, stringsAsFactors = FALSE)
    } else {
      data$peptide = peptide
    }
  }
  
  quant_col = grep('\\.Quantity$', colnames(report), value = TRUE)
  quantity = sapply(quant_col, function(col) {
    quant = report[[col]]
    if (!is.numeric(quant)) {
      quant = as.numeric(ifelse(quant == 'Filtered', NA, quant))
    }
    quant[quant == 0] = NA
    quant
  })
  colnames(quantity) = paste0('quantity.', gsub(
    '^\\[[0-9]+\\] |\\.(PEP|PG)\\.Quantity$', '', quant_col
  ))
  
  data = cbind(data, quantity)
  data = subset(data, rowSums(!is.na(data[grep(
    '^quantity\\.', colnames(data), value = TRUE
  )])) > 0)
  data
}


check.diann.report = function(report, level) {
  if (level == 'protein') {
    any(colnames(report) == 'Protein.Group')
  } else if (level == 'peptide') {
    any(colnames(report) == 'Stripped.Sequence')
  } else {
    FALSE
  } &&
    any(grepl('\\.raw$', colnames(report)))
}

arrange.diann.report = function(report, level) {
  data = NULL
  
  pg_col = which(colnames(report) == 'Protein.Group')
  if (any(pg_col)) {
    protein = report[[pg_col]]
    protein = gsub('(^|;)([a-z]{2})\\|([^|;]+)\\|([^;]*)', '\\1\\3', protein)
    if (is.null(data)) {
      data = data.frame(protein = protein, stringsAsFactors = FALSE)
    } else {
      data$protein = protein
    }
  }
  
  if (level == 'peptide') {
    pep_col = which(colnames(report) == 'Stripped.Sequence')
    peptide = report[[pep_col]]
    if (is.null(data)) {
      data = data.frame(peptide = peptide, stringsAsFactors = FALSE)
    } else {
      data$peptide = peptide
    }
  }
  
  quant_col = grep('\\.raw$', colnames(report), value = TRUE)
  quantity = sapply(quant_col, function(col) {
    quant = report[[col]]
    if (!is.numeric(quant)) {
      quant = as.numeric(ifelse(quant == '', NA, quant))
    }
    quant[quant == 0] = NA
    quant
  })
  colnames(quantity) = paste0('quantity.', basename(quant_col))
  
  if (level == 'peptide') {
    # precursor to peptide
    quantity = aggregate(quantity, by = list(peptide), FUN = sum, na.rm = TRUE)
    data = subset(data, !duplicated(peptide))
    quantity = quantity[match(data$peptide, quantity[[1]]), -1, drop = FALSE]
    quantity[quantity == 0] = NA
  }
  
  data = cbind(data, quantity)
  data = subset(data, rowSums(!is.na(data[grep(
    '^quantity\\.', colnames(data), value = TRUE
  )])) > 0)
  data
}


check.peaks.report = function(report, level) {
  if (level == 'protein') {
    any(colnames(report) =='Accession') && (
      any(grepl('Sample [0-9]+ Area \\(top-3 peptides\\)', colnames(report))) ||
        any(grepl('^Intensity TMT', colnames(report)))
    )
  } else if (level == 'peptide') {
    any(colnames(report) == 'Peptide') && (
      any(grepl('Sample [0-9]+$', colnames(report))) ||
        any(grepl('^Intensity TMT', colnames(report)))
    )
  } else {
    FALSE
  }
}

arrange.peaks.report = function(report, level) {
  data = NULL
  
  pg_col = head(which(colnames(report) %in% c('Accession', 'Protein Accession')), 1)
  if (any(pg_col)) {
    protein = report[[pg_col]]
    protein = ifelse(
      grepl('^([a-z]{2})\\|([^|]+)\\|(.*)', protein),
      sub('^([a-z]{2})\\|([^|]+)\\|(.*)', '\\2', protein),
      sub('^([^|]+)\\|(.*)', '\\1', protein)
    )
    
    if (level == 'protein') {
      group_col = which(colnames(report) == 'Protein Group')
      if (any(group_col)) {
        leading = !duplicated(report[[group_col]])
        protein = aggregate(
          protein, by = list(report[[group_col]]), 
          FUN = function(x) {
            if (length(unique(x)) == 1) x[1] else paste(x, collapse = ';')
          }
        )
        protein = protein[match(report[[group_col]][leading], protein[, 1]), 2]
        
        report = subset(report, leading)
      }
    }
    
    if (is.null(data)) {
      data = data.frame(protein = protein, stringsAsFactors = FALSE)
    } else {
      data$protein = protein
    }
  }
  
  if (level == 'peptide') {
    pep_col = which(colnames(report) == 'Peptide')
    peptide = report[[pep_col]]
    peptide = gsub(
      '\\([+\\-][0-9.]+\\)', '',
      gsub('(^[A-Z]\\.)|(\\.[A-Z]$)', '', peptide)
    )
    if (is.null(data)) {
      data = data.frame(peptide = peptide, stringsAsFactors = FALSE)
    } else {
      data$peptide = peptide
    }
    
    if ('protein' %in% colnames(data)) {
      leading = !duplicated(peptide)
      protein = aggregate(
        data$protein, by = list(peptide), 
        FUN = function(x) {
          if (length(unique(x)) == 1) x[1] else paste(x, collapse = ';')
        }
      )
      data = subset(data, leading)
      data$protein = protein[match(peptide[leading], protein[, 1]), 2]
      report = subset(report, leading)
    }
  }
  
  if (level == 'protein') {
    quant_col = grep('Sample [0-9]+ Area \\(top-3 peptides\\)', colnames(report), value = TRUE)
  } else {
    quant_col = grep('Sample [0-9]+$', colnames(report), value = TRUE)
  }
  if (length(quant_col) == 0) {
    quant_col = grep('^Intensity TMT', colnames(report), value = TRUE)
  }
  quantity = sapply(quant_col, function(col) {
    quant = report[[col]]
    if (!is.numeric(quant)) {
      quant = as.numeric(ifelse(quant == '', NA, quant))
    }
    quant[quant == 0] = NA
    quant
  })
  colnames(quantity) = paste0('quantity.', sub(' Area \\(top-3 peptides\\)$| $|^Intensity ', '', quant_col))
  
  data = cbind(data, quantity)
  data = subset(data, rowSums(!is.na(data[grep(
    '^quantity\\.', colnames(data), value = TRUE
  )])) > 0)
  data
}


check.maxquant.report = function(report, level) {
  if (level == 'protein') {
    any(colnames(report) == 'Protein IDs')
  } else if (level == 'peptide') {
    any(colnames(report) == 'Sequence')
  } else {
    FALSE
  } && (
    any(grepl('^LFQ intensity ', colnames(report))) ||
      any(grepl('^Reporter intensity corrected ', colnames(report)))
  )
}

arrange.maxquant.report = function(report, level) {
  data = NULL
  
  report = subset(report, is.na(Reverse) & is.na(`Potential contaminant`))
  
  # duplicated columns 'Reporter intensity corrected '
  report = report[, !endsWith(colnames(report), ' 1') | !(sub(' 1$', '', colnames(report)) %in% colnames(report))]
  
  pg_col = head(which(colnames(report) %in% c('Protein IDs', 'Proteins')), 1)
  if (any(pg_col)) {
    protein = report[[pg_col]]
    protein = gsub('(^|;)([a-z]{2})\\|([^|;]+)\\|([^;]*)', '\\1\\3', protein)
    if (is.null(data)) {
      data = data.frame(protein = protein, stringsAsFactors = FALSE)
    } else {
      data$protein = protein
    }
  }
  
  if (level == 'peptide') {
    pep_col = which(colnames(report) == 'Sequence')
    peptide = report[[pep_col]]
    if (is.null(data)) {
      data = data.frame(peptide = peptide, stringsAsFactors = FALSE)
    } else {
      data$peptide = peptide
    }
  }
  
  quant_col = grep('^LFQ intensity ', colnames(report), value = TRUE)
  if (length(quant_col) == 0) {
    quant_col = grep('^Reporter intensity corrected ', colnames(report), value = TRUE)
  }
  quantity = sapply(quant_col, function(col) {
    quant = report[[col]]
    if (!is.numeric(quant)) {
      quant = as.numeric(ifelse(quant == '', NA, quant))
    }
    quant[quant == 0] = NA
    quant
  })
  colnames(quantity) = paste0('quantity.', sub('^LFQ intensity |^Reporter intensity corrected ', '', quant_col))
  
  data = cbind(data, quantity)
  data = subset(data, rowSums(!is.na(data[grep(
    '^quantity\\.', colnames(data), value = TRUE
  )])) > 0)
  data
}


check.fragpipe.report = function(report, level) {
  if (level == 'protein') {
    any(colnames(report) %in% c('Protein ID', 'Index'))
  } else if (level == 'peptide') {
    any(colnames(report) %in% c('Peptide Sequence', 'Peptide'))
  } else {
    FALSE
  } && (
    any(grepl('MaxLFQ Intensity$', colnames(report))) ||
      any(grepl('^sample-', colnames(report)))
  )
}

arrange.fragpipe.report = function(report, level) {
  data = NULL
  
  pg_col = head(which(colnames(report) %in% c('ProteinID', 'Protein ID')), 1)
  if (!any(pg_col) && level == 'protein') {
    pg_col = which(colnames(report) == 'Index')
  }
  if (any(pg_col)) {
    report = subset(report, !grepl('contam', report[[pg_col]]))
    protein = report[[pg_col]]
    protein = ifelse(
      grepl('^([a-z]{2})\\|([^|]+)\\|(.*)', protein),
      sub('^([a-z]{2})\\|([^|]+)\\|(.*)', '\\2', protein),
      sub('^([^|]+)\\|(.*)', '\\1', protein)
    )
    if (is.null(data)) {
      data = data.frame(protein = protein, stringsAsFactors = FALSE)
    } else {
      data$protein = protein
    }
  }
  
  if (level == 'peptide') {
    pep_col = head(which(colnames(report) %in% c('Peptide Sequence', 'Peptide')), 1)
    peptide = report[[pep_col]]
    if (is.null(data)) {
      data = data.frame(peptide = peptide, stringsAsFactors = FALSE)
    } else {
      data$peptide = peptide
    }
  }
  
  quant_col = grep('^sample-', colnames(report), value = TRUE)
  if (length(quant_col) == 0) {
    quant_col = grep('MaxLFQ Intensity$', colnames(report), value = TRUE)
  }
  quantity = sapply(quant_col, function(col) {
    quant = report[[col]]
    if (!is.numeric(quant)) {
      quant = as.numeric(ifelse(quant == '', NA, quant))
    }
    quant[quant == 0] = NA
    quant
  })
  if (grepl('^sample-', quant_col[1]) && max(quantity, na.rm = TRUE) < 100) {
    # TMT report: log2 intensity
    quantity = 2 ^ quantity
  }
  colnames(quantity) = paste0('quantity.', sub(' MaxLFQ Intensity$', '', quant_col))
  
  data = cbind(data, quantity)
  data = subset(data, rowSums(!is.na(data[grep(
    '^quantity\\.', colnames(data), value = TRUE
  )])) > 0)
  data
}


find.organism = function(report, fasta, use_abbr = TRUE) {
  accession = sub('^([a-z]{2})\\|([^|]+)\\|(.*)', '\\2', attr(fasta, 'name'))
  
  organism = sapply(strsplit(report$protein, ';'), function(x) {
    idx = match(x, accession)
    organism = sapply(idx, function(i) {
      if (!is.na(i)) {
        desc = attr(fasta[[i]], 'Annot')
        os = sub('^.* OS=([^=]+) [A-Z]+=.*', '\\1', desc)
        os = ifelse(os == desc, 'Unknown', os)
      } else {
        'Unknown'
      }
    })
    
    if (use_abbr) {
      organism = sub('^(\\S)\\S+\\s(\\S+).*', '\\1. \\2', organism)
    }
    
    if (length(unique(organism)) == 1) {
      organism[1] 
    } else {
      paste(organism, collapse = ';')
    } 
  })
  
  organism
}


append.xlsx = function(reports, file, sheets = names(reports)) {
  library(openxlsx)
  
  if (file.exists(file)) {
    wb = loadWorkbook(file)
  }
  else {
    wb = createWorkbook()
  }
  
  if (is.data.frame(reports)) {
    reports = list(reports)
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

