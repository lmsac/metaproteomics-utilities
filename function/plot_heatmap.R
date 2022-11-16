library(ggplot2)

plot.heatmap = function(protein_report, protein_annotation,
                        filename = 'heatmap.svg', 
                        width = 12, height = 8, unit = 'cm',
                        hide_text = FALSE,
                        return_data = FALSE) {
  data = data.frame(
    protein_annotation,
    protein_report[
      match(protein_annotation$protein, sapply(strsplit(protein_report$PG.ProteinAccessions, ';'), head, 1)),
      grep('Quantity', colnames(protein_report)),
      drop = FALSE
    ],
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
    geom_tile(aes(x = sample, y = protein, fill = log10(quantity))) +
    scale_fill_gradientn(
      name = 'log10 quantity',
      rescaler = function(x, ...) scales::rescale_mid(
        x, ..., 
        mid = mean(log10(df$quantity))
      ),
      values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
      colors = c(
        '#323898', '#4474B6', '#75B0D0', '#AADAE7', '#E4F2F8',
        '#FFFEC1', 
        '#F7DE9A', '#FDAD62', '#F57246', '#D63026', '#A50024'
      ),
      limits = mean(log10(df$quantity)) + 
        c(-1, 1) * abs(max(log10(df$quantity) - mean(log10(df$quantity))))
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
      axis.text.y = if (hide_text) element_blank() else element_text(color = 'black'),
      legend.position = 'bottom',
      legend.key.size = unit(0.4, 'cm'), 
      legend.title = if (hide_text) element_blank() else element_text(size = 7), 
      legend.text = if (hide_text) element_blank() else element_text(size = 6)
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


plot.heatmap.kegg = function(protein_report, protein_annotation, ko_list, pathway = NULL, ...) {
  if (!is.null(pathway) && !is.na(pathway[1])) {
    protein_annotation = subset(
      protein_annotation, 
      grepl(pathway, KEGG_Pathway)
    )
  }
  
  if (is.null(ko_list) || is.na(ko_list[1])) {
    ko_list = sub('ko:', '', unique(unlist(strsplit(protein_annotation$KEGG_ko, ','))))
  }
  
  data = do.call(rbind, lapply(ko_list, function(k) {
    data = subset(
      protein_annotation, 
      grepl(k, KEGG_ko),
      select = c('protein', 'eggNOG_OGs', 'KEGG_ko')
    )
    data$KEGG_ko = gsub('ko:', '', data$KEGG_ko)
    data$taxonomy = sapply(strsplit(data$eggNOG_OGs, ','), function(s) sub('^[^|]+\\|', '',tail(s, 1)))
    data = data[order(data$taxonomy), ]
  }))
  
  data = subset(data, !duplicated(protein))
  data$eggNOG_OGs = NULL
  data$description = paste0(data$KEGG_ko, ' | ', data$taxonomy)
  
  plot.heatmap(protein_report, data, ...)
}


plot.heatmap.cog = function(protein_report, protein_annotation, cog_list, cog_category = NULL, ...) {
  if (!is.null(cog_category) && !is.na(cog_category[1])) {
    protein_annotation = subset(
      protein_annotation, 
      grepl(cog_category, COG_category)
    )
  }
  
  if (is.null(cog_list) || is.na(cog_list[1])) {
    cog_list = unique(sub('^(COG[0-9]+).*', '\\1', protein_annotation$eggNOG_OGs))
  }
  
  data = do.call(rbind, lapply(cog_list, function(cog) {
    data = subset(
      protein_annotation, 
      grepl(cog, eggNOG_OGs),
      select = c('protein', 'eggNOG_OGs')
    )
    data$COG = cog
    data$taxonomy = sapply(strsplit(data$eggNOG_OGs, ','), function(s) sub('^[^|]+\\|', '',tail(s, 1)))
    data = data[order(data$taxonomy), ]
  }))
  
  data = subset(data, !duplicated(protein))
  data$eggNOG_OGs = NULL
  data$description = paste0(data$COG, ' | ', data$taxonomy)

  plot.heatmap(protein_report, data, ...)
}



library(readr)


differential_protein_annotation = read_csv('differential_protein_annotation.csv')
protein_report = read_csv('protein_report.csv')

protein_report$PG.ProteinAccessions = sub('^(\\S+)\\s.*', '\\1', protein_report$PG.ProteinAccessions)

cog_list = c(
  'COG0443',
  'COG0484',
  'COG0542',
  'COG0466',
  'COG0740',
  'COG1220',
  'COG1225',
  'COG2077'
)

data_heatmap = plot.heatmap.cog(
  protein_report,
  differential_protein_annotation,
  cog_list = cog_list, cog_category = 'O',
  filename = 'heatmap_cog_O_text.svg',
  width = 8, height = 8, unit = 'cm',
  return_data = TRUE
)

write.csv(data_heatmap, 'heatmap_cog_O.csv', row.names = F)



ko_list = c(
  # 'K01067', # acetyl-CoA hydrolase [EC:3.1.2.1]
  # 'K00925', # ackA, acetate kinase [EC:2.7.2.1]
  # 'K13788', 'K00625', 'K15024', # pta, phosphate acetyltransferase [EC:2.3.1.8]
  'K11262', 'K01946', 'K01962', 'K02160', 'K01961', 'K01963', 'K18472',
  # ACACA/ACACB/accA/accB/accC/accD/accD6, acetyl-CoA carboxylase [EC:6.4.1.2]
  'K00645', # fabD, S-malonyltransferase [EC:2.3.1.39]
  'K00648', 'K18473',  
  # fabH, fabY 3-oxoacyl-[acyl-carrier-protein] synthase III [EC:2.3.1.180]
  # fabY, acetoacetyl-[acyl-carrier protein] synthase
  'K00647', # fabB, 3-oxoacyl-[acyl-carrier-protein] synthase I [EC:2.3.1.41]
  'K09458', # fabF, 3-oxoacyl-[acyl-carrier-protein] synthase II [EC:2.3.1.179]
  'K00059', # fabG, 3-oxoacyl-[acyl-carrier protein] reductase [EC:1.1.1.100]
  'K01716', 'K02372', 'K16363',
  # fabA/fabZ/lpxC-fabZ, 3-hydroxyacyl-[acyl-carrier protein] dehydratase [EC:4.2.1.59]
  'K10780', # fabL, enoyl-[acyl-carrier protein] reductase III [EC:1.3.1.104]
  'K00208', 'K02371', 'K00209',
  # fabI/fabK/fabV, enoyl-[acyl-carrier protein] reductase I/II/ [EC:1.3.1.9]
  'K10781', 'K10782',
  # FATB/FATA, fatty acyl-ACP thioesterase B/A [EC:3.1.2.21 3.1.2.14]
  'K01897' # ACSL, long-chain acyl-CoA synthetase [EC:6.2.1.3]
)

data_heatmap = plot.heatmap.kegg(
  protein_report,
  differential_protein_annotation,
  ko_list = ko_list,
  filename = 'heatmap_kegg_fatty_acid.svg',
  width = 8, height = 8, unit = 'cm',
  # hide_text = TRUE,
  return_data = TRUE
)

write.csv(data_heatmap, 'heatmap_kegg_fatty_acid.csv', row.names = F)



data_heatmap = plot.heatmap.kegg(
  protein_report,
  differential_protein_annotation,
  ko_list = NULL, pathway = 'ko00290',
  filename = 'heatmap_kegg_BCAA.svg',
  width = 8, height = 8, unit = 'cm',
  # hide_text = TRUE,
  return_data = TRUE
)

write.csv(data_heatmap, 'heatmap_kegg_BCAA.csv', row.names = F)
