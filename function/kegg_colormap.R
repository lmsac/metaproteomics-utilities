get.kegg_colormap = function(protein_annotation, taxonomy,
                             color = c('up' = 'red', 'down' = 'blue', 'shared' = 'green')) {
  kegg_colormap = subset(
    protein_annotation,
    grepl(taxonomy, eggNOG_OGs) & KEGG_ko != '-',
    select = c('KEGG_ko', 'log2fc')
  )
  
  kegg_colormap = cbind(
    ko = gsub('ko:', '', kegg_colormap$KEGG_ko),
    color = ifelse(kegg_colormap$log2fc > 0, color['up'], color['down'])
  )
  
  kegg_colormap = do.call(rbind, apply(kegg_colormap, 1, function(x) {
    if (grepl(',', x[1])) {
      s = strsplit(x[1], ',')[[1]]
      cbind(ko = s, color = rep(x[2], length(s)))
    }
    else {
      x
    }
  }))
  
  kegg_colormap = aggregate(color ~ ko, kegg_colormap, function(x) {
    if (length(unique(x)) > 1) {
      color['shared']
    }
    else {
      x[1]
    }
  })
  
  kegg_colormap
}



library(readr)


differential_protein_annotation = read_csv('differential_protein_annotation.csv')


taxonomy = 'Bacteroidetes'

kegg_colormap = get.kegg_colormap(differential_protein_annotation, taxonomy)

write.table(
  kegg_colormap, 
  paste0('kegg_colormap_', taxonomy, '.txt'), 
  row.names = FALSE, quote = FALSE, col.names = FALSE
)
