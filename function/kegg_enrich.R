library(httr)

get_http_content = function (url, ...) {
  res = GET(url, ...)
  content(res)
}

get_kegg_table = function (api, columns) {
  cfg = user_agent('Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0')
  s = get_http_content(paste0('http://rest.kegg.jp/', api), config = cfg)
  data = read.delim(text = s, sep = '\t', header = FALSE, quote = '')
  colnames(data) = columns
  data
}


KEGG_organism = get_kegg_table('list/organism', columns = c('entry', 'code', 'name', 'lineage'))

KEGG_pathway = get_kegg_table('list/pathway/ko', columns = c('entry', 'name'))
KEGG_pathway$number = sub('path:ko', '', KEGG_pathway$entry, fixed = TRUE)

KEGG_ko = get_kegg_table('list/ko', columns = c('entry', 'name'))

KEGG_ko_pathway = get_kegg_table('link/pathway/ko', columns = c('ko', 'pathway'))


filter_pathway_by_organism = function(pathway, organism) {
  pathway$found_in_org = Reduce('|', lapply(organism, function(org) {
    message(paste0(org, '\t'), appendLF = FALSE)
    pathway_org = get_kegg_table(paste0('list/pathway/', org), columns = c('entry', 'name'))
    pathway_org$number = sub(paste0('path:', org), '', pathway_org$entry, fixed = TRUE)
    !is.na(match(pathway$number, pathway_org$number))
  }))
  
  pathway
}

filter_ko_by_organism = function(ko, organism) {
  ko$found_in_org = Reduce('|', lapply(organism, function(org) {
    message(paste0(org, '\t'), appendLF = FALSE)
    ko_org = get_kegg_table(paste0('link/ko/', org), columns = c(org, 'ko'))
    !is.na(match(ko$entry, ko_org$ko))
  }))
  
  ko
}


get_kegg_by_taxonomy = function(taxonomy) {
  selected_organism = KEGG_organism$code[grep(taxonomy, KEGG_organism$lineage)]
  
  pathway = filter_pathway_by_organism(KEGG_pathway, selected_organism)
  ko = filter_ko_by_organism(KEGG_ko, selected_organism)
  
  ko_pathway = cbind(
    KEGG_ko_pathway,
    found_in_org = KEGG_ko_pathway$ko %in% ko$entry[ko$found_in_org] & 
      KEGG_ko_pathway$pathway %in% pathway$entry[pathway$found_in_org]
  )
  
  list(
    organism = selected_organism,
    pathway = pathway,
    ko = ko,
    ko_pathway = ko_pathway
  )
}



taxonomy_list = c(
  'Actinobacteria',
  'Bacteroidetes',
  'Bacilli',
  'Clostridia',
  'Negativicutes',
  'Gammaproteobacteria'
)

invisible(lapply(taxonomy_list, function(taxonomy) {
  message(taxonomy)
  kegg_tables = get_kegg_by_taxonomy(taxonomy)
  
  pathway = kegg_tables$pathway
  ko = kegg_tables$ko
  ko_pathway = kegg_tables$ko_pathway
  
  write.csv(pathway, paste0('pathway_', taxonomy, '.csv'), row.names = FALSE)
  write.csv(ko, paste0('ko_', taxonomy, '.csv'), row.names = FALSE)
  write.csv(ko_pathway, paste0('ko_pathway_', taxonomy, '.csv'), row.names = FALSE)
}))





library(readr)
library(clusterProfiler)
library(cowplot)
library(ggplot2)


differential_protein_annotation = read_csv('differential_protein_annotation.csv')


enrich_results = lapply(taxonomy_list, function(taxonomy) {
  ko_pathway = read_csv(paste0('ko_pathway_', taxonomy, '.csv'))
  pathway = read_csv(paste0('pathway_', taxonomy, '.csv'))
  
  selected_ko = differential_protein_annotation$KEGG_ko[grep(
    taxonomy, 
    differential_protein_annotation$eggNOG_OGs, 
    fixed = TRUE
  )]
  selected_ko = setdiff(unlist(strsplit(selected_ko, ';')), '-')
  
  enrich_result = enricher(
    gene = sub('ko:', '', selected_ko),
    TERM2GENE = cbind(
      sub('path:', '', ko_pathway$pathway[ko_pathway$found_in_org]),
      sub('ko:', '', ko_pathway$ko[ko_pathway$found_in_org])
    ),
    TERM2NAME = cbind(
      sub('path:', '', pathway$entry[pathway$found_in_org]),
      pathway$name[pathway$found_in_org]
    )
  )
})
names(enrich_results) = taxonomy_list


plots = lapply(enrich_results, function(x) {
  pl = dotplot(x)
  pl = pl + theme(
    panel.background = element_blank(), 
    panel.border = element_blank(), 
    legend.position = 'top', 
    legend.title = element_blank(), 
    axis.title = element_blank(), 
    axis.line.x = element_line(color = 'black'), 
    axis.ticks.y = element_blank()
  )
})
plots = align_plots(plotlist = plots, align = "v")
names(plots) = taxonomy_list


invisible(lapply(names(plots), function(taxonomy) {
  svg(
    filename = paste0('kegg_enrich_', taxonomy, '.svg'), 
    width = 8, height = 4, 
    bg = 'transparent'
  )
  plot_grid(plots[[taxonomy]])
  dev.off()
}))


invisible(lapply(names(enrich_results), function(taxonomy) {
  write.csv(
    enrich_results[[taxonomy]]@result, 
    paste0('kegg_enrich_', taxonomy, '.csv'), 
    row.names = FALSE
  )
}))



taxonomy_color = c(
  "Actinobacteria" = "#ff9933",
  "Bacteroidetes" = "#009926",
  "Bacilli;Clostridia" = "#9933cc",
  "Bacilli" = "#0000ee",
  "Clostridia" = "#ff3366",
  "Negativicutes" = "#cccc00"
)

# taxonomy_color = c(
#   "Actinobacteria" = "#ff9933",
#   "Bacteroidetes" = "#009926",  
#   "Clostridia" = "#ff3366",
#   "Gammaproteobacteria" = "#0000ee"
# )


taxonomy_pathway = lapply(names(taxonomy_color), function(taxonomy) {
  data = read_csv(paste0('kegg_enrich_', taxonomy, '.csv'))
  data = enrich_results[[taxonomy]]@result
  data = subset(data, p.adjust < 0.05 & Count >= 3)
  data$taxonomy = taxonomy
  data
})
taxonomy_pathway = taxonomy_pathway[sapply(taxonomy_pathway, nrow) > 0]

taxonomy_ko = do.call(rbind, lapply(taxonomy_pathway, function(data) {
  cbind(
    taxonomy = data$taxonomy[1],
    ko = unique(unlist(strsplit(data$geneID, '/')))
  )
}))
taxonomy_ko = aggregate(taxonomy ~ ko, taxonomy_ko, function(x) {
  if(length(x) == 1) {
    x
  }
  else {
    paste0(x, collapse = ';')
  }
})

taxonomy_ko$color = taxonomy_color[taxonomy_ko$taxonomy]
taxonomy_ko$color[is.na(taxonomy_ko$color)] = "#000000"

write.table(
  subset(taxonomy_ko, select = c('ko', 'color')), 
  'kegg_colormap_taxonomy.txt', 
  row.names = FALSE, quote = FALSE, col.names = FALSE
)

