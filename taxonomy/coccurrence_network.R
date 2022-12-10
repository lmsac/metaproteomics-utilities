library(igraph)

build.cooccurrence.network = function(
  taxonomy_quantity, 
  run_groups,
  level = 'genus',
  min_percentage = 0.05,
  r_threshold = 0.8,
  p_threshold = 0.05,
  filename = 'cooccurrence_network.graphml',
  return_data = FALSE
) {
  data = cbind(
    taxonomy_quantity[, c('level', 'name')],
    do.call(cbind, lapply(run_groups, function(group) {
      taxonomy_quantity[, grep(group, colnames(taxonomy_quantity))]
    }))
  )
  
  data = do.call(cbind, lapply(run_groups, function(group) {
    data_level = data[data$level == level & data$name != '-', ]
    
    data_level = subset(data_level, apply(data_level[, -1:-2], 1, function(x) {
      percentage = x / colSums(data[match('root', data$level), -1:-2])
      max(percentage) >= min_percentage * 0.01
    }))
    
    data = as.matrix(data_level[, -1:-2])
    rownames(data) = data_level$name
    data
  }))
  
  correlations = Hmisc::rcorr(t(data), type = 'spearman')
  correlations$P.adj = matrix(
    p.adjust(correlations$P, method = 'bonferroni'),
    nrow = nrow(correlations$P),
    ncol = ncol(correlations$P),
    dimnames = dimnames(correlations$P)
  )
  
  r_filtered = correlations$r
  r_filtered[abs(correlations$r) < r_threshold | 
               correlations$P.adj > p_threshold] = 0
  r_filtered = r_filtered[rowSums(r_filtered) > 1, colSums(r_filtered) > 1]
  
  correlation_graph = graph.adjacency(
    r_filtered, 
    mode = 'undirected',
    weighted = TRUE
  )
  correlation_graph = simplify(correlation_graph)
  V(correlation_graph)$label = V(correlation_graph)$name
  V(correlation_graph)$degree = degree(correlation_graph)
  
  cluster = cluster_walktrap(correlation_graph)
  V(correlation_graph)$cluster = cluster$membership + 1
  
  if (!is.null(filename)) {
    write_graph(correlation_graph, filename, format = 'graphml')
    if (return_data) {
      data
    }
  }
  else {
    if (return_data) {
      list(graph = correlation_graph, data = data)
    }
    else {
      correlation_graph
    }
  }
}



library(readr)

taxonomy_quantity = read_csv('taxonomy_quantity.csv')


build.cooccurrence.network(
  taxonomy_quantity, 
  run_groups = c(
    'W Sp' = '^value1',
    'W Su' = '^value2',
    'W Au' = '^value3'
  ),
  filename = 'cooccurrence_network.graphml'
)

# plot.igraph(
#   correlation_graph,
#   vertex.size = 4,
#   vertex.color = V(correlation_graph)$cluster,
#   vertex.label.cex = 0.75,
#   edge.curved = TRUE,
#   edge.size = 1.5,
#   layout = layout.fruchterman.reingold
# )


