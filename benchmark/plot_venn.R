library(VennDiagram)

venn.diagram.no.bg = function(x, filename, 
                              height = 3000, width = 3000, resolution = 500, 
                              imagetype = "png", units = "px", 
                              compression = "lzw", 
                              ...) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  grob.list = venn.diagram(x, filename = NULL, ...)
  
  if (!is.null(filename)) {
    current.type <- getOption("bitmapType")
    if (length(grep("Darwin", Sys.info()["sysname"]))) {
      options(bitmapType = "quartz")
    }
    else {
      options(bitmapType = "cairo")
    }
    if ("tiff" == imagetype) {
      tiff(filename = filename, height = height, width = width, 
           units = units, res = resolution, compression = compression)
    }
    else if ("png" == imagetype) {
      png(filename = filename, height = height, width = width, 
          units = units, res = resolution, bg = "transparent")
    }
    else if ("svg" == imagetype) {
      svg(filename = filename, height = height, width = width, bg = "transparent")
    }
    else {
      flog.error("You have misspelled your 'imagetype', please try again", 
                 name = "VennDiagramLogger")
      stop("You have misspelled your 'imagetype', please try again")
    }
    grid.draw(grob.list)
    dev.off()
    options(bitmapType = current.type)
    return(1)
  }
  return(grob.list)
}


plot.venn = function(reports, id = 'protein', 
                     filename = 'venn.svg', 
                     fill = c('#339dff', '#ff3355', '#65c3ba'), 
                     height = 2.5, width = 2.5, units = 'cm', margin = 0.1,
                     hide.text = FALSE,
                     ...) {
  id_list = lapply(reports, function(report) { 
    unique(report[[id]])
  })
  
  arg_list = list(
    id_list,
    fill = array(fill, dim = length(id_list)), 
    alpha = rep(0.25, length(id_list)),
    fontfamily = "sans", cat.fontfamily = "sans",
    col = 'white', imagetype = 'svg',
    filename = filename, 
    ext.text = FALSE,
    height = height, width = width, units = units, margin = margin,
    ...
  )
  
  if (hide.text) {
    arg_list = c(arg_list, list(cex = 0, cat.cex = 0))
  }
  
  do.call(venn.diagram.no.bg, args = arg_list)
}

