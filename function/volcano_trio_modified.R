library(ggpubr)

volcano_plot = function (pvalues_df, comparison, p_cutoff = 0.05, fc_cutoff = 1, 
          label_rows = NULL, p_col_suffix = "pvalue", padj_col_suffix = "padj", 
          fc_col_suffix = "logFC", cutoff_criteria = "pvalue", label_col = "label", 
          label_size = 3, text_size = 10, marker_alpha = 0.7, marker_size = 3, 
          marker_outline_colour = "white", marker_outline_width = 0.5, 
          sig_names = NULL, colour_col = FALSE, 
          colours = c("salmon", "steelblue", "limegreen", "grey60"), drop_unused_cols = TRUE, 
          fc_line = TRUE, p_line = TRUE, line_colours = c("black", 
                                                          "black")) 
{
  if (colour_col & !"col" %in% colnames(pvalues_df)) {
    stop("pvalues_df must have a 'col' column for the polar colour scheme")
  }
  if (colour_col) {
    toptable <- pvalues_df[, c(colnames(pvalues_df)[grepl(comparison, 
                                                          colnames(pvalues_df))], label_col, "col")]
  }
  else {
    toptable <- pvalues_df[, c(colnames(pvalues_df)[grepl(comparison, 
                                                          colnames(pvalues_df))], label_col)]
  }
  colnames(toptable)[colnames(toptable) == label_col] <- "label"
  if (!length(colours) %in% c(0, 3, 4, 7, 8)) {
    stop(paste("colours must be of length 3, 4, 7, or 8 according to", 
               "colour-coding. See `?volcano_plot` for more info"))
  }
  if (is.null(label_col)) {
    toptable$label <- rownames(toptable)
  }
  use_cols <- c(pvalue = p_col_suffix, padj = padj_col_suffix, 
                logFC = fc_col_suffix)
  colnames(toptable)[unlist(lapply(use_cols, function(x) {
    which(grepl(x, colnames(toptable)))
  }))] <- names(use_cols)
  if (!is.null(sig_names)) {
    if (length(sig_names) != length(colours)) {
      stop("length(sig_names) must equal length(colours)")
    }
  }
  else {
    if (length(colours) == 3) {
      sig_names <- c(paste(cutoff_criteria, "<=", p_cutoff, 
                           "and log(FC) <", -1 * fc_cutoff), paste(cutoff_criteria, 
                                                                   "<=", p_cutoff, "and log(FC) >=", fc_cutoff), 
                     "Not Significant")
    }
    else if (length(colours) == 4) {
      sig_names <- c(paste(cutoff_criteria, "<=", p_cutoff, 
                           "and |log(FC)| >=", fc_cutoff), paste(cutoff_criteria, 
                                                                 "<=", p_cutoff, "and |log(FC)| <", fc_cutoff), 
                     paste(cutoff_criteria, ">", p_cutoff, "and |log(FC)| >=", 
                           fc_cutoff), paste(cutoff_criteria, ">", p_cutoff, 
                                             "and |log(FC)| <", fc_cutoff))
    }
    else {
      sig_names <- c(paste(cutoff_criteria, "<=", p_cutoff, 
                           "and log(FC) >=", fc_cutoff), paste(cutoff_criteria, 
                                                               "<=", p_cutoff, "and", -1 * fc_cutoff, "< log(FC) <= 0"), 
                     paste(cutoff_criteria, "<=", p_cutoff, "and 0 < log(FC) <", 
                           fc_cutoff), paste(cutoff_criteria, "<=", p_cutoff, 
                                             "and log(FC) <=", -1 * fc_cutoff), paste(cutoff_criteria, 
                                                                                      ">", p_cutoff, "and log(FC) >=", fc_cutoff), 
                     paste(cutoff_criteria, ">", p_cutoff, "and", 
                           -1 * fc_cutoff, "< log(FC) <= 0"), paste(cutoff_criteria, 
                                                                    ">", p_cutoff, "and 0 < log(FC) <", fc_cutoff), 
                     paste(cutoff_criteria, ">", p_cutoff, "and log(FC) <=", 
                           -1 * fc_cutoff))
    }
  }
  toptable$cols <- NA
  toptable$cols[toptable[, cutoff_criteria] <= p_cutoff & 
                  toptable$logFC <= -1 * fc_cutoff] <- "a"
  toptable$cols[toptable[, cutoff_criteria] <= p_cutoff & 
                  toptable$logFC > -1 * fc_cutoff & toptable$logFC <= 
                  0] <- "b"
  toptable$cols[toptable[, cutoff_criteria] <= p_cutoff & 
                  toptable$logFC < fc_cutoff & toptable$logFC > 0] <- "c"
  toptable$cols[toptable[, cutoff_criteria] <= p_cutoff & 
                  toptable$logFC >= fc_cutoff] <- "d"
  toptable$cols[toptable[, cutoff_criteria] > p_cutoff & toptable$logFC <= 
                  -1 * fc_cutoff] <- "e"
  toptable$cols[toptable[, cutoff_criteria] > p_cutoff & toptable$logFC > 
                  -1 * fc_cutoff & toptable$logFC <= 0] <- "f"
  toptable$cols[toptable[, cutoff_criteria] > p_cutoff & toptable$logFC < 
                  fc_cutoff & toptable$logFC > 0] <- "g"
  toptable$cols[toptable[, cutoff_criteria] > p_cutoff & toptable$logFC >= 
                  fc_cutoff] <- "h"
  if (length(colours) == 3) {
    toptable$cols[toptable$cols %in% letters[c(2:3, 5:8)]] <- "c"
    toptable$cols[toptable$cols %in% "d"] <- "b"
  }
  if (length(colours) == 4) {
    toptable$cols[toptable$cols %in% "c"] <- "b"
    toptable$cols[toptable$cols %in% "d"] <- "a"
    toptable$cols[toptable$cols %in% c("e", "h")] <- "c"
    toptable$cols[toptable$cols %in% c("f", "g")] <- "d"
  }
  toptable$cols <- factor(toptable$cols, labels = sig_names[letters[seq_len(length(colours))] %in% 
                                                              unique(toptable$cols)])
  if (colour_col) {
    mapping <- setNames(levels(toptable$col), colours)
    toptable$cols <- toptable$col
    sig_names <- levels(toptable$cols)
  }
  else {
    mapping <- setNames(sig_names, colours)
  }
  toptable <- toptable[!is.na(toptable$logFC), ]
  toptable <- toptable[!is.na(-log10(toptable$pvalue)), ]
  if (cutoff_criteria == "padj") {
    toptable <- toptable[!is.na(-log10(toptable$padj)), 
    ]
  }
  if (cutoff_criteria == "padj") {
    toptable$lp <- -log10(toptable$padj)
  }
  else {
    toptable$lp <- -log10(toptable$pvalue)
  }
  p <- ggplot(toptable, aes_string(x = "logFC", y = "lp")) + 
    geom_point(data = toptable, size = marker_size, alpha = marker_alpha, 
               color = marker_outline_colour, stroke = marker_outline_width, 
               aes_string(fill = "cols"), shape = 21) + 
    scale_fill_manual(limits = sig_names, 
                      values = colours, drop = drop_unused_cols) + 
    labs(y = if (cutoff_criteria == "padj") expression(-log["10"] * `Adjusted p-value`) else expression(-log["10"] * p), 
         x = expression(log["2"] * FC), 
         title = gsub("_", " vs ", comparison), 
         fill = "Significance") + 
    theme_classic() + 
    theme(text = element_text(size = text_size), 
          legend.background = element_rect(fill = "transparent", colour = NA), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.background = element_rect(fill = "transparent", colour = NA))
  if (fc_line) {
    p <- p + geom_vline(xintercept = fc_cutoff, linetype = "dashed", 
                        color = line_colours[1], size = 0.5) + geom_vline(xintercept = -fc_cutoff, 
                                                                          linetype = "dashed", color = line_colours[1], size = 0.5)
  }
  if (any(toptable[, cutoff_criteria] <= p_cutoff) & p_line == 
      TRUE) {
    top <- max(toptable$padj[toptable[, cutoff_criteria] <= 
                               p_cutoff], na.rm = TRUE)
    p <- p + geom_hline(yintercept = -log10(top), linetype = "dashed", 
                        color = line_colours[2], size = 0.5)
  }
  if (!is.null(label_rows)) {
    if (!all(is.numeric(label_rows))) {
      if (!all(label_rows %in% rownames(toptable))) {
        stop("label_rows must be in rownames(pvalues)")
      }
    }
    if (all(is.numeric(label_rows))) {
      if (!all(label_rows < nrow(toptable))) {
        stop("label_rows not in 1:nrow(pvalues)")
      }
    }
    label_df <- toptable[label_rows, ]
    if (cutoff_criteria == "padj") {
      label_df$lp <- -log10(label_df$padj)
    }
    else {
      label_df$lp <- -log10(label_df$pvalue)
    }
    p <- p + geom_text_repel(data = label_df, aes_string(x = "logFC", 
                                                         y = "lp", label = "label"), size = label_size, box.padding = unit(0.35, 
                                                                                                                           "lines"), point.padding = unit(0.3, "lines"))
  }
  return(p)
}

volcano_trio = function (polar, p_cutoff = 0.05, cutoff_criteria = "pvalue", 
                         fc_cutoff = 1, label_rows = NULL, label_size = 3, text_size = 10, 
                         marker_alpha = 0.7, marker_size = 3, marker_outline_colour = "white", 
                         marker_outline_width = 0.5, shared_legend_size = 1, sig_names = NULL, 
                         colour_scheme = "none", colours = c("salmon", "steelblue", 
                                                             "limegreen", "grey60"), fc_line = TRUE, p_line = TRUE, 
                         line_colours = c("black", "black"), share_axes = TRUE) 
{
  if (class(polar) != "polar") 
    stop("polar must be an object of class polar")
  pvalues <- polar@pvalues
  pvalues$col <- polar@polar$sig
  if (!colour_scheme %in% c("none", "polar", "upregulated")) {
    stop("colour_scheme must be in c(\"none\", \"polar\", \"upregulated\")")
  }
  if (colour_scheme %in% c("polar", "upregulated") & length(colours) != 
      7) {
    stop(paste("For polar or upregulated colour schemes, colours must be of", 
               "length 7 to match c(levels(polar@sampledata[,", 
               "polar@contrast]), polar@non_sig_name)"))
  }
  groups <- levels(polar@sampledata[, polar@contrast])
  sig_groups <- c(paste0(groups[1], "+"), paste0(groups[1], 
                                                 "+", groups[2], "+"), paste0(groups[2], "+"), paste0(groups[2], 
                                                                                                      "+", groups[3], "+"), paste0(groups[3], "+"), paste0(groups[1], 
                                                                                                                                                           "+", groups[3], "+"))
  pvalues$col <- factor(pvalues$col, levels = c(sig_groups, 
                                                polar@non_sig_name))
  if (!length(colours) %in% c(0, 3, 4, 7, 8)) {
    stop(paste("colours must be of length 3, 4, 7 or 8 according to", 
               "colour-coding. See `?volcano_trio` for more info"))
  }
  if (!any(grepl(cutoff_criteria, colnames(pvalues)))) {
    stop("There must be cutoff_criteria columns in pvalues")
  }
  if (any(unlist(lapply(colours, function(x) {
    class(try(col2rgb(x), silent = TRUE)) == "try-error"
  })))) {
    stop(paste(paste(colours[unlist(lapply(colours, function(x) {
      class(try(col2rgb(x), silent = TRUE)) == "try-error"
    }))], collapse = ", "), "is not a valid colour"))
  }
  if (length(line_colours) != 2) {
    stop(paste("The line_colour vector must be of length 2 for fc and p", 
               "cutoff lines respectively"))
  }
  if (any(unlist(lapply(line_colours, function(x) {
    class(try(col2rgb(x), silent = TRUE)) == "try-error"
  })))) {
    stop(paste(paste(line_colours[unlist(lapply(line_colours, function(x) {
      class(try(col2rgb(x), silent = TRUE)) == "try-error"
    }))], collapse = ", "), "is not a valid colour"))
  }
  if (!is.logical(fc_line)) {
    stop(paste("fc_line must be a logical dictating whether to add vertical", 
               "lines for fold change cutoff"))
  }
  if (!is.logical(p_line)) {
    stop(paste("p_line must be a logical dictating whether to add horizontal", 
               "lines for p-value cutoff"))
  }
  if (!is.logical(share_axes)) {
    stop(paste("share_axes must be a logical dictating whether combined plots", 
               "should share axes"))
  }
  v <- c("p_cutoff", "label_size", "fc_cutoff", "text_size", 
         "marker_size", "shared_legend_size")
  if (any(unlist(lapply(v, function(x) !is.numeric(get(x)))))) {
    stop(paste(paste(v[unlist(lapply(v, function(x) !is.numeric(get(x))))], 
                     collapse = ", "), "must be numeric"))
  }
  if (!(0 <= p_cutoff & p_cutoff <= 1)) 
    stop("p_cutoff must be between 0 and 1")
  v <- c("label_size", "fc_cutoff", "text_size", "marker_size", 
         "shared_legend_size")
  if (any(unlist(lapply(v, function(x) get(x) <= 0)))) {
    stop(paste(v[unlist(lapply(v, function(x) get(x) <= 
                                 0))], "must be greater than 0"))
  }
  groups <- gsub("_pvalue", "", colnames(pvalues)[grepl("pvalue", 
                                                        colnames(pvalues))])
  if (!is.null(polar@multi_group_test)) {
    groups <- groups[groups != polar@multi_group_test]
  }
  if (!all(paste0(groups, "_logFC") %in% colnames(pvalues))) {
    stop(paste("No", paste(groups[!paste0(groups, "_logFC") %in% colnames(pvalues)]
                           , collapse = ", "), 
               "fold-change columns"))
  }
  plot_outputs <- lapply(groups, function(x) {
    if (colour_scheme == "upregulated") {
      options <- levels(polar@sampledata[, polar@contrast])
      colours_use <- colours[c(1, 3, 5, 7)]
      colours_use <- colours_use[c(rev(match(unlist(strsplit(x, 
                                                             "_")), options)), 4)]
      sig_names <- c(paste("Up in", rev(unlist(strsplit(x, 
                                                        "_")))), polar@non_sig_name)
    }
    else {
      colours_use <- colours
    }
    if (cutoff_criteria == "padj") 
      pa <- "padj"
    else pa <- NULL
    volcano_plot(pvalues_df = pvalues, comparison = x, p_cutoff = p_cutoff, 
                 padj_col_suffix = pa, fc_cutoff = fc_cutoff, label_rows = label_rows, 
                 label_size = label_size, text_size = text_size, 
                 marker_size = marker_size, marker_alpha = marker_alpha, 
                 marker_outline_colour = marker_outline_colour, marker_outline_width = marker_outline_width, 
                 sig_names = sig_names, colour_col = colour_scheme %in% 
                   c("polar"), colours = colours_use, fc_line = fc_line, 
                 p_line = p_line, drop_unused_cols = FALSE, line_colours = line_colours, 
                 cutoff_criteria = cutoff_criteria)
  })
  names(plot_outputs) <- groups
  out <- plot_outputs
  if (colour_scheme != "upregulated") {
    out[[1]] <- out[[1]] + theme(legend.key.size = unit(shared_legend_size, "cm"), 
                                 legend.key.width = unit(shared_legend_size, "cm"))
  }
  if (share_axes) {
    yRange <- unlist(lapply(out, function(x) layer_scales(x)$y$range$range))
    xRange <- unlist(lapply(out, function(x) layer_scales(x)$x$range$range))
    yMax <- max(yRange, na.rm = TRUE)
    xMax <- max(abs(xRange), na.rm = TRUE)
    out[1:3] <- lapply(out[1:3], function(p) {
      p + lims(x = c(-1 * xMax, xMax), y = c(0, yMax))
    })
  }
  plot_outputs[["All"]] <- ggarrange(plotlist = out, ncol = length(out), 
                                     nrow = 1, legend = "right", common.legend = colour_scheme != 
                                       "upregulated")
  return(plot_outputs)
}

