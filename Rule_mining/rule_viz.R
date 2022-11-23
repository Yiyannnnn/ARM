#rule viz

.get_parameters <- function(parameter, defaults) {
  
  defaults <- as.list(defaults)
  parameter <- as.list(parameter)
  
  ## add verbose
  if (is.null(defaults$verbose))
    defaults$verbose <- FALSE
  
  if (length(parameter) != 0) {
    o <- pmatch(names(parameter), names(defaults))
    defaults[o[!is.na(o)]] <- parameter[!is.na(o)]
  }
  if (defaults$verbose) {
    cat("Used control parameters:\n")
    #print(defaults)
    cat(rbind(names(defaults), " = ", gsub("\n", " ", as.character(defaults))),
        sep = c("\t", " ", "\n"))
  }
  
  defaults
}

.col_picker <- function(level, palette, alpha = NULL) {
  col <- palette[floor(level * (length(palette) - 1)) + 1]
  if (!is.null(alpha)) {
    col <- apply(sapply(col, grDevices::col2rgb) / 255, 2,
                 function(x)
                   grDevices::rgb(x[1], x[2], x[3], alpha = alpha))
  }
  col
}

default_colors <- function(n , alpha = 1,type = "taxa"){
  if (type == "rule") {
    grDevices::colorRampPalette(c("#668DC0", "#C0D0EF", "#D9E0E6"), alpha = alpha)(n)
  }
  else{
    colorspace::rainbow_hcl(n,alpha = alpha)
    #grDevices::colorRampPalette(c("#FD5993", "#FFC3D8", "#FFFFFF"), alpha = alpha)(n)
  }
}


map <- function(x, range = c(0, 1), from.range = NA) {
  if (is.null(from.range) || any(is.na(from.range)))
    from.range <- range(x, na.rm = TRUE)
  
  ## check if all values are the same
  if (!diff(from.range)) {
    if (is.matrix(x))
      return(matrix(
        mean(range),
        ncol = ncol(x),
        nrow = nrow(x),
        dimnames = dimnames(x)
      ))
    else
      return(structure(rep(mean(range), length(x)), names = names(x)))
    
  }
  
  ## map to [0,1]
  x <- (x - from.range[1])
  x <- x / diff(from.range)
  ## handle single values
  if (diff(from.range) == 0)
    x <- 0
  
  ## map from [0,1] to [range]
  if (range[1] > range[2])
    x <- 1 - x
  x <- x * (abs(diff(range))) + min(range)
  
  x[x < min(range) | x > max(range)] <- NA
  
  x
}

limit <- function(x, max, shading = NULL, measure = NULL, quiet = FALSE) {
  if (is.null(max)) return(x)
  if (length(x) <= max) return(x)
  
  l <- if (!is.null(shading)) shading[1] else measure[1]
  
  if (!quiet)
    warning(
      "Too many ", class(x), " supplied. Only plotting the best ",
      max, " using ",
      if (!is.null(l)) sQuote(l) else "the original order",
      " (change control parameter max if needed).",
      call. = FALSE
    )
  
  head(x,
       n = max,
       by = l,
       decreasing = TRUE)
}

graph_visNetwork <- function(x, measure = "support", shading = "lift", control = NULL, ...) {
  
  control <- c(control, list(...))
  control <- .get_parameters(
    control,
    list(
      #main = paste("Graph for", length(x), "rules"),
      itemCol = grDevices::hcl(h = 0),
      nodeColrule = default_colors(100, alpha = 0,type = "rule"),
      nodeColtaxa = default_colors(100, alpha = 0),
      precision = 3,
      igraphLayout = "layout_nicely",
      interactive = TRUE,
      engine = "visNetwork",
      max = 100,
      selection_menu = TRUE,
      degree_highlight = 1
    )
  )
  
  x <- limit(x, control$max, shading, measure)
  g <- associations2igraph(x)
  
  va <- igraph::get.vertex.attribute(g)
  n_items <- sum(va$type == 1)
  
  titleAssoc <- paste0('<B>[',1:length(x),']</B><BR>',
                       labels(x, itemSep = ',<BR>&nbsp;&nbsp;', ruleSep = '<BR>&nbsp;&nbsp; => ', 
                              setStart = '<B>{', setEnd = '}</B>'),"<BR><BR>",
                       apply(quality(x), MARGIN = 1,
                             FUN = function(x) paste(names(x), "=", signif(x, control$precision), collapse = "<BR>")))
  
  title <- c(va$label[va$type == 1], titleAssoc)
  
  size <- map(present_proportion[va$label][va$type == 1], c(10, 20))
  size[which(va$type == 2)] <- 0.25
  
  phylum <- 
    factor(stool_phyloseq@tax_table[which(rownames(stool_phyloseq@tax_table) %in% va$label),"phylum"])
  family <- 
    factor(stool_phyloseq@tax_table[which(rownames(stool_phyloseq@tax_table) %in% va$label),"family"])
  genus <- 
    factor(stool_phyloseq@tax_table[which(rownames(stool_phyloseq@tax_table) %in% va$label),"genus"])
  color <- c(c("olivedrab","lightsalmon","paleturquoise","plum","lightcoral","khaki","beige")[as.numeric(genus)],
             #.col_picker(map(present_proportion[va$label][va$type == 1], c(0.9, 0.1)), control$nodeColtaxa)
             .col_picker(map(va[[shading]][va$type == 2], c(0.9, 0.1)), control$nodeColrule))
  #rep(control$itemCol[1], n_items),
  
  
  label <- str_to_title(gsub('\\b(\\pL)\\pL{2,}|.','\\L\\1',va$label,perl = TRUE))
  label[va$type == 2] <- NA
  #substr(class(x), 1, nchar(class(x)) - 1L)#seq_along(x))
  
  nodes <- data.frame(
    id = seq_along(va$name),
    label = label,
    group = va$type,
    #value = size,
    color = color,
    title = title,
    shape = ifelse(va$type == 2, "box", "circle"),
    phylum = c(phylum,1:length(x))
  )
  
  e <- igraph::as_data_frame(g, what = "edges")
  edges <- data.frame(
    from = as.integer(factor(e$from, levels = va$name)),
    to = as.integer(factor(e$to, levels = va$name)),
    arrows = "to", arrow.size=0.01
  )
  
  lnodes <- data.frame(label = c("Faecali- \n bacterium", "Fusicateni- \n bacter","Anaero- \n stipes","Dorea","Anaero- \n butyricum","Blautia","Eubacterium"),
                       shape = c("circle"), color = c("lightcoral", "khaki","lightsalmon","plum","olivedrab","paleturquoise",NA),
                       title = "Genus",icon = list(size = 10))
  
  
  visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
    visNetwork::visIgraphLayout(layout = control$igraphLayout) %>%
    visNetwork::visOptions(
      highlightNearest =
        list(
          enabled = TRUE,
          degree = control$degree_highlight,
          hover = TRUE
        ),
      nodesIdSelection = control$selection_menu
    ) %>% visNetwork :: visNodes(scaling = list(min = 1, max = 2))%>%
    #visGroups(groupname = "Faecalibacterium", color = "lightcoral",icon = list(size = 10)) %>%
    #visGroups(groupname = "Fusicatenibacter", color = "khaki",icon = list(size = 10)) %>%
    #visGroups(groupname = "Anaerostipes", color = "lightsalmon", icon = list(size = 10)) %>%
    #visGroups(groupname = "Dorea", color = "plum", icon = list(size = 10)) %>%
    #visGroups(groupname = "Anaerobutyricum", color = "olivedrab", icon = list(size = 10)) %>%
    #visGroups(groupname = "Blautia", color = "paleturquoise", icon = list(size = 10)) %>%
    #visGroups(groupname = "Eubacterium", color = NA, icon = list(size = 10)) %>%
    addFontAwesome() %>%
    visLegend(addNodes = lnodes, useGroups = FALSE,width = 0.5, position = "right", main = "Genus")
}



