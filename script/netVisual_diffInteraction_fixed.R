# ============================================================
# Fixed netVisual_diffInteraction()
#
# Patches two bugs reported in sqjin/CellChat issue #708:
#
#  BUG 1 – "non-conformable arrays" (obj2 - obj1)
#    When the two CellChat objects contain different cell-type
#    labels, the net matrices have different dimensions and the
#    subtraction fails.
#    FIX: intersect row/col labels before subtracting so only
#    shared cell types are compared.
#
#  BUG 2 – igraph:::i_set_edge_attr length mismatch
#    "Length of new attribute value must be 1 or <N edges>, not <M>"
#    When self-loops (source == target) are present, the graph
#    has more edges than the non-loop `edge.start` rows, so
#    assigning loop.angle only to the looped subset crashes
#    because igraph expects a vector the length of ALL edges.
#    FIX: pre-initialise igraph::E(g)$loop.angle to a zero
#    vector of length nrow(edge.start) before writing the
#    loop-specific angles.
# ============================================================

netVisual_diffInteraction <- function(
    object,
    comparison     = c(1, 2),
    measure        = c("count", "weight", "count.merged", "weight.merged"),
    color.use      = NULL,
    color.edge     = c('#b2182b', '#2166ac'),
    title.name     = NULL,
    sources.use    = NULL,
    targets.use    = NULL,
    remove.isolate = FALSE,
    top            = 1,
    weight.scale   = FALSE,
    vertex.weight  = 20,
    vertex.weight.max  = NULL,
    vertex.size.max    = 15,
    vertex.label.cex   = 1,
    vertex.label.color = "black",
    edge.weight.max    = NULL,
    edge.width.max     = 8,
    alpha.edge         = 0.6,
    label.edge         = FALSE,
    edge.label.color   = 'black',
    edge.label.cex     = 0.8,
    edge.curved        = 0.2,
    shape              = 'circle',
    layout             = in_circle(),
    margin             = 0.2,
    arrow.width        = 1,
    arrow.size         = 0.2
) {
  options(warn = -1)
  measure <- match.arg(measure)

  obj1_raw <- object@net[[comparison[1]]][[measure]]
  obj2_raw <- object@net[[comparison[2]]][[measure]]

  # ── BUG 1 FIX ──────────────────────────────────────────────
  # Restrict to cell types present in BOTH datasets so matrix
  # dimensions match for the subtraction.
  shared_label <- intersect(rownames(obj1_raw), rownames(obj2_raw))
  obj1 <- obj1_raw[shared_label, shared_label]
  obj2 <- obj2_raw[shared_label, shared_label]
  # ───────────────────────────────────────────────────────────

  net.diff <- obj2 - obj1

  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name <- "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name <- "Differential interaction strength"
    }
  }

  color.use <- tryCatch(
    scPalette(length(shared_label)),
    error = function(e) scPalette(nrow(net.diff))
  )

  vertex.weight <- tryCatch(
    {
      if (length(vertex.weight) > 1) {
        # keep only shared labels if a per-cell vector was supplied
        if (!is.null(names(vertex.weight))) {
          vertex.weight[shared_label]
        } else {
          vertex.weight
        }
      } else {
        vertex.weight
      }
    },
    error = function(e) 20
  )

  # ── prepare the net matrix (filter by top / sources / targets) ──
  net <- net.diff
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      net[-sources.use, ] <- 0
    } else {
      net[!rownames(net) %in% sources.use, ] <- 0
    }
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      net[, -targets.use] <- 0
    } else {
      net[, !colnames(net) %in% targets.use] <- 0
    }
  }

  # keep only the top fraction of absolute weights
  net.abs <- abs(net)
  cutoff  <- stats::quantile(net.abs[net.abs > 0], probs = 1 - top,
                              na.rm = TRUE)
  if (length(cutoff) > 0) net[net.abs < cutoff] <- 0

  if (remove.isolate) {
    idx <- which(rowSums(abs(net)) == 0 & colSums(abs(net)) == 0)
    if (length(idx) > 0) {
      net       <- net[-idx, -idx]
      color.use <- color.use[-idx]
      if (length(vertex.weight) > 1) vertex.weight <- vertex.weight[-idx]
    }
  }

  # ── build igraph object ──────────────────────────────────────
  g <- igraph::graph_from_adjacency_matrix(
    net, mode = "directed", weighted = TRUE, diag = TRUE
  )

  edge.start <- igraph::ends(g, igraph::E(g), names = FALSE)

  # Resolve layout to a numeric matrix of vertex coordinates.
  # `layout` may arrive as: a pre-computed matrix, a layout function,
  # or the result of in_circle() which is itself a function.
  coords <- layout
  if (is.function(coords)) {
    coords <- coords(g)
  }
  # Final guard: if still not a 2-column numeric matrix, fall back.
  if (!is.matrix(coords) || ncol(coords) < 2) {
    coords <- igraph::layout_in_circle(g)
  }

  # vertex colours
  if (is.null(color.use)) color.use <- scPalette(igraph::vcount(g))
  igraph::V(g)$color       <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]

  # edge colours: red = increased, blue = decreased
  edge.color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(edge.color, alpha.f = alpha.edge)

  # edge widths
  edge.weight.abs <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) edge.weight.max <- max(edge.weight.abs)
  if (weight.scale) {
    igraph::E(g)$width <- 0.5 + edge.width.max * edge.weight.abs / edge.weight.max
  } else {
    igraph::E(g)$width <- 0.5 + edge.width.max * edge.weight.abs / max(edge.weight.abs)
  }

  # edge labels
  if (label.edge) {
    igraph::E(g)$label       <- round(igraph::E(g)$weight, 2)
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex   <- edge.label.cex
  }

  # ── BUG 2 FIX ──────────────────────────────────────────────
  # Self-loops need a loop.angle attribute.  igraph adds a
  # self-loop as an extra edge, making E(g) longer than
  # edge.start.  Pre-fill the whole edge attribute with 0 first,
  # then overwrite just the looped positions so the vector is
  # always exactly length(E(g)).
  loop.idx <- which(edge.start[, 1] == edge.start[, 2])
  if (length(loop.idx) > 0) {
    loop.angle <- base::ifelse(
      coords[edge.start[loop.idx, 1], 1] > 0, -pi / 2, pi / 2
    )
    # initialise attribute to 0 for ALL edges first
    igraph::E(g)$loop.angle                <- rep(0, igraph::ecount(g))
    igraph::E(g)$loop.angle[loop.idx]      <- loop.angle
  }
  # ───────────────────────────────────────────────────────────

  # vertex sizes
  if (is.null(vertex.weight.max)) vertex.weight.max <- max(vertex.weight)
  vertex.size <- vertex.size.max * vertex.weight / vertex.weight.max

  # ── plot ─────────────────────────────────────────────────────
  igraph::plot.igraph(
    g,
    layout             = coords,
    vertex.color       = igraph::V(g)$color,
    vertex.frame.color = igraph::V(g)$frame.color,
    vertex.size        = vertex.size,
    vertex.label.color = vertex.label.color,
    vertex.label.cex   = vertex.label.cex,
    edge.color         = igraph::E(g)$color,
    edge.width         = igraph::E(g)$width,
    edge.arrow.width   = arrow.width,
    edge.arrow.size    = arrow.size,
    edge.curved        = edge.curved,
    edge.loop.angle    = if (!is.null(igraph::edge_attr(g, "loop.angle")))
                           igraph::E(g)$loop.angle else 0,
    margin             = margin,
    vertex.shape       = shape,
    main               = title.name
  )

  invisible(g)
}
