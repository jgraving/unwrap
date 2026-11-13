# Bats (Goldshtein): figure construction.
#
# Sourced by notebooks_r/goldshtein_bats.ipynb (no standalone output). R
# twin of notebooks/goldshtein_bats_plots.py: the notebook keeps the data,
# both fits, the diagnostic values and the post-hoc tables; only the
# figure construction lives here (plus the basemap tile helpers it
# needs).

CARTO_SUBDOMAINS <- c("a", "b", "c", "d")
R_EARTH <- 20037508.342789244
TILE_PX <- 256L
TILE_HANDLE <- curl::new_handle()
curl::handle_setheaders(TILE_HANDLE,
  "User-Agent" = "unwrap-r-notebooks-tiles (bake; jgraving/unwrap)")

trace_plot <- function(fit, parameters, angular = character()) {
  d <- draw_matrix(fit, parameters)
  long <- do.call(rbind, lapply(parameters, function(p) {
    x <- d[[p]]
    if (p %in% angular) x <- wrap(x - circ_mean(x))
    data.frame(parameter = p, value = x, iteration = d[[".iteration"]],
               chain = factor(d[[".chain"]]))
  }))
  ggplot(long, aes(iteration, value, colour = chain, group = chain)) +
    geom_line(linewidth = 0.25, alpha = 0.6) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(x = "retained iteration", y = "draw (angular traces centered)") +
    theme(legend.position = "bottom")
}

# The Python _format_treatment_name rule: lowercase, split on commas,
# drop "control" when it is one of several parts, sort the parts,
# capitalise each part's first letter, rejoin with ", ". A coefficient
# figure additionally names its intercept row "Control (Intercept)".
fmt_treat_label <- function(name) {
  s <- tolower(trimws(as.character(name)))
  parts <- trimws(unlist(strsplit(s, ",", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) > 1L && "control" %in% parts) parts <- parts[parts != "control"]
  if (!length(parts)) parts <- "control"
  parts <- sort(parts)
  parts <- paste0(toupper(substr(parts, 1, 1L)), substr(parts, 2, 1000L))
  paste(parts, collapse = ", ")
}

fetch_carto_tile <- function(z, x, y, variant) {
  cache_dir <- file.path(
    Sys.getenv("UNWRAP_TILE_CACHE",
      file.path(Sys.getenv("XDG_CACHE_HOME",
        file.path(Sys.getenv("HOME", "."), ".cache")), "unwrap-carto-tiles")),
    variant)
  f <- file.path(cache_dir, sprintf("%d_%d_%d.png", z, x, y))
  if (file.exists(f)) {
    img <- try(png::readPNG(f), silent = TRUE)
    if (!inherits(img, "try-error")) return(img)
  }
  url <- sprintf("https://%s.basemaps.cartocdn.com/%s/%d/%d/%d.png",
    CARTO_SUBDOMAINS[((x + y) %% 4L) + 1L], variant, z, x, y)
  if (nzchar(CARTO_KEY)) url <- paste0(url, "?key=", CARTO_KEY)
  resp <- curl::curl_fetch_memory(url, handle = TILE_HANDLE)
  # CartoDB answers 200 with a small watermarked placeholder for a missing
  # or invalid key; treat that as a failure too.
  if (resp$status_code != 200L || length(resp$content) < 5000L) return(NULL)
  img <- png::readPNG(resp$content)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  writePNG(img, f)
  img
}

# Compose the tiles covering [xmin, xmax] x [ymin, ymax] (EPSG:3857,
# meters) at the given zoom into one raster. The extent is quantized to
# the tile-pixel grid first, so image pixels line up 1:1 with tile pixels
# and the copy is exact. Returns list(img, xmin, ymin, xmax, ymax) where
# the limits are the (quantized) range the raster spans.
basemap_panel_img <- function(xmin, ymin, xmax, ymax, variant, zoom,
                              fallback_hex, warn = TRUE) {
  if (!nzchar(CARTO_KEY)) {
    if (warn) cat("Basemap: flat-panel fallback (no UNWRAP_CARTO_API_KEY set).\n")
    return(list(flat = TRUE, hex = fallback_hex, xmin = xmin, ymin = ymin,
      xmax = xmax, ymax = ymax))
  }
  n <- 2 ^ zoom
  ts <- 2 * R_EARTH / n
  px_m <- ts / TILE_PX
  U0 <- R_EARTH - ymax       # u measures distance south of the north pole
  U1 <- R_EARTH - ymin
  U0 <- floor(U0 / px_m) * px_m
  U1 <- ceiling(U1 / px_m) * px_m
  xmin <- floor(xmin / px_m) * px_m
  xmax <- ceiling(xmax / px_m) * px_m
  # The panel draws the tiles as a fine vector mosaic (~64 cells across the
  # long edge): each fetched tile contributes the average color of every
  # mosaic cell it covers, so the full-resolution tile grid never has to
  # exist as one big array.
  n_mosaic <- 64L
  cell_m <- max(xmax - xmin, (R_EARTH - U0) - (R_EARTH - U1)) / n_mosaic
  mx0 <- floor(xmin / cell_m) * cell_m
  my0 <- floor((R_EARTH - U1) / cell_m) * cell_m
  mx0i <- mx0 / cell_m
  my0i <- my0 / cell_m
  mnx <- as.integer(ceiling((xmax - mx0) / cell_m))
  mny <- as.integer(ceiling((R_EARTH - U0 - my0) / cell_m))
  acc <- array(0, c(mny, mnx, 3L))
  cnt <- matrix(0L, nrow = mny, ncol = mnx)
  tx <- floor(c(xmin, xmax) / ts)
  ty <- floor(c(U0, U1) / ts)
  failed <- FALSE
  for (txi in as.integer(tx[1]):as.integer(tx[2])) {
    for (tyi in as.integer(ty[1]):as.integer(ty[2])) {
      tile <- fetch_carto_tile(zoom, txi, tyi, variant)
      if (is.null(tile)) { failed <- TRUE; next }
      if (dim(tile)[3] == 4L) tile <- tile[, , 1:3]
      tlx <- txi * ts
      tby <- R_EARTH - (tyi + 1) * ts
      cxa <- as.integer(floor(pmax(xmin, tlx) / cell_m))
      cxb <- as.integer(ceiling(pmin(xmax, tlx + ts) / cell_m))
      cya <- as.integer(floor(pmax(R_EARTH - U1, tby) / cell_m))
      cyb <- as.integer(ceiling(pmin(R_EARTH - U0, tby + ts) / cell_m))
      for (cxi in cxa:cxb) {
        ix <- as.integer(cxi - mx0i) + 1L
        if (ix < 1L || ix > mnx) next
        for (cyi in cya:cyb) {
          iy <- as.integer(cyi - my0i) + 1L
          if (iy < 1L || iy > mny) next
          p0 <- pmax(1L, as.integer(floor((cxi * cell_m - tlx) / px_m)) + 1L)
          p1 <- pmin(TILE_PX, as.integer(ceiling(((cxi + 1) * cell_m - tlx) / px_m)))
          q0 <- pmax(1L, as.integer(floor((cyi * cell_m - tby) / px_m)) + 1L)
          q1 <- pmin(TILE_PX, as.integer(ceiling(((cyi + 1) * cell_m - tby) / px_m)))
          if (p1 < p0 || q1 < q0) next
          blk <- tile[q0:q1, p0:p1, ]
          acc[iy, ix, ] <- acc[iy, ix, ] + colSums(blk)
          cnt[iy, ix] <- cnt[iy, ix] + prod(dim(blk)[1:2])
        }
      }
    }
  }
  if (failed || !any(cnt > 0L)) {
    if (warn) cat("Basemap failed (continuing without tiles): some CartoDB tiles were unavailable.\n")
    return(list(flat = TRUE, hex = fallback_hex, xmin = xmin,
      ymin = R_EARTH - U1, xmax = xmax, ymax = R_EARTH - U0))
  }
  have <- which(cnt > 0L, arr.ind = TRUE)
  tdf <- data.frame(
    x = (mx0i + have[, 2] - 1L + 0.5) * cell_m,
    y = (my0i + have[, 1] - 1L + 0.5) * cell_m,
    col = grDevices::rgb(
      acc[cbind(have[, 1], have[, 2], rep(1L, nrow(have)))] / cnt[have],
      acc[cbind(have[, 1], have[, 2], rep(2L, nrow(have)))] / cnt[have],
      acc[cbind(have[, 1], have[, 2], rep(3L, nrow(have)))] / cnt[have]))
  list(flat = FALSE, tiles = tdf, cell = cell_m, xmin = xmin,
    ymin = R_EARTH - U1, xmax = xmax, ymax = R_EARTH - U0)
}

# Consecutive-point track segments in 3857 (solid steps + dashed gap
# connectors), the R twin of the Python segment/gap layer. With pcol the
# per-segment probability is the mean of the two endpoint values.
segments_3857 <- function(track, pcol = NULL) {
  trk_xy <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(track, coords = c("lon", "lat"), crs = 4326, remove = FALSE), 3857))
  track$px <- trk_xy[, 1]
  track$py <- trk_xy[, 2]
  do.call(rbind, lapply(split(track, track$tag), function(g) {
    g <- g[order(g$time), ]
    if (nrow(g) < 2L) return(NULL)
    n <- nrow(g)
    out <- data.frame(
      x = g$px[-n], y = g$py[-n], xend = g$px[-1], yend = g$py[-1],
      treatment = factor(g$treatment[-1], levels = levels(df_track$treatment)),
      destination = g$dest[-1], gap = as.integer(g$gap_before[-1]) == 1L,
      time = g$time[-1], row.names = NULL)
    if (!is.null(pcol)) out$p <- 0.5 * (g[[pcol]][-n] + g[[pcol]][-1])
    out
  }))
}

# The unified data extent the Python figures compute (min/max over all
# layers + a small pad).
map_extent <- function(layer_bounds, pad = 0.06) {
  minx <- min(vapply(layer_bounds, function(b) b["xmin"], numeric(1)))
  miny <- min(vapply(layer_bounds, function(b) b["ymin"], numeric(1)))
  maxx <- max(vapply(layer_bounds, function(b) b["xmax"], numeric(1)))
  maxy <- max(vapply(layer_bounds, function(b) b["ymax"], numeric(1)))
  c(xmin = minx - (maxx - minx) * pad, ymin = miny - (maxy - miny) * pad,
    xmax = maxx + (maxx - minx) * pad, ymax = maxy + (maxy - miny) * pad)
}

seg_bounds <- function(seg)
  c(xmin = min(seg$x, seg$xend), ymin = min(seg$y, seg$yend),
    xmax = max(seg$x, seg$xend), ymax = max(seg$y, seg$yend))

# The Python panels use set_aspect("equal", adjustable = "box"): the axis
# box is forced to the panel aspect and the data is centered inside it.
# Emulate that by padding the data extent to the panel aspect (the extra
# space shows the white panel background, as in the references).
aspect_extent <- function(bm, aspect, xlim = NULL, ylim = NULL) {
  if (is.null(xlim)) xlim <- c(bm$xmin, bm$xmax)
  if (is.null(ylim)) ylim <- c(bm$ymin, bm$ymax)
  ew <- diff(xlim); eh <- diff(ylim)
  if (ew / eh > aspect) ew <- eh * aspect else eh <- ew / aspect
  list(xlim = c(mean(xlim) - ew / 2, mean(xlim) + ew / 2),
    ylim = c(mean(ylim) - eh / 2, mean(ylim) + eh / 2))
}

# One basemap panel: the composed tile raster (or flat panel) across its
# range with no axes; meters per unit stay equal on both axes. The optional
# lake outline is drawn as part of the panel base, because its sf layers
# must be added before the coordinate system: in this ggplot2 build adding a
# geom_sf afterwards silently replaces the coordinate (and with it the fixed
# limits and the forced equal-aspect box), and the panel collapses.
map_panel <- function(bm, title = NULL, title_size = 16, border_col = NULL,
                      aspect = 1, lake = FALSE,
                      lake_colour = "#17BECF", lake_alpha = 0.35,
                      lake_lw = 1.0) {
  ae <- aspect_extent(bm, aspect)
  p <- ggplot()
  if (isTRUE(bm$flat)) {
    p <- p + annotate("rect", xmin = bm$xmin, xmax = bm$xmax,
      ymin = bm$ymin, ymax = bm$ymax, fill = bm$hex, colour = NA)
  } else {
    p <- p + geom_tile(data = bm$tiles, aes(x = x, y = y, fill = col),
      inherit.aes = FALSE, width = bm$cell, height = bm$cell) +
      scale_fill_identity() + guides(fill = "none")
  }
  if (lake) {
    # Teal lake overlay (Agamon panels + p_goal panels); in the Python
    # z-order the lake sits below the track lines.
    p <- p +
      geom_sf(data = sf::st_sf(geometry = lake_union_3857), inherit.aes = FALSE,
        fill = lake_colour, colour = lake_colour, alpha = lake_alpha,
        linewidth = lake_lw * 0.3528) +
      geom_sf(data = sf::st_sf(geometry = sf::st_boundary(lake_union_3857)),
        inherit.aes = FALSE, colour = lake_colour, fill = NA,
        linewidth = (lake_lw + 0.2) * 0.3528)
  }
  # theme_void() carries its own plot.title (relative size, left-just);
  # compose the explicit title/border themes AFTER it so they win.
  th <- theme_void()
  if (!is.null(border_col)) {
    th <- th + theme(panel.border = element_rect(
      colour = border_col, fill = NA, linewidth = 0.35))
  }
  if (!is.null(title)) {
    p <- p + labs(title = title)
    th <- th + theme(plot.title = element_text(size = title_size, hjust = 0.5,
      margin = margin(0, 0, 2, 0)))
  }
  p + th +
    coord_sf(crs = sf::st_crs(3857), datum = NA, expand = FALSE,
      xlim = ae$xlim, ylim = ae$ylim)
}

# Marker radius in meters for a matplotlib marker-size s (pt^2) drawn on
# a panel of panel_in inches that shows extent_m meters wide.
marker_radius_m <- function(s_pt2, panel_in, extent_m)
  sqrt(s_pt2) / 72 / panel_in * extent_m / 2

star_poly_pts <- function(cx, cy, r) {
  k <- 0:9
  ang <- (90 - 36 * k) * pi / 180
  rad <- ifelse(k %% 2 == 0L, r, r * 0.38196601)
  cbind(cx + rad * cos(ang), cy + rad * sin(ang))
}

add_roost_marker <- function(r) {
  st <- star_poly_pts(ROOST_XY[1], ROOST_XY[2], r)
  annotate("polygon", x = st[, 1], y = st[, 2],
    fill = "white", colour = "black", linewidth = 0.35)
}

add_x_markers <- function(points_xy, r) {
  do.call(c, lapply(seq_len(nrow(points_xy)), function(i) {
    cx <- points_xy[i, 1]; cy <- points_xy[i, 2]
    c(
      annotate("segment", x = cx - r, y = cy - r, xend = cx + r, yend = cy + r,
        linewidth = 1.0, colour = "white"),
      annotate("segment", x = cx + r, y = cy - r, xend = cx - r, yend = cy + r,
        linewidth = 1.0, colour = "white"),
      annotate("segment", x = cx - r, y = cy - r, xend = cx + r, yend = cy + r,
        linewidth = 0.45, colour = "black"),
      annotate("segment", x = cx + r, y = cy - r, xend = cx - r, yend = cy + r,
        linewidth = 0.45, colour = "black")
    )
  }))
}

# Scale bar + north arrow in a panel corner (data coordinates), the R
# twin of the Python helper: a two-tone bar with an outline, the length
# label below it (nice 1-2-5 length ~20% of the panel width) and a
# compass arrow centered above the bar. Returns the plot with the layers.
nice_scale_len <- function(width_m) {
  target <- max(width_m * 0.20, 1)
  e <- floor(log10(target))
  b <- target / 10 ^ e
  nice <- if (b < 1.5) 1 else if (b < 3.5) 2 else if (b < 7.5) 5 else 10
  nice * 10 ^ e
}

fmt_scale_label <- function(m)
  paste0(as.integer(round(m)), " m")

add_scale_compass <- function(p, xlim, ylim, side = "right",
                              panel_w_in = 2.5, panel_h_in = 5.85,
                              variant = "classic") {
  w_m <- diff(xlim); h_m <- diff(ylim)
  pt2m <- h_m / (panel_h_in * 72)
  L <- nice_scale_len(w_m)
  xr <- if (side == "right") xlim[2] else xlim[1]
  if (identical(variant, "classic")) {
    # 16-px corner pad (100 dpi), 7-pt bar, 9-pt label, 30-pt compass -
    # the Agamon map geometry.
    pad_pt <- 16 * 72 / 100
    seg_h_pt <- 7
    breath_pt <- 7
    lab_pt <- 9
    bar_edge_pt <- 0.9
    tail_pt <- 1.9
    noff <- 2
    ytxt <- ylim[1] + (pad_pt + 10.5) * pt2m
    yb0 <- ytxt + breath_pt * pt2m
    yb1 <- yb0 + seg_h_pt * pt2m
    cbase <- (yb0 + yb1) / 2 + (seg_h_pt / 2 + 8) * pt2m
    arrow_pt <- 30; hw_pt <- 9; hl_pt <- 10
  } else {
    # p_goal geometry: equal 22-px (full) / 14-px (compact) corner pads at
    # 80 dpi, 8-pt bar, 10-pt label, compass with its base on the bar top
    # (20-pt arrow full / 18-pt compact) and the N 10-pt.
    small <- identical(variant, "pgoal_small")
    pad_pt <- if (small) 14 * 72 / 80 else 22 * 72 / 80
    seg_h_pt <- 8
    breath_pt <- if (small) 5 else 6
    lab_pt <- 10
    bar_edge_pt <- if (small) 0.8 else 0.9
    tail_pt <- if (small) 1.8 else 2.0
    noff <- if (small) 2 else 3
    # 10-pt font bbox is ~13 pt tall: the label's bottom edge sits on the
    # pad line, the bar breathes above it.
    ytxt <- ylim[1] + pad_pt * pt2m
    yb0 <- ytxt + (13 + breath_pt) * pt2m
    yb1 <- yb0 + seg_h_pt * pt2m
    cbase <- yb1 + 2 * pt2m
    arrow_pt <- if (small) 18 else 20
    hw_pt <- if (small) 8 else 9
    hl_pt <- if (small) 9 else 10
  }
  ctip <- cbase + arrow_pt * pt2m
  hw <- hw_pt / 2 * pt2m; hl <- hl_pt * pt2m
  x1 <- xr - pad_pt * pt2m
  x0 <- x1 - L
  xc <- (x0 + x1) / 2
  # vjust in this build: 1 pins the text TOP at the anchor, 0 pins the
  # BOTTOM - classic keeps its top-pinned label, p_goal pins the label
  # bottom on the pad line.
  lab_vjust <- if (identical(variant, "classic")) 1 else 0
  p +
    geom_rect(inherit.aes = FALSE,
      data = data.frame(xmin = x0, xmax = xc, ymin = yb0, ymax = yb1),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#ffffff", colour = NA) +
    geom_rect(inherit.aes = FALSE,
      data = data.frame(xmin = xc, xmax = x1, ymin = yb0, ymax = yb1),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#000000", colour = NA) +
    geom_rect(inherit.aes = FALSE,
      data = data.frame(xmin = x0, xmax = x1, ymin = yb0, ymax = yb1),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, colour = "#000000",
      linewidth = bar_edge_pt * 0.3528) +
    annotate("text", x = xc, y = ytxt,
      label = fmt_scale_label(L), size = lab_pt * 0.3528, colour = "#000000",
      hjust = 0.5, vjust = lab_vjust) +
    annotate("segment", x = xc, y = cbase,
      xend = xc, yend = ctip - hl + 2 * pt2m, linewidth = tail_pt * 0.3528,
      colour = "#000000") +
    annotate("polygon", x = c(xc, xc - hw, xc + hw),
      y = c(ctip, ctip - hl, ctip - hl), fill = "#000000", colour = "#000000",
      linewidth = 0.9 * 0.3528) +
    annotate("text", x = xc, y = ctip + noff * pt2m, label = "N",
      size = lab_pt * 0.3528, colour = "#000000", hjust = 0.5, vjust = 0)
}

# Hand-drawn figure-level legend panel (the Python fig.legend items): a
# framed white box with marker/line glyphs and labels, no axes.
# Side/bottom legend: a white box with one stacked column of keys (the
# destination/flow map columns) or one packed horizontal row of keys
# (the p_goal figure legend). Data coordinates are inches of the cell the
# panel is measured to occupy, so keys stay round and the box keeps its
# fixed shape however patchwork allocates the cell.
legend_glyph <- function(p, kind, gx, gy) {
  if (kind == "star") {
    st <- star_poly_pts(gx, gy, 0.07)
    p + annotate("polygon", x = st[, 1], y = st[, 2],
      fill = "white", colour = "black")
  } else if (kind == "x") {
    r <- 0.065
    p +
      annotate("segment", x = gx - r, y = gy - r, xend = gx + r, yend = gy + r,
        linewidth = 1.2, colour = "white") +
      annotate("segment", x = gx + r, y = gy - r, xend = gx - r, yend = gy + r,
        linewidth = 1.2, colour = "white") +
      annotate("segment", x = gx - r, y = gy - r, xend = gx + r, yend = gy + r,
        linewidth = 0.5, colour = "black") +
      annotate("segment", x = gx + r, y = gy - r, xend = gx - r, yend = gy + r,
        linewidth = 0.5, colour = "black")
  } else if (kind == "dash") {
    p + annotate("segment", x = gx - 0.09, y = gy, xend = gx + 0.09, yend = gy,
      linewidth = 0.5, colour = "black", linetype = "dashed")
  } else {  # swatch (Lake): pale cyan chip
    p + annotate("rect", xmin = gx - 0.07, xmax = gx + 0.07,
      ymin = gy - 0.05, ymax = gy + 0.05, fill = "#17BECF", alpha = 0.35,
      colour = "#17BECF", linewidth = 0.3)
  }
}

legend_panel <- function(items, cell_w = 1.8, cell_h = 5.85, font_mm = 3.0) {
  if (is.null(items$n)) items$n <- rep(1, nrow(items))
  nrow_i <- length(unique(items$n))
  # rows 0.34 in apart; the box hugs them with a 0.15 in pad. A single
  # row keeps the bare 0.15 in pad.
  box_half <- if (nrow_i > 1L) ((nrow_i - 1L) / 2) * 0.34 + 0.15 else 0.15
  p <- ggplot() +
    annotate("rect", xmin = 0.03, xmax = cell_w - 0.03,
      ymin = cell_h / 2 - box_half, ymax = cell_h / 2 + box_half,
      fill = "#ffffff", colour = "#000000", linewidth = 0.25)
  if (nrow_i > 1L) {
    for (j in seq_len(nrow_i)) {
      it <- items[items$n == j, , drop = FALSE]
      gy <- cell_h / 2 + ((nrow_i - 1L) / 2 - (j - 1L)) * 0.34
      p <- legend_glyph(p, it$kind, 0.32, gy)
      # annotate("text") sizes are millimetres of glyph height here
      # (3.0 mm ~ 8.5 pt); theme element_text sizes are points.
      p <- p + annotate("text", x = 0.55, y = gy, label = it$label,
        size = font_mm, hjust = 0, vjust = 0.5)
    }
  } else {
    # one packed, centred row: key, label, gap, key, label, ...
    lw_in <- strwidth(items$label, units = "inches") * (8 / 12)
    key_in <- 0.18; gk_in <- 0.06; gi_in <- 0.16
    tot_in <- sum(lw_in) + nrow(items) * (key_in + gk_in) +
      (nrow(items) - 1) * gi_in
    starts_in <- (cell_w - tot_in) / 2 +
      cumsum(c(0, key_in + gk_in + lw_in[-nrow(items)] + gi_in))
    for (i in seq_len(nrow(items))) {
      p <- legend_glyph(p, items$kind[i], starts_in[i] + key_in / 2, cell_h / 2)
      p <- p + annotate("text", x = starts_in[i] + key_in + gk_in,
        y = cell_h / 2, label = items$label[i], size = font_mm,
        hjust = 0, vjust = 0.5)
    }
  }
  p + coord_cartesian(xlim = c(0, cell_w), ylim = c(0, cell_h), expand = FALSE) +
    theme_void()
}

destination_panels <- function(segs, bm, xlim, ylim, r_roost, r_trans) {
  panel_dest <- function(dest, title, lake = FALSE) {
    sd <- segs[segs$destination == dest, , drop = FALSE]
    p <- map_panel(bm, title, title_size = 15, aspect = 0.43,
      lake = lake, lake_colour = "#56B4E9", lake_alpha = 0.35, lake_lw = 1.2)
    p +
      geom_segment(data = sd[sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = treatment),
        linetype = "dashed", linewidth = 1.44 * 0.3528, alpha = 0.85) +
      geom_segment(data = sd[!sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = treatment),
        linewidth = 1.6 * 0.3528, alpha = 0.95) +
      scale_colour_manual(values = treatment_colors_cb) +
      guides(colour = "none") +
      add_roost_marker(r_roost) +
      add_x_markers(TRANS_XY, r_trans)
  }
  p_dest_home <- panel_dest("Home", "Home")
  p_dest_agam <- panel_dest("Agamon", "Agamon Hula Lake", lake = TRUE)
  p_dest_agam <- add_scale_compass(p_dest_agam, xlim, ylim,
    panel_w_in = 2.53, panel_h_in = 5.88)
  p_dest_legend <- legend_panel(data.frame(
    n = c(1, 2, 3), kind = c("star", "x", "dash"),
    label = c("Roost", "Translocation", "Signal Dropout")),
    cell_w = 1.8, cell_h = 5.85)
  (p_dest_home | p_dest_agam | p_dest_legend) +
    plot_layout(widths = c(0.37, 0.37, 0.26)) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

destination_panels_compact <- function(segs, bm, xlim, ylim, r_roost, r_trans) {
  panel_compact <- function(dest, title, lake = FALSE) {
    sd <- segs[segs$destination == dest, , drop = FALSE]
    p <- map_panel(bm, title, title_size = 9, border_col = "black",
      aspect = 1.22 / 3.10, lake = lake, lake_alpha = 0.35, lake_lw = 1.0)
    p +
      geom_segment(data = sd[sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = treatment),
        linetype = "dashed", linewidth = 1.22 * 0.3528, alpha = 0.85) +
      geom_segment(data = sd[!sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = treatment),
        linewidth = 1.35 * 0.3528, alpha = 0.95) +
      scale_colour_manual(values = treatment_colors_cb) +
      guides(colour = "none") +
      add_roost_marker(r_roost) +
      add_x_markers(TRANS_XY, r_trans)
  }
  # Compact legend inside the Agamon panel, top-right (the Python twin
  # draws it at (0.985, 0.985) of that panel).
  add_compact_legend <- function(p) {
    # Point scale from the panel's rendered height (2.5x3.8 page with this
    # layout draws the equal-aspect panel ~2.855in tall), so keys, frame pad
    # and 6pt labels land at true point sizes.
    m_per_in <- diff(ylim) / 2.855
    pp <- function(pts) pts * m_per_in / 72
    fs_pt <- 6
    labels <- c("Roost", "Translocation", "Signal dropout")
    lw_in <- strwidth(labels, units = "inches") * fs_pt / 12
    pad <- pp(0.35 * fs_pt)          # borderpad
    mcol <- pp(7.2)                  # marker column (markersize 7.2)
    hpad <- pp(0.45 * fs_pt)         # handletextpad
    pitch <- pp(1.35 * fs_pt)        # row pitch (labelspacing 0.35)
    box_w <- 2 * pad + mcol + hpad + max(lw_in) * 1.05 * m_per_in
    box_h <- 2 * pad + pp(fs_pt) + (length(labels) - 1) * pitch
    x1 <- xlim[1] + 0.985 * diff(xlim)
    x0 <- x1 - box_w
    y1 <- ylim[2] - (1 - 0.985) * diff(ylim)
    y0 <- y1 - box_h
    gx <- x0 + pad + mcol / 2
    tx <- x0 + pad + mcol + hpad
    ys <- y1 - (pad + pp(fs_pt) / 2 + (0:(length(labels) - 1)) * pitch)
    r_s <- pp(7.2 * 0.9 / 2)         # Roost star radius (scale 0.9)
    r_x <- pp(6.2 * 0.9 / 2)         # X half-diagonal
    d_h <- pp(1.4 * fs_pt / 2)       # dashed-handle half-length
    st <- star_poly_pts(gx, ys[1], r_s)
    p +
      annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1,
        fill = "#fbfbfb", colour = "#cccccc", linewidth = 0.3528) +
      annotate("polygon", x = st[, 1], y = st[, 2],
        fill = "white", colour = "black", linewidth = 0.8 * 0.3528) +
      annotate("segment", x = gx - r_x, y = ys[2] - r_x,
        xend = gx + r_x, yend = ys[2] + r_x, linewidth = 0.6, colour = "white") +
      annotate("segment", x = gx + r_x, y = ys[2] - r_x,
        xend = gx - r_x, yend = ys[2] + r_x, linewidth = 0.6, colour = "white") +
      annotate("segment", x = gx - r_x, y = ys[2] - r_x,
        xend = gx + r_x, yend = ys[2] + r_x, linewidth = 0.8 * 0.3528, colour = "black") +
      annotate("segment", x = gx + r_x, y = ys[2] - r_x,
        xend = gx - r_x, yend = ys[2] + r_x, linewidth = 0.8 * 0.3528, colour = "black") +
      annotate("text", x = tx, y = ys[1], label = "Roost",
        size = fs_pt * 0.3528, hjust = 0, vjust = 0.5) +
      annotate("text", x = tx, y = ys[2], label = "Translocation",
        size = fs_pt * 0.3528, hjust = 0, vjust = 0.5) +
      annotate("segment", x = gx - d_h, y = ys[3], xend = gx + d_h, yend = ys[3],
        linewidth = 0.3528, colour = "black", linetype = "dashed") +
      annotate("text", x = tx, y = ys[3], label = "Signal dropout",
        size = fs_pt * 0.3528, hjust = 0, vjust = 0.5)
  }
  p_compact_home <- panel_compact("Home", "Home")
  p_compact_home <- add_scale_compass(p_compact_home, xlim, ylim,
    panel_w_in = 1.22, panel_h_in = 3.10)
  p_compact_agam <- panel_compact("Agamon", "Agamon Hula Lake", lake = TRUE)
  p_compact_agam <- add_compact_legend(p_compact_agam)
  (p_compact_home | p_compact_agam) +
    theme(plot.margin = margin(1, 3, 1, 3))
}

flow_field_panels <- function(bm, xlim, ylim, r_roost, r_trans) {
  # Flow grids (200 m, the Python resolution) over the full panel extent.
  flow_res <- 200
  gx <- seq(ceiling(xlim[1] / flow_res) * flow_res,
    floor(xlim[2] / flow_res) * flow_res, by = flow_res)
  gy <- seq(ceiling(ylim[1] / flow_res) * flow_res,
    floor(ylim[2] / flow_res) * flow_res, by = flow_res)
  gpts <- expand.grid(x = gx, y = gy)
  lake_mask <- lengths(sf::st_within(
    sf::st_as_sf(gpts, coords = c("x", "y"), crs = 3857), lake_union_3857)) > 0L
  nxf <- length(gx); nyf <- length(gy)
  # Field matrices follow the twin's meshgrid layout: row = gy index,
  # column = gx index. expand.grid cycles x fastest, so point k =
  # (b - 1) * nxf + a sits at row b, column a.
  M <- matrix(FALSE, nrow = nyf, ncol = nxf)
  M[cbind((seq_len(nrow(gpts)) - 1L) %/% nxf + 1L,
    (seq_len(nrow(gpts)) - 1L) %% nxf + 1L)] <- lake_mask
  # squared distance (in grid steps) to the nearest lake cell: the exact
  # per-pixel minimum, the twin of scipy distance_transform_edt(~mask)
  d2 <- matrix(Inf, nyf, nxf)
  li <- which(M, arr.ind = TRUE)
  for (i in seq_len(nrow(li))) {
    d2 <- pmin(d2, outer((seq_len(nyf) - li[i, 1])^2,
      (seq_len(nxf) - li[i, 2])^2, "+"))
  }
  dist_f <- sqrt(d2) * flow_res
  dist_f[M] <- 0
  # Gaussian pass on each axis with the twin's sigma (3 px) and nearest
  # edge handling (scipy gaussian_filter mode="nearest"): each axis is
  # pre-padded with its edge value before the pass. filter() is not used
  # because it returns NAs over the 9-tap edge on a plain vector, which on
  # this panel's 26-column axis would erase 18 of 26 grid columns.
  gkern <- dnorm(seq(-9L, 9L), 0, 3); gkern <- gkern / sum(gkern)
  half_k <- length(gkern) %/% 2L
  gpass1d <- function(v) {
    vq <- c(rep(v[1], half_k), v, rep(v[length(v)], half_k))
    as.numeric(stats::convolve(vq, gkern, type = "filter"))
  }
  dist_f <- apply(dist_f, 2, gpass1d)
  dist_f <- t(apply(dist_f, 1, gpass1d))
  # unit vectors pointing at the lake (negative gradient of distance)
  txf <- matrix(0, nyf, nxf); tyf <- matrix(0, nyf, nxf)
  txf[, 2:(nxf - 1L)] <- -(dist_f[, 3:nxf] - dist_f[, 1:(nxf - 2L)]) / (2 * flow_res)
  txf[, 1L] <- if (nxf > 1L) -(dist_f[, 2L] - dist_f[, 1L]) / flow_res else 0
  txf[, nxf] <- if (nxf > 1L) -(dist_f[, nxf] - dist_f[, nxf - 1L]) / flow_res else 0
  tyf[2:(nyf - 1L), ] <- -(dist_f[3:nyf, ] - dist_f[1:(nyf - 2L), ]) / (2 * flow_res)
  tyf[1L, ] <- if (nyf > 1L) -(dist_f[2L, ] - dist_f[1L, ]) / flow_res else 0
  tyf[nyf, ] <- if (nyf > 1L) -(dist_f[nyf, ] - dist_f[nyf - 1L, ]) / flow_res else 0
  normf <- sqrt(txf^2 + tyf^2)
  txf <- txf / normf; tyf <- tyf / normf
  txf[normf <= 1e-12] <- NA_real_
  tyf[normf <= 1e-12] <- NA_real_
  txf[M > 0] <- NA_real_
  tyf[M > 0] <- NA_real_
  # radial field toward the roost (no mask)
  rx <- gpts$x - ROOST_XY[1]; ry <- gpts$y - ROOST_XY[2]
  nr <- sqrt(rx^2 + ry^2)
  txr <- -rx / nr; tyr <- -ry / nr
  # Arrow length: the Python quiver with unit vectors and scale 30
  # (default scale_units = "width") draws each unit vector 1/30 of
  # the panel width, so the whole arrow (shaft + head) is this long.
  arrow_len_f <- diff(xlim) / 30
  # expand.grid order: x index fastest, so as.vector(M) matches the point
  # rows; the fields must be transposed first for the same reason.
  lake_dir <- data.frame(x = gpts$x, y = gpts$y,
    xend = gpts$x + arrow_len_f * as.vector(t(txf)), yend = gpts$y + arrow_len_f * as.vector(t(tyf)),
    in_lake = lake_mask)
  lake_dir <- lake_dir[!lake_dir$in_lake & is.finite(lake_dir$xend), , drop = FALSE]
  roost_dir <- data.frame(x = gpts$x, y = gpts$y,
    xend = gpts$x + arrow_len_f * txr, yend = gpts$y + arrow_len_f * tyr)
  roost_dir <- roost_dir[is.finite(roost_dir$xend), , drop = FALSE]

  p_flow_home <- map_panel(bm, "Home (roost flow)", title_size = 14,
    aspect = 0.41) +
    geom_segment(data = roost_dir, inherit.aes = FALSE,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(0.04, "inches"), type = "closed"),
      linewidth = 0.35, alpha = 0.95, colour = "#ffffff") +
    add_roost_marker(r_roost) +
    add_x_markers(TRANS_XY, r_trans)
  p_flow_agam <- map_panel(bm, "Agamon Hula Lake", title_size = 14,
    aspect = 0.41, lake = TRUE, lake_alpha = 0.35, lake_lw = 1.2) +
    geom_segment(data = lake_dir, inherit.aes = FALSE,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(0.04, "inches"), type = "closed"),
      linewidth = 0.35, alpha = 0.95, colour = "#ffffff") +
    add_roost_marker(r_roost) +
    add_x_markers(TRANS_XY, r_trans)
  p_flow_legend <- legend_panel(data.frame(
    n = c(1, 2), kind = c("star", "x"),
    label = c("Roost", "Translocation")),
    cell_w = 1.9, cell_h = 7.3)
  (p_flow_home | p_flow_agam | p_flow_legend) +
    plot_layout(widths = c(0.38, 0.38, 0.24)) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# Per-tag polar KDE panels: 0 deg at right with the counter-clockwise
# orientation, one curve per tag in the treatment colour, the two density
# rings (r = 1 at density 0, r = 2 at density 1), four grey direction
# spokes, and a darker 0-deg segment from r = 1 to r = 2.
polar_data_panels <- function(kde_data, tag_meta) {
  xmap_polar <- function(a) (pi / 2 - a) %% (2 * pi)
  kde_data$theta_x <- xmap_polar(kde_data$theta)
  kde_data$taglabel <- factor(paste0("tag ", kde_data$tag),
    levels = paste0("tag ", as.character(tag_meta$tag)))
  rmax_polar <- 1 + max(1, max(kde_data$density)) + 0.2
  pg_a <- seq(0, 2 * pi, length.out = 361L)
  polar_rings <- data.frame(theta = rep(xmap_polar(pg_a), 2),
    y = c(rep(1, 361L), rep(2, 361L)), ring = rep(1:2, each = 361L))
  spoke_a <- c(-pi / 2, 0, pi / 2, pi)
  polar_spokes <- data.frame(theta = rep(xmap_polar(spoke_a), each = 2),
    y = rep(c(1, rmax_polar), 4), sp = rep(seq_along(spoke_a), each = 2))
  polar_zero <- data.frame(theta = c(xmap_polar(0), xmap_polar(0)), y = c(1, 2))
  nc_polar <- ceiling(sqrt(nrow(tag_meta)))
  ggplot(kde_data, aes(theta_x, 1 + density, colour = treatment,
      group = interaction(treatment, tag))) +
    geom_line(linewidth = 0.7, alpha = 0.95) +
    geom_path(data = polar_rings, inherit.aes = FALSE, aes(theta, y, group = ring),
      colour = "#999999", linewidth = 0.32, alpha = 0.65) +
    geom_path(data = polar_spokes, inherit.aes = FALSE, aes(theta, y, group = sp),
      colour = "#bfbfbf", linewidth = 0.28, alpha = 0.7) +
    geom_path(data = polar_zero, inherit.aes = FALSE, aes(theta, y),
      colour = "#a6a6a6", linewidth = 0.57, alpha = 0.95) +
    facet_wrap(~ taglabel, ncol = nc_polar) +
    scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0),
      breaks = xmap_polar(c(0, pi / 2, pi, -pi / 2)), labels = c("0°", "90°",
      "±180°", "-90°")) +
    coord_polar() +
    scale_y_continuous(limits = c(0, rmax_polar), expand = c(0, 0)) +
    scale_colour_manual(values = treatment_colors, name = NULL) +
    guides(colour = "none") +
    theme_minimal(base_size = 9) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6, colour = "grey25"),
      panel.grid = element_blank(), panel.border = element_blank(),
      strip.background = element_blank(), strip.text = element_text(size = 9)) +
    labs(x = "", y = "")
}

# Axis labels: comma parts of the raw treatment name, sorted, "control"
# dropped when other components are present, first letter of each part
# upper-cased ("Anosmia, Eye cover, Magnet").
fmt_treat_name <- function(s) {
  parts <- trimws(strsplit(tolower(s), ",", fixed = TRUE)[[1]])
  parts <- parts[nchar(parts) > 0]
  if (length(parts) > 1L && "control" %in% parts) parts <- parts[parts != "control"]
  if (!length(parts)) parts <- "control"
  parts <- paste0(toupper(substring(sort(parts), 1, 1)),
    substring(sort(parts), 2))
  paste(parts, collapse = ", ")
}

# One panel per treatment: the per-tag estimates as jittered dots, the
# median as a filled circle, and the 25-75 percentile range as a bar.
desc_plot <- function(variable, label, ylim, descriptive, treat_name_fmt, grid_y = FALSE) {
  set.seed(0)
  blocks <- split(descriptive, descriptive$treatment)
  jit <- lapply(blocks, function(g) runif(nrow(g), -0.10, 0.10))
  dots <- do.call(rbind, lapply(seq_along(blocks), function(i) {
    data.frame(x = (i - 1) + jit[[i]], y = blocks[[i]][[variable]],
      treatment = blocks[[i]]$treatment[1])
  }))
  # Dots beyond the fixed axis limits are clipped at the edge, as in the
  # source; dropping them keeps the render identical and quiet.
  dots <- dots[is.finite(dots$y) & dots$y >= ylim[1] & dots$y <= ylim[2], ]
  stats <- do.call(rbind, lapply(seq_along(blocks), function(i) {
    g <- blocks[[i]]
    q <- quantile(g[[variable]], c(0.25, 0.50, 0.75), names = FALSE)
    data.frame(x = i - 1, treatment = g$treatment[1],
      med = q[2], lo = q[1], hi = q[3])
  }))
  p <- ggplot(dots, aes(x, y, colour = treatment, group = treatment)) +
    geom_point(shape = 16, size = 1.23, alpha = 0.9) +
    geom_errorbar(data = stats, aes(y = med, ymin = lo, ymax = hi),
      colour = "black", linewidth = 0.35, width = 0.14) +
    geom_point(data = stats, aes(y = med, fill = treatment),
      shape = 21, size = 1.76, colour = "black") +
    scale_colour_manual(values = treatment_colors, name = NULL) +
    scale_fill_manual(values = treatment_colors, name = NULL) +
    guides(colour = "none", fill = "none") +
    scale_x_continuous(breaks = 0:(length(blocks) - 1L),
      labels = treat_name_fmt, expand = c(0.2, 0.2)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0),
      breaks = if (variable == "kappa_hat") c(0, 5, 10) else c(0, 0.5, 1)) +
    labs(x = "", y = label) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  if (grid_y) {
    p <- p + theme(panel.grid.major.y = element_line(linewidth = 0.2,
      colour = "grey80"))
  }
  p
}

total_elpd_panel <- function(res_total, vline = TRUE, xlabel = "PSIS-LOO ELPD") {
  best <- res_total$elpd[which.max(res_total$elpd)]
  p <- ggplot(res_total, aes(elpd, reorder(model, elpd))) +
    geom_errorbar(aes(xmin = elpd - se_elpd, xmax = elpd + se_elpd),
      orientation = "y", width = 0.15, linewidth = 0.35) +
    geom_point(size = 2)
  if (vline) p <- p + geom_vline(xintercept = best, linetype = "dashed", linewidth = 0.3)
  p + theme_classic(base_size = 9) +
    labs(x = xlabel, y = "")
}

# Total + per-treatment deltas: Total row on top (black), treatments in
# tab10 below, named y ticks, dashed line at zero.
elpd_delta_plot <- function(delta_df) {
  delta_df$rcol <- c("black", treatment_colors[as.character(levels(df_track$treatment))])
  ggplot(delta_df, aes(delta_elpd, ypos)) +
    geom_errorbar(aes(xmin = delta_elpd - se_delta_elpd, xmax = delta_elpd + se_delta_elpd,
      colour = rcol), orientation = "y", width = 0.18, linewidth = 0.35) +
    geom_point(aes(colour = rcol), size = 2.2) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_y_reverse(limits = c(-0.41, 4 * 1.35 + 0.41),
      breaks = delta_df$ypos, labels = delta_df$label, expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
    scale_colour_identity(guide = "none") +
    theme_classic(base_size = 9) +
    labs(x = "ELPD difference", y = "")
}

# Bracket panel for the combined summaries: colorblind treatment rows,
# only the "Total" tick, a left-hand bracket over the treatment rows and
# a left outward tick + "Within-treatment" label sitting in the left
# margin (plot.margin below), with clipping off so they can render there.
elpd_bracket_plot <- function(delta_df) {
  brk <- delta_df
  brk$rcol <- c("black", treatment_colors_cb[as.character(levels(df_track$treatment))])
  gx <- range(c(brk$delta_elpd - brk$se_delta_elpd, brk$delta_elpd + brk$se_delta_elpd))
  padx <- 0.05 * diff(gx)
  xlim_b <- c(gx[1] - padx, gx[2] + padx)
  # Bracket geometry mirrors the Python twin, which draws in axes-fraction
  # coordinates with clipping off: the bracket sits on the left axis line,
  # the caps go into the plot, and an outward tick reaches a label that sits
  # entirely in the reserved left margin (plot.margin below, clip off above).
  xb <- xlim_b[1]                    # bracket on the y-axis line
  cap_w <- 0.020 * diff(xlim_b)      # caps into the plot (right)
  tick_w <- 0.055 * diff(xlim_b)     # tick out to the label (left)
  ymid_b <- (1.35 + 4 * 1.35) / 2    # middle of the treatment rows
  ggplot(brk, aes(delta_elpd, ypos)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(xmin = delta_elpd - se_delta_elpd, xmax = delta_elpd + se_delta_elpd,
      colour = rcol), orientation = "y", width = 0.18, linewidth = 0.35) +
    geom_point(aes(colour = rcol), size = 2.2) +
    annotate("segment", x = xb, y = 1.35, xend = xb, yend = 4 * 1.35,
      linewidth = 0.35, colour = "black") +
    annotate("segment", x = xb, y = 1.35, xend = xb + cap_w, yend = 1.35,
      linewidth = 0.35, colour = "black") +
    annotate("segment", x = xb, y = 4 * 1.35, xend = xb + cap_w, yend = 4 * 1.35,
      linewidth = 0.35, colour = "black") +
    annotate("segment", x = xb, y = ymid_b, xend = xb - tick_w, yend = ymid_b,
      linewidth = 0.35, colour = "black") +
    annotate("text", x = xb - tick_w - 0.012 * diff(xlim_b), y = ymid_b,
      label = "Within-treatment", size = 2.8, colour = "black",
      hjust = 1, vjust = 0.5) +
    scale_y_reverse(limits = c(-0.54, 4 * 1.35 + 0.54), breaks = 0,
      labels = "Total", expand = c(0, 0)) +
    scale_x_continuous(limits = xlim_b, oob = scales::oob_keep) +
    coord_cartesian(clip = "off") +
    scale_colour_identity(guide = "none") +
    theme_classic(base_size = 9) +
    theme(plot.margin = margin(5, 5, 20, 0.95 * 72)) +
    labs(x = "ELPD difference", y = "")
}

elpd_individual_plot <- function(delta_individual) {
  delta_individual$tlabel <- factor(treatment_nice[as.character(delta_individual$treatment)],
    levels = treatment_nice[levels(df_track$treatment)])
  di_grid <- range(c(delta_individual$delta_elpd - delta_individual$se_delta_elpd,
    delta_individual$delta_elpd + delta_individual$se_delta_elpd))
  di_pad <- 0.06 * diff(di_grid)
  ggplot(delta_individual,
    aes(delta_elpd, reorder(tag, delta_elpd), colour = treatment)) +
    geom_errorbar(aes(xmin = delta_elpd - se_delta_elpd, xmax = delta_elpd + se_delta_elpd),
      orientation = "y", width = 0.15, linewidth = 0.3) +
    geom_point(size = 1.4) +
    facet_wrap(~ tlabel, ncol = 4, scales = "free_y") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_colour_manual(values = treatment_colors_cb[tlabels_all], guide = "none") +
    scale_x_continuous(limits = di_grid - c(di_pad, -di_pad), expand = c(0, 0)) +
    labs(x = "ELPD difference (vMU - vM)", y = "") +
    theme_classic(base_size = 8) + theme(strip.background = element_blank(),
      strip.text = element_text(size = 9))
}

# Twin look: colorblind rows, hidden treatment ticks/labels, twin x label
# and reference line. (Python plot_mean_hdi_horizontal, no-labels variant.)
mean_hdi_plot <- function(table, quantity, label, reference = NULL, xlimits = NULL,
  flip = FALSE) {
  x <- table[table$quantity == quantity, ]
  p <- ggplot(x, aes(mean, treatment, colour = treatment)) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.18,
      linewidth = 0.4) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = treatment_colors_cb[tlabels_all], guide = "none") +
    scale_y_discrete(
      limits = if (flip) levels(df_track$treatment)
               else rev(levels(df_track$treatment)), labels = "") +
    labs(x = label, y = "") +
    theme_classic(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 12, 5.5, 5.5))
  if (!is.null(reference)) p <- p + geom_vline(xintercept = reference, linetype = "dashed")
  if (!is.null(xlimits)) p <- p + scale_x_continuous(limits = xlimits, expand = c(0, 0))
  p
}

# Twin look: one row of four named panels, animals sorted by value within
# each, a light-gray population HDI band with a solid population-mean
# line, optional dashed reference, one shared x-axis + figure label.
tag_hdi_plot <- function(quantity, label, reference = NULL, xlimits = NULL) {
  x <- tag_summary[tag_summary$quantity == quantity, ]
  pop <- pop_summary[pop_summary$quantity == quantity, ]
  lv <- tlabels_all
  if (is.null(xlimits)) {
    gr <- range(c(min(x$low, pop$low), max(x$high, pop$high)))
    if (!is.null(reference)) gr <- range(c(gr, reference))
    pd <- 0.06 * diff(gr)
    xlimits <- gr - c(pd, -pd)
  }
  panels <- lapply(lv, function(t) {
    sub <- x[x$treatment == t, ]
    sub <- sub[order(sub$mean, decreasing = TRUE), ]
    ny <- nrow(sub)
    poprow <- pop[pop$treatment == t, ]
    panel <- ggplot(sub, aes(mean, y = seq_len(ny))) +
      annotate("rect", xmin = poprow$low, xmax = poprow$high,
        ymin = 0.5, ymax = ny + 0.5,
        fill = "#d9d9d9", colour = NA, inherit.aes = FALSE) +
      geom_vline(xintercept = poprow$mean, linewidth = 0.3, colour = "#000000") +
      geom_errorbar(aes(xmin = low, xmax = high, colour = treatment),
        orientation = "y", width = 0.22, linewidth = 0.3) +
      geom_point(aes(colour = treatment), size = 1.5) +
      scale_colour_manual(values = treatment_colors_cb[tlabels_all], guide = "none") +
      scale_y_reverse(limits = c(0.5, ny + 0.5), expand = c(0, 0),
        breaks = seq_len(ny), labels = as.character(sub$tag)) +
      scale_x_continuous(limits = xlimits, expand = c(0, 0)) +
      labs(title = treatment_nice[[t]], x = label, y = "") +
      theme_classic(base_size = 8) +
      theme(plot.title = element_text(size = 9, hjust = 0.5))
    if (!is.null(reference)) {
      panel <- panel + geom_vline(xintercept = reference,
        linetype = "dashed", linewidth = 0.3)
    }
    panel
  })
  patchwork::wrap_plots(panels[[1]], panels[[2]], panels[[3]], panels[[4]],
    ncol = 4, axis_titles = "collect_x")
}

# Twin layout: one panel; "Control (Intercept)" + formatted rows colored
# tab10; mu diamonds above each row (circular wrap) and kappa triangles
# below; dashed zero line; mu/kappa legend at the lower right.
base_coef_plot <- function(coef_pt, coef_lab) {
  # Explicit x limits (ggplot2's default 5% multiplicative padding) give the
  # hand-drawn legend box stable data coordinates.
  cb_raw <- range(c(coef_pt$low, coef_pt$high))
  cspan_c <- diff(cb_raw)
  cband_lim <- cb_raw + c(-0.05, 0.05) * cspan_c
  ggplot(coef_pt, aes(x, y)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = low, xmax = high, colour = treatment),
      height = 0.18, linewidth = 0.3) +
    geom_point(aes(colour = treatment, shape = comp), size = 1.9) +
    scale_shape_manual(values = c(mu = 18, kappa = 17), guide = "none") +
    scale_colour_manual(values = treatment_colors[tlabels_all], guide = "none") +
    scale_y_reverse(limits = c(0.5, length(coef_lab) + 0.5),
      breaks = seq_along(coef_lab), labels = coef_lab, expand = c(0, 0)) +
    scale_x_continuous(limits = cband_lim, expand = c(0, 0)) +
    # Framed legend, drawn by hand (this ggplot2 build drops
    # legend.position = "bottomright" without an error).
    annotate("rect", xmin = cband_lim[1] + 0.55 * cspan_c,
      xmax = cband_lim[2] - 0.04 * cspan_c, ymin = 2.95, ymax = 4.35,
      fill = "white", colour = "grey40", linewidth = 0.4) +
    annotate("point", x = cband_lim[1] + 0.63 * cspan_c, y = 3.35,
      shape = 18, size = 2.6) +
    annotate("text", x = cband_lim[1] + 0.82 * cspan_c, y = 3.35,
      label = "\u03bc", hjust = 0.5, vjust = 0.5, size = 3.0) +
    annotate("point", x = cband_lim[1] + 0.63 * cspan_c, y = 3.98,
      shape = 17, size = 2.6) +
    annotate("text", x = cband_lim[1] + 0.82 * cspan_c, y = 3.98,
      label = "\u03ba", hjust = 0.5, vjust = 0.5, size = 3.0) +
    theme_classic(base_size = 9) +
    labs(x = "Coefficient value", y = "")
}

# Polar mixture panels follow the Python twin: 0 deg at the top, rotating
# clockwise; the radius is 1 + pdf/peak_max, where peak_max is the shared
# peak of every panel's 95% band upper bound, so density zero sits on the
# ring at radius 1 and the shared peak on the ring at radius 2.
# Colourblind treatment colours, 95% band at the twin's 0.35 opacity, no
# angle labels; the rings and the 0-deg spoke are drawn explicitly.
mixture_polar_plot <- function(polar_pdf) {
  peak_max_mix <- max(polar_pdf$high)
  polar_pdf$r_med <- 1 + polar_pdf$median / peak_max_mix
  polar_pdf$r_low <- 1 + polar_pdf$low / peak_max_mix
  polar_pdf$r_high <- 1 + polar_pdf$high / peak_max_mix
  rmax_mix <- 2.05
  mix_ring <- function(r) data.frame(theta = seq(0, 2 * pi, length.out = 361L), r = r)
  ggplot(polar_pdf, aes(theta, r_med, colour = treatment, group = treatment)) +
    geom_ribbon(aes(ymin = r_low, ymax = r_high, fill = treatment), alpha = 0.35, colour = NA) +
    geom_path(linewidth = 0.5) +
    geom_path(data = mix_ring(1), inherit.aes = FALSE, aes(theta, r),
      colour = "grey55", linewidth = 0.25) +
    geom_path(data = mix_ring(2), inherit.aes = FALSE, aes(theta, r),
      colour = "grey75", linewidth = 0.25) +
    geom_path(data = data.frame(theta = c(0, 0), r = c(1, rmax_mix)),
      inherit.aes = FALSE, aes(theta, r), colour = "grey70", linewidth = 0.3) +
    facet_wrap(~ treatment) +
    scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0), breaks = NULL) +
    coord_polar() +
    scale_y_continuous(limits = c(0, rmax_mix), expand = c(0, 0)) +
    scale_colour_manual(values = treatment_colors_cb) +
    scale_fill_manual(values = treatment_colors_cb) +
    labs(x = "", y = "") +
    theme(legend.position = "none", strip.text = element_blank(),
      strip.background = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
}

polar_ring360 <- function(r) {
  a <- seq(0, 2 * pi, length.out = 361L)
  data.frame(x = a, y = r)
}

# Polar population PPC, one panel per treatment (Python twin look): 0 deg
# at top, clockwise; edge ring at r = 0.55 (density 0); observed headings
# black dashed with a 95% HDI band; posterior-predictive headings in the
# treatment's colorblind colour with a band; Data/Model legend at the
# bottom. The histogram rows stay in the table; the figure shows KDE
# bands only, as in the Python twin.
population_ppc_plot <- function(population_ppc) {
  ppc_plot <- population_ppc[population_ppc$kind == "KDE", , drop = FALSE]
  ppc_plot$key <- ifelse(ppc_plot$source == "observed bootstrap", "Data", "Model")
  # Facet titles mirror the Python twin's format_treatment_title: comma
  # parts stay in source order, Control is dropped from mixed names, and
  # every word is initial-cased.
  treat_title_part <- function(pp) {
    ws <- strsplit(tolower(pp), " ", fixed = TRUE)[[1]]
    ws <- ws[ws != ""]
    paste0(toupper(substring(ws, 1L, 1L)), substring(ws, 2L), collapse = " ")
  }
  treat_title <- function(code) {
    nm <- treatment_labels[as.character(code)]
    if (is.na(nm)) nm <- as.character(code)
    nm <- trimws(nm)
    if (!nzchar(nm)) return(as.character(code))
    parts <- trimws(unlist(strsplit(nm, ",", fixed = TRUE)))
    parts <- parts[nzchar(parts)]
    if (length(parts) > 1L) parts <- parts[tolower(parts) != "control"]
    if (!length(parts)) return(as.character(code))
    paste(vapply(parts, treat_title_part, character(1)), collapse = ", ")
  }
  tlabel_vec <- vapply(as.character(levels(df_track$treatment)), treat_title,
    character(1))
  ppc_plot$tlabel <- factor(tlabel_vec[as.character(ppc_plot$treatment)],
    levels = tlabel_vec)
  edge_ppc <- 0.55
  amp_ppc <- 0.5
  rmax_ppc <- edge_ppc + amp_ppc * 1.15
  ppc_plot$ymed <- pmin(edge_ppc + ppc_plot$median, rmax_ppc)
  ppc_plot$ylow <- pmin(edge_ppc + ppc_plot$low, rmax_ppc)
  ppc_plot$yhi <- pmin(edge_ppc + ppc_plot$high, rmax_ppc)
  # The model rows carry their treatment's colorblind colour through the
  # treatment code; the data rows carry "Data" (black). The colour and
  # fill maps are therefore named identically to the data, stay hidden
  # from the legend, and the legend shows the Data/Model line styles.
  ppc_plot$cmap <- ifelse(ppc_plot$key == "Data", "Data",
    as.character(ppc_plot$treatment))
  tcol_pop <- c(Data = "black", setNames(treatment_colors_cb[
    levels(df_track$treatment)], as.character(levels(df_track$treatment))))
  ppc_obs <- ppc_plot[ppc_plot$source == "observed bootstrap", , drop = FALSE]
  ppc_mod <- ppc_plot[ppc_plot$source == "posterior predictive", , drop = FALSE]
  ggplot(ppc_plot, aes(x, ymed)) +
    geom_ribbon(aes(ymin = ylow, ymax = yhi, group = interaction(key, tlabel),
      fill = cmap), alpha = 0.25, colour = NA) +
    geom_line(data = ppc_obs, aes(group = interaction(key, tlabel),
      linetype = key), colour = "black", linewidth = 0.3) +
    geom_line(data = ppc_mod, aes(group = interaction(key, tlabel),
      linetype = key, colour = cmap), linewidth = 0.3) +
    geom_path(data = polar_ring360(edge_ppc), inherit.aes = FALSE, aes(x, y),
      colour = "#595959", linewidth = 0.25) +
    geom_path(data = polar_ring360(edge_ppc + amp_ppc), inherit.aes = FALSE, aes(x, y),
      colour = "#d9d9d9", linewidth = 0.25) +
    geom_path(data = data.frame(x = c(0, 0), y = c(edge_ppc, edge_ppc + amp_ppc)),
      inherit.aes = FALSE, aes(x, y), colour = "#d9d9d9", linewidth = 0.3) +
    annotate("text", x = 0, y = edge_ppc + amp_ppc + 0.035,
      label = "0\u00b0", size = 2.0) +
    facet_wrap(~ tlabel, ncol = 4) +
    scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0),
      breaks = 0, labels = "0°") +
    coord_polar() +
    scale_y_continuous(limits = c(0, rmax_ppc), expand = c(0, 0)) +
    scale_colour_manual(values = tcol_pop, name = NULL, guide = "none") +
    scale_fill_manual(values = tcol_pop, name = NULL, guide = "none") +
    scale_linetype_manual(values = c(Data = "dashed", Model = "solid"), name = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6, colour = "grey25"),
      panel.grid = element_blank(), panel.border = element_blank(),
      strip.background = element_blank(), strip.text = element_text(size = 9),
      legend.position = "bottom", legend.background = element_rect(
        fill = "white", colour = "grey60"),
      legend.key = element_blank(), legend.text = element_text(size = 8)) +
    labs(x = "", y = "")
}

# Per-tag polar PPC grid (Python twin look): 8 columns, 0 deg at top and
# clockwise, only the 0-deg label, two grey rings and a 0-deg spoke;
# observed bootstrap black dashed + grey band, posterior-predictive in the
# tag's treatment colour (colorblind) + band; bold tag titles; legend with
# Data/Model plus the named treatments.
individual_ppc_plot <- function(individual_ppc) {
  ind_edge <- 0.55
  ind_rmax <- ind_edge + 1.01
  ind <- individual_ppc
  # Twin panel order: treatment blocks (sorted treatment code), tags
  # ascending within each block; the stable-sorted twin order is the one the
  # grid reads off.
  ind$taglabel <- as.character(ind$tag)
  tags_by_treatment <- unlist(Map(
    function(t) {
      z <- as.character(tag_grid$tag[as.character(tag_grid$treatment) == t])
      z[order(as.numeric(z))]
    },
    levels(df_track$treatment)))
  ind$taglabel <- factor(ind$taglabel,
    levels = intersect(unique(tags_by_treatment), unique(ind$taglabel)))
  titlecase_t <- function(s) {
    parts <- strsplit(tolower(s), " ", fixed = TRUE)[[1]]
    paste0(toupper(substring(parts, 1L, 1L)), substring(parts, 2L),
      collapse = " ")
  }
  # The factor level names double as the display labels: manual scales
  # match the data by level name, so the values and the scale keys must
  # be byte-identical.
  key_treat <- vapply(tnames, titlecase_t, character(1))
  ind$key <- ifelse(ind$source == "observed bootstrap", "Data",
    key_treat[match(as.character(ind$treatment), levels(df_track$treatment))])
  ind$key <- factor(ind$key, levels = c("Data", "Model", key_treat))
  ind$fcol <- ifelse(ind$source == "observed bootstrap", "grey20",
    treatment_colors_cb[as.character(ind$treatment)])
  # Peak normalization follows the Python twin: radius = 0.55 + 1 *
  # d/peak_max, where peak_max is the maximum of all per-tag 95% band
  # peaks, so the twin's 1.56 outer limit is preserved.
  peak_max_ind <- max(ind$high)
  ind$ymed <- pmin(ind_edge + ind$median / peak_max_ind, ind_rmax)
  ind$ylow <- pmin(ind_edge + ind$low / peak_max_ind, ind_rmax)
  ind$yhi <- pmin(ind_edge + ind$high / peak_max_ind, ind_rmax)
  ind_obs <- ind[ind$source == "observed bootstrap", , drop = FALSE]
  ind_ppc <- ind[ind$source == "posterior predictive", , drop = FALSE]
  p_individual_ppc <- ggplot(ind, aes(x, ymed)) +
    geom_ribbon(data = ind_obs, aes(ymin = ylow, ymax = yhi, group = tag),
      fill = "grey20", alpha = 0.18, colour = NA) +
    geom_ribbon(data = ind_ppc, aes(ymin = ylow, ymax = yhi, group = tag),
      fill = ind_ppc$fcol, alpha = 0.28, colour = NA) +
    geom_line(data = ind_obs, aes(group = tag), colour = "grey5",
      linetype = "dashed", linewidth = 0.3) +
    geom_line(data = ind_ppc, aes(group = tag, colour = key, linetype = key),
      linewidth = 0.4) +
    geom_path(data = polar_ring360(ind_edge), inherit.aes = FALSE, aes(x, y),
      colour = "#cccccc", linewidth = 0.2) +
    geom_path(data = polar_ring360(ind_edge + 0.5), inherit.aes = FALSE, aes(x, y),
      colour = "#dbdbdb", linewidth = 0.2) +
    geom_path(data = data.frame(x = c(0, 0), y = c(ind_edge, ind_rmax)),
      inherit.aes = FALSE, aes(x, y), colour = "#d9d9d9", linewidth = 0.3) +
    facet_wrap(~ taglabel, ncol = 8) +
    scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0),
      breaks = 0, labels = "0°") +
    coord_polar() +
    scale_y_continuous(limits = c(0, ind_rmax), expand = c(0, 0)) +
    scale_colour_manual(values = c(Data = "black", Model = "black",
      setNames(treatment_colors_cb[levels(df_track$treatment)], key_treat)), name = NULL) +
    scale_linetype_manual(values = c(Data = "dashed", Model = "solid",
      setNames(rep("solid", length(tnames)), key_treat))) +
    theme_minimal(base_size = 9) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6, colour = "grey25"),
      panel.grid = element_blank(), panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      legend.position = "none") +
    labs(x = "", y = "")

  # One framed legend row under the grid (the twin of the shared Python
  # fig.legend): Data (dashed), Model (solid) and the named treatments in
  # one row, sized against the 6.9 in page width.
  # Legend labels mirror the Python twin's format_treatment_nicely:
  # every word initial-cased, empty names read as Control.
  leg_treat <- vapply(tnames, function(nm) {
    nm <- trimws(as.character(nm))
    if (!nzchar(nm)) return("Control")
    titlecase_t(nm)
  }, character(1))
  leg_items <- data.frame(
    label = c("Data", "Model", leg_treat),
    colour = c("black", "black",
      as.character(treatment_colors_cb[levels(df_track$treatment)])),
    dashed = c(TRUE, FALSE, rep(FALSE, length(tnames))),
    stringsAsFactors = FALSE)
  # Hand-drawn one-row legend (the Python twin's fig.legend: six framed
  # keys). A one-row gtable with explicit inch units has a deterministic
  # width that centres itself in the full-width legend row; a gTree of
  # absolute-inch children does not (its extent is viewport-dependent),
  # which clipped the last key on the narrower 67-tag grid.
  tw_in <- strwidth(leg_items$label, units = "inches") * (8 / 12)  # 8 pt glyphs
  lbl_pad <- 0.04  # room for glyph anti-aliasing past the measured width
  key_in <- 0.18; gkt_in <- 0.05; gi_in <- 0.10; leg_pad <- 0.08
  row_in <- sum(tw_in) + nrow(leg_items) * (key_in + gkt_in) +
    (nrow(leg_items) - 1) * gi_in + nrow(leg_items) * lbl_pad
  if (2 * leg_pad + row_in > 5.7) {
    key_in <- 0.14; gi_in <- 0.06
    row_in <- sum(tw_in) + nrow(leg_items) * (key_in + gkt_in) +
      (nrow(leg_items) - 1) * gi_in + nrow(leg_items) * lbl_pad
  }
  ws <- leg_pad
  for (i in seq_len(nrow(leg_items))) {
    ws <- c(ws, key_in, gkt_in, tw_in[i] + lbl_pad)
    if (i < nrow(leg_items)) ws <- c(ws, gi_in)
  }
  ws <- c(ws, leg_pad)
  leg_gtable <- gtable::gtable(widths = unit(ws, "in"), heights = unit(0.34, "in"))
  leg_gtable <- gtable::gtable_add_grob(leg_gtable,
    grid::rectGrob(gp = grid::gpar(fill = "white", col = "grey60", lwd = 0.3)),
    t = 1, l = 1, r = length(ws))
  for (i in seq_len(nrow(leg_items))) {
    it <- leg_items[i, ]
    ckey <- 1 + (i - 1) * 4 + 1  # columns: pad + per item (key, gap, label, gap)
    leg_gtable <- gtable::gtable_add_grob(leg_gtable,
      grid::segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"), default.units = "npc",
        gp = grid::gpar(col = it$colour, lwd = 1.1,
          lty = if (it$dashed) "dashed" else "solid")), t = 1, l = ckey)
    leg_gtable <- gtable::gtable_add_grob(leg_gtable,
      grid::textGrob(it$label, x = unit(0, "npc"), y = unit(0.5, "npc"),
        just = "left", gp = grid::gpar(fontsize = 8)), t = 1, l = ckey + 2)
  }
  base_ppc_grob <- ggplotGrob(p_individual_ppc)
  base_rows <- max(base_ppc_grob$layout$b)
  gt_ind <- gtable::gtable_add_rows(base_ppc_grob,
    heights = unit(0.5, "in"), pos = base_rows + 1)
  gt_ind <- gtable::gtable_add_grob(gt_ind, leg_gtable, t = base_rows + 1,
    l = 1, r = max(base_ppc_grob$layout$r))
  gt_ind
}

directedness_map <- function(segs, bm, xlim, ylim, r_roost, r_trans) {
  panel_pg <- function(tt, border_col) {
    sd <- segs[as.character(segs$treatment) == tt, , drop = FALSE]
    # Border as a data-level rect added LAST: the flat basemap cover would
    # otherwise swallow half of a theme-drawn panel border.
    p <- map_panel(bm, aspect = 1.67 / 5.58,
      lake = TRUE, lake_alpha = 0.25, lake_lw = 1.0)
    p +
      geom_segment(data = sd[sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = p),
        linetype = "dashed", linewidth = 0.95 * 0.3528, alpha = 0.80) +
      geom_segment(data = sd[!sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = p),
        linewidth = 1.0 * 0.3528, alpha = 0.92) +
      scale_colour_viridis_c(limits = c(0, 1)) +
      guides(colour = "none") +
      add_roost_marker(r_roost) +
      add_x_markers(TRANS_XY, r_trans) +
      geom_rect(inherit.aes = FALSE,
        data = data.frame(xmin = xlim[1], xmax = xlim[2],
          ymin = ylim[1], ymax = ylim[2]),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = NA, colour = border_col, linewidth = 0.78)
  }
  p_pg_1 <- panel_pg("1", treatment_colors_cb["1"])
  p_pg_2 <- panel_pg("2", treatment_colors_cb["2"])
  p_pg_3 <- panel_pg("3", treatment_colors_cb["3"])
  p_pg_4 <- panel_pg("4", treatment_colors_cb["4"])
  p_pg_4 <- add_scale_compass(p_pg_4, xlim, ylim,
    variant = "pgoal", panel_w_in = 1.48, panel_h_in = 4.95)
  # Thin vertical colourbar (Python: figure axes [0.922, 0.20, 0.012,
  # 0.62], label EV_LABEL).
  pb <- seq_len(240)
  cbar_rows <- data.frame(x = 0.30, y = 0.19 + 0.62 * (pb - 0.5) / 240,
    v = (pb - 0.5) / 240)
  p_cbar <- ggplot(cbar_rows, aes(x, y, fill = v)) +
    geom_tile(width = 0.13, height = 0.62 / 240) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1.35), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  for (v in seq(0, 1, by = 0.2)) {
    yv <- 0.19 + 0.62 * v
    p_cbar <- p_cbar +
      annotate("segment", x = 0.375, y = yv, xend = 0.40, yend = yv,
        linewidth = 0.15, colour = "#000000") +
      annotate("text", x = 0.42, y = yv, label = sprintf("%.1f", v),
        size = 2.2, hjust = 0, vjust = 0.5)
  }
  p_cbar <- p_cbar +
    annotate("text", x = 1.10, y = 0.5,
      label = EV_LABEL, size = 2.8, angle = 90,
      hjust = 0.5, vjust = 0.5) +
    theme_void() +
    guides(fill = "none")
  p_pg_legend <- legend_panel(data.frame(
    kind = c("swatch", "star", "x", "dash"),
    label = c("Lake", "Roost", "Translocation", "Signal Dropout")),
    cell_w = 7.4, cell_h = 0.84)
  p_prob_top <- (p_pg_1 | p_pg_2 | p_pg_3 | p_pg_4) / p_pg_legend
  p_prob_top <- p_prob_top + plot_layout(heights = c(1, 0.17))
  p_probability_map <- p_prob_top | p_cbar
  p_probability_map +
    plot_layout(widths = c(1, 0.09)) +
    theme(plot.margin = margin(1, 1, 1, 1))
}

directedness_map_compact <- function(segs, bm, xlim, ylim, r_roost, r_trans) {
  panel_pg_s <- function(tt, border_col) {
    sd <- segs[as.character(segs$treatment) == tt, , drop = FALSE]
    # Border as a data-level rect added LAST (see panel_pg above).
    map_panel(bm, aspect = 1.34 / 3.59,
      lake = TRUE, lake_alpha = 0.25, lake_lw = 0.8) +
      geom_segment(data = sd[sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = p),
        linetype = "dashed", linewidth = 0.76 * 0.3528, alpha = 0.75) +
      geom_segment(data = sd[!sd$gap, ], inherit.aes = FALSE,
        aes(x = x, y = y, xend = xend, yend = yend, colour = p),
        linewidth = 0.8 * 0.3528, alpha = 0.90) +
      scale_colour_viridis_c(limits = c(0, 1)) +
      guides(colour = "none") +
      add_roost_marker(r_roost) +
      add_x_markers(TRANS_XY, r_trans) +
      geom_rect(inherit.aes = FALSE,
        data = data.frame(xmin = xlim[1], xmax = xlim[2],
          ymin = ylim[1], ymax = ylim[2]),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = NA, colour = border_col, linewidth = 0.53)
  }
  p_pgs_1 <- panel_pg_s("1", treatment_colors_cb["1"])
  p_pgs_2 <- panel_pg_s("2", treatment_colors_cb["2"])
  p_pgs_3 <- panel_pg_s("3", treatment_colors_cb["3"])
  p_pgs_4 <- panel_pg_s("4", treatment_colors_cb["4"])
  p_pgs_4 <- add_scale_compass(p_pgs_4, xlim, ylim,
    variant = "pgoal_small", panel_w_in = 1.33, panel_h_in = 3.49)
  # Wide horizontal colourbar (Python: figure axes [0.20, 0.060, 0.60,
  # 0.034]).
  pb <- seq_len(240)
  cbar_h_rows <- data.frame(x = 0.20 + 0.60 * (pb - 0.5) / 240, y = 0.62,
    v = (pb - 0.5) / 240)
  p_cbar_h <- ggplot(cbar_h_rows, aes(x, y, fill = v)) +
    geom_tile(width = 0.60 / 240, height = 0.24) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  for (v in seq(0, 1, by = 0.2)) {
    xv <- 0.20 + 0.60 * v
    p_cbar_h <- p_cbar_h +
      annotate("segment", x = xv, y = 0.49, xend = xv, yend = 0.44,
        linewidth = 0.15, colour = "#000000") +
      annotate("text", x = xv, y = 0.40, label = sprintf("%.1f", v),
        size = 2.2, hjust = 0.5, vjust = 1)
  }
  p_cbar_h <- p_cbar_h +
    annotate("text", x = 0.5, y = 0.05,
      label = EV_LABEL, size = 2.8, hjust = 0.5,
      vjust = 0) +
    theme_void() +
    guides(fill = "none")
  p_prob_top_s <- (p_pgs_1 | p_pgs_2 | p_pgs_3 | p_pgs_4) / p_cbar_h
  p_prob_top_s <- p_prob_top_s + plot_layout(heights = c(1, 0.15))
  p_prob_top_s +
    theme(plot.margin = margin(1, 1, 1, 1))
}

# Full-size twin (Python plot_population_curve_smooth): colorblind bands
# + lines, dashed 0.5 reference, boxed right legend with formatted
# treatment names, y on [0, 1].
curve_plot_full <- function(curves, reverse = FALSE) {
  x <- curves$population
  x$treatment <- factor(x$treatment, levels = 1:4)
  lv <- levels(x$treatment)
  p <- ggplot(x, aes(x, mean, colour = treatment, fill = treatment,
    group = treatment)) +
    geom_ribbon(aes(ymin = low, ymax = high),
      alpha = 0.20, colour = NA, linewidth = 0) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0.5, colour = "grey90", linetype = "dashed", linewidth = 0.3) +
    scale_colour_manual(values = treatment_colors_cb[lv],
      labels = treatment_nice[lv], name = NULL) +
    scale_fill_manual(values = treatment_colors_cb[lv], name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      labels = function(v) sprintf("%.1f", v), expand = c(0, 0.02)) +
    guides(fill = "none") +
    labs(x = if (reverse) "Distance to destination (m)" else "Cumulative path distance from end (m)",
      y = "Mean posterior directedness") +
    theme_classic(base_size = 9) +
    theme(legend.position = "right",
      legend.background = element_rect(fill = "white", colour = "grey30"),
      legend.key = element_blank(),
      legend.text = element_text(size = 8, colour = "black"),
      legend.margin = margin(2, 2, 2, 2, "mm"),
      legend.spacing.y = unit(1.3, "lines"))
  if (reverse) {
    p <- p + scale_x_reverse(labels = function(v) format(v, trim = TRUE, scientific = FALSE))
  } else {
    p <- p + scale_x_continuous(limits = c(-110 * 60, 0),
      labels = function(v) format(abs(v), big.mark = ","))
  }
  p
}

# Compact (paper-width) twin: same curves, no legend, smaller text.
curve_plot_compact <- function(curves, reverse = FALSE) {
  x <- curves$population
  x$treatment <- factor(x$treatment, levels = 1:4)
  lv <- levels(x$treatment)
  p <- ggplot(x, aes(x, mean, colour = treatment, fill = treatment,
    group = treatment)) +
    geom_ribbon(aes(ymin = low, ymax = high),
      alpha = 0.18, colour = NA, linewidth = 0) +
    geom_line(linewidth = 0.45) +
    geom_hline(yintercept = 0.5, colour = "grey88", linetype = "dashed", linewidth = 0.3) +
    scale_colour_manual(values = treatment_colors_cb[lv], name = NULL) +
    scale_fill_manual(values = treatment_colors_cb[lv], name = NULL) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
    guides(colour = "none", fill = "none") +
    labs(x = if (reverse) "Distance to destination (m)" else "Cumulative path distance from end (m)",
      y = "Posterior\ndirectedness") +
    theme_classic(base_size = 7) +
    theme(legend.position = "none")
  if (reverse) {
    p <- p + scale_x_reverse(labels = function(v) format(v, trim = TRUE, scientific = FALSE))
  } else {
    p <- p + scale_x_continuous(limits = c(-110 * 60, 0),
      labels = function(v) format(abs(v), big.mark = ","))
  }
  p
}

# Per-animal 2x2 grid (Python plot_tag_curves_faceted_4panel): one panel
# per treatment, animals drawn in the treatment color, 0.5 reference,
# absolute comma-formatted ticks, y label only on the left column.
tag_curve_plot <- function(curves, reverse = FALSE) {
  x <- curves$individual
  x$treatment <- factor(x$treatment, levels = 1:4)
  lv <- levels(x$treatment)
  xlabs <- function(v) format(abs(v), big.mark = ",")
  panels <- lapply(seq_along(lv), function(k) {
    t <- lv[k]
    sub <- x[x$treatment == t, ]
    p <- ggplot(sub, aes(x, mean, group = tag)) +
      geom_ribbon(aes(ymin = low, ymax = high),
        fill = treatment_colors_cb[t], alpha = 0.10, colour = NA, linewidth = 0) +
      geom_line(colour = treatment_colors_cb[t], alpha = 0.75, linewidth = 0.35) +
      geom_hline(yintercept = 0.5, colour = "grey90", linetype = "dashed", linewidth = 0.3) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
        labels = function(v) sprintf("%.1f", v), expand = c(0, 0)) +
      labs(title = treatment_nice[[t]],
        x = if (reverse) "Distance to destination (m)" else "Cumulative path distance from end (m)",
        y = if (k == 1L || k == 3L) "Posterior directedness" else "") +
      theme_classic(base_size = 8) +
      theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.text.y = if (k %% 2L == 0L) element_blank() else element_text())
    if (reverse) {
      p <- p + scale_x_reverse(labels = function(v) format(v, trim = TRUE, scientific = FALSE))
    } else {
      p <- p + scale_x_continuous(limits = c(-110 * 60, 0), labels = xlabs)
    }
    p
  })
  patchwork::wrap_plots(panels[[1]], panels[[2]], panels[[3]], panels[[4]],
    ncol = 2)
}
