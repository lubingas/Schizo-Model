# app/charts.R
# ---------------------------
# Charts for the Shiny app (plotly). One fixed colour per drug; hairline grids;
# every chart has a table or direct labels so colour never carries meaning alone.
# ---------------------------

axis_style <- function(title = "", ...) {
  list(title = list(text = title, font = list(size = 13, color = ui_ink$secondary), standoff = 12),
       gridcolor = ui_ink$grid, gridwidth = 1, zerolinecolor = ui_ink$axis, zerolinewidth = 1,
       linecolor = ui_ink$axis, tickfont = list(size = 12, color = ui_ink$secondary),
       automargin = TRUE, ...)
}

style_chart <- function(p, filename = "chart", legend = FALSE, ...) {
  p %>%
    layout(font = list(family = app_font, size = 13, color = ui_ink$secondary),
           paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
           hoverlabel = list(bgcolor = "#fcfcfb", bordercolor = "#c5ccc8",
                             font = list(family = app_font, size = 13, color = ui_ink$primary), align = "left"),
           showlegend = legend, separators = ".,", ...) %>%
    config(displaylogo = FALSE, responsive = TRUE,
           modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "hoverClosestCartesian",
                                      "hoverCompareCartesian", "toggleSpikelines"),
           toImageButtonOptions = list(format = "png", scale = 2, filename = filename))
}

pad_range <- function(x, frac = 0.12, include = NULL) {
  r <- range(c(x, include), na.rm = TRUE)
  d <- diff(r); if (d == 0) d <- max(abs(r), 1) * 0.2
  c(r[1] - frac * d, r[2] + frac * d)
}

# ── Base case: cost-effectiveness plane with the efficiency frontier ───────────
plane_chart <- function(model, width = 1000) {
  ia  <- IncrementalAnalysis(model)
  pts <- ia$ce_points
  fr  <- ia$ce_lines
  ref <- pts$intervention[1]          # the drug with the most DALYs

  p <- plot_ly()
  if (nrow(fr) > 1) {
    p <- add_trace(p, x = fr$delta_outcomes, y = fr$delta_costs, type = "scatter", mode = "lines",
                   line = list(color = ui_ink$primary, width = 2), hoverinfo = "skip", showlegend = FALSE)
  }
  # label offsets: alternate below/above in order of cost so close points don't collide
  ord <- rank(pts$delta_costs, ties.method = "first")
  ay  <- (if (width < 600) c(40, -40, 84, -84, 40) else c(34, -34, 64, -64, 34))[ord]
  ann <- list()
  for (i in seq_len(nrow(pts))) {
    d <- pts$intervention[i]; on <- pts$status[i] == "On frontier"
    p <- add_trace(p, x = pts$delta_outcomes[i], y = pts$delta_costs[i], type = "scatter", mode = "markers",
                   marker = list(size = 15, color = if (on) drug_cols[[d]] else ui_ink$surface,
                                 line = list(color = drug_cols[[d]], width = 3)),
                   hovertemplate = paste0("<b>", d, "</b><br>Cost per patient: ", fmt_usd(pts$costs[i]),
                                          "<br>DALYs per patient: ", fmt_num(pts$outcomes[i], 3),
                                          "<br>", pts$status[i], "<extra></extra>"),
                   showlegend = FALSE)
    ann[[i]] <- list(x = pts$delta_outcomes[i], y = pts$delta_costs[i], text = d, showarrow = TRUE,
                     arrowhead = 0, arrowwidth = 1, arrowcolor = ui_ink$axis, ax = 0, ay = ay[i],
                     font = list(size = 13, color = ui_ink$primary), standoff = 9)
  }
  p %>%
    layout(xaxis = axis_style(if (width < 600) paste("DALYs averted vs", tolower(ref)) else paste("DALYs averted per patient, compared with", tolower(ref)),
                              range = { r <- pad_range(pts$delta_outcomes, 0.15); if (width < 600) r[1] <- r[1] - 0.25 * diff(r); r },
                              tickformat = ".3f", nticks = if (width < 600) 4 else 8),
           yaxis = axis_style("Incremental cost per patient (US$)", range = pad_range(pts$delta_costs, 0.3),
                              tickformat = "$,.0f"),
           annotations = ann, margin = list(t = 10, r = 20, b = 10, l = 10)) %>%
    style_chart("cost-effectiveness-plane")
}

# ── One-way analysis: tornado diagram ──────────────────────────────────────────
owsa_value_fmt <- function(outcome) {
  switch(outcome,
         cost = function(x) fmt_usd(x, 2),
         outcome = function(x) fmt_num(x, 4),
         ICER = function(x) paste0(fmt_usd(x), " per DALY averted"),
         NMB = function(x) fmt_usd(x, 2))
}
owsa_axis_title <- function(outcome, drug, comp) {
  switch(outcome,
         cost = sprintf("Incremental cost per patient, %s minus %s (US$)", tolower(drug), tolower(comp)),
         outcome = sprintf("DALYs averted per patient by %s, compared with %s", tolower(drug), tolower(comp)),
         ICER = sprintf("Cost per DALY averted, %s compared with %s (US$)", tolower(drug), tolower(comp)),
         NMB = sprintf("Net monetary benefit per patient, %s compared with %s (US$)", tolower(drug), tolower(comp)))
}
# what a run with no cost per DALY averted means, for hover text and the table
owsa_status_text <- c(dominant = "Dominant: costs less and averts at least as many DALYs",
                      dominated = "Dominated: costs more and averts no more DALYs",
                      "cheaper, less effective" = "Costs less but averts fewer DALYs")

# ow: output of owsaPair(); info: table with code, label, low_in, high_in for every row
tornado_data <- function(ow, info, n, max_chars = 46) {
  d <- ow$owsa_data                # sorted by effect, largest last
  d$code <- sub(" \\(.*$", "", d$variable)
  d <- dplyr::left_join(d, info, by = "code")
  d <- d[rev(seq_len(nrow(d))), ][seq_len(min(n, nrow(d))), ]
  lab <- d$label
  d$label_short <- ifelse(nchar(lab) > max_chars, paste0(substr(lab, 1, max_chars - 2), "…"), lab)
  d$label_short <- make.unique(d$label_short, sep = " ")
  d[rev(seq_len(nrow(d))), ]       # largest effect at the top
}
tornado_chart <- function(ow, info, n, outcome, drug, comp, width = 1000) {
  narrow <- width < 600
  d <- tornado_data(ow, info, n, if (narrow) 24 else 46)
  base <- ow$owsa_basecase
  f <- owsa_value_fmt(outcome)
  col <- c(low = "#86b6ef", high = "#1c5cab")
  hov <- function(which, inval, res) paste0("<b>", d$label, "</b><br>", which, " input: ", inval,
                                            "<br>Result: ", f(res), "<extra></extra>")
  # draw each row's longer bar first, so the shorter one stays visible when both fall on the same side
  low_behind <- is.na(d$high) | (!is.na(d$low) & abs(d$low - base) >= abs(d$high - base))
  layer <- function(is_low) list(x = as.numeric(ifelse(is_low, d$low, d$high)) - base,
                                 col = ifelse(is_low, col[["low"]], col[["high"]]),
                                 hover = ifelse(is_low, hov("Low", d$low_in, d$low), hov("High", d$high_in, d$high)))
  back <- layer(low_behind); front <- layer(!low_behind)

  p <- plot_ly(y = d$label_short) %>%
    add_trace(type = "bar", orientation = "h", x = back$x, base = base, marker = list(color = back$col),
              width = 0.62, hovertemplate = back$hover, showlegend = FALSE) %>%
    add_trace(type = "bar", orientation = "h", x = front$x, base = base, marker = list(color = front$col),
              width = 0.62, hovertemplate = front$hover, showlegend = FALSE) %>%
    # legend entries only (zero-length bars)
    add_trace(type = "bar", orientation = "h", x = 0, y = d$label_short[1], base = base, hoverinfo = "skip",
              name = "Low value of the input", marker = list(color = col[["low"]])) %>%
    add_trace(type = "bar", orientation = "h", x = 0, y = d$label_short[1], base = base, hoverinfo = "skip",
              name = "High value of the input", marker = list(color = col[["high"]]))

  # cost per DALY averted only: runs with no ratio are marked in a strip at each edge,
  # dominant on the left and dominated (or cheaper but less effective) on the right
  xaxis_extra <- list(); shapes <- list(); notes <- list()
  if (outcome == "ICER") {
    xr <- range(c(base, d$low, d$high), na.rm = TRUE)
    strip <- 0.3 * max(diff(xr), abs(base))
    lim <- c(xr[1] - 1.4 * strip, xr[2] + 1.4 * strip)
    ticks <- pretty(xr, n = if (narrow) 3 else 6); ticks <- ticks[ticks >= xr[1] - 0.2 * strip & ticks <= xr[2] + 0.2 * strip]
    xaxis_extra <- list(range = lim, tickvals = ticks, ticktext = fmt_usd(ticks))
    marks <- function(side_low, status, x_at, symbol) {
      keep <- (if (side_low) d$low_status else d$high_status) %in% status
      if (!any(keep)) return(NULL)
      st <- (if (side_low) d$low_status else d$high_status)[keep]
      list(x = rep(x_at, sum(keep)), y = d$label_short[keep], col = col[[if (side_low) "low" else "high"]],
           symbol = ifelse(st == "cheaper, less effective", paste0(symbol, "-open"), symbol),
           hover = paste0("<b>", d$label[keep], "</b><br>", if (side_low) "Low" else "High", " input: ",
                          (if (side_low) d$low_in else d$high_in)[keep], "<br>", owsa_status_text[st], "<extra></extra>"))
    }
    gap <- 0.22 * strip
    for (m in list(marks(TRUE, "dominant", xr[1] - 0.7 * strip - gap, "triangle-left"),
                   marks(FALSE, "dominant", xr[1] - 0.7 * strip + gap, "triangle-left"),
                   marks(TRUE, c("dominated", "cheaper, less effective"), xr[2] + 0.7 * strip - gap, "triangle-right"),
                   marks(FALSE, c("dominated", "cheaper, less effective"), xr[2] + 0.7 * strip + gap, "triangle-right"))) {
      if (is.null(m)) next
      p <- p %>% add_trace(type = "scatter", mode = "markers", x = m$x, y = m$y, showlegend = FALSE,
                           marker = list(color = m$col, size = 13, symbol = m$symbol, line = list(color = m$col, width = 2)),
                           hovertemplate = m$hover)
    }
    edge <- function(x) list(type = "line", xref = "x", yref = "paper", x0 = x, x1 = x, y0 = 0, y1 = 1,
                             line = list(color = ui_ink$axis, width = 1, dash = "dot"))
    shapes <- list(edge(xr[1] - 0.35 * strip), edge(xr[2] + 0.35 * strip))
    # strip names sit in the row of tick labels, under each strip
    head <- function(x, text) list(x = x, y = 0, xref = "x", yref = "paper", yanchor = "top", yshift = -6, showarrow = FALSE,
                                   text = text, font = list(size = 12, color = ui_ink$secondary))
    notes <- list(head(xr[1] - 0.7 * strip, "Dominant"), head(xr[2] + 0.7 * strip, "Dominated"))
  }

  p %>%
    layout(barmode = "overlay", bargap = 0.3,
           xaxis = do.call(axis_style, c(list(if (narrow) switch(outcome, cost = "Incremental cost (US$)", outcome = "DALYs averted",
                                                                  ICER = "Cost per DALY averted (US$)", NMB = "Net monetary benefit (US$)")
                                              else owsa_axis_title(outcome, drug, comp),
                                              tickformat = if (outcome == "outcome") ".4f" else "$,.2~f", zeroline = TRUE,
                                              nticks = if (narrow) 4 else 8, tickangle = 0), xaxis_extra)),
           yaxis = axis_style("", categoryorder = "array", categoryarray = d$label_short,
                              tickfont = list(size = if (narrow) 11.5 else 12.5, color = ui_ink$primary), gridcolor = "rgba(0,0,0,0)"),
           shapes = c(list(list(type = "line", xref = "x", yref = "paper", x0 = base, x1 = base, y0 = 0, y1 = 1,
                                line = list(color = ui_ink$primary, width = 1.5))), shapes),
           annotations = notes,
           legend = list(orientation = "h", x = 0, y = 1, yanchor = "bottom", xanchor = "left", itemclick = FALSE,
                         itemdoubleclick = FALSE, traceorder = "normal", font = list(size = 12.5, color = ui_ink$secondary)),
           margin = list(t = 36, r = if (narrow) 8 else 20, b = 10, l = 6)) %>%
    style_chart("tornado", legend = TRUE)
}

# ── Probabilistic analysis: acceptability curves ───────────────────────────────
ceac_chart <- function(ce, wtp) {
  p <- plot_ly()
  for (d in drug_order) {
    p <- add_trace(p, x = ce$WTP, y = ce[[d]], type = "scatter", mode = "lines", name = d,
                   line = list(color = drug_cols[[d]], width = 2.5, shape = "linear"),
                   hovertemplate = paste0(d, ": %{y:.0%}<extra></extra>"))
  }
  p %>%
    layout(hovermode = "x unified",
           xaxis = axis_style("Willingness to pay per DALY averted (US$)", tickformat = "$,.0f",
                              hoverformat = ",.0f", range = c(0, max(ce$WTP))),
           yaxis = axis_style("Probability of being the most cost-effective drug", tickformat = ".0%", range = c(0, 1.02)),
           shapes = list(list(type = "line", xref = "x", yref = "paper", x0 = wtp, x1 = wtp, y0 = 0, y1 = 1,
                              line = list(color = ui_ink$secondary, width = 1, dash = "dot"))),
           annotations = list(list(x = wtp, y = 0.97, xref = "x", yref = "paper", yanchor = "top", xanchor = "left",
                                   text = paste0("threshold ", fmt_usd(wtp)), showarrow = FALSE, xshift = 4,
                                   font = list(size = 12, color = ui_ink$secondary), bgcolor = "rgba(252,252,251,0.85)")),
           legend = list(orientation = "h", x = 0, y = 1.02, yanchor = "bottom", xanchor = "left",
                         font = list(size = 13, color = ui_ink$primary)),
           margin = list(t = 40, r = 20, b = 10, l = 10)) %>%
    style_chart("acceptability-curves", legend = TRUE)
}

# ── Probabilistic analysis: one cost-effectiveness plane per drug ──────────────
psa_plane_chart <- function(psa, ref, wtp) {
  others <- setdiff(drug_order, ref)
  dx <- lapply(others, function(d) psa$outcomes[[ref]] - psa$outcomes[[d]])
  dy <- lapply(others, function(d) psa$costs[[d]] - psa$costs[[ref]])
  xr <- pad_range(unlist(dx), 0.06, 0); yr <- pad_range(unlist(dy), 0.06, 0)
  wx <- c(xr[1], xr[2]) * 1.5

  panels <- lapply(seq_along(others), function(i) {
    d <- others[i]
    plot_ly() %>%
      add_trace(x = wx, y = wtp * wx, type = "scatter", mode = "lines", hoverinfo = "skip",
                line = list(color = ui_ink$muted, width = 1, dash = "dot"), showlegend = FALSE) %>%
      add_trace(x = dx[[i]], y = dy[[i]], type = "scatter", mode = "markers", showlegend = FALSE,
                marker = list(color = drug_cols[[d]], size = 6, opacity = 0.45),
                hovertemplate = paste0("<b>", d, "</b><br>DALYs averted: %{x:.3f}<br>Incremental cost: %{y:$,.0f}<extra></extra>")) %>%
      add_trace(x = mean(dx[[i]]), y = mean(dy[[i]]), type = "scatter", mode = "markers", showlegend = FALSE,
                marker = list(color = drug_cols[[d]], size = 13, line = list(color = ui_ink$surface, width = 2)),
                hovertemplate = paste0("<b>", d, ", mean</b><br>DALYs averted: %{x:.3f}<br>Incremental cost: %{y:$,.0f}<extra></extra>"))
  })
  p <- subplot(panels, nrows = 2, shareX = TRUE, shareY = TRUE, margin = c(0.03, 0.03, 0.07, 0.07))
  # with shared axes, the 2 x 2 grid has x and x2 (one per column) and y and y2 (one per row)
  axn <- function(k) if (k == 1) "" else as.character(k)
  colk <- function(i) (i - 1) %% 2 + 1
  rowk <- function(i) (i - 1) %/% 2 + 1
  titles <- lapply(seq_along(others), function(i) {
    list(text = paste0("<b>", others[i], "</b>"), xref = paste0("x", axn(colk(i)), " domain"), yref = paste0("y", axn(rowk(i)), " domain"),
         x = 0.02, y = 0.98, xanchor = "left", yanchor = "top", showarrow = FALSE,
         font = list(size = 13, color = ui_ink$primary), bgcolor = "rgba(252,252,251,0.85)")
  })
  lay <- list(annotations = titles, margin = list(t = 10, r = 20, b = 10, l = 10))
  for (k in 1:2) {
    lay[[paste0("xaxis", axn(k))]] <- axis_style(paste("DALYs averted, compared with", tolower(ref)), range = xr, tickformat = ".2f")
    lay[[paste0("yaxis", axn(k))]] <- axis_style("Incremental cost (US$)", range = yr, tickformat = "$,.0f")
  }
  do.call(layout, c(list(p), lay)) %>% style_chart("psa-planes")
}
