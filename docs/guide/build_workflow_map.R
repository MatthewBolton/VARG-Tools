script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/")
guide_dir <- dirname(script_path)
output_file <- file.path(guide_dir, "workflow-map.png")

grDevices::png(
  output_file,
  width = 2200,
  height = 920,
  res = 200,
  bg = "white",
  type = "cairo"
)

grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "#FFFFFF", col = NA))

font_family <- "Arial"
ink <- "#1C2528"
muted <- "#526167"
panel <- "#F4F6F6"

grid::grid.text(
  "From data to a defensible decision",
  x = 0.04,
  y = 0.92,
  just = "left",
  gp = grid::gpar(fontfamily = font_family, fontsize = 23, fontface = "bold", col = ink)
)
grid::grid.text(
  "Each VARG-Tools stage produces an auditable output and a scientific checkpoint.",
  x = 0.04,
  y = 0.865,
  just = "left",
  gp = grid::gpar(fontfamily = font_family, fontsize = 13, col = muted)
)

stages <- list(
  list(
    number = "1",
    title = "PROCESS",
    colour = "#2A6F79",
    actions = c("Clean and map columns", "Transform compositions", "Cluster or project"),
    output = "Prepared analyses",
    check = "Do groups remain coherent\nin oxide space?"
  ),
  list(
    number = "2",
    title = "COMPARE",
    colour = "#C96A3D",
    actions = c("Rank variable pairs", "Plot unknowns and references", "Inspect density fields"),
    output = "Diagnostic evidence",
    check = "Which variables distinguish\nthe candidates?"
  ),
  list(
    number = "3",
    title = "CORRELATE",
    colour = "#3E6FA8",
    actions = c("Set coordinate directions", "Test candidate ties", "Inspect interval slopes"),
    output = "Accepted tie table",
    check = "Which links survive geochemistry\nand stratigraphy?"
  ),
  list(
    number = "4",
    title = "MODEL TIME",
    colour = "#6A7B3D",
    actions = c("Choose model structure", "Generate OxCal code", "Run sensitivity checks"),
    output = "Reviewed model input",
    check = "How do assumptions and ties\naffect the ages?"
  )
)

centres <- c(0.14, 0.38, 0.62, 0.86)
box_width <- 0.205
box_height <- 0.56
box_y <- 0.53

for (i in seq_along(stages)) {
  stage <- stages[[i]]
  x <- centres[[i]]

  grid::grid.roundrect(
    x = x,
    y = box_y,
    width = box_width,
    height = box_height,
    r = grid::unit(5, "pt"),
    gp = grid::gpar(fill = panel, col = "#CDD5D7", lwd = 1.2)
  )
  grid::grid.rect(
    x = x,
    y = 0.765,
    width = box_width,
    height = 0.09,
    gp = grid::gpar(fill = stage$colour, col = stage$colour)
  )
  grid::grid.text(
    paste(stage$number, stage$title),
    x = x,
    y = 0.765,
    gp = grid::gpar(fontfamily = font_family, fontsize = 14, fontface = "bold", col = "white")
  )

  action_y <- c(0.665, 0.605, 0.545)
  for (j in seq_along(stage$actions)) {
    grid::grid.points(
      x = grid::unit(x - 0.077, "npc"),
      y = grid::unit(action_y[[j]], "npc"),
      pch = 16,
      size = grid::unit(2.1, "mm"),
      gp = grid::gpar(col = stage$colour)
    )
    grid::grid.text(
      stage$actions[[j]],
      x = x - 0.064,
      y = action_y[[j]],
      just = "left",
      gp = grid::gpar(fontfamily = font_family, fontsize = 10.5, col = ink)
    )
  }

  grid::grid.text(
    "OUTPUT",
    x = x - 0.078,
    y = 0.455,
    just = "left",
    gp = grid::gpar(fontfamily = font_family, fontsize = 8.5, fontface = "bold", col = muted)
  )
  grid::grid.text(
    stage$output,
    x = x - 0.078,
    y = 0.417,
    just = "left",
    gp = grid::gpar(fontfamily = font_family, fontsize = 11, fontface = "bold", col = stage$colour)
  )

  grid::grid.rect(
    x = x,
    y = 0.315,
    width = box_width - 0.025,
    height = 0.13,
    gp = grid::gpar(fill = "white", col = "#D8DEDF", lwd = 1)
  )
  grid::grid.text(
    "CHECK",
    x = x - 0.068,
    y = 0.347,
    just = "left",
    gp = grid::gpar(fontfamily = font_family, fontsize = 8.5, fontface = "bold", col = muted)
  )
  grid::grid.text(
    stage$check,
    x = x,
    y = 0.293,
    just = "centre",
    gp = grid::gpar(fontfamily = font_family, fontsize = 9.2, col = ink)
  )
}

for (i in 1:3) {
  grid::grid.lines(
    x = grid::unit(c(centres[[i]] + box_width / 2 + 0.006, centres[[i + 1]] - box_width / 2 - 0.006), "npc"),
    y = grid::unit(c(0.53, 0.53), "npc"),
    arrow = grid::arrow(length = grid::unit(3.5, "mm"), type = "closed"),
    gp = grid::gpar(col = "#7E8B8F", lwd = 1.8)
  )
}

grid::grid.text(
  "Expert judgement connects the stages: document exclusions, settings, rejected alternatives, and uncertainty.",
  x = 0.5,
  y = 0.105,
  gp = grid::gpar(fontfamily = font_family, fontsize = 12.5, fontface = "bold", col = ink)
)

grDevices::dev.off()
cat(normalizePath(output_file, winslash = "/"), "\n")
