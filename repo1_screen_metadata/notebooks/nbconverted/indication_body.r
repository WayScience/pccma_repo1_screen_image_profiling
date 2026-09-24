suppressPackageStartupMessages(library(gganatogram))
source(file.path("..", "plotting_helpers", "repo1_figure_utils.R"))

output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

map_file <- file.path("..", "plotting_helpers", "indication_organ_map.csv")
max_other <- 10        # most common indications listed in the "Other" block (the rest are summed in one line)
n_body <- 50           # indications shown on the body
# organs are filled on a brown gradient (no brown in the site or plate colours of the other anatogram figure)
gradient_low <- "#E3BD8B"; gradient_high <- "#5A3418"
organ_labels <- c(urinary_bladder = "Bladder", cartilage = "Joints", skeletal_muscle = "Muscle")
# the package has no eye: the two eyes are drawn as ellipses (adult body coordinates), and the right eye is pointed to
eye_points <- data.frame(x = c(49.4, 55.8), y = c(-12.8, -12.8))
eye_radius <- c(2.1, 1.1)
other_anchor <- c(52.9, -80)   # where the "Other" block points on the body: the middle of the abdomen

compounds <- read_repo1_compounds()
organ_map <- read.csv(map_file, stringsAsFactors = FALSE)

# compounds per indication (a compound with several indications is counted in each), with the organ each maps to
indication_rows <- explode_annotation(compounds, "Indication")   # one row per compound and indication
counts <- indication_rows %>% count(value, name = "n") %>% arrange(desc(n), value) %>%
  left_join(organ_map, by = c(value = "indication"))

# the n_body most common indications that map to an organ (ties at the cut are broken alphabetically)
body_items <- counts %>% filter(!is.na(organ), organ != "none") %>% slice_head(n = n_body)
cutoff <- min(body_items$n)
# every indication at least as common as the last one shown must be in the mapping table
stopifnot(!anyNA(counts$organ[counts$n >= cutoff]))
# the most common indications with no single organ, listed in the "Other" block
other_items <- counts %>% filter(organ == "none", n >= cutoff) %>% slice_head(n = max_other)
# everything else: the remaining indications, and the compounds with at least one of them (each compound counted once)
shown <- c(body_items$value, other_items$value)
rest <- indication_rows %>% filter(!value %in% shown)
rest_line <- sprintf("All %d other indications (n = %s)", n_distinct(rest$value), format(n_distinct(rest$BROAD_CPD_ID), big.mark = ","))
# compounds per organ: those with at least one of the indications listed for it
organ_compounds <- indication_rows %>% inner_join(body_items %>% select(value, organ), by = "value") %>%
  group_by(organ) %>% summarise(n_compounds = n_distinct(BROAD_CPD_ID), .groups = "drop")
cat(rest_line, "\n")
cat(sprintf("indications: %d | on the body: %d (down to n = %d) | not placed: %d\n", nrow(counts), nrow(body_items), cutoff, nrow(other_items)))
print(as.data.frame(body_items %>% select(display_name, organ, n)))

# where each labelled organ is pointed to: the centre of its right-hand part (the image right, the body's left).
# The lung and breast shapes in the package are on the left, so those two are pointed to on the right of the chest by hand.
anchor_overrides <- list(lung = c(59.5, -36), breast = c(61.5, -47.5), eye = c(55.8, -12.8))
male <- gganatogram::hgMale_list
organ_anchor <- function(organ) {
  d <- male[[organ]]
  if (organ == "skin") {   # the skin is the outline of the body: point at its right shoulder
    shoulder <- !is.na(d$x) & -d$y > -38 & -d$y < -28
    return(data.frame(organ = organ, x = max(d$x[shoulder]), y = -d$y[shoulder][which.max(d$x[shoulder])]))
  }
  if (organ %in% names(anchor_overrides)) return(data.frame(organ = organ, x = anchor_overrides[[organ]][1], y = anchor_overrides[[organ]][2]))
  parts <- split(d[!is.na(d$x), ], cumsum(is.na(d$x))[!is.na(d$x)])   # organs made of several shapes are separated by NA rows
  centres <- bind_rows(lapply(parts, function(p) data.frame(x = mean(p$x), y = -mean(p$y))))
  cbind(organ = organ, centres[which.max(centres$x), ])
}
anchors <- bind_rows(lapply(unique(body_items$organ), organ_anchor))

outline <- male[["human_male_outline"]]
body_x <- range(outline$x, na.rm = TRUE)
body_y <- range(-outline$y, na.rm = TRUE)   # the drawing has y pointing down; the plot has it pointing up
midline <- 52.9   # x of the middle of the body
force_side <- c(breast = "right", liver = "right")   # text blocks kept on one side to balance the two columns
bilateral <- c("lung", "breast", "eye", "skin", "cartilage")   # organs that exist on both sides: pointed to on the side of their text block

# text blocks: one per organ, on the left or right of the body
header_h <- 6.6        # plot units for an organ name
line_h <- 4.15         # plot units per indication line
block_gap <- 2.8       # plot units between blocks
gap_to_text <- 9       # plot units between the body and the text columns
units_per_inch <- 21   # plot units per inch: sets the text size relative to the body
balance_weight <- 1    # how strongly blocks avoid the longer text column

blocks <- body_items %>%
  mutate(line = sprintf("%s (n = %d)", display_name, n)) %>%
  group_by(organ) %>%
  summarise(text = paste(line, collapse = "\n"), n_lines = n(), .groups = "drop") %>%
  left_join(anchors, by = "organ") %>%
  mutate(is_other = FALSE)
other_block <- data.frame(organ = "other", is_other = TRUE, x = other_anchor[1], y = other_anchor[2], n_lines = nrow(other_items) + 1,
                          text = paste(c(sprintf("%s (n = %d)", other_items$display_name, other_items$n), rest_line), collapse = "\n"))
blocks <- bind_rows(blocks, other_block) %>%
  left_join(organ_compounds, by = "organ") %>%
  mutate(header = coalesce(unname(organ_labels[organ]), str_to_sentence(str_replace_all(organ, "_", " "))),
         # each organ's header carries the number of compounds with one of its listed indications
         header = ifelse(is_other, header, sprintf("%s (n = %d)", header, n_compounds)),
         height = header_h + line_h * n_lines,
         # organs clearly on one side of the body keep their text on that side; the others can go on either
         natural = ifelse(organ %in% bilateral | is_other, "either", ifelse(x < midline - 3, "left", ifelse(x > midline + 3, "right", "either"))),
         natural = coalesce(unname(force_side[organ]), natural)) %>%
  arrange(desc(y))

# place the blocks from the top down: each goes on the side where it can sit closest to its organ (centred on it, or just below the block above)
lowest <- c(left = body_y[2] + block_gap, right = body_y[2] + block_gap)   # lowest point used so far on each side
blocks$side <- NA_character_
blocks$top <- NA_real_
for (i in seq_len(nrow(blocks))) {
  desired <- min(blocks$y[i] + blocks$height[i] / 2, body_y[2])
  sides <- if (blocks$natural[i] == "either") c("right", "left") else blocks$natural[i]
  tops <- sapply(sides, function(s) min(desired, lowest[[s]] - block_gap))
  # cost of a side: how far the block ends up from its organ, plus a penalty for a side that is already long (keeps the columns balanced)
  cost <- (desired - tops) + balance_weight * (body_y[2] - sapply(sides, function(s) lowest[[s]]))
  side <- sides[which.min(cost)]
  blocks$side[i] <- side
  blocks$top[i] <- tops[[side]]
  lowest[[side]] <- blocks$top[i] - blocks$height[i]
}
# a paired organ pointed to from the left is pointed to on its left-hand side
blocks <- blocks %>%
  mutate(x = ifelse(side == "left" & organ %in% bilateral, 2 * midline - x, x),
         header_y = top - header_h / 2, text_y = top - header_h,
         right_side = side == "right")

text_x_right <- body_x[2] + gap_to_text
text_x_left <- body_x[1] - gap_to_text
blocks <- blocks %>% mutate(text_x = ifelse(right_side, text_x_right, text_x_left))
longest <- function(side) max(nchar(strsplit(paste(blocks$text[blocks$side == side], collapse = "\n"), "\n")[[1]]))
legend_room <- 14   # plot units under the feet for the colour bar
x_lim <- c(text_x_left - longest("left") * 1.9, text_x_right + longest("right") * 1.9)
y_lim <- c(min(body_y[1] - legend_room, min(blocks$top - blocks$height)) - 2, body_y[2] + 2)
fig_width <- diff(x_lim) / units_per_inch
fig_height <- diff(y_lim) / units_per_inch
cat(sprintf("figure size: %.1f x %.1f in | blocks on the left: %d, on the right: %d\n", fig_width, fig_height, sum(!blocks$right_side), sum(blocks$right_side)))

ellipse <- function(cx, cy, rx, ry, n = 40) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + rx * cos(t), y = cy + ry * sin(t))
}
eye_value <- organ_compounds$n_compounds[organ_compounds$organ == "eye"]
eyes <- bind_rows(lapply(seq_len(nrow(eye_points)), function(i) cbind(ellipse(eye_points$x[i], eye_points$y[i], eye_radius[1], eye_radius[2]), eye = i, value = eye_value)))

# the organs to fill (the eye is drawn separately)
highlight <- organ_compounds %>% filter(organ != "eye") %>%
  transmute(organ, type = "labelled", colour = "grey", value = n_compounds) %>% as.data.frame()

fig <- gganatogram(data = highlight, fillOutline = "#f2f2f2", organism = "human", sex = "male", fill = "value") +
  {if ("eye" %in% body_items$organ) geom_polygon(data = eyes, aes(x, y, group = eye, fill = value), colour = "black", linewidth = 0.25, inherit.aes = FALSE)} +
  geom_segment(data = blocks, aes(x = x, y = y, xend = ifelse(right_side, text_x - 2, text_x + 2), yend = header_y, linetype = is_other),
               colour = "grey45", linewidth = 0.35,
               show.legend = FALSE, inherit.aes = FALSE) +
  scale_linetype_manual(values = c("solid", "22")) +
  geom_point(data = blocks %>% filter(!is_other), aes(x = x, y = y), size = 1.5, colour = "black", inherit.aes = FALSE) +
  geom_point(data = blocks %>% filter(is_other), aes(x = x, y = y), shape = 21, size = 2.6, fill = "white", colour = "black", stroke = 0.6, inherit.aes = FALSE) +
  geom_text(data = blocks, aes(x = text_x, y = header_y, label = header, hjust = ifelse(right_side, 0, 1)), fontface = "bold", size = 5.2, inherit.aes = FALSE) +
  geom_text(data = blocks, aes(x = text_x, y = text_y, label = text, hjust = ifelse(right_side, 0, 1)), vjust = 1, size = 4.1, lineheight = 1.05, inherit.aes = FALSE) +
  scale_fill_gradient(low = gradient_low, high = gradient_high, name = "Compounds per organ",
                      guide = guide_colourbar(direction = "horizontal", title.position = "top", barwidth = unit(8, "lines"), barheight = unit(0.8, "lines"),
                                              ticks.colour = "grey30", frame.colour = "grey30")) +
  coord_fixed(xlim = x_lim, ylim = y_lim, expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(6, 6, 6, 6),
        legend.position = "inside", legend.position.inside = c((mean(body_x) - x_lim[1]) / diff(x_lim), 0.005), legend.justification = c(0.5, 0),
        legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10))

print(fig)

ggsave(file.path(output_dir, "indication_body.png"), fig, width = fig_width, height = fig_height, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "indication_body.pdf"), fig, width = fig_width, height = fig_height)
