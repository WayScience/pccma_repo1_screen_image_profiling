source(file.path("..", "plotting_helpers", "repo1_figure_utils.R"))

output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

target_file <- file.path("..", "raw_data", "repo1_target_metadata_annotation.csv")
positions_file <- file.path("..", "raw_data", "repo1_target_gene_positions.csv")

bar_colour <- "#4C78A8"
oncology_colour <- "#D62728"
n_labelled <- 15   # the most targeted genes are labelled

compounds <- read_repo1_compounds()
targets_all <- read.csv(target_file, colClasses = c(entrez_id = "character"))
targets <- targets_all %>% filter(!is.na(entrez_id))

# compounds per target gene, and the position of each gene
genes <- read.csv(positions_file, colClasses = c(entrez_id = "character")) %>%
  left_join(targets %>% group_by(entrez_id) %>% summarise(n_compounds = n_distinct(compound_id)), by = "entrez_id")

# how many compounds are in the plot: all, with a target annotation, and with a target gene that has a position on the genome
n_total <- nrow(compounds)
n_with_target <- n_distinct(targets_all$compound_id)
n_placed <- targets %>% filter(entrez_id %in% genes$entrez_id) %>% pull(compound_id) %>% n_distinct()
not_placed <- setdiff(targets_all$compound_id, targets$compound_id[targets$entrez_id %in% genes$entrez_id])
not_placed_targets <- targets_all %>% filter(compound_id %in% not_placed) %>% distinct(gene_symbol) %>% pull(gene_symbol)
cat(sprintf("compounds: %d | with a target: %d | placed: %d | not placed: %d (targets: %s)\n",
            n_total, n_with_target, n_placed, length(not_placed), paste(not_placed_targets, collapse = ", ")))

# genes targeted by at least one oncology compound
oncology_compounds <- explode_annotation(compounds, "Disease Area") %>% filter(value == "oncology") %>% pull(BROAD_CPD_ID)
oncology_genes <- targets %>% filter(compound_id %in% oncology_compounds) %>% pull(entrez_id) %>% unique()
genes <- genes %>% mutate(oncology = entrez_id %in% oncology_genes)
cat(sprintf("oncology compounds: %d | with a target gene on the genome: %d | genes they target: %d\n",
            length(oncology_compounds), n_distinct(targets$compound_id[targets$compound_id %in% oncology_compounds & targets$entrez_id %in% genes$entrez_id]),
            sum(genes$oncology)))

# chromosome lengths (GRCh38)
chromosomes <- data.frame(
  chromosome = c(1:22, "X", "Y"),
  length = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
             135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
             46709983, 50818468, 156040895, 57227415))
gap <- 40e6   # empty space between chromosomes, in bases
chromosomes <- chromosomes %>%
  mutate(offset = gap / 2 + cumsum(lag(length + gap, default = 0))) %>%
  left_join(genes %>% count(chromosome, name = "n_genes"), by = "chromosome") %>%
  mutate(n_genes = coalesce(n_genes, 0L))
genome_size <- sum(chromosomes$length) + nrow(chromosomes) * gap

# angle (radians, clockwise from the top) of a genome position, and its x/y at radius r
angle_of <- function(chromosome, position) 2 * pi * (chromosomes$offset[match(chromosome, chromosomes$chromosome)] + position) / genome_size
xy <- function(theta, r) data.frame(x = r * sin(theta), y = r * cos(theta))

genes <- genes %>% mutate(theta = angle_of(chromosome, (start + end) / 2))

# ring radii
r_ring <- c(0.955, 1.0)   # chromosomes
r_base <- 0.62            # baseline of the bars
r_top <- 0.935            # bar length of the most targeted gene
max_n <- max(genes$n_compounds)
bar_length <- function(n) (r_top - r_base) * n / max_n

# chromosome ring: one arc per chromosome, with its name in the ring and its number of target genes above it
arc <- function(chromosome) {
  info <- chromosomes[chromosomes$chromosome == chromosome, ]
  theta <- seq(angle_of(chromosome, 0), angle_of(chromosome, info$length), length.out = 60)
  rbind(cbind(xy(theta, r_ring[2]), chromosome = chromosome), cbind(xy(rev(theta), r_ring[1]), chromosome = chromosome))
}
ring <- bind_rows(lapply(chromosomes$chromosome, arc)) %>%
  mutate(shade = ifelse(match(chromosome, chromosomes$chromosome) %% 2 == 1, "a", "b"))
chromosome_theta <- angle_of(chromosomes$chromosome, chromosomes$length / 2)
ring_labels <- bind_cols(chromosomes, xy(chromosome_theta, mean(r_ring)))

# bars and spokes
bars <- bind_cols(genes, xy(genes$theta, r_base), xy(genes$theta, r_base + bar_length(genes$n_compounds)) %>% rename(xend = x, yend = y))
spokes <- bind_cols(genes[c("theta", "oncology")], xy(genes$theta, r_base), data.frame(xend = 0, yend = 0))

# scale circles for the bars, labelled in the gap at the top of the circle
scale_counts <- c(25, 50, 75)
scale_circles <- bind_rows(lapply(c(0, scale_counts), function(n) {
  theta <- seq(0, 2 * pi, length.out = 300)
  cbind(xy(theta, r_base + bar_length(n)), n = n)
}))
scale_labels <- data.frame(n = scale_counts) %>% bind_cols(xy(0, r_base + bar_length(scale_counts)))

# labels for the most targeted genes, pushed apart where neighbours would overlap
top <- genes %>% slice_max(n_compounds, n = n_labelled, with_ties = FALSE) %>% arrange(theta)
min_gap <- 3.4 * pi / 180   # smallest angle between two labels
label_theta <- top$theta
for (iteration in 1:200) {
  for (i in seq_along(label_theta)[-1]) {
    shortfall <- min_gap - (label_theta[i] - label_theta[i - 1])
    if (shortfall > 0) {
      label_theta[i - 1] <- label_theta[i - 1] - shortfall / 2
      label_theta[i] <- label_theta[i] + shortfall / 2
    }
  }
}
r_label <- 1.11
degrees <- label_theta * 180 / pi
top <- top %>% mutate(
  label_theta = label_theta,
  label = sprintf("%s (n = %d)", gene_symbol, n_compounds),
  right_side = sin(label_theta) >= 0,
  text_angle = ifelse(right_side, 90 - degrees, 270 - degrees),
  hjust = ifelse(right_side, 0, 1))
top <- bind_cols(top, xy(top$label_theta, r_label) %>% rename(lx = x, ly = y), xy(top$theta, r_ring[2]) %>% rename(cx = x, cy = y),
                 xy(top$label_theta, r_label - 0.005) %>% rename(ex = x, ey = y))

# chromosome counts sit just outside the ring, nudged along their chromosome away from any gene label
half_span <- pi * chromosomes$length / genome_size
count_theta <- sapply(seq_len(nrow(chromosomes)), function(i) {
  offsets <- sort(c(0, seq(0.25, 0.85, by = 0.05) * rep(c(-1, 1), each = 13) * half_span[i]), decreasing = FALSE)
  offsets <- offsets[order(abs(offsets))]
  for (offset in offsets) {
    if (all(abs(chromosome_theta[i] + offset - label_theta) >= 3.3 * pi / 180)) return(chromosome_theta[i] + offset)
  }
  chromosome_theta[i]
})
count_labels <- bind_cols(chromosomes, xy(count_theta, 1.015)) %>%
  mutate(hjust = 0.5 - 0.5 * sin(count_theta), vjust = 0.5 - 0.5 * cos(count_theta))

spoke_labels <- c(
  oncology = sprintf("Targeted by an\noncology compound\n(n = %d genes)", sum(genes$oncology)),
  other = sprintf("Other target\ngenes\n(n = %s genes)", format(sum(!genes$oncology), big.mark = ",")))
spokes <- spokes %>% mutate(group = ifelse(oncology, "oncology", "other"))

fig <- ggplot() +
  geom_path(data = scale_circles, aes(x, y, group = n), colour = "grey85", linewidth = 0.3, linetype = "dashed") +
  geom_segment(data = spokes %>% filter(!oncology), aes(x, y, xend = xend, yend = yend, colour = group), alpha = 0.25, linewidth = 0.12) +
  geom_segment(data = spokes %>% filter(oncology), aes(x, y, xend = xend, yend = yend, colour = group), alpha = 0.55, linewidth = 0.22) +
  geom_segment(data = bars, aes(x, y, xend = xend, yend = yend), colour = bar_colour, linewidth = 0.45, lineend = "butt") +
  geom_polygon(data = ring, aes(x, y, group = chromosome, fill = shade), colour = "white", linewidth = 0.3) +
  geom_text(data = ring_labels, aes(x, y, label = chromosome), colour = "white", size = 2.9, fontface = "bold") +
  geom_label(data = scale_labels, aes(x, y, label = n), size = 2.6, colour = "grey30", fill = "white", label.size = 0,
             label.padding = unit(0.08, "lines")) +
  geom_segment(data = top, aes(cx, cy, xend = ex, yend = ey), colour = "grey40", linewidth = 0.25) +
  geom_label(data = count_labels, aes(x, y, label = sprintf("n = %d", n_genes), hjust = hjust, vjust = vjust), size = 2.7, colour = "grey25", fill = "white", label.size = 0,
             label.padding = unit(0.06, "lines")) +
  geom_text(data = top, aes(lx, ly, label = label, angle = text_angle, hjust = hjust), size = 3.1, colour = "grey10") +
  scale_fill_manual(values = c(a = "#5C6773", b = "#8E98A4"), guide = "none") +
  scale_colour_manual(values = c(oncology = oncology_colour, other = "grey55"), labels = spoke_labels, name = "Lines",
                      breaks = c("oncology", "other"),
                      guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1.1))) +
  coord_fixed(xlim = c(-1.22, 1.42), ylim = c(-1.42, 1.39), clip = "off", expand = FALSE) +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 12, 8, 8),
        legend.position = "right", legend.box.spacing = unit(2, "pt"), legend.key.height = unit(1.6, "lines"), legend.key.spacing.y = unit(8, "pt"),
        legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10),
        legend.background = element_blank(), legend.key.width = unit(2, "lines"))

print(fig)

ggsave(file.path(output_dir, "target_circos.png"), fig, width = 9.8, height = 8.9, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "target_circos.pdf"), fig, width = 9.8, height = 8.9)
