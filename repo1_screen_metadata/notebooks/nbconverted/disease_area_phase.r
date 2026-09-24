suppressPackageStartupMessages(library(forcats))

source(file.path("..", "plotting_helpers", "repo1_figure_utils.R"))

output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# disease areas with fewer compounds than this are grouped into "other"
min_compounds <- 10

compounds <- read_repo1_compounds()
n_compounds <- nrow(compounds)

# one row per compound and disease area (a compound with several disease areas is counted in each)
areas <- explode_annotation(compounds, "Disease Area")

area_totals <- areas %>% count(value, name = "total")
shown_areas <- area_totals %>% filter(total >= min_compounds) %>% pull(value)
n_areas_other <- nrow(area_totals) - length(shown_areas)

no_area_label <- "no disease area"
other_label <- sprintf("other (%d areas)", n_areas_other)

# the bars: each shown area, "other" (compounds in any of the grouped areas, counted once) and compounds with no disease area
bars <- bind_rows(
  areas %>% filter(value %in% shown_areas) %>% transmute(BROAD_CPD_ID, phase, bar = value),
  areas %>% filter(!value %in% shown_areas) %>% distinct(BROAD_CPD_ID, phase) %>% mutate(bar = other_label),
  compounds %>% filter(is.na(`Disease Area`)) %>% transmute(BROAD_CPD_ID, phase, bar = no_area_label)
)

bar_totals <- bars %>% count(bar, name = "total")
# axis order: the two special bars at the top, the disease areas below them from most to least compounds
area_order <- bar_totals %>% filter(!bar %in% c(no_area_label, other_label)) %>% arrange(total) %>% pull(bar)
bar_levels <- c(area_order, other_label, no_area_label)

plot_data <- bars %>% count(bar, phase) %>% mutate(bar = factor(bar, bar_levels))
bar_totals <- bar_totals %>% mutate(bar = factor(bar, bar_levels))

cat(sprintf("compounds: %d | no disease area: %d | disease areas: %d (shown: %d, other: %d)\n",
            n_compounds, bar_totals$total[bar_totals$bar == no_area_label], nrow(area_totals), length(shown_areas), n_areas_other))

# legend labels with the number of compounds in each phase (over all compounds, each counted once)
phase_counts <- table(compounds$phase)
phase_labels <- setNames(sprintf("%s (n = %s)", phase_levels, format(as.integer(phase_counts[phase_levels]), big.mark = ",", trim = TRUE)), phase_levels)

fig <- ggplot(plot_data, aes(n, bar, fill = phase)) +
  geom_col(width = 0.75, position = position_stack(reverse = TRUE), colour = "grey30", linewidth = 0.15) +
  geom_text(data = bar_totals, aes(total, bar, label = total), hjust = -0.25, size = 3.6, colour = "grey15", inherit.aes = FALSE) +
  scale_fill_manual(values = phase_colors, breaks = phase_levels, labels = phase_labels, drop = FALSE, name = "Phase") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_y_discrete(labels = function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))) +
  labs(x = "Compounds", y = NULL, title = "Disease areas of the REPO1 compounds",
       subtitle = sprintf("(n = %s total compounds)", format(n_compounds, big.mark = ","))) +
  theme_repo1() +
  theme(legend.position = "inside", legend.position.inside = c(0.72, 0.5), legend.justification = c(0.5, 0.5),
        legend.background = element_rect(colour = "grey60", linewidth = 0.3))

print(fig)

ggsave(file.path(output_dir, "disease_area_phase.png"), fig, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "disease_area_phase.pdf"), fig, width = 8, height = 7)
