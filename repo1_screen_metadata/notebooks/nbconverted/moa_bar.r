source(file.path("..", "plotting_helpers", "repo1_figure_utils.R"))

output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

moa_file <- file.path("..", "raw_data", "repo1_moa_metadata_annotation.csv")

n_top <- 40             # MOAs shown individually
group_gap <- 0.7        # extra space between the groups of bars, in bar heights

compounds <- read_repo1_compounds()
n_compounds <- nrow(compounds)
moa <- read.csv(moa_file)   # one row per compound and MOA

# for each compound, the number of other MOAs it is annotated with
per_compound <- moa %>% group_by(compound_id) %>% mutate(other_moas = n_distinct(moa) - 1) %>% ungroup()

moa_counts <- moa %>% count(moa, name = "total") %>% arrange(desc(total))
# the n_top most common MOAs (ties at the cut are broken alphabetically)
shown <- moa_counts %>% arrange(desc(total), moa) %>% slice_head(n = n_top)
n_tied_out <- sum(moa_counts$total == min(shown$total)) - sum(shown$total == min(shown$total))
n_moas_other <- nrow(moa_counts) - nrow(shown)

no_moa_label <- "no MOA"
other_label <- sprintf("other (%d MOAs)", n_moas_other)

# the bars: each shown MOA, "other" (compounds with any of the remaining MOAs, counted once) and compounds with no MOA
bars <- bind_rows(
  per_compound %>% filter(moa %in% shown$moa) %>% transmute(compound_id, other_moas, bar = moa),
  per_compound %>% filter(!moa %in% shown$moa) %>% distinct(compound_id, other_moas) %>% mutate(bar = other_label),
  compounds %>% filter(!BROAD_CPD_ID %in% moa$compound_id) %>% transmute(compound_id = BROAD_CPD_ID, other_moas = NA_real_, bar = no_moa_label)
)
cat(sprintf("compounds: %d | with a MOA: %d | MOAs: %d (shown: %d, other: %d; %d more MOAs tie with the last one shown) | most other MOAs on one compound: %d\n",
            n_compounds, n_distinct(moa$compound_id), nrow(moa_counts), nrow(shown), n_moas_other, n_tied_out, max(per_compound$other_moas)))

# group the MOAs by drug action (the last word of the name), largest group first
actions <- c(inhibitor = "Inhibitors", antagonist = "Antagonists", agonist = "Agonists", blocker = "Blockers",
             modulator = "Modulators", activator = "Activators")
action_pattern <- sprintf("(%s)$", paste(names(actions), collapse = "|"))
action_of <- function(x) replace_na(str_extract(x, action_pattern), "other")
special <- c(no_moa_label, other_label)

layout <- bars %>% count(bar, name = "total") %>%
  mutate(action = ifelse(bar %in% special, "special", action_of(bar)))
group_order <- layout %>% filter(action != "special") %>% group_by(action) %>% summarise(size = sum(total)) %>%
  arrange(desc(size)) %>% pull(action)
# groups in which every bar is a receptor MOA (e.g. "adrenergic receptor antagonist") are named "Receptor antagonists"
receptor_groups <- layout %>% filter(!action %in% c("special", "other")) %>% group_by(action) %>%
  summarise(all_receptor = all(str_detect(bar, paste0("receptor\\s+", action, "$")))) %>% filter(all_receptor) %>% pull(action)

layout <- layout %>%
  mutate(action = factor(action, c("special", group_order))) %>%
  arrange(action, match(bar, special), desc(total)) %>%
  # y position from the top down, with extra space where the group changes
  mutate(step = 1 + group_gap * (action != lag(action, default = first(action))),
         y = -cumsum(step),
         # the action word is written once per group, so it is dropped from the bar labels
         stripped = str_remove(bar, paste0("\\s*", action_pattern)),
         # keep the full name where the rest would be too short to make sense (e.g. "DNA inhibitor")
         stripped = ifelse(action %in% receptor_groups, str_remove(stripped, "\\s*receptor$"), stripped),
         label = ifelse(action %in% c("special", "other") | (nchar(stripped) <= 4 & !action %in% receptor_groups), bar, stripped),
         label = paste0(toupper(substring(label, 1, 1)), substring(label, 2)),
         panel = factor(ifelse(action == "special", "special", "moa"), c("special", "moa")))

# a shaded band and a name for each group of bars
group_bands <- layout %>% filter(action != "special") %>% group_by(action) %>%
  summarise(ymin = min(y) - 0.6, ymax = max(y) + 0.6, .groups = "drop") %>%
  mutate(shade = rep(c("grey93", "white"), length.out = n()),
         name = ifelse(action == "other", "Other actions", actions[as.character(action)]),
         name = ifelse(action %in% receptor_groups, paste0("Receptor\n", str_to_lower(name)), name),
         y = (ymin + ymax) / 2, panel = factor("moa", c("special", "moa")))

plot_data <- bars %>% count(bar, other_moas) %>% left_join(layout %>% select(bar, y, panel), by = "bar") %>%
  mutate(stack_order = coalesce(other_moas, -1))
max_other <- max(per_compound$other_moas)

fig <- ggplot(plot_data, aes(n, y, fill = other_moas, group = stack_order)) +
  geom_rect(data = group_bands, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = group_bands$shade, inherit.aes = FALSE) +
  geom_col(orientation = "y", width = 0.75, position = position_stack(reverse = TRUE)) +
  geom_text(data = layout, aes(total, y, label = scales::comma(total)), hjust = -0.25, size = 3.6, colour = "grey15", inherit.aes = FALSE) +
  geom_text(data = group_bands, aes(Inf, y, label = name), hjust = 1.08, fontface = "bold", size = 4.4, colour = "grey40", lineheight = 0.9, inherit.aes = FALSE) +
  facet_wrap(~panel, ncol = 1, scales = "free") +   # each panel gets its own axis
  scale_fill_gradient(low = "#9ECAE1", high = "#08306B", na.value = "#BAB0AC", limits = c(0, max_other), breaks = 0:max_other,
                      name = "Other MOAs\nper compound",
                      guide = guide_colourbar(barwidth = unit(1, "lines"), barheight = unit(9, "lines"),
                                          ticks.colour = "grey30", frame.colour = "grey30")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), labels = scales::comma) +
  scale_y_continuous(breaks = layout$y, labels = layout$label, expand = expansion(add = 0.8)) +
  labs(x = "Number of compounds", y = NULL, title = "Mechanisms of action of the REPO1 compounds",
       subtitle = sprintf("(n = %s total compounds)", format(n_compounds, big.mark = ","))) +
  theme_repo1() +
  theme(strip.text = element_blank(), strip.background = element_blank(), panel.spacing.y = unit(0.6, "lines"),
        legend.position = "right", legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(6, 16, 6, 6))

# facet_wrap gives the panels equal heights: make each proportional to its number of bars (plus the axis padding)
fig_grob <- ggplotGrob(fig)
panel_rows <- fig_grob$layout$t[grepl("^panel", fig_grob$layout$name)]
panel_span <- layout %>% group_by(panel) %>% summarise(span = max(y) - min(y) + 1 + 1.6) %>% arrange(panel) %>% pull(span)
fig_grob$heights[panel_rows] <- grid::unit(panel_span, "null")

grid::grid.newpage()
grid::grid.draw(fig_grob)

ggsave(file.path(output_dir, "moa_bar.png"), fig_grob, width = 10.2, height = 12, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "moa_bar.pdf"), fig_grob, width = 10.2, height = 12)
