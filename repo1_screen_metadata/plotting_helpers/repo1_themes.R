# Colours and the ggplot theme shared by every REPO1 compound figure.
suppressPackageStartupMessages(library(ggplot2))

# clinical phase, in the order shown in legends and stacks (most advanced first)
phase_levels <- c("Launched", "Withdrawn", "Phase 1-3", "Preclinical", "Not annotated")
phase_colors <- c(
  "Launched"      = "#4C78A8",
  "Withdrawn"     = "#E45756",
  "Phase 1-3"     = "#F2A93B",
  "Preclinical"   = "#72B7B2",
  "Not annotated" = "#BAB0AC"
)

theme_repo1 <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "grey15"),
      axis.title = element_text(colour = "grey10"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, colour = "grey25"),
      plot.title.position = "plot"
    )
}
