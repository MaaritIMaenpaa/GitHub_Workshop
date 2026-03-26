rm(list=ls())

# PACKAGES ----
library(tidyverse)
library(eurostat)
library(sf);library(giscoR)
library(pheatmap)
library(viridis)
library(grid)
library(patchwork)
library(ggbiplot)
library(FactoMineR)
library(factoextra)
library(ggcorrplot)
library(reshape2)
library(ggResidpanel)
library(dplyr)
library(purrr)
library(glmmTMB)
library(officer)
library(flextable)

# Theme for all figures ----
theme_set(theme_bw())

# DATA ---- 
field_data <- readRDS("processed_data/field_data.rds")
greenhouse_data <- readRDS("processed_data/greenhouse_data.rds")

# Rescale c and n
greenhouse_data <- greenhouse_data |>
  mutate(c = 10000*c_g_100g_ts,
         n = 10000*n_g_100g_ts) 

## Field ionome ----
field_ionome <- field_data |>
  select(c(id, al:zn, straw_si:straw_n)) |>
  rename("Al" = "al",
         "Ca" = "ca",
         "Cu" = "cu",
         "Fe" = "fe",
         "K" = "k",
         "Mg" = "mg",
         "Mn" = "mn",
         "Na" = "na",
         "P" = "p",
         "S" = "s",
         "Zn" = "zn",
         "Si" ="straw_si",
         "C" = "straw_c",
         "N" = "straw_n") 
field_ionome

## Greenhouse ionome ----
greenhouse_ionome <- greenhouse_data |>
  select(c(pot_id, 
           silicon_mg_kg, 
           c_g_100g_ts:n_g_100g_ts, 
           al:zn)) |>
  rename("Si" = "silicon_mg_kg",
         "C" = "c_g_100g_ts",
         "N" = "n_g_100g_ts",
         "Al" = "al",
         "Ca" = "ca",
         "Cu" = "cu",
         "Fe" = "fe",
         "K" = "k",
         "Mg" = "mg",
         "Mn" = "mn",
         "Na" = "na",
         "P" = "p",
         "S" = "s",
         "Zn" = "zn")
greenhouse_ionome

#
#
#

# FIELD ----

#
## Map ----
field_data <- field_data |>
  mutate(latitude = as.numeric(latitude),
         longitude = as.numeric(longitude))

field_data_sf <- field_data |> 
  st_as_sf(coords = c("longitude", "latitude"), 
           remove = FALSE, crs = 4326)

eu_sh <- 
  get_eurostat_geospatial(resolution = 10, nuts_level = 0)

exp1_map <- 
  ggplot(eu_sh) +
  geom_sf() +
  geom_sf(data=field_data_sf) +
  geom_sf_label(data = field_data_sf,
                aes(label = id),   
                nudge_y = 0.1,     
                alpha = 0.8,
                size = 3) +
  scale_x_continuous(limits = c(7.5, 13)) +
  scale_y_continuous(limits = c(54.5, 58)) +
  labs(x="", y="")
exp1_map

#
## Heatmap ----
field_ionome_hm <- field_ionome |>
  tibble::column_to_rownames("id")

exp1_heatmap <-
  pheatmap(field_ionome_hm,
           scale = "column",
           clustering_method = "ward.D2",
           color = viridis(10),
           angle_col = 0,
           fontsize_row = 12,
           fontsize_col = 10)
exp1_heatmap

### Combined map and heatmap ----

# Extract heatmap object
heat_grob <- exp1_heatmap$gtable

# Combine figures
exp1_combined <- exp1_map + wrap_elements(heat_grob) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(0.75, 1))
exp1_combined

# Save figure
ggsave(plot=exp1_combined,
       filename="figures/ionomes/exp1_field_ionome.png",
       width=20, height=10, units="cm")

#
#
#

## Transformations ----

# Log transformation to stabilise variance
element_cols <- setdiff(names(field_ionome), "id")
field_ionome_log <- field_ionome |>
  mutate(across(all_of(element_cols), ~ log(.x)))

## PCA ----
pca_field <- prcomp(field_ionome_log[, element_cols],
                    center = TRUE,
                    scale. = TRUE)

summary(pca_field)
pca_field

# Biplot
field_biplot <-
  ggbiplot(pca_field,
         obs.scale = 1,
         var.scale = 1, 
         labels = field_ionome_log$id, 
         ellipse = FALSE,
         circle = FALSE, 
         labels.size = 4, 
         varname.adjust = 2, 
         varname.color = "steelblue") +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) 
field_biplot

# Save the plot
ggsave(plot=field_biplot,
       filename="figures/ionomes/exp1_biplot.png",
       width=15, height=15, units="cm")

## Correlations ----
cor_mat <- cor(field_ionome_log[element_cols], 
               use = "pairwise.complete.obs")

field_correlogram <-
  ggcorrplot(cor_mat, 
           method = "square",   
           type = "upper",
           lab = FALSE,
           lab_size = 3, 
           tl.cex = 10, 
           sig.level = 0.05, 
           insig = "blank", 
           hc.order = TRUE,
           colors = c("blue", "white", "red")) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust = 0.5))
field_correlogram

# Save plot
ggsave(plot=field_correlogram,
       filename="figures/ionomes/exp1_correlogram.png",
       width=15, height=15, units="cm")

#
#
#

## Barchart ----

# Long-form data
field_long <- field_ionome |>
  mutate(N = N*10000,
         C = C*10000) |>
  select(-C) |>
  pivot_longer(cols=Al:N,
               names_to="element",
               values_to = "content")

# Select colours (adding to colour blind palette in ggthemes)
my_13_colors <- c("#1b9e77", "#d95f02", "#7570b3",
                  "#e7298a", "#66a61e", "#e6ab02",
                  "#a6761d", "#666666", "#1f78b4",
                  "#33a02c", "#fb9a99", "#e31a1c",
                  "#6a3d9a")

# Figure
ggplot(field_long, aes(x=id, fill=element, y=content)) +
  theme_classic() +
  geom_bar(position=position_dodge(),  stat="identity") +
  scale_fill_manual(values=my_13_colors) +
  scale_y_continuous(limits=c(0,30000), expand=c(0,0)) +
  labs(x="Site", 
       y = "Element concentration (\u00B5g)",
       fill="Element") 

# Save figure
ggsave(plot=last_plot(),
       filename="figures/ionomes/exp1_bars.png",
       units="cm", width=30, height=10)

#
#
#

# GREENHOUSE ----

### Heatmap ----
greenhouse_ionome_hm <- greenhouse_ionome |>
  tibble::column_to_rownames("pot_id")
greenhouse_ionome_hm

# Row annotations: Species and Treatment by pot
row_anno <- greenhouse_data |>
  select(pot_id, treatment, species) |>
  distinct(pot_id, .keep_all = TRUE) |>
  column_to_rownames("pot_id") |>
  rename("Species" = "species",
         "Treatment" = "treatment")

# Example palettes: tweak to your levels
species_lvls   <- sort(unique(row_anno$Species))
treatment_lvls <- sort(unique(row_anno$Treatment))

species_cols <- setNames(c("orange", "maroon", "steelblue1", "darkgreen"), species_lvls)
treat_cols   <- setNames(c("grey", "tomato", "purple", "steelblue"), treatment_lvls)

ann_colors <- list(Species = species_cols, 
                   Treatment = treat_cols)

exp2_heatmap <- 
  pheatmap(greenhouse_ionome_hm,
           scale = "column",
           clustering_method = "ward.D2",
           color = viridis(10),
           angle_col = 90,
           fontsize_row = 10,
           fontsize_col = 10,
           annotation_row = row_anno,         
           annotation_colors = ann_colors, show_rownames = FALSE,
           #border_color = NA, 
           cluster_rows = TRUE)
exp2_heatmap

# Save figure
ggsave(plot=exp2_heatmap,
       bg = "white",
       filename="figures/ionomes/exp2_greenhouse_ionome.png",
       width=12.5, height=15, units="cm")

#
#
#

## Transformations ----

# Log transformation to stabilise variance
element_cols <- setdiff(names(greenhouse_ionome), "pot_id")

# Turn all negative values to 0
greenhouse_ionome_log <- greenhouse_ionome |>
  mutate(across(where(is.numeric), ~ ifelse(.x < 0, 0, .x))) |>
  mutate(across(all_of(element_cols), ~ log(.x + 0.001)))

## PCA ----
pca_greenhouse <- prcomp(greenhouse_ionome_log[, element_cols],
                         center = TRUE,
                         scale. = TRUE)

summary(pca_greenhouse)
pca_greenhouse

# Biplot
greenhouse_biplot <-
  ggbiplot(pca_greenhouse,
         obs.scale = 1,
         var.scale = 1, 
         groups = greenhouse_data$species,
         ellipse = TRUE,
         circle = FALSE, 
         labels.size = 4, 
         varname.adjust = 2, 
         varname.color = "steelblue") +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) 
greenhouse_biplot

# Save the plot
ggsave(plot=greenhouse_biplot,
       filename="figures/ionomes/exp2_biplot.png",
       width=15, height=15, units="cm")

ggbiplot(pca_greenhouse,
         obs.scale = 1,
         var.scale = 1, 
         groups = greenhouse_data$treatment,
         ellipse = TRUE,
         circle = FALSE, 
         labels.size = 4, 
         varname.adjust = 2, 
         varname.color = "steelblue") +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) 

greenhouse_data$treatspec <- paste0(greenhouse_data$species,
                                    "_",
                                    greenhouse_data$treatment)

ggbiplot(pca_greenhouse,
         obs.scale = 1,
         var.scale = 1, 
         groups = greenhouse_data$treatspec,
         ellipse = TRUE,
         circle = FALSE, 
         labels.size = 4, 
         varname.adjust = 2, 
         varname.color = "steelblue") +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) 

#
## Correlations ----
cor_mat <- cor(greenhouse_ionome_log[element_cols], 
               use = "pairwise.complete.obs")

greenhouse_correlogram <- 
  ggcorrplot(cor_mat, 
           method = "square",   
           type = "upper",
           lab = FALSE,
           lab_size = 3, 
           tl.cex = 10, 
           sig.level = 0.05, 
           insig = "blank", 
           hc.order = TRUE,
           colors = c("blue", "white", "red")) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0.5, 
                                   hjust = 0.5))
greenhouse_correlogram

# Save plot
ggsave(plot=greenhouse_correlogram,
       filename="figures/ionomes/exp2_correlogram.png",
       width=15, height=15, units="cm")

#
#
#

## Barchart ----

# Long-form data
greenhouse_long <- greenhouse_data |>
  select(c(pot_id, 
           species, treatment,
           silicon_mg_kg, 
           c_g_100g_ts:n_g_100g_ts, 
           al:zn)) |>
  rename("Si" = "silicon_mg_kg",
         "C" = "c_g_100g_ts",
         "N" = "n_g_100g_ts",
         "Al" = "al",
         "Ca" = "ca",
         "Cu" = "cu",
         "Fe" = "fe",
         "K" = "k",
         "Mg" = "mg",
         "Mn" = "mn",
         "Na" = "na",
         "P" = "p",
         "S" = "s",
         "Zn" = "zn") |>
  mutate(C = C*10000,
         N = N*10000) |>
  select(-c(C)) |>
  pivot_longer(cols=Si:Zn,
               names_to="element",
               values_to = "content") |>
  group_by(species, treatment, element) |>
  summarise(content = mean(content))
greenhouse_long

# Select colours (adding to colour blind palette in ggthemes)
my_13_colors <- c("#1b9e77", "#d95f02", "#7570b3",
                  "#e7298a", "#66a61e", "#e6ab02",
                  "#a6761d", "#666666", "#1f78b4",
                  "#33a02c", "#fb9a99", "#e31a1c",
                  "black",
                  "#6a3d9a")

# Figure
ggplot(greenhouse_long, aes(x=treatment, fill=element, y=content)) +
  theme_classic() +
  facet_grid(~species) +
  geom_bar(position=position_dodge(),  stat="identity") +
  scale_fill_manual(values=my_13_colors) +
  scale_y_continuous(limits=c(0,40000), expand=c(0,0)) +
  labs(x="Treatment", 
       y = "Element concentration (mg/kg)",
       fill="Element") 

# Save figure
ggsave(plot=last_plot(),
       filename="figures/ionomes/exp2_bars.png",
       units="cm", width=30, height=10)

#
## Association models ----

element_list <- 
  greenhouse_data |>
  select(c(al:zn), c, n) |> # select all elements
  select(-c(c, na)) |> # Manually exclude the ones where standard model doesn't fit
  names()

# Object to store ANOVA results
anova_list <- list()

for(i in seq_along(element_list)) {
  
  # Fit a model
  mod <- glmmTMB(formula = as.formula(paste0(element_list[i], " ~ species + treatment + species:treatment + silicon_mg_kg")),
    dispformula = ~ treatment,
    data = greenhouse_data)
  
  # Save ANOVA table
  anova_list[[ element_list[i] ]] <- car::Anova(mod)
  
  # Create and save residual diagnostic plot 
  residuals <- resid_auxpanel(resid(mod, type = "pearson"), fitted(mod))
  ggsave(plot=residuals, 
         filename = file.path("output/exp2/ionomes/model_fit/", paste0("resid_plot_", element_list[i], ".png")),
         width = 15, height = 12, units="cm")
  
  # Create an effect size plot
  # Plot the effect
  pred_plot <- ggpredict(mod, terms = "silicon_mg_kg") |> plot()
  
  ggsave(plot = pred_plot, 
         filename = paste0("figures/ionomes/individual_elements/", element_list[i], "_effect.png"),
         width=15, height=10, units="cm")
  
}

anova_df <- map2_dfr(anova_list,
                     names(anova_list),
                     ~ as.data.frame(.x) |> 
                       mutate(element = .y, term = rownames(.x), .before = 1))
anova_df

#### C ----

# Figures
ggplot(greenhouse_data, aes(x=treatment, y=c)) +
  facet_wrap(~species) +
  geom_point()

ggplot(greenhouse_data, aes(x=silicon_mg_kg, y=c)) +
  facet_grid(treatment~.) +
  geom_point()

# Fit the model
mod_c <- glmmTMB(c ~ species + treatment + species:treatment + silicon_mg_kg,
                 dispformula = ~ treatment,
                 data = greenhouse_data)

# Assumption check
resid_auxpanel(resid(mod_c, type = "pearson"), fitted(mod_c))

# Rescale C
mod_c <- glmmTMB(scale(c) ~ species + treatment + species:treatment + silicon_mg_kg,
                 dispformula = ~ treatment,
                 data = greenhouse_data)

# Residual diagnostics
(residuals <- resid_auxpanel(resid(mod_c, type = "pearson"), fitted(mod_c)))
ggsave(plot=residuals, 
       filename = file.path("output/exp2/ionomes/model_fit/", paste0("resid_plot_c.png")),
       width = 15, height = 12, units="cm")

# Results
C_result <- car::Anova(mod_c) |> 
  as.data.frame() |>
  tibble::rownames_to_column(var = "term") |>
  mutate(element=rep("c", 4), .before=1)

# Add to the full table
anova_df <- rbind(anova_df, C_result)

# Plot the effect
ggpredict(mod_c, terms = "silicon_mg_kg") |> plot()

ggsave(plot = last_plot(), 
       filename = "figures/ionomes/individual_elements/c_effect.png",
       width=15, height=10, units="cm")

#
#### Na ----

# Fit a model
mod_na <- glmmTMB(na ~ species + treatment + species:treatment + silicon_mg_kg,
                 dispformula = ~ treatment,
                 data = greenhouse_data)

# Residual diagnostics 
resid_auxpanel(resid(mod_na, type = "pearson"), fitted(mod_na))

# Change the variance function
mod_na <- glmmTMB(na ~ species + treatment + species:treatment + silicon_mg_kg,
                  dispformula = ~ species,
                  data = greenhouse_data)

# Residual diagnostics
(residuals <- resid_auxpanel(resid(mod_na, type = "pearson"), fitted(mod_na)))
ggsave(plot=residuals, 
       filename = file.path("output/exp2/ionomes/model_fit/", paste0("resid_plot_c.png")),
       width = 15, height = 12, units="cm")

# Results
Na_result <- car::Anova(mod_na) |> 
  as.data.frame() |>
  tibble::rownames_to_column(var = "term") |>
  mutate(element=rep("na", 4), .before=1)

# Add to the full table
anova_df <- rbind(anova_df, Na_result)

# Plot the effect
ggpredict(mod_na, terms = "silicon_mg_kg") |> plot()

ggsave(plot = last_plot(), 
       filename = "figures/ionomes/individual_elements/na_effect.png",
       width=15, height=10, units="cm")

#
## Table as output ----

# Make the table nicer
(anova_pretty <- anova_df |>
    mutate(across(where(is.numeric), ~round(.x, 3))) |>
    mutate(Chisq = round(Chisq, 2)) |>
    dplyr::rename(Element = element,
                  Term = term,
                  ChiSq = Chisq,
                  df = Df,
                  P = `Pr(>Chisq)`) |>
    mutate(P = ifelse(P == 0, "<0.001", P),
           Term = factor(Term)) |>
    mutate(Term = fct_recode(Term, "Species" = "species",
                             "Treatment" = "treatment",
                             "Silicon (mg/kg)" = "silicon_mg_kg",
                             "Species:Treatment" = "species:treatment"),
           Element = fct_recode(Element, 
                                "C" = "c",
                                "N" = "n",
                                "Al" = "al",
                                "Ca" = "ca",
                                "Cu" = "cu",
                                "Fe" = "fe",
                                "K" = "k",
                                "Mg" = "mg",
                                "Mn" = "mn",
                                "Na" = "na",
                                "P" = "p",
                                "S" = "s",
                                "Zn" = "zn")) |>
    arrange(Element)) 

names(anova_pretty)[names(anova_pretty) == "ChiSq"] <- "\u03C7\u00B2"

# Build a flextable
(ft <- flextable(anova_pretty) |>
    autofit() |>
    theme_booktabs() |>
    align(j = c("Element", "Term"), align = "left", part = "all") |>
    align(j = c("χ²", "df", "P"), align = "right", part = "all") |>
    valign(valign = "center", part = "all") |>
    bold(part = "header") |>
    fontsize(size = 10, part = "all") |>
    padding(padding = 3, part = "all") |>
    bg(i = c(5:8, 13:16, 21:24, 
             29:32, 37:40, 45:48),
       bg = "grey90", part = "body"))
  

# Create a Word document and add the table
doc <- read_docx() |>
  body_add_par("ANOVA Results", style = "heading 1") |>
  body_add_flextable(ft)

print(doc, target = "output/exp2/ionomes/exp2_anova_results.docx")

#
#
#

# END ----