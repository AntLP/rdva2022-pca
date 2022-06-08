# Setup ----

library(data.table)
library(haven)
library(tidyverse)
library(tidymodels)
library(tidytext)
library(ggspatial)
library(sf)
library(ggthemr)
library(ggridges)
library(ggalt)

PROMUTUEL_YELLOW = "#FDDB00"
PROMUTUEL_GREY = "#53565A"


# Import data ----

dir.create("temp")
unzip("./data/ad_boundaries.zip", exdir = "temp")
unzip("./data/sdr_boundaries.zip", exdir = "temp")

ad_layer <- read_sf("temp/lda_000b16a_e.shp") %>% 
  janitor::clean_names() %>% 
  filter(pruid == "24") %>% 
  rename(adidu = dauid,
         sdridu = csduid,
         shape = geometry)


sdr_layer <- read_sf("temp/lcsd000b16a_e.shp") %>% 
  janitor::clean_names() %>% 
  filter(pruid == "24") %>% 
  rename(name = csdname)


unlink("temp", recursive=TRUE)


geo_demo_dat <- read_csv(file = "./data/Donnees_GeoDemo_QC_NB_2016.csv",
                         col_types = cols(ADIDU = col_character(), 
                                          PRIDU = col_character(),
                                          SDRIDU = col_character()), 
                         locale = locale(encoding = "latin1")) %>% 
  janitor::clean_names() %>% 
  filter(pridu == "24")


geo_demo_dat <- bind_rows(
  geo_demo_dat,
  ad_layer %>% 
    select(adidu, sdridu) %>% 
    st_drop_geometry() %>% 
    anti_join(geo_demo_dat, by = "adidu") %>% 
    filter(sdridu %in% c("2423027", "2466023"))
)

# Data manipulation ----

geo_demo_dat <- geo_demo_dat %>% 
  mutate(
    prop_immigrants = nb_immigrants / (nb_immigrants + nb_non_immigrants),
    prop_bacc = nb_bacc / population_2016
  ) %>% 
  mutate(
    across(
      .cols = c(densite_pop_km2, taille_moy_menage_priv, revenu_total_median_2015, revenu_total_median_menage, frais_logement_mensuel_moy, taux_chomage),
      .fns = ~ifelse(.x == 0, NA, .x)
    ),
  ) %>% 
  left_join(ad_layer %>% filter(pruid == "24") %>% select(adidu, shape), by = "adidu") %>% 
  mutate(cent_x = st_coordinates(st_centroid(shape))[,1],
         cent_y = st_coordinates(st_centroid(shape))[,2])

rec <- geo_demo_dat %>% 
  recipe(~ adidu + densite_pop_km2 + revenu_total_median_2015 + 
           revenu_total_median_menage + prop_immigrants + frais_logement_mensuel_moy + 
           valeur_moy_resi + age_moy_residences + prop_bacc) %>% 
  update_role("adidu", new_role = "key") %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_impute_mean(all_numeric_predictors()) %>% 
  step_pca(all_numeric(), 
           keep_original_cols = T,
           num_comp = 999)


prep_rec <- prep(rec)
baked_rec <- prep_rec %>% bake(new_data = NULL)


pca_res <- tidy(prep_rec, 3)



# Presentation data ----

plt_geodat_range <- geo_demo_dat %>% 
  na.omit() %>% 
  summarise(
    across(
      .cols = c(densite_pop_km2, revenu_total_median_2015, revenu_total_median_menage,
                prop_immigrants, frais_logement_mensuel_moy, valeur_moy_resi, 
                age_moy_residences, prop_bacc),
      .fns = list(min, max),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(
    fun = ifelse(grepl("_1", name), "min", "max"),
    var = gsub('.{2}$', '', name)
  ) %>% 
  pivot_wider(
    id_cols = var,
    names_from = fun,
    values_from = value
  ) %>% 
  arrange(var) %>% 
  ggplot() + 
  geom_segment(aes(x = fct_rev(as.factor(var)), xend = as.factor(var), y = min, yend = max), 
               color = PROMUTUEL_YELLOW,
               size = 4) +
  coord_flip() + 
  scale_y_continuous(labels = scales::comma,
                     trans = "log10") +
  labs(
    x = "",
    y = "Valeurs"
  )

plt_geodat_boxplot <- geo_demo_dat %>% 
  pivot_longer(cols = c(densite_pop_km2, revenu_total_median_2015, revenu_total_median_menage,
                        prop_immigrants, frais_logement_mensuel_moy, valeur_moy_resi, 
                        age_moy_residences, prop_bacc), 
               names_to = "var", 
               values_to = "val") %>% 
  ggplot(aes(x = fct_rev(as.factor(var)), y = val)) + 
  geom_boxplot(color = "black", 
               fill = PROMUTUEL_YELLOW, 
               size = 0.5, 
               outlier.alpha = 0,
               na.rm = T) + 
  coord_flip() + 
  scale_y_continuous(labels = scales::comma,
                     trans = "log10") +
  labs(
    x = "",
    y = "Valeurs"
  )

# PCA results ----

## pourcent variance ----

pct_var <- tibble(
  PC = paste0(rep("PC", length(prep_rec$steps[[3]]$res$sdev)), 1:length(prep_rec$steps[[3]]$res$sdev)),
  pct_var = prep_rec$steps[[3]]$res$sdev^2 / sum(prep_rec$steps[[3]]$res$sdev^2)
) %>% 
  ggplot(aes(x = PC, y = pct_var)) +
  geom_col(alpha = 0.8) +
  # coord_flip() +
  labs(x = "", y = "", title = "% de la variance")

pct_var_cum <- tibble(
  PC = c("", paste0(rep("PC", length(prep_rec$steps[[3]]$res$sdev)), 1:length(prep_rec$steps[[3]]$res$sdev))),
  pct_var = c(0, cumsum(prep_rec$steps[[3]]$res$sdev^2 / sum(prep_rec$steps[[3]]$res$sdev^2)))
) %>% 
  ggplot(aes(x = PC, y = pct_var, group = 1)) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(x = "", y = "", title = "% cumulatif de la variance")


plt_pct_var <- cowplot::plot_grid(pct_var, pct_var_cum)

## PC 1 et 2 scores ----

plt_pc_scores <- pca_res %>%
  filter(component %in% paste0("PC", 1:2)) %>%
  group_by(component) %>%
  top_n(5, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component),
         value_flip = case_when(
           component %in% c("PC2") ~ -value,
           TRUE ~ value
         )) %>%
  ggplot(aes(value_flip, terms, fill = value_flip > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  scale_fill_manual(palette = function(x){c("#d4716a", "#8bd46a")}) +
  labs(
    x = "Contribution",
    y = NULL, fill = "Positif?"
  ) + 
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = 15))






limit <- 0.05

ad_layer_pc <- ad_layer %>% 
  inner_join(baked_rec, by = "adidu") %>% 
  mutate(
    PC2 = -PC2,
    pc1_rescaled = PC1 %>% pmin(quantile(PC1, c(1-limit))) %>% pmax(quantile(PC1, c(limit))),
    pc2_rescaled = PC2 %>% pmin(quantile(PC2, c(1-limit))) %>% pmax(quantile(PC2, c(limit))),
    score_aisance_financiere = (pc1_rescaled - min(pc1_rescaled)) / (max(pc1_rescaled) - min(pc1_rescaled)) * 100,
    score_urbanisme = (pc2_rescaled - min(pc2_rescaled)) / (max(pc2_rescaled) - min(pc2_rescaled)) * 100,
    mixed_id = ifelse(is.na(ctuid), adidu, ctuid)
  ) %>% 
  group_by(mixed_id) %>% 
  mutate(across(
    c(score_aisance_financiere, score_urbanisme),
    ~ mean(.x, na.rm = T),
    .names = "sr_{.col}"
  )) %>% 
  ungroup()

## Points ----

pc_points <- ad_layer_pc %>% 
  left_join(as_tibble(sdr_layer) %>% select(csduid, sdr_name = name), by = c("sdridu" = "csduid")) %>% 
  mutate(sdr_name = case_when(sdr_name %in% c("Québec", "Montréal") ~ sdr_name,
                              TRUE ~ "Ailleurs")) %>% 
  ggplot(aes(x = PC1, y = PC2, color = sdr_name)) +
  geom_point(aes(x = PC1, y = PC2, color = sdr_name, fill = sdr_name), alpha = 0.03) +
  stat_ellipse(level = 0.98) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-5, 5)) +
  scale_color_brewer(type = "div") +
  scale_fill_brewer(type = "div", guide = "none") +
  labs(
    y = "Axe d'urbanisme",
    x = "Axe d'aisance financière",
    color = ""
  )



## Densité ----

dens_aisance <- ad_layer_pc %>% 
  left_join(as_tibble(sdr_layer) %>% select(csduid, sdr_name = name), by = c("sdridu" = "csduid")) %>% 
  mutate(sdr_name = case_when(sdr_name %in% c("Québec", "Montréal") ~ sdr_name,
                              TRUE ~ "Ailleurs")) %>% 
  ggplot(aes(x = score_aisance_financiere, y = sdr_name, fill = sdr_name, color = sdr_name)) +
  geom_density_ridges(alpha = 0.7) +
  labs(
    title = "Score d'aisance financière",
    y = "",
    x = ""
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_brewer(type = "div") +
  scale_fill_brewer(type = "div", guide = "none") +
  theme(legend.position = "None")

dens_urbain <- ad_layer_pc %>% 
  left_join(as_tibble(sdr_layer) %>% select(csduid, sdr_name = name), by = c("sdridu" = "csduid")) %>% 
  mutate(sdr_name = case_when(sdr_name %in% c("Québec", "Montréal") ~ sdr_name,
                              TRUE ~ "Ailleurs")) %>% 
  ggplot(aes(x = score_urbanisme, y = sdr_name, fill = sdr_name, color = sdr_name)) +
  geom_density_ridges(alpha = 0.7) +
  labs(
    title = "Score d'urbanisme",
    y = "",
    x = ""
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_brewer(type = "div") +
  scale_fill_brewer(type = "div", guide = "none") +
  theme(legend.position = "None")

dens_scores <- cowplot::plot_grid(dens_aisance, dens_urbain)



## Map PCs ----

filter_within <- function(ad_layer, ad_id, dist_km){
  ad_layer <- ad_layer %>% 
    mutate(
      cent = st_centroid(shape)
    )
  
  cent_ad_select <- ad_layer %>% 
    filter(adidu == ad_id) %>% 
    select(cent) %>% 
    pull()
  
  ad_layer %>% 
    mutate(dist_to_ad = st_distance(cent, cent_ad_select)) %>% 
    filter(dist_to_ad <= units::as_units(dist_km, "km")) %>% 
    "$"(adidu)
  
}

make_pc_map <- function(layer, pc, full_layer = layer, center_on_id = NULL, dist_km = NULL, fill_label = ""){
  
  if(!is.null(center_on_id)){
    layer <- layer %>% 
      filter(adidu %in% filter_within(full_layer, center_on_id, dist_km))
  }
  
  
  layer %>% 
    mutate(shape = st_transform(shape, crs = st_crs(3857L))) %>% 
    ggplot() +
    annotation_map_tile(
      zoomin = 1L,
      type = "cartolight",
    ) +
    geom_sf(aes(geometry = shape, fill = {{pc}}), alpha = 0.75, color = NA) +
    scale_fill_gradient2(low = "#5f89e3", 
                         mid = "white", 
                         high = "#FF7F7F", 
                         midpoint = 50) +
    labs(fill = fill_label, x = "", y = "") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0),
          legend.position = "None")
}

mtl_aisance <- make_pc_map(ad_layer_pc, sr_score_aisance_financiere, ad_layer, "24663370", 6) + labs(title = "Score d'aisance financière")
mtl_urbanisme <- make_pc_map(ad_layer_pc, sr_score_urbanisme, ad_layer, "24663370", 6) + labs(title = "Score d'urbanisme")
legend <- cowplot::get_legend(mtl_aisance + theme(legend.box.margin = margin(0, 0, 0, 12),
                                                  legend.position = "right"))

mtl_scores <- cowplot::plot_grid(mtl_aisance, mtl_urbanisme, legend, ncol = 3, rel_widths = c(3, 3, .6))

excl_fleuve <- c("24231037", "24230223", "24250008", "24250259", "24230055", "24250281", "24250195", "24250225", "24250226", "24230066", "24250220")

qc_aisance <- make_pc_map(ad_layer_pc %>% filter(!adidu %in% excl_fleuve), sr_score_aisance_financiere, ad_layer, "24230175", 8) + labs(title = "Score d'aisance financière")
qc_urbanisme <- make_pc_map(ad_layer_pc %>% filter(!adidu %in% excl_fleuve), sr_score_urbanisme, ad_layer, "24230175", 8) + labs(title = "Score d'urbanisme")

qc_scores <- cowplot::plot_grid(qc_aisance, qc_urbanisme, legend, ncol = 3, rel_widths = c(3, 3, .6))








