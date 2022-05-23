library(tidymodels)
library(tidyverse)

set.seed(12)

nsim <- 8

points <- tibble(
  a = runif(nsim, 0, 5),
  b = a * 2 / 5 + rnorm(nsim, 0, 0.5) + 2
)

pca_res <- points %>% 
  recipe(~a + b) %>% 
  step_pca(all_numeric(), 
           options = list(center = TRUE, scale. = TRUE)) %>% 
  prep()

points <- bind_cols(
  points, 
  points %>% 
    recipe(~a + b) %>% 
    step_center(all_numeric()) %>% 
    prep() %>% 
    bake(new_data = NULL) %>% 
    rename(
      a_centered = a,
      b_centered = b
    ),
  points %>% 
    recipe(~a + b) %>% 
    step_normalize(all_numeric()) %>% 
    prep() %>% 
    bake(new_data = NULL) %>% 
    rename(
      a_scaled = a,
      b_scaled = b
    ),
  pca_res %>% 
    bake(new_data = NULL)
) %>% arrange(a_scaled)




write_csv(points, "./data/pca_example.csv")

points %>%
  ggplot(aes(x = a, y = b)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 8))

points %>%
  ggplot(aes(x = a_scaled, y = b_scaled)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  stat_smooth(method = "lm", se = F)
# 
# pca_res$steps[[1]]$res
# 
 prcomp(points %>% select(a, b), center = T, scale. = F)




