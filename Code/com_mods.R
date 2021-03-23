## Purpose: Run community-level models using predicted occupancy from the MSOM
## Project: High_Severity_Birds
## Author: Zack Steel

com_models <- function() 
{
  library(tidyverse)
  library(brms)
  library(rethinking)
  
  load("data/Results/MSOM/z_mat.RData")
  
  ## Get point-wise estimates of richness
  rich_draws <- apply(z_draws, c(1,2), sum)
  rich_mn <- apply(rich_draws, 2, mean)
  rich_md <- apply(rich_draws, 2, median)
  rich_sd <- apply(rich_draws, 2, sd)
  rich_hpdi <- apply(rich_draws, 2, HPDI, prob = 0.9)
  
  ## Of dispersion
  disp_draws <- read.csv("Data/Results/MSOM/dis_c.csv")
  disp_mn <- apply(disp_draws, 2, mean)
  disp_md <- apply(disp_draws, 2, median)
  disp_sd <- apply(disp_draws, 2, sd)
  disp_hpdi <- apply(disp_draws, 2, HPDI, prob = 0.9)
  
  ## Richness by guild
  nest <- read.csv("Data/Reference/MSOM_guilds.csv") %>% 
    dplyr::select(sp = species, common, nest = Nest_Site)
  
  spp <- merge(spp.tab, nest)
  ## Pull out cavity (primary and secondary), 
  ## ground, shrub and tree nesters
  cavid <- filter(spp, nest %in% c("cavity_primary", "cavity_secondary")) %>% 
    pull(s)
  cavid1 <- filter(spp, nest %in% c("cavity_primary")) %>% 
    pull(s)
  cavid2 <- filter(spp, nest %in% c("cavity_secondary")) %>% 
    pull(s)
  groid <- filter(spp, nest == "ground") %>% 
    pull(s)
  shrid <- filter(spp, nest == "shrub") %>% 
    pull(s)
  treid <- filter(spp, nest == "tree") %>% 
    pull(s)
  
  crich1 <- z_draws[,,cavid1] %>% 
    apply(c(1,2), sum)
  crich1_mn <- apply(crich1, 2, mean)
  crich1_md <- apply(crich1, 2, median)
  crich1_sd <- apply(crich1, 2, sd)
  crich1_hpdi <- apply(crich1, 2, HPDI, prob = 0.9)
  
  crich2 <- z_draws[,,cavid2] %>% 
    apply(c(1,2), sum)
  crich2_mn <- apply(crich2, 2, mean)
  crich2_md <- apply(crich2, 2, median)
  crich2_sd <- apply(crich2, 2, sd)
  crich2_hpdi <- apply(crich2, 2, HPDI, prob = 0.9)
  
  grich <- z_draws[,,groid] %>% 
    apply(c(1,2), sum)
  grich_mn <- apply(grich, 2, mean)
  grich_md <- apply(grich, 2, median)
  grich_sd <- apply(grich, 2, sd)
  grich_hpdi <- apply(grich, 2, HPDI, prob = 0.9)
  
  srich <- z_draws[,,shrid] %>% 
    apply(c(1,2), sum)
  srich_mn <- apply(srich, 2, mean)
  srich_md <- apply(srich, 2, median)
  srich_sd <- apply(srich, 2, sd)
  srich_hpdi <- apply(srich, 2, HPDI, prob = 0.9)
  
  trich <- z_draws[,,treid] %>% 
    apply(c(1,2), sum)
  trich_mn <- apply(trich, 2, mean)
  trich_md <- apply(trich, 2, median)
  trich_sd <- apply(trich, 2, sd)
  trich_hpdi <- apply(trich, 2, HPDI, prob = 0.9)
  
  ## Set up dataframe for modeling
  d <- data.frame(rich_mn, rich_md, rich_sd, 
                  rich_low = rich_hpdi[1,], rich_hi = rich_hpdi[2,],
                  
                  crich1_mn, crich1_md, crich1_sd, 
                  crich1_low = crich1_hpdi[1,], crich1_hi = crich1_hpdi[2,],
                  
                  crich2_mn, crich2_md, crich2_sd, 
                  crich2_low = crich2_hpdi[1,], crich2_hi = crich2_hpdi[2,],
                  
                  grich_mn, grich_md, grich_sd, 
                  grich_low = grich_hpdi[1,], grich_hi = grich_hpdi[2,],
                  
                  srich_mn, srich_md, srich_sd, 
                  srich_low = srich_hpdi[1,], srich_hi = srich_hpdi[2,],
                  
                  trich_mn, trich_md, trich_sd, 
                  trich_low = trich_hpdi[1,], trich_hi = trich_hpdi[2,],
                  
                  disp_mn, disp_md, disp_sd,
                  disp_low = disp_hpdi[1,], disp_hi = disp_hpdi[2,],
                  
                  locs = stan_data$locs,
                  patches = stan_data$patches,
                  fires = stan_data$fires,
                  dist = stan_data$dist,
                  ysf = stan_data$ysf,
                  ysfsq = stan_data$ysfsq,
                  dy = stan_data$dy,
                  parea = stan_data$parea,
                  py = stan_data$py,
                  elev = stan_data$elev,
                  elevsq = stan_data$elevsq,
                  lat = stan_data$lat)
  
  m_rich <- brm(rich_mn | se(rich_sd, sigma = TRUE) ~ 
                  dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                  (1|locs) + (1|patches) + (1|fires), 
                family = gaussian(link = "log"),
                data = d, chains = 4, cores = 4, iter = 1000)
  
  mc1_rich <- brm(crich1_mn | se(crich1_sd, sigma = TRUE) ~ 
                   dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                   (1|locs) + (1|patches) + (1|fires), 
                  family = gaussian(link = "log"),
                 data = d, chains = 4, cores = 4, iter = 1000)
  
  mc2_rich <- brm(crich2_mn | se(crich2_sd, sigma = TRUE) ~ 
                    dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                    (1|locs) + (1|patches) + (1|fires), 
                  family = gaussian(link = "log"),
                  data = d, chains = 4, cores = 4, iter = 1000)
  
  mg_rich <- brm(grich_mn | se(grich_sd, sigma = TRUE) ~ 
                   dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                   (1|locs) + (1|patches) + (1|fires), 
                 family = gaussian(link = "log"),
                 data = d, chains = 4, cores = 4, iter = 1000)
  
  ms_rich <- brm(srich_mn | se(srich_sd, sigma = TRUE) ~ 
                   dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                   (1|locs) + (1|patches) + (1|fires), 
                 family = gaussian(link = "log"),
                 data = d, chains = 4, cores = 4, iter = 1000)
  
  mt_rich <- brm(trich_mn | se(trich_sd, sigma = TRUE) ~ 
                   dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                   (1|locs) + (1|patches) + (1|fires), 
                 family = gaussian(link = "log"),
                 data = d, chains = 4, cores = 4, iter = 1000)
  
  ## Dispersion model (not splitting by guild)
  m_disp <- brm(disp_mn | se(disp_sd, sigma = TRUE) ~ 
                  dist + ysf + ysfsq + dy + parea + py + elev + elevsq + lat +
                  (1|locs) + (1|patches),# + (1|fires), 
                family = gaussian(link = "log"),
                data = d, chains = 4, cores = 4, iter = 1000)
  
  ## Save models
  save(d, m_rich, mc1_rich, mc2_rich, mg_rich, ms_rich, mt_rich, m_disp,
       file = "Models/community.RData")
  
}
