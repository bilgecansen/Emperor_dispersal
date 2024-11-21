
library(foreach)
library(tidyverse)
library(party)
library(permimp)
library(caret)
library(pROC)
library(R.matlab)
library(patchwork)
library(pdp)
library(MCMCvis)


index1 <- c(1:18, 20:26, 28, 31:34, 38:39, 41, 43:46, 48:52, 54:55, 57:59, 
            61:62, 64)
index2 <- c(1:18, 19:25, 28, 30:33, 34:35, 37, 38:41, 42:46, 47:48, 49:51, 
            52:54)

empe_sites <- read.csv("data/empe_sitesNewNB.csv")
sites <- empe_sites$site_id[index1]

data_pop <- readRDS("data/data_pop_empe.rds")

data_emg_year <- 
  read.csv("data/Data_emigration_median_per_year.csv", header = F) %>%
  as.matrix()
data_emg_year <- data_emg_year[index2,]

data_emg <- read.csv("data/Data_emigration_median.csv", header = F)
data_emg <- data_emg[index2,]


# Load environmental data -------------------------------------------------

# ESM data
data_esm <- readRDS("data/data_env_empe.rds") %>%
  select(-contains(c("phtc_",
                     "zooc_incubation", "zooc_arrival",
                     "hmxl_incubation", "hmxl_arrival",
                     "aice_"))) %>%
  filter(year >= 2009 & year <= 2013) %>%
  filter(site_id %in% sites)

data_faice <- readRDS("data/data_faice.rds")[[6]] %>%
  filter(year >= 2009 & year <= 2013) %>%
  filter(site_id %in% sites)

data_fdice <- readRDS("data/data_fdice.rds") %>%
  filter(year >= 2009 & year <= 2013) %>%
  filter(site_id %in% sites)

data_ftice <- readRDS("data/data_ftice.rds")[[6]] %>%
  filter(year >= 2009 & year <= 2013) %>%
  filter(site_id %in% sites)

data_env <- left_join(data_esm, data_faice, by = c("site_id", "year")) %>%
  left_join(data_fdice, by = c("site_id", "year"))  %>%
  left_join(data_ftice, by = c("site_id", "year"))


# Random forests with environmental data ----------------------------------

# Environmental data
data_rf <- data_env %>%
  select(-site_id, -year)
data_rf$r <- c(t(data_emg_year))
data_rf$r <- ifelse(data_rf$r > 0, 1, 0) %>%
  as.factor()

nvar <- ncol(data_rf)-1
mtry <- round(nvar*seq(0.1,0.9, 0.1))

auc_all <- foreach(i = 1:length(mtry), .combine = "c") %do% {
  
  foreach(h = 1:10, .combine = "mean") %do% {
    
    idx <- sample(1:250)
    data_rf_test <- data_rf[idx[1:175],]
    data_rf_tr <- data_rf[-idx[1:175],]
    
    rf <- cforest(r ~ ., data = data_rf_tr, 
                  controls = cforest_unbiased(ntree = 2000, mtry = mtry[i]))
    pr <- predict(rf, newdata = data_rf_test, OOB = TRUE, type = "prob")
    pr <- do.call(rbind, pr)
    
    auc(data_rf_test$r, pr[,2])[1]
  }
}

set.seed(101)
rf_main <- cforest(r ~ ., data = data_rf, 
                   controls = cforest_unbiased(ntree = 8000, 
                                               mtry = mtry[which.max(auc_all)]))

varimp_rf <- permimp(rf_main, conditional = T, progressBar = T, AUC = T)
varimp_rf <- varimp_rf$values[order(varimp_rf$values, decreasing = T)]
varimp_rf

# Partial Dependence plots
par1 <- partial(rf_main, pred.var = "zooc_nonbreed", 
                type = "classification", which.class = 2, 
                rug = T, prob = T, progress = "text")

par2 <- partial(rf_main, pred.var = "fdice_rearing", 
                type = "classification", which.class = 2, 
                rug = T, prob = T, progress = "text")


# Random forests with demography ------------------------------------------

# Colony size and growth
data_dem <- readRDS("data/data_dem.rds")

# Blinking
bl <- readMat("data/mean_blinking_10years.mat")[[1]]
data_blink <- data.frame(blink = bl,
                         site_id = sites)

# Environment
data_env_avg <- group_by(data_env, site_id) %>%
  summarise(across(-year, mean)) %>%
  ungroup() 

# Average emigration
dz <- ifelse(data_emg_year > 0 , 1, 0)
z <- apply(dz, 1, sum)
data_emg_avg <- data.frame(e = z/5,
                           site_id = sites)

# Site level fast ice data
data_sites <- readRDS("data/data_icef_sites.rds") %>%
  select(site_id, contains("_100"), -contains(c("icef", "time")))

# Combine all data
data_rf2 <- left_join(data_env_avg, data_sites, by = "site_id") %>%
  left_join(data_blink, by = "site_id")  %>%
  left_join(data_dem, by = "site_id") %>%
  left_join(data_emg_avg, by = "site_id") 

set.seed(102)
rf_main2 <- cforest(e ~ ., data = data_rf2[,-1], 
                    controls = cforest_unbiased(ntree = 8000, mtry = 15))

cforestStats(rf_main2)
varimp_rf2 <- permimp(rf_main2, conditional = T, progressBar = T)
varimp_rf2 <- varimp_rf2$values[order(varimp_rf2$values, decreasing = T)]
varimp_rf2

par_d1 <- partial(rf_main2, pred.var = "size", 
                  rug = T, progress = "text")

rf_results <- list(rf_main = rf_main,
                   rf_main2 = rf_main2,
                   par1 = par1,
                   par2 = par2,
                   par_d1 = par_d1)

saveRDS(rf_results, "results/results_rf_disp.rds")


# Visualizations ----------------------------------------------------------

theme_set(theme_bw())

dat_rug1 <- filter(data_rf, r == 0)
dat_rug2 <- filter(data_rf, r == 1)
g1 <- ggplot() +
  stat_smooth(data = par1, aes(x = zooc_nonbreed, y = yhat), geom = "line", 
              se = F, size = 1,
              col = "darkblue", alpha = 0.8) +
  geom_rug(data = dat_rug1, mapping = aes(x = zooc_nonbreed), size = 0.1) +
  geom_rug(data = dat_rug2, mapping = aes(x = zooc_nonbreed), sides = "top", 
           size = 0.1) +
  labs(y = "Annual Emigration Probability", 
       x = "Zooplankton Biomass (Nonbreeding)") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0.29, 0.55))

g1
ggsave("fig3a.pdf", width = 10, height = 10, units = "cm", dpi = 600)

g2 <- ggplot() +
  stat_smooth(data = par2, aes(x = fdice_rearing, y = yhat), geom = "line", 
              size = 1, alpha = 0.8, col = "darkblue") +
  labs(y = "Annual Emigration Probability", 
       x = "Distance to nearest fast ice edge (Rearing)") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0.29, 0.55))

g2
ggsave("fig3b.pdf", width = 10, height = 10, units = "cm", dpi = 600)

g3 <- ggplot() +
  stat_smooth(data = par_d1, aes(x = size, y = yhat), geom = "line", size = 1, 
              alpha = 0.8, col = "darkblue") +
  labs(y = "Average Emigration Probability", 
       x = "Blinking Probability") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0.29, 0.55))

g3
ggsave("fig3c.pdf", width = 10, height = 10, units = "cm", dpi = 600)

g4 <- ggplot(mapping = aes(x = factor(names(varimp_rf2), 
                                      levels = names(varimp_rf2)),
                           y = varimp_rf2)) +
  geom_col(fill = "dark blue") +
  labs(y = "Variable Importance") +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8))
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 10))

g4
ggsave("fig3d_alt.pdf", width = 10, height = 10, units = "cm", dpi = 600)

g5 <- ggplot(mapping = aes(x = factor(names(varimp_rf), 
                                      levels = names(varimp_rf)),
                           y = varimp_rf)) +
  geom_col(fill = "dark blue") +
  labs(y = "Variable Importance") +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8))
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 10))

g5
ggsave("figS4.pdf", width = 6, height = 6, units = "cm", dpi = 600)

#(g1 + g2 + g3) +
#plot_annotation(tag_levels = "a")

#ggsave("fig4new.pdf", width = 24, height = 8, units = "cm", dpi = 600)

# Variable importance vs cor with Ice area
data_ice <- readRDS("data_env_empe.rds") %>%
  select(site_id, year, contains("aice_")) %>%
  filter(year >= 2009 & year <= 2013) %>%
  filter(site_id %in% sites)

cor(data_ice$aice_nonbreed, select(data_rf, -r)) 
cor(data_ice$aice_laying, select(data_rf, -r))

varimp_rf3 <- varimp_rf[order(names(varimp_rf))]
cor1 <- cor(data_ice$aice_nonbreed, select(data_rf, -r))[1,]
cor1 <- cor1[order(names(cor1))]
cor2 <- cor(data_ice$aice_laying, select(data_rf, -r))[1,]
cor2 <- cor2[order(names(cor2))]
cor3 <- cor(data_ice$aice_incubation, select(data_rf, -r))[1,]
cor3 <- cor3[order(names(cor3))]
cor4 <- cor(data_ice$aice_laying, select(data_rf, -r))[1,]
cor4 <- cor4[order(names(cor4))]

gc1 <- ggplot() +
  geom_point(mapping = aes(x = abs(cor1), y = varimp_rf3), size = 2) +
  labs(x = "Absolute Correlation with Ice Area (Nonbreeding)", 
       y = "Variable Importance") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

gc2 <- ggplot() +
  geom_point(mapping = aes(x = abs(cor2), y = varimp_rf3), size = 2) +
  labs(x = "Absolute Correlation with Ice Area (Laying)", 
       y = "Variable Importance") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

gc3 <- ggplot() +
  geom_point(mapping = aes(x = abs(cor3), y = varimp_rf3), size = 2) +
  labs(x = "Absolute Correlation with Ice Area (Incubation)", 
       y = "Variable Importance") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

gc4 <- ggplot() +
  geom_point(mapping = aes(x = abs(cor4), y = varimp_rf3), size = 2) +
  labs(x = "Absolute Correlation with Ice Area (Rearing)", 
       y = "Variable Importance") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

(gc1 + gc2) / (gc3 + gc4)
ggsave("plot_disp_imp_vs_cor.jpeg", width = 24, height = 24, units = "cm", 
       dpi = 600)
