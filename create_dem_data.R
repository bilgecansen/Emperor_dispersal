
library(tidyverse)
library(foreach)

# Colony size and growth
N_chains <- readRDS("data/N_chains.rds")
sites_update <- read.csv("data/colony_attributes_update.csv")

N_mean_chains <- foreach(i = 1:50, .combine = "cbind") %do% {
  
  z <- t(apply(N_chains[,seq(i, 500, 50)], 1, function(x) {
    mean(log(x[which(x != 0)]))
  }))
  
  as.vector(z)
}

N_mean <- apply(N_mean_chains, 2, mean)

r_mean_chains <- foreach(i = 1:50, .combine = "cbind") %do% {
  
  z1 <- t(apply(N_chains[,seq(i, 500, 50)], 1, function(x) {
    z2 <- log(x[2:10]/x[1:9])
    z2 <- z2[!is.infinite(z2)]
    mean(z2[!is.nan(z2)])
  }))
  
  as.vector(z1)
}

r_mean <- apply(r_mean_chains, 2, mean, na.rm = T)

data_dem <- data.frame(size = N_mean,
                       growth = r_mean,
                       site_id = sites_update$site_id)

saveRDS(data_dem, "data/data_dem.rds")
