setwd(dirname(rstudioapi::getActiveDocumentContext()[["path"]]))
## simulate spatial dataset with R pkg "sp"
library(sp) # https://cran.r-project.org/web/packages/sp/vignettes/intro_sp.pdf
library(data.table)
library(ggplot2)
library(rgdal)
source("supportfunctions.R")
load("data/GyeonggiMap.Rdata")
Sys.setlocale("LC_ALL", "korean")



# Toy data
xy.df = fread("data/GyeonggiFD_coords.csv", encoding = "UTF-8")
coords.all = cbind(xy.df$lon, xy.df$lat)
colnames(coords.all) = c("lon", "lat")

xy.df2 = xy.df
xy.df2$y = rnorm(nrow(xy.df))*10
GyeonggiMap + 
  labs(title = "Measurement sites", x = "", y = "") +
  geom_point(data = xy.df2, aes(x = lon, y = lat, col = y), inherit.aes = F, alpha = 0.8, size = 2)+
  guides(col = guide_legend(title = "y"))+
  theme(legend.position = "right")


# simulate true data with nonstationary GP
Beta = matrix(c(1, 0.1, -0.1), ncol = 1)
sigma = 0.05
rho = 1
tau = 0.5
X = coords.all

std_df = function(X){
  means = colMeans(X)
  sds = apply(X, 2, sd)
  X = sweep(X, 2, means, "-")
  X = sweep(X, 2, sds, "/")
  return(list(std_X = X, means = means, sds = sds))
}
X.std.info = std_df(coords.all)


Y = nonstGP(60, X, coords.all,  Beta, sigma, rho, tau)

y = Y[,1]
xy.df2 = xy.df
xy.df2$y = y
yplot = GyeonggiMap + 
  labs(title = "Nonstationary GP", x = "", y = "") +
  geom_point(data = xy.df2, aes(x = lon, y = lat, col = y), inherit.aes = F, alpha = 0.8, size = 2)+
  guides(col = guide_legend(title = "y"))+
  theme(legend.position = "right")
yplot

# generate train-test split
site_idx = sample.int(nrow(Y), 5)
Y.train = Y[-site_idx,]
coords.train = coords.all[-site_idx, ]
Y.pred = Y[site_idx, ]
coords.pred = coords.all[site_idx, ]

# generate obs-missing split within train
site_idx = sample.int(nrow(Y.train), 20)
time_idx = sample.int(ncol(Y.train), 50)
Y.obs = Y.train
Y.obs[site_idx, time_idx] = NA_real_
sum(is.na(Y.obs))

Y.miss = Y.train[site_idx, time_idx]

# grid to be krigged
grid.df = fread("data/GyeonggiFD_grid.csv", encoding = "UTF-8")

