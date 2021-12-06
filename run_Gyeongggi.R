set.seed(19940827)
setwd(dirname(rstudioapi::getActiveDocumentContext()[["path"]]))
source("supportfunctions.R")
source("MCMC_learn.R")
source("MCMC_predict.R")
library(sp) # https://cran.r-project.org/web/packages/sp/vignettes/intro_sp.pdf
library(data.table)
library(ggplot2)
library(rgdal)
load("data/GyeonggiMap.Rdata")
Sys.setlocale("LC_ALL", "korean")


#### prepare data
DF = as.data.frame(fread("data/GyeonggiFD_pm10.csv", encoding = "UTF-8"))
dim(DF) # replications of responses
Y.obs = as.matrix(DF[, -c(1, ncol(DF)-1, ncol(DF))])
dim(Y.obs)
coords.train = as.matrix(DF[, c(ncol(DF), ncol(DF)-1)])
dim(coords.train)

std_df = function(X){
  means = colMeans(X)
  sds = apply(X, 2, sd)
  X = sweep(X, 2, means, "-")
  X = sweep(X, 2, sds, "/")
  return(list(std_X = X, means = means, sds = sds))
}
X.std.info = std_df(coords.train)
X = matrix(rep(1, nrow(Y.obs)), ncol = 1)
X = cbind(X, X.std.info$std_X)
dim(X)
Xs = array(rep(X, ncol(Y.obs)), dim = c(dim(X), ncol(Y.obs))) # replications of covariates
dim(Xs)

#### MCMC
# subtract site-specific mean
# Y.obs.2 = sweep(Y.obs,1, rowMeans(Y.obs, na.rm = T), "-")

RES = SDPlearn(Y.obs, Xs, coords.train, phigrid_d = 0.1, B = 3e3, burnin = 1e3)
save(RES, file = "Gyg_RES.Rdata")

png("Gyg_traceplots.png", width = 600, height = 300)
B = 3e3; burnin = 1e3
traceidx = 1:(B-burnin)
opar = par(no.readonly = T)
par(mfrow = c(2, 4))
par(mar = c(2.5, 3, 2, 1), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0))
plot(RES$BETA.samp[1, traceidx], type = "l")
abline(h = mean(RES$BETA.samp[1, traceidx]), col = "red")
plot(RES$BETA.samp[2, traceidx], type = "l")
abline(h = mean(RES$BETA.samp[2, traceidx]), col = "red")
plot(RES$BETA.samp[3, traceidx], type = "l")
abline(h = mean(RES$BETA.samp[3, traceidx]), col = "red")
plot(RES$TAU.samp[traceidx], type = "l")
abline(h = mean(RES$TAU.samp[traceidx]), col = "red")
plot(RES$NU.samp[traceidx], type = "l")
abline(h = mean(RES$NU.samp[traceidx]), col = "red")
plot(RES$PHI.samp[traceidx], type = "l")
abline(h = mean(RES$PHI.samp[traceidx]), col = "red")
plot(RES$SIG.samp[traceidx], type = "l")
abline(h = mean(RES$SIG.samp[traceidx]), col = "red")
par(opar)
dev.off()

# grid locations
grid.mat = as.matrix(fread("data/GyeonggiFD_grid.csv", encoding = "UTF-8"))
X_grid = sweep(grid.mat, 2, X.std.info$means, "-")
X_grid = sweep(X_grid, 2, X.std.info$sds, "/")
X_grid = cbind(1, X_grid)
head(X_grid)
X_grids = array(rep(X_grid, ncol(Y.obs)), dim = c(dim(X_grid), ncol(Y.obs))) # replications of covariates
RESgrid = SDPpredict(RES, X_grids, coords.train, grid.mat)
save(RES, file = "Gyg_RESgrid.Rdata")

#### predict
RESPRED = function(RESpred){
  res = matrix(0, nrow = dim(RESpred)[1], dim(RESpred)[2])
  for(i in 1:dim(RESpred)[3]){
    res = res + RESpred[,,i] 
  }
  return(res / dim(RESpred)[3])
}

RESgridmean = RESPRED(RESgrid$Y.pred.samp)
REScompmean = RESPRED(RES$Y.comp.samp)

idx = sample.int(92,1)
plotdf2 = data.frame(lon = coords.train[,1],
                     lat = coords.train[,2],
                     y = REScompmean[, idx])
obsplot = GyeonggiMap + 
  labs(title = "Observed Fine Dust", x = "", y = "") +
  geom_point(data = plotdf2, aes(x = lon, y = lat, col = y), size = 3, inherit.aes = F, alpha = 0.8)+
  guides(col = guide_legend(title = "y"))+
  theme(legend.position = "right")+
  scale_color_gradient(low = "green", high = "red")
obsplot

idx = sample.int(92,1)
plotdf1 = data.frame(lon = grid.mat[,1],
                     lat = grid.mat[,2],
                     y = RESgridmean[,idx])
gridplot1 = GyeonggiMap + 
  labs(title = "Krigged Fine Dust (1)", x = "", y = "") +
  geom_point(data = plotdf1, aes(x = lon, y = lat, col = y), size = 3, inherit.aes = F, alpha = 0.8)+
  guides(col = guide_legend(title = "y"))+
  theme(legend.position = "right")+
  scale_color_gradient(low = "green", high = "red")+
  geom_point(data=plotdf2, aes(x=lon, y= lat), size = 3, inherit.aes=F, shape=3)
gridplot1

idx = sample.int(92,1)
plotdf1 = data.frame(lon = grid.mat[,1],
                     lat = grid.mat[,2],
                     y = RESgridmean[,idx])
gridplot2 = GyeonggiMap + 
  labs(title = "Krigged Fine Dust (2)", x = "", y = "") +
  geom_point(data = plotdf1, aes(x = lon, y = lat, col = y), size = 3, inherit.aes = F, alpha = 0.8)+
  guides(col = guide_legend(title = "y"))+
  theme(legend.position = "right")+
  scale_color_gradient(low = "green", high = "red")+
  geom_point(data=plotdf2, aes(x=lon, y= lat), size = 3, inherit.aes=F, shape=3)
gridplot2

save(obsplot, gridplot1, gridplot2, file = "Gyg_plots.Rdata")