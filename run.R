set.seed(19940827)
setwd(dirname(rstudioapi::getActiveDocumentContext()[["path"]]))
source("simdat.R")

#### prepare data
dim(Y.obs) # replications of responses
X = matrix(rep(1, nrow(Y.obs)), ncol = 1)
X = cbind(X, coords.train)
head(X)
Xs = array(rep(X, ncol(Y.obs)), dim = c(dim(X), ncol(Y.obs))) # replications of covariates

# test locations
dim(Y.pred)
X_pred = matrix(rep(1, nrow(Y.pred)), ncol = 1)
X_pred = cbind(X_pred, coords.pred)
head(X_pred)
X_preds = array(rep(X_pred, ncol(Y)), dim = c(dim(X_pred), ncol(Y))) # replications of covariates

#### MCMC
source("MCMC_learn.R")
B = 3e3; burnin = 1e3
RES = SDPlearn(Y.obs, Xs, coords.train, phigrid_d = 0.1, B = B, burnin = burnin)

png("Sim_traceplots.png", width = 600, height = 300)
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

#### predict
source("MCMC_predict.R")
RESpred = SDPpredict(RES, X_preds, coords.train, coords.pred)
RESpred = RESpred$Y.pred.samp

png("Sim_imputations.png", width = 800, height = 300)
opar = par(no.readonly = T)
par(mfrow = c(1, 3))
par(mar = c(2.5, 3, 2, 1), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0))
hist(Y.pred[1,], breaks = 20, prob = T)
hist(RESpred[1,,1], add =T, breaks = 20, col = "blue", alpha = 0.3, prob = T)
hist(Y.pred[2,], breaks = 20, prob = T)
hist(RESpred[2,,1], add =T, breaks = 20, col = "blue", alpha = 0.3, prob = T)
hist(Y.pred[3,], breaks = 20, prob = T)
hist(RESpred[3,,1], add =T, breaks = 20, col = "blue", alpha = 0.3, prob = T)
par(opar)
dev.off()
