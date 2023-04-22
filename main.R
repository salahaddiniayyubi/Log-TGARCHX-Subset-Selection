# Libraries:

suppressPackageStartupMessages({
library(xtable)
library(latex2exp)
library(Boruta)
library(ranger)
library(readxl)
library(lgarch)
library(xts)
library(readr)
library(zoo)
library(fpp2)
library(glmnet)
library(caret)
library(doParallel)
library(abess)
library(scales)
library(MCS)
})


############################################################################################
#                             1. Data:                                                     #
############################################################################################

options(future.globals.maxSize = 7680000 * 1024^2)
# setwd("Log-GARCH_Subset-Selection"
# )
rm(list = ls()
)
start_time <- Sys.time()
# Data and Settings:
# setwd("/share/castor/home/orujov/NewLgarchx/Log-TGARCHX-Subset-Selection")
data <- read_csv("data/data_stationary.csv",
                 col_types = cols(...1 = col_date(format = "%Y-%m-%d")))
colnames(data)[1] <- "Date"
dates <- as.Date(data$Date, format = "%m-%d-%Y")
data <- data[, !(names(data) %in% c("Date", "dl_BTCUSD"))]
data <- as.xts(data, order.by = dates)
data_scaled <- data.frame(row.names = dates)
names(data) <- gsub(x = names(data), pattern = "=|-|\\^", replacement = "")

############################################################################################
#                             2. LASSO:                                                    #
############################################################################################


test_len <- val_len <- round(0.2 * NROW(data))
data_scaled <- scale(data)
data_scaled$asym <- scale(data$returns^2) * (data$returns < 0)

# CV data:
cv_dat <- xts(glag(data_scaled[seq(1, NROW(data_scaled) - test_len), ],
                  pad.value = colMeans(data_scaled)))
names(cv_dat) <- names(data_scaled)
# Create error series:
cv_dat$errs <- xts(arima(log(scale(data$returns[index(cv_dat)],
      scale = FALSE, center = FALSE)^2),
      order = c(1, 0, 1))$residuals, order.by = index(cv_dat))

## Start a parallel cluster:
doParallel::registerDoParallel(cores = detectCores())

## Define the grid of lambda values of lasso and set alpha = 1:
set.seed(777)
trainx <- cv_dat[, !(names(cv_dat) %in% c("returns", "errs", "dl_BTCUSD"))]
trainy <- cv_dat$errs
mod <- glmnet(trainx, trainy, alpha = 1)
grid_reg <- expand.grid(alpha = 1,
                          lambda = mod$lambda)

## Define the control settings of validation:
myTimeControl <- trainControl(method = "timeslice",
                              initialWindow = NROW(data) - val_len - test_len,
                              horizon = 1,
                              fixedWindow = FALSE,
                              allowParallel = TRUE)
glmnet_mod <- train(errs ~ . - returns,
                    intercept = TRUE,
                    data = cv_dat,
                    method = "glmnet",
                    family = "gaussian",
                    trControl = myTimeControl,
                    tuneGrid = grid_reg,
                    metric = 'RMSE')

#getting the coefficients of the final model
coefficients <- coef(glmnet_mod$finalModel, glmnet_mod$bestTune$lambda)
# plot(coefficients)

#create a list of the selected coefficients
vars_lasso <- names(coefficients[which(coefficients != 0),])
vars_lasso <- print(vars_lasso[vars_lasso != "(Intercept)"])
best_lambda <- glmnet_mod$bestTune$lambda

# Predictions
print("start cluster for predictions")
test_set <- createTimeSlices(data$returns,
      initialWindow = NROW(data) - test_len,
            fixedWindow = FALSE
            )
cl <- makeCluster(detectCores())
registerDoParallel(cl)
preds_lasso <- foreach(i = iter(test_set$train),
      .packages = c("lgarch", "xts", "zoo"), .combine = "rbind") %dopar% {
        obs <- index(data)[i]
        ret <- data$returns[i]
        exot <- data_scaled[obs, vars_lasso]
        nams <- names(exot)
        exo <- glag(exot, pad.value = colMeans(exot))
        names(exo) <- nams
        mod <- lgarch(ret, xreg = exo, solve.tol = 1e-21, method = "ls")
        myhat <- fitted(mod, verbose = TRUE)
        a <- tail(cbind(rep(1, length(obs)),
              log(ret^2), myhat[, "lnsd2"], exot) %*%
              coef(mod)[-length(coef(mod))], 1)
        rownames(a) <- as.character(tail(obs, 1) + 1)
        colnames(a) <- "lnsd2"
        print(a)
        return(a)
      }

pred <- as.vector(sqrt(exp(preds_lasso)))
act <- as.vector(sqrt(data$rv_d[rownames(preds_lasso)]))
lasso_acc <- print(forecast::accuracy(act, pred))
rownames(lasso_acc) <- "LASSO"
stopCluster(cl)

############################################################################################
#                             2. ABESS                                                     #
############################################################################################

cv_dat <- as.data.frame(cv_dat
)

cv_set <- createTimeSlices(cv_dat$errs, initialWindow =
          NROW(cv_dat) - val_len, fixedWindow = TRUE
          )
cv_dat <- cv_dat[, !(names(cv_dat) %in% c("returns"))
          ]
ncores <- parallel::detectCores()
mycl <- parallel::makeCluster(ncores)
# i <- cv_set$train$Training1630
doParallel::registerDoParallel(mycl)
x <- foreach(size = 1:(ncol(cv_dat) - 1),
          .combine = "cbind", .packages = c("abess", "xts", "zoo")) %:%
          foreach(i = iter(cv_set$train), .combine = "rbind") %dopar% {
            obs <- index(cv_dat)[i]
            dat <- cv_dat[obs, ]
            fit <- abess(errs ~ ., data = dat, support.size = size)
            test_idx <- tail(index(dat), 1) + 1
            new <- subset(cv_dat[test_idx, ], select = -errs)
            a <- (cv_dat[test_idx, "errs"] - predict(fit, newx = new))^2
            return(a)
          }
stopCluster(mycl)

## Average cross-validation scores for each sparsity
scores <- colMeans(x)
best_size <- print(as.numeric(which.min(scores)))

## Fit the final model with the best sparsity score
fit <- abess(as.matrix(errs) ~ ., data = cv_dat, support.size = best_size)

## Selected variables excluding the lag of log(y^2):
print("vars abess:")
vars_abess <- print(abess::extract(fit, support.size = best_size)$support.vars)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
preds_abess <- foreach(i = iter(test_set$train),
      .packages = c("lgarch", "xts", "zoo"), .combine = "rbind") %dopar% {
        obs <- index(data)[i]
        ret <- data$returns[i]
        exot <- data_scaled[obs, vars_abess]
        nams <- names(exot)
        exo <- glag(exot, pad.value = colMeans(exot))
        names(exo) <- nams
        mod <- lgarch(ret, xreg = exo, solve.tol = 1e-21, method = "ls")
        myhat <- fitted(mod, verbose = TRUE)
        a <- tail(cbind(rep(1, length(obs)),
              log(ret^2), myhat[, "lnsd2"], exot) %*%
              coef(mod)[-length(coef(mod))], 1)
        rownames(a) <- as.character(tail(obs, 1) + 1)
        colnames(a) <- "lnsd2"
        print(a)
        return(a)
      }

pred <- as.vector(sqrt(exp(preds_abess)))
act <- as.vector(sqrt(data$rv_d[rownames(preds_abess)]))
fabess_acc <- print(forecast::accuracy(act, pred))
rownames(fabess_acc) <- "ABESS"
stopCluster(cl)
end_time <- Sys.time()
print(end_time - start_time)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################################################################
#                             3. BORUTA                                                    #
############################################################################################

cv_dat <- as.data.frame(cv_dat
)

# Parameter Grid for random forest:
sensitivity_df <- expand.grid(
  num.trees = c(500, 750, 1000),
  mtry = 1:10,
  block.size = c(7, 14, 21, 28),
  min.node.size = c(5, 10, 20)
)

# Seed:
set.seed(777)

#create and register cluster
mycl <- parallel::makeCluster(detectCores()
)
doParallel::registerDoParallel(cl = mycl
)

#fitting each rf model with different hyperparameters
prediction_error <- foreach(
  num.trees = sensitivity_df$num.trees,
  mtry = sensitivity_df$mtry,
  min.node.size = sensitivity_df$min.node.size,
  block.size = sensitivity_df$block.size, 
  .combine = "c",
  .packages = "rangerts"
) %dopar% {
  #fit model:
  m_i <- rangerts::rangerts(errs ~ ., data = cv_dat,
                 num.trees = num.trees,
                 mtry = mtry,
                 replace = TRUE, # default = T too
                 seed = 1,
                 bootstrap.ts = "moving",
                 block.size = block.size,
                 min.node.size = min.node.size
                 )

  #returning prediction error as percentage
  return(m_i$prediction.error)
}
stopCluster(mycl)

#adding the prediction error column
sensitivity_df$prediction_error <- prediction_error
print(head(sensitivity_df))

# Plotting hyper-parameters:
# ggplot2::ggplot(data = sensitivity_df) + 
#   ggplot2::aes(
#     x = mtry,
#     y = as.factor(min.node.size),
#     fill = prediction_error
#   ) + 
#   ggplot2::facet_wrap(as.factor(sensitivity_df$num.trees)) +
#   ggplot2::geom_tile() + 
#   ggplot2::scale_y_discrete(breaks = c(1, 10, 20)) +
#   ggplot2::scale_fill_viridis_c() + 
#   ggplot2::ylab("min.node.size")

# Best set of hyper-params:
best_hype <- sensitivity_df %>% 
dplyr::arrange(prediction_error) %>% 
dplyr::slice(1)

set.seed = 1
boruta_sel <- Boruta(errs ~ ., data = cv_dat,
      getImp=getImpRfZ, num.trees = best_hype$num.trees,
      mtry = best_hype$mtry, min.node.size = best_hype$min.node.size
      )
borsel <- getSelectedAttributes(boruta_sel
)
print(borsel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
preds_fboruta <- foreach(i = iter(test_set$train),
      .packages = c("lgarch", "xts", "zoo"), .combine = "rbind") %dopar% {
        obs <- index(data)[i]
        ret <- data$returns[i]
        exot <- data_scaled[obs, borsel]
        nams <- names(exot)
        exo <- glag(exot, pad.value = colMeans(exot))
        names(exo) <- nams
        mod <- lgarch(ret, xreg = exo, solve.tol = 1e-21, method = "ls")
        myhat <- fitted(mod, verbose = TRUE)
        a <- tail(cbind(rep(1, length(obs)),
              log(ret^2), myhat[, "lnsd2"], exot) %*%
              coef(mod)[-length(coef(mod))], 1)
        rownames(a) <- as.character(tail(obs, 1) + 1)
        colnames(a) <- "lnsd2"
        print(a)
        return(a)
      }

pred <- as.vector(sqrt(exp(preds_fboruta)))
act <- as.vector(sqrt(data$rv_d[rownames(preds_fboruta)]))
fboruta_acc <- print(forecast::accuracy(act, pred))
rownames(fboruta_acc) <- "BORUTA"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

############################################################################################
#                             4. BENCHMARKS                                                #
############################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mycl <- makeCluster(detectCores())
registerDoParallel(cl = mycl)
preds_lgarch <- foreach(i = iter(test_set$train),
            .packages = c("xts", "zoo", "lgarch")) %dopar% {
            obs <- index(data)[i]
            ret <- scale(data$returns[obs], scale = FALSE)
            mod <- lgarch(ret, method = "ls")
            myhat <- fitted(mod, verbose = TRUE)
            preds <- cbind(rep(1, length(obs)),
                        log(ret^2), myhat[, "lnsd2"]) %*%
                        coef(mod)[-length(coef(mod))
                        ]
            a <- tail(preds, 1)
            rownames(a) <- as.character(tail(obs, 1) + 1)
            colnames(a) <- "lnsd2"
            return(a)
            }
preds_lgarch <- do.call("rbind", preds_lgarch)
act <- as.vector(sqrt(data$rv_d[rownames(preds_lgarch)]))
pred <- as.vector(sqrt(exp(preds_lgarch)))
lgarch_acc <- print(forecast::accuracy(act, pred))
rownames(lgarch_acc) <- "Log - GARCH"
# Fit lgarchx model to 20 percent of the data:
start_time <- Sys.time()
preds_lgarchx <- foreach(i = iter(test_set$train),
            .packages = c("xts", "zoo", "lgarch")) %dopar% {
            obs <- index(data)[i]
            ret <- scale(data$returns[obs], scale = FALSE)
            exot <- scale(data[obs, !(names(data) %in% "returns")])
            nams <- names(exot)
            exo <- glag(exot, pad.value = colMeans(exot))
            colnames(exo) <- nams
            mod <- lgarch(ret, xreg = exo, method = "ls")
            myhat <- fitted(mod, verbose = TRUE)
            preds <- cbind(rep(1, length(obs)),
                        log(ret^2), myhat[, "lnsd2"], exot) %*%
                        coef(mod)[-length(coef(mod))
                        ]
            a <- tail(preds, 1)
            print(a)
            rownames(a) <- as.character(tail(obs, 1) + 1)
            colnames(a) <- "lnsd2"
            return(a)
            }
stopCluster(mycl)
preds_lgarchx <- do.call("rbind", preds_lgarchx)
act <- as.vector(sqrt(data$rv_d[rownames(preds_lgarchx)]))
pred <- as.vector(sqrt(exp(preds_lgarchx)))
lgarchx_acc <- print(forecast::accuracy(act, pred))
rownames(lgarchx_acc) <- "Log - TGARCHX"
accuracies <- print(rbind(lgarch_acc, lgarchx_acc,
            lasso_acc, fabess_acc, fboruta_acc)
            )
end_time <- Sys.time()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################################################################
#                             5. Model Confidence Set                                      #
############################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for (i in c('lgarch', 'lgarchx', 'lasso', 'abess', 'fboruta')) {
  k <- paste0("preds_", i)
  a <- LossVol(sqrt(data$rv_d[rownames(eval(parse(text = k)))]),
  sqrt(exp(eval(parse(text = k)))), "QLIKE")
  colnames(a) <- i
  b <- LossVol(sqrt(data$rv_d[rownames(eval(parse(text = k)))]),
  sqrt(exp(eval(parse(text = k)))), "SE1")
  colnames(b) <- i
  assign(eval(paste0("qlike_",i)), a)
  assign(eval(paste0("mse_",i)), b)
}

losses_qlike <- cbind(qlike_lgarch, qlike_lgarchx, qlike_lasso,
            qlike_abess, qlike_fboruta)
mcs_qlike <- MCSprocedure(losses_qlike)
losses_mse <- cbind(mse_lgarch, mse_lgarchx, mse_lasso,
            mse_abess, mse_fboruta)
mcs_mse <- MCSprocedure(losses_mse)

accuracies <- cbind(accuracies, colMeans(losses_qlike))
colnames(accuracies)[6] <- "QLIKE"

# MCS without lasso and abess:

mcs_qlike_boruta <- MCSprocedure(losses_qlike[, c(1, 2, 5)])
mcs_mse_boruta <- MCSprocedure(losses_qlike[, c(1, 2, 5)])
print(xtable(accuracies, type = "latex"), file = "accuracies.tex")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################################################################
#                             6. Plot                                                      #
############################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# dev.off()
setEPS()
postscript("./images/benchmark.eps")
par(mfrow = c(2, 1), mar = c(2, 4, 2, 2)
)
plot(tail(dates, test_len), sqrt(tail(data$rv_d, test_len)),
     type = "l", ylim = c(0,20), col = "black",
     main = "Log-GARCH (1,1)", xlab = "Date", ylab = TeX(r"($\sigma_t$)")
     )
lines(tail(dates, test_len), sqrt(exp(preds_lgarch)), col = "red"
)
legend("topright", inset = 0, horiz = TRUE,
       legend = c(TeX(r"($\sigma_t$)"), TeX(r"($\hat{\sigma_t}$)")),
       col = c("black", "red"), lty = 1, cex = 0.8)

plot(tail(dates, test_len), sqrt(tail(data$rv_d, test_len)),
     type = "l", ylim = c(0, 20), col = "black",
     main = "Log - TGARCHX", xlab = "Date", ylab = TeX(r"($\sigma_t$)")
     )
lines(tail(dates, test_len), sqrt(exp(preds_lgarchx)), col = "red"
)
legend("topright", inset = 0, horiz = TRUE, 
       legend = c(TeX(r"($\sigma_t$)"), TeX(r"($\hat{\sigma_t}$)")),
       col = c("black", "red"), lty = 1, cex = 0.8
       )
# dev.off()

setEPS()
postscript("./images/superiors.eps")
par(mfrow = c(3, 1), mar = c(2, 5, 2, 2)
)
plot(tail(dates, test_len), sqrt(tail(data$rv_d, test_len)),
     type = "l", ylim = c(0,22),
     col = "black",
     main = "VS-LTGARCHX (LASSO)", xlab = "Date", ylab = TeX(r"($\sigma_t$)")
     )
lines(tail(dates, test_len), sqrt(exp(preds_lasso)), col = "red"
)
legend("topright", inset = 0, horiz = TRUE,
       legend = c(TeX(r'($\sigma_t$)'), TeX(r'($\hat{\sigma_t}$)')),
       col = c("black", "red"), lty = 1, cex = 0.8
       )
plot(tail(dates, test_len), sqrt(tail(data$rv_d, test_len)),
    type = 'l', ylim = c(0, 22),
     col = 'black',
     main = "VS-LTGARCHX (ABESS)", xlab = "Date", ylab = TeX(r'($\sigma_t$)')
     )
lines(tail(dates,test_len), sqrt(exp(preds_abess)), col = "red"
)
legend("topright", inset = 0, horiz = TRUE,
       legend = c(TeX(r'($\sigma_t$)'), TeX(r'($\hat{\sigma_t}$)')),
       col = c("black", "red"), lty = 1, cex = 0.8
       )
plot(tail(dates, test_len), sqrt(tail(data$rv_d, test_len)), 
     type = 'l', ylim = c(0,22), col = 'black',
     main = "VS-LTGARCHX (Boruta)", xlab = "Date", ylab = TeX(r'($\sigma_t$)')
     )
lines(tail(dates, test_len), sqrt(exp(preds_fboruta)), col = "red"
)
legend("topright", inset = 0.01, horiz = TRUE,
       legend=c(TeX(r'($\sigma_t$)'), TeX(r'($\hat{\sigma_t}$)')),
       col = c("black", "red"), lty = 1, cex = 0.8
       )
# dev.off()
print(end_time - start_time)
save.image(file = "Log_GARCH.Rdata")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
