#Bayesian Optimization from scratch
# starting with https://rpubs.com/Argaadya/bayesian-optimization
library(tidymodels)
library(scales)
library(pso)
library(tidyverse)
install.packages("rBayesianOptimization")
library(rBayesianOptimization)
library(lubridate)
install.packages("GPfit")
library(GPfit)

# For Machine Learning
library(tidytext)
library(keras)
library(RVerbalExpressions)
library(textclean)


# Toy function

f <- function(x) {
  y <- (2 * x - 10)^2 * sin(3.14*2 * x - 4)
  return(y)
}

x <- c(0, 1/3, 1/2, 2/3, 1)
x<- c(0, 0.05, 1/3, 1/2, 2/3, 0.88, 1)
plot(x, f(x))

#Fit GP to f(x) for 0 <= x <= 1
eval <- data.frame(x = x, y = f(x)) %>% as.matrix()
fit <- GP_fit(X = eval[ , "x"],
              Y = eval[ , "y"],
              corr = list(type = "exponential", power = 1.95))

# calculate expected values for each possible x, as well as their uncertainty
x_new <-seq(0, 1, length.out=100)
pred <- predict.GP(fit, xnew = data.frame(x = x_new))
mu <- pred$Y_hat
sigma <- sqrt(pred$MSE)

ggplot(as.data.frame(eval)) +
  geom_line(data = data.frame(x = x_new, y = f(x_new)),
            aes(x = x, y = y), color = "black", linetype = "solid")+
  theme_minimal() +
  labs(title = "F(x) for one variable",
       y = "f(x)")

#Visualize the result
ggplot(as.data.frame(eval))+
  geom_line(data = data.frame(x = x_new, y = mu),
            aes(x = x, y = y), color = "red", linetype = "dashed") +
  geom_line(data = data.frame(x = x_new, y = f(x_new)),
            aes(x = x, y = y), color = "black", linetype = "solid")+
  geom_ribbon(data = data.frame(x = x_new, y_up = mu + sigma, y_low = mu - sigma), 
              aes(x = x_new, ymax = y_up, ymin = y_low), fill = "skyblue", alpha = 0.5) +
  geom_point(aes(x,y), size = 2)+
  theme_minimal() +
  labs(title = "Gaussian Process Posterior of f(x)",
       subtitle = "Blue area indicate the credible intervals",
       y = "f(x)")

# Testing gpfit in two dimension
computer_simulator <- function(x) {
  x1 = 4 * x[, 1] - 2
  x2 = 4 * x[, 2] - 2
  t1 = 1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 
                                6 * x1 *x2 + 3 * x2^2)
  t2 = 30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 
                                     36 * x1 * x2 + 27 * x2^2)
  y = t1 * t2
  return(y)
}
n = 30
d = 2
set.seed(100)
library(lhs)
x = maximinLHS(n, d)
y = computer_simulator(x)
GPmodel = GP_fit(x, y)
print(GPmodel)

Model_pred <- predict(GPmodel, maximinLHS(50, 2))
plot(GPmodel, range = c(0, 1), resolution = 50, surf_check = T, response = F, contour = F)




x2 <- rbind(x, c(0,1))
y2 <- computer_simulator(x2)
GPmodel2 = GP_fit(x2, y2)
plot(GPmodel2, range = c(0, 1), resolution = 50, surf_check = T, response = T, contour = F)




#Modeling on my data.....
predict.exp.res <- readRDS(file = "C:/Users/asolivai/Desktop/R_Files/AutoClustR/New_Predict_exp_params.rds")
zeisel <- predict.exp.res$index.frame$Zeisel
z.res <- zeisel[c("K.Param","Resolution", "Sil")]


z.start <- z.res[sample(1:24, 12, replace = F)*23 + sample(1:24, 12, replace = F) ,]
#Scale to unit cube
z.start[1] <- z.start[1] / 256
z.start[2] <- z.start[2] / 2.56


#Check if scaling the response changes expected improvement
z.def <- z.start
predictors <- z.def[1:2]
response <- z.def[[3]] # Y/response has to be a numeric
z.def.model <- GP_fit(X = predictors, Y = response)
plot(z.def.model, range = c(0, 1), resolution = 50, surf_check = T, response = F)

z.scaled <- z.start
z.scaled[3] <- z.scaled[3] * 10^1.55
predictors <- z.scaled[1:2]
scaled.response <- z.scaled[[3]] # Y/response has to be a numeric
z.scaled.model <- GP_fit(X = predictors, Y = scaled.response)
plot(z.scaled.model, range = c(0, 1), resolution = 50, surf_check = T, response = F)

# construct space on which to explore gaussian process
grid.space <- expand.grid(x = seq(0, 1, 0.02),
                          y = seq(0, 1, 0.02))

predicted.return <- predict(z.def.model, xnew = grid.space, each = 50)$complete_data
mu <- predicted.return[ , 3]
sigma <- predicted.return[ , 4]
y.max <- max(response)

expected_improvement <- map2_dbl(
  mu, sigma,
  function(m, s) {
    if(s == 0) {return(0)}
    gamma <- (m - y.max) / s
    phi <- pnorm(gamma)
    return(s * (gamma * phi + dnorm(gamma)))
  }
)

# same thing but for the scaled responses
scaled.return <- predict(z.scaled.model, xnew = grid.space, each = 50)$complete_data
mu <- scaled.return[ , 3]
sigma <- scaled.return[ , 4]
y.max <- max(response)

scaled_improvement <- map2_dbl(
  mu, sigma,
  function(m, s) {
    if(s == 0) {return(0)}
    gamma <- (m - y.max) / s
    phi <- pnorm(gamma)
    return(s * (gamma * phi + dnorm(gamma)))
  }
)

which.max(scaled_improvement)

#alright, cool, scaling doesn't really matter

#if you take the root of the means squared error, it might not matter too much

sqrt(0.004)
sqrt(3)


z.start[51, ] <- c(0.01, 0.8, 0.7)

#So the actual workflow is going to look like..........
#Get initial parameters (as [0, 1])
params <- maximinLHS(16, 2)
#Vector of space limits
scaled.params <- c(240, 2.4) * t(params)
tab.params <- cbind(round(scaled.params[1,]), scaled.params[2,])

#now that we have params as a latin square, need to grab SIl index values for each param pair



#%%%%%%%%%%%%%%%%%%% Creating the acquistion function   %%%%%%%%%%%%%%%%%%%%%%%%

#idk why the tutorial defines "y_best" as the minimum, but that's whatever
y_best <- min(eval[,2])

#influences expoloration. Higher epsilon = more exploration
epsilon <- 0.01

#function for calculating expected improvement
ei_calc <- function(m_s)  {        # mu and sigma, respectively
  m <- m_s[1]
  s <- m_s[2]
  if (s == 0) {return(0)}            # if deviation is 0, the point has already been explored
  Z <- (m - y_best - epsilon)/s
  expected_imp <- (m - y_best - epsilon) * pnorm(Z) + s * dnorm(Z)
  return(expected_imp)
}

m_s <- data.frame(mu = mu,
                  sigma = sigma)
expected_improvement <- apply(m_s, MARGIN = 1, ei_calc)
max(expected_improvement)


exp_imp <- data.frame(x = x_new,
                      y = expected_improvement)



expected_improvement <- map2_dbl(
  mu, sigma,
  function(m, s) {
    if(s == 0) {return(0)}
    gamma <- (m - y.max) / s
    phi <- pnorm(gamma)
    return(s * (gamma * phi + dnorm(gamma)))
  }
)

#Visualize expected improvement
ggplot(exp_imp, aes(x, y))+
  geom_line()+
  geom_ribbon(aes(ymin = 0, ymax = y), fill = "skyblue", alpha = 0.5, color = "white")+ 
  geom_vline(xintercept = exp_best$x, linetype = "dashed", color = "red")+
  geom_point(data = exp_best, size = 2)+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  scale_x_continuous(breaks = c(seq(0,1,0.25), round(exp_best$x,2)))+
  labs(title = "Expected Improvement",
       subtitle = "x with the highest expected improvement will be evaluated",
       y = "Expected Improvement")

#Get the next value to test
exp_best <- exp_imp %>% filter(y == max(y))

exp_best$x

#Fit GP to f(x) for 0 <= x <= 1
eval.2 <- rbind(eval, c(exp_best$x, f(exp_best$x))) %>% as.matrix()
fit <- GP_fit(X = eval.2[2:6 , "x"],
              Y = eval.2[2:6 , "y"],
              corr = list(type = "exponential", power = 1.95))

# calculate expected values for each possible x, as well as their uncertainty
x_new <-seq(0, 1, length.out=100)
pred <- predict.GP(fit, xnew = data.frame(x = x_new))
mu <- pred$Y_hat
sigma <- sqrt(pred$MSE)

#Visualize the result
ggplot(as.data.frame(eval.2))+
  geom_line(data = data.frame(x = x_new, y = mu),
            aes(x = x, y = y), color = "red", linetype = "dashed")+
  geom_ribbon(data = data.frame(x = x_new, y_up = mu + sigma, y_low = mu - sigma), 
              aes(x = x_new, ymax = y_up, ymin = y_low), fill = "skyblue", alpha = 0.5) +
  geom_point(aes(x,y), size = 2)+
  theme_minimal() +
  labs(title = "Gaussian Process Posterior of f(x)",
       subtitle = "Blue area indicate the credible intervals",
       y = "f(x)")





















