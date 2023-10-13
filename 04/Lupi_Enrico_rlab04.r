# Enrico Lupi, 2090596, 14 May 2023


# LABORATORY 4

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(knitr)
library(plotly)

# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- FALSE


# function to compute quantiles
quant <- function(val, p, post, dp) {
    q <- 0
    sum <- 0
    for (i in seq(1, length(p))) {
        sum <- sum + dp * post[i]
        # take the first value when the integral surpasses the
        # desired value
        if (sum >= val) {
            q <- p[i]
            break
        }
    }
    return(q)
}
# function to compute credibility intervals
cred_interval <- function(val, p, post, dp) {
    limits <- c(0, 0)
    v <- (1 - val) / 2
    limits[1] <- quant(v,     p, post, dp)
    limits[2] <- quant(1 - v, p, post, dp)
    return(limits)
}


# EXERCISE 1

n <- 10
claims <- c(5, 8, 4, 6, 11, 6, 6, 5, 6, 4)
delta <- 0.01
mu <- seq(0, 20, delta) # we will only consider mu values up till 20,
                        # even though formally they could go up to +Inf,
                        # as the probability gets too small for those
                        # values to be interesting

# PART 1 : suppose to use a prior uniform distribution for mu
#          - find the posterior distribution for mu and compute the posterior
#            mean, median and variance
#          - plot the posterior distribution and the 95% credibility interval

# We can use an improper prior g(mu) = 1 for mu > 0
# The posterior thus becomes proportional to the likelihood, a Gamma
# distribution with:
alpha_unif <- sum(claims) + 1
lambda <- n
post_unif <- dgamma(mu, alpha_unif, lambda)
post_unif <- post_unif / (delta * sum(post_unif))

mean_unif <- delta * sum(mu * post_unif)
var_unif  <- delta * sum((mu - mean_unif)^2 * post_unif)
med_unif  <- quant(0.5, mu, post_unif, delta)

# 95% C. I.
lim_unif <- cred_interval(0.95, mu, post_unif, delta)

val_unif <- c(mean_unif, var_unif, med_unif, lim_unif) |> round(digits = 3)

# plot results
df_unif <- data.frame(mu, post_unif)
g_unif  <-  ggplot() +
            geom_line(data = df_unif, mapping = aes(x = mu, y = post_unif,
                      color = "PDF"), linewidth = 1, linetype = 1) +
            geom_vline(aes(xintercept = mean_unif,
                           color = "Mean"), linetype = 2) +
            geom_area(data = subset(df_unif, mu >= lim_unif[1] &
                                             mu <= lim_unif[2]),
                      mapping = aes(x = subset(mu, mu >= lim_unif[1] &
                                                   mu <= lim_unif[2]),
                                    y = post_unif,
                                    fill = "95% C.I."),
                      alpha = 0.7) +
            scale_color_manual(values = c("black",
                                          "blue")) +
            scale_fill_manual(values = "lightblue") +
            labs(x = "\u00b5", y = "P(\u00b5 | y, n, M)",
                 title = "Ex. 1 : Posterior from Uniform Prior",
                 color = "", fill = "")

if (save_images) {
    ggsave("Ex4-1_unif.png")
}


# PART 2 : suppose to use a Jeffrey's prior for mu
#          - find the posterior distribution for mu and compute the posterior
#            mean, median and variance
#          - plot the posterior distribution and the 95% credibility interval

# We can use an improper prior g(mu) = 1/sqrt(mu) for mu > 0
# The posterior thus becomes a Gamma distribution with:
alpha_jeff <- sum(claims) + 0.5
lambda <- n
post_jeff <- dgamma(mu, alpha_jeff, lambda)
post_jeff <- post_jeff / (delta * sum(post_jeff))

mean_jeff <- delta * sum(mu * post_jeff)
var_jeff  <- delta * sum((mu - mean_jeff)^2 * post_jeff)
med_jeff  <- quant(0.5, mu, post_jeff, delta)

# 95% C. I.
lim_jeff <- cred_interval(0.95, mu, post_jeff, delta)

val_jeff <- c(mean_jeff, var_jeff, med_jeff, lim_jeff) |> round(digits = 3)

# plot results
df_jeff <- data.frame(mu, post_jeff)
g_jeff <-   ggplot() +
            geom_line(data = df_jeff, mapping = aes(x = mu, y = post_jeff,
                      color = "PDF"), linewidth = 1, linetype = 1) +
            geom_vline(aes(xintercept = mean_jeff,
                              color = "Mean"), linetype = 2) +
            geom_area(data = subset(df_jeff, mu >= lim_jeff[1] &
                                             mu <= lim_jeff[2]),
                      mapping = aes(x = subset(mu, mu >= lim_jeff[1] &
                                                   mu <= lim_jeff[2]),
                                    y = post_jeff,
                                    fill = "95% C.I."),
                      alpha = 0.7) +
            scale_color_manual(values = c("black",
                                          "red")) +
            scale_fill_manual(values = "#f45c5c") +
            labs(x = "\u00b5", y = "P(\u00b5 | y, n, M)",
                 title = "Ex. 1 : Posterior from Jeffrey's Prior",
                 color = "", fill = "")

if (save_images) {
    ggsave("Ex4-1_jeff.png")
}


# PART 3 : evaluate a 95% credibility interval for the results obtained with
#          both priors. Compare the result with that obtained using a normal
#          approximation for the posterior distribution

# normal approximation of posterior from uniform prior
norm_unif <- dnorm(mu, mean_unif, sqrt(var_unif))
norm_unif <- norm_unif / (delta * sum(norm_unif))

med_normu <- quant(0.5, mu, norm_unif, delta)
lim_normu <- cred_interval(0.95, mu, norm_unif, delta)

val_normu <- c(mean_unif, var_unif, med_normu, lim_normu) |> round(digits = 3)

# normal approximation of posterior from Jeffrey's prior
norm_jeff <- dnorm(mu, mean_jeff, sqrt(var_jeff))
norm_jeff <- norm_jeff / (delta * sum(norm_jeff))

med_normj <- quant(0.5, mu, norm_jeff, delta)
lim_normj <- cred_interval(0.95, mu, norm_jeff, delta)

val_normj <- c(mean_jeff, var_jeff, med_normj, lim_normj) |> round(digits = 3)

# plot results
g_normunif <- ggplot() +
            geom_line(aes(x = mu, y = post_unif, color = "Posterior"),
                      linewidth = 0.8, linetype = 1, alpha = 0.8) +
            geom_line(aes(x = mu, y = norm_unif,
                          color = "Normal Approximation"),
                      linewidth = 0.8, linetype = 1, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_unif[1], color = "Posterior"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_unif[2], color = "Posterior"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_normu[1],
                           color = "Normal Approximation"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_normu[2],
                           color = "Normal Approximation"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            scale_color_manual(values = c("black",
                                          "blue")) +
            labs(x = "\u00b5", y = "P(\u00b5 | y, n, M)",
                 title = "Ex. 1 : Posterior from Uniform Prior",
                 color = "")

g_normjeff <- ggplot() +
            geom_line(aes(x = mu, y = post_jeff, color = "Posterior"),
                      linewidth = 0.8, linetype = 1, alpha = 0.8) +
            geom_line(aes(x = mu, y = norm_jeff,
                          color = "Normal Approximation"),
                      linewidth = 0.8, linetype = 1, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_jeff[1], color = "Posterior"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_jeff[2], color = "Posterior"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_normj[1],
                           color = "Normal Approximation"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = lim_normj[2],
                           color = "Normal Approximation"),
                       linewidth = 0.5, linetype = 2, alpha = 0.8) +
            scale_color_manual(values = c("black",
                                          "red")) +
            labs(x = "\u00b5", y = "P(\u00b5 | y, n, M)",
                 title = "Ex. 1 : Posterior from Jeffrey's Prior",
                 color = "")

g_norm <- ggarrange(g_normunif, g_normjeff,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)

if (save_images) {
    ggsave("Ex4-1_normCI.png")
}

# dataframe containing statistics estimators for the two posteriors
df_ex1 <- data.frame(val_unif, val_normu, val_jeff, val_normj)
rownames(df_ex1) <- c("Mean", "Var", "Median",
                      "95% C.I. Left", "95% C.I. Right")
colnames(df_ex1) <- c("Uniform P.", "Nor. App. Uniform",
                      "Jeffrey's P.", "Nor. App. Jeffrey's")


# ---------------------------------------------------------------------- #


# EXERCISE 2

# In this case the null hypothesis is given by the standard test, which has a
# probability of failure of
p0 <- 0.15

# The new test is performed on
n2 <- 75


# PART 1 : what is the probability distribution of y, the number of times the
#          new method fails to detect the disease?

# It is a binomial distribution with size = n2 and probability p where we
# consider as successes the failures of the test


# PART 2 : on the n2 = 75 patients sample, the new method fails to detect the
#          disease in y = 6 cases. What is the frequentist estimator of the
#          failure probability of the new method?

y <- 6
p_freq <- y / n2

cat(c("Ex. 2.2:\nThe frequentist estimator of the failure probability",
      " of the new method is", p_freq, "\n"))


# PART 3 : setup a bayesian computation of the posterior probability, assuming
#          a beta distribution with mean value 0.15 and standard deviation
#          0.14. Plot the posterior distribution for p, and mark on the plot
#          the mean value and variance

# Probability values
p <- seq(0, 1, length.out = 201)
n_p <- length(p)
delta_p <- p[2] - p[1]

# get prio parameters
mu2 <- 0.15
var <- 0.14
alpha_prior <-      mu2  * (mu2 * (1 - mu2) / var - 1)
beta_prior  <- (1 - mu2) * (mu2 * (1 - mu2) / var - 1)

# compute posterior
alpha_post <- alpha_prior + y
beta_post <- beta_prior + n2 - y

post_beta <- dbeta(x = p, shape1 = alpha_post, shape2 = beta_post)
post_beta <- post_beta / (delta_p * sum(post_beta))

mean_beta <- delta_p * sum(p * post_beta)
var_beta  <- delta_p * sum((p - mean_beta)^2 * post_beta)

df <- data.frame(p, post_beta)
dx <- mean_beta + sqrt(var_beta)
sx <- mean_beta - sqrt(var_beta)
g_posterior <- ggplot() +
               geom_line(data = df, mapping = aes(x = p, y = post_beta,
                         color = "PDF"), linewidth = 1, linetype = 1) +
               geom_vline(aes(xintercept = mean_beta,
                              color = "Mean"), linetype = 2) +
               geom_vline(aes(xintercept = p0,
                              color = "p0"), linetype = 2) +
               geom_area(data = subset(df, p > sx & p < dx),
                         mapping = aes(x = subset(p, p > sx & p < dx),
                                       y = post_beta,
                         fill = "Mean \u00b1 \u03c3"), alpha = 0.7) +
               geom_area(data = subset(df, p >= p0),
                         mapping = aes(x = subset(p, p >= p0),
                                       y = post_beta,
                         fill = "Null H probability"), alpha = 0.7) +
               scale_color_manual(values = c("blue",
                                             "red",
                                             "black")) +
               scale_fill_manual(values = c("lightblue",
                                            "#f45c5c")) +
               labs(x = "p", y = "P(p | y, n, M)",
                    title = "Ex. 2 : Posterior from Beta Prior",
                    color = "", fill = "")

if (save_images) {
    ggsave("Ex4-2_posterior.png")
}


# PART 4 : Perform a test of hypothesis assuming that if the probability of
#          failing to the detect the desease in ill patients is greater or
#          equal than 15%, the new test is no better that the traditional
#          method. Test the sample at a 5% level of significance in the
#          Bayesian way.

# We want to test the null hypothesis H0 p >= p0 = 0.15 against the alternative
# hypothesis H1 p < p0 = 0.15. In order to do so, we evaluate the probability
# of the null hypothesis by integrating the posterior in the required region
# and check wether the value we obtain is lesser than the significance
# or not, and we respectively can or cannot reject the null hypothesis

# integrate right side of posterior
i_max <- which(p == p0)
p_nullh <- 0
for (i in seq(n_p, i_max, -1)){
    p_nullh <- p_nullh + delta_p * post_beta[i]
}

significance <- 0.05
cat(c("Ex 4.2.4:\nThe probability of the null hypotheisis is",
      p_nullh * 100, "%"))
if (p_nullh > significance) {
    cat(c("\nThat is more than the significance", significance * 100, "%",
          "\nso we cannot reject the null hypothesis\n"))
} else {
    cat(c("\nThat is less than the significance", significance * 100, "%",
          "\nso we can reject the null hypothesis\n"))
}


# PART 5 : Perform the same hypothesis test in the classical frequentist way.

ys <- seq(0, 10)
# calculate the p_values for every y
# in this case the p-value is the cumulative distribution
# as we need to integrate from the left side
p_values <- cumsum(dbinom(ys, n2, p0))
g_pdfbar <- ggplot(data.frame(ys, p_values)) +
            geom_bar(aes(x = ys, y = p_values), color = "black",
                     stat = "identity", fill = "lightblue", alpha = 0.5) +
            geom_hline(aes(yintercept = significance,
                            color = "\u03b1 = 5%"), linetype = 2) +
            scale_color_manual(values = c("red")) +
            labs(x = "# test failures", y = "F(y)",
                 title = "Ex. 2 : Frequentist Test", color = "")

if (save_images) {
    ggsave("Ex4-2_freqTest.png")
}

p_val <- p_values[which(ys == y)]
cat(c("Ex 4.2.5:\nThe p-value is",
      p_val * 100, "%"))
if (p_val > significance) {
    cat(c("\nThat is more than the significance", significance * 100, "%",
          "\nso we cannot reject the null hypothesis\n"))
} else {
    cat(c("\nThat is less than the significance", significance * 100, "%",
          "\nso we can reject the null hypothesis\n"))
}


# ---------------------------------------------------------------------- #


# EXERCISE 3

# A lighthouse is located at a position alpha along the shore and at a distance
# beta out at sea. It emits a series of short highly collimated flashes at
# random intervals and at random angles. we detect the pulses on the coast
# using photo-detectors; they record only the position xk of the flash arrival
# on the coast, but not the angle of emission. N flashes have been recorded at
# positions {xk} We want to estimate the position of the lighthouse

# The distribution of the azimuthal angle of emission is uniform in
# [-pi/2; pi/2]. By a change of variable we can derive the distribution of the
# positions xk, a Cauchy distribution.
# We can thus assume true values for alpha and beta and use them to generate
# a sample D of the positions

alpha_true <- 1 #km
beta_true <- 1.5 #km

set.seed(4888)
x <- rcauchy(100, alpha_true, beta_true)

# We assume a uniform prior for both paramaters, so 1/4 for alpha in [-2; 2]km
# and 1/3 for beta in [0.5; 3.5]km. Note that we do not actually care for the
# value of this prior, but only that it limits the parameters to a specific 2D
# region we will investigate in, as all constant factor will be taken care of
# in the end after the normalisation.

# The likelihood of the entire dataset is simply given by the product of the
# single Cauchy probabilities, as the measures are independent one from the
# other. The posterior is proportional to the likelihood since, as we discussed
# before, the prior is just a constant factor. It is easier to work with the
# Log of the posterior: we can thus define (ignoring constant factors)

p_log <- function(a, b, data) {
    logL <- -1 * length(data) * log(b)
    for (x in data) {
        logL <- logL - log(1 + ((x - a) / b)^2)
    }
    return(logL)
}

# sample grid for parameters A and B
alphalim <- c(-2, 2)
betalim <- c(0.5, 3.5)
n3 <- 200
uni_grid <- seq(from = 1 / (2 * n3), to = 1 - 1 / (2 * n3), by = 1 / n3)
delta_alpha <- diff(alphalim) / n3
delta_beta  <- diff(betalim) / n3
alpha <- alphalim[1] + diff(alphalim) * uni_grid
beta  <- betalim[1]  + diff(betalim)  * uni_grid

# Now we compute the value of the Log posterior given the dataset for every
# combination of alpha and beta

logp <- outer(alpha, beta,
              Vectorize(function(alpha, beta) p_log(alpha, beta, x)))
logp <- logp - max(logp)

max_ind <- which(logp == max(logp), arr.ind = TRUE)
alpha_max <- alpha[max_ind[1]]
beta_max  <- beta[max_ind[2]]

# Now we get the actual posterior by taking the exponential...
ab_post <- exp(logp)
# and normalising with a simple 2D Riemann sum
ab_post <- ab_post / (delta_alpha * delta_beta * sum(ab_post))

# plot results

# 3D plot of posterior
fig <- plot_ly(x = alpha, y = beta, z = t(ab_post)) |> add_surface()

g_alpha <-  ggplot() +
            geom_line(aes(x = alpha, y = ab_post[, max_ind[2]],
                          color = "PDF"),
                      linewidth = 1, linetype = 1, alpha = 1) +
            geom_vline(aes(xintercept = alpha_max, color = "Max \u03b1"),
                       linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = alpha_true, color = "True \u03b1"),
                       linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = mean(x), color = "Mean value D"),
                       linetype = 2) +
            scale_color_manual(values = c("blue",
                                          "#00a100",
                                          "black",
                                          "red")) +
            labs(x = "\u03b1 [km]", y = "P(\u03b1, \u03b2 | D)",
                 title = paste("Ex. 3 : Cross Section of Posterior",
                               "Distribution for \u03b2 = \u03b2_max"),
                 color = "")

g_beta  <-  ggplot() +
            geom_line(aes(x = beta, y = as.numeric(ab_post[max_ind[1], ]),
                          color = "PDF"),
                      linewidth = 1, linetype = 1, alpha = 1) +
            geom_vline(aes(xintercept = beta_max, color = "Max \u03b2"),
                       linetype = 2, alpha = 0.8) +
            geom_vline(aes(xintercept = beta_true, color = "True \u03b2"),
                       linetype = 2, alpha = 0.8) +
            scale_color_manual(values = c("blue",
                                          "black",
                                          "red")) +
            labs(x = "\u03b2 [km]", y = "P(\u03b1, \u03b2 | D)",
                 title = paste("Ex. 3 : Cross Section of Posterior",
                               "Distribution for \u03b1 = \u03b1_max"),
                 color = "")

g_lighthouse <- ggarrange(g_alpha, g_beta,
                          labels = c("A", "B"),
                          ncol = 1, nrow = 2)

if (save_images) {
    ggsave("Ex4-3_lighthouseCross.png")
}

dat <- expand_grid(x = alpha, y = beta)
dat <- mutate(dat, z = as.vector(t(ab_post)))
g_contour <- ggplot(dat, aes(x, y)) +
             geom_contour_filled(aes(z = z), alpha = 0.9) +
             theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
             labs(x = "\u03b1 [km]", y = "\u03b2 [km]",
                  title = "Ex. 3 : Posterior Distribution Contour Plot",
                  color = "") +
             geom_hline(aes(yintercept = beta_true,  color = "True values"),
                        linetype = 2) +
             geom_vline(aes(xintercept = alpha_true, color = "True values"),
                        linetype = 2) +
             scale_color_manual(values = c("black",
                                           "black")) +
             xlim(alphalim[1], alphalim[2]) +
             ylim(betalim[1], betalim[2])

if (save_images) {
    ggsave("Ex4-3_contour.png")
}


# ---------------------------------------------------------------------- #


# EXERCISE 4

# set model parameters
X0 <- 0                     # Peak position
W <- 1                      # signal resolution
A_true <- 2                 # signal amplitude
B_true <- 1                 # background amplitude
Dt <- 5                     # exposure time

# sample grid for parameters A and B
alim <- c(0, 4)
blim <- c(0.5, 1.5)
n4 <- 100
uni_grid <- seq(from = 1 / (2 * n4), to = 1 - 1 / (2 * n4), by = 1 / n4)
delta_a <- diff(alim) / n4
delta_b <- diff(blim) / n4
a <- alim[1] + diff(alim) * uni_grid
b <- blim[1] + diff(blim) * uni_grid

# generative model
signal <- function(x, a, b, x0, w, dt) {
    dt * (a * exp(-(x - x0)^2 / (2 * w^2)) + b)
}

# log posterior
p_log <- function(d, x, a, b, x0, w, dt) {
    if (a < 0 || b < 0) {
        return(-Inf)
    }
    sum(dpois(d, lambda = signal(x, a, b, x0, w, dt), log = TRUE))
}

# function to calcolate 2D posterior
computePost <- function(ddat, xdat, a, b, x0, w, dt) {
    # compute log posterior (unnormalised)
    z <- matrix(data = NA, nrow = length(a), ncol = length(b))
    for (i in seq(1, length(a))) {
        for (j in seq(1, length(b))) {
            z[i, j] <- p_log(ddat, xdat, a[i], b[j], x0, w, dt)
        }
    }
    # set maximum to zero
    z <- z - max(z)
    # compute posterior
    z <- exp(z)

    return(z)
}

# marginalise posterior along one dimension
marginalise <- function(z, delta, dim = 1) {
    p <- apply(z, dim, sum)
    p <- p / (delta * sum(p))
    return(p)
}

# compute normalized conditional posteriors P(a|b,D) and P(b|a,D)
# using true values of conditioned parameters
cond_prob <- function(par, ddat, xdat, a, b, x0, w, dt, delta) {
    p <- exp(Vectorize(p_log, par)(ddat, xdat, a, b, x0, w, dt))
    p <- p / (delta * sum(p))
}

# compute covariance
covariance <- function(a, b, mean_a, mean_b, z) {
    cov_ab <- 0
    for (i in seq(1, length(a))) {
        for (j in seq(1, length(b))) {
            cov_ab <- cov_ab + (a[i] - mean_a) * (b[j] - mean_b) * z[i, j]
        }
    }
    cov_ab <- cov_ab / sum(z)
    return(cov_ab)
}

# produce contour plot given the parameters and posterior
contourPlot <- function(a, b, z, a_true, b_true) {

    print(c(a_true, b_true))
    # plot results
    dat <- expand_grid(x = a, y = b)
    dat <- mutate(dat, z = as.vector(t(z)))
    g_contour <- ggplot(dat, aes(x, y)) +
                 geom_contour_filled(aes(z = z), alpha = 0.9) +
                 theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
                 labs(x = "A", y = "B", color = "",
                      title = "") +
                 geom_vline(aes(xintercept = a_true, color = "True values"),
                            linetype = 2) +
                 geom_hline(aes(yintercept = b_true, color = "True values"),
                            linetype = 2) +
                 scale_color_manual(values = c("black",
                                               "black")) +
                 xlim(alim[1], alim[2]) +
                 ylim(blim[1], blim[2])

    return(g_contour)
}

# produce plot showing marginalised and conditional posterior
postPlot <- function(a, b, p_a_D, p_b_D, p_a_Db, p_b_Da, a_true, b_true) {
    g_a <-  ggplot() +
            geom_line(aes(x = a, y = p_a_D, color = "Marginalised"),
                      linewidth = 1, linetype = 1) +
            geom_line(aes(x = a, y = p_a_Db, color = "Conditional"),
                      linewidth = 1, linetype = 5) +
            geom_vline(aes(xintercept = a_true, color = "True Value"),
                       linewidth = 0.8, linetype = 1) +
            scale_color_manual(values = c("blue",
                                          "red",
                                          "grey")) +
            labs(x = "amplitude A", y = "P(A | D)  and  P(A | B,D)",
                 title = "", color = "")

    g_b <-  ggplot() +
            geom_line(aes(x = b, y = p_b_D, color = "Marginalised"),
                      linewidth = 1, linetype = 1) +
            geom_line(aes(x = b, y = p_b_Da, color = "Conditional"),
                      linewidth = 1, linetype = 5) +
            geom_vline(aes(xintercept = b_true, color = "True Value"),
                       linewidth = 0.8, linetype = 1) +
            scale_color_manual(values = c("blue",
                                          "red",
                                          "grey")) +
            labs(x = "background B", y = "P(B | D)  and  P(B | A,D)",
                 title = "", color = "")

    g_ab <- ggarrange(g_a, g_b,
                      labels = c("", ""),
                      ncol = 2, nrow = 1,
                      common.legend = TRUE)
    return(g_ab)
}


# PART 1 : vary the sampling resolution of used to generate the data and
#          check the effect on the result

ws <- c(0.1, 0.25, 1, 2, 3) # signal resolutions

# objects to hold loop results
labels_w <- vector("list", length(ws))
ls_contw <- vector("list", length(ws))
ls_postw <- vector("list", length(ws))
cnames <- c("w", "Mean A", "Std A", "Mean B", "Std B", "Cov AB", "Rho AB")
stat_w <- data.frame(matrix(nrow = length(ws), ncol = length(cnames)))
colnames(stat_w) <- cnames

# Loop over different resolutions w
for (i in seq(1, length(ws))) {
    labels_w[[i]] <- paste("w = ", ws[i])
    # generate data
    xdat <- seq(from = -7 * W, to = 7 * W, by = 0.5 * ws[i])
    true_sig <- signal(xdat, A_true, B_true, X0, ws[i], Dt)
    set.seed(4888)
    ddat <- rpois(length(true_sig), true_sig)

    # compute posterior and plot contour plot
    z <- computePost(ddat, xdat, a, b, X0, ws[i], Dt)

    ls_contw[[i]] <- contourPlot(a, b, z, A_true, B_true)

    # compute marginalised posteriors and plot them
    p_a_D <- marginalise(z, delta_a, 1)
    p_b_D <- marginalise(z, delta_b, 2)
    p_a_Db <- cond_prob("a", ddat, xdat, a, B_true, X0, ws[i], Dt, delta_a)
    p_b_Da <- cond_prob("b", ddat, xdat, A_true, b, X0, ws[i], Dt, delta_b)

    ls_postw[[i]] <- postPlot(a, b, p_a_D, p_b_D,
                              p_a_Db, p_b_Da, A_true, B_true)

    # Compute mean, standard deviation, covariance, correlation of A and B
    mean_a <- delta_a * sum(a * p_a_D)
    mean_b <- delta_b * sum(b * p_b_D)
    sd_a <- sqrt(delta_a * sum((a - mean_a)^2 * p_a_D))
    sd_b <- sqrt(delta_b * sum((b - mean_b)^2 * p_b_D))
    cov_ab <- covariance(a, b, mean_a, mean_b, z)
    rho_ab <- cov_ab / (sd_a * sd_b)

    stat_w[i, 1] <- ws[i]  |> round(digits = 3)
    stat_w[i, 2] <- mean_a |> round(digits = 3)
    stat_w[i, 4] <- sd_a   |> round(digits = 3)
    stat_w[i, 3] <- mean_b |> round(digits = 3)
    stat_w[i, 5] <- sd_b   |> round(digits = 3)
    stat_w[i, 6] <- cov_ab |> round(digits = 3)
    stat_w[i, 7] <- rho_ab |> round(digits = 3)
}

g_contw <- ggarrange(plotlist = ls_contw,
                     labels = labels_w,
                     ncol = 2, nrow = 3,
                     common.legend = TRUE)
g_contw <- annotate_figure(g_contw,
                top = text_grob("Ex. 4 : Posterior Distribution Contour Plot"))
if (save_images) {
    ggsave("Ex4-4_contourW.png")
}

g_postw1 <- ggarrange(plotlist = ls_postw[1:3],
                      labels = labels_w[1:3],
                      ncol = 1, nrow = 3,
                      common.legend = TRUE)
g_postw1 <- annotate_figure(g_postw1,
                    top = text_grob("Ex. 4 : Posteriors varying Resolution w"))
if (save_images) {
    ggsave("Ex4-4_posteriorsW1.png")
}

g_postw2 <- ggarrange(plotlist = ls_postw[4:5],
                      labels = labels_w[4:5],
                      ncol = 1, nrow = 2,
                      common.legend = TRUE)
g_postw2 <- annotate_figure(g_postw2,
                    top = text_grob("Ex. 4 : Posteriors varying Resolution w"))
if (save_images) {
    ggsave("Ex4-4_posteriorsW2.png")
}

# A greater resolution leads to a higher accuracy in determining the values
# of the parameters: in particular, B is highly sensitive to changes in w
# and its standard deviation changes of a order of magnitude across the
# investigated range. Rho also varies, as A and B appear to be more
# anticorrelated the smaller the resolution we use.


# PART 2 : change the ratio A=B used to simulate the data (keeping both
#          positive in accordance with the prior) and
#          check the effect on the result

As <- c(0.1, 0.5,  1,   2, 2.75, 3.5)
Bs <- c(1.4, 1.33, 1.1, 1, 0.8,  0.6)
ratio <- (As / Bs) |> round(digits = 3)

# objects to hold loop results
labels_r <- vector("list", length(ratio))
ls_contr <- vector("list", length(ratio))
ls_postr <- vector("list", length(ratio))
cnames <- c("Ratio", "True A", "Mean A", "Std A", "True B",
            "Mean B", "Std B", "Cov AB", "Rho AB")
stat_r <- data.frame(matrix(nrow = length(ratio), ncol = length(cnames)))
colnames(stat_r) <- cnames

# Loop over different A/B ratios
for (i in seq(1, length(ratio))) {
    # plot labels
    labels_r[[i]] <- paste("R = ", ratio[i])

    # generate data
    xdat <- seq(from = -7 * W, to = 7 * W, by = 0.5 * W)
    true_sig <- signal(xdat, As[i], Bs[i], X0, W, Dt)
    set.seed(4888)
    ddat <- rpois(length(true_sig), true_sig)

    # compute posterior and plot contour plot
    z <- computePost(ddat, xdat, a, b, X0, W, Dt)

    ls_contr[[i]] <- contourPlot(a, b, z, As[i], Bs[i])

    # compute marginalised posteriors and plot them
    p_a_D <- marginalise(z, delta_a, 1)
    p_b_D <- marginalise(z, delta_b, 2)
    p_a_Db <- cond_prob("a", ddat, xdat, a, Bs[i], X0, W, Dt, delta_a)
    p_b_Da <- cond_prob("b", ddat, xdat, As[i], b, X0, W, Dt, delta_b)

    ls_postr[[i]] <- postPlot(a, b, p_a_D, p_b_D,
                              p_a_Db, p_b_Da, As[i], Bs[i])

    # Compute mean, standard deviation, covariance, correlation of A and B
    mean_a <- delta_a * sum(a * p_a_D)
    mean_b <- delta_b * sum(b * p_b_D)
    sd_a <- sqrt(delta_a * sum((a - mean_a)^2 * p_a_D))
    sd_b <- sqrt(delta_b * sum((b - mean_b)^2 * p_b_D))
    cov_ab <- covariance(a, b, mean_a, mean_b, z)
    rho_ab <- cov_ab / (sd_a * sd_b)

    stat_r[i, 1] <- ratio[i]  |> round(digits = 3)
    stat_r[i, 2] <- As[i]     |> round(digits = 3)
    stat_r[i, 3] <- mean_a    |> round(digits = 3)
    stat_r[i, 4] <- sd_a      |> round(digits = 3)
    stat_r[i, 5] <- Bs[i]     |> round(digits = 3)
    stat_r[i, 6] <- mean_b    |> round(digits = 3)
    stat_r[i, 7] <- sd_b      |> round(digits = 3)
    stat_r[i, 8] <- cov_ab    |> round(digits = 3)
    stat_r[i, 9] <- rho_ab    |> round(digits = 3)
}


g_contr <- ggarrange(plotlist = ls_contr,
                     labels = labels_r,
                     ncol = 2, nrow = 3,
                     common.legend = TRUE)
g_contr <- annotate_figure(g_contr,
                top = text_grob("Ex. 4 : Posterior Distribution Contour Plot"))
if (save_images) {
    ggsave("Ex4-4_contourR.png")
}

g_postr1 <- ggarrange(plotlist = ls_postr[1:3],
                      labels = labels_r[1:3],
                      ncol = 1, nrow = 3,
                      common.legend = TRUE)
g_postr1 <- annotate_figure(g_postr1,
                    top = text_grob("Ex. 4 : Posteriors varying A/B Ratio"))
if (save_images) {
    ggsave("Ex4-4_posteriorsR1.png")
}

g_postr2 <- ggarrange(plotlist = ls_postr[4:6],
                      labels = labels_r[4:6],
                      ncol = 1, nrow = 3,
                      common.legend = TRUE)
g_postr2 <- annotate_figure(g_postr2,
                    top = text_grob("Ex. 4 : Posteriors varying A/B Ratio"))

if (save_images) {
    ggsave("Ex4-4_posteriorsR2.png")
}


# The ratio A/B does not seem to have a great effect on the parameter
# estimation, as both the standard deviations and the rho do not show any
# particular trend. The results for the lower and higher ratios, however,
# could be less precise as the true values were located at the edges of the
# parameter space, so further analysis on a bigger range of parameters
# may be necessary.


# ---------------------------------------------------------------------- #


# Save output to pdf file
if (save_pdf) {
    pdf("Lupi_Enrico_rlab04_Outputs.pdf")
    print(g_unif)
    print(g_jeff)
    print(g_norm)
    grid.newpage()
    grid.table(df_ex1)
    print(g_posterior)
    print(g_pdfbar)
    print(g_lighthouse)
    print(g_contour)
    print(g_contw)
    print(g_postw1)
    print(g_postw2)
    grid.newpage()
    grid.table(stat_w)
    print(g_contr)
    print(g_postr1)
    print(g_postr2)
    grid.newpage()
    grid.table(stat_r)
    dev.off()
}