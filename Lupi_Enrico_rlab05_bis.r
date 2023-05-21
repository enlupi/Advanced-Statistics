# Enrico Lupi, 2090596, 28 May 2023


# LABORATORY 4

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(knitr)
library(plotly)
library(rjags)
library(R2jags)


# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- FALSE


# compute quantiles
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
# compute credibility intervals
cred_interval <- function(val, p, post, dp) {
    limits <- c(0, 0)
    v <- (1 - val) / 2
    limits[1] <- quant(v,     p, post, dp)
    limits[2] <- quant(1 - v, p, post, dp)
    return(limits)
}

# plot pdf and credibility intervals
plotGraph <- function(x, y, mean, lim, color, range) {
    fill_color <- "lightblue"
    if (color == "red") {
        fill_color <- "#f45c5c"
    } else if (color == "#039303") {
        fill_color <- "#71de71"
    }
cat(c(mean, "\n"))
    df <- subset(data.frame(x, y), x >= range[1] & x <= range[2])
    g <- ggplot() +
        geom_line(data = df, mapping = aes(x = x, y = y, color = "PDF"),
                  linewidth = 1, linetype = 1) +
        geom_vline(aes(xintercept = mean, color = "Mean"), linetype = 2) +
        geom_area(data = subset(df, x >= lim[1] & x <= lim[2]),
                  mapping = aes(x = subset(x, x >= lim[1] & x <= lim[2]),
                                y = y, fill = "95% C.I."),
                  alpha = 0.5) +
        scale_color_manual(values = c("black",
                                      color)) +
        scale_fill_manual(values = fill_color)

    return(g)
}
colors <- c("blue", "red", "#039303")


# EXERCISE 1

# number of dead soldiers
y1  <- c(0, 1, 2, 3, 4, 5)
# observations
n1 <- vector("list", length(2))
n1[[1]] <- c(109, 65, 22, 3, 1, 0)
n1[[2]] <- c(144, 91, 32, 11, 2, 0)
# priors
priors <- c("Uniform", "Jeffrey's")
p_param <- c(1, 0.5) # parameter to sum to sum(y) to obtain
                     # the alpha parameter of posterior

delta <- 0.001
lambda <- seq(0, 3, delta)


# PART 1 & 2 :
# assuming a uniform and Jeffrey's priors, compute and plot the posterior
# distribution for lambda, the death rate over the measurement time.
# Determine the posterior mean, median and variance, and compute the
# 95% credibility interval.

# As the uniform prior we can use the improper prior
# g(lambda) = 1 for lambda > 0
# The posterior thus becomes proportional to the likelihood, a Gamma
# distribution with: alpha = sum(y) + 1 and beta = n
# Jeffrey's prior is g(lambda) = 1/sqrt(lambda) for lambda > 0
# The posterior thus becomes proportional to the likelihood, a Gamma
# distribution with: alpha = sum(y) + 0.5 and beta = n

l <- length(n1) * length(priors)
post <- vector("list", l)
mean <- vector("list", l)
med  <- vector("list", l)
var  <- vector("list", l)
lim  <- vector("list", l)
ls_g <- vector("list", l)

for (i in seq(1, length(n1))) {
    for (j in seq(1, length(priors))) {
        # list index
        pos <- (i - 1) * length(n1) + j

        # compute posterior
        alpha <- sum(y1 * n1[[i]]) + p_param[j]
        beta <- sum(n1[[i]])
        post[[i]] <- dgamma(lambda, alpha, beta)
        post[[i]] <- post[[i]] / (delta * sum(post[[i]]))

        # mean, median and variance
        mean[[pos]] <- delta * sum(lambda * post[[i]])
        var[[pos]]  <- delta * sum((lambda - mean[[pos]])^2 * post[[i]])
        med[[pos]]  <- quant(0.5, lambda, post[[i]], delta)

        # 95% C. I.
        lim[[pos]] <- cred_interval(0.95, lambda, post[[i]], delta)

        # plot results
        ls_g[[pos]] <- plotGraph(lambda, post[[i]], mean[[pos]],
                                 lim[[pos]], colors[j], c(0.25, 1.25)) +
                       labs(x = "\u03bb", y = "P(\u03bb | y, M)",
                            title = paste("   From ", priors[j], " Prior"),
                            color = "", fill = "")
    }
}

g1 <- ggarrange(ls_g[[1]], ls_g[[3]],
                labels = c("n1", "n2"),
                ncol = 1, nrow = 2,
                common.legend = TRUE)
g2 <- ggarrange(ls_g[[2]], ls_g[[4]],
                labels = c("n1", "n2"),
                ncol = 1, nrow = 2,
                common.legend = TRUE)
g_post <- ggarrange(g1, g2,
                    labels = c("", ""),
                    ncol = 2, nrow = 1)
g_post <- annotate_figure(g_post,
                          top = text_grob("Ex. 4 : Compute Posteriors"))


# ---------------------------------------------------------------------- #


# EXERCISE 3

p <- seq(0, 1, delta)

# observations
n3 <- c(116, 165)
y3 <- c(11, 9)


# PART a-e : frequentist estimator for p

p_freq <- y3 / n3


# PART b-f : using a Beta(1; 10) prior for p, calculate and posterior
#            distribution P(p | y). For the second measurement, assume also
#            the posterior probability of the older measurement as the prior
#            for the new one.

priors3 <- c("Beta(1, 10)", "Beta(1, 10)",
             "Posterior of Previous Measurement")
post3 <- vector("list", 3)
alpha_prior <- 1
beta_prior <- 10
# compute posterior using beta distribution as prior
for (i in seq(1, 2)) {
    alpha <- alpha_prior + y3[i]
    beta  <- beta_prior + n3[i] - y3[i]
    post3[[i]] <- dbeta(x = p, shape1 = alpha, shape2 = beta)
    post3[[i]] <- post3[[i]] / (delta * sum(post3[[i]]))
}
# compute posterior for the second observation using the
# posterior of the first as a prior
post3[[3]] <- dbinom(y3[2], n3[2], p) * post3[[1]]
post3[[3]] <- post3[[3]] / (delta * sum(post3[[3]]))


# PART c-g : find the bayesian estimator for p, the posterior mean and
#            variance, and a 95% credible interval

mean3 <- vector("list", 3)
med3  <- vector("list", 3)
var3  <- vector("list", 3)
lim3  <- vector("list", 3)
for (i in seq(1, 3)) {
    # mean, median and variance
    mean3[[i]] <- delta * sum(p * post3[[i]])
    var3[[i]]  <- delta * sum((p - mean3[[i]])^2 * post3[[i]])
    med3[[i]]  <- quant(0.5, p, post3[[i]], delta)

    # 95% C. I.
    lim3[[i]] <- cred_interval(0.95, p, post3[[i]], delta)
}


# PART d-h : test the hypotesis H0 : p = 0.1 versus H1 : p != 0.1
#            at 5% level of significance with both the frequentist
#            and bayesian approach

twoSidesPval <- function(x, x0, p) {
    n <- length(p)
    p_val <- rep(0, n)
    for (j in seq(1, n)) {
        p_val[j] <- ifelse(x[j] < x0, sum(p[1: j]),
                                      sum(p[j: n]))
    }
    return(p_val)
}

rejectRegion <- function(x, x0, p, a) {
    p_val <- twoSidesPval(x, x0, p)
    idx <- which(x == floor(x0))
    lim <- c(x[idx], x[idx + 1], 1)
    i <- 1
    while ((idx - i) > 0 && (idx + 1 + i) < length(x)) {
        sum <- p_val[idx - i] + p_val[idx + 1 + i]
        if (sum <= a) {
            sum_left  <- p_val[idx - i + 1] + p_val[idx + 1 + i]
            sum_right <- p_val[idx - i] + p_val[idx + i]
            if (sum_left > sum_right && sum_left <= a) {
                lim[1] <- x[idx - i + 1]
                lim[2] <- x[idx + 1 + i]
                lim[3] <- sum_left
            } else if (sum_right >= sum_left && sum_right <= a) {
                lim[1] <- x[idx - i]
                lim[2] <- x[idx + i]
                lim[3] <- sum_right
            } else {
                lim[1] <- x[idx - i]
                lim[2] <- x[idx + 1 + i]
                lim[3] <- sum
            }
            break
        }
        i <- i + 1
    }
    return(lim)
}

# null hypothesis
p0 <- 0.1
# significance
a <- 0.05

# Frequentist test
probs <- vector("list", 2)
rlims <- vector("list", 2)
ls_g_freq  <- vector("list", 2)
for (i in seq(1, 2)) {
    x <- seq(0, n3[i])
    probs[[i]] <- dbinom(x, n3[i], p0)
    rlims[[i]] <- rejectRegion(x, n3[i] * p0, probs[[i]], a)

    # plot results
    df <- subset(data.frame(a = x, b = probs[[i]]), a < 35)
    cat(c(y3[i], "\n"))
    f <- y3[i]
    ls_g_freq[[i]] <- ggplot() +
            geom_bar(data = subset(df, a >  rlims[[i]][1] & a <  rlims[[i]][2]),
                     mapping = aes(x = a, y = b), color = "black",
                     stat = "identity", fill = "lightblue", alpha = 0.5) +
            geom_bar(data = subset(df, a == y3[i]),
                     mapping = aes(x = a, y = b,
                                   fill = "Observed Data"),
                     stat = "identity",  color = "black", alpha = 0.5) +
            geom_bar(data = subset(df, a <= rlims[[i]][1] | a >= rlims[[i]][2]),
                     mapping = aes(x = a, y = b,
                                   fill = "Rejection Region\n\u03b1 = 5%"),
                     stat = "identity",  color = "black", alpha = 0.5) +
            annotate(geom = "text", x = 28, y = 0.01,
                     label = paste("Area = ", round(rlims[[i]][3], digits = 3)
                                    * 100, "%"),
                     color = "red") +
            scale_fill_manual(values = c("#71de71", "#f45c5c")) +
            labs(x = "# High Bacteria Samples", y = "P(y | n, M)",
                 title = "", color = "", fill = "")
}

# Bayesian test
ls_g_bayes <- vector("list", 3)
for (i in seq(1, 3)) {
    ls_g_bayes[[i]] <- plotGraph(p, post3[[i]], mean3[[i]],
                                 lim3[[i]], colors[i], c(0, 0.25)) +
                       geom_vline(aes(xintercept = p0, linetype = "p0"),
                                  color = "black", linewidth = 0.8) +
                       scale_linetype_manual(values = c(6)) +
                       labs(x = "p", y = "P(p | y, n, M)",
                            title = paste("   From ", priors3[i], " Prior"),
                            color = "", fill = "", linetype = "")
}

g_freq <- ggarrange(plotlist = ls_g_freq,
                    labels = c("n1", "n2"),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE)
g_freq <- annotate_figure(g_freq,
                          top = text_grob("Ex. 3 : Frequentist Test"))

g_bayes <- ggarrange(plotlist = ls_g_bayes,
                    labels = c("n1", "n2", "n2"),
                    ncol = 1, nrow = 3)
g_bayes <- annotate_figure(g_bayes,
                           top = text_grob("Ex. 3 : Bayesian Test"))



