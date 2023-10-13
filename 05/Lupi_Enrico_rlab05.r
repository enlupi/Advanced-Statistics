# Enrico Lupi, 2090596, 28 May 2023


# LABORATORY 5

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(knitr)
library(plotly)
library(rjags)
library(R2jags)


set.seed(4848)

# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- FALSE
if (save_pdf) {
    pdf("Lupi_Enrico_rlab05_Outputs.pdf")
}

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
plotGraph <- function(x, y, mean, lim, color, fill_color, range) {
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
fill_colors <- c("lightblue", "#f45c5c", "#71de71")


# ---------------------------------------------------------------------- #


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
ls_g <- vector("list", l)

cnames <- c("Mean", "St Dev", "Median", "95% C.I. sx", "95% C.I. dx")
rnames <- c("n1 Uniform", "n1 Jeffrey's", "n2 Uniform", "n2 Jeffrey's")
stat1 <- data.frame(matrix(nrow = length(rnames), ncol = length(cnames)))
colnames(stat1) <- cnames
rownames(stat1) <- rnames

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
        stat1[pos, 1] <- delta * sum(lambda * post[[i]])
        stat1[pos, 2] <- sqrt(delta * sum((lambda - stat1[pos, 1])^2
                                    * post[[i]]))
        stat1[pos, 3] <- quant(0.5, lambda, post[[i]], delta)

        # 95% C. I.
        stat1[pos, 4:5] <- cred_interval(0.95, lambda, post[[i]], delta)

        # plot results
        ls_g[[pos]] <- plotGraph(lambda, post[[i]], stat1[pos, 1],
                                 c(stat1[pos, 4], stat1[pos, 5]),
                                 colors[j], fill_colors[j], c(0.25, 1.25)) +
                       labs(x = "lambda", y = "P(lambda | y, M)",
                            title = paste("   From ", priors[j], " Prior"),
                            color = "", fill = "")
    }
}

stat1 <- stat1 |> round(digits = 3)

#plots
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
                          top = text_grob("Ex. 1 : Compute Posteriors"))
if (save_images) {
    ggsave("Ex5-1_posteriors.png")
}
if (save_pdf) {
    print(g_post)
    grid.newpage()
    grid.table(stat1)
}


# ---------------------------------------------------------------------- #


# EXERCISE 2

# prepare BUGS models with different priors
model_string <- c("
    model {
        # data likelihood
        for (i in 1:length(X)) {
            X[i] ~ dpois(lambda);
        }

        # a (approximately) uniform prior for lambda
        lambda ~ dexp(0.00001)

        # Predicted data , given lambda
        Y ~ dpois(lambda);
    }
", "
    model {
        # data likelihood
        for (i in 1:length(X)) {
            X[i] ~ dpois(lambda);
        }

        # (approximately) Jeffrey's prior for lambda
        lambda ~ dgamma(0.5, 0.00001)

        # Predicted data , given lambda
        Y ~ dpois(lambda);
    }
")
bugs_files <- c("ex2_jags_poiss_model1.bug", "ex2_jags_poiss_model2.bug")
for (i in seq(1, 2)) {
    writeLines(model_string[i], con = bugs_files[i])
}

ls_glambda <- vector("list", l)
ls_gcounts <- vector("list", l)
ls_gcorr   <- vector("list", l)
stat2 <- data.frame(matrix(nrow = length(rnames), ncol = length(cnames)))
colnames(stat2) <- cnames
rownames(stat2) <- rnames

# run MCMC for different models and observations
for (i in seq(1, length(n1))) {
    for (j in seq(1, length(priors))) {
        # list index
        pos <- (i - 1) * length(n1) + j

        # create JAGS model
        data <- NULL
        data$X <- rep(y1, n1[[i]])
        model <- bugs_files[j]
        jm2 <- jags.model(model, data)

        # Update the Markov chain (Burn-in)
        update(jm2, 1000)

        # run chain
        chain <- coda.samples(jm2, c("lambda", "Y"), n.iter = 10000)
        chain_df <- as.data.frame(as.mcmc(chain))

        # print results
        plot(chain, col = "navy")
        print(summary(chain))
        print(cor(chain_df))

        # mean, variance and median
        stat2[pos, 1] <- summary(chain)[[1]]["lambda", "Mean"]
        stat2[pos, 2] <- summary(chain)[[1]]["lambda", "SD"]
        stat2[pos, 3] <- summary(chain)[[2]]["lambda", "50%"]

        # 95% C. I.
        stat2[pos, 4] <- summary(chain)[[2]]["lambda", "2.5%"]
        stat2[pos, 5] <- summary(chain)[[2]]["lambda", "97.5%"]

        # plots results
        ls_glambda[[pos]] <- chain_df |> ggplot() +
                    geom_histogram(mapping = aes(x = lambda, y = ..density..),
                                   bins = 100, color = "black",
                                   fill = fill_colors[j]) +
                    labs(x = "lambda", y = "P(lambda | y, M)",
                         title = paste("   From ", priors[j], " Prior"))

        ls_gcounts[[pos]] <- chain_df$Y |> table() |> data.frame() |> ggplot() +
                    geom_bar(mapping = aes(x = Var1, y = Freq / sum(Freq)),
                             color = "black", fill = "lightblue",
                             stat = "identity", alpha = 0.5) +
                    labs(x = "Y", y = "f (Y)",
                         title = paste("   From ", priors[j], " Prior"))

        ls_gcorr[[pos]] <- chain_df |> ggplot() +
                    geom_point(aes(x = lambda, y = Y), color = colors[j]) +
                    labs(x = "lambda", y = "Y",
                         title = paste("   From ", priors[j], " Prior"))
    }
}

stat2 <- stat2 |> round(digits = 3)

#plots
g_lambda <- ggarrange(plotlist = ls_glambda,
                      labels = c("n1", "n1", "n2", "n2"),
                      ncol = 2, nrow = 2)
g_lambda <- annotate_figure(g_lambda,
                          top = text_grob("Ex. 2 : Inference on Lambda"))
if (save_images) {
    ggsave("Ex5-2_lambda.png")
}

g_counts <- ggarrange(plotlist = ls_gcounts,
                      labels = c("n1", "n1", "n2", "n2"),
                      ncol = 2, nrow = 2)
g_counts <- annotate_figure(g_counts,
                          top = text_grob("Ex. 2 : Predicted Counts"))
if (save_images) {
    ggsave("Ex5-2_counts.png")
}

g_corr <- ggarrange(plotlist = ls_gcorr,
                      labels = c("n1", "n1", "n2", "n2"),
                      ncol = 2, nrow = 2)
g_corr <- annotate_figure(g_corr,
                          top = text_grob("Ex. 2 : Correlation Plot"))
if (save_images) {
    ggsave("Ex5-2_corr.png")
}
if (save_pdf) {
    print(g_lambda)
    print(g_counts)
    print(g_corr)
    grid.newpage()
    grid.table(stat2)
}


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

rnames3 <- rnames <- c("n1 Beta", "n2 Beta", "n2 Prev. Post.")
stat3 <- data.frame(matrix(nrow = length(rnames3), ncol = length(cnames)))
colnames(stat3) <- cnames
rownames(stat3) <- rnames3
for (i in seq(1, 3)) {
    # mean, median and variance
    stat3[i, 1] <- delta * sum(p * post3[[i]])
    stat3[i, 2] <- sqrt(delta * sum((p - stat3[i, 1])^2
                                * post3[[i]]))
    stat3[i, 3]  <- quant(0.5, p, post3[[i]], delta)

    # 95% C. I.
    stat3[i, 4:5] <- cred_interval(0.95, p, post3[[i]], delta)
}

stat3 <- stat3 |> round(digit = 3)


# PART d-h : test the hypotesis H0 : p = 0.1 versus H1 : p != 0.1
#            at 5% level of significance with both the frequentist
#            and bayesian approach

# compute probability that a random sample from distribution is smaller
# than x if x < x0 or greater than x if x > x0
twoSidesPval <- function(x, x0, p) {
    n <- length(p)
    p_val <- rep(0, n)
    for (j in seq(1, n)) {
        p_val[j] <- ifelse(x[j] < x0, sum(p[1: j]),
                                      sum(p[j: n]))
    }
    return(p_val)
}

# computes rejection region for significance a in two-sides
# hypotheis testing
rejectRegion <- function(x, x0, p, a) {
    # compute left and right side "p-values"
    p_val <- twoSidesPval(x, x0, p)
    idx <- which(x == floor(x0))
    lim <- c(x[idx], x[idx + 1], 1)

    # find the rejection egion by checking if the integral in the proposed
    # rejection range is less than the significance of the test
    # start at the central value x0 and move outwards at each iteration
    i <- 1
    while ((idx - i) > 0 && (idx + 1 + i) < length(x)) {
        sum <- p_val[idx - i] + p_val[idx + 1 + i]
        if (sum <= a) {
            sum_left  <- p_val[idx - i + 1] + p_val[idx + 1 + i]
            sum_right <- p_val[idx - i] + p_val[idx + i]
            # check if adding a bin on the left or right to increase the
            # rejection range is acceptable
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
                     stat = "identity", fill = "lightgrey", alpha = 0.5) +
            geom_bar(data = subset(df, a == y3[i]),
                     mapping = aes(x = a, y = b,
                                   fill = "Observed Value"),
                     stat = "identity",  color = "black", alpha = 0.5) +
            geom_bar(data = subset(df, a <= rlims[[i]][1] | a >= rlims[[i]][2]),
                     mapping = aes(x = a, y = b,
                                   fill = "Rejection Region\n\u03b1 = 5%"),
                     stat = "identity",  color = "black", alpha = 0.5) +
            annotate(geom = "text", x = 30, y = 0.01,
                     label = paste("Area = ", round(rlims[[i]][3], digits = 3)
                                    * 100, "%"),
                     color = "red") +
            scale_fill_manual(values = c("#71de71", "#f45c5c")) +
            labs(x = "# High Bacteria Samples", y = "P(y | p0, n, M)",
                 title = "", color = "", fill = "")
}

# Observed values is always in the accepted region, so we cannot
# reject the null hypothesis in any observations. For the second one,
# though, the observed value is at the very edge of the rejection region
# so the result is not very strong.


# Bayesian test
ls_g_bayes <- vector("list", 3)
for (i in seq(1, 3)) {
    ls_g_bayes[[i]] <- plotGraph(p, post3[[i]], stat3[i, 1],
                                 c(stat3[i, 4], stat3[i, 5]),
                                 colors[i], fill_colors[i], c(0, 0.25)) +
                       geom_vline(aes(xintercept = p0, linetype = "p0"),
                                  color = "black", linewidth = 0.8) +
                       scale_linetype_manual(values = c(6)) +
                       labs(x = "p", y = "P(p | y, n, M)",
                            title = paste("   From ", priors3[i], " Prior"),
                            color = "", fill = "", linetype = "")
}

# For the first measurement p0 is well within the credibility interval
# of the posterior, so we cannot reject the null hypothesis.
# On the other hand, the result changes for the second measurement
# depending on the prior used: using a Beta(1, 10) prior p0 lies outside
# the credibility interval so we can reject the null hypothesis, while
# using the posterior of the first measurement as a prior it is contained
# in the credibility interval once again.


g_freq <- ggarrange(plotlist = ls_g_freq,
                    labels = c("n1", "n2"),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE)
g_freq <- annotate_figure(g_freq,
                          top = text_grob("Ex. 3 : Frequentist Test"))
if (save_images) {
    ggsave("Ex5-3_freq.png")
}

g_bayes <- ggarrange(plotlist = ls_g_bayes,
                    labels = c("n1", "n2", "n2"),
                    ncol = 1, nrow = 3)
g_bayes <- annotate_figure(g_bayes,
                           top = text_grob("Ex. 3 : Bayesian Test"))
if (save_images) {
    ggsave("Ex5-3_bayes.png")
}
if (save_pdf) {
    print(g_freq)
    print(g_bayes)
    grid.newpage()
    grid.table(stat3)
}


# ---------------------------------------------------------------------- #


# EXERCISE 4

# prepare BUGS model
model_string4 <- "
    model {
        # data likelihood
        for (i in 1:length(X)) {
            X[i] ~ dbern(p);
        }

        # a Beta prior for p
        p ~ dbeta(1, 10)

        # Predicted data, given p
        Y ~ dbin(p, n_next);
    }
"
bugs_file4 <- "ex2_jags_binom_model.bug"
writeLines(model_string4, con = bugs_file4)


ls_gp <- vector("list", 2)
ls_gcounts4 <- vector("list", 2)
ls_gcorr4   <- vector("list", 2)
rnames4 <- c("n1 Beta", "n2 Beta")
stat4 <- data.frame(matrix(nrow = length(rnames4), ncol = length(cnames)))
colnames(stat4) <- cnames
rownames(stat4) <- rnames4

# run MCMC for different models and observations
for (i in seq(1, length(n3))) {
    # create JAGS model
    data4 <- NULL
    data4$X <- rep(c(0, 1), c(n3[i] - y3[i], y3[i]))
    data4$n_next <- 10 # Predictions
    model4 <- bugs_file4
    jm4 <- jags.model(model4, data4)

    # Update the Markov chain (Burn-in)
    update(jm4, 1000)

    # run chain
    chain <- coda.samples(jm4, c("p", "Y"), n.iter = 10000)
    chain_df <- as.data.frame(as.mcmc(chain))

    # print results
    plot(chain, col = "navy")
    print(summary(chain))
    print(cor(chain_df))

    # mean, variance and median
    stat4[i, 1] <- summary(chain)[[1]]["p", "Mean"]
    stat4[i, 2] <- summary(chain)[[1]]["p", "SD"]
    stat4[i, 3] <- summary(chain)[[2]]["p", "50%"]

    # 95% C. I.
    stat4[i, 4] <- summary(chain)[[2]]["p", "2.5%"]
    stat4[i, 5] <- summary(chain)[[2]]["p", "97.5%"]

    # plots results
    ls_gp[[i]] <- chain_df |> ggplot() +
                    geom_histogram(mapping = aes(x = p, y = ..density..),
                                   bins = 100, color = "black",
                                   fill = fill_colors[i]) +
                    labs(x = "p", y = "P(p | y, M)",
                         title = paste("   From ", priors3[i], " Prior"))

    ls_gcounts4[[i]] <- chain_df$Y |> table() |> data.frame() |> ggplot() +
                    geom_bar(mapping = aes(x = Var1, y = Freq / sum(Freq)),
                             color = "black", fill = "lightblue",
                             stat = "identity", alpha = 0.5) +
                    labs(x = "Y", y = "f (Y)",
                         title = paste("   From ", priors3[i], " Prior"))

    ls_gcorr4[[i]] <- chain_df |> ggplot() +
                    geom_point(aes(x = p, y = Y), color = colors[i]) +
                    labs(x = "p", y = "Y",
                         title = paste("   From ", priors3[i], " Prior"))
}

stat4 <- stat4 |> round(digits = 3)

#plots
g_p <- ggarrange(plotlist = ls_gp,
                 labels = c("n1", "n2"),
                 ncol = 1, nrow = 2)
g_p <- annotate_figure(g_p,
                          top = text_grob("Ex. 4 : Inference on p"))
if (save_images) {
    ggsave("Ex5-4_p.png")
}

g_counts4 <- ggarrange(plotlist = ls_gcounts4,
                       labels = c("n1", "n2"),
                       ncol = 1, nrow = 2)
g_counts4 <- annotate_figure(g_counts4,
                          top = text_grob("Ex. 4 : Predicted Counts"))
if (save_images) {
    ggsave("Ex5-4_counts.png")
}

g_corr4 <- ggarrange(plotlist = ls_gcorr4,
                     labels = c("n1", "n2"),
                     ncol = 1, nrow = 2)
g_corr4 <- annotate_figure(g_corr4,
                          top = text_grob("Ex. 4 : Correlation Plot"))
if (save_images) {
    ggsave("Ex5-4_corr.png")
}
if (save_pdf) {
    print(g_p)
    print(g_counts4)
    print(g_corr4)
    grid.newpage()
    grid.table(stat4)
}




if (save_pdf) {
    dev.off()
}
