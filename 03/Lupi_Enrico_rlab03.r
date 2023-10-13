# Enrico Lupi, 2090596, 7 May 2023

# LABORATORY 2

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)

# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- FALSE

# Probability values
p <- seq(0, 1, length.out = 201)
n_p <- length(p)
delta_p <- p[2] - p[1]


# EXERCISE 1

# study the binomial inference for a study that reports y = 7
# successes in n = 20 independent trial

y <- 7
n <- 20

# Using a uniform distribution, the posterior pdf is simply
# proportional to the Likelihood, a Binomial distribution with
# n = 20 and y successes
post_unif <-  dbinom(x = y, size = n, prob = p)
# We still have to normalise it as it was normalised in y and not in p
post_unif <- post_unif / (delta_p * sum(post_unif))


# The Jeffrey's prior for Bernoulli trials is the arcsin distribution,
# that is a Beta distribution with
alpha_prior <- 0.5
beta_prior <- 0.5
# Using a Beta prior, the posterior is a Beta distribution with
alpha <- alpha_prior + y
beta <- beta_prior + n - y

post_beta <- dbeta(x = p, shape1 = alpha, shape2 = beta)
post_beta <- post_beta / (delta_p * sum(post_beta))


# Using an arbitrary step function as a prior
prior_g <- function(p) {
    ifelse(p <= 0.2, p,
        ifelse(0.2 < p && p <= 0.3, 0.2,
            ifelse(0.3 < p && p <= 0.5, 0.5 - p, 0)))
}

# To obtain the posterior we need to multiply the prior by the likelihood
# (a Bernoulli distribution) and then normalise

post_g <- rep(0, length(p))
for (i in seq(1, length(post_g))){
    post_g[i] <- dbinom(x = y, size = n, prob = p[i]) * prior_g(p[i])
}
post_g <- post_g / (delta_p * sum(post_g))


# PART 1 : plot the posterior distribution and summarize the results
#          computing the first two moments

# We will compute the mean and variance numerically using their definitions

mean_unif <- delta_p * sum(p * post_unif)
var_unif  <- delta_p * sum((p - mean_unif)^2 * post_unif)

mean_beta <- delta_p * sum(p * post_beta)
var_beta  <- delta_p * sum((p - mean_beta)^2 * post_beta)

mean_g <- delta_p * sum(p * post_g)
var_g  <- delta_p * sum((p - mean_g)^2 * post_g)

# Plot the results

g_posterior <- ggplot() +
               geom_line(aes(x = p, y = post_unif, color = "Uniform Prior"),
                         linewidth = 1, linetype = 1) +
               geom_line(aes(x = p, y = post_beta, color = "Jeffrey's Prior"),
                         linewidth = 1, linetype = 1) +
               geom_line(aes(x = p, y = post_g,    color = "Arbitrary Prior g"),
                         linewidth = 1, linetype = 1) +
               geom_vline(xintercept = mean_unif,
                          linetype = 2, color = "red") +
               geom_vline(xintercept = mean_beta,
                          linetype = 2, color = "blue") +
               geom_vline(xintercept = mean_g,
                          linetype = 2, color = "#00a100") +
               scale_color_manual(values = c("#00a100",
                                             "blue",
                                             "red")) +
               labs(x = "p", y = "P(p|y, n, M)",
                    title = "Posterior", color = "")

if (save_images) {
    ggsave("Ex3-1_posteriors.png", g_posterior)
}


# PART 2 : compute a 95% credibility interval and give
#          the results in a summary table

# We will compute a symmetrical credibility interval with equal probability
# in the two tails, thus centered around the median

cred_interval <- function(val, p, post, dp = delta_p) {
    limits <- c(0, 0)
    # compute the integral in each tail
    v <- (1 - val) / 2
    sum <- 0
    for (i in seq(1, length(p))) {
        sum <- sum + dp * post[i]
        # take the first value when the integral surpasses the
        # desired value
        if (sum >= v) {
            limits[1] <- p[i]
            break
        }
    }
    sum <- 0
    # do the same but for the right tail
    for (i in seq(length(p), 1, -1)) {
        sum <- sum + dp * post[i]
        if (sum >= v) {
            limits[2] <- p[i]
            break
        }
    }
    return(limits)
}

lim_unif <- cred_interval(0.95, p, post_unif)
lim_beta <- cred_interval(0.95, p, post_beta)
lim_g    <- cred_interval(0.95, p, post_g)

val_unif <- c(mean_unif, var_unif, lim_unif) |> round(digits = 3)
val_beta <- c(mean_beta, var_beta, lim_beta) |> round(digits = 3)
val_g    <- c(mean_g,    var_g,    lim_g)    |> round(digits = 3)

df <- data.frame(val_unif, val_beta, val_g)
rownames(df) <- c("Mean", "Var", "95% C.I. Left", "95% C.I. Right")
colnames(df) <- c("Uniform P.", "Jeffrey's P.", "Arbitrary P. g")


# PART 3 : draw the limits on the plot of the posterior distribution

g_unif <- ggplot() +
          geom_line(aes(x = p, y = post_unif), color = "red",
                    linewidth = 1, linetype = 1) +
          geom_vline(xintercept = mean_unif,
                     linetype = "dashed", color = "red") +
          geom_vline(xintercept = lim_unif[1], linewidth = 1,
                     linetype = "dashed", color = "black") +
          geom_vline(xintercept = lim_unif[2], linewidth = 1,
                     linetype = "dashed", color = "black") +
          labs(x = "p", y = "P(p|y, n, M)",
                    title = "Uniform Prior", color = "")

g_beta <- ggplot() +
          geom_line(aes(x = p, y = post_beta), color = "blue",
                    linewidth = 1.5, linetype = 1) +
          geom_vline(xintercept = mean_beta,
                     linetype = "dashed", color = "blue") +
          geom_vline(xintercept = lim_beta[1], linewidth = 1,
                     linetype = "dashed", color = "black") +
          geom_vline(xintercept = lim_beta[2], linewidth = 1,
                     linetype = "dashed", color = "black") +
          labs(x = "p", y = "P(p|y, n, M)",
                    title = "Jeffrey's Prior", color = "")

g_g    <- ggplot() +
          geom_line(aes(x = p, y = post_g),    color = "#00a100",
                    linewidth = 1.5, linetype = 1) +
          geom_vline(xintercept = mean_g,
                     linetype = "dashed", color = "#00a100") +
          geom_vline(xintercept = lim_g[1], linewidth = 1,
                     linetype = "dashed", color = "black") +
          geom_vline(xintercept = lim_g[2], linewidth = 1,
                     linetype = "dashed", color = "black") +
          labs(x = "p", y = "P(p|y, n, M)",
                    title = "Arbitrary Prior", color = "")

g_cred <- ggarrange(g_unif, g_beta, g_g,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)

if (save_images) {
    ggsave("Ex3-1_CIplot.png", width = 6300, units = "px")
}


# ---------------------------------------------------------------------- #


# EXERCISE 2


# This is once again a case of inference for Binomial distribution with n trials
# and y successes (when the sample contains Giard cystis)
n2 <- 116
y2 <- 17

# We can use again the results of the previous exercise
# Using a uniform distribution, the posterior pdf is simply
# proportional to the Likelihood, a Binomial distribution with
# n = 116 and y = 17 successes
post2_unif <- dbinom(x = y2, size = n2, prob = p)
post2_unif <- post2_unif / (delta_p * sum(post2_unif))


# Using a Beta prior with
alpha2_prior <- 1
beta2_prior <- 4
# the posterior is a Beta distribution with
alpha2 <- alpha2_prior + y2
beta2 <- beta2_prior + n2 - y2

post2_beta <- dbeta(x = p, shape1 = alpha2, shape2 = beta2)
post2_beta <- post2_beta / (delta_p * sum(post2_beta))


# PART 1 : plot the posterior distribution and summarize the results
#          computing the first two moments

mean2_unif <- delta_p * sum(p * post2_unif)
var2_unif  <- delta_p * sum((p - mean2_unif)^2 * post2_unif)

mean2_beta <- delta_p * sum(p * post2_beta)
var2_beta  <- delta_p * sum((p - mean2_beta)^2 * post2_beta)

# The resulting plot is
g_posterior2 <- ggplot() +
                geom_line(aes(x = p, y = post2_unif, color = "Uniform Prior"),
                          linewidth = 1, linetype = 1) +
                geom_line(aes(x = p, y = post2_beta,
                              color = "Beta(1, 4) Prior"),
                          linewidth = 1, linetype = 1) +
                geom_vline(xintercept = mean2_unif,
                           linetype = "dashed", color = "red") +
                geom_vline(xintercept = mean2_beta,
                           linetype = "dashed", color = "blue") +
                scale_color_manual(values = c("blue",
                                              "red")) +
                labs(x = "p", y = "P(p|y, n, M)",
                     title = "Posterior", color = "")

if (save_images) {
    ggsave("Ex3-2_posterior.png", width = 4200, units = "px")
}


# PART 2 : find a normal approximation for the posterior g(p|y)

# We can use as a normal approximation for the posterior the gaussian curve
# that has the same mean and variation of the posterior. Since they are very
# similar, we can use the posterior obtained either with a uniform prior or
# with a Beta(1,4) prior without any effective difference

post2_approx <- dnorm(p, mean = mean2_unif, sd = sqrt(var2_unif))

g_approx <- ggplot() +
            geom_line(aes(x = p, y = post2_unif, color = "Uniform Prior"),
                      linewidth = 1, linetype = 1) +
            geom_line(aes(x = p, y = post2_approx,
                          color = "Normal Approximation"),
                      linewidth = 1, linetype = 1) +
            geom_vline(xintercept = mean2_unif,
                       linetype = "dashed", color = "black") +
            scale_color_manual(values = c("#00a100",
                                          "red")) +
            labs(x = "p", y = "P(p|y, n, M)",
                 title = "Posterior", color = "")

if (save_images) {
    ggsave("Ex3-2_normalApprox.png")
}


# PART 3 : compute a 95% credibility interval and give
#          the results in a summary table

lim2_unif   <- cred_interval(0.95, p, post2_unif)
lim2_beta   <- cred_interval(0.95, p, post2_beta)
lim2_approx <- cred_interval(0.95, p, post2_approx)

val2_unif   <- c(mean2_unif, var2_unif,   lim2_unif) |> round(digits = 3)
val2_beta   <- c(mean2_beta, var2_beta,   lim2_beta) |> round(digits = 3)
val2_approx <- c(mean2_unif, var2_unif, lim2_approx) |> round(digits = 3)

df2 <- data.frame(val2_unif, val2_beta, val2_approx)
rownames(df2) <- c("Mean", "Var", "95% C.I. Left", "95% C.I. Right")
colnames(df2) <- c("Uniform P.", "Beta(1,4) P.", "Normal Approximation")


# PART 4 : draw the limits on the plot of the posterior distribution

g_unif2 <- ggplot() +
           geom_line(aes(x = p, y = post2_unif), color = "red",
                     linewidth = 1, linetype = 1) +
           geom_vline(xintercept = mean2_unif,
                      linetype = "dashed", color = "red") +
           geom_vline(xintercept = lim2_unif[1], linewidth = 1,
                      linetype = "dashed", color = "black") +
           geom_vline(xintercept = lim2_unif[2], linewidth = 1,
                      linetype = "dashed", color = "black") +
           labs(x = "p", y = "P(p|y, n, M)",
                title = "Uniform Prior", color = "")

g_beta2 <- ggplot() +
           geom_line(aes(x = p, y = post2_beta), color = "blue",
                     linewidth = 1, linetype = 1) +
           geom_vline(xintercept = mean2_beta,
                      linetype = "dashed", color = "blue") +
           geom_vline(xintercept = lim2_beta[1], linewidth = 1,
                      linetype = "dashed", color = "black") +
           geom_vline(xintercept = lim2_beta[2], linewidth = 1,
                      linetype = "dashed", color = "black") +
           labs(x = "p", y = "P(p|y, n, M)",
                title = "Beta(1,4) Prior", color = "")

g_approx2 <- ggplot() +
             geom_line(aes(x = p, y = post2_approx), color = "#00a100",
                      linewidth = 1, linetype = 1) +
             geom_vline(xintercept = mean2_unif,
                        linetype = "dashed", color = "#00a100") +
             geom_vline(xintercept = lim2_approx[1], linewidth = 1,
                        linetype = "dashed", color = "black") +
             geom_vline(xintercept = lim2_approx[2], linewidth = 1,
                        linetype = "dashed", color = "black") +
             labs(x = "p", y = "P(p|y, n, M)",
                  title = "Normal Approximation", color = "")

g_cred2 <- ggarrange(g_unif2, g_beta2, g_approx2,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3)

if (save_images) {
    ggsave("Ex3-2_CIplot.png", width = 6300, units = "px")
}


# ---------------------------------------------------------------------- #


# EXERCISE 3

# Once again, we are dealing with a Binomial distribution with n trials
# and y successes (Head outcome in the coin toss)
n3 <- 30
y3 <- 15


# PART 1 : assuming a flat prior, and a beta prior, plot the likelihood
#          prior and posterior distributions for the data set

# The likelihood is a Binomial distribution with the number of successes
# as the variable and their probability as a fixed parameter.
# We can also compute what is the likelihood of obtaining y3 successes
# as a function of the probability p. As it is not a distribution in p
# we need to normalise it

like <- dbinom(x = y3, size = n3, prob = p)
like <- like / (delta_p * sum(like))

# We will use a uniform U(0,1) prior
prior_unif <- dunif(x = p, min = 0, max = 1)
# and a Beta prior with
alpha3_prior <- 2
beta3_prior <- 2
# as we assume that a fair coin is more probable
prior_beta <- dbeta(x = p, shape1 = alpha3_prior, shape2 = beta3_prior)

# The posterior distrbution are, respectively, proportional to the likelihood
post3_unif <- like
# and a Beta distribuion with
alpha3 <- alpha3_prior + y3
beta3 <- beta3_prior + n3 - y3
post3_beta <- dbeta(x = p, shape1 = alpha3, shape2 = beta3)
post3_beta <- post3_beta / (delta_p * sum(post3_beta))

# plotting the results we have
g_unif3 <-  ggplot() +
            geom_line(aes(x = p, y = post3_unif, color = "Posterior"),
                      linewidth = 1, linetype = 1, alpha = 0.5) +
            geom_line(aes(x = p, y = prior_unif, color = "Prior"),
                      linewidth = 1, linetype = 2) +
            geom_line(aes(x = p, y = like, color = "Likelihood"),
                      linewidth = 1, linetype = 3) +
            scale_color_manual(values = c("blue",
                                          "red",
                                          "#00a100")) +
            labs(x = "p", y = "Density",
                 title = "Uniform Prior", color = "")

g_beta3 <-  ggplot() +
            geom_line(aes(x = p, y = post3_beta, color = "Posterior"),
                      linewidth = 1, linetype = 1, alpha = 0.5) +
            geom_line(aes(x = p, y = prior_beta, color = "Prior"),
                      linewidth = 1, linetype = 2) +
            geom_line(aes(x = p, y = like, color = "Likelihood"),
                      linewidth = 1, linetype = 3) +
            scale_color_manual(values = c("blue",
                                          "red",
                                          "#00a100")) +
            labs(x = "p", y = "Density",
                 title = "Beta(2, 2) Prior", color = "")

g_post3 <- ggarrange(g_unif3, g_beta3,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE)

if (save_images) {
    ggsave("Ex3-3_posterior.png", width = 4200, units = "px")
}


# PART 2 : evaluate the most probable value for the coin probability p
#          and, integrating the posterior probability distribution,
#          give an estimate for a 95% credibility interval.

# The most probable outcome for the probability p is
mode3_unif <- p[which.max(post3_unif)]
mode3_beta <- p[which.max(post3_beta)]
# both return the same value, 0.5

# The 95% credibility interval limits are
lim3_unif   <- cred_interval(0.95, p, post3_unif)
lim3_beta   <- cred_interval(0.95, p, post3_beta)


val3_unif   <- c(mode3_unif, lim3_unif) |> round(digits = 3)
val3_beta   <- c(mode3_beta, lim3_beta) |> round(digits = 3)

df3 <- data.frame(val3_unif, val3_beta)
rownames(df3) <- c("Mode", "95% C.I. Left", "95% C.I. Right")
colnames(df3) <- c("Uniform P.", "Beta(2,2) P.")


# PART 3 : Repeat the same analysis assuming a sequential analysis of the
#          data. Show how the most probable value and the credibility
#          interval change as a function of the number of coin tosses

# get number of successes after every trial
ns <- seq(1, n3)
ys <- rep(0, n3)
res <- strsplit("TTTTTHTTHHTTHHHTHTHTHHTHTHTHHH", "")[[1]]
for (i in seq(1, length(res))) {
  if (res[i] == "H") {
    ys[i] <- 1
    }
}

modes_unif   <- rep(0, n3)
lims_unif_sx <- rep(0, n3)
lims_unif_dx <- rep(0, n3)
modes_beta   <- rep(0, n3)
lims_beta_sx <- rep(0, n3)
lims_beta_dx <- rep(0, n3)

for (i in ns) {
    # compute likelihood of the i-th extraction
    likelihood <- dbinom(x = ys[i], size = 1, prob = p)

    # compute posteriors
    posts_unif <- likelihood * prior_unif
    posts_unif <- posts_unif / (delta_p * sum(posts_unif))

    posts_beta <- likelihood * prior_beta
    posts_beta <- posts_beta / (delta_p * sum(posts_beta))

    modes_unif[i] <- p[which.max(posts_unif)]
    modes_beta[i] <- p[which.max(posts_beta)]
    lims_unif <- cred_interval(0.95, p, posts_unif)
    lims_unif_sx[i] <- lims_unif[1]
    lims_unif_dx[i] <- lims_unif[2]
    lims_beta <- cred_interval(0.95, p, posts_beta)
    lims_beta_sx[i] <- lims_unif[1]
    lims_beta_dx[i] <- lims_unif[2]

    # update prior for the next loop using current posterior
    prior_unif <- posts_unif
    prior_beta <- posts_beta
}

g_iteru <- ggplot() +
           geom_errorbar(aes(x = ns, y = modes_unif,
                             ymin = lims_unif_sx, ymax = lims_unif_dx,
                             color = "Sequential C.I."),
                         width = 1, ) +
           geom_line(aes(x = ns, y = modes_unif),
                     linewidth = 0.5, linetype = 2) +
           geom_point(aes(x = ns, y = modes_unif,
                          color = "Sequential Mode"), size = 2) +
           geom_hline(aes(color = "One-step Mode",
                          yintercept = mode3_unif)) +
           geom_hline(aes(yintercept = lim3_unif[1],
                          color = "One-step C.I."),
                      linetype = "dashed") +
           geom_hline(yintercept = lim3_unif[2],
                      linetype = "dashed", color = "#ff5664") +
           scale_color_manual(values = c("#ff5664",
                                         "red",
                                         "#767676",
                                         "black")) +
           labs(x = "n", y = "",
                title = "Uniform Prior", color = "")

g_iterb <- ggplot() +
           geom_errorbar(aes(x = ns, y = modes_beta,
                             ymin = lims_beta_sx, ymax = lims_beta_dx,
                             color = "Sequential C.I."),
                         width = 1, ) +
           geom_line(aes(x = ns, y = modes_beta),
                     linewidth = 0.5, linetype = 2) +
           geom_point(aes(x = ns, y = modes_beta,
                          color = "Sequential Mode"), size = 2) +
           geom_hline(aes(color = "One-step Mode",
                          yintercept = mode3_beta)) +
           geom_hline(aes(yintercept = lim3_beta[1],
                          color = "One-step C.I."),
                      linetype = "dashed") +
           geom_hline(yintercept = lim3_beta[2],
                      linetype = "dashed", color = "#46b6dc") +
           scale_color_manual(values = c("#46b6dc",
                                         "blue",
                                         "#767676",
                                         "black")) +
           labs(x = "n", y = "",
                title = "Beta(2,2) Prior", color = "")

g_iter <- ggarrange(g_iteru, g_iterb,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2,
                    common.legend = FALSE)

if (save_images) {
    ggsave("Ex3-3_iter.png", width = 4200, units = "px")
}


# PART 4 : Do you get a different results, by analyzing the data sequentially
#          with respect to a one-step analysis?

# We get the same result, which is expected as we are doing the same operation
# from an analytical point of view. The binomial distribution is, save for
# normalisation factors, the product of the single Bernoulli probabilities,
# so we get the same result by multiplying the prior for the binomial in the
# one-step analysis or by multiplying the prior for the single Bernoulli
# probability at each step in the sequential analysis. The advantage of the
# sequential analysis is that, if we have a particular goal, we can check if we
# have reached it at every step and stop if so, saving resources.



# ---------------------------------------------------------------------- #


# EXERCISE 4

# set the constants of the experiment: we will use 6 boxes
# and perform a total of 30 extractions
n_trial <- 30
trial <- seq(0, n_trial)
n_boxes  <- 6
n_stones <- n_boxes - 1

# select randomly one box
set.seed(4848)
w <- sample(0:n_stones, 1)
# and sample one stone with replacement at each trial
e <- c(NA, rbinom(n_trial, 1, w / n_stones))
res <- "Extractions:"
for (i in seq(1, n_trial)) {
    if (e[i + 1] == 1) {
        res <- paste(res, "W")
    } else {
        res <- paste(res, "B")
    }
}
cat(res)

# Let us create a dataframe where to keep the trial number, the result of the
# extraction (0 for B and 1 for W) and the probability of having chosen
# the w-th box after each trial
H <- data.frame(trial, e)
cnames <- c("Trial", "E")
for (i in seq(1, n_boxes)){
    h <- rep(1 / n_boxes, n_trial + 1)
    H <- cbind(H, h)
    cnames <- c(cnames, paste("H", toString(i - 1), sep = ""))
}
colnames(H) <- cnames

# loop over the trials
for (i in seq(1, n_trial)) {
    sum <- 0
    # loop over the boxes
    for (j in seq(1:n_boxes)) {
        # compute the probability of getting the i-th result
        # for the j-th box
        p_e <- (j - 1) / n_stones
        if (e[i + 1] == 0) {
            p_e <- 1 - p_e
        }
        # update the probability of having selected the j-th box
        H[i + 1, j + 2] <- p_e * H[i, j + 2]
        # update normalisation factor
        sum <- sum + H[i + 1, j + 2]
    }
    # normlise probabilities
    H[i + 1, 3:(3 + n_stones)] <- H[i + 1, 3:(3 + n_stones)] / sum

    cat(c("After trial ", i, ":\n"))
    print(H[i + 1, 2:(3 + n_stones)])
    cat("\n")
}

# plot results
g0 <- ggplot(data = H) +
      geom_line(aes(x = Trial, y = H[, 3]),
                linewidth = 0.5, linetype = 2) +
      geom_point(aes(x = Trial, y = H[, 3]), size = 2) +
      labs(x = "# trial", y = "P(H)", title = colnames(H)[3])
g1 <- ggplot(data = H) +
      geom_line(aes(x = Trial, y = H[, 4]),
                linewidth = 0.5, linetype = 2) +
      geom_point(aes(x = Trial, y = H[, 4]), size = 2) +
      labs(x = "# trial", y = "P(H)", title = colnames(H)[4])
g2 <- ggplot(data = H) +
      geom_line(aes(x = Trial, y = H[, 5]),
                linewidth = 0.5, linetype = 2) +
      geom_point(aes(x = Trial, y = H[, 5]), size = 2) +
      labs(x = "# trial", y = "P(H)", title = colnames(H)[5])
g3 <- ggplot(data = H) +
      geom_line(aes(x = Trial, y = H[, 6]),
                linewidth = 0.5, linetype = 2) +
      geom_point(aes(x = Trial, y = H[, 6]), size = 2) +
      labs(x = "# trial", y = "P(H)", title = colnames(H)[6])
g4 <- ggplot(data = H) +
      geom_line(aes(x = Trial, y = H[, 7]),
                linewidth = 0.5, linetype = 2) +
      geom_point(aes(x = Trial, y = H[, 7]), size = 2) +
      labs(x = "# trial", y = "P(H)", title = colnames(H)[7])
g5 <- ggplot(data = H) +
      geom_line(aes(x = Trial, y = H[, 8]),
                linewidth = 0.5, linetype = 2) +
      geom_point(aes(x = Trial, y = H[, 8]), size = 2) +
      labs(x = "# trial", y = "P(H)", title = colnames(H)[8])

g_H <- ggarrange(g0, g1, g2, g3, g4, g5,
                 labels = "AUTO",
                 ncol = 3, nrow = 2)

if (save_images) {
    ggsave("Ex3-4_box.png", width = 6300, height = 6300, units = "px")
}



# Save output to pdf file
if (save_pdf) {
    pdf("Lupi_Enrico_rlab03_Outputs.pdf")
    print(g_posterior)
    grid.newpage()
    grid.table(df)
    print(g_cred)
    print(g_posterior2)
    print(g_approx)
    grid.newpage()
    grid.table(df2)
    print(g_cred2)
    print(g_post3)
    grid.newpage()
    grid.table(df3)
    print(g_iter)
    print(g_H)
    dev.off()
}