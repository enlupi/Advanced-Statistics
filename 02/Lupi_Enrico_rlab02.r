# Enrico Lupi, 2090596, 23 Apr 2023

# LABORATORY 2

library(FRACTION)
library(tidyverse)
library(ggpubr)
library(GoFKernel)

# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- FALSE


# EXERCISE 1


# PART 1

# Define the probability density function
dk <- function(n) {
    stopifnot(is.wholenumber(n) == TRUE) # checks if the input
                                         # variable is discrete
    return(ifelse(n > 0 & n < 6, floor(n) / 15, 0))
}

# Define the cumulative distribution function
pk <- function(n) {
    stopifnot(is.wholenumber(n) == TRUE)
    # create array with cdf values different than 0 and 1
    v <- seq(1, 4)
    cdf_val <- cumsum(dk(v))

    m <- rep(0, length(n))
    # if value is in [1, 4] gives corresponding cdf value, otherwise 0 or 1
    for (i in seq(1, length(n))){
        m[i] <- ifelse((n[i] > 0 & n[i] < 5), cdf_val[n[i]],
                                              ifelse(n[i] > 4, 1, 0))
    }
    return(m)
}


# PART 2

x <- -1:7
# create bar plots for pdf...

g_pdfbar <- ggplot(data.frame(x), aes(x = x, y = dk(x))) +
            geom_bar(stat = "identity", color = "black",
                     fill = "lightblue", alpha = 0.5) +
            labs(x = "k", y = "P(k)",
                 title = "Probability Distribution")

g_cdfbar <- ggplot(data.frame(x), aes(x = x, y = pk(x))) +
            geom_bar(stat = "identity", color = "black",
                     fill = "lightblue", alpha = 0.5) +
            labs(x = "k", y = "CDF(k)",
                 title = "Cumulative Distribution")

g_bar <- ggarrange(g_pdfbar, g_cdfbar,
                   labels = c("A", "B", "C"),
                   ncol = 2, nrow = 1)

if (save_images) {
    ggsave("Ex2-1_pdf-cdf.png", g_bar)
}


# PART 3

# We will compute the mean value and variance from their definition
k <- 1:5
exp_value <- sum(k * dk(k))
var <- sum((k - exp_value)^2 * dk(k))

cat(c("Exercise 1.3:", "\n",
      "The mean value is", exp_value, "\n",
      "The variance is", var, "\n"))


# PART 4

# Once again we use the definition
exp_value2 <- sum(k * (6 - k) * dk(k))

cat(c("Exercise 1.4:", "\n",
      "The expected value is", exp_value2, "\n"))


# PART 5

# Define the function that allows to sample random numbers from the
# distribution by inverting the cdf
rk <- function(n) {
    # compute cdf values
    x <- 1:5
    cdf_val <- pk(x)

    m <- rep(0, n)
    # create array of length n of uniform variables in [0,1]
    u <- runif(n)
    # get the cdf inverse
    for (i in 1:n) {
        for (j in x) {
            if (u[i] <= cdf_val[j]) {
                m[i] <- j
                break
            }
        }
    }
    return(m)
}


# PART 6

# Sample 100,000 events from the distribution
rnd_val <- rk(100000)
# Build the histogram
g_k <- ggplot(data.frame(rnd_val), aes(rnd_val)) +
       geom_histogram(aes(y = after_stat(density)),
                      breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
                      color = "black", fill = "lightblue", alpha = 0.5) +
       geom_step(data.frame(x = seq(-0.5, 6.5), y = dk(-1:6)),
                 mapping = aes(x = x, y = y, color = "red"),
                 direction = "vh", linewidth = 1, linetype = 2) +
       scale_color_discrete(labels = c("Expected Distribution")) +
       labs(x = "k", y = "P(k)",
            title = "Discrete Distribution", color = element_blank())

if (save_images) {
    ggsave("Ex2-1_HistK.png", g_k)
}


# ---------------------------------------------------------------------- #


# EXERCISE 2


# PART 1

# Define the tringular pdf between [a, b] with a peak in c
dtriang <- function(x, a = 0, b = 1, c = 0.5) {
    m <- rep(0, length(x))
    for (i in seq(1, length(x))){
        if (a <= x[i] && x[i] < c) {
            m[i] <- 2 * (x[i] - a) / ((b - a) * (c - a))
        } else if (c <= x[i] && x[i] <= b) {
            m[i] <- 2 * (b - x[i]) / ((b - a) * (b - c))
        }
    }
    return(m)
}

# plot the function with a = 0, b = 1, c = 0.5
x <- seq(-0.5, 1.5, 0.01)
g_triang <- ggplot(data.frame(x), aes(x)) +
            geom_function(fun = dtriang) +
            labs(x = "X", y = "f(X)",
            title = "Triangular Distribution")

if (save_images) {
    ggsave("Ex2-2_TriangFunc.png", g_triang)
}


# PART 2

# Define the triangular cdf function
ptriang <- function(x, a = 0, b = 1, c = 0.5) {
    m <- rep(0, length(x))
    for (i in seq(1, length(x))){
        if (a <= x[i] && x[i] < c) {
            m[i] <- (x[i] - a)^2 / ((b - a) * (c - a))
        } else if (c <= x[i] && x[i] <= b) {
            m[i] <- (c - a) / (b - a) +
                (2 * b * x[i] - x[i]^2 - 2 * b * c + c^2) / ((b - a) * (b - c))
        }
        # alyernatively we could integrate numerically
        # using m[i] <- integrate(dtriang, lower = a, upper = x[i],
        #                         a = a, b = b, c = c)[[1]]
        }
    return(m)
}

# Define the function that allows to sample random numbers from the
# distribution by inverting the cdf
rtriang <- function(n, a = 0, b = 1, c = 0.5) {
    m <- rep(0, n)
    u <- runif(n)
    # numerically invert the cdf
    inv <- inverse(ptriang, lower = a, upper = b)
    for (i in 1:n) {
        m[i] <- inv(u[i])
    }
    return(m)
}


# PART 3

# Sample 10,000 events from the distribution
rnd_val <- rtriang(10000)
# Build the histogram
g_histtriang <- ggplot() +
                geom_histogram(data = data.frame(rnd_val),
                               mapping = aes(rnd_val, y = after_stat(density)),
                               breaks = seq(-0.05, 1.05, 0.1), color = "black",
                               fill = "lightblue", alpha = 0.5) +
                geom_function(data = data.frame(seq(-0.05, 1.05, 0.01)),
                              fun = dtriang, mapping = aes(color = "red"),
                              linewidth = 1, linetype = 2) +
                scale_color_discrete(labels = c("Expected Distribution")) +
                labs(x = "x", y = "P(x)",
                     title = "Triangular Distribution", color = element_blank())

if (save_images) {
    ggsave("Ex2-2_HistTriang.png", g_histtriang)
}


# ---------------------------------------------------------------------- #


# EXERCISE 3

# The waiting time at the doctor's follows an exponenetial distribution
# with mean
exp_true_mean <- 30 #min
# and rate
rate <- 1 / exp_true_mean

# PART 1

# Simulate the waiting time for 60 people
x <- rexp(60, rate)
# plot the histogram
g_exp <- ggplot() +
         geom_histogram(data = data.frame(x),
                        mapping = aes(x, y = after_stat(density)),
                        color = "black", fill = "lightblue", alpha = 0.5) +
         geom_function(data = data.frame(seq(0, 150, 1)),
                       fun = dexp, args = list(rate = rate),
                       mapping = aes(color = "red"),
                       linewidth = 1, linetype = 2) +
         scale_color_discrete(labels = c("Expected Distribution")) +
         labs(x = "Wait Time [min]", y = "Frequency", color = "",
              title = "Exponential Ditribution")

if (save_images) {
    ggsave("Ex2-3_ExpSim.png", g_exp)
}


# PART 2

# The probability that a person waits for less than 12 minutes is
p_less12 <- pexp(12, rate)
# Alternatively we could integrate the pdf
# p_less12 <- integrate(dexp, lower = 0, upper = 12, rate = 1 / 30)

cat(c("Exercise 3.2:", "\n",
      "The probability that a person waits for less than 12 minutes is",
      p_less12, "\n"))


# PART 3

# The average waiting time from simulated data is
exp_sim_mean <- mean(x)
# The value obtained by manipulating the probability distribution is
exp_calc_mean <- integrate((function(x) x * rate * exp(-x * rate)),
                            lower = 0, upper = Inf)
# The expected value from theory is 1/rate, that is the true_mean we
# defined at the beginning of this ex.
# The relative differences are
exp_diff_simTr <- abs(exp_sim_mean - exp_true_mean) / exp_true_mean
exp_diff_calcTr <- abs(exp_calc_mean[[1]] - exp_true_mean) / exp_true_mean

cat(c("Exercise 3.3:", "\n",
      "The theoretical average is", exp_true_mean, "\n",
      "The average from simulated data is", exp_sim_mean, "\n",
      "Difference of", exp_diff_simTr * 100, "% from theory", "\n",
      "The numerical average by manipulating the pdf is", exp_calc_mean[[1]],
      "\n",
      "Difference of", exp_diff_calcTr * 100, "% from theory", "\n"))


# PART 4

# The probability of waiting more than one hour is
p_more60 <- 1 - pexp(60, rate)
# ALternatively
# p_more60 <- integrate(dexp, lower = 60, upper = Inf, rate = 1 / 30)
# p_more60 <- 1 -  integrate(dexp, lower = 0,  upper = 60, rate = 1 / 30)[[1]]

cat(c("Exercise 3.4:", "\n",
      "The probability that a person waits for more than 1 hour is",
      p_more60, "\n"))


# ---------------------------------------------------------------------- #


# EXERCISE 4


# We will indicate with K (known) the event "the student knows the answer"
# and with NK (not known) the event "the student does not know the answer".
# These events are complementary, so P(NK) = 1 - P(K).
# Similarly, we will use T (true) for the event "the given answer is correct"
# and F (false) for "the given answer is not correct", with P(F) = 1 - P(T).
# From the text we know P(K) and thus P(NK)
p_k <- 0.7
p_nk <- 1 - p_k
# We can also define the conditioned probability of getting a correct/incorrect
# answer in a certain situation. In particular, we can say that there is a 100%
# probability of the answer being correct if the student knows the answer P(T|K)
# (the student has a really good memory and recalls only correct information),
# and a 1/(N. of choices) chance if the student does not know the answer and
# selects randomly from the possible choices P(T|NK). Of course, the parameter
# P(T|K) is not clearly defined in the text, so we could adjust it to express
# how faulty the student's memory is.
p_t_k <- 1
p_t_nk <- 0.2
# We can now exploit Bayes Theorem to find the "reversed" conditioned
# probability of the student knowing the answer when the answer is correct
# using the formula P(K|T) = P(T|K)*P(K)/( P(T|K)*P(K) + P(T|NK)*P(NK) )
p_k_t <- p_t_k * p_k / (p_t_k * p_k + p_t_nk * p_nk)

cat(c("Exercise 4:", "\n",
      "once a correct answer is given, the probability that the student",
      "really knew the correct answer is", p_k_t, "\n"))

# As a bonus, let us see how this probability varies by changing P(K)
x <- seq(0, 1, 0.01)

g_bayes <- ggplot(data.frame(x), aes(x)) +
           geom_function(fun = function(p) {p / (p + 0.2 * (1 - p))}) +
           geom_hline(yintercept = p_k_t, linetype = "dashed", color = "red") +
           geom_vline(xintercept = p_k,   linetype = "dashed", color = "red") +
           geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue") +
           geom_vline(xintercept = 1 / 6, linetype = "dashed", color = "blue") +
           labs(x = "P(K)", y = "P(K|T)",
           title = "Probability of Knowing the Answer given a Correct Answer")

if (save_images) {
    ggsave("Ex2-4_Bayes.png", g_bayes)
}

# The red lines represent the situation we had in the text of the problem,
# while the blue lines are when there is an equal chance that the student
# actually knew the answer or selected randomly. This happens, for an
# arbitrary P(T|K) = a and P(T|NK) = 1/N, for P(K) = 1/(1 + N*a)


# ---------------------------------------------------------------------- #


# EXERCISE 5

# For convenience, let us set 10:30 as t_zero = 0. The probability distribution
# of the arrival time, then, is uniform in the range [10:45, 11:45] -> [15, 75].
# The probability distribution of the waiting time, as a consequence, is also
# uniform in the range [0, 30] (the two extremes when a person arrives right on
# time or as soon as the previous train left). We can check the distribution
# explicitly by defining a "wait" function that gives us the time a person has
# to wait for the next train when arriving at a certain time:

train_time <- c(0, 30, 60, 90)
wait <- function(t, times) {
    w <- Inf
    # find the minumum positive difference between t
    # and the values in the vector
    for (i in seq(1, length(times))) {
        w_i <- times[i] - t
        if (0 <= w_i && w_i < w)
            w <- w_i
    }
    return(w)
}

# and then build a histogram of the values we obtain by extracting randomly a
# starting time in [15, 75]
x <- runif(100000, min = 15, max = 75)
w <- rep(0, length(x))
for (i in seq(1, length(x))){
    w[i] <- wait(x[i], train_time)
}
bin_breaks <- seq(0, 30)
g_wait <- ggplot(data.frame(w), aes(w)) +
          geom_histogram(aes(y = after_stat(density)), breaks = bin_breaks,
                         color = "black", fill = "lightblue", alpha = 0.5) +
          geom_function(fun = function(x) {ifelse(0 <= x & x <= 30, 1 / 30, 0)},
                        aes(color = "red")) +
          scale_color_discrete(labels = c("Expected Distribution")) +
          labs(x = "Wait Time [min]", y = "Probability",
               title = "Wait Time Ditribution", color = element_blank())

if (save_images) {
    ggsave("Ex2-5_HistWait.png", g_wait)
}


# PART 1

# The probability of waiting at most 10 minutes is
p_less10 <- punif(10, min = 0, max = 30)
# or, alternatively, p_less10 <- integrate(dunif, lower = 0, upper = 10,
#                                          min = 0, max = 30)

cat(c("Exercise 5.1:", "\n",
      "The probability of waiting at most 10 minutes is", p_less10, "\n"))


# PART 2

# The probability of waiting at least 15 minutes is
p_more15 <- 1 - punif(15, min = 0, max = 30)
# or, alternatively, p_more15 <- integrate(dunif, lower = 15, upper = 30,
#                                          min = 0, max = 30)

cat(c("Exercise 5.2:", "\n",
      "The probability of waiting at least 15 minutes is", p_more15, "\n"))


# PART 3

# The average time spent waiting is the value in the middle of the interval
wait_true_mean <- (30 - 0) / 2
# We can also check it usingthe simulated values
wait_sim_mean <- mean(w)
# or numerically by using the definition of expected value
# and the wait function
wait_calc_mean <- integrate((function(x) {x * dunif(x, min = 0, max = 30)}),
                             lower = 0, upper = 30)

wait_diff_simTr <- abs(wait_sim_mean - wait_true_mean) / wait_true_mean
wait_diff_calcTr <- abs(wait_calc_mean[[1]] - wait_true_mean) / wait_true_mean

cat(c("Exercise 5.3:", "\n",
      "The theoretical average is", wait_true_mean, "\n",
      "The average from simulated data is", wait_sim_mean, "\n",
      "Difference of", wait_diff_simTr * 100, "% from theory", "\n",
      "The numerical average by manipulating the pdf is", wait_calc_mean[[1]],
      "\n",
      "Difference of", wait_diff_calcTr * 100, "% from theory", "\n"))


# ---------------------------------------------------------------------- #


# EXERCISE 6

# Total investment
tot <- 200 * 85 #euros
# The return rate is a normal variable with
mu <- tot * 0.1
sigma <- tot * 0.12
# The probability that after a year his net profit from the
# investment is at least 800 euros is
p_more800 <- 1 - pnorm(800, mean = mu, sd = sigma)
# or alternatively p_more800 <- integrate(dnorm, lower = 800, upper = Inf,
#                                         mean = mu, sd = sigma)

cat(c("Exercise 6:", "\n",
      "The probability that after a year his net profit",
      "from the investment is at least 800 euros is", p_more800, "\n"))


# Save output to pdf file
if (save_pdf) {
     pdf("Lupi_Enrico_rlab02_Outputs.pdf")
     print(g_bar)
     print(g_k)
     print(g_triang)
     print(g_histtriang)
     print(g_exp)
     print(g_bayes)
     print(g_wait)
     dev.off()
}