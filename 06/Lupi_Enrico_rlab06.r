# Enrico Lupi, 2090596, 4 June 2023


# LABORATORY 6

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(knitr)
library(plotly)
library(rjags)
library(R2jags)
library(tidybayes)
library(runjags)
library(bayestestR)


set.seed(09112001)

# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- TRUE
if (save_pdf) {
    pdf("Lupi_Enrico_rlab06_Outputs.pdf")
}


#  EXERCISE 1

# Define the 1 dimensional Metropolis-Hastings algorithm using an independent
# proposal density (Norm(0; sd))
# Parameters:
# func : a function whose first argument is a real vector of parameters
# func returns a log10 of the likelihood function
# theta_init : the initial value of the Markov Chain (and of func)
# n_sample: number of required samples
# sigma : standar deviation of the gaussian MCMC sampling pdf
metropolis_1d <- function(func, theta_init, n_sample, sigma) {
    theta_curr <- theta_init
    func_curr <- func(theta_curr)
    func_samp <- matrix(data = NA, nrow = n_sample, ncol = 2 + 1)
    n_accept <- 0
    rate_accept <- 0.0
    for (n in 1:n_sample) {
        theta_prop <- rnorm(n = 1, mean = 0, sd = sigma)
        func_prop <- func(theta_prop)
        q_curr <- log10(dnorm(theta_curr, mean = 0, sd = sigma))
        q_prop <- log10(dnorm(theta_prop, mean = 0, sd = sigma))
        # Log10 of the Metropolis ratio
        log_mr <- func_prop - func_curr + q_curr - q_prop
        if (log_mr >= 0 || log_mr > log10(runif(1))) {
            theta_curr <- theta_prop
            func_curr <- func_prop
            n_accept <- n_accept + 1
            rate_accept <- n_accept / n
        }
        func_samp[n, 1] <- func_curr
        func_samp[n, 2] <- theta_curr
        func_samp[n, 3] <- rate_accept
    }
return(func_samp)
}

# posterior to sample from
gfunc <- function(theta) {
    res <- 0.5 * exp(-1 * (theta + 3)**2 / 2) +
           0.5 * exp(-1 * (theta - 3)**2 / 2)
    return(res)
}

# posterior integral
# used only for visualisation and checking purposes
integral <- integrate(gfunc, -Inf, Inf)

# loglikelihood required by metropolis function
logpost <- function(theta) {
    return(log10(gfunc(theta)))
}

# markov Chain initial values
theta_init <- 0
# Markov Chain step size
sample_cov <- 1
# burn-in values
burnin_val <- c(0, 100, 500)
# thinning values
thinning_val <- c(1, 10, 20, 30)
#samples
n_samp <- 1e5

# Let us now produce various chain with different burn-in and thinning settings
index <- 1
for (b in burnin_val) {
    for (t in thinning_val) {
        samples <- metropolis_1d(func = logpost, theta_init = theta_init,
                                 n_sample = b + n_samp, sigma = sample_cov)

        # burn=in and thinning
        selection <- seq(from = b + 1, to = b + n_samp, by = t)
        post_samp <- samples[selection, 2]

        # compute autocorrelation function
        chain <- as.mcmc(post_samp)
        lags <- seq(0, 500, 10)
        acf <- autocorr(chain, lags = lags)

        # plot chain
        g_chain <- ggplot() +
            geom_line(aes(x = selection, y = post_samp),
                      color = "black", size = 0.1) +
            labs(x = "Iteration", y = "Theta", title = "")

        # plot theta distribution
        post_den <- density(post_samp, n = 2^10)
        g_theta <- ggplot() +
            geom_line(mapping = aes(x = post_den$x, y = post_den$y,
                                    color = "Monte Carlo")) +
            geom_line(mapping = aes(x = post_den$x,
                                    y = gfunc(post_den$x) / integral$value,
                                    color = "Theoretical Distribution"),
                      alpha = 0.5) +
            geom_line(mapping = aes(x = post_den$x,
                                    y = dnorm(post_den$x, mean = 0,
                                                          sd = sample_cov),
                                    color = "Proposal Distribution")) +
            scale_color_manual(values = c("black", "grey", "red")) +
            labs(x = "Theta", y = "Density", color = "", title = "")

        # plot autocorrelation
        g_corr <- ggplot() +
            geom_point(aes(x = lags, y = acf),
                       color = "black") +
            geom_line(aes(x = lags, y = acf),
                      color = "black", alpha = 0.8) +
            labs(x = "Lag", y = "ACF", title = "") +
            annotate(geom = "text", x = 200, y = 0.8,
                     label = paste("Effective Size = ",
                                   round(effectiveSize(chain), digits = 3)))

        info <- paste("Ex. 1 : Burn-in = ", b, "  Thinning = ", t)

        g_g <- ggarrange(g_chain, g_theta, g_corr,
                         labels = c("", "", ""),
                         ncol = 1, nrow = 3)
        g_g <- annotate_figure(g_g, top = text_grob(info))
        if (save_pdf) {
            print(g_g)
        }

        index <- index + 1
    }
}

# Let us now analyse the results.
# First of all, by looking at the chain itself and the density plot we can see
# that the proposed independent distribution is not fit for the posterior
# under study: it has a peak at 0 and falls off rather quickly for higher
# (absolute) values, while the posterior has a minimum at 0 and maxima at +-3.
# This means that a lot of proposed thetas will be rejected and the chain has
# a tendency of getting "stuck". A simple solution would be to increase the
# sigma of the proposal distribution.
# Given the poor quality of the chain, it is difficult to analyse the results
# of the autocorrelation function as they heavily depend on random effects.
# In general, we can say that the burn-in does not significantly influence the
# results (as expected, since the proposal distribution is independent), while
# the thinning factor helps in quickly reducing the value of ACF while keeping
# the effective size of the chain acceptable.


# ---------------------------------------------------------------------- #


# EXERCISE 2

# function to plot results
plotRes <- function(res, var, ci, color, bf = 0) {
    res <- as.data.frame(res)

    # prepare credibility interval plot
    delta <- 0.0001
    ndigits <- 3
    if (var == "diff_rate") {
        delta <- 0.01
        ndigits <- 1
    }
    x <- seq(ci[1], ci[2], delta)
    xmin  <- res[, var] |> min()
    xmax  <- res[, var] |> max()
    perc    <- ((ci[1] + ci[2]) / 2 - xmin) / (xmax - xmin)
    perc_sx <- (ci[1]               - xmin) / (xmax - xmin)
    perc_dx <- (ci[2]               - xmin) / (xmax - xmin)
    label_ci <- textGrob(label = "C.I. 95%", x = perc, y = 0.08,
                    just = c("center", "bottom"),
                    gp = gpar(col = "black", size = 5))
    label_sx <- textGrob(label = round(ci[1], digits = ndigits),
                         x = perc_sx, y = 0.05,
                         just = c("right", "bottom"),
                         gp = gpar(col = "black", size = 4))
    label_dx <- textGrob(label = round(ci[2], digits = ndigits),
                         x = perc_dx, y = 0.05,
                         just = c("left", "bottom"),
                         gp = gpar(col = "black", size = 4))

    # prepare mean plot
    xmean <- res[, var] |> mean()
    perc_mean <- (xmean - xmin) / (xmax - xmin)
    label_mean <- textGrob(label = paste("Mean = ",
                                         round(xmean, digits = ndigits)),
                    x = perc_mean, y = 0.7,
                    just = c("left", "bottom"),
                    gp = gpar(col = "black", size = 3))

    # plot histogram of samples
    g <- ggplot() +
        geom_histogram(mapping = aes(x = res[, var], y = ..density..),
                       bins = 100, color = color,
                       fill = color, alpha = 0.5) +
        geom_line(mapping = aes(x = x, y = rep(0, length(x))),
                  color = "black", linewidth = 2) +
        geom_vline(xintercept = xmean, color = "black",
                   linetype = "dashed", alpha = 0.5) +
        annotation_custom(label_ci, xmin = -Inf, xmax = Inf,
                                    ymin = -Inf, ymax = Inf) +
        annotation_custom(label_sx, xmin = -Inf, xmax = Inf,
                                    ymin = -Inf, ymax = Inf) +
        annotation_custom(label_dx, xmin = -Inf, xmax = Inf,
                                    ymin = -Inf, ymax = Inf) +
        annotation_custom(label_mean, xmin = -Inf, xmax = Inf,
                                      ymin = -Inf, ymax = Inf)

    # add Bayes Factor
    if (var == "diff_rate") {
        label_bf <- textGrob(label = paste("Effectiveness > 50%",
                                            "\nlog( BF ) = ", bf),
                    x = 0.1, y = 0.7,
                    just = c("left", "bottom"),
                    gp = gpar(col = "black", size = 4))
        g <- g +
            annotation_custom(label_bf, xmin = -Inf, xmax = Inf,
                                        ymin = -Inf, ymax = Inf) +
            labs(x = "Diff Rate", y = "Density", title = "")
    } else {
        g <- g +
            labs(x = "Infection Rate", y = "Density", title = var)
    }

    return(g)
}

# We will analyse the efficacy of the following vaccines
vaccines <- c("Jcovden (ex Janssen)", "SpikeVax (ex. Moderna)",
              "Vaxzevria (ex AstraZeneca)")

# The studies have been conducted on two groups of people, one was
# administered the vaccine while the other just a placebo
tot_vaccine <- c(19630, 14134, 5258)
tot_placebo <- c(19691, 14073, 5210)
# Number of patients tested postive after RCT:
pos_vaccine <- c(116, 11, 64)
pos_placebo <- c(348, 185, 154)

# JAGS model
model_string <- "
    model {
        for ( i in 1:Ntot ) {
            tested[i] ~ dbern( theta[patient[i]] )
        }
        for ( k in 1:Nclass ) {
            theta[k] ~ dbeta(3 , 100)
        }
    }"

for (i in 1:length(vaccines)) {
    # organise data
    patient <- c(rep("Vaccine", tot_vaccine[i]), rep("Placebo", tot_placebo[i]))
    tested <- c(rep("Pos", pos_vaccine[i]),
                rep("Neg", tot_vaccine[i] - pos_vaccine[i]),
                rep("Pos", pos_placebo[i]),
                rep("Neg", tot_placebo[i] - pos_placebo[i]))
    vaccine_tb <- tibble(tested = tested, patient = patient)
    table(vaccine_tb[[2]], vaccine_tb[[1]])

    # setup data list for JAGS
    data_list <- list(tested = ifelse(vaccine_tb$tested == "Neg", 0, 1),
                      patient = as.integer(factor(vaccine_tb$patient)),
                      Ntot = nrow(vaccine_tb),
                      Nclass = nlevels(factor(vaccine_tb$patient)))

    # run JAGS model
    chain <- run.jags(model_string,
                      sample = 15000,
                      n.chains = 4,
                      method = "parallel",
                      monitor = "theta",
                      data = data_list)

    # from the MC chain, we can compute the efficacy
    res <- tidybayes::tidy_draws(chain) %>%
           select("theta[1]":"theta[2]") %>%
           rename(Placebo = "theta[1]", Vaccine = "theta[2]") %>%
           mutate(diff_rate = (Placebo - Vaccine) / Placebo * 100,
                  Placebo_perc = Placebo * 100, Vaccine_perc = Vaccine * 100)

    mcmc <- as.mcmc(res, vars = "diff_rate")

    # compute 95% credibility interval
    summary <- summary(mcmc)
    ci_95 <- list(c(summary[[2]]["Placebo",    "2.5%"],
                    summary[[2]]["Placebo",   "97.5%"]),
                  c(summary[[2]]["Vaccine",    "2.5%"],
                    summary[[2]]["Vaccine",   "97.5%"]),
                  c(summary[[2]]["diff_rate",  "2.5%"],
                    summary[[2]]["diff_rate", "97.5%"]))

    # assuming as a prior a Normal distribution with mean = 50% and sd = 15%
    # compute the odds that with our prior and posterior the Vaccine is more
    # that 50% effective
    prior  <- bayestestR::distribution_normal(60000, mean = 50, sd = 15)
    bf <- bayestestR::bayesfactor_parameters(res$diff_rate, prior,
                            direction = "two-sided", null = 50)[[1]] |>
          round(digits = 2)

    # plot results
    g_placebo <- plotRes(res, "Placebo",   ci_95[[1]], "lightpink")
    g_vaccine <- plotRes(res, "Vaccine",   ci_95[[2]], "lightblue")
    g_rate    <- plotRes(res, "diff_rate", ci_95[[3]], "yellow", bf)

    g_pv <- ggarrange(g_placebo, g_vaccine,
                      labels = c("", ""),
                      ncol = 2, nrow = 1)
    g_pvdr <- ggarrange(g_pv, g_rate,
                        labels = c("", ""),
                        ncol = 1, nrow = 2)
    g_pvdr <- annotate_figure(g_pvdr,
                              top = text_grob(paste("Ex 2 : ", vaccines[i])))
    if (save_pdf) {
        print(g_pvdr)
    }
}


# ---------------------------------------------------------------------- #


# EXERCISE 3

# load dataset
owid_df <- read_csv("06/owid-covid-data.csv")


# Vaccinated People :

# Cumulative :
# we will use "new_people_vaccinated_smoothed", the daily number of people
# receiving their first vaccine dose (7-day smoothed), so as to find the
# cumulative number of people that have received at least one dose
cumvax <- owid_df |>
    select(date, continent,
           new_people_vaccinated_smoothed) |>
    mutate(vax = new_people_vaccinated_smoothed) |>
    group_by(date) |>
    summarize(tot_vaxxed    = sum(new_people_vaccinated_smoothed,
                                  na.rm = TRUE),
              tot_vaxxed_eu = sum(new_people_vaccinated_smoothed
                                  [continent == "Europe"],
                                  na.rm = TRUE)) |>
    mutate(cum_vaxxed    = cumsum(tot_vaxxed),
           cum_vaxxed_eu = cumsum(tot_vaxxed_eu))

g_cumworld <- cumvax |> ggplot() +
    geom_line(aes(x = date, y = cum_vaxxed), color = "red") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "3 months",
                 limits = as.Date(c("2020-12-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 2.1e10),
                       breaks = seq(0, 2.1e10, by = 0.25e10)) +
    labs(x = "Date", y = "Total Vaccinated People", title = "World")

g_cumeu <- cumvax |> ggplot() +
    geom_line(aes(x = date, y = cum_vaxxed_eu), color = "blue") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "3 months",
                 limits = as.Date(c("2020-12-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 6e8),
                       breaks = seq(0, 6e8, by = 1e8)) +
    labs(x = "Date", y = "Total Vaccinated People", title = "Europe")

g_cum <- ggarrange(g_cumworld, g_cumeu,
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
g_cum <- annotate_figure(g_cum,
            top = text_grob("Ex 3 : Cumulative Sum of Vaccinated People"))
if (save_pdf) {
    print(g_cum)
}


# Daily number :
# we will use "new_vaccinations", the new COVID-19 vaccination doses
# administered, so as to find the daily amount of doses administered
# regardless of wether it is the first dose or not
dailyvax <- owid_df |>
    select(date, continent, new_vaccinations) |>
    mutate(vax = new_vaccinations) |>
    group_by(date) |>
    summarize(daily    = sum(new_vaccinations,
                             na.rm = TRUE),
              daily_eu = sum(new_vaccinations[continent == "Europe"],
                             na.rm = TRUE))

g_dailyworld <- dailyvax |> ggplot() +
    geom_line(aes(x = date, y = daily), color = "red") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "3 months",
                 limits = as.Date(c("2020-12-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 1.8e8),
                       breaks = seq(0, 1.8e8, by = 0.2e8)) +
    labs(x = "Date", y = "New Daily Doses", title = "World", color = "")

g_dailyeu <- dailyvax |> ggplot() +
    geom_line(aes(x = date, y = daily_eu), color = "blue") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "3 months",
                 limits = as.Date(c("2020-12-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 7e6),
                       breaks = seq(0, 7e6, by = 1e6)) +
    labs(x = "Date", y = "New Daily Doses", title = "Europe", color = "")

g_daily <- ggarrange(g_dailyworld, g_dailyeu,
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
g_daily <- annotate_figure(g_daily,
                top = text_grob("Ex 3 : Daily Number of Vaccinated People"))
if (save_pdf) {
    print(g_daily)
}


# Weekly average :
# we will use "new_vaccinations" once again
weeklyvax <- owid_df |>
    select(date,  continent, new_vaccinations) |>
    mutate(vax = new_vaccinations) |>
    mutate(date = strftime(date, format = "%Y-%W")) |>
    group_by(date) |>
    summarise(weekly    = sum(new_vaccinations,
                              na.rm = TRUE) / 7,
              weekly_eu = sum(new_vaccinations[continent == "Europe"],
                              na.rm = TRUE) / 7) |>
    mutate(date = as_date(paste(date, "-1", sep = ""), format = "%Y-%W-%u"))

g_weeklyworld <- weeklyvax |> ggplot() +
    geom_line(aes(x = date, y = weekly), color = "red") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "3 months",
                 limits = as.Date(c("2020-12-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 1.8e8),
                       breaks = seq(0, 1.8e8, by = 0.2e8)) +
    labs(x = "Date", y = "New Doses (Weekly Average)",
         title = "World", color = "")

g_weeklyeu <- weeklyvax |> ggplot() +
    geom_line(aes(x = date, y = weekly_eu), color = "blue") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "3 months",
                 limits = as.Date(c("2020-12-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 7e6),
                       breaks = seq(0, 7e6, by = 1e6)) +
    labs(x = "Date", y = "New Doses (Weekly Average)",
         title = "Europe", color = "")

g_weekly <- ggarrange(g_weeklyworld, g_weeklyeu,
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
g_weekly <- annotate_figure(g_weekly,
                top = text_grob("Ex 3 : Weekly Average of Vaccinated People"))
if (save_pdf) {
    print(g_weekly)
}


# Number of Deaths :

# Cumulative :
cumdeath <- owid_df |>
    select(date, continent, new_deaths) |>
    group_by(date) |>
    summarize(tot_deaths    = sum(new_deaths, na.rm = TRUE),
              tot_deaths_eu = sum(new_deaths[continent == "Europe"],
                                  na.rm = TRUE)) |>
    mutate(cum_deaths    = cumsum(tot_deaths),
           cum_deaths_eu = cumsum(tot_deaths_eu))

g_deathworld <- cumdeath |> ggplot() +
    geom_line(aes(x = date, y = cum_deaths), color = "red") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "4 months",
                 limits = as.Date(c("2020-3-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 3e7),
                       breaks = seq(0, 3e7, by = 0.5e7)) +
    labs(x = "Date", y = "Total Deaths", title = "World")

g_deatheu <- cumdeath |> ggplot() +
    geom_line(aes(x = date, y = cum_deaths_eu), color = "blue") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "4 months",
                 limits = as.Date(c("2020-3-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 2.2e6),
                       breaks = seq(0, 2e6, by = 0.25e6)) +
    labs(x = "Date", y = "Total Deaths", title = "Europe")

g_deaths <- ggarrange(g_deathworld, g_deatheu,
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
g_deaths <- annotate_figure(g_deaths,
                top = text_grob("Ex 3 : Cumulative Sum of Deaths"))
if (save_pdf) {
    print(g_deaths)
}


# Weekly average :
weeklydeath <- owid_df |>
    select(date,  continent, new_deaths) |>
    mutate(date = strftime(date, format = "%Y-%W")) |>
    group_by(date) |>
    summarise(weekly    = sum(new_deaths,
                              na.rm = TRUE) / 7,
              weekly_eu = sum(new_deaths[continent == "Europe"],
                              na.rm = TRUE) / 7) |>
    mutate(date = as_date(paste(date, "-1", sep = ""), format = "%Y-%W-%u"))

g_weeklydeathsworld <- weeklydeath |> ggplot() +
    geom_line(aes(x = date, y = weekly), color = "red") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "4 months",
                 limits = as.Date(c("2020-3-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 7e4),
                       breaks = seq(0, 7e4, by = 0.5e4)) +
    labs(x = "Date", y = "New Deaths (Weekly Average)",
         title = "World", color = "")

g_weeklydeathseu <- weeklydeath |> ggplot() +
    geom_line(aes(x = date, y = weekly_eu), color = "blue") +
    scale_x_date(date_labels = "%b %Y",
                 date_breaks = "4 months",
                 limits = as.Date(c("2020-3-1", "2023-7-1"))) +
    scale_y_continuous(limits = c(0, 6e3),
                       breaks = seq(0, 6e3, by = 0.5e3)) +
    labs(x = "Date", y = "New Deaths (Weekly Average)",
         title = "Europe", color = "")

g_weeklydeaths <- ggarrange(g_weeklydeathsworld, g_weeklydeathseu,
                   labels = c("", ""),
                   ncol = 1, nrow = 2)
g_weeklydeaths <- annotate_figure(g_weeklydeaths,
                top = text_grob("Ex 3 : Weekly Average of Deaths"))
if (save_pdf) {
    print(g_weeklydeaths)
}


if (save_pdf) {
    dev.off()
}
