# Enrico Lupi, 2090596, 16 Apr 2023

# Auxiliary variables to save plots
save_images <- FALSE
save_pdf <- TRUE
# Avoid warning messages
options(warn = -1)

# EXERCISE 1.1

library(tidyverse)
library(ggpubr)
#library(gridExtra)


# PART 1 : read the data and import them in a data.frame or tibble structure

aa_empl <- read_delim("Data/american_airline_empl.txt", delim = "\t",
                      locale = locale(grouping_mark = ",", decimal_mark = "."),
                      col_types = "cc???")
da_empl <- read_delim("Data/delta_airline_empl.txt",    delim = "\t",
                      locale = locale(grouping_mark = ",", decimal_mark = "."),
                      col_types = "cc???")
fe_empl <- read_delim("Data/federal_express_empl.txt",  delim = "\t",
                      locale = locale(grouping_mark = ",", decimal_mark = "."),
                      col_types = "cc???")
ua_empl <- read_delim("Data/united_airline_empl.txt",   delim = "\t",
                      locale = locale(grouping_mark = ",", decimal_mark = "."),
                      col_types = "cc???")

# Add a column for the date
aa_empl <- aa_empl |> mutate(Date = as.Date(paste(Year, Month, "1", sep = "-")))
da_empl <- da_empl |> mutate(Date = as.Date(paste(Year, Month, "1", sep = "-")))
fe_empl <- fe_empl |> mutate(Date = as.Date(paste(Year, Month, "1", sep = "-")))
ua_empl <- ua_empl |> mutate(Date = as.Date(paste(Year, Month, "1", sep = "-")))


# PART 2 : merge the four data tibble in a common tibble

tot_empl <- bind_rows(aa_empl, da_empl,
                      fe_empl, ua_empl)

# Add a new column with the name of the corresponding company
airlines <- c(rep("American Airlines", times = nrow(aa_empl)),
              rep("Delta Airlines",    times = nrow(da_empl)),
              rep("Federal Express",   times = nrow(fe_empl)),
              rep("United Airlines",   times = nrow(ua_empl)))

tot_empl <- tot_empl |> mutate(Company = airlines)


# PART 3 : produce a plot of the behaviour of the employees
#          as a function of time for all four companies,
#          separately for the number of full-time and part-time employees

g_parttime <- ggplot(data = tot_empl) +
              geom_line(aes(x = Date, y = tot_empl$'Part-time',
                        colour = Company)) +
              labs(x = "Date", y = "Number of Employees",
                   title = "Part time", colour = "Company")


g_fulltime <- ggplot(data = tot_empl) +
              geom_line(aes(x = Date, y = tot_empl$'Full-time',
                        colour = Company)) +
              labs(x = "Date", y = "Number of Employees",
                   title = "Full time", colour = "Company")

g_empl <- ggarrange(g_parttime, g_fulltime,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)


# PART 4 : when did each company reach the minimum and
#	       maximum number of employees?

aa_max <-  aa_empl[which.max(aa_empl$'Grand Total'), "Date"] |> pull()
aa_min <-  aa_empl[which.min(aa_empl$'Grand Total'), "Date"] |> pull()
da_max <-  da_empl[which.max(da_empl$'Grand Total'), "Date"] |> pull()
da_min <-  da_empl[which.min(da_empl$'Grand Total'), "Date"] |> pull()
fe_max <-  fe_empl[which.max(fe_empl$'Grand Total'), "Date"] |> pull()
fe_min <-  fe_empl[which.min(fe_empl$'Grand Total'), "Date"] |> pull()
ua_max <-  ua_empl[which.max(ua_empl$'Grand Total'), "Date"] |> pull()
ua_min <-  ua_empl[which.min(ua_empl$'Grand Total'), "Date"] |> pull()


# PART 5 : plot the fraction of part-time worker over
#	       the total employess as a function of time

tot_empl <- tot_empl |> mutate(Ratio = tot_empl$'Part-time' / tot_empl$'Grand Total')

g_ratio <- ggplot(data = tot_empl) +
           geom_line(aes(x = Date, y = Ratio, colour = Company)) +
           labs(x = "Date", y = "Ratio",
                title = "Part Time Employees over Total", colour = "Company")


# PART 6 : did the COVID-19 pandemic have any influence
#		   in the employed workers of the airline companies?
#		   Can you see a trend in the years 2019-2023?

tot_empl_2019 <- tot_empl |> filter(Year >= 2019)
g_total_2019    <- ggplot(data = tot_empl_2019) +
                   geom_line(aes(x = Date, y = tot_empl_2019$'Grand Total',
                             colour = Company)) +
                   labs(x = "Date", y = "Number of Employees",
                   title = "Total (after 2019)", colour = "Company")

g_parttime_2019 <- ggplot(data = tot_empl_2019) +
                   geom_line(aes(x = Date, y = tot_empl_2019$'Part-time',
                             colour = Company)) +
                   labs(x = "Date", y = "Number of Employees",
                   title = "Part time (after 2019)", colour = "Company")


g_fulltime_2019 <- ggplot(data = tot_empl_2019) +
                   geom_line(aes(x = Date, y = tot_empl_2019$'Full-time',
                             colour = Company)) +
                   labs(x = "Date", y = "Number of Employees",
                   title = "Full time (after 2019)", colour = "Company")

g_empl_2019 <- ggarrange(g_total_2019, g_parttime_2019, g_fulltime_2019,
                         labels = c("A", "B", "C"),
                         ncol = 1, nrow = 3)


# We can see that there was a sharp decline in the number of full- time
# workers for Delta Airlines at the start of 2020, right as the pandemic
# started, while for American and United it happened later on, in the second
# half of the year; on the other hand, it does not seem that Federal
# Express was hit as the others as, afeter a very small initial dip, the
# number of full-time employees rose even higher than the pre-pandemic period.
# The number of part-time employees shows similar trends for United and Delta
# Airlines, while it differs for the other companies: American Airlines does
# not seem to be influenced by the pandemic, while Federal Express has a more
# oscillating trend that resulted in a net negative, although in a longer
# period of time so it may be not directly linked to the pandemic.


# Save plots and print outputs
if (save_images) {
    ggsave("Ex1-1_Employement.png", g_empl)
    ggsave("Ex1-1_Ratio.png", g_ratio)
    ggsave("Ex1-1_Employement_2019.png", g_empl_2019)
}
print("Ex 1.1.4 - Number of Employees:")
print(paste("American Airlines    max:", aa_max, "  min:", aa_min))
print(paste("Delta Airlines       max:", da_max, "  min:", da_min))
print(paste("Federal Express      max:", fe_max, "  min:", fe_min))
print(paste("United Airlines      max:", ua_max, "  min:", ua_min))



#--------------------------------------------------------------------------#



# EXERCISE 1.2

library(nycflights13)

# PART 1.1 : Plot the total number of flights departed from each
#            of the three NYC airports as a function of time
#            (one entry for each of the 365 days of the year).

# Keep the total amount of flights for each day and airport
f_per_day <- flights |> transmute(Date = as_date(time_hour), origin, flight) |>
                        group_by(origin, Date) |>
                        summarise(tot = n_distinct(flight))

g_fpd <- ggplot(data = f_per_day) +
         geom_line(aes(x = Date, y = tot, colour = origin)) +
         labs(x = "Date", y = "Flights per Day",  colour = "Airport")


# PART 1.2 : Plot the average number of flights computed over the first
#            five working days of each week as a function of the week number
#            of the year. Produce the same plot for the flights departing
#            over the weekend (Saturdays and Sundays).

# add a column with week number (starting with the first Monday as the first
# day of week 01) and day of the week (with Sunday as day 1)
f_weekly <- f_per_day |> mutate(Week = as.integer(strftime(Date, format = "%W")),
                                Day = wday(Date)) |>
                         group_by(Week, weekend = (Day == 1 | Day == 7)) |>
                         summarise(avg = mean(tot))

g_fpw <- ggplot(data = f_weekly) +
         geom_line(aes(x = Week, y = avg, colour = weekend)) +
         scale_color_discrete(labels = c("Weekday", "Weekend")) +
         labs(x = "Number of the Week", y = "Average Number of Flights",
              colour = "")


# PART 2.1 : For each flight in the data frame, compute the departure delay
#            and extract the following pieces of information (separately
#            for each NYC airport):
#            - min, max and average delay for each day of the year
#             (show the data in corresponding plots)

# "%j" retruns the day of the year as a decimal number (range 001 to 366).
f_delay <- flights |> transmute(Day = as.integer(strftime(time_hour, format = "%j")),
                                origin, dep_delay) |>
                      group_by(origin, Day) |>
                      summarise(min = min(dep_delay,  na.rm = TRUE),
                                max = max(dep_delay,  na.rm = TRUE),
                                avg = mean(dep_delay, na.rm = TRUE))

g_min <- ggplot(data = f_delay) +
         geom_line(aes(x = Day, y = min, colour = origin)) +
         labs(x = "Day of the Year", y = "Min Dep. Delay [min]",
              title = "", colour = "Airport")

g_max <- ggplot(data = f_delay) +
         geom_line(aes(x = Day, y = max, colour = origin)) +
         labs(x = "Day of the Year", y = "Max Dep. Delay [min]",
              title = "", colour = "Airport")

g_avg <- ggplot(data = f_delay) +
         geom_line(aes(x = Day, y = avg, colour = origin)) +
         labs(x = "Day of the Year", y = "Average Dep. Delay [min]",
              title = "", colour = "Airport")

g_delay <- ggarrange(g_min, g_max, g_avg,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3)


# PART 3 : assuming the distance flew by the plane is, at first approximation,
#          the distance between the two connecting airports (as given in
#          the data frame), compute the average speed of each plane. Produce
#          a plot of the average plane speed as a function of departure
#          day of the year

f_speed <- flights |> transmute(Day = as.integer(strftime(time_hour, format = "%j")),
                                speed = 60 * distance / air_time) |>
                      group_by(Day) |>
                      summarise(avg_speed = mean(speed, na.rm = TRUE))

g_speed <- ggplot(data = f_speed) +
           geom_line(aes(x = Day, y = avg_speed)) +
           labs(x = "Day of the Year", y = "Average Speed [mph]", title = "")


# PART 4 : analyze the flights offered by each airline company and determine:
#          - the airline companies offering the largest two numbers of
#            flights per day and per week;
#          - the airline company offering the smallest number of flights per month
#          - the airline company offering the longest distance flight per month.
#          (you can produce plots, if you like, to visualize the results
#           of the analysis)

f_data <- flights |> transmute(Day  = as.integer(strftime(time_hour, format = "%j")),
                               Week = as.integer(strftime(time_hour, format = "%W")),
                               month, distance, carrier, flight)

company_max_fpw <- f_data |> group_by(Week, carrier) |>
                             summarise(tot = n_distinct(flight)) |>
                             filter(tot == sort(tot, decreasing = TRUE)[1] |
                                    tot == sort(tot, decreasing = TRUE)[2])

company_max_fpd <- f_data |> group_by(Day, carrier) |>
                             summarise(tot = n_distinct(flight)) |>
                             filter(tot == sort(tot, decreasing = TRUE)[1] |
                                    tot == sort(tot, decreasing = TRUE)[2])

company_min_fpm <- f_data |> group_by(month, carrier) |>
                             summarise(tot = n_distinct(flight)) |>
                             filter(tot == min(tot))

company_max_dpm <- f_data |> group_by(month) |>
                             filter(distance == max(distance)) |>
                             transmute(carrier, distance) |>
                             distinct() |>
                             arrange(month)


# Save plots and print outputs
if (save_images) {
    ggsave("Ex1-2_FlightsPerDay.png", g_fpd)
    ggsave("Ex1-2_FlightsWeekly.png", g_fpw)
    ggsave("Ex1-2_FlightsDelay.png", g_delay)
    ggsave("Ex1-2_FlightsSpeed.png", g_speed)
}

print("Ex 1.2.4:")
print("Ailines offering argest two numbers of flights per day:")
print(company_max_fpd, n = 20)
print("...and per week:")
print(company_max_fpw, n = 20)
print("Airlines offering the smallest number of flights per month:")
print(company_min_fpm)
print("Airlines offering the longest distance flight per month:")
print(company_max_dpm)

# Save output to pdf file
if (save_pdf) {
     pdf("Lupi_Enrico_rlab01_Outputs.pdf")
     print(g_empl)
     print(g_ratio)
     print(g_empl_2019)
     print(g_fpd)
     print(g_fpw)
     print(g_delay)
     print(g_speed)
     dev.off()
}