rm(list=ls())

require(tidyverse)
require(readxl)
require(lubridate)
library(zoo)
require(MMWRweek)
library(RColorBrewer)

setwd("~/GitHub/thresholds_leptospirosis")


################################################################################
################################################################################
# ISLANDWIDE: load confirmed/probable case data 1999 - present
################################################################################
################################################################################

#previous dataset is 1999-2023
ds_1999_2018 <- left_join(data.frame(year=rep(1999:2018,each=52),
                                     week=rep(1:52,length(1999:2018))) %>%
                            mutate(weekdate = MMWRweek2Date(MMWRyear = year,MMWRweek = week)),
                          read.csv("~/Lepto/data/PRThresholdModelling_MasterLineList_2024March.csv")%>%
                            rename(week = Investigation_MMWR_Week,
                                   year = Investigation_MMWR_Year) %>%
                            filter(
                              week != 53,
                              calc_case_classification != "Suspected") %>%
                            group_by(year,week) %>%
                            summarize(cases = n()) %>%
                            ungroup() %>%
                            mutate(weekdate = MMWRweek2Date(MMWRyear = year,MMWRweek = week))) %>%
  replace(is.na(.), 0)

#newest dataset is 2019 - present
ds_2019_2025 <- full_join((data.frame(year=rep(2019:2025,each=52),
                            week=rep(1:52,length(2019:2025))) %>%
                   mutate(weekdate = MMWRweek2Date(MMWRyear = year,MMWRweek = week))),
                (read.csv("~/Lepto/data/LeptoCases_2019-20250715.csv")%>%
                   rename(week = Investigation.MMWR.Week,
                          year = Investigation.MMWR.Year) %>%
                   filter(
                     week != 53,
                     Investigation.Case.Status  != "Suspected") %>%
                   group_by(year,week) %>%
                   summarize(cases = n()) %>%
                   ungroup() %>%
                   mutate(weekdate = MMWRweek2Date(MMWRyear = year,MMWRweek = week)))) %>%
  replace(is.na(.), 0) %>%
  mutate(
    cases = replace(cases, 341:364, NA) #replace missing / future data with NA
  )

# combine datasets
df_all <- bind_rows(ds_1999_2018,ds_2019_2025)%>%
  mutate(ind2022on = ifelse(year>=2022,1,0)) #add indicator for pre/post 2022

# Note dates for Hurricanes Maria and Fiona
fiona.dt <- as.Date("2022-09-17")
irma.dt <- as.Date("2017-09-06")
maria.dt <- as.Date("2017-09-20")

# plot reported cases
ggplot() +
  geom_col(data = df_all, aes(x = weekdate, y = cases, fill = "Reported cases"), na.rm = TRUE) +
  scale_fill_manual(name = NULL, values = c("Reported cases" = "grey")) +
  # scale_x_date(date_breaks = "3 months", labels = scales::date_format("%b-%Y")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.key = element_blank()
  ) +
  labs(x = "Onset week", y = "Cases")


################################################################################
################################################################################
# ISLANDWIDE: Conduct time-series bootstrapping on case data
################################################################################
################################################################################

#temporal bootstrapping
sample_window = 5
n_samples = 5

#reformat
df_bs <- cbind(df_all,
               data.frame(rollapply( df_all$cases , width = sample_window , FUN=sample, 
                                     size=n_samples, replace=T, fill = NA, align = "center" )))
  
#merge with original dataset
df <- full_join(df_all,
                (df_bs %>%
                   pivot_longer(
                     cols = starts_with("X"),
                     names_to = "sim",
                     names_prefix = "X",
                     values_to = "bootstrap_val",
                     values_drop_na = T
                   ))) %>%
  mutate(week_factor = as.factor(week))


################################################################################
################################################################################
# ISLANDWIDE: calculate historical thresholds 2013-present
# allow for change in reporting definition in 2022
################################################################################
################################################################################

threshold_evol <- tibble()
for (yr in 2015:2025) {
  for (wk in 1:52) {
    this_nb <- MASS:::glm.nb(bootstrap_val ~ week_factor*ind2022on,
                             data=filter(df, weekdate < MMWRweek2Date(MMWRyear = yr-1,MMWRweek = wk)))
    if (wk == 1){
      threshold_evol <- bind_rows(threshold_evol,
                                  tibble(
                                    year = yr,
                                    week = wk,
                                    q50 = qnbinom(0.5, size=this_nb$theta, 
                                                  mu=exp(this_nb$coefficients[[1]])), # median
                                    q75 = qnbinom(0.75, size=this_nb$theta, 
                                                  mu=exp(this_nb$coefficients[[1]])),
                                    q80 = qnbinom(0.80, size=this_nb$theta, 
                                                  mu=exp(this_nb$coefficients[[1]])),
                                    q85 = qnbinom(0.85, size=this_nb$theta, 
                                                  mu=exp(this_nb$coefficients[[1]])),
                                    q90 = qnbinom(0.90, size=this_nb$theta, 
                                                  mu=exp(this_nb$coefficients[[1]])),
                                    q95 = qnbinom(0.95, size=this_nb$theta, 
                                                  mu=exp(this_nb$coefficients[[1]]))
                                  ))
    } else {
      threshold_evol <- bind_rows(threshold_evol,
                                  tibble(
                                    year = yr,
                                    week = wk,
                                    q50 = qnbinom(0.5, size=this_nb$theta, 
                                                  mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]],na.rm=T))), # median
                                    q95 = qnbinom(0.95, size=this_nb$theta, 
                                                  mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]],na.rm=T))),
                                    q90 = qnbinom(0.90, size=this_nb$theta, 
                                                  mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]],na.rm=T))),
                                    q85 = qnbinom(0.85, size=this_nb$theta, 
                                                  mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]],na.rm=T))),
                                    q80 = qnbinom(0.80, size=this_nb$theta, 
                                                  mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]],na.rm=T))),
                                    q75 = qnbinom(0.75, size=this_nb$theta, 
                                                  mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]],na.rm=T)))
                                
                                  ))
    }
    
    print(c(yr,wk))
  }
}


# Join datasets and plot
conf_thresh <- df_all %>% 
  left_join(threshold_evol)

# Apply smoothing filter to all thresholds
smooth.filt <- c(1, 2, 3, 4, 5, 4, 3, 2, 1)/25

conf_thresh$q50[829:832] <- (conf_thresh$q50)[1193:1196] #for smoothing purposes (to avoid missing vals), wrap last observations of 2021 to end of 2014 
conf_thresh$q75[829:832] <- (conf_thresh$q75)[1193:1196]
conf_thresh$q80[829:832] <- (conf_thresh$q80)[1193:1196]
conf_thresh$q85[829:832] <- (conf_thresh$q85)[1193:1196]
conf_thresh$q90[829:832] <- (conf_thresh$q90)[1193:1196]
conf_thresh$q95[829:832] <- (conf_thresh$q95)[1193:1196]
N = nrow(conf_thresh)

conf_thresh$q50.sm <- (stats::filter(c(conf_thresh$q50, conf_thresh$q50[(which(conf_thresh$year==2025))[1:4]]), smooth.filt, circular=T))[-c((N+1):(N+4))]
conf_thresh$q75.sm <- (stats::filter(c(conf_thresh$q75, conf_thresh$q75[(which(conf_thresh$year==2025))[1:4]]), smooth.filt, circular=T))[-c((N+1):(N+4))]
conf_thresh$q80.sm <- (stats::filter(c(conf_thresh$q80, conf_thresh$q80[(which(conf_thresh$year==2025))[1:4]]), smooth.filt, circular=T))[-c((N+1):(N+4))]
conf_thresh$q85.sm <- (stats::filter(c(conf_thresh$q85, conf_thresh$q85[(which(conf_thresh$year==2025))[1:4]]), smooth.filt, circular=T))[-c((N+1):(N+4))]
conf_thresh$q90.sm <- (stats::filter(c(conf_thresh$q90, conf_thresh$q90[(which(conf_thresh$year==2025))[1:4]]), smooth.filt, circular=T))[-c((N+1):(N+4))]
conf_thresh$q95.sm <- (stats::filter(c(conf_thresh$q95, conf_thresh$q95[(which(conf_thresh$year==2025))[1:4]]), smooth.filt, circular=T))[-c((N+1):(N+4))]


########################################
# plot all options
conf_thresh_options <- conf_thresh %>%
  pivot_longer(cols = c(q75.sm, q80.sm, q85.sm, q90.sm, q95.sm),
               names_to = "type", values_to = "value") %>%
  mutate(type = recode(type,
                       q75.sm = "Real-time 75% outbreak threshold",
                       q80.sm = "Real-time 80% outbreak threshold",
                       q85.sm = "Real-time 85% outbreak threshold",
                       q90.sm = "Real-time 90% outbreak threshold",
                       q95.sm = "Real-time 95% outbreak threshold"))

# Create a new column for plotting categories
conf_thresh$plot_type <- "Reported cases"
conf_thresh_options$plot_type <- conf_thresh_options$type

# Combine all plot labels into one factor to ensure order
plot_levels <- c("Reported cases", "Hurricane Fiona", "Hurricane Maria",
                 "Real-time 75% outbreak threshold", "Real-time 80% outbreak threshold", "Real-time 85% outbreak threshold", "Real-time 90% outbreak threshold", "Real-time 95% outbreak threshold")

hurricane_dates <- data.frame(
  hurricane = c("Hurricane Fiona", "Hurricane Maria"),
  date = as.Date(c(fiona.dt, maria.dt))
)

ggplot() +
  
  # Reported cases bars
  geom_col(data = conf_thresh, aes(x = weekdate, y = cases, fill = "Reported cases"), color = NA, na.rm = TRUE) +
  
  # Vertical lines as segments (y from bottom to top of plot)
  geom_segment(data = hurricane_dates,
               aes(x = date, xend = date, y = 0, yend = Inf, color = hurricane, linetype = hurricane),
               linewidth = 1.1, show.legend = TRUE) +
  
  # Real-time lines
  geom_line(data = conf_thresh_options,
            aes(x = weekdate, y = value, color = type, linetype = type),
            linewidth = 1, show.legend = TRUE) +
  
  # Manual fills
  scale_fill_manual(name = NULL, values = c("Reported cases" = "grey")) +
  
  # Manual colors (all)
  scale_color_manual(name = NULL, values = c(
    # "Reported cases" = "grey",
    "Hurricane Fiona" = "seagreen",
    "Hurricane Maria" = "limegreen",
    "Real-time 75% outbreak threshold" = "purple",
    "Real-time 80% outbreak threshold" = "hotpink",
    "Real-time 85% outbreak threshold" = "firebrick",
    "Real-time 90% outbreak threshold" = "orangered",
    "Real-time 95% outbreak threshold" = "orange"
  )) +
  
  # Manual linetypes (all)
  scale_linetype_manual(name = NULL, values = c(
    "Hurricane Fiona" = "dotted",
    "Hurricane Maria" = "dotted",
    "Real-time 75% outbreak threshold" = "solid",
    "Real-time 80% outbreak threshold" = "solid",
    "Real-time 85% outbreak threshold" = "solid",
    "Real-time 90% outbreak threshold" = "solid",
    "Real-time 95% outbreak threshold" = "solid"
  )) +
  
  # Legends: override to get horizontal lines with no points
  guides(
    fill = guide_legend(order = 1, override.aes = list(
      fill = "grey", color = NA, linetype = "blank", shape = NA)),
    color = guide_legend(order = 2, override.aes = list(
      linetype = c("dotted", "dotted", "solid", "solid", "solid", "solid", "solid"),
      color = c("seagreen", "limegreen", "purple", "hotpink","firebrick","orangered","orange"),
      linewidth = 1.5,
      shape = NA
    )),
    linetype = "none"
  ) +
  
  scale_x_date(date_breaks = "3 months", labels = scales::date_format("%b-%Y")) +
  scale_y_continuous(limits = c(0, NA)) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.key = element_blank()
  ) +
  labs(x = "Onset week", y = "Cases")

#####################################################
# full time series with smoothed 90% thresholds

conf_thresh_realtime <- conf_thresh %>%
  pivot_longer(cols = c(q50.sm, q90.sm),
               names_to = "type", values_to = "value") %>%
  mutate(type = recode(type,
                       q50.sm = "Real-time historical median",
                       q90.sm = "Real-time outbreak threshold"))

# Create a new column for plotting categories
conf_thresh$plot_type <- "Reported cases"
conf_thresh_realtime$plot_type <- conf_thresh_realtime$type

# Combine all plot labels into one factor to ensure order
plot_levels <- c("Reported cases", "Hurricane Fiona", "Hurricane Maria",
                 "Real-time historical median", "Real-time outbreak threshold")

hurricane_dates <- data.frame(
  hurricane = c("Hurricane Fiona", "Hurricane Maria"),
  date = as.Date(c(fiona.dt, maria.dt))
)

ggplot() +
  
  # Reported cases bars
  geom_col(data = conf_thresh, aes(x = weekdate, y = cases, fill = "Reported cases"), color = NA, na.rm = TRUE) +
  
  # Vertical lines as segments (y from bottom to top of plot)
  geom_segment(data = hurricane_dates,
               aes(x = date, xend = date, y = 0, yend = Inf, color = hurricane, linetype = hurricane),
               linewidth = 1.1, show.legend = TRUE) +
  
  # Real-time lines
  geom_line(data = conf_thresh_realtime,
            aes(x = weekdate, y = value, color = type, linetype = type),
            linewidth = 1, show.legend = TRUE) +
  
  # Manual fills
  scale_fill_manual(name = NULL, values = c("Reported cases" = "grey")) +
  
  # Manual colors (all)
  scale_color_manual(name = NULL, values = c(
    # "Reported cases" = "grey",
    "Hurricane Fiona" = "seagreen",
    "Hurricane Maria" = "limegreen",
    "Real-time historical median" = "steelblue",
    "Real-time outbreak threshold" = "firebrick"
  )) +
  
  # Manual linetypes (all)
  scale_linetype_manual(name = NULL, values = c(
    "Hurricane Fiona" = "dotted",
    "Hurricane Maria" = "dotted",
    "Real-time historical median" = "solid",
    "Real-time outbreak threshold" = "solid"
  )) +
  
  # Legends: override to get horizontal lines with no points
  guides(
    fill = guide_legend(order = 1, override.aes = list(
      fill = "grey", color = NA, linetype = "blank", shape = NA)),
    color = guide_legend(order = 2, override.aes = list(
      linetype = c("dotted", "dotted", "solid", "solid"),
      color = c("seagreen", "limegreen", "steelblue", "firebrick"),
      linewidth = 1.5,
      shape = NA
    )),
    linetype = "none"
  ) +
  
  scale_x_date(date_breaks = "3 months", labels = scales::date_format("%b-%Y")) +
  scale_y_continuous(limits = c(0, NA)) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.key = element_blank()
  ) +
  labs(x = "Onset week", y = "Cases")

##########################################
# 2022-present smoothed with retrospective 2025 threshold (what we would use now looking back)

conf_thresh$thresh2025.sm <- NA
conf_thresh$med2025.sm <- NA

foo <- conf_thresh %>%
  filter(conf_thresh$year == 2025)

#apply most recent year's threshold to 2022 onward (what we would now consider the threshold)
conf_thresh$thresh2025.sm[which(conf_thresh$year %in% c(2022,2023, 2024, 2025))] <- foo$q90.sm
conf_thresh$med2025.sm[which(conf_thresh$year %in% c(2022,2023, 2024, 2025))] <- foo$q50.sm

conf_thresh_long <- conf_thresh %>%
  pivot_longer(cols = c(med2025.sm, thresh2025.sm),
               names_to = "type", values_to = "value") %>%
  mutate(type = recode(type,
                       med2025.sm = "Retrospective 2025 median",
                       thresh2025.sm = "Retrospective 2025 threshold"
  ))

legend_data <- tibble(
  weekdate = as.Date("2022-01-01"),
  value = 0,  # irrelevant
  type = factor(c("Retrospective 2025 median", "Retrospective 2025 threshold", "Hurricane Fiona"))
)


ggplot() +
  # Reported cases as grey bars
  geom_col(data = conf_thresh, aes(x = weekdate, y = cases, fill = "Reported cases"), color = "white", na.rm = TRUE) +
  
  # Actual threshold lines
  geom_line(data = conf_thresh_long,
            aes(x = weekdate, y = value, color = type, linetype = type),
            linewidth = 1) +
  
  # Actual vertical vline for Fiona (no legend)
  geom_vline(xintercept = as.Date(fiona.dt), linetype = "dotted", color = "black", linewidth = 1.1) +
  
  # Fake horizontal line for Hurricane Fiona to get it in the legend
  geom_line(data = legend_data %>% filter(type == "Hurricane Fiona"),
            aes(x = weekdate, y = value, color = type, linetype = type),
            linewidth = 1) +
  
  # Manual fill for bars
  scale_fill_manual(name = NULL, values = c("Reported cases" = "grey")) +
  
  # Manual color and linetype for lines
  scale_color_manual(name = NULL, values = c(
    "Retrospective 2025 median" = "steelblue",
    "Retrospective 2025 threshold" = "firebrick",
    "Hurricane Fiona" = "black"
  )) +
  scale_linetype_manual(name = NULL, values = c(
    "Retrospective 2025 median" = "solid",
    "Retrospective 2025 threshold" = "solid",
    "Hurricane Fiona" = "dotted"
  )) +
  
  # Legend appearance overrides
  guides(
    fill = guide_legend(order = 1, override.aes = list(color = NA)),
    color = guide_legend(order = 2, override.aes = list(
      shape = NA,
      linetype = c("solid", "solid", "dotted"),
      color = c("steelblue", "firebrick", "black"),
      linewidth = 1.2
    )),
    linetype = "none"  # suppress separate linetype legend
  ) +
  
  # Axis and theme
  scale_x_date(date_breaks = "1 month",
               labels = scales::date_format("%b-%Y"),
               limits = as.Date(c('2022-01-01','2025-12-31'))) +
  scale_y_continuous(limits = c(0, 20)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.key = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = "Onset week", y = "Cases")




##############################################################################
##############################################################################
# ISLANDWIDE: Sensitivity analysis with different levels of thresholds and consecutive week rules
##############################################################################
##############################################################################
df_thresh <- conf_thresh #save in case we need unmodified version later

conf_thresh <- df_thresh %>%
  select(weekdate,year,week,cases,q50.sm,q75.sm,q80.sm,q85.sm,q90.sm,q95.sm) %>%
  pivot_longer(
    cols = starts_with("q"),
    names_to = "percentile",
    values_to = "thresh_realt",
    values_drop_na = F
  ) %>%
  mutate(percentile = as.numeric(as.character(fct_recode(percentile,
                                 "50" = "q50.sm",
                                 "75" = "q75.sm",
                                 "80" = "q80.sm",
                                 "85" = "q85.sm",
                                 "90" = "q90.sm",
                                 "95" = "q95.sm"))),
         cases_nonmissing = cases, 
         cases_nonmissing = replace_na(cases_nonmissing, 0)) 

conf_thresh_tmp <- conf_thresh %>%
   filter(conf_thresh$year == 2025)

excess.table.retro <- tibble()
excess.table.realt <- tibble()
for (wks in 2:4) {
  outbr.wks = wks #apply one consecutive week rule at a time
  
  for (perc in seq(75,95,5)) {
    df_tmp <- conf_thresh %>% #filter to one threshold level
      filter(percentile == perc)
    df_tmp_retro <- conf_thresh %>% #filter to current year
      filter(conf_thresh$year == 2025,
             percentile == perc)
    df_tmp$thresh_retro <- NA
    df_tmp$thresh_retro[which(df_tmp$year %in% c(2022:2025))] <- df_tmp_retro$thresh_realt #retrospective thresholds 
    
    ###########################################################
    ##### calculate outbreaks with retrospective thresholds 2022-2025
    ###########################################################
    df_tmp$abovethresh <- (df_tmp$cases_nonmissing>df_tmp$thresh_retro) #indicator for whether cases are above retrospective threshold
    # count consecutive weeks above threshold
    df_tmp$count <-sequence(rle(as.character(df_tmp$abovethresh))$lengths)
    df_tmp$count_cum <- with(rle(df_tmp$abovethresh), rep(values * cumsum(values & lengths >= outbr.wks ),lengths))
    
    init_y <- 1999 #starting year
    thresh_init_y <- 2022 #first year of thresholds we are concerned with
    curr_y <- 2025 #last observed year in dataset (ok if incomplete but can't be nonmissing - fill 0s)
    
    df_tmp$count_y <- NA
    df_tmp$count_cum_y <- NA
    for (y in 1:length(table(df_tmp$year))){
      df_tmp$count_y[which(df_tmp$year==(y+(init_y - 1)))] <-sequence(rle(as.character(df_tmp$abovethresh[which(df_tmp$year==(y+(init_y - 1)))]))$lengths)
      df_tmp$count_cum_y[which(df_tmp$year==(y+(init_y - 1)))] <- with(rle(df_tmp$abovethresh[which(df_tmp$year==(y+(init_y - 1)))]), rep(values * cumsum(values & lengths >= outbr.wks ),lengths))
    }
    
    outbrk_y <- list()
    # for (y in (thresh_init_y ):((init_y - 1)+length(table(df_tmp$year)))){
    for (y in (thresh_init_y ):(curr_y)){
      outbrk_wk_tbl <- (table(df_tmp$count_cum_y[which(df_tmp$year==y)]))
      if (dimnames(outbrk_wk_tbl)[[1]][1]=="0") {outbrk_wk_tbl <- outbrk_wk_tbl[-1]}
      outbrk_num <- length(outbrk_wk_tbl)
      
      num_wks <- NULL
      if (outbrk_num>1) {
        for(i in 1:(outbrk_num)){
          num_wks <- c(num_wks,max(df_tmp$count_y[which((df_tmp$year==y) & 
                                                               df_tmp$count_cum_y == i)]))
        }} 
      if (outbrk_num == 1) {
        num_wks <- c(num_wks,max(df_tmp$count_y[which((df_tmp$year==y) & 
                                                             df_tmp$count_cum_y == 1)]))
      }
      if (outbrk_num == 0){
        num_wks <- c(num_wks,0)
        
      } 
      
      outbrk_y[[y-(thresh_init_y - 1)]] <- num_wks
      
      # print(c(y,num_wks))
    }
    
    ds_y <- data.frame(year=(thresh_init_y ):curr_y)
    ds_y$num_outbrk_wks <- NA
    ds_y$num_outbrks <- NA
    for (i in 1:length((thresh_init_y ):curr_y)){
      ds_y$num_outbrk_wks[i] <- sum(outbrk_y[[i]])
      ds_y$num_outbrks[i] <- length(outbrk_y[[i]])
    }
    
    ds_y$num_outbrks[which(ds_y$num_outbrk_wks==0)] <- NA
    r <- rle(df_tmp$abovethresh)
    r$values <- r$lengths >= outbr.wks & r$values 
    df_tmp$outbrk_flag <- inverse.rle(r)
    
    # calculate excess weeks and cases
    df_tmp$ex.all <- (df_tmp$cases_nonmissing - df_tmp$thresh_retro)*(df_tmp$outbrk_flag)
    df_tmp$ex.all <- round((df_tmp$ex.all), 0)
    df_tmp$ex.all[which(df_tmp$ex.all<0)] <- 0
    df_tmp$ex.all[which(is.na(df_tmp$ex.all))] <- 0
    
    excess.table_retro_tmp <- df_tmp %>%
      group_by(year) %>%
      summarise(no.consec.epidemic.weeks = sum(outbrk_flag),
                excess.all = sum(ex.all)) %>%
      full_join(ds_y) %>%
      mutate(no.consec.epidemic.weeks = if_else(no.consec.epidemic.weeks < outbr.wks, 0, no.consec.epidemic.weeks), 
             threshold_perc = perc,
             weekrule = wks,
             threshold_label = paste0(perc,"th %ile")) %>%
      filter(year >= 2022) %>%
      select (weekrule,threshold_perc,threshold_label,year,num_outbrk_wks,num_outbrks)
    
    
    excess.table.retro <- bind_rows(excess.table.retro,excess.table_retro_tmp)
    
    
    ###########################################################
    ##### calculate outbreaks with realtime thresholds pre-2022
    ###########################################################
    df_tmp <- conf_thresh %>% #filter to one threshold level
      filter(percentile == perc)
    
    df_tmp$abovethresh <- (df_tmp$cases_nonmissing>df_tmp$thresh_realt) #indicator for whether cases are above realtime threshold
    # count consecutive weeks above threshold
    df_tmp$count <-sequence(rle(as.character(df_tmp$abovethresh))$lengths)
    df_tmp$count_cum <- with(rle(df_tmp$abovethresh), rep(values * cumsum(values & lengths >= outbr.wks ),lengths))
    
    init_y <- 1999 #starting year
    thresh_init_y <- 2015 #first year of thresholds we are concerned with
    curr_y <- 2025 #last observed year in dataset (ok if incomplete but can't be nonmissing - fill 0s)
    
    df_tmp$count_y <- NA
    df_tmp$count_cum_y <- NA
    for (y in 1:length(table(df_tmp$year))){
      df_tmp$count_y[which(df_tmp$year==(y+(init_y - 1)))] <-sequence(rle(as.character(df_tmp$abovethresh[which(df_tmp$year==(y+(init_y - 1)))]))$lengths)
      df_tmp$count_cum_y[which(df_tmp$year==(y+(init_y - 1)))] <- with(rle(df_tmp$abovethresh[which(df_tmp$year==(y+(init_y - 1)))]), rep(values * cumsum(values & lengths >= outbr.wks ),lengths))
    }
    
    outbrk_y <- list()
    for (y in (thresh_init_y ):(curr_y)){
      outbrk_wk_tbl <- (table(df_tmp$count_cum_y[which(df_tmp$year==y)]))
      if (dimnames(outbrk_wk_tbl)[[1]][1]=="0") {outbrk_wk_tbl <- outbrk_wk_tbl[-1]}
      outbrk_num <- length(outbrk_wk_tbl)
      
      num_wks <- NULL
      if (outbrk_num>1) {
        for(i in 1:(outbrk_num)){
          num_wks <- c(num_wks,max(df_tmp$count_y[which((df_tmp$year==y) & 
                                                          df_tmp$count_cum_y == i)]))
        }} 
      if (outbrk_num == 1) {
        num_wks <- c(num_wks,max(df_tmp$count_y[which((df_tmp$year==y) & 
                                                        df_tmp$count_cum_y == 1)]))
      }
      if (outbrk_num == 0){
        num_wks <- c(num_wks,0)
        
      } 
      
      outbrk_y[[y-(thresh_init_y - 1)]] <- num_wks
    }
    
    ds_y <- data.frame(year=(thresh_init_y ):curr_y)
    ds_y$num_outbrk_wks <- NA
    ds_y$num_outbrks <- NA
    for (i in 1:length((thresh_init_y ):curr_y)){
      ds_y$num_outbrk_wks[i] <- sum(outbrk_y[[i]])
      ds_y$num_outbrks[i] <- length(outbrk_y[[i]])
    }
    
    ds_y$num_outbrks[which(ds_y$num_outbrk_wks==0)] <- NA
    r <- rle(df_tmp$abovethresh)
    r$values <- r$lengths >= outbr.wks & r$values 
    df_tmp$outbrk_flag <- inverse.rle(r)
    
    # calculate excess weeks and cases
    df_tmp$ex.all <- (df_tmp$cases_nonmissing - df_tmp$thresh_realt)*(df_tmp$outbrk_flag)
    df_tmp$ex.all <- round((df_tmp$ex.all), 0)
    df_tmp$ex.all[which(df_tmp$ex.all<0)] <- 0
    df_tmp$ex.all[which(is.na(df_tmp$ex.all))] <- 0
    
    excess.table_realt_tmp <- df_tmp %>%
      group_by(year) %>%
      summarise(no.consec.epidemic.weeks = sum(outbrk_flag),
                excess.all = sum(ex.all)) %>%
      full_join(ds_y) %>%
      mutate(no.consec.epidemic.weeks = if_else(no.consec.epidemic.weeks < outbr.wks, 0, no.consec.epidemic.weeks), 
             threshold_perc = perc,
             weekrule = wks,
             threshold_label = paste0(perc,"th %ile")) %>%
      filter(year < 2022, year >= 2015) %>%
      select(weekrule,threshold_perc,threshold_label,year,num_outbrk_wks,num_outbrks) 
    
    
    excess.table.realt <- bind_rows(excess.table.realt,excess.table_realt_tmp)
    
    print(c(wks,perc))
  }
}

excess.table.retro <- excess.table.retro %>%
  filter(threshold_perc != 50) %>%
  mutate(threshold_label = factor(threshold_label, levels = c("75th %ile", "80th %ile", "85th %ile", "90th %ile", "95th %ile")),
         weekrule = factor(weekrule))

excess.table.realt <- excess.table.realt %>%
  filter(threshold_perc != 50) %>%
  mutate(threshold_label = factor(threshold_label, levels = c("75th %ile", "80th %ile", "85th %ile", "90th %ile", "95th %ile")),
         weekrule = factor(weekrule))

# Plot retrospective threshold levels and week rules
ggplot(excess.table.retro, aes(x = factor(year), y = num_outbrk_wks, fill = threshold_label)) +
  geom_col(position = "dodge") +
  facet_wrap(
    ~ weekrule,
    ncol = 1,
    labeller = labeller(
      weekrule = c(
        "2" = "2+ consecutive weeks",
        "3" = "3+ consecutive weeks",
        "4" = "4+ consecutive weeks"
      )
    )
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Threshold") +
  labs(
    x = "Year",
    y = "Consecutive outbreak weeks",
    title = "Annual outbreak duration by threshold definition: \nRetrospective thresholds 2022-present"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot real-time threshold levels and week rules
ggplot(excess.table.realt, aes(x = factor(year), y = num_outbrk_wks, fill = threshold_label)) +
  geom_col(position = "dodge") +
  facet_wrap(
    ~ weekrule,
    ncol = 1,
    labeller = labeller(
      weekrule = c(
        "2" = "2+ consecutive weeks",
        "3" = "3+ consecutive weeks",
        "4" = "4+ consecutive weeks"
      )
    )
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Threshold") +
  labs(
    x = "Year",
    y = "Consecutive outbreak weeks",
    title = "Annual outbreak duration by threshold definition: \nReal-time thresholds 2015-2021"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
