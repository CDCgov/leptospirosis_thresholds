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
# HEALTH REGION: load data, time series bootstrapping, thresholds
################################################################################
################################################################################

hreg <- c("Aguadilla Region","Arecibo Region","Bayamon Region","Caguas Region","Fajardo Region","Mayaguez Region","Metro Region","Ponce Region")


#used 90th percentile and counted each week above thresh here, but can change
conf_thresh_list <- tibble()
excess.table.retro <- tibble()
for (j in 1:length(hreg)){
  ds <- full_join((data.frame(year=rep(2019:2025,each=52),
                             week=rep(1:52,length(2019:2025))) %>%
                    mutate(weekdate = MMWRweek2Date(MMWRyear = year,MMWRweek = week))),
                  (read.csv("~/Lepto/data/LeptoCases_2019-20250715.csv")%>%
                     filter(Epi.Health.Region == hreg[j]) %>% #only selected health region
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

  ##############s
  
  #temporal bootstrapping
  sample_window = 5
  n_samples = 5
  df_bs <- cbind(ds,
                 (data.frame(rollapply( ds$cases , width = sample_window , FUN=sample, 
                                        size=n_samples, replace=T, fill = NA, align = "center" ))))
  
  df <- full_join(ds,
                  (df_bs %>%
                     pivot_longer(
                       cols = starts_with("X"),
                       names_to = "sim",
                       names_prefix = "X",
                       values_to = "bootstrap_val",
                       values_drop_na = T
                     ))) %>%
    mutate(week_factor = as.factor(week),
           ind2022on = ifelse(year>=2022,1,0))
  
  
  ### calculate thresholds in the usual way
  threshold_evol <- tibble()
  for (yr in 2022:2025) {
    for (wk in 1:52) {
      this_nb <- MASS:::glm.nb(bootstrap_val ~ week_factor*ind2022on,
                               data=filter(df, weekdate < MMWRweek2Date(MMWRyear = yr-1,MMWRweek = wk)),init.theta = 3)
      
      
      if (wk == 1){
        threshold_evol <- bind_rows(threshold_evol,
                                    tibble(
                                      year = yr,
                                      week = wk,
                                      q50 = qnbinom(0.5, size=this_nb$theta, 
                                                    mu=exp(this_nb$coefficients[[1]])), # median
                                      q90 = qnbinom(0.90, size=this_nb$theta, 
                                                    mu=exp(this_nb$coefficients[[1]]))
                                    ))
      } else {
        threshold_evol <- bind_rows(threshold_evol,
                                    tibble(
                                      year = yr,
                                      week = wk,
                                      q50 = qnbinom(0.5, size=this_nb$theta, 
                                                    mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]], this_nb$coefficients[[52+wk]],na.rm=T))), # median
                                      q90 = qnbinom(0.90, size=this_nb$theta, 
                                                    mu=exp(sum(this_nb$coefficients[[1]], this_nb$coefficients[[wk]], this_nb$coefficients[[53]], this_nb$coefficients[[52+wk]],na.rm=T)))
                                      
                                    ))
      }
    }
  }
  
  
  conf_thresh <- ds %>% 
    left_join(threshold_evol)  %>%
    mutate(cases_nonmissing = cases, 
           cases_nonmissing = replace_na(cases_nonmissing, 0)) 
  
  
  #smooth thresholds and median
  filt_length <- 11
  smooth.filt <-c(1:filt_length, (filt_length-1):1)/sum(c(1:filt_length, (filt_length-1):1))
  
  conf_thresh$q50[147:156] <- (conf_thresh$q50)[355:364] #for smoothing purposes (to avoid missing vals), wrap last observations of 2025 to end of 2021 
  conf_thresh$q90[147:156] <- (conf_thresh$q90)[355:364]
  N = nrow(conf_thresh)
  
  conf_thresh$q50.sm <- (stats::filter(c(conf_thresh$q50, conf_thresh$q50[(which(conf_thresh$year==2025))[1:11]]), smooth.filt, circular=T))[-c((N+1):(N+11))]
   conf_thresh$q90.sm <- (stats::filter(c(conf_thresh$q90, conf_thresh$q90[(which(conf_thresh$year==2025))[1:11]]), smooth.filt, circular=T))[-c((N+1):(N+11))]
  
  
  ##########################################################
  ######## COUNT WEEKS ABOVE THRESHOLD#############
   # RETROSPECTIVE
  
  conf_thresh$thresh_retro <- NA
  conf_thresh$med_retro <- NA
  conf_thresh_tmp <- conf_thresh %>%
    filter(conf_thresh$year == 2025)
  conf_thresh$thresh_retro[which(conf_thresh$year %in% c(2022:2025))] <- conf_thresh_tmp$q90.sm
  conf_thresh$med_retro[which(conf_thresh$year %in% c(2022:2025))] <- conf_thresh_tmp$q50.sm
  
  
  
  outbr.wks = 1 #3
  
  conf_thresh$abovethresh <- (conf_thresh$cases_nonmissing>conf_thresh$thresh_retro)
  conf_thresh$count <-sequence(rle(as.character(conf_thresh$abovethresh))$lengths)
  conf_thresh$count_cum <- with(rle(conf_thresh$abovethresh), rep(values * cumsum(values & lengths >= outbr.wks ),lengths))
  
  init_y <- 2019 #starting year
  thresh_init_y <- 2022 #first year of thresholds
  curr_y <- 2025 #last observed year in dataset (ok if incomplete)
  
  conf_thresh$count_y <- NA
  conf_thresh$count_cum_y <- NA
  for (y in 1:length(table(conf_thresh$year))){
    conf_thresh$count_y[which(conf_thresh$year==(y+(init_y - 1)))] <-sequence(rle(as.character(conf_thresh$abovethresh[which(conf_thresh$year==(y+(init_y - 1)))]))$lengths)
    conf_thresh$count_cum_y[which(conf_thresh$year==(y+(init_y - 1)))] <- with(rle(conf_thresh$abovethresh[which(conf_thresh$year==(y+(init_y - 1)))]), rep(values * cumsum(values & lengths >= outbr.wks ),lengths))
  }
  
  outbrk_y <- list()
  # for (y in (thresh_init_y ):((init_y - 1)+length(table(conf_thresh$year)))){
  for (y in (thresh_init_y ):(curr_y)){
    outbrk_wk_tbl <- (table(conf_thresh$count_cum_y[which(conf_thresh$year==y)]))
    if (dimnames(outbrk_wk_tbl)[[1]][1]=="0") {outbrk_wk_tbl <- outbrk_wk_tbl[-1]}
    outbrk_num <- length(outbrk_wk_tbl)
    
    num_wks <- NULL
    if (outbrk_num>1) {
      for(i in 1:(outbrk_num)){
        num_wks <- c(num_wks,max(conf_thresh$count_y[which((conf_thresh$year==y) & 
                                                             conf_thresh$count_cum_y == i)]))
      }} 
    if (outbrk_num == 1) {
      # if (names(outbrk_wk_tbl)=="1"){
      num_wks <- c(num_wks,max(conf_thresh$count_y[which((conf_thresh$year==y) & 
                                                           conf_thresh$count_cum_y == 1)]))
    }#}
    if (outbrk_num == 0){
      num_wks <- c(num_wks,0)
      
    } 
    
    outbrk_y[[y-(thresh_init_y - 1)]] <- num_wks
    
  }
  
  
  ds_y <- data.frame(year=(thresh_init_y ):curr_y)
  ds_y$num_outbrk_wks <- NA
  ds_y$num_outbrks <- NA
  ds_y$same_outbrk_y_b4 <- 0
  ds_y$same_outbrk_y_after <- 0
  for (i in 1:length((thresh_init_y ):curr_y)){
    ds_y$num_outbrk_wks[i] <- sum(outbrk_y[[i]])
    ds_y$num_outbrks[i] <- length(outbrk_y[[i]])
    if (i<length((thresh_init_y ):curr_y)){
      if ((conf_thresh$count_y[which((conf_thresh$year==(i+(thresh_init_y ))) & 
                                     (conf_thresh$week==(52)))] > 0) &
          (conf_thresh$count_y[which((conf_thresh$year==(i+(thresh_init_y ))) & 
                                     (conf_thresh$week==(1)))] > 0) & 
          (conf_thresh$count_cum_y[which((conf_thresh$year==(i+(thresh_init_y -1))) & 
                                         (conf_thresh$week==(52)))] > 0 ) & 
          (conf_thresh$count_cum_y[which((conf_thresh$year==(i+(thresh_init_y -1 ))) & 
                                         (conf_thresh$week==(1)))] > 0 )) {
        ds_y$same_outbrk_y_b4[i] <- ds_y$same_outbrk_y_after[i+1] <- 1
      }}
  }
  
  ds_y$num_outbrks[which(ds_y$num_outbrk_wks==0)] <- NA
  r <- rle(conf_thresh$abovethresh)
  r$values <- r$lengths >= outbr.wks & r$values 
  conf_thresh$outbrk_flag <- inverse.rle(r)
  
  # calculate excess weeks and cases
  conf_thresh$ex.all <- (conf_thresh$cases_nonmissing - conf_thresh$thresh_retro)*(conf_thresh$outbrk_flag) 
  
  conf_thresh$ex.all <- round((conf_thresh$ex.all), 0)
  conf_thresh$ex.all[which(conf_thresh$ex.all<0)] <- 0
  conf_thresh$ex.all[which(is.na(conf_thresh$ex.all))] <- 0

  excess.table_retro_tmp <- df_tmp %>%
    group_by(year) %>%
    summarise(no.consec.epidemic.weeks = sum(outbrk_flag),
              excess.all = sum(ex.all)) %>%
    full_join(ds_y) %>%
    mutate(no.consec.epidemic.weeks = if_else(no.consec.epidemic.weeks < outbr.wks, 0, no.consec.epidemic.weeks), 
           threshold_perc = perc,
           threshold_label = paste0(perc,"th %ile"),
           hreg = hreg[j]) %>%
    filter(year >= 2022) %>%
    select (hreg,threshold_perc,threshold_label,year,num_outbrk_wks,num_outbrks)
  
  
  excess.table.retro <- bind_rows(excess.table.retro,excess.table_retro_tmp)
  
  
  conf_thresh$health_region <- hreg[j]
  
  conf_thresh_list <-  bind_rows(conf_thresh_list,
                                    tibble(conf_thresh))
  

  print(j)
  
}

#################################

# Plot each health region from 2022-present with median and retrospective thresholds
legend_data <- tibble(
  weekdate = as.Date("2022-01-01"),
  value = 0,  # irrelevant
  type = factor(c(#"Real-time historical median", "Real-time outbreak threshold", 
                  "Retrospective 2025 median", "Retrospective 2025 threshold", "Hurricane Fiona"), levels = c(#"Real-time historical median", "Real-time outbreak threshold", 
                    "Retrospective 2025 median", "Retrospective 2025 threshold", "Hurricane Fiona"))
)




ggplot() +
  # Reported cases as grey bars
  geom_col(
    data = conf_thresh_list,
    aes(x = weekdate, y = cases, fill = "Reported cases"),
    # color = "white",      # Optional: outlines for clarity
    width = 7,            # Wider bars for week-based x-axis
    na.rm = TRUE
  ) +
  
  # Threshold and median lines
  geom_line(data = conf_thresh_long,
            aes(x = weekdate, y = value, color = type, linetype = type),
            size = 1) +
  
  # Actual vertical line for Hurricane Fiona (not in legend)
  geom_vline(xintercept = as.Date(fiona.dt),
             linetype = "dotted", color = "black", size = 1.1) +
  
  # Dummy line for Fiona to include in legend
  geom_line(data = legend_data %>% filter(type == "Hurricane Fiona"),
            aes(x = weekdate, y = value, color = type, linetype = type),
            size = 1) +
  
  # Facet by health region
  facet_wrap(~health_region, ncol = 2, scales = "free_y") +
  
  # Manual colors and linetypes
  scale_fill_manual(name = NULL, values = c("Reported cases" = "grey")) +
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
  
  # Legend formatting
  guides(
    fill = guide_legend(order = 1, override.aes = list(color = NA)),
    color = guide_legend(order = 2, override.aes = list(
      shape = NA,
      linetype = c("solid", "solid", "dotted"),
      color = c("steelblue", "firebrick", "black"),
      size = 1.2
    )),
    linetype = "none"
  ) +
  
  # Axes and theme
  scale_x_date(date_breaks = "1 month",
               labels = scales::date_format("%b-%Y"),
               limits = as.Date(c('2022-01-01', '2025-12-31'))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Onset week", y = "Cases") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    legend.key = element_blank(),
    legend.position = "bottom"
  )


# plot # outbreak weeks by year and health region, sorted from overall highest to lowest
outbrk_df <- excess.table.retro %>%
  group_by(hreg) %>%
  mutate(total_outbreak_weeks = sum(num_outbrk_wks, na.rm = TRUE),
         health_region = factor(hreg)  ) %>%
  ungroup() %>%
  mutate(health_region = reorder(health_region, -total_outbreak_weeks))

ggplot(outbrk_df, aes(x = factor(year), y = num_outbrk_wks, fill = health_region)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Dark2", name = "Health Region") +
  labs(
    x = "Year",
    y = "Total outbreak weeks",
    title = "Annual outbreak weeks by health region"
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
