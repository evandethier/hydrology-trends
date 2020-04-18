#### LIBRARIES ####
library(rnoaa)
library(modelr)
library(scales)
library(kdensity)
library(EnvCpt)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(grid)
library(gridExtra)
library(egg)
library(ggplotify)
library(dplyr)
library(zoo)
library(changepoint)
library(tidyverse)
library(readr)
library(readxl)
library(tidyquant)
library(RColorBrewer)
library(matrixStats)
library(smoother)
library(reshape2)
library(boot)
library(kernelboot)
library(np)
library(ggpubr)
library(gstat)
library(automap)
library(tidyr)
library(modelr)
library(sp)
library(USAboundaries)
library(sf)
library(rgeos)
library(raster)
library(rgdal)
library(maptools)
library(PBSmapping)
library(svglite)
library(ecp)
library(dataRetrieval)
library(scatterpie)
library(maps)
library(ggpubr)
library(devtools)
library(data.table)
library(splus2R)
library(Kendall)


#### THEMES ####
theme_evan <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed',color = 'grey70'),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    text = element_text(size=8),
    axis.text = element_text(size = 8), 
    plot.title = element_text(size = 9)
  )

theme_evan_facet <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    # legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    strip.background = element_rect(fill = 'white'),
    text = element_text(size=8),
    axis.text = element_text(size = 8), 
    plot.title = element_text(size = 9)
  )
season_facet <- theme_evan_facet + theme(
  legend.position = 'none', 
  strip.background = element_blank(),
  strip.text = element_text(hjust = 0, margin = margin(0,0,0,0, unit = 'pt'))
)

# Converts log10 axis values to format 10^x
fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- log10(l)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(parse(text = paste("10^",as.character(l),sep = "")))
} 
#### SET DIRECTORIES ####
# Set root directory
wd_root <- "~/hydrology-trends"

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root,"/hydro-imports/")
# Exports folder (save all figures, tables here)
wd_exports <- paste0(wd_root,"/hydro-exports/")

wd_figures <- paste0(wd_exports, "hydro-figures/")

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_exports, wd_figures)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}
#### IMPORT AND CLEAN CLUSTER DATA, INITIALIZE ####
setwd(wd_imports)
#separate main dataframe into different regions

    # Replace cluster numbers with new numbers for better ordering
    replace_cl_nos <- c(4,3,9,10,11,6,NA,5,12,7,2,8,NA,1,NA) # ordered list of replacement
    hydroregion_names <- c('Hawaii', 'Pacific NW','Pacific Coast','Appalachians','Mid-A Lowlands','NE/upper-Midwest','Rocky Highlnd', 
                           'High Plains', 'Rocky Mtns','Midwest','CA Rockies','Southeast')
    # Import HCDN sites 
    station_info <- data.table(read_excel("hcdn-huc08.xlsx"))
    
    # Import stations with cluster attached
    cluster_file <- fread("cluster_data_51419.dat")[Cluster_15_dist < 8][, ':='(
      agency_cd = ifelse(nchar(site_number) == 7, 'WSC', 'USGS'))
    ][,cluster_clean:=replace_cl_nos[Cluster_15]][
      !is.na(cluster_clean)][
        , ':='(site_no = site_number,
               cluster = cluster_clean)
        ]
    
    # [, .(agency_cd, site_number, station_name, state, lat, long, drainage_area_sqkm, eco03region, elevation, Cluster_15)]

cluster_file_export <- cluster_file # subset station characteristic file by those columns
cluster_list <- sort(unique(cluster_file$cluster_clean)) # sort by cluster
n_clusters <- length(cluster_list) # number of clusters that we will use

#### IMPORT ALL DATA  ####

# Initialize a list to store Discharge timeseries data.tables
# Each hydro-region (cluster) is assigned to a list item
cluster_files_master <- vector("list", length(min(cluster_list):max(cluster_list)))

# Generate import names for files from the WSC and the USGS, which have different naming conventions
cluster_files <- cluster_file[, .(file_name = ifelse(agency_cd == 'WSC', 
                                                    paste0('Q_daily/',site_number,'.csv'),
                                                    paste0('Q_daily/', site_number,'.dat')))
                                                    ]
# Make a vector of all file names
file_names_all <- list.files(paste0(wd_imports,'/Q_daily/')) 

# Initialize columns and categories
# Select timeframes (by season, Annual)
timeframe <- c('SON','DJF','MAM','JJA','Annual')
timeframes <- list(c(5:7), c(9:11), c(12,1,2), c(3:5), c(1:12))

# Select recurrence interval. 
# May want to change the naming convention here, since n_peaks is misleading.
n_peaks <- c(2,3,5,10,25)

# Functions for importing and cleaning discharge timeseries data
  # Import and standardize Discharge timeseries from USGS and WSC
  # TO DO: make sure site_no is correct 8 digit string for USGS
  getQ_ts <- function(file_name){
    if(grepl('.dat',file_name)){
      Q_ts_import <- fread(file_name)[,':='(
        Q_cms = q_cfs * 0.0283168, # Convert cfs to cms for USGS sites
        Date = ymd(Date),
        q_cfs = NULL,
        V1 = NULL,
        `qual code` = NULL
      )][,.(agency_cd,site_no, Date, water_year,Q_cms)]
    }else{
      Q_ts_import <- fread(file_name)[,':='(
        Date = ymd(Date)
      )
      ][,.(agency_cd,site_no, Date, water_year,Q_cms)]
    }
    return(Q_ts_import)
  }
  
  # Identify peaks and minima for discharge timeseries
  getPeaks <- function(dt){
    set.seed(0)
    # Get all 7-day maxima and minima. Will need to take the extreme X% later
    return(dt[,month:=month(Date)][
      peaks(dt$Q_cms, span = 7), Q_peaks := Q_cms][ # high-flows: maximum Q in 7-day period
        peaks(-dt$Q_cms, span = 7), Q_lows := Q_cms]) # low-flows: minimum Q in a 7 day period
  }
  
  # Get dates of continuous daily discharge sampling at each site
  # Return start year, end year, number of years
  getContinuous <- function(dt){
    last_missing <- max(which(diff(dt$Date)>(365*2)),1) # find timeseries gaps of 2+ years
    startYear <- year(dt$Date[last_missing]) #end year of gap
    endYear <- year(max(dt$Date)) # end year of timeseries
    record_yrs <- endYear-startYear # years span of timeseries record
    # Return site number, start year, end year, and years of record in data.table format
    return(data.table(site_no = dt$site_no[1], startYear = startYear, endYear = endYear, record_years = record_yrs))
  }
  
# Import and standardize the cluster data
# TO DO: update usgs and hydat files to use latest data (up to october 2019)
for(n_cluster in 1:length(cluster_list)){
  # Get cluster number
  cluster_sel <- cluster_list[n_cluster]
  # Get site info for sites in cluster from data.table of all sites
  cluster_sites <- cluster_file[cluster_clean==cluster_sel]
  # Get file names for stations in cluster
  cluster_files <- cluster_sites[, .(file_name = ifelse(agency_cd == 'WSC', 
                                                       paste0(site_number,'.csv'),
                                                       paste0(site_number,'.dat')))]$file_name
  
  # Match file names of sites in cluster to timeseries filenames in folder
  # There may be some sites that are missing a file, or the file was incomplete and discarded, etc.
  file_list <- paste0('Q_daily/',cluster_files[cluster_files %in% file_names_all])
  
  # Make sure that there are sufficient sites for analysis
  if(length(file_list)>3){
    print(cluster)
    print(length(file_list))
  # Apply import functions to each station discharge timeseries
    # First import timeseries and standardize USGS and WSC data
    cluster_files_temp <- lapply(file_list, getQ_ts)
    # Get peaks (minimum and maximum in 7-day window)
    cluster_files_temp <- lapply(cluster_files_temp,getPeaks)
  # Copy station timeseries data to cluster master list for subsequent analysis
    cluster_files_master[[cluster_sel]] <- cluster_files_temp
  }}
#   
 
#### CALCULATE TRENDS USING COX-LEWIS Z-SCORES #### 
  # Return top N event rows for all data.tables in a list
  # Iterates over all given month categories
  # Also iterates over all n event categories
  getExtremes_high <- function(dt,n_peaks, months, start_year, end_year){
    
    dt[,month:=month(Date)]
    return(rbindlist(lapply(months, function(months){ # Select top n events by given months
      # Make variable with given month three-letter tag/'Annual'
      month_diagnostic <- ifelse(length(months) == 12, 5, 
                                 ifelse(max(months) == 11, 1, 
                                        ifelse(max(months) == 12,2, 
                                               ifelse(max(months) == 5, 3,4))))
      month_sel <- timeframe[month_diagnostic]
      # Get peak events for each of selected number of peaks (n_peaks)
      return(rbindlist(lapply(n_peaks, function(n_peaks, dt){
        # For each number of n_peaks, get number of peaks for each starting year
        return(dt[,':='(
          n_peaks_sel = n_peaks, # add column with selected number of peaks
          recur = round((end_year - start_year)/n_peaks), # add column with approximate recurrence interval
          timeframe_sel = month_sel)][ # add column with selected timeframe (month range or annual)
            # select only rows that are within year and month ranges
            water_year > start_year & water_year <= end_year & month %in% months][ 
              order(Q_peaks*-1)][ # order by discharge (multiply by -1 so list is descending, with top events at the top)
                1:round((end_year - start_year)/n_peaks)])}, dt = dt)))}))) # select top n_peaks events
  }
  
  # Return minimum N event rows for all data.tables in a list
  # Iterates over all given month categories
  # Also iterates over all n event categories
  # Limited to minimum event per year due to memory effects
  getExtremes_low <- function(dt,n_peaks, months, start_year, end_year){
    
    dt[,month:=month(Date)]
    return(rbindlist(lapply(months, function(months){ # Select top n events by given months
      # Make variable with given month three-letter tag/'Annual'
      month_diagnostic <- ifelse(length(months) == 12, 5, 
                                 ifelse(max(months) == 11, 1, 
                                        ifelse(max(months) == 12,2, 
                                               ifelse(max(months) == 5, 3,4))))
      month_sel <- timeframe[month_diagnostic]
      # Get peak events for each of selected number of peaks (n_peaks)
      return(rbindlist(lapply(n_peaks, function(n_peaks, dt){
        # For each number of n_peaks, get number of peaks for each starting year
        return(dt[,':='(
          n_peaks_sel = n_peaks, # add column with selected number of peaks
          recur = round((end_year - start_year)/n_peaks), # add column with approximate recurrence interval
          timeframe_sel = month_sel)][ # add column with selected timeframe (month range or annual)
            # select only rows that are within year and month ranges
            water_year > start_year & water_year <= end_year & month %in% months][ 
              , Q_lows_sel := ifelse(Q_lows == min(Q_lows, na.rm = T), Q_lows, NA), by = .(water_year)][
                order(Q_lows_sel)][ # order by discharge (multiply by -1 so list is descending, with top events at the top)
                  1:round((end_year - start_year)/n_peaks)])}, dt = dt)))}))) # select top n_peaks events
  }
  
  
  # Make a data.table of Z-Score significance thresholds for plotting
  significance_dt <- data.table(timeframe_sel = rep(timeframe,5), 
                                significance = sort(rep(c(-2.575,-1.96,0,1.96,2.575),5)))     
    # For each cluster loop through all sites
    
    # For each site loop through each season for both floods and droughts
    # Each result returns a data.table of dates and magnitudes
    # Bind all 
    # TO DO: also bind N stations to output tables 
for(i in 1:n_clusters){
  n_cluster <- cluster_list[i]
    # dt_peaks <- lapply(cluster_files_master[[i]], getPeaks)
    dt_peaks <- cluster_files_master[[i]]
    
    # Generate list of start and end dates for stations within hydro-region
    dt_startDate <- rbindlist(lapply(dt_peaks, getContinuous))
    
    # Generate list of start years (first year in each decade with observations)
    start_years <- c((min(dt_startDate[,startYear - startYear%%10])/10 + 1):198)*10
    
    for(k in 1:length(start_years)){
      # Set record length parameters
      # Test a range of start years
      start_year <- start_years[k]
      # TO DO: determine end year based on average end year in the cluster
      # Will allow for more canadian data
      end_year = 2016
      record_length_sel <- end_year - start_year
      
      
      # Eliminate sites with incomplete records
      # Select rows in data.table of first and last measurement that begin before selected start year
      start_year_rows <- which(dt_startDate$startYear <= start_year & dt_startDate$endYear >= 2016)
      if(length(start_year_rows)>0){
      # Count number of stations with adequate record length
      n_stations <- length(start_year_rows)
      
      # Subset list of observation data.tables with start dates before selected start year
      dt_peaks_byStartYr <- dt_peaks[start_year_rows]

      # Find the year (and magnitude) of extreme events for selected time range
      # Accomplished by sorting data.table by high and low flow magnitude, respectively, and selecting top-N events
      # High-flow events
      dt_extremes <- rbindlist(lapply(dt_peaks_byStartYr, getExtremes_high, n_peaks = n_peaks, 
                                      months = timeframes, 
                                      start_year = start_year, end_year = end_year))
      # Low-flow events
      dt_extremes_low <- rbindlist(lapply(dt_peaks_byStartYr, getExtremes_low, n_peaks = n_peaks, 
                                      months = timeframes, 
                                      start_year = start_year, end_year = end_year))
        
      # For entire cluster:
      # Get number of events/year for each N event recurrence interval
      # High-flow events (by cluster)
      dt_count <- dt_extremes[, .(count = .N), by = .(water_year, n_peaks_sel, recur, timeframe_sel)]
      # Low-flow events (by cluster)
      dt_count_low <- dt_extremes_low[, .(count = .N), by = .(water_year, n_peaks_sel, recur, timeframe_sel)]
      
      # for each station
      # Get number of events/year for each N event recurrence interval
      # High-flow events (individual stations)
      dt_count_indiv <- dt_extremes[, .(count = .N), by = .(site_no, water_year, n_peaks_sel, recur, timeframe_sel)]
      # Low-flow events (individual stations)
      dt_count_low_indiv <- dt_extremes_low[, .(count = .N), by = .(site_no, water_year, n_peaks_sel, recur, timeframe_sel)]
      
      # Calculate Cox-Lewis statistic (Z-Score) for each N event recurrence interval, each time period
      # High-flow events (by cluster)
      dt_z <- dt_count[,.(
                z = (sum(count*water_year/sum(count)) - (end_year + start_year)/2)/
                  ((end_year - start_year)/sqrt(12*sum(count)))), by = .(n_peaks_sel, timeframe_sel)][
                    , ':='(cluster = n_cluster,
                           start_year = start_year)
                  ]
      # Low flow events (by cluster)
      dt_z_low <- dt_count_low[,.(
                z = (sum(count*water_year/sum(count)) - (end_year + start_year)/2)/
                  ((end_year - start_year)/sqrt(12*sum(count)))), by = .(n_peaks_sel, timeframe_sel)][
                    , ':='(cluster = n_cluster,
                           start_year = start_year)
                  ]
      
      # High-flow events (individual stations)
      dt_z_indiv <- dt_count_indiv[,.(
                z = (sum(count*water_year/sum(count)) - (end_year + start_year)/2)/
                  ((end_year - start_year)/sqrt(12*sum(count)))), by = .(n_peaks_sel, timeframe_sel, site_no)][
                    , ':='(cluster = n_cluster,
                           start_year = start_year)
                  ]
      # Low-flow events (individual stations)
      dt_z_low_indiv <- dt_count_low_indiv[,.(
                z = (sum(count*water_year/sum(count)) - (end_year + start_year)/2)/
                  ((end_year - start_year)/sqrt(12*sum(count)))), by = .(n_peaks_sel, timeframe_sel, site_no)][
                    , ':='(cluster = n_cluster,
                           start_year = start_year)
                  ]
      
      ## PLOT SECTION
      # Plot Z-Scores for inidividual stations in a cluster
      # One row plotted for each approximate recurrence interval
      # Faceted by season/annual: SON,DJF,MAM,JJA, Annual
          # For high-flow events
          dt_z_high_indiv_plot <- ggplot(dt_z_indiv, aes(y = z, x = reorder(as.character(n_peaks_sel),n_peaks_sel), color = factor(n_peaks_sel))) + 
            geom_point(position = 'jitter') +
            geom_hline(data = significance_dt, 
                       aes(yintercept = significance), lwd = 0.25) +
            # geom_boxplot(fill = NA) +
            scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
            season_facet + 
            scale_y_continuous(limits = c(-6,6)) +
            facet_wrap(.~factor(timeframe_sel, levels = rev(timeframe)), nrow = 5) +
            labs(
              x = 'N high-flow events in record',
              y = 'Z-Score'
            ) +
            rotate()
          
          # For low-flow events
          dt_z_low_indiv_plot <- ggplot(dt_z_low_indiv, aes(y = z, x = reorder(as.character(n_peaks_sel),n_peaks_sel), color = factor(n_peaks_sel))) + 
            geom_point(position = 'jitter') +
            geom_hline(data = significance_dt, 
                       aes(yintercept = significance), lwd = 0.25) +
            # geom_boxplot(fill = NA) +
            scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
            season_facet + 
            scale_y_continuous(limits = c(-6,6)) +
            facet_wrap(.~factor(timeframe_sel, levels = rev(timeframe)), nrow = 5) +
            labs(
              x = 'N low-flow events in record',
              y = 'Z-Score'
            ) +
            rotate()
          
          # Combine high- and low-flow Z-Scores plot for a given cluster, start year
          indiv_z_combined_plot <- ggpubr::ggarrange(dt_z_low_indiv_plot, dt_z_high_indiv_plot, ncol = 1)
          
          # Save combined high- and low-flow Z-Scores plot
          ggsave(indiv_z_combined_plot, width = 3, height = 5.5,
                 filename = paste0(wd_exports,'Q_z_indiv_cl',n_cluster,'_',start_year,'.pdf'), useDingbats = F)
          ggsave(indiv_z_combined_plot, width = 3, height = 5.5,
                 filename = paste0(wd_exports,'Q_z_indiv_cl',n_cluster,'_',start_year,'.png'))
          
          # Write aggregated cluster z-scores to file 
              # (only write once/cluster: after final start year calculation completed)
          if(k == length(start_years)){
            # High-flow events
            fwrite(dt_z_cluster, paste0(wd_exports,'Q_trends_high_sensitivity_cl',n_cluster,'.csv')) 
            # Low-flow events
            fwrite(dt_z_low_cluster, paste0(wd_exports,'Q_trends_low_sensitivity_cl',n_cluster,'.csv')) 
          }
          
          
      # Plot smoothed trends, colored by approximate recurrence interval, facet by month category
          # Annual trends, high-flow
          Q_high_trends_ann <- ggplot(dt_count[timeframe_sel == 'Annual'], 
                 aes(x = water_year, y = count, color = factor(n_peaks_sel))) + 
            # geom_point() +
            facet_wrap(~factor(timeframe_sel, levels = timeframe)) +
            geom_smooth(se = F, lwd = 0.25) +
            season_facet + 
            scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
            theme(legend.position = 'right') +
            labs(x = 'Water year',
                 y = 'Total N hydro-region events/yr',
                 color = paste0('Avg. years/event\nat-a-station\n(N = ', n_stations,' stns)'))
          
          # Monthly trends, low-flow
          Q_high_trends_month <- ggplot(dt_count[timeframe_sel != 'Annual'], 
                 aes(x = water_year, y = count, color = factor(n_peaks_sel))) + 
            # geom_point() +
            facet_wrap(~factor(timeframe_sel, levels = timeframe)) +
            geom_smooth(se = F, lwd = 0.25) +
            season_facet + 
            scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
            theme(legend.position = 'right') +
            labs(x = 'Water year',
                 y = 'Total N hydro-region events/yr',
                 color = paste0('Avg. years/event\nat-a-station\n(N = ', n_stations,' stns)'))
          
          # Combine Annual and monthly trend plots, high-flow
          Q_high_trends_comb <- ggpubr::ggarrange(Q_high_trends_ann, Q_high_trends_month, nrow = 2, common.legend = T)
          
          # Export high-flow event trend plots
          ggsave(Q_high_trends_comb, width = 3.5, height = 5, 
                 filename = paste0(wd_exports,'Q_high_trends_sensitivity_cl',n_cluster,'_',start_year,'.pdf'), useDingbats = F)
          ggsave(Q_high_trends_comb, width = 3.5, height = 5, 
                 filename = paste0(wd_exports,'Q_high_trends_sensitivity_cl',n_cluster,'_',start_year,'.png'))
          
          # Annual trends, low flow
          Q_low_trends_ann <- ggplot(dt_count_low[timeframe_sel == 'Annual'], 
                 aes(x = water_year, y = count, color = factor(n_peaks_sel))) + 
            # geom_point() +
            facet_wrap(~factor(timeframe_sel, levels = timeframe)) +
            geom_smooth(se = F, lwd = 0.25) +
            season_facet + 
            scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
            theme(legend.position = 'right') +
            labs(x = 'Water year',
                 y = 'Total N hydro-region events/yr',
                 color = paste0('Avg. years/event\nat-a-station\n(N = ', n_stations,' stns)'))
          
          # Monthly trends, low flow
          Q_low_trends_month <- ggplot(dt_count_low[timeframe_sel != 'Annual'], 
                 aes(x = water_year, y = count, color = factor(n_peaks_sel))) + 
            # geom_point() +
            facet_wrap(~factor(timeframe_sel, levels = timeframe)) +
            geom_smooth(se = F, lwd = 0.25) +
            season_facet + 
            scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
            theme(legend.position = 'right') +
            labs(x = 'Water year',
                 y = 'Total N hydro-region events/yr',
                 color = paste0('Avg. years/event\nat-a-station\n(N = ', n_stations,' stns)'))
          
          # Combine Annual and monthly trend plots, low-flow
          Q_low_trends_comb <- ggpubr::ggarrange(Q_low_trends_ann, Q_low_trends_month, nrow = 2, common.legend = T)
          
          # Export low-flow event trend plots
          ggsave(Q_low_trends_comb, width = 3.5, height = 5, 
                 filename = paste0(wd_exports,'Q_low_trends_sensitivity_cl',n_cluster,'_',start_year,'.pdf'), useDingbats = F)
          ggsave(Q_low_trends_comb, width = 3.5, height = 5, 
                 filename = paste0(wd_exports,'Q_low_trends_sensitivity_cl',n_cluster,'_',start_year,'.png'))
      
          # For each start date, and for individual and cluster-aggregate calculations, save Z-Scores, event counts per water-year
          if(k == 1){
            dt_z_cluster <- dt_z # high-flow Z-Scores, cluster
            dt_z_low_cluster <- dt_z_low # low-flow Z-Scores, cluster
            dt_z_indiv_cluster <- dt_z_indiv # high-flow Z-Scores, individual
            dt_z_low_indiv_cluster <- dt_z_low_indiv # low-flow Z-Scores, individual
            dt_count_cluster <- dt_count # high-flow event counts, cluster
            dt_count_low_cluster <- dt_count_low # low-flow event counts, cluster
            dt_count_indiv_cluster <- dt_count_indiv # high-flow event counts, individual
            dt_count_low_indiv_cluster <- dt_count_low_indiv # low-flow event counts, individual
            
          }else{
            dt_z_cluster <- rbind(dt_z_cluster, dt_z) # high-flow Z-Scores, cluster
            dt_z_low_cluster <- rbind(dt_z_low_cluster, dt_z_low) # low-flow Z-Scores, cluster
            dt_z_indiv_cluster <- rbind(dt_z_indiv_cluster, dt_z_indiv) # high-flow Z-Scores, individual
            dt_z_low_indiv_cluster <- rbind(dt_z_low_indiv_cluster, dt_z_low_indiv) # low-flow Z-Scores, individual
            dt_count_cluster <- rbind(dt_count_cluster,dt_count) # high-flow event counts, cluster
            dt_count_low_cluster <- rbind(dt_count_low_cluster,dt_count_low) # low-flow event counts, cluster
            dt_count_indiv_cluster <- rbind(dt_count_indiv_cluster,dt_count_indiv) # high-flow event counts, individual
            dt_count_low_indiv_cluster <- rbind(dt_count_low_indiv_cluster,dt_count_low_indiv) # low-flow event counts, individual
          }}
      
    }
    
    # After cluster calculations are complete
    # Bind together Z-Score, event counts per water-year
    # Each cluster identified in cluster tables
    # Each cluster and site identified in individual tables
      if(i == 1){
        dt_z_cluster_master <- dt_z_cluster # high-flow Z-Scores, cluster
        dt_z_low_cluster_master <- dt_z_low_cluster # low-flow Z-Scores, cluster
        dt_z_indiv_cluster_master <- dt_z_indiv_cluster # high-flow Z-Scores, individual
        dt_z_low_indiv_cluster_master <- dt_z_low_indiv_cluster # low-flow Z-Scores, individual
        dt_count_cluster_master <- dt_count_cluster # high-flow event counts, cluster
        dt_count_low_cluster_master <- dt_count_low_cluster # low-flow event counts, cluster
        dt_count_indiv_cluster_master <- dt_count_indiv_cluster # high-flow event counts, individual
        dt_count_low_indiv_cluster_master <- dt_count_low_indiv_cluster # low-flow event counts, individual
      }else{
        dt_z_cluster_master <- rbind(dt_z_cluster_master,dt_z_cluster) # high-flow Z-Scores, cluster
        dt_z_low_cluster_master <- rbind(dt_z_low_cluster_master,dt_z_low_cluster) # low-flow Z-Scores, cluster
        dt_z_indiv_cluster_master <- rbind(dt_z_indiv_cluster_master,dt_z_indiv_cluster) # high-flow Z-Scores, individual
        dt_z_low_indiv_cluster_master <- rbind(dt_z_low_indiv_cluster_master, dt_z_low_indiv_cluster) # low-flow Z-Scores, individual
        dt_count_cluster_master <- rbind(dt_count_cluster_master,dt_count_cluster) # high-flow event counts, cluster
        dt_count_low_cluster_master <- rbind(dt_count_low_cluster_master,dt_count_low_cluster) # low-flow event counts, cluster
        dt_count_indiv_cluster_master <- rbind(dt_count_indiv_cluster_master,dt_count_indiv_cluster) # high-flow event counts, individual
        dt_count_low_indiv_cluster_master <- rbind(dt_count_low_indiv_cluster_master,dt_count_low_indiv_cluster) # low-flow event counts, individual
      }
    
}

#### AGGREGATE CLUSTER Z-SCORE PLOTS ####

    # CREATE THIS PLOT:
      # For a given recurrence interval
      # x = Z-Score, y = season, color = timeframe
    # Plot Z-Scores for inidividual stations in a cluster
    # Facet by high-flow/low-flow events
    # Arrange on grid by hydro-region (cluster)
  
  # First set labeling and color for trends of different significance
  cluster_color_categ <- c('Drier (p < 0.01)','Drier (p < 0.05)','Not sign. (p > 0.05)','Wetter (p < 0.05)','Wetter (p < 0.01)')
  cluster_color_values <- c('#a50f15','#de2d26','grey70','#3182bd','#08519c')
  names(cluster_color_values) <- c('Drier (p < 0.01)','Drier (p < 0.05)','Not sign. (p > 0.05)','Wetter (p < 0.05)','Wetter (p < 0.01)')
  
    # Combine high- and low-flow event tables
      # For cluster aggregate
      # Assign cluster trend significance for each row
      dt_z_comb_cluster_master <- na.omit(rbind(dt_z_cluster_master[,':='(event_type ='High-flow', event_sign = 1)],
                                        dt_z_low_cluster_master[,':='(event_type ='Low-flow', event_sign = -1)])[
                                          , hydroregion_name := hydroregion_names[cluster]
                                        ][,':='(cluster_color = ifelse(z * event_sign < -2.575, cluster_color_categ[1], # code to dark red
                                                                ifelse(z * event_sign < -1.96, cluster_color_categ[2], # light red
                                                                ifelse(z * event_sign > 2.575, cluster_color_categ[5], # dark blue
                                                                ifelse(z * event_sign > 1.96, cluster_color_categ[4], # light blue
                                                                                              cluster_color_categ[3]))))) # grey
                                          ], cols = c('z'))
      
      # For individual stations
      dt_z_comb_indiv_cluster_master <- na.omit(rbind(dt_z_indiv_cluster_master[,':='(event_type ='High-flow', event_sign = 1)],
                                              dt_z_low_indiv_cluster_master[,':='(event_type ='Low-flow', event_sign = -1)])[
                                                , hydroregion_name := hydroregion_names[cluster]
                                                ], cols = c('z'))
      
      # Add color based on cluster, timeframe, start_year signficance. 
      # Color for Z-Score plots
      # dark red (1): z_cluster < -2.575, event_type == 'High-flow'
      #           z_cluster > 2.575, event_type == 'Low-flow'
      # red (2): z_cluster < -1.96, event_type == 'High-flow'
      #           z_cluster > 1.96, event_type == 'Low-flow'
      # dark blue (3): z_cluster > 2.575, event_type == 'High-flow'
      #           z_cluster < -2.575, event_type == 'Low-flow'
      # blue (4): z_cluster > 1.96, event_type == 'High-flow'
      #           z_cluster < 1.96, event_type == 'Low-flow'
      # else (5): grey40
      
      # Requires joining dt_z_comb_cluster_master and dt_z_comb_indiv_cluster_master
      setkey(dt_z_comb_cluster_master, cluster, timeframe_sel, n_peaks_sel, start_year, event_type, event_sign, hydroregion_name)
      setkey(dt_z_comb_indiv_cluster_master, cluster, timeframe_sel, n_peaks_sel, start_year, event_type, event_sign, hydroregion_name)
      # For some reason some rows in cluster 10 are duplicated. This is a temporary fix. TO DO: remove duplicated rows. 
      # I think duplicated rows are gone
      
      dt_z_comb_indiv_cluster_master <- dt_z_comb_indiv_cluster_master[
            dt_z_comb_cluster_master][
              , ':='(z_cluster = i.z, i.z = NULL)
              # Assign trend and signficance to a wetter/drier trend and p-value threshold < 0.05 vs 0.01
              ]
      
      # Make a data.table of Z-Score significance thresholds for plotting
      significance_dt_high_low <- data.table(event_typ = rep(c('High-flow','Low-flow'),5), 
                                    significance = sort(rep(c(-2.575,-1.96,0,1.96,2.575),2)))
    
    
    # Function to generate a Z-Score plot for a selected cluster, given an input data.table
    getZplot <- function(cluster_sel, dt, start_year_sel, n_peaks){
      # return(ggplot(dt[cluster == cluster_sel & start_year < 1970 & start_year > 1930], # for a range of start years
      dt_plot <- dt[cluster == cluster_sel & start_year == start_year_sel & n_peaks_sel == n_peaks]
      hydroregion_name <- dt_plot$hydroregion_name[1]
      return(
          ggplot(dt_plot, # for a single start year
                                                 aes(y = z, x = factor(timeframe_sel, levels = rev(timeframe)),
                                                     color = factor(cluster_color, levels = cluster_color_categ, ordered = T))) + 
        # geom_boxplot() +
        geom_hline(data = significance_dt_high_low,
                   aes(yintercept = significance), lwd = 0.25) +
        geom_point(position = position_jitter(0.2), size = 0.75) +
        # geom_boxplot(fill = NA) +
        scale_color_manual(values = cluster_color_values) +
        # scale_color_manual(values = rev(c('#A6170A','grey80','#5890A6', '#1D4859', 'grey40'))) +
        season_facet +
        scale_y_continuous(limits = c(-6,6)) +
        facet_wrap(.~event_type, ncol = 1, strip.position = 'left') +
        theme(strip.placement = 'outside',
              strip.text = element_text(size = 10, hjust = 0.5),
              axis.ticks.y = element_blank(),
              axis.title = element_blank(),
              plot.title = element_text()
              ) +
        labs(
          # x = '',
          # y = 'Z-Score',
          title = paste0(cluster_sel, '-',hydroregion_name)
        ) +
        rotate())
    }
    
    start_year_sel <- 1960
    # Get individual Z-Score plots for high and low flows for each cluster. Output is a list of ggplots.
    zplots_indiv_byCluster <- lapply(cluster_list, getZplot, dt = dt_z_comb_indiv_cluster_master, 
                                     start_year_sel = start_year_sel, n_peaks = 3)
    # Example Z-Score plot
    # getZplot(2, dt = dt_z_comb_indiv_cluster_master, start_year_sel = start_year_sel, n_peaks = 3)
    
    # Arrange list of individual Z-Score plots in a grid by hydro-region (cluster)
    zplots_indiv_byCluster_comb <- annotate_figure(
                  ggpubr::ggarrange(plotlist = zplots_indiv_byCluster, ncol = 3, nrow = 4),
                  'bottom' = text_grob('Z-Score'))
                  
    # Save Z-Score combination plot with all clusters
    ggsave(zplots_indiv_byCluster_comb, width = 7, height = 8, 
           filename = paste0(wd_figures,'Q_ZScore_indiv_combined_',start_year_sel,'.pdf'), useDingbats = F)
    ggsave(zplots_indiv_byCluster_comb, width = 7, height = 8, 
           filename = paste0(wd_figures,'Q_ZScore_indiv_combined_',start_year_sel,'.png'))
  
#### VISUALIZE TREND SENSITIVITY TO START YEAR ####

    
    # Make a grid of plots for each cluster, season
    # Each grid element is a plot with x-axis start year, 
    # y-axis corresponding to high/low-flow event and magnitude (10 possibilities)
    # Color is given by trend significance
    z_score_start_year_sensitivity_grid_plot <- ggplot(dt_z_comb_cluster_master[start_year < 1970], 
                                                  aes(x = factor(start_year), 
                                                      y = factor(paste0(event_type, n_peaks_sel), 
                                                                 levels = event_type_levels,
                                                                 labels = event_type_labels), 
                                                      fill = factor(cluster_color, levels = cluster_color_categ))) + 
      geom_tile() +
      facet_grid(reorder(paste0(cluster, '-', hydroregion_name), cluster)~factor(timeframe_sel, levels = timeframe)) + 
      scale_fill_manual(values = cluster_color_values) +
      season_facet + 
      theme(legend.position = 'right',
            axis.text.x = element_text(angle = 90)) +
      labs(
        x = 'Start year',
        y = 'Event recurrence interval (yr)',
        fill = 'Trend (signif.)'
      )
    
    ggsave(z_score_start_year_sensitivity_grid_plot, width = 8, height = 10, 
           filename = paste0(wd_exports,'z_score_start_year_n',n_peaks_sel,'_sensitivity_grid_plot.pdf'), useDingbats = F)
    ggsave(z_score_start_year_sensitivity_grid_plot, width = 8, height = 10, 
           filename = paste0(wd_exports,'z_score_start_year_n',n_peaks_sel,'_sensitivity_grid_plot.png')) 
    
    # Slight modification on above plot
    # Make two grid plots: each element of gird corresponds to a unique cluster, season, event type
    # Each grid element is a plot with x-axis start year, 
    # y-axis corresponding to high OR low-flow event and magnitude (5 possibilities for each plot)
    # Color is given by trend significance
    # Change size of facet label text
    strip_text_size = 8.5
    z_score_start_year_sensitivity_HF_grid_plot <- ggplot(dt_z_comb_cluster_master[start_year < 1970 & event_type == 'High-flow'], 
                                                  aes(x = start_year, 
                                                      y = factor(paste0(event_type, n_peaks_sel), 
                                                                 levels = event_type_levels,
                                                                 labels = event_type_labels), 
                                                      fill = factor(cluster_color, levels = cluster_color_categ))) + 
      geom_tile() +
      facet_grid(reorder(paste0(cluster, '\n', gsub(' ', '\n',
                                                    gsub('upper-Mid','upper- Mid',hydroregion_name))), cluster)~
                   factor(timeframe_sel, levels = timeframe)) +  
      scale_fill_manual(values = cluster_color_values) +
      # scale_x_continuous(breaks = c(1920,1940,1960)) +
      season_facet + 
      theme(legend.position = 'right',
            axis.text.x = element_text(angle = 90, size = strip_text_size, vjust = 0.5),
            axis.text.y = element_text(size = strip_text_size),
            strip.text = element_text(size = strip_text_size),
            axis.title = element_text(size = strip_text_size + 1),
            legend.text = element_text(size = strip_text_size),
            legend.title = element_text(size = strip_text_size + 1)) +
      labs(
        x = 'First year in analysis',
        y = 'Record-averaged recurrence interval (yr)',
        fill = 'Trend (signif.)'
      )
    
    z_score_start_year_sensitivity_LF_grid_plot <- ggplot(dt_z_comb_cluster_master[start_year < 1970 & event_type == 'Low-flow'], 
                                                  aes(x = start_year, 
                                                      y = factor(paste0(event_type, n_peaks_sel), 
                                                                 levels = event_type_levels,
                                                                 labels = event_type_labels), 
                                                      fill = factor(cluster_color, levels = cluster_color_categ))) + 
      geom_tile() +
      facet_grid(reorder(paste0(cluster, '\n', gsub(' ', '\n',
                                               gsub('upper-Mid','upper- Mid',hydroregion_name))), cluster)~
                   factor(timeframe_sel, levels = timeframe)) + 
      # labeller = label_wrap_gen(width=10)
      scale_fill_manual(values = cluster_color_values) +
      # scale_x_continuous(breaks = c(1920,1940,1960)) +
      season_facet + 
      theme(legend.position = 'right',
            axis.text.x = element_text(angle = 90, size = strip_text_size, vjust = 0.5),
            axis.text.y = element_text(size = strip_text_size),
            strip.text = element_text(size = strip_text_size),
            axis.title = element_text(size = strip_text_size + 2),
            legend.text = element_text(size = strip_text_size),
            legend.title = element_text(size = strip_text_size + 1)) +
      labs(
        x = 'First year in analysis',
        y = 'Record-averaged recurrence interval (yr)',
        fill = 'Trend (signif.)'
      )
    
    z_score_start_year_sensitivity_comb_grid_plot <- ggpubr::ggarrange(z_score_start_year_sensitivity_HF_grid_plot,
                                                               z_score_start_year_sensitivity_LF_grid_plot,
                                                               nrow = 1, common.legend = T, labels = c('A','B'))
    
    ggsave(z_score_start_year_sensitivity_comb_grid_plot, width = 8, height = 11, 
           filename = paste0(wd_exports,'z_score_start_year_sensitivity_comb_grid_plot.pdf'), useDingbats = F)
    ggsave(z_score_start_year_sensitivity_comb_grid_plot, width = 8, height = 11, 
           filename = paste0(wd_exports,'z_score_start_year_sensitivity_comb_grid_plot.png'))
    




#### COLOR PALETTE ####
# set color palette for remaining plots
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  
pal_length <- length(cluster_list)
cl_colors <- getPalette(pal_length)
names(cl_colors) <- cluster_list
cl_colors['12'] <- '#cab2d6'
cl_colors['8'] <- '#a6cee3'
cl_colors['5'] <- '#1f78b4'
cl_colors['6'] <- '#419486'
cl_colors['3'] <- "#D96D3B"
cl_colors['9'] <- "#999999"
cl_colors['10'] <- "#6a3d9a"
#### INITIALIZE MAP DATA FOR N.AMERICA ####
# get states shapefile for clipping/display
us_states <-  us_boundaries(map_date = NULL, type = c("state"), resolution = c("low"), states = NULL)
# us_states <- us_states[us_states$state_abbr != "AK" & us_states$state_abbr != "HI" & us_states$state_abbr != "PR",]

# convert shapefile to Spatial class
us_states <- as(us_states, 'Spatial')

# plot(us_states)
# plot(us_states_merge)

# set projection for states shapefile
projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# import canadian provinces
canada_prov <- read_sf(dsn = "Canada", layer = "Canada")
# plot(canada_prov)
canada_geom <- st_geometry(canada_prov)
# attributes(canada_geom)

# do conversions and projections for canadian provinces to match us states
canada_prov <- as(canada_prov, 'Spatial')

canada_prov <- spTransform(canada_prov, projection)
# coordinates(canada_prov) <- ~long+lat
proj4string(canada_prov) <- proj4string(us_states)

names(canada_prov) <- c('name','frenchName') # rename columns labeling canadian provinces
canada_prov <- canada_prov[,c('name')] # only select column with province name
us_states <- us_states[,c('name')] # only select column with province name
us_states <- rbind(us_states,canada_prov) # combine canadian and us shapefiles

# fortify state shapefile
us_ca <- fortify(us_states)
#### GENERATE POLYGONS OUTLINING CLUSTERS ON MAP ####
# Join cluster change data with station info
# Join dt_z_comb_indiv_cluster_master and cluster_file

# Need to make dt_z_comb_indiv_cluster_master site numbers 8 digits for USGS (TO DO fix during import, see above)
dt_z_comb_indiv_cluster_master[!grepl(pattern = '[A-z]', site_no) & nchar(site_no) == 7, site_no:=paste0(0,site_no)]
# Set keys
setkey(dt_z_comb_indiv_cluster_master, site_no, cluster)
setkey(cluster_file, site_no, cluster)

station_change_master <- dt_z_comb_indiv_cluster_master[
    cluster_file[site_no %in% unique(dt_z_comb_indiv_cluster_master$site_no)]] # select only sites used in analysis

# Get coordinates of convex hull of each cluster group
# Result will be polygon outline of cluster group with all cluster points bounded
cl_coords <- station_change_master[n_peaks_sel == 2 & timeframe_sel == 'Annual' & 
                                     start_year == 1960 &
                                      event_type == 'High-flow'][,.(cluster,lat,long)]
getCluster_polygon <- function(cluster_sel,cl_coords){
  cl_coords_sel <- cl_coords[cluster == cluster_sel]
  chull_sel <- chull(cl_coords_sel[,.(lat,long)])
  cl_hull <- cl_coords_sel[chull_sel]
  # cl_hull_coords <- cl_out_sel[c(cl_out_hull,cl_out_hull[1]),]
  return(cl_hull)
}

cluster_polygons <- rbindlist(lapply(cluster_list,getCluster_polygon, cl_coords = cl_coords))
cluster_centroids <- cl_coords[,.(long = mean(long), 
                                  lat = mean(lat)), by = cluster]

cluster_map_raw <- ggplot(cluster_polygons, aes(x = long, y = lat, 
                                                color = as.character(cluster), fill = as.character(cluster), 
                                                group = cluster)) + 
  geom_map(data = us_ca, map = us_ca, aes(map_id = id, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) +
  # geom_polygon(fill = NA, lty = 'dashed') + # Add bounding polygon for each cluster
  geom_point(data = cl_coords, color = 'grey30', pch = 21, stroke = 0.2) +
  geom_text_repel(data = cluster_centroids, aes(label = cluster), color = 'black', size = 5, point.padding = NA, force = 0.1) +
  season_facet +
  scale_color_manual(values = cl_colors) +
  scale_fill_manual(values = cl_colors) +
  scale_y_continuous(limits = c(23, 60)) + 
  scale_x_continuous(limits = c(-131.2, -54)) +
  # guides(color = guide_legend(keyheight = 1, keywidth = 0.2, title.hjust = 0.5, ncol = 3)) +
  guides(fill = guide_legend(keyheight = 1, keywidth = 0.2, title.hjust = 0.5, ncol = 3)) +
  theme(legend.position = c(0.85,0.3),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  labs(
    x = 'Longitude',
    y = 'Latitude',
    # color = 'Hydro-region',
    fill = 'Hydro-region'
  )

cluster_map_raw_HI <- cluster_map_raw + 
  scale_y_continuous(limits = c(18, 23), breaks = c(19, 22)) + 
  scale_x_continuous(limits = c(-161, -154), breaks = c(-160, -155)) +
  theme(axis.title = element_blank(),
        plot.background = element_blank(),
        legend.position = 'none')

HI_inset <- annotation_custom(grob = ggplotGrob(cluster_map_raw_HI), # joins plot with centroid coordinates
                              xmin = -135.5, ymin = 20.5, 
                              xmax = -113, ymax = 33)
cluster_map_raw_wHI <- cluster_map_raw + HI_inset


#### CREATE NEW GEOM FOR PLOTTING HAWAII INSET ####
# Create new ggplot geom for adding a different inset plot for each facet
# In this case, for adding a Hawaii inset to each season facet plot

GeomCustom <- ggproto(
  "GeomCustom",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    data
  },
  
  draw_group = function(data, panel_scales, coord, xmin, xmax, ymin, ymax) {
    vp <- grid::viewport(x=data$xmin, y=data$ymin, 
                         width = data$xmax - data$xmin, height = data$ymax - data$ymin,
                         just = c(0,0))
    g <- grid::editGrob(data$grob[[1]], vp=vp)
    ggplot2:::ggname("geom_custom", g)
  },
  
  required_aes = c("grob","xmin","xmax","ymin","ymax")
  # required_aes = c("grob","xmin","xmax","ymin","ymax")
  
)

geom_custom <-  function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )}

#### MAP ANNUAL HIGH AND LOW FLOW TRENDS ####

# 1. Select rows for a given start year, recurrence that have significant cluster trend
# This will eliminate need to subset, can use facet instead
# 2. For each timeframe_sel (season) + event_type (high-flow/low-flow) combination:
# a) Identify hydro-regions (clusters) with significant trends
# b) Use those clusters to select cluster polygons, append corresponding timeframe & event_type columns
# 3. Make a function with three inputs: 
#       timeframe_sel (to be changed in one run. options: 1) 'Annual' and 2) c('SON','DJF','MAM','JJA'))
#       event magnitude (n_peaks_sel) and start_year (unchanged for a given run)
# a) Function creates a map of trends for that season, data for ggplot is entire dataframe of Z-Scores
# b) Facet ggplot by timeframe_sel (season) and event_type (high-flow/low-flow)
# Plot annual trends in high-flow and low flow
getQtrend_plot <- function(start_year_sel, recurrence_sel){
  stations_sig <- station_change_master[start_year == start_year_sel & 
                                                   # timeframe_sel == season_sel & 
                                                   n_peaks_sel == recurrence_sel &
                                                   abs(z_cluster) >= 1.96][
                                                   ]
  # Take inputs of selected season, event type
  # Return polygons with that season, event type 
  # With columns identifying season and event type
  getPolys_timeframe <- function(season_sel, event_type_sel){
  cluster_polygons_sel <- cluster_polygons[cluster %in% unique(stations_sig[event_type == event_type_sel &
                                                                            timeframe_sel == season_sel]$cluster)][
                                    , ':='(event_type = event_type_sel,
                                           timeframe_sel = season_sel
                                           )]
  return(cluster_polygons_sel)
  }
  
  cluster_polygons_sel <- rbind(rbindlist(lapply(timeframe, getPolys_timeframe, 'High-flow')),
                          rbindlist(lapply(timeframe, getPolys_timeframe, 'Low-flow')))
  
  getCluster_trend_map <- function(timeframe_sel_plot, event_type_sel){
    stations_sig_plot <- stations_sig[timeframe_sel %chin% timeframe_sel_plot & event_type == event_type_sel][
      z*z_cluster < 0, cluster_color:= 'Not sign. (p > 0.05)' # set stations with trend opposite of cluster to grey
    ][
      ,alpha:=ifelse(abs(z) > 1.96, 1, 0.5)
    ][order(alpha)]
    
    cluster_trend_map <- ggplot(cluster_polygons_sel[timeframe_sel %chin% timeframe_sel_plot & event_type == event_type_sel], 
                                aes(x = long, y = lat, 
                                    # color = as.character(cluster), fill = as.character(cluster), 
                                    group = cluster)) + 
      geom_map(data = us_ca, map = us_ca, aes(map_id = id, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) +
      geom_polygon(fill = NA, color = 'black', lty = 'dashed', size = 0.25) + # Add bounding polygon for each cluster
      geom_point(data = stations_sig_plot, 
                 color = 'grey30', aes(fill = cluster_color, alpha = alpha), pch = 21, stroke = 0.2) +
      season_facet +
      scale_alpha_identity() +
      # scale_color_manual(values = cl_colors) +
      scale_fill_manual(values = cluster_color_values) +
      scale_y_continuous(limits = c(23, 60)) + 
      scale_x_continuous(limits = c(-131.2, -54)) +
      # guides(color = guide_legend(keyheight = 1, keywidth = 0.2, title.hjust = 0.5, ncol = 3)) +
      guides(fill = guide_legend(keyheight = 1, keywidth = 0.2, title.hjust = 0.5, ncol = 1)) +
      theme(
            # legend.position = c(0.95,0.3),
            strip.text = element_text(size = 10),
            legend.position = 'none',
            legend.background = element_blank(),
            legend.key = element_blank()
            ) +
      labs(
        x = 'Longitude',
        y = 'Latitude',
        # color = 'Hydro-region',
        fill = 'Trend'
      )
    
    cluster_trend_map_HI <- lapply(timeframe_sel_plot, 
      function(timeframe_sel_indiv){
        return(
          ggplotGrob(ggplot(cluster_polygons_sel[timeframe_sel == timeframe_sel_indiv & event_type == event_type_sel], 
                      aes(x = long, y = lat, 
                          # color = as.character(cluster), fill = as.character(cluster), 
                          group = cluster)) + 
                 geom_map(data = us_ca, map = us_ca, aes(map_id = id, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) +
                 geom_polygon(fill = NA, color = 'black', lty = 'dashed', size = 0.25) + # Add bounding polygon for each cluster
                 geom_point(data = stations_sig_plot[timeframe_sel %chin% timeframe_sel_indiv], 
                            color = 'grey30', aes(fill = cluster_color, alpha = alpha), pch = 21, stroke = 0.2) +
                 season_facet +
                   scale_alpha_identity() +
                 # scale_color_manual(values = cl_colors) +
                 scale_fill_manual(values = cluster_color_values) +
                 scale_y_continuous(limits = c(18, 23), breaks = c(19, 22)) + 
                 scale_x_continuous(limits = c(-161, -154), breaks = c(-160, -155)) +
                 theme(axis.title = element_blank(),
                       plot.background = element_blank(),
                       legend.position = 'none',
                       strip.text = element_blank(),
                       strip.background = element_blank()
                 ))
        )
      })
          
     HI_inset_list <- tibble(timeframe_sel=timeframe_sel_plot, grob = cluster_trend_map_HI)

    cluster_map_raw_wHI <- cluster_trend_map + 
      facet_wrap(~timeframe_sel, nrow = 2) +
      geom_custom(data = HI_inset_list, aes(grob = grob),
                  xmin = 0, ymin = 0, xmax = 0.27, ymax = 0.29)

    return(cluster_map_raw_wHI)
  }
  
  cluster_high_trend_map_ann <- getCluster_trend_map('Annual','High-flow')
  cluster_high_trend_map_month <- getCluster_trend_map(c('SON','DJF','MAM','JJA'), 'High-flow')
  cluster_low_trend_map_ann <- getCluster_trend_map('Annual','Low-flow')
  cluster_low_trend_map_month <- getCluster_trend_map(c('SON','DJF','MAM','JJA'), 'Low-flow')
  cluster_high_trend_map <- ggpubr::ggarrange(cluster_high_trend_map_ann,
                                              cluster_high_trend_map_month,
                                      nrow = 2, common.legend = F, heights = c(1,1.2), labels = c('A','B'))
  
  cluster_low_trend_map <- ggpubr::ggarrange(cluster_low_trend_map_ann,
                                             cluster_low_trend_map_month, 
                                      nrow = 2, common.legend = F, heights = c(1,1.2), labels = c('C','D'))
  map_width <- 5.75
  map_height <- 7.5
  ggsave(cluster_high_trend_map, width = map_width, height = map_height, 
         filename = paste0(wd_figures,'cluster_high_trend_map_', start_year_sel, '_n', recurrence_sel,'.pdf'), useDingbats = F)
  ggsave(cluster_high_trend_map, width = map_width, height = map_height, 
         filename = paste0(wd_figures,'cluster_high_trend_map_', start_year_sel, '_n', recurrence_sel,'.png'))
  ggsave(cluster_low_trend_map, width = map_width, height = map_height, 
         filename = paste0(wd_figures,'cluster_low_trend_map_', start_year_sel, '_n', recurrence_sel,'.pdf'), useDingbats = F)
  ggsave(cluster_low_trend_map, width = map_width, height = map_height, 
         filename = paste0(wd_figures,'cluster_low_trend_map_', start_year_sel, '_n', recurrence_sel,'.png'))
}

# Save maps to disk
lapply(c(1940,1950,1960), getQtrend_plot, 2)
lapply(c(1940,1950,1960), getQtrend_plot, 3)
lapply(c(1940,1950,1960), getQtrend_plot, 5)
lapply(c(1940,1950,1960), getQtrend_plot, 10)
lapply(c(1940,1950,1960), getQtrend_plot, 25)


#### CLUSTER PARALLEL COORDINATE PLOTS ####
# select normalized max flows
parallel_coord_vars <- colnames(cluster_file)[which(grepl('_max',colnames(cluster_file)))]

# melt cluster file to select only max flow columns
parallel_coord_1 <- melt(cluster_file[site_no %in% unique(dt_z_comb_indiv_cluster_master$site_no)], 
                         measure.vars = parallel_coord_vars, 
                         id.vars = c('station_name','state','elevation','lat','long','cluster'))
# create month variables
parallel_coord_vars_1 <- gsub(x = parallel_coord_vars, pattern = '_normalized_peak_max', replacement = '')
# level month variables for water year
parallel_coord_vars_2 <- gsub(x = parallel_coord_vars_1, pattern = 'max', replacement = 'mar')[c(10:12,1:9)]

parallel_coord_1$variable <- gsub(x = parallel_coord_1$variable, pattern = '_normalized_peak_max', replacement = '')
parallel_coord_1$variable <- gsub(x = parallel_coord_1$variable, pattern = 'max', replacement = 'mar')


parallel_coord_1$variable <- factor(parallel_coord_1$variable, levels = parallel_coord_vars_2)
parallel_coord_2 <- parallel_coord_1 %>% group_by(cluster, variable) %>% summarise_at(c('value'), list(mean = mean))
parallel_coord_plot <- ggplot(parallel_coord_1 %>% filter(cluster %in% 1:12), 
                              aes(x = variable, y = value, color = as.factor(cluster))) + 
         geom_line(size = 0.25, aes(group = station_name)) +
         geom_line(data = parallel_coord_2 %>% filter(cluster %in% 1:12), # add cluster group mean line
                   aes(x = variable, y = mean, group = cluster), color = 'black', size = 0.5) +
         facet_wrap(.~paste0(cluster,"-",hydroregion_names[cluster])) + 
  theme_evan_facet + 
  # scale_x_discrete(labels = c('O','N','D','J','F','M','A','M','J','J','A','S')) +
  scale_x_discrete(breaks = c('oct','jan','apr','jul'),labels = c('Oct','Jan','Apr','Jul')) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_color_manual(values = cl_colors) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(hjust = 0, vjust = -2, size = 8),
    text = element_text(size = 9), 
    axis.text = element_text(size = 8), 
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.spacing.y = unit(-0.5, 'lines'), 
    legend.position = 'none',
    plot.margin = unit(c(0,5,0,5), "pt")
  ) +
  labs(
    x = "",
    y = "Discharge, normalized \n by annual maximum"
  )

# Tukey HSD box plots for visualization of seasonality
# and identification of typical high- and low-flow months
# parallel_coord_tukey_plot <- ggplot(parallel_coord_1 %>% filter(Cluster_15 %in% z_cl_sel),
parallel_coord_tukey_plot <- ggplot(parallel_coord_1,
                              aes(x = variable, y = value, fill = as.factor(cluster))) +
  geom_boxplot(aes(x = variable, y = value, group = variable), outlier.shape = NA, color = 'black', size = 0.2) +
  # geom_point(size = 0.25, aes(group = station_name)) +
  # geom_point(data = parallel_coord_2 %>% filter(cluster %in% z_cl_sel), # add cluster group mean line
  #           aes(x = variable, y = mean, group = cluster), color = 'black', size = 1) +
  facet_wrap(.~cluster) + theme_evan_facet +
  # scale_x_discrete(labels = c('O','N','D','J','F','M','A','M','J','J','A','S')) +
  scale_x_discrete(breaks = c('oct','jan','apr','jul'),labels = c('Oct','Jan','Apr','Jul')) +
  scale_y_continuous(breaks = c(0,0.5,1), expand = expand_scale(add = c(0.25, 0.25))) +
  # scale_fill_manual(values = cl_colors) + 
  scale_fill_manual(values = cl_colors[replace_cl_nos]) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(hjust = 0, vjust = -2, size = 8),
    text = element_text(size = 8), 
    axis.text = element_text(size = 8), 
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.spacing.y = unit(-0.5, 'lines'), 
    legend.position = 'none'
  ) +
  labs(
    x = "",
    y = "Discharge, normalized \n by annual maximum"
  )

# Arrange map of stations by cluster with average annual peak hydrograph plot
cluster_map_wParallels <- ggpubr::ggarrange(plotlist = list(cluster_map_raw_wHI,parallel_coord_plot), 
                                    heights = c(1,0.6), align = 'h',
                                    ncol = 1)
ggsave(cluster_map_wParallels, width = 5, height = 6, 
       filename = paste0(wd_figures,'cluster_map_wParallels.pdf'), useDingbats = F)
ggsave(cluster_map_wParallels, width = 5, height = 6, 
       filename = paste0(wd_figures,'cluster_map_wParallels.png'))

#### GENERATE PIE PLOTS ####

pie_size <- 4
# Change centroid position for cluster centroids that are too close to each other
changeDistance <- function(cluster_sel){
  cluster_centroid_sel <- cluster_centroids[cluster == cluster_sel]
  long_sel <- cluster_centroid_sel$long
  lat_sel <- cluster_centroid_sel$lat
  nearest <- cluster_centroids[cluster != cluster_sel][
    , distance := sqrt((lat-lat_sel)^2 + (long-long_sel)^2)
    ][order(distance)][1]
  long_near <- nearest$long
  lat_near <- nearest$lat
  
  pie_size_sel <- nearest$distance*(nearest$distance < pie_size*1.2)
  # Get magnitude of shift
  change_mag_long <- max(0,pie_size_sel - abs(long_sel - long_near))
  change_mag_lat <- max(0,pie_size_sel - abs(lat_sel - lat_near))
  # Get direction of shift
  change_dir_long <- sign(long_sel - long_near)
  change_dir_lat <- sign(lat_sel - lat_near)
  
  # Shift coordinates
  long_new <- long_sel + change_mag_long*change_dir_long
  lat_new <- lat_sel + change_mag_lat*change_dir_lat
  
  return(cluster_centroid_sel[,':='(long = long_new,
                                    lat = lat_new,
                                    distance = nearest$distance,
                                    pie_size_sel = pie_size_sel)])
}
cluster_centroids_adj <- rbindlist(lapply(cluster_list,changeDistance))[
  cluster == 1,':='(long = -126,
                    lat = 35) # replace Hawaii coordinates
  ]
cluster_centroids_adj[]

getCluster_pies <- function(cluster_sel, event_type_sel, start_year_sel, n_peaks){
  long <- cluster_centroids_adj[cluster == cluster_sel]$long
  lat <- cluster_centroids_adj[cluster == cluster_sel]$lat
  pie_plot_sel <- dt_z_comb_cluster_master[
    cluster == cluster_sel &
      n_peaks_sel == n_peaks &
      start_year == start_year_sel & 
      event_type == event_type_sel &
      timeframe_sel != 'Annual']
  return(annotation_custom(grob = ggplotGrob(
    ggplot(pie_plot_sel, 
           aes(x = timeframe_sel, fill = cluster_color, color = cluster_color)) + 
      geom_bar(width = 1, 
               color = "black",
               # color = .$color_sel,
               size = 0.25) +
      coord_polar(start = -pi/4) +
      scale_fill_manual(values = cluster_color_values) +
      # scale_color_manual(values = cl_colors)) %>% 
      scale_color_manual(values = cluster_color_values) +
      theme_void() + guides(fill = F) + theme(legend.position = "none")),
    xmin = long - pie_size, xmax = long + pie_size,
    ymin = lat - pie_size, ymax = lat + pie_size
  )
  )
}

getPie_combined_maps <- function(event_type_sel){
  cluster_pie_sel <- lapply(cluster_list, getCluster_pies,event_type_sel,start_year_sel,recurrence_sel)
  
  cluster_map_pie <- ggplot() +
    geom_map(data = us_ca, map = us_ca, aes(map_id = id, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) +
    # geom_polygon(fill = NA, lty = 'dashed') + # Add bounding polygon for each cluster
    season_facet +
    scale_y_continuous(limits = c(23, 60)) + 
    scale_x_continuous(limits = c(-131.2, -54)) +
    theme(
      # legend.position = c(0.85,0.3),
      # legend.background = element_blank(),
      # legend.key = element_blank(),
      plot.title = element_text(size = strip_text_size* 1.4)) +
    labs(
      x = 'Longitude',
      y = 'Latitude',
      title = paste0(event_type_sel, ' frequency change'),
      # color = 'Hydro-region',
      fill = 'Trend'
    )
  
  cluster_map_pie_HI <- cluster_map_pie + 
    scale_y_continuous(limits = c(18, 23), breaks = c(19, 22)) + 
    scale_x_continuous(limits = c(-161, -154), breaks = c(-160, -155)) +
    theme(axis.title = element_blank(),
          plot.background = element_blank(),
          legend.position = 'none',
          plot.title = element_blank())
  
  HI_pie_inset <- annotation_custom(grob = ggplotGrob(cluster_map_pie_HI), # joins plot with centroid coordinates
                                    xmin = -135.5, ymin = 20.5, 
                                    xmax = -113, ymax = 33)
  cluster_map_pie_wHI <- cluster_map_pie + HI_pie_inset +
    cluster_pie_sel +
    geom_text_repel(data = cluster_centroids_adj, aes(x = long, y = lat, label = cluster), 
                    color = 'black', size = 5, point.padding = NA, force = 0.1, nudge_x = -3, nudge_y = 3, segment.alpha = 0)
  
  return(cluster_map_pie_wHI)
}

pie_trend_combined_maps <- lapply(c('High-flow','Low-flow'),getPie_combined_maps)

pie_trend_combined_plot <- ggpubr::ggarrange(plotlist = pie_trend_combined_maps, ncol = 1,
                                             labels = c('A','B'))

ggsave(pie_trend_combined_plot, width = map_width, height = map_height, 
       filename = paste0(wd_figures,'pie_trend_combined_plot', start_year_sel, '_n', recurrence_sel,'.pdf'), useDingbats = F)
ggsave(pie_trend_combined_plot, width = map_width, height = map_height, 
       filename = paste0(wd_figures,'pie_trend_combined_plot', start_year_sel, '_n', recurrence_sel,'.png'))



#####
#### IN DEVELOPMENT ####

#### DOWNLOAD CLIMATE DATA (TEMPERATURE) ####

# Import all station metadata for ghcn stations (takes a few minutes, about 700,000 stations)
  ghcnd_stations_all <- ghcnd_stations()
# Filter by stations with long records extending to the present
  ghcnd_stations_long <- ghcnd_stations_all %>% subset(last_year == 2019 & first_year < 1961)

# Map long-record stations on US map
  ghcn_stations_long_map <- ggplot() +
  geom_map(data = us_ca, map = us_ca, 
           aes(map_id = id, x = long, y = lat),
           color = "grey30", fill = 'white', lwd = 0.25) +
    geom_point(data = ghcnd_stations_long, aes(x = longitude, y = latitude, color = elevation), size = 0.25) +
    season_facet + 
    scale_color_distiller(palette = 'Blues') +
    scale_y_continuous(limits = c(25, 51)) + scale_x_continuous(limits = c(-131, -68)) + 
    # guides(fill = guide_legend(ncol = 2)) +
    theme(
      legend.position = c(0.1, 0.25),
      legend.spacing.x = unit(0, 'in'),
      legend.spacing.y = unit(0, 'in'),
      legend.key = element_blank(),
      legend.key.size = unit(0.75, 'picas'),
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1),
      text = element_text(size=9),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.title = element_text(hjust = 0.5, size = 7),
      strip.text = element_text(margin = margin(0,0,1.5,0)),
      legend.background = element_rect(size = 0.25, color = 'grey20'),
      legend.margin = margin(0.02,0.02,0.02,0.02, unit = 'inches')
    ) + 
    labs(
      x = "Longitude",
      y = "Latitude",
      color = 'Elevation'
    )  

  # Save map to PDF
    setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency/temp_plots")
    ggsave(ghcn_stations_long_map, filename = 'ghcn_stations_long_map.pdf', 
           width = 3.4, height = 2.7, useDingbats = F)
    setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency")

# Select cluster stations from cluster list
    # Select only Lat/Long and site info columns
      cluster_stations_for_ghcn <- cluster_file_export[,c('lat','long','site_number')]
      colnames(cluster_stations_for_ghcn) <- c('lat','long','id')
# Import closest station within 100 km of each station in hydro-regions list
      ghcn_nearby_stns <- meteo_nearby_stations(lat_lon_df = cluster_stations_for_ghcn, 
                      lat_colname = 'lat', lon_colname = 'long',
                      station_data = ghcnd_stations_long,
                      var = 'all', year_min = 1960, year_max = 2020, radius = 100, limit = 1)
# Convert confusing nested list to a useful dataframe
      ghcn_stns_sel <- do.call(rbind, lapply(names(ghcn_nearby_stns), 
                                       function(x) data.frame(c(site_number=x, ghcn_nearby_stns[[x]]))))
# Remove NAs and duplicate stations (some gages are closest to the same meteorological station)
      ghcn_stns_sel_2 <- merge(ghcn_stns_sel, cluster_file_export, by = c('site_number')) %>% subset(!is.na(id))
      ghcn_stns_sel_3 <- ghcn_stns_sel_2[which(!duplicated(ghcn_stns_sel_2$id)),]

# Find unique cluster numbers -- need to change the cluster number later to align with hcluster analysis
      ghcn_cl_sel <- unique(ghcn_stns_sel_3$Cluster_15) %>% sort()

# Initialize list for storing imported station observation data for each cluster
ghcn_import <- list(rep(NA, 15))

# For each cluster, import all station data, add it to a dataframe. This takes a long time (~10 seconds/station)
for(i in ghcn_cl_sel[1:15]){
  # select only met stations corresponding to cluster of interest
  monitor_stations <- ghcn_stns_sel_3 %>% subset(Cluster_15 == i)
  monitor_ids <- monitor_stations$id # make vector of met station names for batch import
  ghcn_import[[i]] <- meteo_pull_monitors(monitors = monitor_ids, 
                    date_min = '1959-10-01', var = c('TMIN','TMAX','TAVG','PRCP','SNOW','SNWD')) %>% data.frame()

}


# For each cluster, make some charts/analysis of trends in temperature
for(i in ghcn_cl_sel){
  if(i %in% z_cl_sel){
    z_sel_pos <- z_cl_sel[which(i == z_cl_sel)]
    renamed_cl <- replace_cl_nos[z_sel_pos]
  # filter only observations 1960-present, add date, month, year, decade columns
  cluster_tmax_grouped <- ghcn_import[[i]] %>% subset(year(date) > 1959) %>% 
    mutate(doy = yday(date), month = month(date), year = year(date), 
           half_decade = (year(date)-year(date)%%5), decade = (year(date)-year(date)%%10)) 
  
  # Summarise by year -- average monthly temp
  monthly_tmax_ann_ts <- cluster_tmax_grouped %>% 
    group_by(month, year, decade) %>%
    summarise_at(c('tmax'), mean_na)
  # Summarise by half decade -- average monthly temp
  monthly_tmax_hdec_ts <- cluster_tmax_grouped %>% 
    group_by(month, half_decade) %>%
    summarise_at(c('tmax'), mean_na)
  # Summarise by decade -- average monthly temp
  monthly_tmax_dec_ts <- cluster_tmax_grouped %>% 
    group_by(month, decade) %>%
    summarise_at(c('tmax'), mean_na)
  # # Plot monthly average trends over years
  # ggplot(monthly_tmax_ann_ts, aes(x = year, y = tmax/10)) + 
  #   geom_point() +
  #   geom_smooth(method = 'lm') +
  #   facet_wrap(.~month, scales = 'free_y') +
  #   season_facet +
  #   labs(
  #     x = 'Year',
  #     y = 'Maximum temperature (C)'
  #   )
  # Plot monthly average trends over half decades
  half_dec_trend_plot <- ggplot(monthly_tmax_hdec_ts, aes(x = half_decade, y = tmax/10)) + 
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(.~month, scales = 'free_y') +
    season_facet +
    labs(
      x = 'Year',
      y = 'Maximum temperature (C)'
    )
  # Plot average monthly temperature graph, color by decade
  avg_daily_tmax_plot <- ggplot(monthly_tmax_dec_ts, 
         aes(x = month, y = tmax/10, color = decade, group = decade)) + 
    geom_point(size = 0.5) +
    geom_line(lwd = 0.25) +
    # facet_wrap(.~paste0(decade,'s'), nrow = 1) +
    scale_color_distiller(palette = 'Reds', direction = -1) +
    season_facet +
    theme(
      legend.position = 'right'
    ) +
    labs(
      x = 'Month',
      y = 'Average daily maximum temperature (C)',
      color = 'Decade'
    )
  # Save plots
  setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency/temp_plots")
  ggsave(half_dec_trend_plot, filename = paste0('half_dec_trend_plot_cluster', ghcn_cl_sel[renamed_cl], '.pdf'), 
         width = 5, height = 4.2, useDingbats = F)
  ggsave(avg_daily_tmax_plot, filename = paste0('avg_daily_tmax_plot_cluster', ghcn_cl_sel[renamed_cl], '.pdf'), 
         width = 3.4, height = 3, useDingbats = F)
  setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency")
}}


#### 
#####
