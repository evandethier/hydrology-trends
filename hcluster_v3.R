#### LIBRARIES ####
library(rnoaa)
library(modelr)
library(scales)
library(kdensity)
library(EnvCpt)
library(ggplot2)
library(ggthemes)
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

    # Import HCDN sites 
    station_info <- data.table(read_excel("hcdn-huc08.xlsx"))
    
    # Import stations with cluster attached
    cluster_file <- fread("cluster_data_51419.dat")[Cluster_15_dist < 8][, ':='(
      agency_cd = ifelse(nchar(site_number) == 7, 'WSC', 'USGS'))
    ]
    # [, .(agency_cd, site_number, station_name, state, lat, long, drainage_area_sqkm, eco03region, elevation, Cluster_15)]

cluster_file_export <- cluster_file # subset station characteristic file by those columns
cluster_list <- sort(unique(cluster_file$Cluster_15)) # sort by cluster
n_clusters <- length(cluster_list) # number of clusters that we will use

master_column <- 2
timeframe_select <- c("djf","mam","jja","son","annual")

n_events <- 12
cluster_z <- data.frame(matrix(ncol = 10,nrow = n_clusters)) # initialize output data file for z-scores
  colnames(cluster_z) <- c("djf_flood","mam_flood","jja_flood","son_flood","annual_flood",
                          "djf_drought","mam_drought","jja_drought","son_drought","annual_drought")
  rownames(cluster_z) <- cluster_list
smooths <- c(1,2,5,10,15,20) # smoothing windows (in years) for gaussian kernel smoothing frequency analysis

indiv_z <- data.frame(matrix(ncol = 10,nrow = nrow(cluster_file_export)))
indiv_z_cl <- data.frame(matrix(ncol = 10,nrow = nrow(cluster_file_export)))
colnames(indiv_z) <- colnames(cluster_z)
colnames(indiv_z_cl) <- paste0(colnames(cluster_z),"_cl")
station_info_export <- cbind(cluster_file_export,indiv_z,indiv_z_cl)
infocols <- ncol(cluster_file_export)

#### IMPORT ALL DATA -- RUN ONCE ONLY ####
# commented out to avoid running
cluster_files_master <- vector("list", length(min(cluster_list):max(cluster_list)))

cluster_files <- cluster_file[, .(file_name = ifelse(agency_cd == 'WSC', 
                                                    paste0('Q_daily/',site_number,'_Daily_Flow_ts.csv'),
                                                    paste0('Q_daily/', site_number,'.dat')))
                                                    ]

file_names_all <- list.files(paste0(wd_imports,'/Q_daily/')) 

for(n_cluster in 1:length(cluster_list)){
  cluster <- cluster_list[n_cluster]
  cluster_sites <- cluster_file[Cluster_15==cluster]
  cluster_stations <- cluster_sites$site_number
  cluster_files <- cluster_sites[, .(file_name = ifelse(agency_cd == 'WSC', 
                                                       paste0(site_number,'_Daily_Flow_ts.csv'),
                                                       paste0(site_number,'.dat')))]$file_name
  
#   cluster_files <- createFileList_cluster(cluster_sites)
# 
  file_list <- paste0('Q_daily/',cluster_files[cluster_files %in% file_names_all])
  print(length(file_list))
  if(length(file_list)>2){
    print(cluster)
#     # df_chosen <- cluster_sites
# 
    cluster_files_temp <- lapply(file_list, getQ_ts)
#     cluster_files_temp <- lapply(cluster_files_temp,getDates)
    cluster_files_temp <- lapply(cluster_files_temp,getPeaks)
# 
    cluster_files_master[[cluster]] <- cluster_files_temp
  }}
#   

    # Import and standardize Discharge timeseries from USGS and WSC
    getQ_ts <- function(file_name){
      if(grepl('.dat',file_name)){
        Q_ts_import <- fread(file_name)[,':='(
          Q_cms = q_cfs * 0.0283168,
          Date = ymd(Date),
          q_cfs = NULL,
          V1 = NULL,
          `qual code` = NULL
        )][,.(agency_cd,site_no, Date, water_year,Q_cms)]
      }else{
        Q_ts_import <- fread(file_name)[,':='(
          agency_cd = 'WSC',
          site_no = ID,
          Date = ymd(Date),
          Q_cms = Flow,
          ID = NULL,
          SYM = NULL,
          PARAM = NULL,
          water_year = ifelse(month(Date) > 9, as.numeric(year(Date) + 1), as.numeric(year(Date)))
        )
        ][,.(agency_cd,site_no, Date, water_year,Q_cms)]
      }
      return(Q_ts_import)
    }
    
    dt <- cluster_files_temp[[2]]
    
    getPeaks <- function(dt){
      set.seed(0)
      # Get all 7-day maxima and minima. Will need to take the extreme X% later
      return(dt[,month:=month(Date)][
                peaks(dt$Q_cms, span = 7), Q_peaks := Q_cms][
                peaks(-dt$Q_cms, span = 7), Q_lows := Q_cms]) 
    }
    
    dt_peaks <- getPeaks(dt)
    
    # Return top N event rows for all data.tables in a list
    getExtremes_high <- function(dt,n_peaks, months, start_year, end_year){
      # Select top n events by given months
      return(dt[,month:=month(Date)][
        water_year > start_year & water_year <= end_year & month %in% months][
                order(Q_peaks*-1)][
                  1:n_peaks]) # select top n_peaks events
    }
    
    # Return top N event rows for all data.tables in a list
    # Iterates over all given month categories
    getExtremes_high <- function(dt,n_peaks, months, start_year, end_year){
      # Select top n events by given months
      dt[,month:=month(Date)]
      return(rbindlist(lapply(months, function(months, dt){
        return(dt[
        water_year > start_year & water_year <= end_year & month %in% months][
                order(Q_peaks*-1)][
                  1:n_peaks])}, dt = dt))) # select top n_peaks events
    }
    
    # Return top N event rows for all data.tables in a list
    # Iterates over all given month categories
    # Also iterates over all n event categories
    timeframe <- c('DJF','MAM','JJA','SON','Annual')
    getExtremes_high <- function(dt,n_peaks, months, start_year, end_year){
      # Select top n events by given months
      dt[,month:=month(Date)]
      return(rbindlist(lapply(months, function(months){
        month_diagnostic <- ifelse(length(months) == 12, 5, round(min(months)/3,0)+1)
        month_sel <- timeframe[month_diagnostic]
        return(rbindlist(lapply(n_peaks, function(n_peaks, dt){
          return(dt[,':='(
                      n_high_peaks = round((end_year - start_year)/n_peaks),
                      timeframe_high = month_sel)][
        water_year > start_year & water_year <= end_year & month %in% months][
                order(Q_peaks*-1)][
                  1:n_peaks])}, dt = dt)))}))) # select top n_peaks events
    }
    
    getExtremes_low <- function(dt, start_year){
      return(dt[,month:=month(Date)][
        water_year > start_year & month == 12])
    }
    # For each cluster loop through all sites
    # Eliminate sites with incomplete records
    # For each site loop through each season for both floods and droughts
    # Each result returns a data.table of dates and magnitudes
    # Bind all 
    
    dt_peaks <- lapply(cluster_files_master[[i]], getPeaks)
    
    
    dt <- dt_peaks[[1]]
    
    # Get dates of continuous daily discharge sampling at each site
    # Return start year, end year, number of years
    getContinuous <- function(dt){
      last_missing <- max(which(diff(dt$Date)>(365*2)),1)
      startYear <- year(dt$Date[last_missing])
      endYear <- year(max(dt$Date))
      record_yrs <- endYear-startYear
      return(data.table(site_no = dt$site_no[1], startYear = startYear, endYear = endYear, record_years = record_yrs))
    }
    
    dt_startDate <- rbindlist(lapply(dt_peaks, getContinuous))
    
    # Set record length parameters
    start_years <- c(192:198)*10
    for(k in 1:7){
      # Test a range of start years
      start_year <- start_years[k]
      end_year = 2016
      start_year_rows <- which(dt_startDate$startYear <= 1940)
      n_stations <- length(start_year_rows)
      timeframes <- list(c(5:7), c(9:11), c(12,1,2), c(3:5), c(1:12))
      n_peaks <- round((end_year - start_year)/c(2,3,5,10,25))
      
      dt_peaks_byStartYr <- dt_peaks[start_year_rows]
      
      dt_extremes <- rbindlist(lapply(dt_peaks_byStartYr, getExtremes_high, n_peaks = n_peaks, 
                                      months = timeframes, 
                                      start_year = start_year, end_year = end_year))
        
      # Get peaks over threshold for a given N, but vary the N
      dt_count <- dt_extremes[, .(count = .N), by = .(water_year, n_high_peaks, timeframe_high)]
    
      Q_high_trends_ann <- ggplot(dt_count[timeframe_high == 'Annual'], 
             aes(x = water_year, y = count, color = factor(n_high_peaks))) + 
        # geom_point() +
        facet_wrap(~factor(timeframe_high, levels = timeframe)) +
        geom_smooth(se = F, lwd = 0.25) +
        season_facet + 
        scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
        theme(legend.position = 'right') +
        labs(x = 'Water year',
             y = 'N events/yr',
             color = paste0('Avg. years/event\nat-a-station\n(N = ', n_stations,' stns)'))
      
      Q_high_trends_month <- ggplot(dt_count[timeframe_high != 'Annual'], 
             aes(x = water_year, y = count, color = factor(n_high_peaks))) + 
        # geom_point() +
        facet_wrap(~factor(timeframe_high, levels = timeframe)) +
        geom_smooth(se = F, lwd = 0.25) +
        season_facet + 
        scale_color_manual(values = rev(c('#025159','#03A696','#F28705','#F25D27','#F20505'))) +
        theme(legend.position = 'right') +
        labs(x = 'Water year',
             y = 'N events/yr',
             color = paste0('Avg. years/event\nat-a-station\n(N = ', n_stations,' stns)'))
      
      Q_high_trends_comb <- ggarrange(Q_high_trends_ann, Q_high_trends_month, nrow = 2, common.legend = T)
      ggsave(Q_high_trends_comb, width = 3.5, height = 5, 
             filename = paste0(wd_exports,'Q_trends_sensitivity_cl',n_cluster,'_',start_year,'.pdf'), useDingbats = F)
      ggsave(Q_high_trends_comb, width = 3.5, height = 5, 
             filename = paste0(wd_exports,'Q_trends_sensitivity_cl',n_cluster,'_',start_year,'.png'))
      
    }
    
  


# # Import all files at once from a folder (makes a single data.table)
# importFiles <- function(path, pattern) {
#   files = list.files(path, pattern, full.names = TRUE) 
#   rbindlist(lapply(files, function(x) fread(x)))
# }
# 
# list.files(paste0(wd_imports,'/Q_daily/'), pattern = "*.dat", full.names = TRUE)
# 
# allUSGS_files <- importFiles(paste0(wd_imports,'/Q_daily/'), pattern = "*.dat")[,':='(
#   Q_cms = q_cfs * 0.0283168,
#   q_cfs = NULL,
#   V1 = NULL,
#   `qual code` = NULL
# )]
# 
# allWSC_files <- importFiles(paste0(wd_imports,'/Q_daily/'), pattern = "*.csv")[,':='(
#   agency_cd = 'WSC',
#   site_no = ID,
#   Q_cms = Flow,
#   ID = NULL,
#   SYM = NULL,
#   PARAM = NULL,
#   water_year = ifelse(month(Date) > 9, as.numeric(year(Date) + 1), as.numeric(year(Date)))
# )
# ]


#### GENERATE INDIV. AND CLUSTER SIGNFICANCE ####
for(cluster in cluster_list[12:length(cluster_list)]){
# for(cluster in cluster_list[2:length(cluster_list)]){
  # # for(cluster in 6:6){
    # cluster <- 6
  
  cluster_sites <- subset(cluster_file, cluster_file$Cluster_15==cluster)
  cluster_stations <- cluster_sites$site_number
  cluster_files <- createFileList_cluster(cluster_sites)
  
  file_list <- cluster_files
  print(length(file_list))
  if(length(file_list)>2){
    print(cluster)
  df_chosen <- cluster_sites
  
  cluster_files_test <- cluster_files_master[[cluster]]
  if(!is.null(cluster_files_test)){  
  # cluster_files_test <- lapply(file_list, read.csv)
  # cluster_files_test <- lapply(cluster_files_test,getDates)
  # cluster_files_test <- lapply(cluster_files_test,getPeaks)
    ###### CHOOSE PLOT TO INITIALIZE (MAX, TOP5 AVG) #########
  
  ### CURRENT USEFUL PARAMETERS ###
  year_start <- 1960
  year_end <- 2019
  bounds <- (c(year_start,year_end))
  region_events <- data.frame(matrix(ncol = (length(file_list)+1), nrow = length(year_start:year_end)))
  colnames(region_events) <- c("water_year",matrix(df_chosen$site_number))
  
  region_events$water_year <- c(year_start:year_end)
  
  t_list_fl <- replicate(n = 5,
            expr = {region_events},
            simplify = F)
  names(t_list_fl) <- timeframe_select
  t_list_dr <- replicate(n = 5,
            expr = {region_events},
            simplify = F)
  names(t_list_dr) <- timeframe_select
  
  good_data <- rep(0,length(file_list))
  
  #initialize for-loop
  
  # setwd("~/Documents/Dartmouth/Discharge Project/R files/station_data_all")
  setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency")
  
  #### DO I MAKE GET PEAKS RETURN NAS??????######
  # for(file in file_list){
  for(i in 1:length(file_list)){
  
    print(file_list[i])
    # #extract full data file for specific site
    # station_df <- read.csv(file)
    # 
    # #adjust dates and get normalized discharge for each station
    # station_df <- getDates(station_df)
    # # station_df <- getNormDischarge(station_df, station_info)
    # station_df <- getPeaks(station_df) # set Q to 1 cfs for all but peak days
    # 
    if(grepl("csv",file_list[i])){ # determines whether a station is canadian or usgs by filename
      canada_yes_no <- "canada"
    }else{
      canada_yes_no <- "usa"
    }
    station_df <- cluster_files_test[[i]]
    print(typeof(station_df[1]))
  # }}}}#all curly brackets in this line are test
    # there was some issue in importing data that caused conversion to list, need to correct by switching back
    if(typeof(station_df$site_no)=="integer"){station_df$site_no <- as.character(station_df$site_no)} 
    if(typeof(station_df$datetime)=="list"){station_df$datetime <- as.POSIXct(as.Date(station_df$datetime))} 
    if(typeof(station_df$Date)=="list"){station_df$Date <- as.POSIXct(as.Date(station_df$Date))} 
    
    station_df <- station_df[which(!is.na(station_df$q_cfs)),] # remove NA rows
    if(min(diff(station_df$datetime))==1){
    last_missing <- which(diff(station_df$datetime)>(365*2))}else{ # find gaps in the timeseries (>2 yrs missing)
      last_missing <- which(diff(station_df$datetime)>(365*2*24))
      }
    if(!is_empty(last_missing)){ # if there is a 2-yr gap, start timeseries after last major gap
       station_df <- station_df[(max(last_missing)+1):nrow(station_df),]}
    
    # start timeseries at 1960 and only use if there is at least some data for every month
    if(any(station_df$water_year==1960 & length(unique(station_df$month))==12)){
       station_df <- station_df[min(which(station_df$water_year == 1960)):nrow(station_df),]
    
    
    # optional filtering by season
    
    seasonality_options <- c(c(11,0,1),c(2,3,4),c(5,6,7),c(8,9,10), c(1,1,1))
    for(season_select in 1:5){
      # print(i)
      timeframe <- timeframe_select[season_select]
      # print(timeframe)
      row_add <- season_select*length(1899:2019)-121
      month_select <- seasonality_options[(season_select*3-2):((season_select*3-2)+2)]
      # select season
      if(season_select<5){station_df_season <- subset(station_df,(station_df$month==month_select[1] | station_df$month==month_select[2] | station_df$month==month_select[3]))}else{ # if we want to filter by month
        station_df_season <- station_df}
      
      # find flood and drought years
      station_t_event <- getT(station_df_season,"timeframe",n_events)
      # don't run test if data have gaps, or too many floods or droughts
      n_too_many <- n_events*1.5
      if(nrow(station_t_event)>=(2015-1961) & sum(station_t_event$drought_count)<n_too_many & sum(station_t_event$flood_count)<n_too_many){
      # print(station_t_flood)
        # print(i)
      good_data[i] <- good_data[i]+i+1
      t_list_fl[[timeframe]][which(region_events$water_year %in% station_t_event$water_year),(i+1)] <- station_t_event$flood_count
      t_list_dr[[timeframe]][which(region_events$water_year %in% station_t_event$water_year),(i+1)] <- station_t_event$drought_count
      }
      
      # compute individual site density
       station_df_kernel <- getKernelDensity(station_df_season,smooths,timeframe,"individual")
       station_df_plot <- station_df_kernel[[1]]
       if(nchar(station_df$site_no[1])==7 & canada_yes_no == "usa"){station_no_str <- paste0(0,station_df$site_no[1])}else{
         station_no_str <- station_df$site_no[1]
       }
       ggexport(station_df_plot, filename = paste(toString(station_df$site_no[1]), "_cluster_", toString(cluster), "_", timeframe_select[season_select],"_n",toString(n_events), "_drought_flood_recur_ts.pdf",sep = ""), res = 250)
       station_info_export[which(station_info_export$site_number == station_no_str),(season_select+infocols)] <- station_df_kernel[[2]]
       station_info_export[which(station_info_export$site_number == station_no_str),(season_select+infocols+5)] <- station_df_kernel[[3]]
    }
    print(station_info_export[which(station_info_export$site_number == station_no_str),]) 
     
    }
    }
  
  good_data <- c(1,(good_data[which(good_data%%5==0 & good_data!=0)]/5)) # selects only series with reliable data (>60 yrs)
  
  for(i in 1:5){
    cluster_fl <- t_list_fl[[i]][,good_data]
    cluster_dr <- t_list_dr[[i]][,good_data]
    fl_plot <- getAreaPlot(cluster_fl,"flood",bounds,timeframe_select[i])
    dr_plot <- getAreaPlot(cluster_dr,"drought",bounds,timeframe_select[i])
    
    if(!is_empty(2:ncol(cluster_fl))&ncol(cluster_fl)>2){
    fl_sum <- rowSums(cluster_fl[,2:ncol(cluster_fl)], na.rm = T)/ncol(cluster_fl[,2:ncol(cluster_fl)])
    dr_sum <- rowSums(cluster_dr[,2:ncol(cluster_dr)], na.rm = T)/ncol(cluster_dr[,2:ncol(cluster_dr)])
    
    fl_dr <- data.frame(cbind(cluster_fl$water_year,fl_sum,dr_sum))
    colnames(fl_dr) <- c("water_year","flood","drought")
    # print(timeframe_select[i])
    kernel_output <- getKernelDensity(fl_dr,smooths,timeframe_select[i],"cluster")
    kernel_plot <- kernel_output[[1]]
    cluster_z[which(as.numeric(rownames(cluster_z)) == cluster), i] <- kernel_output[[2]]
    cluster_z[which(as.numeric(rownames(cluster_z)) == cluster), (i+5)] <- kernel_output[[3]]
    # plot(kernel_plot)
    # area_plot_combine <- ggarrange(fl_plot, dr_plot,
    #                              heights = c(1,1),
    #                              ncol = 1, nrow = 2, align = "v", labels = c("A","B"))
    # plot(area_plot_combine)}
    # 
  
  # for plotting kernel density per cluster
  kernel_area_plot <- ggarrange(kernel_plot,
                                 # theme(axis.title.x=element_blank(),axis.text.x.bottom = element_blank()), 
                                 fl_plot, dr_plot, heights = c(1,1,1),
                                 ncol = 1, nrow = 3, align = "v", labels = c("A","B","C"))
  ggexport(kernel_area_plot, filename = paste("cluster_", toString(cluster), "_", timeframe_select[i],"_n",toString(n_events), "_drought_flood_recur_ts.pdf",sep = ""), res = 250)
  # plot(kernel_area_plot)
}
  }
  
  station_info_export[which(as.numeric(station_info_export$Cluster_15) == cluster),paste0(colnames(cluster_z),"_cl")]<- cluster_z[which(as.numeric(rownames(cluster_z)) == cluster),]
  print(cluster_z)
  }}
}

#### CLEAN UP AND EXPORT CLUSTER TREND FILES ####
# after-the-fact matching of each station to the cluster change info
for(cluster in cluster_list){
  station_info_export[which(as.numeric(station_info_export$Cluster_15) == cluster),
                    paste0(colnames(cluster_z),"_cl")] <- cluster_z[
                      which(as.numeric(rownames(cluster_z)) == cluster),]}

  write.csv(cluster_z,file = paste("cluster_z_n",toString(n_events), "_us_ca_drought_flood_recur.csv",sep = ""))
  write.csv(station_info_export,file = paste("indiv_z_n",toString(n_events), "_us_ca_drought_flood_recur.csv",sep = ""))

  n_events
  setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency")
   # manually import station and cluster change master files
  cluster_z <- read_csv('cluster_z_n12_us_ca_drought_flood_recur.csv')
  cluster_z <- cluster_z[,c("djf_flood","mam_flood","jja_flood","son_flood","annual_flood",
    "djf_drought","mam_drought","jja_drought","son_drought","annual_drought")]
  station_change_master <- read_csv('indiv_z_n12_us_ca_drought_flood_recur.csv')
  station_change_master <- station_change_master[,colnames(cbind(cluster_file_export,indiv_z,indiv_z_cl))]
# station_change_master_n20 <- station_change_master  # write n20 to more permanent df
  
station_change_master <- station_info_export
station_change_master$Cluster_15 <- as.factor(station_change_master$Cluster_15)

# CHECK HERE MAKE SURE CLUSTER REPLACEMENT IS CONSISTENT!!! #
# replace cluster numbers with new numbers for better ordering
replace_cl_nos <- c(4,3,9,10,11,6,5,12,7,2,8,1) # ordered list of replacement

z_cl_sel <- c(1,2,3,4,5,6,8,9,10,11,12,14)
# Remove NA rows and change row order of cluster_z
cluster_z <- cluster_z[z_cl_sel,][replace_cl_nos,]

station_change_master$Cluster_12 <- NA # make a new cluster no. column.
station_change_master <- station_change_master[which(station_change_master$Cluster_15 %in% z_cl_sel),]
for(i in 1:12){
  cl_no <- c(1,2,3,4,5,6,8,9,10,11,12,14)[i]
  rows <- which(station_change_master$Cluster_15 == cl_no)
  station_change_master$Cluster_12[rows] <- replace_cl_nos[i] # replace rows with new cl. number
}
station_change_master$Cluster_15 <- station_change_master$Cluster_12
station_change_master <- station_change_master[,colnames(station_change_master)!='Cluster_12']

cluster_file_2 <- cluster_file[which(cluster_file$Cluster_15 %in% z_cl_sel),]
cluster_file_2$Cluster_12 <- NA
for(i in 1:12){
  cl_no <- c(1,2,3,4,5,6,8,9,10,11,12,14)[i]
  print(i)
  rows <- which(cluster_file_2$Cluster_15 == cl_no)
  print(rows)
  cluster_file_2$Cluster_12[rows] <- replace_cl_nos[i] # replace rows with new cl. number
}
cluster_file_2$Cluster_15 <- cluster_file_2$Cluster_12
cluster_file_2 <- cluster_file_2[,colnames(cluster_file_2)!='Cluster_12']
unique(cluster_file_2$Cluster_15)

# new z_cl_sel
z_cl_sel_new <- 1:12


#### INITIALIZE MAP DATA FOR N.AMERICA ####
# get states shapefile for clipping/display
us_states <-  us_boundaries(map_date = NULL, type = c("state"), resolution = c("low"), states = NULL)
us_states <- us_states[us_states$state_abbr != "AK" & us_states$state_abbr != "HI" & us_states$state_abbr != "PR",]

# convert shapefile to Spatial class
us_states <- as(us_states, 'Spatial')

# generate merged shapefile of all states--each state part of one group rather than split up
us_states_merge <- gSimplify(gUnaryUnion(us_states),0.25)
IDs <- data.frame(ID=sapply(slot(us_states_merge, "polygons"), function(x) slot(x, "ID")))
rownames(IDs)  <- IDs$ID
us_states_merge <- SpatialPolygonsDataFrame(us_states_merge,IDs)

# plot(us_states)

# set projection for states shapefile
projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# generate coordinates for each station in cluster analysis
coordinates(station_change_master) = ~long+lat
# project to match states shapefile
proj4string(station_change_master) <- proj4string(us_states_merge)

setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency")
# import canadian provinces
canada_prov <- read_sf(dsn = "Canada", layer = "Canada")
# plot(canada_prov)
canada_geom <- st_geometry(canada_prov)
# attributes(canada_geom)

# do conversions and projections for canadian provinces to match us states
canada_prov <- as(canada_prov, 'Spatial')

canada_prov <- spTransform(canada_prov, projection)
# coordinates(canada_prov) <- ~long+lat
proj4string(canada_prov) <- proj4string(us_states_merge)

names(canada_prov) <- c('name','frenchName') # rename columns labeling canadian provinces
canada_prov <- canada_prov[,c('name')] # only select column with province name
us_states <- us_states[,c('name')] # only select column with province name
us_states <- rbind(us_states,canada_prov) # combine canadian and us shapefiles

plot(us_states)
# initialize plotting variable to avoid overwriting station_change_master
plt <- station_change_master

ca_us_number_stations <- plt %>% data.frame() %>% summarise_at(c('canada'), list(ca = function(x) length(which(x == 7)),
                                              us = function(x) length(which(x == 8))))
#### OUTLINE CLUSTERS ON MAP ####
poly_list = list(NA)
cl_names <- c()
for(cl in 1:length(unique(plt$Cluster_15))){
  cluster_poly <- sort(unique(plt$Cluster_15))[cl]
  cl_out_sel <- plt[plt$Cluster_15 == cluster_poly,]
  if(nrow(cl_out_sel)>2){
    cl_names[cl] <- cl
    cl_out_hull <- chull(cl_out_sel@coords)
    cl_hull_coords <- cl_out_sel[c(cl_out_hull,cl_out_hull[1]),]
    poly_list[cl] <- Polygons(list(Polygon(as.matrix(cl_hull_coords@coords))),
                              # paste0("cluster_",cluster_poly,"_poly"))
                              # cluster_poly)
                              ID = cl)
  }}
# poly_list <- poly_list[-which(sapply(poly_list,is.null))] this doesn't work anymore for some reason
cl_names <- cl_names[!is.na(cl_names)]
cl_poly_list <- SpatialPolygons(poly_list, proj4string = projection)

#### COLOR PALETTE ####
# set color palette for remaining plots
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  
pal_length <- length(cl_names)
cl_colors <- getPalette(pal_length)
names(cl_colors) <- as.numeric(names(cl_poly_list))
cl_colors['12'] <- '#cab2d6'
cl_colors['8'] <- '#a6cee3'
cl_colors['5'] <- '#1f78b4'
cl_colors['6'] <- '#419486'
cl_colors['3'] <- "#D96D3B"
cl_colors['9'] <- "#999999"
cl_colors['10'] <- "#6a3d9a"

#### CLUSTER PARALLEL COORDINATE PLOTS ####
# select normalized max flows
parallel_coord_vars <- colnames(cluster_file_2)[which(grepl('max',colnames(cluster_file_2)))][-1]

# melt cluster file to select only max flow columns
parallel_coord_1 <- melt(cluster_file_2, measure.vars = parallel_coord_vars, 
                         id.vars = c('station_name','state','elevation','lat','long','Cluster_15'))
# create month variables
parallel_coord_vars_1 <- gsub(x = parallel_coord_vars, pattern = '_normalized_peak_max', replacement = '')
# level month variables for water year
parallel_coord_vars_2 <- gsub(x = parallel_coord_vars_1, pattern = 'max', replacement = 'mar')[c(10:12,1:9)]

parallel_coord_1$variable <- gsub(x = parallel_coord_1$variable, pattern = '_normalized_peak_max', replacement = '')
parallel_coord_1$variable <- gsub(x = parallel_coord_1$variable, pattern = 'max', replacement = 'mar')


parallel_coord_1$variable <- factor(parallel_coord_1$variable, levels = parallel_coord_vars_2)
parallel_coord_2 <- parallel_coord_1 %>% group_by(Cluster_15, variable) %>% summarise_at(c('value'), list(mean = mean))
parallel_coord_plot <- ggplot(parallel_coord_1 %>% filter(Cluster_15 %in% 1:12), 
                              aes(x = variable, y = value, color = as.factor(Cluster_15))) + 
         geom_line(size = 0.25, aes(group = station_name)) +
         geom_line(data = parallel_coord_2 %>% filter(Cluster_15 %in% 1:12), # add cluster group mean line
                   aes(x = variable, y = mean, group = Cluster_15), color = 'black', size = 0.5) +
         facet_wrap(.~Cluster_15) + theme_evan_facet + 
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
    legend.position = 'none'
  ) +
  labs(
    x = "",
    y = "Discharge, normalized \n by annual maximum"
  )

# Tukey HSD box plots for visualization of seasonality
# and identification of typical high- and low-flow months
# parallel_coord_tukey_plot <- ggplot(parallel_coord_1 %>% filter(Cluster_15 %in% z_cl_sel),
parallel_coord_tukey_plot <- ggplot(parallel_coord_1,
                              aes(x = variable, y = value, fill = as.factor(Cluster_15))) +
  geom_boxplot(aes(x = variable, y = value, group = variable), outlier.shape = NA, color = 'black', size = 0.2) +
  # geom_point(size = 0.25, aes(group = station_name)) +
  # geom_point(data = parallel_coord_2 %>% filter(Cluster_15 %in% z_cl_sel), # add cluster group mean line
  #           aes(x = variable, y = mean, group = Cluster_15), color = 'black', size = 1) +
  facet_wrap(.~Cluster_15) + theme_evan_facet +
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
ggsave(parallel_coord_tukey_plot, filename = 'parallel_coord_tukey.pdf', width = 3.4, height = 3)
ggsave(parallel_coord_tukey_plot, filename = 'parallel_coord_tukey.png', width = 3.4, height = 3)

#### INDIV Z ANALYSIS ####
season_sig_summary <- data.frame('Cluster_15' = sort(unique(plt$Cluster_15)))
plt %>% data.frame() %>% group_by(Cluster_15) %>% summarise_at(c('annual_flood_cl'),list(length))
for(plot_num in 1:length(colnames(cluster_z))){
  # plot_num <- 10
  plot_timeframe <- colnames(cluster_z)[plot_num] # timeframe name selected for plot
  plot_timeframe_cl <- paste0(plot_timeframe, "_cl") # timeframe for plot, with cluster tag

  season_sig_station_indiv <- plt %>% data.frame() %>% group_by(Cluster_15) %>% 
    select(c('lat', 'long','station_name','state','drainage_area_sqkm','elevation'), # select ID rows
           value = plot_timeframe) %>% # select season & rename
    filter(abs(value) > 1.645 & !is.na(value)) # select significant rows for seasons
  
  rename_season <- c('Cluster_15',paste0(plot_timeframe,'_inc'),paste0(plot_timeframe, '_dec'))
  season_sig_station_indiv_summary <- season_sig_station_indiv %>% 
    summarise_at(c('value'), list(increase = function(x) length(which(x > 0)),
                                  decrease = function(x) length(which(x < 0))))
  colnames(season_sig_station_indiv_summary) <- rename_season
  season_sig_summary <- merge(season_sig_summary, season_sig_station_indiv_summary, all = T)
}

write_csv(season_sig_summary,path = 'season_sig_summary.csv')
#### INDIV Z CORRESPOND TO CLUSTER Z ####

  for(plot_num in 1:length(colnames(cluster_z))){
    # plot_num <- 10
    plot_timeframe <- colnames(cluster_z)[plot_num] # timeframe name selected for plot
    plot_timeframe_cl <- paste0(plot_timeframe, "_cl") # timeframe for plot, with cluster tag
    
    season_sig_stations <- plt %>% data.frame() %>% group_by(Cluster_15) %>% 
      select(c('lat', 'long','station_name','state','drainage_area_sqkm','elevation'), # select ID rows
             value = plot_timeframe, value_cl = plot_timeframe_cl) %>% # select season & rename
      filter(abs(value_cl) > 1.645 & !is.na(value)) # select significant rows for seasons
    
    same_sign_byCluster <- season_sig_stations %>% 
        mutate(same = value*value_cl) %>% # column w/ sign for same/oppos. change as cluster
        mutate(sig = as.numeric(same > 0 & abs(value) > 1.645), # column w/ 1 = significant change in same dir. as cluster
               sig_opp = as.numeric(same <= 0 & abs(value) > 1.645)) %>% # column w/ 1 = signif. change in opp. dir. as cluster
      summarise_at(c('same', 'sig','sig_opp'), list(n_true = function(x) length(which(x > 0)), # count same sign, significance, signif. reverse
                                  n_false = function(x) length(which(x <= 0)), # count oppos. sign, significance
                                  n_stations = length)) %>% # count total n
        select(c('Cluster_15','same_n_true','sig_n_true','same_n_false','sig_opp_n_true','same_n_stations')) %>% # select relevant stations
                                  mutate(percent_same = same_n_true/same_n_stations, percent_same_signif = sig_n_true/same_n_stations, 
                                         percent_opp_signif = sig_opp_n_true/same_n_stations)
    write_csv(same_sign_byCluster, path = paste0(plot_timeframe, '_same_sign.csv'))
  }
                                                   

plt %>% data.frame() %>% group_by(Cluster_15) %>% 
  filter(abs(annual_drought_cl) > 1.645 & !is.na(annual_drought)) %>% nrow()
#### PLOT INDIVIDUAL Z-SCORES BY CLUSTER ####
# melt plt dataframe for plotting
# make all z-scores as one column, with columns for site and season
plt_melt <- melt(data.frame(plt), id.vars = c('site_number','station_name','Cluster_15'), 
     measure.vars = c('djf_flood','mam_flood','jja_flood','son_flood','annual_flood',
                      'djf_drought','mam_drought','jja_drought','son_drought','annual_drought'))
# create same columns for cluster z-scores
plt_melt_cl <- melt(data.frame(plt), id.vars = c('site_number','station_name','Cluster_15'), 
                    measure.vars = c('djf_flood_cl','mam_flood_cl','jja_flood_cl','son_flood_cl','annual_flood_cl',
                                     'djf_drought_cl','mam_drought_cl','jja_drought_cl','son_drought_cl','annual_drought_cl'),
                    value.name = 'cluster_sig')

# combine indiv. and cluster z-score columns into one dataframe
plt_melt_comb <- cbind(plt_melt, plt_melt_cl[,c('cluster_sig')])

plt_melt_comb$col <- NA # new column for color
plt_levels <- c('navy','blue','black','pink','red') # set color options

# rename columns
colnames(plt_melt_comb) <- c(colnames(plt_melt), 'cluster_sig', 'col')

# assign a color value to each station based on drought/flood and z-score (drier = red, wetter = blue)
# couldn't figure out a better way to do this with cut
plt_melt_comb$col[which(grepl('drought', plt_melt_comb$variable) == T & plt_melt_comb$cluster_sig >= 1.645|
                    grepl('drought', plt_melt_comb$variable) == F & plt_melt_comb$cluster_sig <= -1.645)] <- 'pink'
plt_melt_comb$col[which(grepl('drought', plt_melt_comb$variable) == T & plt_melt_comb$cluster_sig >= 1.96|
                    grepl('drought', plt_melt_comb$variable) == F & plt_melt_comb$cluster_sig <= -1.96)] <- 'red'
plt_melt_comb$col[which(grepl('drought', plt_melt_comb$variable) == T & plt_melt_comb$cluster_sig <= -1.645|
                    grepl('drought', plt_melt_comb$variable) == F & plt_melt_comb$cluster_sig >= -1.645)] <- 'light blue'
plt_melt_comb$col[which(grepl('drought', plt_melt_comb$variable) == T & plt_melt_comb$cluster_sig <= -1.96|
                    grepl('drought', plt_melt_comb$variable) == F & plt_melt_comb$cluster_sig >= 1.96)] <- 'blue'

plt_melt_comb$col[which(abs(plt_melt_comb$cluster_sig)<1.645)] <- 'black' # insig. change = black

# categorize by flood and drought in new column
plt_melt_comb$flood_drought[which(grepl('drought',plt_melt_comb$variable))] <- 'Drought'
plt_melt_comb$flood_drought[which(grepl('flood',plt_melt_comb$variable))] <- 'Flood'
plt_melt_comb$flood_drought <- as.factor(plt_melt_comb$flood_drought)

# levels with spaces to distinguish flood from drought categories
levels(plt_melt_comb$variable) <- c(' DJF',' MAM',' JJA',' SON',' Annual',
                                    'DJF','MAM','JJA','SON','Annual')

# levels without spaces
levels(plt_melt_comb$variable) <- c('DJF','MAM','JJA','SON','Annual',
                                    'DJF','MAM','JJA','SON','Annual')

# set colors and color names for matching with flood_drought column values
plt_cols <- c('#a50f15','#de2d26','grey40','#3182bd','#08519c')
names(plt_cols) <- c('red','pink','black','light blue','blue')
# # select clusters with sufficient data for analysis (remove cluster 7,13,15)
# z_cl_sel <- sort(unique(data.frame(plt_melt_comb)[which(!is.na(plt_melt_comb$cluster_sig) & plt_melt_comb$Cluster_15 != 13),]$Cluster_15))

z_cl_sel <- 1:12
# function for getting z-score plots for all individual stations, divided by clusters
getZplot <- function(z_cl_sel){
indiv_z_plot <- ggplot(data = data.frame(plt_melt_comb) %>% filter(Cluster_15 %in% z_cl_sel)) + 
  geom_vline(aes(xintercept = -1.96), lwd = 0.2) + # add significance lines
  geom_vline(aes(xintercept = 1.96), lwd = 0.2) + # add significance lines
  geom_vline(aes(xintercept = 1.645), lwd = 0.2) + # add significance lines
  geom_vline(aes(xintercept = -1.645), lwd = 0.2) + # add significance lines
  geom_vline(aes(xintercept = 0), lwd = 0.2) + # add zero line
  geom_jitter(aes(y = variable, x = value, # jitter points to avoid overplotting
                  color = col),
                  height = 0.15, size = 0.2) + 
  scale_color_manual(values = plt_cols) + # color by flood/drought and signficance
  scale_x_continuous(limits = c(-6,6)) + 
  # facet_wrap(~Cluster_15, drop = T) + theme_evan_facet + 
  facet_grid(flood_drought~Cluster_15, drop = T, switch = 'y') + # grid by flood/drought, cluster
  theme_evan_facet + 
  theme(
    legend.position = 'none',
    strip.placement = 'outside',
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.background = element_rect(size = 0.5),
    strip.text.x = element_text(hjust = 0, vjust = -2, size = 6,margin = margin(0,1,2,1)),
    strip.text.y = element_text(size = 6, margin = margin(0,1,2,1))
    # strip.text.x = element_blank()
  ) +
  labs(
    x = 'Z-Score',
    # y = 'Season'
    y = ''
  )
return(indiv_z_plot)
}

# collect plots in stacked grid, annotate axis
indiv_z_plots <- annotate_figure(
  ggarrange(plotlist = list(getZplot(z_cl_sel[1:4])+theme(axis.text.x = element_blank()),
                            getZplot(z_cl_sel[5:8])+theme(axis.text.x = element_blank()),
                            getZplot(z_cl_sel[9:12])), 
            ncol = 1, heights = c(1,1,1.1)),
  bottom = 'Z-Score', fig.lab.size = 8)

ggsave(indiv_z_plots, filename = 'indiv_z_plots.pdf', useDingbats = F, height = 5, width = 3.4)
ggsave(indiv_z_plots, filename = 'indiv_z_plots.png', height = 5, width = 3.4)
ggsave(indiv_z_plots, filename = 'indiv_z_plots_nolab.pdf', useDingbats = F, height = 5, width = 3.4)
ggsave(indiv_z_plots, filename = 'indiv_z_plots_nolab.png', height = 5, width = 3.4)

#### PLOT INDIV Z-SCORES ON N.AMERICA MAP; COLOR BY SIGNIFICANCE ####
setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/lydia_frequency/hcluster_figures")

sp_names <- colnames(plt@data)
station_elev <- plt[,colnames(plt@data)[6]] # elevation specific plot

timeframe_plots <- list()
# generate maps of north america for high- and low-flow changes for each timeframe (seasons & annual)
for(plot_num in 1:length(colnames(cluster_z))){
  # plot_num <- 10
  plot_timeframe <- colnames(cluster_z)[plot_num] # timeframe name selected for plot
  plot_timeframe_cl <- paste0(plot_timeframe, "_cl") # timeframe for plot, with cluster tag
  stn_plt <- plt[,c(plot_timeframe,plot_timeframe_cl, "Cluster_15")] # select timeframe columns
  if(!is_empty(stn_plt[which(abs(stn_plt@data[2])>1.645),])){
  stn_plt <- stn_plt[which(abs(stn_plt@data[2])>1.645),] # select only significant cluster change rows
  
  if(grepl("flood", plot_timeframe, fixed = T)){
    color1 <- "red" # if flood, decrease is red (color1), increase is blue (color2)
    color2 <- "blue"
    color_levels <- c(color2, color1)}else{ 
    color1 <- "blue" # if drought, decrease is blue (color1), increase is red (color2)
    color2 <- "red"
    color_levels <- c(color1, color2)
    }
   
  stn_plt$color <- "white"
  stn_plt$color[stn_plt@data[,2]<0] = color1 # if z-score negative, decreasing color (color1)
  stn_plt$color[stn_plt@data[,2]>0] = color2 # if z-score negative, increasing color (color2)

  # if z-score individual is different sign than z-score cluster, color is grey
  stn_plt$color[stn_plt@data[,1]*stn_plt@data[,2]<0] <- "grey"
  
  sig_cls <- unique(stn_plt$Cluster_15)
  
  sig_yes_up_subset <- stn_plt[which(stn_plt@data[,1] > 1.645 & is.na(stn_plt@data[,1])==FALSE & stn_plt@data[,2] > 1.645),]
  sig_yes_down_subset <- stn_plt[which(stn_plt@data[,1] < -1.645 & is.na(stn_plt@data[,1])==FALSE & stn_plt@data[,2] < -1.645),]
  sig_no_sign <- stn_plt[which(abs(stn_plt@data[,1]) < 1.645 & is.na(stn_plt@data[,1])==FALSE & stn_plt@data[,1]*stn_plt@data[,2]>0),]
  sig_no <- stn_plt[which(stn_plt@data[,1]*stn_plt@data[,2]<0 & is.na(stn_plt@data[,1])==FALSE),]
  sig_na <- stn_plt[which(is.na(stn_plt@data[,1])==TRUE),]
 
  cl_poly_sel <- tidy(cl_poly_list[as.character(sig_cls)])

  # determine plotting parameters (whether to plot, color, opacity) by cluster signficance, drought/flood
  # this set has outlines -- commented out for easier illustrator import
  # pl <- list(NA,NA,NA,NA,NA)
  # if(!is_empty(sig_na)){
  #   pl[[1]] <- geom_point(data = data.frame(sig_na), aes(x = long, y = lat), pch = 21, color = 'grey50', fill = 'grey70')}
  # if(!is_empty(sig_no)){
  #   pl[[2]] <- geom_point(data = data.frame(sig_no), aes(x = long, y = lat), pch = 21, color = 'grey50', fill ='grey70')}
  # if(!is_empty(sig_no_sign)){
  #   pl[[3]] <- geom_point(data = data.frame(sig_no_sign), aes(x = long, y = lat, fill = color), alpha = 0.6, pch = 21, color = 'black')}
  # if(!is_empty(sig_yes_down_subset)){
  #   pl[[4]] <- geom_point(data = data.frame(sig_yes_down_subset), aes(x = long, y = lat), pch = 21, color = 'black', fill = color1)}
  # if(!is_empty(sig_yes_up_subset)){
  #   pl[[5]] <- geom_point(data = data.frame(sig_yes_up_subset), aes(x = long, y = lat), pch = 21, color = 'black', fill = color2)}
pt_size <- 0.75
pl <- list(NA,NA,NA,NA,NA)
  if(!is_empty(sig_na)){
    pl[[1]] <- geom_point(data = data.frame(sig_na), aes(x = long, y = lat), pch = 16, color = 'grey70', size = pt_size)}
  if(!is_empty(sig_no)){
    pl[[2]] <- geom_point(data = data.frame(sig_no), aes(x = long, y = lat), pch = 16, color = 'grey70', size = pt_size)}
  if(!is_empty(sig_no_sign)){
    pl[[3]] <- geom_point(data = data.frame(sig_no_sign), aes(x = long, y = lat, color = color), alpha = 0.6, pch = 16, size = pt_size)}
  if(!is_empty(sig_yes_down_subset)){
    pl[[4]] <- geom_point(data = data.frame(sig_yes_down_subset), aes(x = long, y = lat), pch = 16, color = color1, size = pt_size)}
  if(!is_empty(sig_yes_up_subset)){
    pl[[5]] <- geom_point(data = data.frame(sig_yes_up_subset), aes(x = long, y = lat), pch = 16, color = color2, size = pt_size)}
  
  # only select plot grobs that have data in them
  pl_add <- pl[which(!c(is_empty(sig_na),is_empty(sig_no),is_empty(sig_no_sign),is_empty(sig_yes_down_subset),is_empty(sig_yes_up_subset)))]
    
  
  # timeframe_plot <- ggplot(data = data.frame(stn_plt), aes(x = long, y = lat)) + 
  timeframe_plot <- ggplot() + 
    scale_color_manual(values = color_levels) + # if pch = 16
    # scale_fill_manual(values = color_levels) + # if pch = 21
    geom_map(data = tidy(us_states, region = 'name'), map = tidy(us_states, region = 'name'), 
             aes(x = long, y = lat, map_id = id),color = 'grey40', fill = NA, lwd = 0.25) +
    pl_add + 
    theme_evan_facet +
    geom_map(data = cl_poly_sel,map = cl_poly_sel,
             aes(map_id = id, x = long, y = lat), lwd = 0.25, linetype = 'dashed', fill = NA, color = 'black') +
    scale_x_continuous(limits = c(-131.2, -54), breaks = c(-60,-80,-100,-120)) + scale_y_continuous(limits = c(25, 61)) + scale_alpha_identity() + 
    theme(legend.position = 'none') + 
    labs(
      x = 'Longitude',
      y = 'Latitude'
    )
  
  ggsave(timeframe_plot, filename = paste0(plot_timeframe,"_n_us_ca", n_events, "_map.pdf"), width = 4, height = 3)
  
  }else{
    # if there were no significant stations, plot just the state/province outlines (blank map)
  timeframe_plot <- ggplot() + 
    geom_map(data = tidy(us_states, region = 'name'), map = tidy(us_states, region = 'name'), 
             aes(x = long, y = lat, map_id = id),color = 'grey40', fill = NA, lwd = 0.25) +
    theme_evan_facet +
    scale_x_continuous(limits = c(-131.2, -54), breaks = c(-60,-80,-100,-120)) + scale_y_continuous(limits = c(25, 61)) + scale_alpha_identity() + 
    theme(legend.position = 'none') + 
    labs(
      x = 'Longitude',
      y = 'Latitude'
    )
  
  ggsave(timeframe_plot, filename = paste0(plot_timeframe,"_n_us_ca", n_events, "_map.pdf"), width = 4, height = 3)
  
  
  }
  timeframe_plots[[plot_num]] <- timeframe_plot # save each timeframe as a new plot within list
  }

plot(timeframe_plots[[1]])

# arrange seasonal plots in a 2x2 grid
flood_seasons <- ggarrange(timeframe_plots[[1]]+theme(axis.title = element_blank()),
                           timeframe_plots[[2]]+theme(axis.title = element_blank()),
                           timeframe_plots[[3]]+theme(axis.title = element_blank()),
                           timeframe_plots[[4]]+theme(axis.title = element_blank()))
drought_seasons <- ggarrange(timeframe_plots[[6]]+theme(axis.title = element_blank()),
                             timeframe_plots[[7]]+theme(axis.title = element_blank()),
                             timeframe_plots[[8]]+theme(axis.title = element_blank()),
                             timeframe_plots[[9]]+theme(axis.title = element_blank()))
                             # labels = c('b','c','d','e'), label.y = 0.23, label.x = 0.18, hjust = 0) 

# annotate seasonal plot grids with axis labels (drought)
drought_seasons <- annotate_figure(drought_seasons, 
                                   bottom = text_grob('Longitude', size = 8), 
                                   left = text_grob('Latitude', size = 8, rot = 90))

# arrange seasonal plots with larger annual plot (drought)
drought_seasons_ann <- ggarrange(plotlist=list(timeframe_plots[[10]],
                                               drought_seasons), 
                                 # align = 'h', 
                                 ncol = 1)
# labels = c('a'), hjust = 0, label.y = 0.29, label.x = 0.14)

# annotate seasonal plot grids with axis labels (flood)
flood_seasons <- annotate_figure(flood_seasons, 
                                   bottom = text_grob('Longitude', size = 8), 
                                   left = text_grob('Latitude', size = 8, rot = 90))

# arrange seasonal plots with larger annual plot (flood)
flood_seasons_ann <- ggarrange(plotlist=list(timeframe_plots[[5]],
                                               flood_seasons), 
                                 # align = 'h', 
                               ncol = 1, heights = c(1, 1))
# labels = c('a'), hjust = 0, label.y = 0.29, label.x = 0.14)

# save drought and flood plots to pdf and png
ggsave(drought_seasons_ann, filename = 'drought_seasons_annual_plot.pdf', useDingbats = F, width = 6.5, height = 8)
ggsave(flood_seasons_ann, filename = 'flood_seasons_annual_plot.pdf', useDingbats = F, width = 6.5, height = 8)
ggsave(drought_seasons_ann, filename = 'drought_seasons_annual_plot.png', width = 6.5, height = 8)
ggsave(flood_seasons_ann, filename = 'flood_seasons_annual_plot.png', width = 6.5, height = 8)

drought_flood_season_combine <- ggarrange(flood_seasons_ann,drought_seasons_ann, ncol = 2)
ggsave(drought_flood_season_combine, filename = 'drought_flood_seasons_annual_plot.pdf',useDingbats = F, width = 6.8, height = 4.5)
ggsave(drought_flood_season_combine, filename = 'drought_flood_seasons_annual_plot.png', width = 6.8, height = 4.5)

#### GENERATE PIE PLOTS ####
n_am <- subset(map_data("world"), region == "USA" | region == "Canada")
unique(n_am$region)

# get cluster centroids
cluster_lat_centroid <- aggregate(cluster_file_2$lat,list(cluster_file_2$Cluster_15), function(x){mean(x)})
cluster_lon_centroid <- aggregate(cluster_file_2$long,list(cluster_file_2$Cluster_15), function(x){mean(x)})

cluster_lon_lat <- data.frame(cbind(cluster_lon_centroid$Group.1,cluster_lon_centroid$x, cluster_lat_centroid$x))
colnames(cluster_lon_lat) <- c("cluster","lon","lat")

# select only applicable cluster centroids
cluster_lon_lat <- cluster_lon_lat[which(cluster_lon_lat$cluster %in% cl_names),]
cluster_file_plot <- cluster_file_2[which(cluster_file_2$Cluster_15 %in% cl_names),]

cluster_map_raw <- ggplot(cluster_lon_lat %>% filter(cluster %in% z_cl_sel), aes(x = lon, y = lat)) + 
  geom_map(data = n_am, map = n_am, aes(map_id=region, x = long, y = lat), color = "grey30", fill = 'white', lwd = 0.5) +
  geom_point(aes(x = long, y = lat, 
                 # fill = as.factor(Cluster_15),
                 color = as.factor(Cluster_15)),
             # color = "grey10", pch = 21, size = 1.25, stroke = 0.25, # remove for smaller size
             # , # for large size
             pch = 16, size = 1.25, # for small size
             data = cluster_file_plot %>% filter(Cluster_15 %in% z_cl_sel)) +
  # geom_point(aes(fill = as.factor(cluster)), color = 'black', pch = 21, size = 4) + # for displaying circle at centroid
  # geom_point(aes(color = as.factor(cluster)), pch = 3, size = 4) + # for displaying cross at centroid, but messes up for some reason
  # geom_polygon(data = fortify(cl_poly_list), aes(x = long, y = lat, color = as.factor(id)), fill = NA) + # add polygon lines
                                                 # fill = as.factor(id), alpha = 0.4)) + # test for adding fill to polygons
  # theme_bw() + scale_fill_manual(values = getPalette(pal_length)) + scale_color_manual(values = getPalette(pal_length)) +
  theme_bw() + scale_fill_manual(values = cl_colors) + scale_color_manual(values = cl_colors) +
  # coordinates need to extend to full extent of polygon data or it gets cut off
  scale_y_continuous(limits = c(23, 60)) + scale_x_continuous(limits = c(-131.2, -54)) + 
  geom_text(aes(label=cluster),hjust=0, vjust=0) +
                    theme_bw() +
    # guides(fill = guide_legend(ncol = 2)) +
    guides(color = guide_legend(ncol = 3)) +
                    theme(
                          # legend.position = "none",
                          legend.position = c(0.86, 0.21),
                          # legend.box = element_rect(size = 0.5),
                          # legend.position = 'top',
                          legend.spacing.x = unit(0, 'in'),
                          # legend.key = element_rect(size = 0.1),
                          legend.key = element_blank(),
                          legend.key.size = unit(1, 'picas'),
                          # legend.justification = c(1, 0),
                          panel.grid = element_blank(),
                          panel.border = element_rect(size = 1),
                          text = element_text(size=9),
                          axis.text = element_text(size = 8),
                          legend.text = element_text(size = 6),
                          legend.title = element_text(hjust = 0.5, size = 7),
                          legend.background = element_blank()
                          # axis.title = element_blank(),
                          # axis.text = element_blank(),
                          # axis.ticks = element_blank()
                          ) + 
        labs(
          x = "Longitude",
          y = "Latitude",
          color = 'Hydro-region'
          # fill = 'Hydro-region'
          # fill = ""
        )
 
cluster_map_raw

cluster_map_raw_HI <- ggplot(cluster_lon_lat %>% filter(cluster %in% z_cl_sel), aes(x = lon, y = lat)) + 
  geom_map(data = n_am, map = n_am, aes(map_id=region, x = long, y = lat), color = "grey30", fill = 'white', lwd = 0.5) +
  geom_point(aes(x = long, y = lat, 
                 # fill = as.factor(Cluster_15)), 
                 color = as.factor(Cluster_15)), 
             # color = "grey10", pch = 21, size = 1.75, stroke = 0.25, # remove for smaller size
                  pch = 16, size = 1.25, # for small size
                        data = cluster_file_plot %>% filter(Cluster_15 %in% z_cl_sel)) +
   theme_bw() + scale_fill_manual(values = cl_colors) + scale_color_manual(values = cl_colors) +
  # coordinates need to extend to full extent of polygon data or it gets cut off
  scale_y_continuous(limits = c(18, 23), breaks = c(19, 22)) + scale_x_continuous(limits = c(-161, -154), breaks = c(-160, -155)) +
  geom_text(aes(label=cluster),hjust=0, vjust=0) +
  theme_bw() +
  # guides(fill = guide_legend(ncol = 2)) +
  # guides(color = guide_legend(ncol = 3)) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.75),
    # text = element_text(size=10),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.background = element_blank(),
    axis.title = element_blank()
    # axis.text = element_blank(),
    # axis.ticks = element_blank()
  ) + 
  labs(
    x = "Longitude",
    y = "Latitude"
    # fill = 'Hydro-region',
    # fill = ""
  )

HI_inset <- annotation_custom(grob = ggplotGrob(cluster_map_raw_HI), # joins plot with centroid coordinates
                  xmin = -135.5, ymin = 20.5, 
                  xmax = -113, ymax = 33)
cluster_map_raw_wHI <- cluster_map_raw + HI_inset
cluster_map_wParallels <- ggarrange(plotlist = list(cluster_map_raw_wHI,parallel_coord_plot), 
                                    heights = c(1,0.8),
                                    ncol = 1)
ggsave(cluster_map_wParallels, 
       filename = 'cluster_map_parallel.pdf', useDingbats = F, width = 3.4, height = 4.5)


ggsave(cluster_map_raw,filename = "cluster_map_raw_polys.pdf",
       useDingbats=FALSE
) 

geom_map(data = cl_poly_list,map = cl_poly_list,aes(region = ID, x = long, y = lat))
  

geom_polygon(data = fortify(cl_poly_list), aes(x = long, y = lat, color = as.factor(group),fill = NA)) + # test out polygon plotting
 scale_color_manual(values = getPalette(15))


# add a column with cluster number to cluster_z
cluster_pie <- cbind(data.frame('cluster_list' = 1:12),cluster_z)
cluster_pie_coords <- cbind(cluster_pie,cluster_lon_lat)
cluster_pie_coords_melt <- melt(cluster_pie_coords,'cluster_list')

cluster_pie_melt <- melt(cluster_pie,'cluster_list')
cluster_pie_melt$sig <- cut(cluster_pie_melt$value,c(-5,-1.96,-1.65,1.65,1.96,5),
                       labels = c("sig95_d","sig90_d","not_sig","sig90_u","sig95_u"))
cluster_pie_drought <- cluster_pie_melt[which(grepl("drought",cluster_pie_melt$variable, fixed = T)==T & 
                                        grepl("annual",cluster_pie_melt$variable, fixed = T)==F),]
cluster_pie_flood <- cluster_pie_melt[which(grepl("flood",cluster_pie_melt$variable, fixed = T)==T & 
                                           grepl("annual",cluster_pie_melt$variable, fixed = T)==F),]

# these plots generate a pie chart for each cluster
# each pie chart is a polar coordinate bar chart, which is a conventional pie chart
cluster_pie_drought_plot <- ggplot(data = cluster_pie_drought, aes(x = variable, fill = as.factor(sig))) + geom_bar(width = 1, color = "black") + 
  coord_polar(start = -pi/4) +
  facet_wrap(.~as.factor(cluster_list)) + 
  scale_fill_manual(values = rev(c("red","pink","white","blue","navy"))) +
  theme_void()

cluster_pie_flood_plot <- ggplot(data = cluster_pie_flood, aes(x = variable, fill = as.factor(sig))) + geom_bar(width = 1, color = "black") + 
  coord_polar(start = -pi/4) +
  facet_wrap(.~as.factor(cluster_list)) + 
  scale_fill_manual(values = c("red","pink","white","blue","navy")) +
  theme_void()

# save the pie charts
ggsave(cluster_pie_drought_plot, filename = "cluster_pie_drought.pdf",
       useDingbats=FALSE
)
ggsave(cluster_pie_flood_plot, filename = "cluster_pie_flood.pdf",
       useDingbats=FALSE
)


# adds consistency to coloring even when not all possibilities are present
pie_colors_f <- c("red","pink","white","blue","navy")
names(pie_colors_f) <- c("sig95_d","sig90_d","not_sig","sig90_u","sig95_u")
pie_colors_d <- rev(c("red","pink","white","blue","navy"))
names(pie_colors_d) <- c("sig95_d","sig90_d","not_sig","sig90_u","sig95_u")


# add pies to map, potentially
# https://stackoverflow.com/questions/43984614/rggplot2geom-points-how-to-swap-points-with-pie-charts/44125392#44125392

# join cluster trend data with cluster coordinates
# for floods
cluster_pie_flood_coords <- merge(cluster_pie_flood, cluster_pie_coords[,c('lat','lon','cluster_list')]) %>%
                            group_by(cluster_list)
cluster_pie_flood_coords[which(grepl("flood",cluster_pie_flood_coords$variable, fixed = T)==T & 
                              grepl("annual",cluster_pie_flood_coords$variable, fixed = T)==F) &
                              !is.na(cluster_pie_flood_coords$value),] %>%
            as.data.frame(cluster_pie_flood_coords) %>% group_by(cluster_list, lat, lon)

cluster_pie_flood_coords$sig <- cut(cluster_pie_flood_coords$value,c(-5,-1.96,-1.65,1.65,1.96,5),
                                    labels = c("sig95_d","sig90_d","not_sig","sig90_u","sig95_u"))
cluster_pie_flood_coords$color_sel <- as.factor(cluster_pie_flood_coords$cluster_list)

#for droughts
cluster_pie_drought_coords <- merge(cluster_pie_drought, cluster_pie_coords[,c('lat','lon','cluster_list')]) %>%
                            group_by(cluster_list)
cluster_pie_drought_coords[which(grepl("drought",cluster_pie_drought_coords$variable, fixed = T)==T & 
                              grepl("annual",cluster_pie_drought_coords$variable, fixed = T)==F) &
                              !is.na(cluster_pie_drought_coords$value),]
  cluster_pie_drought_coords <- as.data.frame(cluster_pie_drought_coords) %>% group_by(cluster_list, lat, lon)

cluster_pie_drought_coords$sig <- cut(cluster_pie_drought_coords$value,c(-5,-1.96,-1.65,1.65,1.96,5),
                                    labels = c("sig95_d","sig90_d","not_sig","sig90_u","sig95_u"))
cluster_pie_drought_coords$color_sel <- as.factor(cluster_pie_drought_coords$cluster_list)
# labels = c("red","pink","white","blue","navy"))

# creates a list of pie plots; each element of list corresponds to a cluster
# each plot is joined with the coordinates of the cluster centroid

pie_size <- 4
flood_pie.grobs <- cluster_pie_flood_coords %>% # this symbol passes the item to the left to the function on right
  do(subplots = 
       # ggplot(., aes(x = variable, fill = as.factor(sig))) + geom_bar(width = 1, color = "black") +
       ggplot(., aes(x = variable, fill = as.factor(sig), color = as.factor(color_sel))) + 
       geom_bar(width = 1, 
                # color = "black", size = 1) +
                # color = .$color_sel,
                size = 0.75) +
       coord_polar(start = -pi/4) +
       scale_fill_manual(values = pie_colors_f) +
       # scale_color_manual(values = cl_colors)) %>% 
       scale_color_manual(values = cl_colors) +
       theme_void() + guides(fill = F) + theme(legend.position = "none")) %>%
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), # joins plot with centroid coordinates
                                           x = lon-pie_size, y = lat-pie_size, 
                                           xmax = lon+pie_size, ymax = lat+pie_size))) # sets resulting size of plot

drought_pie.grobs <- cluster_pie_drought_coords %>% # this symbol passes the item to the left to the function on right
  do(subplots = 
       # ggplot(., aes(x = variable, fill = as.factor(sig))) + geom_bar(width = 1, color = "black") +
       ggplot(., aes(x = variable, fill = as.factor(sig), color = as.factor(color_sel))) + 
       geom_bar(width = 1, 
                # color = "black", size = 1) +
                # color = .$color_sel,
                size = 0.75) +
       coord_polar(start = -pi/4) +
       scale_fill_manual(values = pie_colors_d) +
       # scale_color_manual(values = cl_colors)) %>% 
       scale_color_manual(values = cl_colors) +
       theme_void() + guides(fill = F) + theme(legend.position = "none")) %>%
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), # joins plot with centroid coordinates
                                           x = lon-pie_size, y = lat-pie_size, 
                                           xmax = lon+pie_size, ymax = lat+pie_size))) # sets resulting size of plot
# test out the individual pie plots
drought_pie.grobs$subplots[1]
# flood_pie.grobs$subplots[2]
# flood_pie.grobs$subplots[3]

# adds pies to cluster map
cl_pie_fl_all <- flood_pie.grobs%>% 
{cluster_map_raw + .$subgrobs}

ggsave(cl_pie_fl_all, filename = "cl_pie_fl_all_sm.pdf",
       useDingbats=FALSE
)
cl_pie_dr_all <- drought_pie.grobs%>% 
{cluster_map_raw + .$subgrobs}

ggsave(cl_pie_dr_all, filename = "cl_pie_dr_all_sm.pdf",
       useDingbats=FALSE
)

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

#### FUNCTIONS ####
getKernelDensity <- function(station_df, smooths, timeframe, individual){
  
  # number of events to find (like choosing a percentile, but more stable)
  rank_no <- 20
  
  if(individual == "individual"){
    drought_years <- unique(station_df$water_year)
    drought_dates <- data.frame(matrix(nrow = length(drought_years)))
    drought_dates <- rep(NA,length(drought_years))
    for(i in 1:length(drought_years)){
      drought_year <- drought_years[i]
      start_row <- min(which(is.na(drought_dates)))
      drought_yr_sub <- station_df[which(station_df$water_year == drought_year),]
      drought_row <- drought_yr_sub[which.min(drought_yr_sub$q_cfs),]
      drought_date <- drought_row$datetime
      drought_q_cfs <- drought_row$q_cfs
      # print(c(drought_row$water_year,drought_q_cfs))
      drought_dates[i] <- drought_date
    }
    
    station_df[which(!(station_df$datetime %in% drought_dates)),]$q_cfs_droughts <- 1e8
    drought_rank_val <- sort(station_df$q_cfs_droughts, decreasing = F)[rank_no]
    # rank_no <- 20
    # make a vector with years of flood above threshold value based on rank
    peak_rank_val <- sort(station_df$q_cfs_peaks,decreasing = T)[rank_no]
    # make a vector with years of drought above threshold value based on rank
    # drought_rank_val <- sort(station_df$q_cfs_droughts,decreasing = F)[rank_no]
    peak_rank <- station_df[which(station_df$q_cfs_peaks >= peak_rank_val),]
    drought_rank <- station_df[which(station_df$q_cfs_droughts <= drought_rank_val),]
    
    t = c(station_df$water_year[which(station_df$q_cfs_peaks >= peak_rank_val)])
    t_d = c(station_df$water_year[which(station_df$q_cfs_droughts <= drought_rank_val)])
     }else{
    t <- station_df[which(station_df$flood>0),c(1,2)]
    t_d <- station_df[which(station_df$drought>0),c(1,3)]
  } 
  # print(c(t, timeframe))
  # get year bounds of record - beyond bounds of flood dates
  t_obs_l <- min(station_df$water_year) # first year
  t_obs_r <- max(station_df$water_year) # last year
  
  bounds <- c(t_obs_l,t_obs_r)
  # smooths <- c(1,2,5,10,15,20,25,35)
  smooths <- c(1,2,5,10,15,20) # smoothing windows to test
  smoothTest <- getSeries(t, bounds, smooths,T,F) # get occurrence rates for floods for different windows
  smoothTest_d <- getSeries(t_d, bounds,smooths,T,F) # get occurrence rates for droughts for different windows
  # returns a list with two items:
  # 1) the kernel density smoother $recur_df_subset
  # 2) the vector of dates with pseudo-dates appended to each side $t_add
  
  t_final <- smoothTest$recur_df_subset # selects the kernel density smoothed data (entry for every year)
  t_final_d <- smoothTest_d$recur_df_subset # selects the kernel density smoothed data (entry for every year)
  # colSums(t_final[2:ncol(t_final)]) # sum the columns - each column is the smooth timeseries for a given window
  
  # t_final[(round(t_final$year) %in% t),2:ncol(t_final)] # converts years from decimal value back to matching years from t
  
  h_opt <- getCV(t,bounds,t_final,smooths)
  h_opt_d <- getCV(t_d,bounds,t_final_d,smooths)
  
  h_opt <- 15
  h_opt_d <- 15
  seriesFinal <- getSeries(t,bounds,h_opt,T,F)
  t_opt <- seriesFinal$recur_df_subset
  t_add <- data.frame(cbind(seriesFinal$t_add,seriesFinal$weights))
  colnames(t_add) <- c('t_add','weights')
  
  seriesFinal_d <- getSeries(t_d,bounds,h_opt_d,T,F)
  t_opt_d <- seriesFinal_d$recur_df_subset
  t_add_d <- data.frame(cbind(seriesFinal_d$t_add,seriesFinal_d$weights))
  colnames(t_add_d) <- c('t_add','weights')
  # plot(smooths,t_final_sq)
  # plot(smooths,h_cv_sum)
  # plot(t_add)
  
  t_opt_combine <- getBoot(t_opt,bounds,t_add,h_opt)
  t_opt_combine_d <- getBoot(t_opt_d,bounds,t_add_d,h_opt_d)
  
  t_opt_flood_drought <- cbind(t_opt_combine,t_opt_combine_d)
  colnames(t_opt_flood_drought) <- c("year","t_opt","Ely","t_opt_l","t_opt_u","year_d","t_opt_d","Ely_d","t_opt_l_d","t_opt_u_d")
  
  if(individual == "individual"){
  z <- (sum(t)/length(t) - (t_obs_r + t_obs_l)/2)/((t_obs_r - t_obs_l)/sqrt(12*length(t)))
  z_d <- (sum(t_d)/length(t_d) - (t_obs_r + t_obs_l)/2)/((t_obs_r - t_obs_l)/sqrt(12*length(t_d)))
  }else{
    n_stat_t <- nrow(t)
    n_stat_td <- nrow(t_d)
    weight_sum <- t[,2]/(sum(t[,2])/n_stat_t)
    t_centroid <- sum(t[,1]*weight_sum)/n_stat_t
    weight_sum_d <- t_d[,2]/(sum(t_d[,2])/n_stat_td)
    t_centroid_d <- sum(t_d[,1]*weight_sum_d)/n_stat_td
    z <- (t_centroid - (t_obs_r + t_obs_l)/2)/((t_obs_r - t_obs_l)/sqrt(12*n_stat_t))
    z_d <- (t_centroid_d - (t_obs_r + t_obs_l)/2)/((t_obs_r - t_obs_l)/sqrt(12*n_stat_td))
  }
  #### PLOT SETUP ####
  
  y_label = parse(text=paste("'occurrence [yr'", "^-1 ", "*']'", sep=""))
  
  t_opt_flood_drought <- cbind(t_opt_combine,t_opt_combine_d)
  colnames(t_opt_flood_drought) <- c("year","t_opt","Ely","t_opt_l","t_opt_u","year_d","t_opt_d","Ely_d","t_opt_l_d","t_opt_u_d")
  y_label = parse(text=paste("'Occurrence [yr'", "^-1 ", "*']'", sep=""))
  cf_plot <- ggplot(data = t_opt_flood_drought) + theme_bw() + theme(
         # panel.border = element_rect(colour = "black", size = 1, fill = NA),
         axis.line.x.top = element_line(color = "black", size = 0.75),
         axis.line.y.left = element_line(color = "black", size = 0.75),
         axis.line.y.right = element_line(color = "black", size = 0.75),
         axis.ticks = element_line(color = "black", size = 0.75),
         axis.ticks.length = unit(4,"points"),
         panel.border = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank()) + 
          scale_x_continuous(limits = c(t_obs_l,t_obs_r), position = 'top', expand = c(0,0))
          # scale_y_continuous(trans = "reverse")
  cf_plot <- cf_plot + geom_ribbon(aes(x = year,ymin = t_opt_l, ymax = t_opt_u), fill = "navy", alpha = 0.25)
  cf_plot <- cf_plot + geom_line(aes(x = year, y = Ely), color = 'navy', lwd = 1)
  cf_plot <- cf_plot + geom_ribbon(aes(x = year_d,ymin = t_opt_l_d, ymax = t_opt_u_d), fill = "red", alpha = 0.25)
  cf_plot <- cf_plot + geom_line(aes(x = year_d, y = Ely_d), color = 'red', lwd = 1) +
    # theme(axis.title.x.top = element_text('Year')) +
    labs(
      x = paste("Year", " [", timeframe, "]", sep = ""),
      y = y_label
      # shape = "Drainage Area"
  ) 

  # plot(cf_plot)
  
  if(individual == "individual"){
  cf_plot_events <- ggplot(data = data.frame(t)) + geom_linerange(data = data.frame(t), 
                  aes(x = t, ymin = rep(0,length(t)), ymax = rep(1,length(t))), 
                  color = 'navy', lwd = 0.75)
  cf_plot_events <- cf_plot_events + geom_linerange(data = data.frame(t_d), 
                   aes(x = t_d+0.5, ymin = rep(0,length(t_d)), ymax = rep(1,length(t_d))), 
                   color = 'red', lwd = 0.75)
  cf_plot_events <- cf_plot_events + theme_void() +
    scale_x_continuous(limits = c(t_obs_l,t_obs_r),expand = c(0,0))
    
  
  # plot(cf_plot_events)
  
  
  
  y_label_q = parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep=""))
  station_df$frac_date <- station_df$year + station_df$doy/365
  cf_plot_ts <- ggplot(data = station_df) + geom_line(aes(frac_date,q_cfs),lwd = 0.25, color = 'navy') +
                              geom_hline(yintercept = peak_rank_val, color = 'grey50', linetype = "dashed") + 
                              geom_hline(yintercept = drought_rank_val, color = 'grey50', linetype = "dashed") + 
                              theme_bw() + 
                              theme(
                             # panel.border = element_rect(colour = "black", size = 1, fill = NA),
                             axis.line.x.bottom = element_line(color = "black", size = 0.75),
                             axis.line.y.left = element_line(color = "black", size = 0.75),
                             axis.ticks = element_line(color = "black", size = 0.75),
                             axis.ticks.length = unit(4,"points"),
                             panel.border = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank()) + 
    scale_x_continuous(limits = c(t_obs_l,t_obs_r),expand = c(0,0)) +
    scale_y_log10(labels = fancy_scientific)+ 
    labs(
      x = paste("Year", " [", timeframe, "]", sep = ""),
      y = y_label_q
      # shape = "Drainage Area"
    ) 
  # plot(cf_plot_ts)
  
  # ggarrange details http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/  
  plot_export <- ggarrange(cf_plot,
              # theme(axis.title.x=element_blank(),axis.text.x.bottom = element_blank()), 
            cf_plot_events, cf_plot_ts, heights = c(2, 0.1, 2),
            ncol = 1, nrow = 3, align = "v", labels = c("A","","B"))
  # plot(plot_export)
  # ggexport(plot_export, filename = paste(toString(station_df$site_no[1]), "_", timeframe, "_drought_flood_recur_ts.pdf",sep = ""))
  }else{
    plot_export <- cf_plot
  }
  
  #### NEED TO MAKE THIS AN IF STATEMENT BASED ON INDIVIDUAL PASS #####
  return(list(plot_export,z,z_d))
  # return(c(z,z_d,plot_export))
}

getT <- function(station_df, timeframe, rank_no){
  # # make a vector with years of flood above threshold percentile
  # t = c(station_df$water_year[which(station_df$q_cfs_peaks>quantile(station_df$q_cfs_peaks,peak_thresh))])
  # # make a vector with years of drought below threshold percentile
  # t_d = c(station_df$water_year[which(station_df$q_cfs_droughts<quantile(station_df$q_cfs_droughts,drought_thresh))])
  
  # select only lowest discharge of each season or year
  
  drought_years <- unique(station_df$water_year)
  drought_dates <- data.frame(matrix(nrow = length(drought_years)))
  drought_dates <- rep(NA,length(drought_years))
  for(i in 1:length(drought_years)){
    drought_year <- drought_years[i]
    # start_row <- min(which(is.na(drought_dates)))
    drought_yr_sub <- station_df[which(station_df$water_year == drought_year),]
    drought_row <- drought_yr_sub[which.min(drought_yr_sub$q_cfs),]
    drought_date <- drought_row$datetime
    drought_q_cfs <- drought_row$q_cfs
    # print(c(drought_row$water_year,drought_q_cfs))
    drought_dates[i] <- drought_date
  }
  
  station_df[which(!(station_df$datetime %in% drought_dates)),]$q_cfs_droughts <- 1e8
  drought_rank_val <- sort(station_df$q_cfs_droughts, decreasing = F)[rank_no]
  # rank_no <- 20
  # make a vector with years of flood above threshold value based on rank
  peak_rank_val <- sort(station_df$q_cfs_peaks,decreasing = T)[rank_no]
  # make a vector with years of drought above threshold value based on rank
  # drought_rank_val <- sort(station_df$q_cfs_droughts,decreasing = F)[rank_no]
  peak_rank <- station_df[which(station_df$q_cfs_peaks >= peak_rank_val),]
  drought_rank <- station_df[which(station_df$q_cfs_droughts <= drought_rank_val),]
  
  t = c(station_df$water_year[which(station_df$q_cfs_peaks >= peak_rank_val)])
  t_d = c(station_df$water_year[which(station_df$q_cfs_droughts <= drought_rank_val)])
  
  # get year bounds of record - beyond bounds of flood dates
  # t_obs_l <- min(station_df$water_year) # first year, but uses full record
  t_obs_l <- 1960 # first year, but just going back to 1960
  t_obs_r <- max(station_df$water_year) # last year
  
  # adding zeros between years with positive values, getting count and cumulative count per year
  t_all <- data.frame(t_obs_l:t_obs_r)
  t_all$drought_count <- matrix(nrow = length(t_obs_l:t_obs_r))
  t_all$flood_count <- matrix(nrow = length(t_obs_l:t_obs_r))
  colnames(t_all) <- c("water_year", "drought_count","flood_count")
  for(i in 1:length(t_all$water_year)){
    # t_all$cumu[i] <- length(which(t<=t_all$water_year[i]))
    t_all$drought_count[i] <- length(which(t_d==t_all$water_year[i]))
    t_all$flood_count[i] <- length(which(t==t_all$water_year[i]))
  }
  
  return(t_all)}
  
getSeries <- function(t, bounds, smooths, plot_binary, pseudodata){
    # repress plot function for batch processing
  plot_binary <- F
  
  t_obs_l <- bounds[1]
  t_obs_r <- bounds[2]
  
  
  
  if(is.null(ncol(t)) == F){
    t_weights <- t[,2]
    t <- t[,1]
    t_add <- c(rev((min(t)-(t-min(t)))), t, rev(max(t)+(max(t)-t))) # reflection pseudodata - the best for this application 
    t_weights_add <- c(rev(t_weights), t_weights, rev(t_weights)) # reflection pseudodata - the best for this application 
    }else if(pseudodata == F){
  # Add pseudo-data to beginning of flood/drought years...
     # but do not add pseudodata if it's already a part of t
      t_add <- c(rev((min(t)-(t-min(t)))), t, rev(max(t)+(max(t)-t))) # reflection pseudodata - the best for this application
      t_weights_add <- c(rep(1, length(t_add)))
    }else{
      t_add <- t
      t_weights_add <- c(rep(1, length(t_add)))}
  
  
  
  
    y_label = parse(text=paste("'flood occurrence [yr'", "^-1 ", "*']'", sep=""))
    
    smooths <- smooths
    recur_df <- data.frame(matrix(ncol = length(smooths)+1, nrow = 512))
    colnames(recur_df) <- (c("year",smooths))
    
    series_length <- length(t_obs_l:t_obs_r)
    hrelmax <- 0.5
    
    design_points <- round(t_obs_l - 3 * series_length * hrelmax):round(t_obs_r + 3 * series_length * hrelmax)
    
    
    
    for(j in 1:length(smooths)){
      bw_select = smooths[j]
      flood_density <- density(x=t_add,bw=bw_select,kernel="gaussian", from = min(design_points), to = max(design_points), weights = t_weights_add)
      # could instead use built-in bandwidth estimator 'ucv' or 'bcv' (unbiased or biased cross-validation)
      # flood_density <- density(x=t_add,bw='ucv',kernel="gaussian", from = min(design_points), to = max(design_points), weights = c(rep(1, length(t_add))))
      recur_df[,(j+1)] <- flood_density[["y"]]
      if(j==1){recur_df$year <- flood_density[["x"]]}
    }
    
    recur_df_subset <- subset(recur_df,recur_df$year>=t_obs_l & recur_df$year<=t_obs_r)
    
    if(plot_binary==T){recur_plot_data <- melt(recur_df_subset,"year")
    recur_plot <- ggplot(data = recur_plot_data,aes(x = year, y = value, group = variable, color = variable)) + theme_bw() 
    recur_plot <- recur_plot + geom_line() + scale_color_viridis_d() + theme(panel.border = element_rect(colour = "black", size = 1, fill = NA),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
    recur_plot <- recur_plot + xlim(t_obs_l, t_obs_r) + labs(
      x = "Year",
      y = y_label,
      colour = "smoothing [yrs]"
      # shape = "Drainage Area"
    )# par(new = T)
    print(recur_plot)}
    
    return(list("recur_df_subset" = recur_df_subset, "t_add" = t_add, "weights" = t_weights_add))
}

getCV <- function(t,bounds,t_final,smooths){ 
  # CROSS VALIDATION
  # this function does cross validation by comparing the area under the curve of the smoothed series
  # to the average area under the curve that would result from smoothing the series with one original t values removed
  # it removes each t value successively, then computes the remaining integral
  # the result checks the sensitivity to removing a value. when it doesn't make much difference to remove a value
  # that is the optimal smoothing window
  if(is.null(ncol(t)) == F){
    t <- t[,1]
      }
  
  h_cv <- data.frame(matrix(ncol=length(smooths),nrow=length(t))) # initialize a cross validation data frame
  colnames(h_cv) <- smooths # adds smoothing window lengths as column names
  
  for(i in 1:length(t)){
    t_mod <- t[-i] # removes one element of timeseries
    t_cv <- getSeries(t_mod,bounds,smooths,F,F)$recur_df_subset # computes recurrence at each year for all the smoothing windows in "smooths" var
    h_cv[i,] <- colSums(t_cv[2:ncol(t_cv)]) # writes sum of "take one out" sum 
  }
  
  
  h_cv_sum <- 2*colSums(h_cv)/length(t) # takes average sum for each smoothing window/2
  t_final_sq <- colSums(t_final^2)[2:(length(smooths)+1)] # integrates the original smoothed data
  h_cv_diff <- t_final_sq-h_cv_sum # difference of original smooth data with the average of the take-one-out data
  # REMOVE PLOT for batch
  # plot(x = smooths,y = h_cv_diff, type = "l")
  
  h_opt <- min(20,smooths[which.min(h_cv_diff)]) # finds optimal smoothing window by x-validation - minimizes value
  h_opt <- max(5,h_opt) # finds optimal smoothing window by x-validation - keeps from no smoothing
  # h_opt <- min(20,smooths[(which.min(h_cv_diff)+1)]) # shifts optimal smoothing window to interval smoother
  return(h_opt)
}

getBoot <- function(t_opt,bounds,t_add,h_opt){
  boot_n = 2000
  
  # year_range <- round(min(t_opt$year):max(t_opt$year))
  # boot_ly <- data.frame(matrix(nrow=length(year_range),ncol=boot_n))
  
  boot_ly <- data.frame(matrix(nrow=nrow(t_opt),ncol=boot_n))
  # bootstrap - resample then repeat process of finding recurrence; use that to get confidence bands
  
  # NEED TO MODIFY RESAMPLE TO ACCOMODATE TADD BEING A DATA FRAME, SO IT JUST SORTS BASED ON
  # THE YEAR COLUMN, BUT KEEPS THE WEIGHTS ASSOCIATED WITH THAT ROW OF TADD
  
  for(i in 1:boot_n){
    resample <- t_add[sort(as.integer(resample_bootstrap(data.frame(t_add$t_add)))),]
    boot_iter <- getSeries(resample, bounds,h_opt, F, T)$recur_df_subset
    # print(c((max(resample)-min(resample)),length(unique(resample)),nrow(boot_iter)))
    # print(nrow(boot_iter))
    boot_ly[1:nrow(boot_iter),i] <- boot_iter[,2]
    # boot_ly[which(year_range %in% round(boot_iter$year)),i] <- boot_iter[,2]
  }
  
  boot_ly[which(boot_ly==0,arr.ind = T)] <- NA
  # max_index <- min(which(is.na(boot_ly[1])))-1
  # boot_ly <- boot_ly[1:max_index,]
  Ely <- rowMeans(boot_ly, na.rm=T) # row average for all bootstrap samples
  Tly <- (boot_ly - Ely)/sqrt(boot_ly) # calculate T-score for each element in boot_ly
  # Tly[!is.finite(Tly)] <- NA
  Ta <- rowQuantiles(data.matrix(abs(Tly)),probs = c(0.9,0.95),na.rm=T)
  # t_opt$year <- round(t_opt$year)
  # t_opt <- t_opt[!duplicated(t_opt$year),]
  t_opt_l <- Ely-Ta[,1]*sqrt(t_opt[,2])
  t_opt_l[t_opt_l<0] <- 0
  t_opt_u <- Ely+Ta[,1]*sqrt(t_opt[,2])
  
  colnames(t_opt) <- c('year','t_opt')
  
  # t_opt_combine <- data.frame(cbind(year_range,Ely,t_opt_l,t_opt_u))
  # colnames(t_opt_combine) <- c('year','Ely','t_opt_l','t_opt_u')
  t_opt_combine <- cbind(t_opt,Ely,t_opt_l,t_opt_u) # for some reason discrepancy btw. t_opt and boot_ly years
  return(t_opt_combine)
}

getFileList <- function(dt){
  dt
}
createFileList_cluster <- function(df_region){
  if(any(which(nchar(as.character(df_region$site_number))==7))){
   file_list_ca <- paste0(df_region$site_number[which(nchar(as.character(df_region$site_number))==7)],"_Daily_Flow_ts.csv")
   }else{
    file_list_ca <- c()   
   }
  if(any(which(nchar(as.character(df_region$site_number))==8))){
   file_list_us <- paste0(df_region$site_number[which(nchar(as.character(df_region$site_number))==8)],".dat")
   }else{
     file_list_us <- c()
   }
  return(c(file_list_ca,file_list_us))
}


getDates <- function(station_df){
  #put date columns in data frame
  if(colnames(station_df)[2]=="agency_cd"){
    dates <- as.POSIXlt(as.character(station_df$Date),format = '%Y-%m-%d')
    
    station_df$q_cfs[station_df$q_cfs < 0] <- NA
    station_df$datetime <- as.POSIXct(as.character(station_df$Date),format = '%Y-%m-%d')
    
    station_df$year <- dates$year # changes year column to date format
    station_df$year <- station_df$year + 1900
    station_df$doy <- dates$yday 
    station_df$month <- dates$mon # makes new column with month for each date (stored as 0 - 11)
    station_df$month_nm <- strftime(dates,'%b')
    station_df <- station_df[c("site_no","Date","water_year","q_cfs","datetime","year","doy","month","month_nm")]
  }else{
    station_df <- station_df[1:(nrow(station_df)-2),]
    dates <- as.POSIXlt(as.character(station_df$Date),format = '%Y/%m/%d')
    
    station_df$Date <- dates
    station_df$site_no <- station_df$ID
    station_df <- addWaterYear(station_df)
    station_df$q_cfs <- station_df$Flow / 0.0283
    station_df$q_cfs[station_df$Flow < 0] <- NA
    station_df$datetime <- station_df$Date
    
    station_df$year <- dates$year # changes year column to date format
    station_df$year <- station_df$year + 1900
    station_df$doy <- dates$yday 
    station_df$month <- dates$mon # makes new column with month for each date (stored as 0 - 11)
    station_df$month_nm <- strftime(dates,'%b') 
    
    station_df <- station_df[c("site_no","Date","waterYear","q_cfs","datetime","year","doy","month","month_nm")]
    colnames(station_df) <- c("site_no","Date","water_year","q_cfs","datetime","year","doy","month","month_nm")
  }
  
  return(station_df)
}

getNormDischarge <- function(station_df, station_info){
  #pull USGS info for that site
  
  site_no <- toString(station_df$site_no[1]) 
  if (nchar(site_no) == 7){
    site_no <- paste("0", site_no, sep = "")
  } else { 
    site_no <- site_no
  }
  rownum <- which(site_no == station_info$site_no)
  
  #add drainage area to station_df
  station_df$drainage_area <- station_info$DRAIN_SQKM[rownum] 
  
  #add normalized discharge column
  station_df$q_normalized <- 100 * station_df$q_cfs / station_df$drainage_area
  
  #add elevation (****CHOOSE ELEV PARAMETER***) and corresponding color
  station_df$elev <- station_info$max[rownum]
  station_df$elev_scale <- station_info$max_scale[rownum]
  station_df$elev_color <- station_info$color[rownum]
  
  return(station_df)
}

getPeaks <- function(station_df){
  end_row <- nrow(station_df)
  station_df$q_cfs_random <- station_df$q_cfs + runif(n = nrow(station_df), min = 0.001, max = 0.01)
  station_df$q_cfs_peaks <- station_df$q_cfs #initialize: assume every day is peak, then replace non-peaks
  station_df$q_cfs_droughts <- station_df$q_cfs #initialize: assume every day is drought, then replace non-droughts
  for(select_row in 4:(end_row-3)){
    row_range <- c((select_row-3):(select_row+3)) # 7-day window of rows
    if(station_df$q_cfs_random[select_row]!=max(station_df$q_cfs_random[row_range],na.rm=T)| # is this day the peak day?
       is.na(station_df$q_cfs[select_row])){ # or is it NA?
      station_df$q_cfs_peaks[select_row] <- 1}
  }
  # for(select_row in 16:(end_row-15)){
  #   row_range <- c((select_row-15):(select_row+15)) # 30-day window of rows
  #   if(station_df$q_cfs_random[select_row]!=min(station_df$q_cfs_random[row_range],na.rm=T)| # is this day the peak day?
  #      is.na(station_df$q_cfs[select_row])){ # or is it NA?
  #     station_df$q_cfs_droughts[select_row] <- 1e8} # if it is not peak day, set value = 1
  #   # print(select_row) # error check row
  # }
  
  
  return(station_df)
}



getAreaPlot <- function(t_all, type_fd, bounds, timeframe){
  t_obs_l <- bounds[1]
  t_obs_r <- bounds[2]
  if(type_fd=="flood"){
    palette_select = "Blues"
    yaxis_title = "Floods"
  }else{
    palette_select = "Reds"
    yaxis_title = "Droughts"
  }
  
  
  
  t_all_melt_plot <- melt(t_all,"water_year", value.name = "flood_count")
  area_plot <- ggplot(data = t_all_melt_plot, aes(water_year, variable)) + geom_raster(aes(fill = as.factor(flood_count))) +
  scale_fill_brewer(palette = palette_select)+
  # guides(fill = guide_legend(title = "n events"))+
  guides(fill = F)+
  theme_bw() + 
  theme(
    # panel.border = element_rect(colour = "black", size = 1, fill = NA),
    axis.line.x.bottom = element_blank(),
    axis.line.y.left = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.75),
    axis.ticks.length = unit(4,"points"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(color = "black", size = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) + 
  scale_x_continuous(limits = c(t_obs_l,t_obs_r), expand=c(0,0)) + # expand eliminates buffer zone around plot
  scale_y_discrete(expand=c(0,0))+ 
    labs(
      x = paste("Year", " [", timeframe, "]", sep = ""),
      y = yaxis_title
      # shape = "Drainage Area"
    ) 
return(area_plot)
}

# converts log10 axis values to format 10^x
fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- log10(l)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(parse(text = paste("10^",as.character(l),sep = "")))
} 

# LEFTOVERS - PROBABLY DELETE

drought_rows <- which(station_df$q_cfs_droughts!=1e8)
drought_dates <- station_df$datetime[drought_rows]
drought_q_cfs <- station_df$q_cfs[drought_rows]
drought_plot <- data.frame(cbind(drought_dates,drought_q_cfs))
ggplot(data = drought_plot,aes(x=drought_dates,y=drought_q_cfs)) + geom_smooth()

# Example: Rio Negro at Manaus - example from 2018 Science Advances Paper
# Flood years (t4)
t = c(1922,1953,1976,1989,1999,2009,2012,2013,2014,2015) # years of flood events on Rio Negro, Manaus
t = c(1908,1922,1953,1974,1976,1989,1994,1999,2009,2012,2013,2014,2015) # years of flood events on Rio Negro, Manaus; from 2018 paper
t = c(1928, 1936,1938,1973,2002,2011)
t = c(1922,1925,1927,1936,1938,1947,1948,1952,1973,1976,1987,1998,2000,2002,2011)
peak_thresh = 0.9995
t_obs_l <-  1903 # first year of observations
t_obs_r <- 2016 # last year of observations
t_obs_l <-  1916 # first year of observations, white r


t_dens <- density(x = t_refl_d, bw = 'bcv', adjust = 0.75, from = 1915, to = 2019, kernel = 'gaussian', weights = c(rep(1, length(t_refl_d))), main = NA)
t_dens <- density(x = t_refl_d, bw = 15, from = 1915, to = 2019, kernel = 'gaussian', weights = c(rep(1, length(t_refl_d))))

t_refl <- c(rev(min(t)-(t-min(t))), t, rev(max(t)+(max(t)-t))) # reflection pseudodata - the best for this application
plot(t_refl)

plot(kdensity(x=(t-min(t))/(max(t)-min(t)),  bw = 1/102, kernel = "gamma",normalized = T))
plot(kdensity(x=t,bw = 15/512, kernel = "gamma",normalized = F))
plot(density(x = t, bw = 5, from = 1915, to = 2017, weights = c(rep(1, length(t)))))
plot(density(x = t_refl, bw = 15, from = 1915, to = 2019, kernel = 'gaussian', weights = c(rep(1, length(t_refl)))))

getSeriesPseudo <- function(t, smooths, plot_binary){
  # Add pseudo-data to beginning of flood/drought years...
  t_p = rep(NA,length(t)) # initialize pseudo-data dataframe
  t_pr = rep(NA,length(t)) # initialize pseudo-data dataframe
  
  for(i in 1:(length(t))){
    if(i<3){floor_t <- t_obs_l}else{floor_t <- t[floor(i/3)]}
    if(i<2){floor_t2 <- t_obs_l}else{floor_t2 <- t[floor((2*i)/3)]}
    if(i<3){ceiling_t <- t[1]}else{ceiling_t <- t[ceiling(i/3)]}
    if(i<2){ceiling_t2 <- t[1]}else{ceiling_t2 <- t[ceiling((2*i)/3)]}
    interp1 <- i%%3/3
    interp2 <- (2*i)%%3/3
    t_i_3_interp <- round((ceiling_t - floor_t)*interp1+floor_t)
    t_i2_3_interp <- round((ceiling_t2 - floor_t2)*interp2+floor_t2)
    
    
    t_i_3 <- 5*round(t_i_3_interp-t_obs_l)
    t_i2_3 <- 4*round(t_i2_3_interp-t_obs_l)
    t_i = round((10/3)*(t[i]-t_obs_l))
    
    t_p[i] = t_obs_l - t_i_3 - t_i2_3 + t_i
    
    # print(c(i,t[i],floor_t, ceiling_t, t_i_3_interp, t_i2_3_interp, t_i_3,t_i2_3,t_i,t_p[i]))
  }
  t_inv = sort(t, decreasing = T)
  for(i in 1:(length(t_inv))){
    if(i<3){floor_t <- t_obs_r}else{floor_t <- t_inv[floor(i/3)]}
    if(i<2){floor_t2 <- t_obs_r}else{floor_t2 <- t_inv[floor((2*i)/3)]}
    if(i<3){ceiling_t <- t_inv[1]}else{ceiling_t <- t_inv[ceiling(i/3)]}
    if(i<2){ceiling_t2 <- t_inv[1]}else{ceiling_t2 <- t_inv[ceiling((2*i)/3)]}
    interp1 <- i%%3/3
    interp2 <- (2*i)%%3/3
    t_i_3_interp <- round((ceiling_t - floor_t)*interp1+floor_t)
    t_i2_3_interp <- round((ceiling_t2 - floor_t2)*interp2+floor_t2)
    
    t_i_3 <- 5*round(t_i_3_interp-t_obs_r)
    t_i2_3 <- 4*round(t_i2_3_interp-t_obs_r)
    t_i = round((10/3)*(t_inv[i]-t_obs_r))
    
    t_pr[i] = t_obs_r - t_i_3 - t_i2_3 + t_i
    
    # print(c(i,t_inv[i],floor_t, ceiling_t, t_i_3_interp, t_i2_3_interp, t_i_3,t_i2_3,t_i,t_pr[i]))
  }
  
  t_extended <- c(t_obs_l,t,t_obs_r)
  t_space <- rep(1,3*length(t_extended))
  
  for(i in 1:length(t_space)){
    if((i+2)%%3==0){
      t_space[i] <- t_extended[(i+2)/3]
    }else if((i+1)%%3==0){
      t_space[i] <- t_extended[(i+1)/3]*(2/3)+t_extended[(i+4)/3]*(1/3)
    }else{
      t_space[i] <- t_extended[(i)/3]*(1/3)+t_extended[(i+3)/3]*(2/3)}
  }
  
  # t(i) = tobs−l − 5 [t4(i/3) − tobs−l] − 4 [t4(2i/3) − tobs−l] + 10/3 [t4(i) − tobs−l] 
  
  
  g <- rep(NA,length(t_space))
  for(index in 1:40){
    i <- (index+3)*3
    i2 <- (index+3)*2
    i3 <- index+3
    yi <- t_space[i]
    y2 <- t_space[i2]
    y3 <- t_space[i3]
    # ti <- (10/3)*(yi-t_obs_l)
    # ti2 <- -4*(y2-t_obs_l)
    # ti3 <- -5*(y3-t_obs_l)
    # g[i3] <- t_obs_l + ti3 + ti2 + ti
    ti <- (10/3)*(yi-1929)
    ti2 <- -4*(y2-1929)
    ti3 <- -5*(y3-1929)
    g[i3] <- 1929 + ti3 + ti2 + ti
    print(round(c(yi,y2,y3,g[i3])))
  }
  plot(rev(g))
  plot(c(rev(g),t))
  
  # find if/where pseudo-dates are greater than first observation
  end_t_p <- if(!is_empty(which(t_p>t_obs_l))){end_t_p <- min(which(t_p>t_obs_l))-1}else{end_t_p <- length(t_p)}
  start_t_p <- if(!is_empty(which(t_pr<t_obs_r))){start_t_p <- min(which(t_pr<t_obs_r))-1}else{start_t_p <- length(t_pr)}
  # add pseudo-data to observations, generating full dataset
  t_p_trunc <- t_p[1:end_t_p]
  t_pr_trunc <- t_pr[1:start_t_p]
  t_add <- NA
  t_add <- sort(c(t_p_trunc,t,t_pr_trunc))
  
  # return(t_add)
  
  
  
  
  # generate design points
  # design points cover the interval:
  # [tobs−l − 3 · [t1(n1) − t1(1)] · hrelmax; tobs−r + 3 · [t1(n1) − t1(1)] · hrelmax]
  # with constant spacing
  
  series_length <- length(t_obs_l:2018)
  hrelmax <- 0.5
  
  design_points <- round(t_obs_l - 3 * series_length * hrelmax):round(t_obs_r + 3 * series_length * hrelmax)
  
  # plot(design_points - t_add[10])
  # The occurrence rate, ly, is estimated using a kernel function (Diggle 1985):
  
  # ly(lx)=sum from 1 to n (􏰋Kh(lx −t(i)), Kh(·)=h−1K(·/h)
  
  # n <- length(design_points)
  # nt <- length(t_add)
  # ly <- data.frame(matrix(ncol=nt,nrow = n))
  # h <- smooth
  
  # flood_density <- density(x=t_add,bw=h,kernel="gaussian", from = min(design_points), to = max(design_points))
  y_label = parse(text=paste("'flood occurrence [yr'", "^-1 ", "*']'", sep=""))
  # recur_plot <- plot(flood_density[["x"]], h*flood_density[["y"]], xlim = c(t_obs_l,2020), "l", xlab = "year", 
  #                    ylab = y_label
  #                    # ylim=c(0, 0.01)
  #                   )
  # print(recur_plot)
  smooths <- smooths
  recur_df <- data.frame(matrix(ncol = length(smooths)+1, nrow = 512))
  colnames(recur_df) <- (c("year",smooths))
  
  for(j in 1:length(smooths)){
    bw_select = smooths[j]
    flood_density <- density(x=t_add,bw=bw_select,kernel="gaussian", from = min(design_points), to = max(design_points), weights = c(rep(1, length(t_add))))
    # could instead use built-in bandwidth estimator 'ucv' or 'bcv' (unbiased or biased cross-validation)
    # flood_density <- density(x=t_add,bw='ucv',kernel="gaussian", from = min(design_points), to = max(design_points), weights = c(rep(1, length(t_add))))
    recur_df[,(j+1)] <- flood_density[["y"]]
    if(j==1){recur_df$year <- flood_density[["x"]]}
  }
  
  
  
  recur_df_subset <- subset(recur_df,recur_df$year>t_obs_l & recur_df$year<t_obs_r)
  
  if(plot_binary==T){recur_plot_data <- melt(recur_df_subset,"year")
  recur_plot <- ggplot(data = recur_plot_data,aes(x = year, y = value, group = variable, color = variable)) + theme_bw() 
  recur_plot <- recur_plot + geom_line() + scale_color_viridis_d() + theme(panel.border = element_rect(colour = "black", size = 1, fill = NA),
                                                                           panel.grid.major = element_blank(),
                                                                           panel.grid.minor = element_blank())
  recur_plot <- recur_plot + xlim(t_obs_l, t_obs_r) + labs(
    x = "Year",
    y = y_label,
    colour = "smoothing [yrs]"
    # shape = "Drainage Area"
  )# par(new = T)
  print(recur_plot)}
  
  return(list("recur_df_subset" = recur_df_subset, "t_add" = t_add))
}
t_refl <- c(t)

# makes a data frame of false flood years
t_all_test <- data.frame(t_all$flood_count)
t_all_plot <- data.frame(matrix(nrow = nrow(t_all_test),ncol = 20))
for(i in 1:20){
  t_all_test_plot[,i] <- t_all_test[as.integer(resample_bootstrap(t_all_test)),]
}
t_all_test_plot$water_year <- t_all$water_year
plot(getAreaPlot(t_all_test_plot,"flood",bounds,"annual"))

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