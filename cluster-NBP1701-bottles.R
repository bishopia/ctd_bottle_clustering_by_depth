# library(naniar)
library(tidyverse)
library(purrr)
library(cluster)

rm(list=ls())

# setwd("/Volumes/GoogleDrive/My Drive/LAB PROJECTS/NBP17_cruise_data/STATION_METADATA")
setwd("~/Library/CloudStorage/GoogleDrive-bishopia@uri.edu/My Drive/LAB PROJECTS/NBP17_cruise_data/STATION_METADATA")

#import data
raw <- read_csv("raw/NBP1701_320620161224.exc.csv", skip=50) #import data with headers and units row
units <- raw[1,]
tmp <- raw[-1,] #delete units row
tmp2 <- tmp[-nrow(tmp),] #delete last row, which 

#fix character to numeric
tmp3 <- tmp2 %>% 
  mutate_at(c("CTDPRS", "CTDSAL", "CTDTMP", "CTDOXY"), as.numeric)

#replace various values that signify NA with NA
na_strings <- c("-999", "-999.0", "-999.00", "-999.000", "-999.0000")
df <- tmp3 %>%
  replace_with_na_all(condition = ~.x %in% na_strings)
rm(tmp, tmp2, tmp3)

#import station region mapping file
station_regions <- read_csv("raw/station_region_mapping.csv")

#prepare df
df2 <- df %>% left_join(station_regions) %>% #add regions
  rename(station=STNNBR, cast=CASTNO) %>% #rename vars
  filter(cast==1) %>% # FOR THIS CRUISE, WE ONLY TOOK FROM FIRST CAST ONE STATION FOR ISOLATES, AND CAST TWO ON STATION WAS GRADIENT, NOT AMENABLE FOR CLUSTERING BOTTLES
  mutate(ID=paste(station, cast, sep="-")) %>% #make ID variable
  arrange(station, cast, BTLNBR) %>% #arrange
  drop_na(CTDPRS, CTDOXY, CTDTMP, CTDPRS) #remove na values from these variables, messes up clustering
rm(df, station_regions)

#list of station/cast groups; 
#IMPORTANT: just use a single grouping variable here. if you have two, make a new ID variable that combines both, that way you only need one for loop
groups <- unique(df2$ID)

#how many ks to iterate over?
max_k <- 5
#max silhouette width?
max_sil_width=.5

#initialize depth cluster variable
depth_cluster <- NULL

#initialize
df_list <- list()

for (i in seq_along(groups)) {
# for (i in 1:9) {

    #pull one group via ID
    tmp <- df2[df2$ID==groups[i],]
    
    # tmp <- df2 %>% filter(station==stations[i], cast==casts[j])
    
    #initialize silhouette width vector
    sil_widths <- numeric(max_k)

    #loop over k values, use CTD data to find clusters
    for(k in 2:max_k) {
      kmeans_cluster <- kmeans(tmp[,c("CTDPRS", "CTDOXY", "CTDSAL", "CTDTMP")], centers=k)
      sil_width <- silhouette(kmeans_cluster$cluster, dist(tmp[,c("CTDPRS", "CTDOXY", "CTDSAL", "CTDTMP")]))
      sil_widths[k] <- mean(sil_width[,3])
    }
    
    #choose max sil
    num_clusters <- which.max(sil_widths)
    if (sil_widths[num_clusters] < max_sil_width) {
      num_clusters <- num_clusters - 1
    }
    
    #run kmeans with that k value
    kmeans_cluster <- kmeans(tmp$CTDPRS, centers=num_clusters)
    
    #pull group assignment vector
    DEPTH_GROUP <- as.factor(kmeans_cluster$cluster)
    
    #add group assignment to tmp df
    data <- data.frame(tmp, DEPTH_GROUP)
    
    #add df to list of lists
    df_list <- append(df_list, list(data))
    
    rm(tmp, sil_widths, kmeans_cluster, sil_width, num_clusters, DEPTH_GROUP, data)
}

#combine/collapse list of dfs into one df
combined_df <- bind_rows(df_list)
rm(df_list)

#filter what you want to look at
short_list <- combined_df %>% dplyr::select(station, cast, CTDPRS, CTDOXY, CTDSAL, CTDTMP, DEPTH_GROUP)

#export csv
write_csv(combined_df, "outputs/first_cast_per_station_CTD_data_20230217.csv")
