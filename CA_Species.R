# CA CPCE Species Distribution Modelling

library(dplyr)
library(rgbif)
library(raster)
library(scales)
library(dismo)
library(maptools)

# Set up map of California & neighboring states to trim to ----
ca.map <- getData('GADM', country="USA", level=1)
mex.map <- getData('GADM', country="MEX", level=1)
ca.neigh <- c("Oregon","California","Nevada") 
ca.map.ne <- ca.map[match(toupper(ca.neigh),toupper(ca.map$NAME_1)),]
mex.neigh <- c("Baja California")
mex.map.ne <- mex.map[match(toupper(mex.neigh),toupper(mex.map$NAME_1)),]
plot(mex.map.ne)
plot(ca.map.ne)
cfp_proxy <- spRbind(ca.map.ne,mex.map.ne)
plot(cfp_proxy)

# Load species list from aggregated data file
t_means <- read.csv("CA CPCE Functional Traits Taxon Means.csv")
t_means_angio <- subset(t_means, group == "Angiosperm")

# GBIF interface funtions ----------

remove_all_issues_except<-function(out,issue_codes_were_ok_with){
  for (i in 1:length(issue_codes_were_ok_with)){
    out$data$issues<-sub(issue_codes_were_ok_with[i],"",out$data$issues)
  }
  out$data<-filter(out$data,issues=="")
  return(out$data)
}

get_species_observations<-function(species="Arctostaphylos densiflora",
                                   issue_codes_were_ok_with=c("gass84"),
                                   max_records=5000){
  require(dplyr)
  key <- name_suggest(q=species, rank='species')$key[1]
  out <- occ_search(taxonKey=key, limit=max_records)
  if (out$meta$count == 0){
    # If the first name doesn't match any records, try the second
    key <- name_suggest(q=species, rank='species')$key[2]
    # If the second name doesn't match any records, give up
    out <- occ_search(taxonKey=key, limit=max_records)
    if (out$meta$count == 0){
    stop(paste("No records were found for",species,sep=" "))
  }}
  filtered.data<-remove_all_issues_except(out,issue_codes_were_ok_with)
  return(filtered.data)
}
# Loop: download all species from file ----
# Ok, so we have our basic set-up, now need to make it work in a loop.

name_data <- t_means_angio$binomial
occur_list <- list(length=length(name_data))
for (i in 1:length(name_data)){
  occur_raw <- get_species_observations(paste(name_data[i]))
  if(any(occur_raw$kingdom != "Plantae")){
    stop(paste("Name matching has failed, found taxon from kingdom",
               occur_raw$kingdom[which(occur_raw$kingdom != "Plantae")[1]],
               sep = " "))
  }
  print(paste(i,"of",length(name_data),occur_raw$name[1], sep=" "))
  occur_geo <- subset(occur_raw, !is.na(decimalLongitude) & !is.na(decimalLatitude))
  dups <- duplicated(occur_geo[, c('decimalLongitude', 'decimalLatitude')])
  occur_unigeo <- occur_geo[!dups, ]
  coordinates(occur_unigeo) <- ~decimalLongitude+decimalLatitude
  proj4string(occur_unigeo) <- CRS("+proj=longlat")
  occur_unigeo <- spTransform(occur_unigeo, CRS(proj4string(cfp_proxy))) # transform CRS
  occur_unicfp <- occur_unigeo[cfp_proxy,]
  occur_list[i] <- occur_unicfp  
}

# Show locations of all the species
spp_cols <- rainbow(length(occur_list))
png(filename = "130_spp.jpg",width = 500,height = 500)
plot(cfp_proxy)
for(i in 1:length(occur_list)){
  points(occur_list[[i]]$decimalLatitude~occur_list[[i]]$decimalLongitude, col=alpha(spp_cols[i], 0.1), pch=20, cex=0.5)  
}
dev.off()


# Now, download some economic data to compare this with
library(choroplethr)
library(acs)

# initialize the Amercian Community Survey (Census) API key
api.key.install("93dd8dd250c93e7df313a771e7b1ab0cb19f1485")
png(filename = "CANV_median_incomes.jpg",width = 500,height = 500)
  county_choropleth_acs("B19301", state_zoom=c("california", "nevada"))
dev.off()




# Plot locations of species ----
spp_cols <- black
png(filename = "130_spp.jpg",width = 500,height = 500)
plot(cfp_proxy)
for(i in 1:10){##length(occur_list)){
  points(occur_list[[i]]$decimalLatitude~occur_list[[i]]$decimalLongitude,
         col=alpha(spp_cols[i], 0.1), pch=20, cex=0.5)  
}
dev.off()

ca.counties <- getData('GADM', country="USA", level=2)
plot(ca.counties)


