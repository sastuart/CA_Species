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

# Test: download occurence data for a few test species ----

# Set up string to search for species name
t_means_angio[1:5,1:3]
i <- 1
genus <- gsub(" .*", "\\1", t_means_angio$binomial[i])
spepi <- paste(gsub(".* ", "\\1", t_means_angio$binomial[i]), "*", sep="")

# Download occurence data from gbif

# * gets species names that inc author and subspecies
test1_occ <- gbif(genus, spepi, geo=FALSE) 
# remove data without geolocations
# not sure why Hijmans & Elith do this; might be unecessary 
test1_geo <- subset(test1_occ, !is.na(lon) & !is.na(lat))
## compare the above results to gbif(genus, spepi, geo=TRUE) to see if they are different ---------
comp_test <- gbif(genus, spepi, geo=TRUE) 
comp_test == test1_geo
comp_test == test1_occ
length(which(is.na(comp_test$lon) & is.na(comp_test$lat)))
length(which(is.na(comp_test$lon)))
length(which(is.na(comp_test$lat)))
length(which(is.na(test1_occ$lon) & is.na(test1_occ$lat)))
length(which(is.na(test1_geo$lon) & is.na(test1_geo$lat)))
# Hmm, so the geo = TRUE option seems to not really be working.

# find duplicated records
dups2 <- duplicated(test1_geo[, c('lon', 'lat')])
# keep the records that are _not_ duplicated
test1_unigeo <- test1_geo[!dups2, ]

## Remove outlier records ----
# turn the records into a spatial points data frame
library(sp)
coordinates(test1_unigeo) <- ~lon+lat
# give the data some projection information CAUTION: Should double-check that this
# is the right way to assign a projection
# proj4string(test1_unigeo) <- CRS("+proj=longlat +ellps=WGS84")
proj4string(test1_unigeo) <- CRS("+proj=longlat")
# re-project to match the cfp_proxy
test1_unigeo <- spTransform(test1_unigeo, CRS(proj4string(cfp_proxy))) # transform CRS
# this is designed for a batch process, so
# we will simply crop to the cfp and/or cfp_proxy map
test1_unicfp <- test1_unigeo[cfp_proxy,]
plot(cfp_proxy)
points(test1_unicfp, col=alpha('red', 0.1), pch=20, cex=0.5)

i <- 10
genus <- gsub(" .*", "\\1", t_means_angio$binomial[i])
spepi <- paste(gsub(".* ", "\\1", t_means_angio$binomial[i]), "*", sep="")
test2_occ <- gbif(genus, spepi, geo=FALSE) 
test2_geo <- subset(test2_occ, !is.na(lon) & !is.na(lat))
dups2 <- duplicated(test2_geo[, c('lon', 'lat')])
test2_unigeo <- test2_geo[!dups2, ]
coordinates(test2_unigeo) <- ~lon+lat
proj4string(test2_unigeo) <- CRS("+proj=longlat")
test2_unigeo <- spTransform(test2_unigeo, CRS(proj4string(cfp_proxy))) # transform CRS
test2_unicfp <- test2_unigeo[cfp_proxy,]
points(test2_unicfp, col=alpha('blue', 0.10), pch=20, cex=0.5)

i <- length(t_means_angio$binomial)
genus <- gsub(" .*", "\\1", t_means_angio$binomial[i])
spepi <- paste(gsub(".* ", "\\1", t_means_angio$binomial[i]), "*", sep="")
test3_occ <- gbif(genus, spepi, geo=FALSE) 
test3_geo <- subset(test3_occ, !is.na(lon) & !is.na(lat))
dups2 <- duplicated(test3_geo[, c('lon', 'lat')])
test3_unigeo <- test3_geo[!dups2, ]
coordinates(test3_unigeo) <- ~lon+lat
proj4string(test3_unigeo) <- CRS("+proj=longlat")
test3_unigeo <- spTransform(test3_unigeo, CRS(proj4string(cfp_proxy))) # transform CRS
test3_unicfp <- test3_unigeo[cfp_proxy,]
points(test3_unicfp, col=alpha('darkgreen', 0.10), pch=20, cex=0.5)

i <- length(t_means_angio$binomial) - 20
genus <- gsub(" .*", "\\1", t_means_angio$binomial[i])
spepi <- paste(gsub(".* ", "\\1", t_means_angio$binomial[i]), "*", sep="")
test4_occ <- gbif(genus, spepi, geo=FALSE) 
test4_geo <- subset(test4_occ, !is.na(lon) & !is.na(lat))
dups2 <- duplicated(test4_geo[, c('lon', 'lat')])
test4_unigeo <- test4_geo[!dups2, ]
coordinates(test4_unigeo) <- ~lon+lat
proj4string(test4_unigeo) <- CRS("+proj=longlat")
test4_unigeo <- spTransform(test4_unigeo, CRS(proj4string(cfp_proxy))) # transform CRS
test4_unicfp <- test4_unigeo[cfp_proxy,]
points(test4_unicfp, col=alpha('pink', 0.10), pch=20, cex=0.5)

plot(cfp_proxy)
points(test1_unicfp, col=alpha('red', 0.1), pch=20, cex=0.5)
points(test2_unicfp, col=alpha('blue', 0.10), pch=20, cex=0.5)
points(test3_unicfp, col=alpha('darkgreen', 0.10), pch=20, cex=0.5)
points(test4_unicfp, col=alpha('purple', 0.10), pch=20, cex=0.5)
legend(x="topleft", legend=t_means_angio$binomial[1,10,130,110], fill=c()

store1<- list(test1_unicfp, test2_unicfp, test3_unicfp)
store1[[3]]$lat
store1[[3]]$lon
# works

store_occ <- list(test1_unicfp)
store_occ <- list(store_occ, test2_unicfp)
store_occ <- list(store_occ, test3_unicfp)
store_occ <- list(store_occ, test4_unicfp)
# this doesn't work beacuse each time the list is just two elements

store_occ <- list(length=4)
store_occ[1] <- test1_unicfp
store_occ[2] <- test2_unicfp
store_occ[3] <- test3_unicfp
store_occ[4] <- test4_unicfp
store_occ[5] <- test4_unicfp #Interesting, can continue to load data beyond index

store_occ[[4]]$lat == test4_unicfp$lat
# works :)

# Will's rgbif funtions ----------
library(dplyr)
library(rgbif)

remove_all_issues_except<-function(out,issue_codes_were_ok_with){
  for (i in 1:length(issue_codes_were_ok_with)){
    out$data$issues<-sub(issue_codes_were_ok_with[i],"",out$data$issues)
  }
  out$data<-filter(out$data,issues=="")
  return(out$data)
}

species <- "Quercus douglasii"
species <- "Quercus durata"

str(out_err)
str(out_ok)
out_err$meta$count
out_ok$meta$count

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
## Let's see how my data is different from Will's ----

gbif_issues()

test1_will<-get_species_observations(paste(t_means_angio$binomial[1]))
names(test1_occ) == names(test1_will)
match(names(test1_will), names(test1_occ))
# in test1_will, not in dismo-function data
names(test1_will)[which(is.na(match(names(test1_will), names(test1_occ))))]
# in dismo-function data, not in test1_will
names(test1_occ)[which(is.na(match(names(test1_occ), names(test1_will))))]

plot(cfp_proxy)
points(test1_unicfp, col=alpha('red', 0.1), pch=20, cex=0.5)
points(test1_will$decimalLatitude~test1_will$decimalLongitude, col=alpha('orange', 0.1), pch=20, cex=0.5)

test1_will_geo <- subset(test1_will, !is.na(decimalLongitude) & !is.na(decimalLatitude))
dups2 <- duplicated(test1_will_geo[, c('decimalLongitude', 'decimalLatitude')])
test1_will_unigeo <- test1_will_geo[!dups2, ]
coordinates(test1_will_unigeo) <- ~decimalLongitude+decimalLatitude
proj4string(test1_will_unigeo) <- CRS("+proj=longlat")
test1_will_unigeo <- spTransform(test1_will_unigeo, CRS(proj4string(cfp_proxy))) # transform CRS
test1_will_unicfp <- test1_will_unigeo[cfp_proxy,]
points(test1_will_unicfp, col=alpha('darkblue', 0.10), pch=20, cex=0.5)

plot(cfp_proxy)
points(test1_unicfp, col='red', pch=20, cex=0.5)
points(test1_will_unicfp, col='green', pch=20, cex=0.5)

duplicated(test1_unigeo[,c("lon","lat")])

which(test1_will$name != "Acacia greggii")
test2_will<-get_species_observations(paste(t_means_angio$binomial[10]))
test2_will$name
testacer_will<-get_species_observations("Acer saccharum")
testacer_will$name
points(testacer_will$decimalLatitude~testacer_will$decimalLongitude, col=alpha('purple', 0.1), pch=20, cex=0.5)

testam_will<-get_species_observations("Acer macrophyllum")
testam_will$name
points(testam_will$decimalLatitude~testam_will$decimalLongitude, col=alpha('yellow', 0.1), pch=20, cex=0.5)

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

# check out how it's working
spp_cols <- rainbow(length(occur_list))
png(filename = "130_spp.jpg",width = 500,height = 500)
plot(cfp_proxy)
for(i in 1:length(occur_list)){
  points(occur_list[[i]]$decimalLatitude~occur_list[[i]]$decimalLongitude, col=alpha(spp_cols[i], 0.1), pch=20, cex=0.5)  
}
dev.off()

# Test: generate background data for a few test species ----

# Start by generating a background dataset that combines the occurrence
# records for all species

all_occurs <- data.frame(decimalLongitude=NULL, decimalLatitude=NULL)
for(i in 1:length(occur_list)) {
  all_occurs <- rbind(all_occurs,coordinates(occur_list[[i]]))
}

# Turn all_occurs into a spatial points dataframe and set its projection
coordinates(all_occurs) <- ~decimalLongitude+decimalLatitude
proj4string(all_occurs) <- CRS("+proj=longlat")

# Lets take a look
plot(cfp_proxy)
points(all_occurs, pch=20, cex=0.1)

# Ok, CA is quite well sampled, OR is well-sampled, NV and Baja are not well sampled
# But there is definitely spatial bias in all four regions.

# Let's test the code on a smaller subset
some_occurs <- data.frame(decimalLongitude=NULL, decimalLatitude=NULL)
for(i in 1:2) {
  some_occurs <- rbind(some_occurs,coordinates(occur_list[[i]]))
}

# Turn some_occurs into a spatial points dataframe and set its projection
coordinates(some_occurs) <- ~decimalLongitude+decimalLatitude
proj4string(some_occurs) <- CRS("+proj=longlat")

# The first two species have light coverage but do reflect the overall pattern
# of spatial bias in the sample
plot(cfp_proxy)
points(some_occurs, pch=20, cex=0.1)

# Next we will generate a circular model around all of these datapoints then
# melt the circles together, so that all the space around the datapoints is
# represented equally.
library(rgeos)
target_bg_circles <- circles(some_occurs, d=25000, lonlat=TRUE)
# What size are these circles?
target_bg_shape <- gUnaryUnion(target_bg_circles@polygons)

# Trim out the edges of circles that are outside the CA floristic province
target_bg_cfp <- intersect(target_bg_shape, cfp_proxy)
# Warning that CRS are non-identical, even though:
identicalCRS(target_bg_shape,some_occurs) #TRUE
proj4string(target_bg_shape) == proj4string(some_occurs) #TRUE
# Hmmm... ---------


# What did we get?
plot(cfp_proxy)
plot(target_bg_cfp, add=TRUE, col='red')
# So far so good.

# Create a grid to define sampling within the cfp_proxy area
# I think what H&E are doing us using the resoultion
# of the BioClim rasters to set the grid size
# So, we need to get the Bioclim rasters
b_clim <- getData("worldclim", var="bio", res=2.5, path="data/")

cfp_grid <- raster(cfp_proxy)
bg_random <- randomPoints(cfp_grid, 500)

plot(!is.na(cfp_grid), legend=FALSE)

target_bg_1 <- spsample(target_bg_shape, 250, type='random', iter=25)

# Loop: generate background data for all species ----
# Maxent model: precip, max temp, min temp ----
# Extract margins from maxent model ----