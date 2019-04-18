# Libraries
packs = c("leaflet","sf")
success = suppressWarnings(sapply(packs, require, character.only = TRUE, quietly = TRUE))
if(any(!success)){
  print(paste("The packages", paste(names(success)[!success], collapse=", "), "are required for this process, but are not installed"))
  for(i in names(success)[!success]){
    y = readline(prompt = paste("Would you like to install the package ", i, "? Please enter 'yes' or 'no': ", sep=""))
    if(y == "yes"){
      install.packages(i, quiet = TRUE, verbose = FALSE)
      suppressWarnings(require(i, character.only = TRUE, quietly = TRUE))
    } else{
      stop(paste("The simulation cannot run if you do not install ", i, "! We hope you reconsider", sep=""))
    }
  }
}

# LTER Site Bounding Boxes

lter_sfc = st_sfc(
  st_polygon(
    list(
      cbind(
        c(-149.8455917, -149.829821,  -149.829821,  -149.8455917, -149.8455917),
        c(-17.47185366, -17.47185366, -17.48641792, -17.48641792, -17.47185366)
      )
    )
  ),
  st_polygon(
    list(
      cbind(
        c(-149.8116849, -149.7961685, -149.7961685, -149.8116849, -149.8116849),
        c(-17.46576169, -17.46576169, -17.48131958, -17.48131958, -17.46576169)
      )
    )
  ),
  st_polygon(
    list(
      cbind(
        c(-149.7708619, -149.7519968, -149.7519968, -149.7708619, -149.7708619),
        c(-17.50382025, -17.50382025, -17.52087158, -17.52087158, -17.50382025)
      )
    )
  ),
  st_polygon(
    list(
      cbind(
        c(-149.7772857, -149.7566866, -149.7566866, -149.7772857, -149.7772857),
        c(-17.53305021, -17.53305021, -17.55064263, -17.55064263, -17.53305021)
      )
    )
  ),
  st_polygon(
    list(
      cbind(
        c(-149.8869755, -149.8561009, -149.8561009, -149.8869755, -149.8869755),
        c(-17.56818162, -17.56818162, -17.59182383, -17.59182383, -17.56818162)
      )
    )
  ),
  st_polygon(
    list(
      cbind(
        c(-149.934537,  -149.9115336, -149.9115336, -149.934537,  -149.934537),
        c(-17.50735955, -17.50735955, -17.52839766, -17.52839766, -17.50735955)
      )
    )
  ),
  crs = 4326
)

lter_sites = st_sf(
  lter_sfc,
  site = sapply(1:6, function(x){paste("LTER", x)})
)

st_write(lter_sites, "Data/Shapefiles/LTER_Site_Outlines.shp")

# Plotting
#
# Plot LTER Site Bounding Boxes (base)
# plot(lter_sites)
# 
# Plot LTER Site Bounding Boxes (leaflet)
# m = leaflet(lter_sites) %>% 
#   addTiles() %>% 
#   addPolygons(popup = lter_sites$site)
# m
