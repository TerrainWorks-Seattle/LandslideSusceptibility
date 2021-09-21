referenceRaster <- terra::rast("C:/Work/netmapdata/Scottsburg/elev_scottsburg.flt")

inputRasters <- list(
  terra::rast("C:/Work/netmapdata/Scottsburg/grad_30.tif"),
  terra::rast("C:/Work/netmapdata/Scottsburg/plan_15.tif")
)

inputRasters <- TerrainWorksUtils::alignRasters(referenceRaster, inputRasters)

initiationPoints <- terra::vect("C:/Work/netmapdata/Scottsburg/Scottsburg_Upslope.shp")
initiationPoints <- terra::project(initiationPoints, referenceRaster)

initiationAreas <- terra::buffer(initiationPoints, width=30)

inputRaster <- c(inputRasters[[1]], inputRasters[[2]])

initiationValues <- terra::extract(inputRaster, initiationAreas)
