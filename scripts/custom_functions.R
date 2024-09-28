# ---------------------------------------------------------------------------- #
# Script Name: <name_of_the_script>.R
# Author: Eleftherios Zormpas
# Date: <date_of_creation_or_modification>
#
# Description:
#   - Custom functions to use in the analysis scripts.
#   - Saves some repetition in coding.
#
# Dependencies:
#   - <List of required R packages (e.g. Seurat, scran, spatialDE)>
#   - <Any external scripts or tools used>
#
# Instructions:
#   - How to run the script, including any required parameters or settings.
#   - Example command or execution method if applicable.
#
# Version History:
#   - <Version updates, e.g. added new functionality, fixed bugs>
# ---------------------------------------------------------------------------- #

# The preprocessing here

# ---------------------------------------------------------------------------- #
# ----Custom functions to use --------------------------------------------------
# A custom function to save plots
savePlot <- function(plot,
                     out,
                     prfx = "",
                     main = "",
                     sfx = "",
                     other = "",
                     device = "png",
                     width = grDevices::dev.size(units = "in")[1],
                     height = grDevices::dev.size(units = "in")[2],
                     units = "in",
                     dpi = 300) {
  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".", device),
         plot = plot,
         device = device,
         width = width,
         height = height,
         units = units,
         dpi = dpi)
}

# A custom function to reload images when re-loading an SFE object from disk
reloadSFEimage <- function(msfe, sampleNames) {
  for (id in sampleNames) {
    sfe <- getSFE(msfe, id)
    sfe@int_metadata$imgData <- NULL
    ## Get scale factors
    scaleF <- jsonlite::fromJSON(txt = file.path(sampleDir[[id]],
                                                 "outs/spatial",
                                                 "scalefactors_json.json"))
    ## Add image
    sfe <- SpatialFeatureExperiment::addImg(sfe,
                                            file.path(sampleDir[[id]],
                                                      "outs/spatial/tissue_lowres_image.png"),
                                            sample_id = id,
                                            image_id = "lowres",
                                            scale_fct = scaleF[["tissue_lowres_scalef"]])
    ## Mirror the image
    sfe <- SpatialFeatureExperiment::mirrorImg(sfe,
                                               sample_id = id,
                                               image_id = "lowres")
    ## Add back into the MSFE
    msfe <- addSFE(msfe, sfe, id)

    ## Housekeeping
    rm(sfe)
  }
}
