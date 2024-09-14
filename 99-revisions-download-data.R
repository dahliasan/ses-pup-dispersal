# Revisions for reviewers comments 


# Download subsurface data from copernicus --------------------------------

# check if recticulate installed else install and load it
library(reticulate)

# Create virtual env
virtualenv_create(envname = "CopernicusMarine")

virtualenv_install("CopernicusMarine", packages = c("copernicusmarine"))

reticulate::use_virtualenv("CopernicusMarine", required = TRUE)

# Import copernicusmarine package
cmt <- import("copernicusmarine")

# Set login credentials
# cmt$login("dfoo1", "soAFO9UK1HllL22%fa%j")

cmt$subset(
  dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
  variables=list("uo", "vo"),
  minimum_longitude=110,
  maximum_longitude=234,
  minimum_latitude=-76.943,
  maximum_latitude=-29.64961,
  start_datetime="1995-11-15T00:00:00",
  end_datetime="2001-10-01T00:00:00",
  minimum_depth=0.5,
  maximum_depth=0.5,
)

# import copernicusmarine
# 
# copernicusmarine.subset(
#   dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
#   variables=["uo", "vo"],
#   minimum_longitude=110,
#   maximum_longitude=234,
#   minimum_latitude=-76.943,
#   maximum_latitude=-29.64961,
#   start_datetime="1995-11-15T00:00:00",
#   end_datetime="2001-10-01T00:00:00",
#   minimum_depth=200,
#   maximum_depth=200,
# )



