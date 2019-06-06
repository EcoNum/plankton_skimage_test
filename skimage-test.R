# Note: with SciViews Box 2019, one needs to do first in a terminal:
# sudo pip3 install scikit-image==0.14.2
# R -e 'install.packages("piggyback")'

# Example sample for zooimage
smp <- here::here("zooimage-set", "STDCE.2005-01-18.H1.zidb")
# Zooscan dataset from https://www.seanoe.org/data/00446/55741/ en CC-BY-NC
features_skimage <- here::here("zooscan-set", "features_skimage.csv.gz")
features_native <- here::here("zooscan-set", "features_native.csv.gz")
taxa <- here::here("zooscan-set", "taxa.csv.gz")

repo <- "EcoNum/skimage-test"
data_tag <- "v1.0.0"

# Manage datasets that are stored as releases in the repository with piggyback
library(piggyback)
# Code in comment is required only once to upload the datasets into the release
#pb_track(c("*.csv.gz", # Zooscan datasets
#           "*.zidb"))  # Zooimage datasets
#usethis::browse_github_pat(scopes = c("repo", "gist"),
#  description = "R:GITHUB_PAT", host = "https://github.com")
#usethis::edit_r_environ()
#Restart R session
#Create a first release manually inside Github
#pb_new_release(repo, tag = data_tag, name = "Zooscan/Zooimage datasets")
#pb_upload(smp, repo = repo, tag = data_tag,
#  name = basename(smp), overwrite = FALSE)
#pb_upload(features_skimage, repo = repo, tag = data_tag,
#  name = basename(features_skimage), overwrite = FALSE)
#pb_upload(features_native, repo = repo, tag = data_tag,
#  name = basename(features_native), overwrite = FALSE)
#pb_upload(taxa, repo = repo, tag = data_tag,
#  name = basename(taxa), overwrite = FALSE)

# Download datasets into local repo
pb_download("STDCE.2005-01-18.H1.zidb",
  dest = here::here("zooimage-set"),
  repo = repo, tag = data_tag, overwrite = FALSE)
pb_download(
  c("features_skimage.csv.gz", "features_native.csv.gz", "taxa.csv.gz"),
  dest = here::here("zooscan-set"),
  repo = repo, tag = data_tag, overwrite = FALSE)


# Customized skimageVars(), initially from zooimage 5.6.0
skimageVars <- function(zidbfile) {
  library(zooimage)
  # Use reticulate to access Python and scikit-image
  library(reticulate)
  use_python("/usr/bin/python3")
  #use_virtualenv("skimage")
  skimage <- import("skimage")
  np <- import("numpy", convert = FALSE)
  
  # Read dataset
  #dat1 <- zidbDatRead(zidbfile)
  #head(dat1)
  #attr(dat1, "metadata")
  #summary(dat1[, c("Area", "Perim.", "Skew", "Kurt")])
  #plot(dat1$Area, dat1$Perim., xlab = "Area", ylab = "Perimeter")
  
  # Lazy loading data from one ZIDB file in R
  db1 <- zidbLink(zidbfile)
  
  # Get the list of all vignettes in this dataset
  items1 <- ls(db1) # Contains data in *_dat1 and vignettes in *_nn
  vigs1 <- items1[-grep("_dat", items1)]
  lvigs1 <- length(vigs1)
  if (!lvigs1)
    stop("No vignettes found in the file '", zidbfile, "'. Are you sure it is correct?")
  
  # Display a 5*5 thumbnail of the first 25 vignettes
  #zidbPlotNew("The 25 first vignettes in MTPS.2004-10-20.H1")
  #for (i in 1:25)
  #  zidbDrawVignette(db1[[vigs1[i]]], item = i, nx = 5, ny = 5)
  
  # Read one vignette at a time into R, and calculate skimage attributes
  for (i in 1:lvigs1) {
    # Note: (try 1 for good vignette -only one item- and
    # 18 for a bad vignette with two items)
    vig <- vigs1[i]
    png <- db1[[vig]]
    img <- png::readPNG(png, native = FALSE)
    #png::writePNG(img, target = "~/shared/projects/skimage/vig1.png")
  
    # Extract the three channels: mask, OD and visual
    # Red channel: OD image is inverted
    #png::writePNG(1 - img[ , , 1], target = "~/shared/projects/skimage/vig1_od.png")
    # Green channel: visual
    #png::writePNG(img[ , , 2], target = "~/shared/projects/skimage/vig1_visu.png")
    # Blue channel: the mask is greylevel 50, and the rest is visual again
    mask <- (img[ , , 3] != 50/255) + 0
    #png::writePNG(mask, target = "~/shared/projects/skimage/vig1_mask.png")
  
    # We convert the [0, 1] scale of png into [0, 255] and get it as Numpy array
    # But we first need integer inside our object
    mask2 <- as.integer(mask * 255L)
    dim(mask2) <- dim(mask)
    mask3 <- skimage$util$img_as_ubyte(np$array(mask2)) < 128
    # Note: connectivity is supposed to be mask3.ndim in Python, with value = 2L here
    label_mask3 <- skimage$measure$label(mask3, connectivity = 2L)
    # Again, in R, automatic conversion into numeric!
    dim_lm3 <- dim(label_mask3)
    label_mask3 <- as.integer(label_mask3)
    dim(label_mask3) <- dim_lm3
    # Get the inverted OD image (note: Zooscan images are 10 pixel higher => add 10 and limit to 255L)
    invod <- as.integer(pmin(255L, ((1 - img[ , , 1]) * 255L) + 10L))
    dim(invod) <- dim(img[ , , 1])
  
    # Measure the particles
    # This does not work for convex_area and solidity because an array
    # passed from R to Python is F-contiguous while a C-contiguous array is needed!
    #props <- skimage$measure$regionprops(label_image = label_mask3, intensity_image = invod)
    py$mask <- label_mask3
    py$invod <- invod
    py_run_string("import skimage; props = skimage.measure.regionprops(label_image = skimage.util.img_as_uint(mask.copy(order = 'C')), intensity_image = invod.copy(order = 'C'))", convert = FALSE)
    props <- py$props
    # In case we got several blobs, check which one is the right one
    item <- 1
    l <- length(props)
    # In case we got several items, check which one is better filling the area
    # (ZooImage increases the bbox by 150%, but sometimes, it fails because the
    # object is too close to the border(s)!)
    if (l > 1) {
      idim <- dim(mask)
      deltas <- numeric(0)
      for (j in 1:l) {
        bbox <- unlist(props[[j]]$bbox)
        deltas[[j]] <- (idim[1] - (bbox[3] - bbox[1]) * 1.5) +
          (idim[2] - (bbox[4] - bbox[2]) * 1.5)
      }
      item <- which.min(deltas)
    }
    prop <- props[[item]]
    
    # Now, get items (note: same columns as for the ZOoscan dataset)
    inertia_tensor <- as.numeric(prop$inertia_tensor)
    inertia_tensor_eigvals <- as.numeric(prop$inertia_tensor_eigvals)
    # Why are Zooscan moment 1e9 times higher???
    moments_hu <- as.numeric(prop$moments_hu) * 1e9
    moments_normalized <- as.numeric(prop$moments_normalized) * 1e9
    # Replace NaN by NA to be fully compatible with Zooimage dataset
    moments_normalized[is.nan(moments_normalized)] <- NA
    weighted_moments_hu <- as.numeric(prop$weighted_moments_hu) * 1e9
    weighted_moments_normalized <- as.numeric(prop$weighted_moments_normalized) * 1e9
    res1 <- data.frame(
      objid = vig,
      area = prop$area,
      convex_area = prop$convex_area,
      eccentricity = prop$eccentricity,
      equivalent_diameter = prop$equivalent_diameter,
      euler_number = prop$euler_number,
      filled_area = prop$filled_area,
      inertia_tensor0 = inertia_tensor[1],
      inertia_tensor1 = inertia_tensor[2],
      inertia_tensor2 = inertia_tensor[3],
      inertia_tensor3 = inertia_tensor[4],
      inertia_tensor_eigvals0 = inertia_tensor_eigvals[1],
      inertia_tensor_eigvals1 = inertia_tensor_eigvals[2],
      major_axis_length = prop$major_axis_length,
      max_intensity = prop$max_intensity,
      mean_intensity = prop$mean_intensity,
      min_intensity = prop$min_intensity,
      minor_axis_length = prop$minor_axis_length,
      moments_hu0 = moments_hu[1],
      moments_hu1 = moments_hu[2],
      moments_hu2 = moments_hu[3],
      moments_hu3 = moments_hu[4],
      moments_hu4 = moments_hu[5],
      moments_hu5 = moments_hu[6],
      moments_hu6 = moments_hu[7],
      moments_normalized0 = moments_normalized[1],
      moments_normalized1 = moments_normalized[2],
      moments_normalized2 = moments_normalized[3],
      moments_normalized3 = moments_normalized[4],
      moments_normalized4 = moments_normalized[5],
      moments_normalized5 = moments_normalized[6],
      moments_normalized6 = moments_normalized[7],
      moments_normalized7 = moments_normalized[8],
      moments_normalized8 = moments_normalized[9],
      moments_normalized9 = moments_normalized[10],
      moments_normalized10 = moments_normalized[11],
      moments_normalized11 = moments_normalized[12],
      moments_normalized12 = moments_normalized[13],
      moments_normalized13 = moments_normalized[14],
      moments_normalized14 = moments_normalized[15],
      moments_normalized15 = moments_normalized[16],
      perimeter = prop$perimeter,
      solidity = prop$solidity,
      weighted_moments_hu0 = weighted_moments_hu[1],
      weighted_moments_hu1 = moments_hu[2],
      weighted_moments_hu2 = moments_hu[3],
      weighted_moments_hu3 = moments_hu[4],
      weighted_moments_hu4 = moments_hu[5],
      weighted_moments_hu5 = moments_hu[6],
      weighted_moments_hu6 = moments_hu[7],
      weighted_moments_normalized0 = weighted_moments_normalized[1],
      weighted_moments_normalized1 = weighted_moments_normalized[2],
      weighted_moments_normalized2 = weighted_moments_normalized[3],
      weighted_moments_normalized3 = weighted_moments_normalized[4],
      weighted_moments_normalized4 = weighted_moments_normalized[5],
      weighted_moments_normalized5 = weighted_moments_normalized[6],
      weighted_moments_normalized6 = weighted_moments_normalized[7],
      weighted_moments_normalized7 = weighted_moments_normalized[8],
      weighted_moments_normalized8 = weighted_moments_normalized[9],
      weighted_moments_normalized9 = weighted_moments_normalized[10],
      weighted_moments_normalized10 = weighted_moments_normalized[11],
      weighted_moments_normalized11 = weighted_moments_normalized[12],
      weighted_moments_normalized12 = weighted_moments_normalized[13],
      weighted_moments_normalized13 = weighted_moments_normalized[14],
      weighted_moments_normalized14 = weighted_moments_normalized[15],
      weighted_moments_normalized15 = weighted_moments_normalized[16]
    )
    # What about centroid, weighted_centroid, moments, moments_central, orientation,   and weighted_moments?
    if (i == 1) res <- res1 else res <- rbind(res, res1)
  }
  res
}

# Calculate skimage attributes for smp:
skim <- skimageVars(smp)

# Read 10 first items of skimage attributes for the Zooscan dataset as a comparison
#library(data.io)
#skim_ex <- read(features_skimage, n_max = 10, guess_max = 10)
skim_ex <- readr::read_csv(features_skimage, n_max = 10, guess_max = 10)

# - Intensities are very different!
# - Moments and Weighted moments are also very different!

# Look for most similar items in both sets
skim_all <- readr::read_csv(features_skimage)

# Search for the neirest neightbourg is done using FNN::get.knnx() or RANN::nn2()
# We need to prepare two matrices and keep only similar columns
keep <- c("area", "convex_area", "eccentricity", "equivalent_diameter",
  "euler_number", "filled_area", "inertia_tensor0", "inertia_tensor1",
  "inertia_tensor2", "inertia_tensor3", "inertia_tensor_eigvals0",
  "inertia_tensor_eigvals1", "perimeter", "solidity")

skim1 <- skim[, keep]
skim1b <- skim1[1:3, ]
skim2 <- skim_all[, keep]

# Both functions are equivalent in speed et reutrn same results
system.time(RANN::nn2(skim2, skim1b, k = 10))
system.time(FNN::get.knnx(skim2, skim1b, k = 10))
# Warning: FNN seems to be UNMAINTAINED or now!
near <- RANN::nn2(skim2, skim1, k = 10)
near_idx <- near$nn.idx
near_dist <- near$nn.dist
nearest <- (1:nrow(skim1))[order(near_dist[, 1])]
nearest[1:30]

# OK, so, let's take the best ones... and look at them!
# I have tried #1 (background artifact), #2 (egg), and #21 (copepod) so far
pos <- 21
item1_1 <- nearest[pos]
items1_2 <- near_idx[item1_1, ]
best1 <- rbind(skim1[item1_1, ], skim2[items1_2, ])
# Get item from zidb, get item from the zoosan dataset and get full attributes
best1full <- rbind(skim[item1_1, ], skim_all[items1_2, ])
# Inspect best1full table...
View(best1full)
# Everything seems OK except weighted_moments_hu

# This is the corresponding vignette in the zooimage dataset
db1 <- zidbLink(smp)
vig1 <- db1[[best1full[1, 1]]]
img <- png::readPNG(vig1, native = FALSE)
png::writePNG(img, target = here::here("zooimage-set",
  paste0("vig", pos, ".png")))


# What next?
# 1) Consolidate training/testing datasets with common native + common skimage
#    (for skimage, need to take weighter_moments_huX out)
# 2) rescale moments: sign(X) * log(abs(X/1e9)), except for moment_hu6 where
#    the sign is not relevant for us => log(abs(moments_hu6/1e6))
# 3) Calculate meaningful descriptors. The idea is to consider representative
#    items for various shapes (egg, bubble, copepod1, chaetognath, ...)
#    and to use centroid (mean of all
#    7 hu moments in the training set for these items). One could do much
#    better later on by considering the codebook of LVQ as the basis of
#    calculation, but let's keep it simple for now
#    See https://www.learnopencv.com/shape-matching-using-hu-moments-c-python/
# 4) With all attrobutes + these egg-like, bubble-like, etc. shape attributes
#    let's do a classification with random forest and let's study the importance
#    of the variables, so that we could see if skimage improves classification.
