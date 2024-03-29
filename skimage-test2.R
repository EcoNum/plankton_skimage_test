# Hu moments

huVars <- function(df) {
  df$moments_hu0 <- sign(df$moments_hu0) * log(abs(df$moments_hu0))
  df$moments_hu1 <- sign(df$moments_hu1) * log(abs(df$moments_hu1))
  df$moments_hu2 <- sign(df$moments_hu2) * log(abs(df$moments_hu2))
  df$moments_hu3 <- sign(df$moments_hu3) * log(abs(df$moments_hu3))
  df$moments_hu4 <- sign(df$moments_hu4) * log(abs(df$moments_hu4))
  df$moments_hu5 <- sign(df$moments_hu5) * log(abs(df$moments_hu5))
  df$moments_hu6 <- log(abs(df$moments_hu6))
  df
}

shapeVars <- function(df, ref) {
  # df: a dataframe with attributes, including Hu moments
  # ref: a reference table for additional variables (mean Hu moments per shape)
  for (i in 1:nrow(ref)) {
    ref0 <- ref[i, ]
    df[[paste0(ref0$taxon, "_like")]] <- sqrt(
      (df$moments_hu0 - ref0$mean_hu0)^2 +
      (df$moments_hu1 - ref0$mean_hu1)^2 +
      (df$moments_hu2 - ref0$mean_hu2)^2 +
      (df$moments_hu3 - ref0$mean_hu3)^2 +
      (df$moments_hu4 - ref0$mean_hu4)^2 +
      (df$moments_hu5 - ref0$mean_hu5)^2 +
      (df$moments_hu6 - ref0$mean_hu6)^2
    )
  }
  df
}


# Zooimage example dataset
smp <- here::here("zooimage-set", "STDCE.2005-01-18.H1.zidb")
# Zooscan dataset from https://www.seanoe.org/data/00446/55741/ en CC-BY-NC
features_skimage <- here::here("zooscan-set", "features_skimage.csv.gz")
features_native <- here::here("zooscan-set", "features_native.csv.gz")
taxa <- here::here("zooscan-set", "taxa.csv.gz")

skim_all <- readr::read_csv(features_skimage)
taxa_all <- data.io::read(taxa)
library(dplyr)
library(flow)
left_join(skim_all, taxa_all) %>.%
  huVars(.) %>.%
  select(., objid, area, moments_hu0:moments_hu6, taxon) -> zooscan_all

# Calculating mean values for moments_hu0 -> 6 by taxon
zooscan_all %>.%
  group_by(., taxon) %>.%
  summarise(.,
    mean_hu0 = mean(moments_hu0),
    mean_hu1 = mean(moments_hu1),
    mean_hu2 = mean(moments_hu2),
    mean_hu3 = mean(moments_hu3),
    mean_hu4 = mean(moments_hu4),
    mean_hu5 = mean(moments_hu5),
    mean_hu6 = mean(moments_hu6)
  ) -> zooscan_ref

zooscan_ref_df <- as.data.frame(zooscan_ref)
rownames(zooscan_ref_df) <- zooscan_ref_df$taxon
zooscan_ref_df$taxon <- NULL
zooscan_dist <- vegan::vegdist(zooscan_ref_df, method = "euclidean")
plot(hclust(zooscan_dist, method = "complete"))

# Create new variables for some groups
zooscan_all <- shapeVars(zooscan_all, zooscan_ref[1:3, ])

# What does it give?
library(chart)
zooscan_all %>.%
  group_by(., taxon) %>.%
  sample_n(., 100, replace =  TRUE) -> zooscan_all50
#zooscan_all50$taxon[!zooscan_all50$taxon %in% c("Acantharea", "Acartiidae", "Actinopterygii")] <- "other"
zooscan_all50 <- zooscan_all50[nrow(zooscan_all50):1, ]
zooscan_all50$objid <- NULL
zooscan_all50$taxon <- as.factor(zooscan_all50$taxon)
table(zooscan_all50$taxon)

# Note: these charts take too long to construct!!!
#chart(zooscan_all50, Acantharea_like ~ Actinopterygii_like %col=% taxon) +
#  geom_point()

#chart(zooscan_all50, Acartiidae_like ~ Actinopterygii_like %col=% taxon) +
#  geom_point()

#chart(zooscan_all50, log(area) ~ Actinopterygii_like %col=% taxon) +
#  geom_point()

library(randomForest)
rf <- randomForest(taxon ~ ., data = zooscan_all50, importance = TRUE)
rf
varImpPlot(rf)
# Conclusion: area is by far the most discriminant variable!

zooscan_all50$area <- NULL
rf <- randomForest(taxon ~ ., data = zooscan_all50, importance = TRUE)
rf
varImpPlot(rf)
# Now, without area, error rate is even larger... and Acantharea, Acartiidae,
# and Actinopterygii are not very well classified >80% error!

# A final comparison with only Hu moments
zooscan_all50$Acantharea_like <- NULL
zooscan_all50$Acartiidae_like <- NULL
zooscan_all50$Actinopterygii_like <- NULL
rf <- randomForest(taxon ~ ., data = zooscan_all50, importance = TRUE)
rf
varImpPlot(rf)

# Error rate is about the same, including for our 3 target classes => useless!



