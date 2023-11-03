library(ggplot2)
library(dplyr)
library(vegan)
library(emmeans) # used for tukey test for anovas with error block consideration
library(tidyverse) # used to transform data into long format for abiotic comparisons following PCA
library(gridExtra) # for ggplot grids

# Ensure you alter the path or set your working directory accordingly
setwd("~/University Files/Year 4/Semester 2/BVB304/Assessments/Manuscript/Statistical Analysis")

# Research Questions
# How does tree basal area, volume, trees/he, and total ground cover vary depending on the site?
#################################################
# Read the data from CSV file
dat1 <- read.csv("csvs/mtsum.csv", header=TRUE)
str(dat1)
# Convert forest_type to a factor
dat1$forest_type <- as.factor(dat1$forest_type)

# Cleaning the data (some rows reading as characters)
dat1$Average.of.basal_area <- as.numeric(gsub(",", "", dat1$Average.of.basal_area))
dat1$Average.of.trees_ha <- as.numeric(gsub(",", "", dat1$Average.of.trees_ha))

# Fit nested ANOVA models for basal area, volume, trees per hectare, 
model1 <- aov(Average.of.basal_area ~ forest_type + Error(site/block), data = dat1)
model2 <- aov(Average.of.trees_ha ~ forest_type + Error(site/block), data = dat1)
model3 <- aov(Average.of.tree_volume ~ forest_type + Error(site/block), data = dat1)
model4 <- aov(Average.of.total_groundcover ~ forest_type + Error(site/block), data = dat1)


# Compute emmeans and run tukey tests
BAem <- emmeans(model1, "forest_type", type = "response")
pairs(BAem) # No significance seen for basal area
THem <- emmeans(model2, "forest_type", type = "response")
pairs(THem) # No significance seen for trees per hectare
TVem <- emmeans(model3, "forest_type", type = "response")
pairs(TVem) # No significance seen for volume
GCem <- emmeans(model4, "forest_type", type = "response")
pairs(GCem) # No significance seen for average ground cover


# Average of Basal Area by Forest Type
p1 <- ggplot(dat1, aes(x = forest_type, y = Average.of.basal_area)) + 
  geom_boxplot(alpha=0.7) +
  labs(title = "Basal Area Across Sites", 
       x = "Forest Type", 
       y = expression("Average Basal Area (m"^2*")")) +
  theme_bw()

# Average of Trees per hectare by Forest Type
p2 <- ggplot(dat1, aes(x = forest_type, y = Average.of.trees_ha)) + 
  geom_boxplot(alpha=0.7) +
  labs(title = "Trees per hectare Across Sites", 
       x = "Forest Type", 
       y = "Average Trees per hectare") +
  theme_bw()

# Average of Tree Volume by Forest Type
p3 <- ggplot(dat1, aes(x = forest_type, y = Average.of.tree_volume)) + 
  geom_boxplot(alpha=0.7) +
  labs(title = "Tree Volume Across Sites", 
       x = "Forest Type", 
       y = expression("Average Tree Volume (m"^3*")")) +
  theme_bw()

# Average of Ground Cover Percentage by Forest Type
p4 <- ggplot(dat1, aes(x = forest_type, y = Average.of.total_groundcover)) + 
  geom_boxplot(alpha=0.7) +
  labs(title = "Ground Cover Across Sites", 
       x = "Forest Type", 
       y = "Average Ground Cover (%)") +
  theme_bw()

# Arrange plots in 2x2 grid
grid.arrange(p1, p2, p3, p4, ncol=2)

#################################################################################################

abiotic_csv <- read.csv("csvs/newpca.csv", header=T)
abiotic_csv$mean_light_flux <- as.numeric(gsub(",", "", abiotic_csv$mean_light_flux))
abiotic_only <- abiotic_csv[, c("litter_depth_mm", "mean_soil_pH_field", "mean_soil_pH_lab", "average_soil_water_percent", "canopy_cover", "mean_light_flux", "mean_rel_humidity", "mean_temp", "total_groundcover_percent")]
pairs(abiotic_only)

grouped_data <- abiotic_csv %>%
  group_by(site, block) %>%
  summarise(across(starts_with("litter"):total_groundcover_percent, mean, na.rm=TRUE))

abiotic_only <- as.data.frame(grouped_data[, 3:ncol(grouped_data)])
my.rda_new <- rda(abiotic_only, scale=TRUE)

biplot(my.rda_new, display = "sites", type = "points")
title("Principal Component Analysis (PCA) of Abiotic Conditions Across Four Sites")

# Colors for each site
site_colors <- c("green", "blue", "red", "purple")

# Loop to plot ellipses for each site with respective color
unique_sites <- unique(grouped_data$site)
for (i in 1:length(unique_sites)) {
  ordiellipse(my.rda_new, group = grouped_data$site, kind="sd", label=FALSE, conf=0.95, col=site_colors[i], lty=2, draw="polygon", show.groups=unique_sites[i])
}

text(-0.8, 1.4,cex = 1.3, "Site 1")
text(2.6, -0.5,cex = 1.3, "Site 2")
text(-0.22, -1.2,cex = 1.3, "Site 3")
text(-0.7, 0.13, cex = 1.3,"Site 4")

scores <- scores(my.rda_new, display="sites")
points(scores[,1], scores[,2], pch=16, col="black")
text(scores[,1], scores[,2], labels = grouped_data$block, pos=3, cex=1, col="black")

sp.scores <- scores(my.rda_new, display = "species")
text(sp.scores, labels = gsub("_", " ", rownames(sp.scores)), pos=3, cex=1)
arrows(0, 0, sp.scores[,1], sp.scores[,2], length = 0.1)

# View eigenvalues
eigenvals <- my.rda_new$CA$eig
print(eigenvals)

explained_var <- eigenvals / sum(eigenvals)
print(explained_var)

explained_var_percent <- explained_var * 100
print(explained_var_percent)


#################################################################################################

# Read the data from CSV file
abiotic_orig <- read.csv("csvs/mabiotic.csv", header=TRUE)
str(abiotic_orig)

# Need long data for light and humidity
data_long_light <- abiotic_orig %>%
  pivot_longer(
    cols = starts_with("Light"), 
    names_to = "Measurement",
    values_to = "Value")



# Fit nested ANOVA models, compute emmeans and run tukey tests for light, humidity, leaf litter, & field pH
model5 <- aov(Average.Light..Lux. ~ Site.Type + Error(Site), data = data_long_light)
LSem <- emmeans(model5, "Site.Type", type = "response")
pairs(LSem)
model6 <- aov(Average.Rel..Humidity.... ~ Site.Type + Error(Site), data = data_long_light)
HSem <- emmeans(model6, "Site.Type", type = "response")
pairs(HSem)
model7 <- aov(litter_depth_mm ~ forest_type + Error(site), data = abiotic_csv)
LSem <- emmeans(model7, "forest_type", type = "response")
pairs(LSem)
model8 <- aov(mean_soil_pH_field ~ forest_type + Error(site), data = abiotic_csv)
PSem <- emmeans(model8, "forest_type", type = "response")
pairs(PSem)

# Create each plot separately then feed to grid extra
p1 <- ggplot(data_long_light, aes(x=`Site.Type`, y=`Average.Light..Lux.`)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Light Intensity Across Sites",
       x = "Forest Type",
       y = "Light Intensity (Lux)") +
  theme_bw()

p2 <- ggplot(data_long_light, aes(x=`Site.Type`, y=`Average.Rel..Humidity....`)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Relative Humidity Across Sites",
       x = "Forest Type",
       y = "Relative Humidity (%)") +
  theme_bw()

p3 <- ggplot(abiotic_csv, aes(x=forest_type, y=`litter_depth_mm`)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Litter Depth Across Sites",
       x = "Forest Type",
       y = "Litter Depth (mm)") +
  theme_bw()

p4 <- ggplot(abiotic_csv, aes(x=forest_type, y=mean_soil_pH_field)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Field Soil pH Across Sites",
       x = "Forest Type",
       y = "Soil pH (Field)") +
  theme_bw()

# Arrange the plots in a 2x2 grid using grid.arrange()
grid.arrange(p1, p2, p3, p4, ncol=2)


#################################################################################################

## nMDS
dat3<-read.csv("csvs/mgcover.csv")

## Obtain only the numeric columns for the community data matrix
com <- dat3[, sapply(dat3, is.numeric)]

## Check for NA, NaN, or Inf values and handle them column-wise 
com <- data.frame(lapply(com, function(x) {
  if(any(is.na(x)) || any(is.nan(x)) || any(is.infinite(x))) {
    cat("Column contains NA, NaN, or Inf values. Handling them...\n")
    x[is.na(x)] <- 0  ## This is a simplistic way, consider better imputation methods as needed
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
  }
  return(x)
}))

m_com = as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds    
#extract NMDS scores (x and y coordinates)
## Extract NMDS scores
data.scores <- as.data.frame(scores(nmds)$sites)
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores$forest_type = dat3$forest_type

figureX = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( colour = forest_type))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "forest_type", y = "NMDS2", title="Non-Metric Multidimensional Scaling (NMDS) of Ground Cover Across Different Forest Types")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00", "blue", "red")) 

figureX

## permanova code for testing significance of the % cover--it is essential an anova but for a multivariate matrix of y response variables.

## Filter only numeric columns for distance matrix computation
numeric_data <- dat3[, sapply(dat3, is.numeric)]

## Compute the distance matrix
dist_dml <- vegdist(x = as.matrix(numeric_data), method = "bray", binary=FALSE, diag=TRUE, upper=TRUE, na.rm=FALSE)

## Perform PERMANOVA test
set.seed(12345)
y_permanova <- adonis2(dist_dml ~ forest_type, data = dat3, permutations = 999)
print(y_permanova)

# The p-value is less than 0.05 (p=0.001), indicating that there is a significant difference in % cover among the different forest types.
# You reject the null hypothesis of no difference among groups.

# Boxplots

# Fit nested ANOVA models, compute emmeans and run tukey tests for light, humidity, leaf litter, & field pH
model9 <- aov(leaf_litter_. ~ forest_type, data = dat3)
LSem <- emmeans(model9, "forest_type", type = "response")
pairs(LSem)
model10 <- aov(woody_debris_. ~ forest_type, data = dat3)
WSem <- emmeans(model10, "forest_type", type = "response")
pairs(WSem)
model11 <- aov(fern_cover_. ~ forest_type, data = dat3)
FSem <- emmeans(model11, "forest_type", type = "response")
pairs(FSem)
model12 <- aov(vine_cover_. ~ forest_type, data = dat3)
DSem <- emmeans(model12, "forest_type", type = "response")
pairs(DSem)

# Create each plot separately then feed to grid extra
p1 <- ggplot(dat3, aes(x=forest_type, y=leaf_litter_.)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Leaf Litter Cover Across Sites",
       x = "Forest Type",
       y = "Leaf Litter Cover (%)") +
  theme_bw()

p2 <- ggplot(dat3, aes(x=forest_type, y=woody_debris_.)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Woody Debris Cover Across Sites",
       x = "Forest Type",
       y = "Woody Debris Cover (%)") +
  theme_bw()

p3 <- ggplot(dat3, aes(x=forest_type, y=fern_cover_.)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Fern Cover Across Sites",
       x = "Forest Type",
       y = "Fern Cover (%)") +
  theme_bw()

p4 <- ggplot(dat3, aes(x=forest_type, y=vine_cover_.)) + 
  geom_boxplot(alpha=0.7) +
  labs(title="Vine Cover Across Sites",
       x = "Forest Type",
       y = "Vine Cover (%)") +
  theme_bw()

# Arrange the plots in a 2x2 grid using grid.arrange()
grid.arrange(p1, p2, p3, p4, ncol=2)

#########################################################################################################



### Shannons Diversity code

# Create dataframe for seedlings providing calculations already performed
seedling_data <- data.frame(
  Site = c(1, 2, 3, 4),
  Diversity_Index = c(0.6931, 0.8951, 1.2059, 1.6094))
seedling_data$Site <- factor(seedling_data$Site, levels = c(1, 2, 3, 4), 
                             labels = c("21 year old", "10 year old", "5 year old", "secondary"))

# Calculate shannons DI for using tree families data set
tree_family_data <- read.csv("csvs/families.csv", header=T)
# Calculate proportions for each species
tree_family_data %>%
  group_by(forest_type) %>%
  summarise_at(vars(-block), ~./sum(.)) -> prop_data

print(prop_data)
# Calculate Shannon Index
shannon_index <- function(row) {
  species <- row[!is.nan(row) & row > 0] # Filter out NaN and zero proportions
  if(length(species) == 0) return(NA) # Return NA if no species present
  return(-sum(species * log(species))) # Calculate Shannon Index using natural log
}


shannon_values <- apply(prop_data[,-1], 1, shannon_index)

final_data <- data.frame(forest_type = prop_data$forest_type, shannon_index = shannon_values)

final_data_grouped <- final_data %>%
                        group_by(forest_type) %>%
                          summarise(avg_shannon_index = mean(shannon_index, na.rm = TRUE)) # Use na.rm to remove NA values if any

# combine the data and create a nice single plot with both indices
tree_family_data <- data.frame(
  Site = c("21 year old", "10 year old", "5 year old", "secondary"),
  Diversity_Index = c(1.03, 0.963, 1.36, 1.50),
  Type = "Tree Family Diversity"
)

# Update seedling_data to include a column indicating it's seedling
seedling_data$Type <- "Seedling Species Diversity"

# Combine data
combined_data <- rbind(seedling_data, tree_family_data)

# Create the plot
ggplot(combined_data, aes(x = Site, y = Diversity_Index, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.9) +
  labs(
    title = "Shannon's Diversity Indices across Forest Sites",
    x = "Forest Type",
    y = "Shannon Diversity Index",
    fill = "" # Remove legend title
  ) +
  scale_fill_manual(values = c("Tree Family Diversity" = "sienna", "Seedling Species Diversity" = "olivedrab")) + # Adding the alpha parameter to make colors paler
  theme_minimal() +
  theme(legend.position = "bottom")

# Omg its over. FREEEEEDDOOMMMMMMM.