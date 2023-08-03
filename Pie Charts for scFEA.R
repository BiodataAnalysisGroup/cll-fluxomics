library(Seurat)
library(tidyseurat)
library(scales)
library(ggplot2)

# Aim: To visualize in pies the various metabolic clusters from scFEA in the CLL PBMCs

# 1. Get the proper Seurat object ready:

# Day 0
DimPlot(CLL_d0_scFEA, label = T)
Idents(CLL_d0_scFEA) <- CLL_d0_scFEA$FLUX_snn_res.0.3
new.cluster.ids <- c("Low metabolism", "High metabolism", "Urea", "Serine", "High glycolysis", "High glutamine",
                     "High Serine", "High pyrimidine", "High aa & pyrimidine")
names(new.cluster.ids) <- levels(CLL_d0_scFEA)
CLL_d0_scFEA<- RenameIdents(CLL_d0_scFEA, new.cluster.ids)
DimPlot(CLL_d0_scFEA, label = F)
CLL_d0_scFEA$Idents <- Idents(CLL_d0_scFEA)
# Get neccesary metadata
SO <- CLL_d0_scFEA
df <- SO$FLUX_snn_res.0.3
subset_df <- df[, c("labels", "FLUX_snn_res.0.3")]

#Start subsetting to useful PBMC populations (B cells, CD4+ T cells, CD8+ T cells)
Bcell <- subset_df[subset_df$labels == "B cells",]
CD4T <- subset_df[subset_df$labels == "CD4+ T cells",]
CD8T <- subset_df[subset_df$labels == "CD8+ T cells",]
Monocytes <- subset_df[subset_df$labels == "Monocytes",]

# CLL B cells
# Get the distinct values in the column
distinct_values <- unique(Bcell$FLUX_snn_res.0.3)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(Bcell$FLUX_snn_res.0.3 == value)
  total_observations <- length(Bcell$FLUX_snn_res.0.3)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("Low metabolism", "High metabolism", "Urea", "Serine", "High glycolysis", "High glutamine",
            "High Serine", "High pyrimidine", "High aa & pyrimidine")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("greenyellow", "deepskyblue", "turquoise", "indianred", "yellow", "tan1",
            "maroon1", "midnightblue", "red")

# Title for the pie chart
title <- "CLL B cells, Day 0, Pre-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

#-----------------------------------------

# CD4 T cells
# Get the distinct values in the column
distinct_values <- unique(CD4T$FLUX_snn_res.0.3)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(CD4T$FLUX_snn_res.0.3 == value)
  total_observations <- length(CD4T$FLUX_snn_res.0.3)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("Low metabolism", "High metabolism", "Urea", "Serine", "High glycolysis", 
            #"High glutamine",
            "High Serine")
#"High pyrimidine", "High aa & pyrimidine")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("greenyellow", "deepskyblue", "turquoise", "indianred", "yellow",
            "maroon1")

# Title for the pie chart
title <- "CD4 T cells, Day 0, Pre-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

# Convert the pie chart to a ggplot2 object
pie_ggplot <- ggplot(data.frame(data$values, labels), aes(x = "", y = data$values, fill = labels)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void()
#-------------------Monocytes------------------------------------

# Monocytes
# Get the distinct values in the column
distinct_values <- unique(Monocytes$FLUX_snn_res.0.3)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(Monocytes$FLUX_snn_res.0.3 == value)
  total_observations <- length(Monocytes$FLUX_snn_res.0.3)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("Low metabolism", "High metabolism", "Urea", "Serine", "High glycolysis", "High glutamine",
            "High Serine", "High pyrimidine", "High aa & pyrimidine")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("greenyellow", "deepskyblue", "turquoise", "indianred", "yellow", "tan1",
            "maroon1", "midnightblue", "red")

# Title for the pie chart
title <- "Monocytes, Day 0, Pre-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

#-------------------------CD8 T cells

# CD8 T cells
# Get the distinct values in the column
distinct_values <- unique(CD8T$FLUX_snn_res.0.3)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(CD8T$FLUX_snn_res.0.3 == value)
  total_observations <- length(CD8T$FLUX_snn_res.0.3)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("Low metabolism", "High metabolism", "Urea", "Serine", "High glycolysis", "High glutamine",
            "High Serine", "High pyrimidine", "High aa & pyrimidine")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("greenyellow", "deepskyblue", "turquoise", "indianred", "yellow", "tan1",
            "maroon1", "midnightblue", "red")

# Title for the pie chart
title <- "CD8 T cells, Day 0, Pre-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Day 30 after Ibrutinib

# Day 30
# DimPlot(CLL_d0_scFEA, label = T)
# Idents(CLL_d0_scFEA) <- CLL_d0_scFEA$FLUX_snn_res.0.3
# new.cluster.ids <- c("Low metabolism", "High metabolism", "Urea", "Serine", "High glycolysis", "High glutamine",
#                      "High Serine", "High pyrimidine", "High aa & pyrimidine")
# names(new.cluster.ids) <- levels(CLL_d0_scFEA)
# CLL_d0_scFEA<- RenameIdents(CLL_d0_scFEA, new.cluster.ids)
# DimPlot(CLL_d0_scFEA, label = F)
# CLL_d0_scFEA$Idents <- Idents(CLL_d0_scFEA)
# Get neccesary metadata
SO <- CLL_d30_scFEA
df <- SO@meta.data
subset_df <- df[, c("labels", "FLUX_snn_res.0.2")]

#Start subsetting to useful PBMC populations (B cells, CD4+ T cells, CD8+ T cells)
Bcell <- subset_df[subset_df$labels == "B cells",]
CD4T <- subset_df[subset_df$labels == "CD4+ T cells",]
CD8T <- subset_df[subset_df$labels == "CD8+ T cells",]
Monocytes <- subset_df[subset_df$labels == "Monocytes",]

# CLL B cells
# Get the distinct values in the column
distinct_values <- unique(Bcell$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(Bcell$FLUX_snn_res.0.2 == value)
  total_observations <- length(Bcell$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", "Low metabolism",  "High metabolism & 2OG/malate", "High Aspartate", "High Glutamine/Glutamate", "High pyrimidine", "High Aspartate & purine synthesis")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow", "deepskyblue4", "hotpink", "tan1", "midnightblue", "purple")

# Title for the pie chart
title <- "CLL B cells, Day 30, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

#-----------------------------------------

# CD4 T cells
# Get the distinct values in the column
distinct_values <- unique(CD4T$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(CD4T$FLUX_snn_res.0.2 == value)
  total_observations <- length(CD4T$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", "Low metabolism",  "High metabolism & 2OG/malate", "High Aspartate", "High Glutamine/Glutamate", "High pyrimidine", "High Aspartate & purine synthesis")

#"High pyrimidine", "High aa & pyrimidine")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow", "deepskyblue4", "hotpink", "tan1", "midnightblue", "purple")


# Title for the pie chart
title <- "CD4 T cells, Day 30, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

# Convert the pie chart to a ggplot2 object
pie_ggplot <- ggplot(data.frame(data$values, labels), aes(x = "", y = data$values, fill = labels)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void()
#-------------------Monocytes------------------------------------

# Monocytes
# Get the distinct values in the column
distinct_values <- unique(Monocytes$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(Monocytes$FLUX_snn_res.0.2 == value)
  total_observations <- length(Monocytes$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", "Low metabolism",  "High metabolism & 2OG/malate", "High Aspartate", "High Glutamine/Glutamate", "High pyrimidine", "High Aspartate & purine synthesis")

#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow", "deepskyblue4", "hotpink", "tan1", "midnightblue", "purple")


# Title for the pie chart
title <- "Monocytes, Day 30, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

#-------------------------CD8 T cells

# CD8 T cells
# Get the distinct values in the column
distinct_values <- unique(CD8T$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(CD8T$FLUX_snn_res.0.2 == value)
  total_observations <- length(CD8T$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", "Low metabolism",  "High metabolism & 2OG/malate", "High Aspartate", "High Glutamine/Glutamate", "High pyrimidine", "High Aspartate & purine synthesis")

#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow", "deepskyblue4", "hotpink", "tan1", "midnightblue", "purple")


# Title for the pie chart
title <- "CD8 T cells, Day 30, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))




#####------------------------120 days post Ibrutinib----------------------------

SO <- CLL_d120_scFEA
df <- SO@meta.data
subset_df <- df[, c("labels", "FLUX_snn_res.0.2")]

#Start subsetting to useful PBMC populations (B cells, CD4+ T cells, CD8+ T cells)
Bcell <- subset_df[subset_df$labels == "B cells",]
CD4T <- subset_df[subset_df$labels == "CD4+ T cells",]
CD8T <- subset_df[subset_df$labels == "CD8+ T cells",]
Monocytes <- subset_df[subset_df$labels == "Monocytes",]

# CLL B cells
# Get the distinct values in the column
distinct_values <- unique(Bcell$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(Bcell$FLUX_snn_res.0.2 == value)
  total_observations <- length(Bcell$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", "Low metabolism",  "High nucleotide synthesis", 
            "TCA/2OG/malate_glycolysis_down", "High Aspartate/Purine synthesis", 
            "High hydroxyproline", "High Glutamine", "High Serine metabolism")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow",  "plum4", "deepskyblue4","hotpink", "khaki", "tan1",  "lightgreen")

# Title for the pie chart
title <- "CLL B cells, Day 120, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

#-----------------------------------------

# CD4 T cells
# Get the distinct values in the column
distinct_values <- unique(CD4T$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(CD4T$FLUX_snn_res.0.2 == value)
  total_observations <- length(CD4T$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", 
            #"Low metabolism",  
            "High nucleotide synthesis", 
            #"TCA/2OG/malate_glycolysis_down", 
            "High Aspartate/Purine synthesis") 
            #"High hydroxyproline", "High Glutamine", "High Serine metabolism")

#"High pyrimidine", "High aa & pyrimidine")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1",  "plum4", "hotpink")


# Title for the pie chart
title <- "CD4 T cells, Day 120, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

# Convert the pie chart to a ggplot2 object
pie_ggplot <- ggplot(data.frame(data$values, labels), aes(x = "", y = data$values, fill = labels)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void()
#-------------------Monocytes------------------------------------

# Monocytes
# Get the distinct values in the column
distinct_values <- unique(Monocytes$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(Monocytes$FLUX_snn_res.0.2 == value)
  total_observations <- length(Monocytes$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", 
            "Low metabolism",  
            #"High nucleotide synthesis", 
            "TCA/2OG/malate_glycolysis_down", 
            "High Aspartate/Purine synthesis", 
#"High hydroxyproline", "High Glutamine", 
"High Serine metabolism")

#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow",   "deepskyblue4","hotpink",   "lightgreen")


# Title for the pie chart
title <- "Monocytes, Day 120, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))

#-------------------------CD8 T cells

# CD8 T cells
# Get the distinct values in the column
distinct_values <- unique(CD8T$FLUX_snn_res.0.2)

# Initialize an empty list to store the percentages
percentages <- list()

# Loop through each distinct value and calculate its percentage
for (value in distinct_values) {
  value_count <- sum(CD8T$FLUX_snn_res.0.2 == value)
  total_observations <- length(CD8T$FLUX_snn_res.0.2)
  percentage_value <- (value_count / total_observations) * 100
  percentages[[as.character(value)]] <- percentage_value
}

percent <- as.data.frame(percentages)


# Pie plotting
data <- percent[1,]
data <- t(data)
data <- as.data.frame(data)
colnames(data) <- "values"
data$clusters <- rownames(data)
data <- data[order(data$clusters), ]


labels <- c("High metabolism with Glycolysis", "Low metabolism",  "High nucleotide synthesis", 
            "TCA/2OG/malate_glycolysis_down", "High Aspartate/Purine synthesis", 
            "High hydroxyproline", "High Glutamine", "High Serine metabolism")
#labels <- rownames(data)

# Remove dots from the labels using gsub()
#labels <- gsub("\\.", " ", labels)

# Custom colors for each slice
colors <- c("dodgerblue1", "greenyellow",  "plum4", "deepskyblue4","hotpink", "khaki", "tan1",  "lightgreen")


# Title for the pie chart
title <- "CD8 T cells, Day 120, Post-Ibrutinib"

# Creating the pie chart with labels and custom colors

pie(data$values, main = title, col = colors, cex = 0.9, cex.lab = 5, labels = sprintf("%s (%s)", labels, percent(data$values/sum(data$values))))
